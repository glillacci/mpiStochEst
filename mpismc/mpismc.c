/*
 *  mpismc.c
 *
 *  SMC sampler for INSIGHT.
 *  This program propagates a population of particles from one tolerance to the next.
 *
 *
 *  This file is part of mpiStochEst.
 *  Copyright 2011-2017 Gabriele Lillacci.
 *
 *  mpiStochEst is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  mpiStochEst is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with mpiStochEst.  If not, see <http://www.gnu.org/licenses/>.
 */


// Standard C library includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Mpi include
#include <mpi.h>

// Libxml2 includes
#include <libxml2/libxml/xmlversion.h>
#include <libxml2/libxml/xmlstring.h>
#include <libxml2/libxml/tree.h>
#include <libxml2/libxml/parser.h>
#include <libxml2/libxml/xpath.h>

// Custom library includes
#include <parestlib.h>
#include <stochmod.h>

// Global variable declarations
stochmod * model = NULL;			// The model to estimate
gsl_matrix * data = NULL;			// The flow cytometry data
gsl_matrix * out = NULL; 			// The output matrix of the model
gsl_vector * times = NULL;			// The time points
gsl_matrix * bg = NULL;				// The background distribution
gsl_matrix * prev = NULL;			// The previous population of particles
gsl_vector * weights = NULL;		// The weights of the previous population
gsl_matrix * cov = NULL; 			// The covariance matrix of the previous population
double * lb, * ub;					// Box constraints on the parameters
size_t M;							// Number of experimental samples to use for identification
size_t P;							// Number of measured species
size_t K;							// Number of time points
size_t R;							// Number of parameters
size_t Z; 							// Number of inputs
double epsilon;						// Distance tolerance

// Tags for message passing
#define WORKTAG 1
#define DIETAG 0

// Enumeration for the data mode
enum DATA_MODE_ENUM
{
	DATA_MODE_FCS = 0,
	DATA_MODE_TEXT = 1,
} DATA_MODE;

// Enumeration for the data scale
enum DATA_SCALE_ENUM
{
	DATA_SCALE_LOG = 0,
	DATA_SCALE_LINEAR = 1,
} DATA_SCALE;

// Functions for parallel simulation
static void master (int argc, char * argv[]);
static void slave ();

// Utility functions
double dunif (double x, double a, double b);
double dmvnorm (const gsl_vector * x, const gsl_vector * mean, const gsl_matrix * var);
double rnorm (const gsl_rng * r, double m, double v);
double runif (const gsl_rng * r, double a, double b);
int rmvnorm (const gsl_rng * r, const gsl_vector * mean, const gsl_matrix * var, double result[]);
double rperturb (const gsl_rng * r, double theta, double wd);
double dperturb (double x, double a, double b);
int mvrperturb (const gsl_rng * r, const double theta[], const gsl_matrix * Sig, double result[]);
int covmat (const gsl_matrix * dat, const gsl_vector * weights, gsl_matrix * Sig);

double sample_fl (const gsl_rng * r, double z, double mufp, double sigfp, gsl_matrix * bg);

// Functions for problem definition processing
int processPriorData (xmlDocPtr doc, double lb[], double ub[]);
int processDataSet (xmlDocPtr doc, gsl_vector * times, gsl_matrix * data);
int processBackground (xmlDocPtr doc, gsl_matrix * bg);
int processModel (xmlDocPtr doc, stochmod * model);


int main (int argc, char * argv[])
{
	// ==================================== INITIALIZE LIBRARIES =============================================
	// Initialize MPI
	MPI_Init (&argc, &argv);
	// Initialize libxml2
	LIBXML_TEST_VERSION

	// =============================== PROCESS COMMAND LINE ARGUMENTS ========================================
	// Find out my identity in the default communicator
	int myrank;
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

	// Check correctness of command line arguments
	if (argc != 7)
	{
		if (myrank == 0)
			fprintf (stderr, "usage: %s problem-file-path old-population-file-path new-population-file-path nSMC S print-every\n\n", argv[0]);
		return EXIT_FAILURE;
	}

	// Put command line arguments into correct variables
	char * problemfile = argv[1];
	char * oldpopfile = argv[2];
	char * newpopfile = argv[3];
	size_t nsim = atoi (argv[4]);


	// =============================== EXTRACT INFO FROM PROBLEM FILE ========================================
	// Parse the problem file
	xmlDocPtr problem = xmlReadFile (problemfile, NULL, 0);
	if (problem == NULL)
	{
		fprintf (stderr, "error reading problem file\n");
		return EXIT_FAILURE;
	}

	// Set up the stochasitc model to estimate
	model = (stochmod *) malloc (sizeof (stochmod));
	if (!model || (processModel (problem, model) != GSL_SUCCESS))
	{
		fprintf (stderr, "error in setting up the stochastic model to estimate\n");
		return EXIT_FAILURE;
	}

	// Initialize XPath evaluation context
	xmlXPathContextPtr ctxt = xmlXPathNewContext (problem);

	// Initialize a xmlXPathObjectPtr object
	xmlXPathObjectPtr xobj;

	// Use XPath to find the number of samples to use for identification
	xobj = xmlXPathEvalExpression ((xmlChar *) "/problem/samples", ctxt);
	M = (size_t) xmlXPathCastToNumber (xobj);
	xmlXPathFreeObject (xobj);

	// Use XPath to find the number of measured species
	xobj = xmlXPathEvalExpression ((xmlChar *) "/problem/outputs", ctxt);
	P = (size_t) xmlXPathCastToNumber (xobj);
	xmlXPathFreeObject (xobj);

	// Use XPath to find the number of time points
	xobj = xmlXPathEvalExpression ((xmlChar *) "/problem/timepoints", ctxt);
	K = (size_t) xmlXPathCastToNumber (xobj);
	xmlXPathFreeObject (xobj);

	// Use XPath to find the number of parameters
	xobj = xmlXPathEvalExpression ((xmlChar *) "/problem/parameters", ctxt);
	R = (size_t) xmlXPathCastToNumber (xobj);
	xmlXPathFreeObject (xobj);

	// Use the results to initialize some useful variables
	lb = (double *) malloc ((2*P + R) * sizeof(double));
	ub = (double *) malloc ((2*P + R) * sizeof(double));
	if (lb == NULL || ub == NULL)
	{
		fprintf (stderr, "error in memory allocation\n");
		return EXIT_FAILURE;
	}
	data = gsl_matrix_alloc (M, K*P);
	times = gsl_vector_alloc (K);
	bg = gsl_matrix_alloc (M, P);

	// Process the prior information
	if (processPriorData (problem, lb, ub) != GSL_SUCCESS)
	{
		fprintf (stderr, "error: the problem definition file does not contain valid prior information\n");
		return EXIT_FAILURE;
	}

	// Load the data for the identification
	if (processDataSet (problem, times, data) != GSL_SUCCESS)
	{
		fprintf (stderr, "error: could not process data set information\n");
		return EXIT_FAILURE;
	}

	// Load the background distribution
	if (processBackground (problem, bg) != GSL_SUCCESS)
	{
		fprintf (stderr, "error: could not process background information\n");
		return EXIT_FAILURE;
	}
	// Check whether the bg was actually loaded or it is all zero
	if (gsl_matrix_isnull (bg))
	{
		gsl_matrix_free (bg);
		bg = NULL;
	}

	// We are now done with the problem definition file, so we can free it
	xmlFreeDoc (problem);

	// Set up output matrix
	out = gsl_matrix_calloc (model->nout, model->nspecies);
	model->output (out);

	// Load previous population of particles
	prev = gsl_matrix_alloc (nsim, 2*model->nout + model->nparams);
	weights = gsl_vector_alloc (nsim);
	cov = gsl_matrix_alloc (2*model->nout + model->nparams, 2*model->nout + model->nparams);
	if (strstr (oldpopfile, "NULL"))
	{
		gsl_matrix_set_zero (prev);
		gsl_vector_set_zero (weights);
		gsl_matrix_set_zero (cov);
		if (myrank == 0)
			printf ("This is the sampling of the initial population.\n");
	}
	else
	{
		// Load the previous population of particles
		FILE * pfile = fopen (oldpopfile, "r");
		if (gsl_matrix_fscanf (pfile, prev) != GSL_SUCCESS)
		{
			fprintf (stderr, "error: could not load previous particle population\n");
			return EXIT_FAILURE;
		}
		fclose (pfile);
		// Load the previous weights
		char wfilename[500];
		sprintf (wfilename, "%s_weights", oldpopfile);
		FILE * wfile = fopen (wfilename, "r");
		if (gsl_vector_fscanf (wfile, weights) != GSL_SUCCESS)
		{
			fprintf (stderr, "error: could not load previous weights\n");
			return EXIT_FAILURE;
		}
		fclose (wfile);
		// Renormalize the weights
		double sow = gsl_blas_dasum (weights);
		gsl_vector_scale (weights, 1/sow);
		// Compute the covariance matrix
		covmat (prev, weights, cov);
	}

	// ============================ RUN THE PARALLEL COMPUTATION ==================================
	// Inform the user that we are ready to start
	if (myrank == 0)
	{
		printf ("\nThis is mpiSmc v1.0\nI am the master process\n\n");
		printf ("Problem file path:\n%s\n", problemfile);
		printf ("The new particle population will be saved in:\n%s\n", newpopfile);
		printf ("Will estimate %s using %d SMC samples.\n\n", model->name, (int) nsim);
	}

	// Start a new clock to keep track of elapsed time
	clock_t tic = clock ();

	// Decide what I have to do
	if (myrank == 0)
	{
		// I am in charge
		master (argc, argv);
	}
	else
	{
		// I have to do what I am told
		slave ();
	}

	// Check how much time has passed
	clock_t toc = clock ();
	if (myrank == 0)
		printf ("\n\nSMC inference complete.\nElapsed time: %g sec.\n\n", (((double) (toc - tic)) / CLOCKS_PER_SEC));


	// ===================================== CLEAN UP =============================================
	// Free objects
	free (model);
	free (lb);
	free (ub);
	gsl_matrix_free (out);
	gsl_matrix_free (data);
	gsl_vector_free (times);
	gsl_matrix_free (prev);
	gsl_vector_free (weights);
	gsl_matrix_free (cov);
	gsl_matrix_free (bg);
	xmlXPathFreeContext (ctxt);

	// Shut down MPI
	MPI_Finalize ();

	// Shut down libxml2
	xmlCleanupParser ();
	xmlMemoryDump ();

	// End program
	return EXIT_SUCCESS;
}


// Master process tasks
static void master (int argc, char * argv[])
{
	// Declare some necessary variables
	int ntasks, rank, count = 0;
	double work;
	MPI_Status status;
	size_t res_no = (model->nspecies)*(times->size);
	double results[res_no];
	size_t nsam = atoi (argv[4]);
	size_t print_every = atoi (argv[6]);

	// Find out how many processes there are in the default communicator
	MPI_Comm_size (MPI_COMM_WORLD, &ntasks);

	// Set up the random number generator
	gsl_rng_env_setup ();
	gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (r, (unsigned long int) time (NULL));

	// Set metric parameters
	double beta = 0.05; 						// "External" confidence level
	double alpha = 1 - sqrt (1 - beta);			// "Internal" condifence level
	size_t S = atoi(argv[5]);					// Number of simulations


	// Tell the user about settings related to the data
	switch (DATA_MODE)
	{
	case DATA_MODE_FCS:
		printf ("\nThe data mode is: FCS.\n");
		break;
	case DATA_MODE_TEXT:
		printf ("\nThe data mode is: TEXT.\n");
		break;
	}
	switch (DATA_SCALE)
	{
	case DATA_SCALE_LOG:
		printf ("The data has been transformed to LOG SCALE.\n");
		break;
	case DATA_SCALE_LINEAR:
		printf ("The data will be treated as LINEAR SCALE.\n");
		break;
	}
	if (bg)
		printf ("A BACKGROUND distribution has been specified.\n");
	else
		printf ("No BACKGROUND distribution has been supplied.\n");


	// DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
	// Save the data matrix that was read to file for debugging purposes
	FILE * datadebugfile = fopen ("data_debug.txt", "w");
	if ((datadebugfile == NULL) || (gsl_matrix_fprintf (datadebugfile, data, "%e") != GSL_SUCCESS))
	{
		printf (">> error saving data debug information.. program exiting\n");
		return;
	}
	fclose (datadebugfile);
	printf ("Information for data matrix debugging saved to file.\n");

	// Determine if we have too many tasks for the number of simulations
	if (ntasks > (S+1))
	{
		// The extra processes need to be killed to avoid messes in the GRID
		for (size_t rank = (S+1); rank < ntasks; rank++)
		{
			// Send the "DIE" message to the slave
			MPI_Send (0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
		}

		// Update ntasks
		ntasks = S+1;

		printf ("\n\nWARNING: some extra tasks were killed. There are %d processes left.\n\n", (int) ntasks);
	}

	// Compute critical epsilon
	epsilon = kolmogorov_cdf_inverse (S, 1 - alpha) + 0.005224; // kolmogorov_cdf_inverse (80000, 1-alpha);

	// Visualize metric parameters
	printf ("\nThe metric parameters are as follows\n");
	printf ("Confidence level beta = %f\n", beta);
	printf ("Confidence level alpha = %f\n", alpha);
	printf ("Number of simulations S = %d\n", (int) S);
	printf ("Number of experimental samples M = %d\n", (int) M);
	printf ("Distance tolerance epsilon = %f\n\n", epsilon);

	// Scale covariance matrix by epsilon
	gsl_matrix_scale (cov, epsilon);

	// Initialize parameter matrix
	gsl_matrix * x = gsl_matrix_calloc (nsam, 2*model->nout + model->nparams);
	// Initialize weights vector
	gsl_vector * wt = gsl_vector_alloc (nsam);

	// Initialize distances
	double dc[data->size2];
	int dcg;

	// Initialize parameter array
	double thetac[2*model->nout + model->nparams + model->nin];

	/*
	 * thetac has the following structure
	 * [mufp_1 sigfp_1 ... mufp_P sigfp_P theta_1 ... theta_R u_1 ... u_Z]
	 */

	// Set the value of the input
	thetac[2*P + R + 0] = 20.0;

	// Allocate rxn_ensemble object to store SCRN simulations
	rxn_ensemble * res = rxn_ensemble_alloc (S, model->nspecies, times->size);

	// Allocate matrix to store SCRN output
	gsl_matrix * counts = gsl_matrix_calloc (res->nreplic, model->nout*res->ntimes);

	// Allocate drow and crow vectors
	gsl_vector_view drow, crow;

	// Start a new clock to provide estimate of remaining time
	clock_t tic2 = clock ();

	// Main SMC loop
	size_t i = 0;
	size_t ct = 0;
	while (i < nsam)
	{
		// Provide info on expected completion time
		if ((ct % 500 == 0) && (i != 0))
		{
			printf ("\n\nParticle acceptance rate: %f\n", ((double) i)/((double) ct));
			// Compute estimate of remaining time
			clock_t toc2 = clock ();
			double sec = ((double) (toc2 - tic2) / CLOCKS_PER_SEC);
			double eta = sec/((double) i)*(nsam - i);
			int hrs = floor (eta/3600);
			eta = eta - hrs*3600;
			int min = floor (eta/60);
			eta = eta - min*60;
			printf ("Estimated time remaining: %d hrs %d min %f sec\n\n", hrs, min, eta);
		}

		// Warn the user if we have a very low acceptance rate
		if ((ct % 1000 == 0) && (i == 0) && (ct > 0))
		{
			printf ("\nWarning: %d iterations without an acceptance!\n", (int) ct);
		}

		// Resample a particle from the previous population or the prior
		if (gsl_matrix_isnull (prev))
		{
			// Sample from the prior
			for (size_t j = 0; j < 2*model->nout + model->nparams; j++)
			{
				thetac[j] = runif (r, lb[j], ub[j]);
			}
		}
		else
		{
			// Sample from the previous population
			// Generate a uniform random number
			double u = gsl_rng_uniform (r);

			// Accumulate the weights until they are larger than u
			double mu = 0.0;
			size_t ix = 0;
			while ((mu < u) && (ix < (nsam-1)))
			{
				mu += gsl_vector_get (weights, ix);
				ix++;
			}

			// Take the particle corresponding to ix
			for (size_t j = 0; j < 2*model->nout + model->nparams; j++)
			{
				thetac[j] = gsl_matrix_get (prev, ix, j);
			}

			// Introduce perturbation but without violating the prior bounds
			/*
			// Multivariate normal perturbation
			double thetass[2*model->nout + model->nparams];
			int ok = 0;
			while (!ok)
			{
				mvrperturb (r, thetac, cov, thetass);
				ok = 1;
				for (size_t j = 0; j < 2*model->nout + model->nparams; j++)
				{
					ok = ok * ((thetass[j] > lb[j] && thetass[j] < ub [j]) ? 1 : 0);
				}
			}
			*/

			// Independent uniform perturbation
			double thetass[2*model->nout + model->nparams];
			for (size_t j = 0; j < 2*model->nout + model->nparams; j++)
			{
				// If the lower and upper bound are both zero...
				if (lb[j] == 0.0 && ub[j] == 0.0)
				{
					// ... then the user wants this parameter to be fixed to zero.
					thetass[j] = 0.0;
				}
				else
				{
					// Otherwise proceed and introduce perturbation
					do
					{
						thetass[j] = rperturb (r, thetac[j], 2*epsilon);
					}
					while (thetass[j] <= lb[j] || thetass[j] >= ub[j]);
				}
			}

			// Thetass is now ok, replace thetac with thetass
			for (size_t j = 0; j < 2*model->nout + model->nparams; j++)
			{
				thetac[j] = thetass[j];
			}
		}

		/*
		// DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
		// Visualize the particle
		for (size_t k = 0; k < 2*model->nout + model->nparams + model->nin; k++)
		{
			printf("%f,\n", thetac[k]);
		}
		printf("\n");
		*/

		// Simulate SCRN with resampled values
		{
			// Pack the message for the slaves that contains the model parameters and inputs
			double msg[model->nparams + model->nin];
			for (size_t k = 0; k < model->nparams + model->nin; k++)
			{
				msg[k] = thetac[2*P+k];
			}

/*
			// DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
			// Visualize the message
			for (size_t k = 0; k < model->nparams + model->nin; k++)
			{
				printf("%f,\n", msg[k]);
			}
			printf("\n");
*/

			// Request simulations from the slaves by sending the parameters to simulate
			for (rank = 1; rank < ntasks; rank++)
			{
				MPI_Send (msg, 								// current parameters
						model->nparams + model->nin,			// # of data items
						MPI_DOUBLE,							// data items are of type double
						rank,									// destination process rank
						WORKTAG + count,						// user chosen message tag
						MPI_COMM_WORLD);						// default communicator

				count++;
			}

			// Keep requesting simulations until we have enough
			while (count < S)
			{
				// Receive results from a slave
				MPI_Recv (results,      	  					// message buffer
						res_no,             					// # of data items
						MPI_DOUBLE,        				 	// of type double
						MPI_ANY_SOURCE,      					// receive from any sender
						MPI_ANY_TAG,     	   					// any type of message
						MPI_COMM_WORLD,      					// default communicator
						&status);            					// info about the received message

				// Process the received results
				for (size_t k = 0; k < res_no; k++)
				{
					// The TAG from the sender is simulation #
					res->data[status.MPI_TAG]->counts->data[k] = results[k];
				}

				// Send the slave a new work unit
				MPI_Send (msg,								// message buffer
						model->nparams + model->nin,			// # of data items
						MPI_DOUBLE,							// data items are of type double
						status.MPI_SOURCE,					// destination process rank (the slave we just received from)
						WORKTAG + count,						// user chosen message tag
						MPI_COMM_WORLD);						// default communicator

				// Get the next unit of work to be done
				count++;
			}

			// We now have enough simulations, so receive all the outstanding results from the slaves
			for (rank = 1; rank < ntasks; rank++)
			{
				MPI_Recv (results, res_no, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

				// Process the results
				for (size_t k = 0; k < res_no; k++)
				{
					res->data[status.MPI_TAG]->counts->data[k] = results[k];
				}
			}

			// Reset count
			count = 0;
		}

		// Evaluate DCMS metric
		{
			// Extract the output
			rxn_ensemble_counts (res, counts, out);

			// Compute the fluorescence levels
			for (size_t k = 0; k < counts->size1; k++)
			{
				for (size_t l = 0; l < counts->size2; l++)
				{
					// Current count
					double z = gsl_matrix_get (counts, k, l);

					// Generate a fluorescence level
					switch (DATA_SCALE)
					{
					case DATA_SCALE_LOG:
						gsl_matrix_set (counts, k, l, log10 (sample_fl (r, z, thetac[0], thetac[1], bg)));
						break;
					case DATA_SCALE_LINEAR:
						gsl_matrix_set (counts, k, l, sample_fl (r, z, thetac[0], thetac[1], bg));
						break;
					}
				}
			}

			dcg = 1;

			// Compute the Kolmogorov distances
			for (size_t k = 0; k < data->size2; k++)
			{
				// Extract the appropriate column in the data matrix
				// ATT!! This vector is ALREADY SORTED in ascending order!
				drow = gsl_matrix_column (data, k);
				// Extract the appropriate column in the simulation matrix and sort it
				crow = gsl_matrix_column (counts, k);
				gsl_sort_vector (&crow.vector);

				// Compute the distance
				dc[k] = ksdist_two_sample_nr_presort (&crow.vector, &drow.vector);

				//printf("%f, ", dc[k]);

				if (dc[k] > epsilon)
				{
					dcg = 0;
					break; // If one time point fails, the whole particle fails.
				}
			}
			//printf ("\n\n");
		}

		// Decide what to do
		if (dcg != 0)
		{
			// Accept - The resampled particle satisfies the new metric condition
			for (size_t j = 0; j < 2*model->nout + model->nparams; j++)
			{
				// Keep the particle
				gsl_matrix_set (x, i, j, thetac[j]);
			}

			/*
			// DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
			// Visualize the particle
			for (size_t k = 0; k < 2*model->nout + model->nparams + model->nin; k++)
			{
				printf("%f,\n", thetac[k]);
			}
			printf("\n");
			*/


			// Calculate particle weight
			if (gsl_matrix_isnull (prev))
			{
				gsl_vector_set (wt, i, 1);
			}
			else
			{
				/*
				// Independent univariate prior & multivariate perturbation
				double lnum = 0.0;
				for (size_t k = 0; k < 2*model->nout + model->nparams; k++)
				{
					lnum += log (dunif (thetac[k], lb[k], ub[k]));
				}
				double den = 0.0;
				gsl_vector_view thetacv = gsl_vector_view_array (thetac, 2*model->nout + model->nparams);
				for (size_t l = 0; l < nsam; l++)
				{
					gsl_vector_view thetap = gsl_matrix_row (prev, l);
					den += gsl_vector_get (weights, l) * dmvnorm (&thetacv.vector, &thetap.vector, cov);
				}
				gsl_vector_set (wt, i, exp (lnum)/den);
				*/


				// Independent univariate prior & independent univariate perturbation
				double lnum = 0.0;
				for (size_t k = 0; k < 2*model->nout + model->nparams; k++)
				{
					if (lb[k] == 0 && ub[k]== 0)
						continue;
					lnum += log (dunif (thetac[k], lb[k], ub[k]));
				}
				double den = 0.0;
				for (size_t l = 0; l < nsam; l++)
				{
					double lden = 0.0;
					for (size_t k = 0; k < 2*model->nout + model->nparams; k++)
					{
						if (lb[k] == 0 && ub[k]== 0)
							continue;
						lden += log (dperturb (thetac[k], gsl_matrix_get (prev, l, k), 2*epsilon));
						//lden += log (dunif (gsl_matrix_get (prev, l, k), thetac[k]/1.1, thetac[k]*1.1));
					}
					den += gsl_vector_get (weights, l) * exp (lden);
				}
				gsl_vector_set (wt, i, exp (lnum)/den);


				// printf ("debug info: lnum = %f , den = %f, weight = %e\n", lnum, den, exp(lnum)/den);
			}

			// Move on
			i++;
			// Diagnostics
			if (i % print_every == 0)
				printf ("i = %d, S = %d, epsilon = %f, particle accepted, weight = %e\n", (int) i, (int) S, epsilon, gsl_vector_get (wt, i-1));
		}
		// Advance the ct counter
		ct++;
	}

	// The SMC loop is done

	// Tell all the slaves to exit by sending an empty message with the DIETAG
	for (rank = 1; rank < ntasks; rank++)
	{
		MPI_Send (0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
	}

	// Save the results

	// Normalize the weights
	double sow = gsl_blas_dasum (wt);
	gsl_vector_scale (wt, 1/sow);

	// Write particles to file
	FILE * datafile = fopen (argv[3], "w");
	if (datafile == NULL)
	{
		printf (">> error opening file for writing.. program exiting\n");
		return;
	}
	if (gsl_matrix_fprintf (datafile, x, "%f")) {
		printf (">> error writing simulation data.. program exiting\n");
		return;
	}
	fclose (datafile);

	// Write weights to file
	char wtfilename[500];
	sprintf (wtfilename, "%s_weights", argv[3]);
	FILE * wtfile = fopen (wtfilename, "w");
	if (wtfile == NULL)
	{
		printf (">> error opening weight file for writing.. program exiting\n");
		return;
	}
	if (gsl_vector_fprintf (wtfile, wt, "%e")) {
		printf (">> error writing weights data.. program exiting\n");
		return;
	}
	fclose (wtfile);

	printf("\nWritten data to file:\n%s\n", argv[3]);
	printf("The distance tolerance was: %f\n", epsilon);
	printf ("The particle acceptance rate was: %f\n", ((double) i)/((double) ct));

	// Free manually allocated resources
	rxn_ensemble_free (res);
	gsl_matrix_free (x);
	gsl_vector_free (wt);
	gsl_matrix_free (counts);
	gsl_rng_free (r);

	return;
}


// Slave process tasks
static void slave ()
{
	// Declare some variables
	int myrank;
	size_t res_no = (model->nspecies)*(times->size);
	double pars[model->nparams + model->nin], results[res_no];
	MPI_Status status;

	// Find my rank
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

	// Set up the random number generator
	gsl_rng_env_setup ();
	gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (r, (unsigned long int) time (NULL) + myrank);

	// Allocate necessary objects
	gsl_vector * X0 = gsl_vector_alloc (model->nspecies);
	rxn_sample_path * rsp = rxn_sample_path_alloc (model->nspecies, times->size);
	gsl_vector_view params;

	while (1)
	{
		// Receive a message from the master
		MPI_Recv (pars, model->nparams + model->nin, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		// Check the tag of the received message
		if (status.MPI_TAG == DIETAG)
		{
			printf ("This is process %d, I am exiting.\n", myrank);
			break;
		}

		// Create the vector view with the received parameters
		params = gsl_vector_view_array (pars, model->nparams + model->nin);

		// Sample a new initial state
		model->initial (X0, r);

		// Do the work
		ssa_direct_trajectory (model, &params.vector, times, X0, rsp, r);

		// Communicate that the simulation has been done
		//printf ("This is slave process %d, completed simulation %d.\n", myrank, status.MPI_TAG);

		// Send the result back
		MPI_Send (rsp->counts->data, res_no, MPI_DOUBLE, 0, status.MPI_TAG - WORKTAG, MPI_COMM_WORLD);
	}

	// Free manually allocated resources
	gsl_rng_free (r);
	gsl_vector_free (X0);
	rxn_sample_path_free (rsp);

	return;
}


// Evaluate the uniform density
double dunif (double x, double a, double b)
{
	return (x >= a && x <=b) ? 1/(b-a) : 0;
}


// Multivariate normal PDF
double dmvnorm (const gsl_vector * x, const gsl_vector * mean, const gsl_matrix * var)
{
	// Detect problem size
	size_t n = x->size;

	// Allocate necessary variables
	int s;
	double ax,ay;
	gsl_matrix * work = gsl_matrix_alloc (n, n);
	gsl_matrix * winv = gsl_matrix_alloc (n, n);
	gsl_permutation * p = gsl_permutation_alloc (n);

	// Compute the determinant of the inverse of the covariance matrix
	gsl_matrix_memcpy (work, var);
	gsl_linalg_LU_decomp (work, p, &s);
	gsl_linalg_LU_invert (work, p, winv);
	ax = gsl_linalg_LU_det (work, s);

	// Evaluate the density
	gsl_vector * xm = gsl_vector_alloc (n);
	gsl_vector_memcpy (xm, x);
	gsl_vector_sub (xm, mean);
	gsl_vector * ym = gsl_vector_alloc (n);
	gsl_blas_dsymv (CblasUpper, 1.0, winv, xm, 0.0, ym);
	gsl_blas_ddot (xm, ym, &ay);
	ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax);

	// Free the manually allocated resources
	gsl_vector_free (xm);
	gsl_vector_free (ym);
	gsl_matrix_free (work);
	gsl_matrix_free (winv);
	gsl_permutation_free (p);

	return ay;
}


// Generate a normal random number with mean m and ***variance*** v
double rnorm (const gsl_rng * r, double m, double v)
{
	return m + gsl_ran_gaussian (r, sqrt (v));
}


// Multivariate normal random numbers
int rmvnorm (const gsl_rng * r, const gsl_vector * mean, const gsl_matrix * var, double result[])
{
	// Detect size of the vector to be generated
	size_t n = mean->size;

	// Allocate a vector to store the results
	gsl_vector * res = gsl_vector_alloc (n);

	// Compute Cholesky factorization of the covariance matrix
	gsl_matrix * work = gsl_matrix_alloc (n,n);
	gsl_matrix_memcpy (work, var);
	gsl_linalg_cholesky_decomp (work);

	// Fill the result vector with standard normal random numbers
	for (size_t k = 0; k < n; k++)
	  gsl_vector_set (res, k, gsl_ran_ugaussian (r));

	// Transform the result vector using the Cholesky factorization
	gsl_blas_dtrmv (CblasLower, CblasNoTrans, CblasNonUnit, work, res);
	gsl_vector_add (res, mean);

	// Transfer results into the double array
	for (size_t k = 0; k < n; k++)
		  result[k] = gsl_vector_get (res, k);

	// Free manually allocated resources
	gsl_matrix_free (work);
	gsl_vector_free (res);

	// Signal that the computation was done correctly
	return GSL_SUCCESS;
}


// Generate a uniform random number between a and b
double runif (const gsl_rng * r, double a, double b)
{
	return a + (b-a) * gsl_rng_uniform (r);
}


// Sample a fluorescence level given the background and the count
double sample_fl (const gsl_rng * r, double z, double mufp, double sigfp, gsl_matrix * bg)
{
	double sambg = 0.0, samfl = 0.0;

	// Sample a background value if a bg is present
	if (bg)
	{
		sambg = gsl_matrix_get (bg, gsl_rng_uniform_int (r, M), 0);
	}

	// Sample a fluorescence level if the count is not zero
	if (z > 0)
	{
		samfl = mufp*z + gsl_ran_gaussian (r, sqrt (sigfp*sigfp*z));
	}

	// Return
	return sambg + samfl;
}


// Generate a randomly perturbed resampled particle
double rperturb (const gsl_rng * r, double theta, double wd)
{
	return runif (r, theta * (1 - wd/2), theta * (1 + wd/2));
}


// Evaluate the perturbation density
double dperturb (double x, double theta, double wd)
{
	return dunif (x, theta * (1 - wd/2), theta * (1 + wd/2));
}


// Multivariate normal perturbation
int mvrperturb (const gsl_rng * r, const double theta[], const gsl_matrix * Sig, double result[])
{
	gsl_vector_const_view mean = gsl_vector_const_view_array (theta, 2*P + R);
	rmvnorm (r, &mean.vector, Sig, result);
	return GSL_SUCCESS;
}


// Weighted covariance matrix
int covmat (const gsl_matrix * dat, const gsl_vector * weights, gsl_matrix * Sig)
{
	// Detect problem dimensions
	size_t nsam = dat->size1;
	size_t nvar = dat->size2;

	// Compute some sums
	double s = 0, ssq = 0;
	for (size_t l = 0; l < nsam; l++)
	{
		// Sum of the weights. This should be 1, but the user might be dumb so...
		s += gsl_vector_get (weights, l);
		// Sum of the squared weights
		ssq += pow (gsl_vector_get (weights, l), 2);
	}
	// Squared sum. This should be 1 as well, but see above!
	double sqs = pow (s, 2);

	// Compute the weighted sample means
	double m[nvar];
	for (size_t i = 0; i < nvar; i++)
	{
		m[i] = 0.0;
		for (size_t l = 0; l < nsam; l++)
		{
			m[i] += 1/s * gsl_vector_get (weights, i) * gsl_matrix_get (dat, l, i);
		}
	}

	// Compute the covariance matrix
	for (size_t i = 0; i < nvar; i++)
	{
		// Compute the diagonal element
		double qii = 0.0;
		for (size_t l = 0; l < nsam; l++)
		{
			qii += s/(sqs-ssq) * gsl_vector_get (weights, l) * (gsl_matrix_get (dat, l, i) - m[i])*(gsl_matrix_get (dat, l, i) - m[i]);
		}
		gsl_matrix_set (Sig, i, i, qii);

		for (size_t j = i + 1; j < nvar; j++)
		{
			// Compute the off-diagonal elements
			double qij = 0.0;
			for (size_t l = 0; l < nsam; l++)
			{
				qij += s/(sqs-ssq) * gsl_vector_get (weights, l) * (gsl_matrix_get (dat, l, i) - m[i])*(gsl_matrix_get (dat, l, j) - m[j]);
			}
			gsl_matrix_set (Sig, i, j, qij);
			gsl_matrix_set (Sig, j, i, qij);
		}
	}

	return GSL_SUCCESS;
}


/**
 This function processes the <priors> node of the XML problem definition file. The lower bound and upper
 bound information for each parameter are written in the lb and ub arrays. The arrays need to be preallocated
 to the correct dimension.
 */
int processPriorData (xmlDocPtr doc, double lb[], double ub[])
{
	// Initialize a new XPath evaluation context
	xmlXPathContextPtr ctxt = xmlXPathNewContext (doc);

	// Initialize a xmlXPathObjectPtr object
	xmlXPathObjectPtr result;

	// Apply XPath to find the set of <prior> nodes
	result = xmlXPathEvalExpression ((xmlChar *) "/problem/priors/prior", ctxt);
	// Find out how many nodes were found (if any)
	size_t size = (result->nodesetval) ? result->nodesetval->nodeNr : 0;
	// Check that we have a non-empty nodeset
	if (size != R)
	{
		fprintf (stderr, "error in processPriorData: the number of priors is incorrect\n");
		return GSL_FAILURE;
	}

	// Iterate through the resulting node set
	for (size_t i = 0; i < size; i++)
	{
		// Get first element child of current node -- this is the node <index>
		xmlNode * cur = xmlFirstElementChild (result->nodesetval->nodeTab[i]);

		// Extract the index of current parameter
		size_t j = (size_t) xmlXPathCastNodeToNumber (cur);

		// Advance to the next element sibling of the current node
		cur = xmlNextElementSibling (cur);

		// This is the node <dist>, which is currently not used.

		// Advance to the next element sibling of the current node
		cur = xmlNextElementSibling (cur);

		// Extract the lower bound encoded in the node <lb>
		double lbc = xmlXPathCastNodeToNumber (cur);

		// Advance to the next element sibling of the current node
		cur = xmlNextElementSibling (cur);

		// Extract the upper bound encoded in the node <ub>
		double ubc = xmlXPathCastNodeToNumber (cur);

		// Put stuff in place
		lb[2*P+j-1] = lbc;
		ub[2*P+j-1] = ubc;
	}

	// Free XML resources
	xmlXPathFreeObject (result);


	// Apply XPath to find the <mean> node for the scaling factors
	result = xmlXPathEvalExpression ((xmlChar *) "/problem/scaling/mean", ctxt);
	// Find out how many nodes were found (if any)
	size = (result->nodesetval) ? result->nodesetval->nodeNr : 0;
	// Check that we have a non-empty nodeset
	if (size != P)
	{
		fprintf (stderr, "error in processPriorData: the number of scaling means is incorrect\n");
		return GSL_FAILURE;
	}

	// Iterate through the resulting node set
	for (size_t i = 0; i < size; i+=2)
	{
		// Get first element child of current node -- this is the node <index>
		xmlNode * cur = xmlFirstElementChild (result->nodesetval->nodeTab[i]);

		// Extract the index of current parameter
		size_t j = (size_t) xmlXPathCastNodeToNumber (cur);

		// Advance to the next element sibling of the current node
		cur = xmlNextElementSibling (cur);

		// This is the node <dist>, which is currently not used.

		// Advance to the next element sibling of the current node
		cur = xmlNextElementSibling (cur);

		// Extract the lower bound encoded in the node <lb>
		double lbc = xmlXPathCastNodeToNumber (cur);

		// Advance to the next element sibling of the current node
		cur = xmlNextElementSibling (cur);

		// Extract the upper bound encoded in the node <ub>
		double ubc = xmlXPathCastNodeToNumber (cur);

		// Put stuff in place
		lb[i] = lbc;
		ub[i] = ubc;
	}

	// Free XML resources
	xmlXPathFreeObject (result);


	// Apply XPath to find the <mean> node for the scaling factors
	result = xmlXPathEvalExpression ((xmlChar *) "/problem/scaling/std", ctxt);
	// Find out how many nodes were found (if any)
	size = (result->nodesetval) ? result->nodesetval->nodeNr : 0;
	// Check that we have a non-empty nodeset
	if (size != P)
	{
		fprintf (stderr, "error in processPriorData: the number of scaling standard deviations is incorrect\n");
		return GSL_FAILURE;
	}

	// Iterate through the resulting node set
	for (size_t i = 0; i < size; i+=2)
	{
		// Get first element child of current node -- this is the node <index>
		xmlNode * cur = xmlFirstElementChild (result->nodesetval->nodeTab[i]);

		// Extract the index of current parameter
		size_t j = (size_t) xmlXPathCastNodeToNumber (cur);

		// Advance to the next element sibling of the current node
		cur = xmlNextElementSibling (cur);

		// This is the node <dist>, which is currently not used.

		// Advance to the next element sibling of the current node
		cur = xmlNextElementSibling (cur);

		// Extract the lower bound encoded in the node <lb>
		double lbc = xmlXPathCastNodeToNumber (cur);

		// Advance to the next element sibling of the current node
		cur = xmlNextElementSibling (cur);

		// Extract the upper bound encoded in the node <ub>
		double ubc = xmlXPathCastNodeToNumber (cur);

		// Put stuff in place
		lb[i+1] = lbc;
		ub[i+1] = ubc;
	}

	// Free XML resources
	xmlXPathFreeObject (result);
	xmlXPathFreeContext (ctxt);

	return GSL_SUCCESS;
}


/**
 This function processes the <data> node of the XML problem definition file and reads
 the specified data set.
 */
int processDataSet (xmlDocPtr doc, gsl_vector * times, gsl_matrix * data)
{
	// Initialize a new XPath evaluation context
	xmlXPathContextPtr ctxt = xmlXPathNewContext (doc);

	// Initialize a xmlXPathObjectPtr object
	xmlXPathObjectPtr result;

	// Apply XPath to find the <data> node and read its attributes to determine the source of the data and the scale
	result = xmlXPathEvalExpression ((xmlChar *) "/problem/data", ctxt);
	xmlChar * mode = xmlGetProp (result->nodesetval->nodeTab[0], (xmlChar *) "mode");
	xmlChar * scale = xmlGetProp (result->nodesetval->nodeTab[0], (xmlChar *) "scale");

	// Set the data mode
	if (strstr (mode, "fcs"))
	{
		DATA_MODE = DATA_MODE_FCS;
	}
	else if (strstr (mode, "text"))
	{
		DATA_MODE = DATA_MODE_TEXT;
	}
	else
	{
		fprintf (stderr, "error in processDataSet: the data mode is not valid\n");
		return GSL_FAILURE;
	}

	// Set the data scale
	if (strstr (scale, "log"))
	{
		DATA_SCALE = DATA_SCALE_LOG;
	}
	else if (strstr (scale, "linear"))
	{
		DATA_SCALE = DATA_SCALE_LINEAR;
	}
	else
	{
		fprintf (stderr, "error in processDataSet: the data scale is not valid\n");
		return GSL_FAILURE;
	}

	// Free the objects used so far and move to the data loading section
	xmlXPathFreeObject (result);
	xmlFree (mode);
	xmlFree (scale);

	// Determine where the data should be read from
	switch (DATA_MODE)
	{
	case DATA_MODE_FCS:
	{
		// The data is contained in a bunch of FCS files

		// Apply XPath to find the set of <timept> nodes
		result = xmlXPathEvalExpression ((xmlChar *) "/problem/data/timept", ctxt);

		// Find out how many nodes were found (if any)
		size_t size = (result->nodesetval) ? result->nodesetval->nodeNr : 0;

		// Check that we have a non-empty nodeset
		if (size != K)
		{
			fprintf (stderr, "error in processDataSet: size mismatch in time points\n");
			return GSL_FAILURE;
		}

		// Iterate through the resulting node set
		for (size_t i = 0; i < size; i++)
		{
			// Get first element child of current node -- this is the node <index>
			xmlNode * cur = xmlFirstElementChild (result->nodesetval->nodeTab[i]);
			// Extract the index of current time point
			size_t j = (size_t) xmlXPathCastNodeToNumber (cur);

			// Advance to the next element sibling of the current node -- this is <time>
			cur = xmlNextElementSibling (cur);
			// Extract the time of the measurement
			double t = xmlXPathCastNodeToNumber (cur);
			// Put the value in the appropriate place of the vector times
			gsl_vector_set (times, j-1, t);

			// Advance to the next element sibling of the current node -- this is <unit> -- currently unused
			cur = xmlNextElementSibling (cur);

			// Advance to the next element sibling of the current node -- this is <fcsfile>
			cur = xmlNextElementSibling (cur);
			// Extract the path to the FCS file containing the data
			xmlChar * path = xmlXPathCastNodeToString (cur);
			// Open the FCS file
			FILE * fcs = fopen (path, "rb");
			if (fcs == NULL)
			{
				fprintf (stderr, "error in processDataSet: could not open FCS file %d\n", (int) j);
				return GSL_FAILURE;
			}

			// Read the FCS file header
			float ver;
			long int dos;
			fscanf (fcs, "FCS%f    %*d%*d%ld%*d%*d%*d", &ver, &dos);

			// printf ("%f, %ld\n", ver, dos);


			// Initialize a new xmlXPathObjectPtr object
			xmlXPathObjectPtr result2;
			// Create the XPath expression to retrieve the node set of <measurement> for the current time point
			char xpath[60];
			sprintf (xpath, "/problem/data/timept[%d]/measuredspecies/measurement", (int) i+1);
			// Apply the XPath expression
			result2 = xmlXPathEvalExpression ((xmlChar *) xpath, ctxt);
			// Find out how many nodes were found (if any)
			size_t size2 = (result2->nodesetval) ? result2->nodesetval->nodeNr : 0;
			// Check that we have a non-empty nodeset
			if (size2 != P)
			{
				fprintf (stderr, "error in processDataSet: size mismatch in measured species\n");
				return GSL_FAILURE;
			}

			// Initialize arrays to contain indices and FCS parameter of the measured species
			size_t l[size2], p[size2];

			// Iterate through the resulting node set
			for (size_t k = 0; k < size2 ; k++)
			{
				// Get first element child of current node -- this is the node <index>
				xmlNode * cur2 = xmlFirstElementChild (result2->nodesetval->nodeTab[k]);
				// Extract the index of current measured species
				l[k] = (size_t) xmlXPathCastNodeToNumber (cur2);

				// Advance to the next element sibling of the current node -- this is <description> -- currently unused
				cur2 = xmlNextElementSibling (cur2);

				// Advance to the next element sibling of the current node -- this is <fcsparam>
				cur2 = xmlNextElementSibling (cur2);
				// Extract the parameter of current measured species
				p[k] = (size_t) xmlXPathCastNodeToNumber (cur2);
			}

			// Find out how many total events are stored in the FCS file
			int tot = 0;
			if (ver < 3.0)
				tot = fcs2_read_int_kw (fcs, "$TOT");
			else
				tot = fcs3_read_int_kw (fcs, "$TOT");
			// Check that there are enough events
			if (tot < M)
			{
				fprintf (stderr, "error in processDataSet: the FCS file %d does not contain enough events\n", (int) i+1);
				return GSL_FAILURE;
			}

			//printf ("$TOT = %d\n", tot);

			// Find out how many parameters are there in each event
			int par = 0;
			if (ver < 3.0)
				par = fcs2_read_int_kw (fcs, "$PAR");
			else
				par = fcs3_read_int_kw (fcs, "$PAR");

			// Initialize read buffer
			float fbuf[par];

			//printf ("$PAR = %d\n", par);

			// Iterate through the parameters that have to be read
			for (size_t k = 0; k < size2; k++)
			{
				// Position the file cursor at the beginning of the data segment
				fseek (fcs, dos, SEEK_SET);

				// Read the first M events greater than 1.1 of each parameter
				for (size_t m = 0; m < M; m++)
				{
					// Get a chunk of floats
					do
					{
						fread_floats_swap (fbuf, par, fcs);
						//printf("read number: %f\n", *fbuf);
					}
					while (fbuf[p[k]-1] <= 1.1);

					// Place the data in the proper spot
					//printf ("number read\n");
					gsl_matrix_set (data, m, P*(j-1)+l[k]-1, fbuf[p[k]-1]);
				}
			}

			// Close the FCS file
			fclose (fcs);

			// Free XML resources
			xmlXPathFreeObject (result2);
			xmlFree (path);
		}

		// Free more XML resources
		xmlXPathFreeObject (result);

		break;
	}

	case DATA_MODE_TEXT:
	{
		// The data is contained in a single text file

		// Apply XPath to find the set of <timept> nodes
		result = xmlXPathEvalExpression ((xmlChar *) "/problem/data/timept", ctxt);

		// Find out how many nodes were found (if any)
		size_t size = (result->nodesetval) ? result->nodesetval->nodeNr : 0;

		// Check that we have a non-empty nodeset
		if (size != K)
		{
			fprintf (stderr, "error in processDataSet: size mismatch in time points\n");
			return GSL_FAILURE;
		}

		// Iterate through the resulting node set
		for (size_t i = 0; i < size; i++)
		{
			// Get first element child of current node -- this is the node <index>
			xmlNode * cur = xmlFirstElementChild (result->nodesetval->nodeTab[i]);
			// Extract the index of current time point
			size_t j = (size_t) xmlXPathCastNodeToNumber (cur);

			// Advance to the next element sibling of the current node -- this is <time>
			cur = xmlNextElementSibling (cur);
			// Extract the time of the measurement
			double t = xmlXPathCastNodeToNumber (cur);
			// Put the value in the appropriate place of the vector times
			gsl_vector_set (times, j-1, t);
		}

		// Free the result object before using it again
		xmlXPathFreeObject (result);

		// Apply XPath to find the path of the data file, contained in the node <datafile>
		result = xmlXPathEvalExpression ((xmlChar *) "/problem/data/datafile", ctxt);
		// Get the string value of the node
		xmlChar * path = xmlXPathCastNodeToString (result->nodesetval->nodeTab[0]);

		// Override FCS file processing to load text data
		FILE * synfile = fopen (path, "r");
		if (gsl_matrix_fscanf (synfile, data) != GSL_SUCCESS)
		{
			fprintf (stderr, "error: could not load text-mode data\n");
				return EXIT_FAILURE;
		}
		fclose (synfile);

		// Free XML resources
		xmlXPathFreeObject (result);
		xmlFree (path);

		break;
	}
	}

	// Free XML resources
	xmlXPathFreeContext (ctxt);

	// Determine whether the data should be transformed to log-scale
	switch (DATA_SCALE)
	{
	case (DATA_SCALE_LOG):
	{
		// The data is to be converted to log scale

		for (size_t i = 0; i < data->size1; i++)
		{
			for (size_t j = 0; j < data->size2; j++)
			{
				gsl_matrix_set (data, i, j, log10 (gsl_matrix_get (data, i, j)));
			}
		}
		break;
	}

	case (DATA_SCALE_LINEAR):
	{
		// The data will be used as supplied in the files
		break;
	}
	}

	// Pre-sort the columns of the data matrix
	gsl_matrix * data2 = gsl_matrix_calloc (data->size1, data->size2);

	// Scan the columns of data2
	for (size_t j = 0; j < data2->size2; j++)
	{
		// Read the j-th column of data2
		gsl_vector_view ccol = gsl_matrix_column (data, j);
		// Sort this vector
		gsl_sort_vector (&ccol.vector);
		// Copy the sorted vector into the matrix data2
		gsl_matrix_set_col (data2, j, &ccol.vector);
	}

	// Overwrite the unsorted matrix with the sorted one
	gsl_matrix_memcpy (data, data2);

	// Free data2
	gsl_matrix_free (data2);

	// Return
	return GSL_SUCCESS;
}


/**
 This function processes the <background> node of the XML problem definition file and reads
 the specified FCS file.
 */
int processBackground (xmlDocPtr doc, gsl_matrix * bg)
{
	// Initialize a new XPath evaluation context
	xmlXPathContextPtr ctxt = xmlXPathNewContext (doc);

	// Initialize a xmlXPathObjectPtr object
	xmlXPathObjectPtr result = NULL;

	// Check if the user specified a <background> node
	result = xmlXPathEval ((xmlChar *) "/problem/data/background/fcsfile", ctxt);

	// If the resulting node set is non-empty (i.e. there is a bg distribution)
	if (result->nodesetval->nodeNr)
	{
		// Extract the path of the background distribution
		xmlChar * path = xmlXPathCastNodeToString (result->nodesetval->nodeTab[0]);

		// Open the FCS file
		FILE * fcs = fopen (path, "rb");
		if (fcs == NULL)
		{
			fprintf (stderr, "error in processDataSet: could not open background FCS file\n");
			return GSL_FAILURE;
		}
		// Read the FCS file header
		float ver;
		long int dos;
		fscanf (fcs, "FCS%f    %*d%*d%ld%*d%*d%*d", &ver, &dos);

		// Initialize a new xmlXPathObjectPtr object
		xmlXPathObjectPtr result2;
		// Apply XPath expression to find the background measurement specification
		result2 = xmlXPathEvalExpression ((xmlChar *) "/problem/data/background/measuredspecies/measurement", ctxt);
		// Find out how many nodes were found (if any)
		size_t size2 = (result2->nodesetval) ? result2->nodesetval->nodeNr : 0;
		// Check that we have a non-empty nodeset
		if (size2 != P)
		{
			fprintf (stderr, "error in processDataSet: size mismatch in background data\n");
			return GSL_FAILURE;
		}
		// Find out how many total events are stored in the FCS file
		int tot = 0;
		if (ver < 3.0)
			tot = fcs2_read_int_kw (fcs, "$TOT");
		else
			tot = fcs3_read_int_kw (fcs, "$TOT");
		// Check that there are enough events
		if (tot < M)
		{
			fprintf (stderr, "error in processDataSet: the background FCS file does not contain enough events\n");
			return GSL_FAILURE;
		}

		// Find out how many parameters are there in each event
		int par = 0;
		if (ver < 3.0)
			par = fcs2_read_int_kw (fcs, "$PAR");
		else
			par = fcs3_read_int_kw (fcs, "$PAR");

		// Initialize read buffer
		float fbuf[par];

		// Iterate through the resulting node set and read the events from the FCS file (keeping only the ones > 1.0)
		for (size_t k = 0; k < size2 ; k++)
		{
			// Get first element child of current node -- this is the node <index>
			xmlNode * cur2 = xmlFirstElementChild (result2->nodesetval->nodeTab[k]);
			// Extract the index of current measured species
			size_t l = (size_t) xmlXPathCastNodeToNumber (cur2);

			// Advance to the next element sibling of the current node -- this is <description> -- currently unused
			cur2 = xmlNextElementSibling (cur2);

			// Advance to the next element sibling of the current node -- this is <fcsparam>
			cur2 = xmlNextElementSibling (cur2);
			// Extract the parameter of current measured species
			size_t p = (size_t) xmlXPathCastNodeToNumber (cur2);

			// Position the file cursor at the beginning of the data segment
			fseek (fcs, dos, SEEK_SET);

			// Read M events
			for (size_t m = 0; m < M; m++)
			{
				// Keep reading until the desired parameter is > 1.1
				do
				{
					fread_floats_swap (fbuf, par, fcs);
				}
				while (fbuf[p-1]<=1.1);

				// Place the data in the proper spot
				gsl_matrix_set (bg, m, l-1, fbuf[p-1]);
			}
		}

		// Close the FCS file
		fclose (fcs);
	}
	else
	{
		gsl_matrix_set_zero (bg);
	}

	// Return
	return GSL_SUCCESS;
}


/**
 This function processes the <model> node of the XML problem definition file and looks
 for the specified model in the library StochMod.
 */
int processModel (xmlDocPtr doc, stochmod * model)
{
	// Initialize return value
	int retval = GSL_SUCCESS;

	// Initialize a new XPath evaluation context
	xmlXPathContextPtr ctxt = xmlXPathNewContext (doc);

	// Initialize a xmlXPathObjectPtr object
	xmlXPathObjectPtr result = NULL;

	// Apply XPath to retrieve the <model> node
	result = xmlXPathEval ((xmlChar *) "/problem/model", ctxt);

	// If the resulting node set is non-empty (i.e. model has been specified)
	if (result->nodesetval->nodeNr)
	{
		// Retrieve the model number
		STOCHASTIC_MODEL modnr = (STOCHASTIC_MODEL) xmlXPathCastNodeToNumber (result->nodesetval->nodeTab[0]);

		// Match the model number to the correct model from the library StochMod
		switch (modnr)
		{
		case (MODEL_SYNCIRC):
		{
			syncirc_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_STOCHREP):
		{
			stochrep_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_AUTOREG):
		{
			autoreg_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_LACGFP):
		{
			lacgfp_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_LACGFP2):
		{
			lacgfp2_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_LACGFP3):
		{
			lacgfp3_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_LACGFP4):
		{
			lacgfp4_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_LACGFP5):
		{
			lacgfp5_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_BIRTHDEATH):
		{
			birthdeath_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_LACGFP6):
		{
			lacgfp6_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_LACGFP7):
		{
			lacgfp7_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_LACGFP8):
		{
			lacgfp8_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_IFF):
		{
			iff_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_FBK):
		{
			fbk_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_LACGFP9):
		{
			lacgfp9_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_LACGFP10):
		{
			lacgfp10_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		case (MODEL_SYNPI1):
		{
			synpi1_mod_setup (model);
			retval = GSL_SUCCESS;
			break;
		}
		default:
		{
			fprintf (stderr, "error in processModel: no match for model code found in library StochMod\n");
			retval = GSL_FAILURE;
			break;
		}
		}
	}
	else
	{
		fprintf (stderr, "error in processModel: <model> node not found in the problem definition file\n");
		retval = GSL_FAILURE;
	}

	// Free XML resources
	xmlXPathFreeObject (result);
	xmlXPathFreeContext (ctxt);

	// Return
	return retval;
}
