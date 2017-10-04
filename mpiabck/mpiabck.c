/*
 *  mpiabck.c
 *  mpiStochEst
 *
 *  Parameter inference for discretely observed stochastic chemical reaction networks
 *  using Approximate Bayesian Computation with the DCMS metric
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
stochmod * model;			// The model to estimate
gsl_matrix * data;			// The flow cytometry data
gsl_matrix * out; 			// The output matrix of the model
gsl_vector * times;			// The time points
gsl_matrix * bg;			// The background distribution
double * lb, * ub;			// Box constraints on the parameters
size_t M;					// Number of experimental samples to use for identification
size_t P;					// Number of measured species
size_t K;					// Number of time points
size_t R;					// Number of parameters
double epsilon;				// Distance tolerance
double mufp = 220.0;		// Mean fluorescence level
double sigfp = 390.0;		// Sd of fluorescence level

// Tags for message passing
#define WORKTAG 1
#define DIETAG 0

// Functions for parallel simulation
static void master (int argc, char * argv[]);
static void slave ();

// Utility functions
static double unif (double x, double a, double b);
static double ran_gaussian (const gsl_rng * r, double m, double v);
static double ran_uniform (const gsl_rng * r, double a, double b);

// Functions for problem definition processing
static int processPriorData (xmlDocPtr doc, double lb[], double ub[]);
static int processDataSet (xmlDocPtr doc, gsl_vector * times, gsl_matrix * data);
static int processBackground (xmlDocPtr doc, gsl_matrix * bg);


int main (int argc, char * argv[])
{
	// ================================== USER-PROVIDED INFORMATION ==========================================
	// Choose the stochasitc model to estimate
	model = (stochmod *) malloc (sizeof (stochmod));
	if (model == NULL)
	{
		fprintf (stderr, "error in memory allocation\n");
		return EXIT_FAILURE;
	}
	lacgfp2_mod_setup (model);

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
	if (argc != 4)
	{
		if (myrank == 0)
			fprintf (stderr, "usage: %s problem-file-path results-file-path nsim\n\n", argv[0]);
		return EXIT_FAILURE;
	}

	// Put command line arguments into correct variables
	char * problemfile = argv[1];
	char * resultsfile = argv[2];
	size_t nsim = atoi (argv[3]);

	// Inform the user that we are ready to start
	if (myrank == 0)
	{
		printf ("\nThis is mpiAbcK v1.0\nI am the master process\n\n");
		printf ("Problem file path:\n%s\n", problemfile);
		printf ("Results will be saved in:\n%s\n", resultsfile);
		printf ("Will estimate %s using %d MCMC samples.\n\n", model->name, (int) nsim);
	}

	// =============================== EXTRACT INFO FROM PROBLEM FILE ========================================
	// Parse the problem file
	xmlDocPtr problem = xmlReadFile (problemfile, NULL, 0);
	if (problem == NULL)
	{
		fprintf (stderr, "error reading problem file\n");
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
	lb = (double *) malloc (R * sizeof(double));
	ub = (double *) malloc (R * sizeof(double));
	if (lb == NULL || ub == NULL)
	{
		fprintf (stderr, "error in memory allocation\n");
		return EXIT_FAILURE;
	}
	data = gsl_matrix_alloc (M, K*P);
	times = gsl_vector_alloc (K);

	// Process the prior information
	if (processPriorData (problem, lb, ub) != GSL_SUCCESS)
	{
		fprintf (stderr, "error: the problem definition file does not contain valid prior information\n");
		return EXIT_FAILURE;
	}

	// Load the data from the FCS files
	if (processDataSet (problem, times, data) != GSL_SUCCESS)
	{
		fprintf (stderr, "error: could not process data set information\n");
		return EXIT_FAILURE;
	}

	// Load the background distribution
	bg = gsl_matrix_alloc (M, P);
	if (processBackground (problem, bg) != GSL_SUCCESS)
	{
		fprintf (stderr, "error: could not process background information\n");
		return EXIT_FAILURE;
	}

	// Set up output matrix
	out = gsl_matrix_calloc (model->nout, model->nspecies);
	model->output (out);


	// ============================ RUN THE PARALLEL COMPUTATION ==================================
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
		printf ("\n\nMCMC inference complete.\nElapsed time: %g sec.\n\n", (((double) (toc - tic)) / CLOCKS_PER_SEC));


	// ===================================== CLEAN UP =============================================
	// Free objects
	free (model);
	free (lb);
	free (ub);
	gsl_matrix_free (out);
	gsl_matrix_free (data);
	gsl_vector_free (times);
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
	size_t nsam = atoi (argv[3]);

	// Find out how many processes there are in the default communicator
	MPI_Comm_size (MPI_COMM_WORLD, &ntasks);

	// Set up the random number generator
	gsl_rng_env_setup ();
	gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (r, (unsigned long int) time (NULL));

	// Set metric parameters
	double beta = 0.05; 						// "External" confidence level
	double alpha = 1 - sqrt (1 - beta);			// "Internal" condifence level
	size_t S = 15;								// Number of simulations

	// Set the batch size
	size_t bs = 50;

	// Compute critical epsilon
	epsilon = kolmogorov_cdf_inverse (S, 1 - alpha) + 0.0049250; // kolmogorov_cdf_inverse (90000, 1-alpha);

	// Visualize metric parameters
	printf("\nThe metric parameters are as follows\n");
	printf("Confidence level beta = %f\n", beta);
	printf("Confidence level alpha = %f\n", alpha);
	printf("Number of simulations S = %d\n", (int) S);
	printf("Number of experimental samples M = %d\n", (int) M);
	printf("Distance tolerance epsilon = %f\n\n", epsilon);

	// Initialize parameter matrix
	gsl_matrix * x = gsl_matrix_calloc (nsam, model->nparams);
	// Sample initial guesses from the priors
	for (size_t j = 0; j < x->size2; j++)
	{
		double sam = ran_uniform (r, lb[j], ub[j]);
		gsl_matrix_set (x, 0, j, sam);
	}

	// Initialize distances
	double dc[data->size2];
	int dcg;

	// Initialize acceptance probability
	double A = 0;

	// Initialize parameter, log standard deviation, acceptance counter and batch acceptance counter arrays
	double thetac[model->nparams + model->nin], thetap[model->nparams + model->nin];
	double ls[model->nparams];
	size_t ac[model->nparams], bac[model->nparams];

	// Fill up the arrays with appropriate values
	for (size_t j = 0; j < model->nparams; j++)
	{
		thetac[j] = gsl_matrix_get (x, 0, j);
		thetap[j] = gsl_matrix_get (x, 0, j);

		ls[j] = 0.0;

		ac[j] = 0;
		bac[j] = 0;
	}

	// The extra entries in the thetac and thetap arrays are the model inputs
	thetac[model->nparams] = 0.0;
	thetap[model->nparams] = 0.0;

	// Allocate rxn_ensemble object to store SCRN simulations
	rxn_ensemble * res = rxn_ensemble_alloc (S, model->nspecies, times->size);

	// Allocate matrix to store SCRN output
	gsl_matrix * counts = gsl_matrix_calloc (res->nreplic, model->nout*res->ntimes);

	// Declare necessary views
	gsl_vector_view drow, crow;


	// Main METROPOLIS-HASTINGS (Metropolis-within-Gibbs ABC) loop
	for (size_t i = 1; i < nsam; i++)
	{
		printf ("%d\n", (int) i+1);

		// Check if this is the first iteration of a new batch
		if (i % bs == 0)
		{
			// Compute the batch number
			double bn = ((double) i)/((double) bs);

			// Loop over parameters
			for (size_t j = 0; j < model->nparams; j++)
			{
				// Compute the batch acceptance rate
				double bar = ((double) bac[j])/((double) bs);

				// Perform the adaptation
				if (bar < 0.44)
				{
					// Decrease the log proposal standard dev
					ls[j] -= (1/sqrt(bn) <= 0.01) ? 1/sqrt(bn) : 0.01;
				}
				else
				{
					// Increase the log proposal standard dev
					ls[j] += (1/sqrt(bn) <= 0.01) ? 1/sqrt(bn) : 0.01;
				}

				// Reset the batch acceptance counter
				bac[j] = 0;
			}
		}

		// Loop over parameters
		for (size_t j = 0; j < model->nparams; j++)
		{
			// If zero do not change!
			if (thetap[j] == 0)
			{
				thetac[j] = 0;
				thetap[j] = thetac[j];
				gsl_matrix_set (x, i, j, thetac[j]);
				continue;
			}

			// Propose a new parameter value
			thetac[j] = ran_gaussian (r, thetap[j], exp (ls[j]));

			// Evaluate the prior contribution to the acceptance probability
			A = unif (thetac[j], lb[j], ub[j])/unif (thetap[j], lb[j], ub[j]);

			// Check if the proposed value is feasible
			if (A == 0)
			{
				// Proposed point has no hope
				// Reject - The chain stays at current value
				thetac[j] = thetap[j];
				gsl_matrix_set (x, i, j, thetac[j]);

				// Diagnostics
				printf ("\t par %d - unfeasible proposal\n", (int) j);
				// Move on
				continue;
			}

			// Simulate SCRN with proposed values
			{
				// Request simulations from the slaves by sending the parameters to simulate
				for (rank = 1; rank < ntasks; rank++)
				{
					MPI_Send (thetac, 								// current parameters
							  model->nparams + model->nin,			// # of data items
							  MPI_DOUBLE,							// data items are of type double
							  rank,									// destination process rank
							  WORKTAG + count,						// user chosen message tag
							  MPI_COMM_WORLD);						// default communicator

					count++;
				}

				// Keep requesting simulations untils we have enough
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
					MPI_Send (thetac,								// message buffer
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
						// Resample a background level
						double sambg = gsl_matrix_get (bg, gsl_rng_uniform_int (r, M), 0);

						// Current count
						double z = gsl_matrix_get (counts, k, l);

						// Sample a fluorescence level
						double samfl = mufp*z + gsl_ran_gaussian (r, sqrt (sigfp*sigfp*z));

						gsl_matrix_set (counts, k, l, samfl + sambg);
					}
				}

				dcg = 1;

				// Compute the Kolmogorov distances
				for (size_t k = 1; k < data->size2; k+=2)
				{
					 drow = gsl_matrix_column (data, k);
					 crow = gsl_matrix_column (counts, k);
					 dc[k] = ksdist_two_sample_2 (&crow.vector, &drow.vector);

					 if (dc[k] > epsilon)
						 dcg = 0;
				}
			}

			// Decide what to do
			if (dcg != 0)
			{
				// Accept - The chain moves to the proposed value
				thetap[j] = thetac[j];
				gsl_matrix_set (x, i, j, thetac[j]);
			}
			else
			{
				// Reject - The chain stays at current value
				thetac[j] = thetap[j];
				gsl_matrix_set (x, i, j, thetac[j]);
			}

			// Diagnostics
			printf ("\t par %d - p = %f - c = %f - ls = %f \n", (int) j,
					gsl_matrix_get (x, i-1, j), gsl_matrix_get (x, i, j), ls[j]);
		}
	}

	// The MCMC loop is done

	// Tell all the slaves to exit by sending an empty message with the DIETAG
	for (rank = 1; rank < ntasks; rank++)
	{
		MPI_Send (0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
	}

	// Save the results

	// Write to file
	FILE * datafile = fopen (argv[2], "w");
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

	printf("Written data to file:\n%s\n", argv[2]);

	// Free manually allocated resources
	rxn_ensemble_free (res);
	gsl_matrix_free (x);
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
static double unif (double x, double a, double b)
{
	return (x >= a && x <=b) ? 1/(b-a) : 0;
}


// Generate a normal random number with mean m and ***variance*** v
static double ran_gaussian (const gsl_rng * r, double m, double v)
{
	return m + gsl_ran_gaussian (r, sqrt (v));
}


// Generate a uniform random number between a and b
static double ran_uniform (const gsl_rng * r, double a, double b)
{
	return a + (b-a) * gsl_rng_uniform (r);
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
		fprintf (stderr, "error in processPriorData: size mismatch\n");
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
		lb[j-1] = lbc;
		ub[j-1] = ubc;
	}

	// Free XML resources
	xmlXPathFreeObject (result);
	xmlXPathFreeContext (ctxt);

	return GSL_SUCCESS;
}


/**
 This function processes the <data> node of the XML problem definition file and reads
 the specified FCS files.
 */
int processDataSet (xmlDocPtr doc, gsl_vector * times, gsl_matrix * data)
{
	// Initialize a new XPath evaluation context
	xmlXPathContextPtr ctxt = xmlXPathNewContext (doc);

	// Initialize a xmlXPathObjectPtr object
	xmlXPathObjectPtr result;

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

		// Find out how many parameters are there in each event
		int par = 0;
		if (ver < 3.0)
			par = fcs2_read_int_kw (fcs, "$PAR");
		else
			par = fcs3_read_int_kw (fcs, "$PAR");

		// Position the file cursor at the beginning of the data segment
		fseek (fcs, dos, SEEK_SET);
		// Initialize read buffer
		float fbuf[par];
		// Read the first M events
		for (size_t m = 0; m < M; m++)
		{
			// Get a chunk of floats
			fread_floats_swap (fbuf, par, fcs);
			// Iterate through the parameters that have to be read
			for (size_t k = 0; k < size2; k++)
			{
				// Place the data in the proper spot
				gsl_matrix_set (data, m, P*(j-1)+l[k]-1, fbuf[p[k]-1]);
			}
		}

		// Close the FCS file
		fclose (fcs);

		// Free XML resources
		xmlXPathFreeObject (result2);
		xmlFree (path);
	}

	// Free XML resources
	xmlXPathFreeObject (result);
	xmlXPathFreeContext (ctxt);

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
	xmlXPathObjectPtr result;

	// Check if the user specified a <background> node
	result = xmlXPathEvalExpression ((xmlChar *) "/problem/data/background/fcsfile", ctxt);

	// If the resulting node set is non-empty (i.e. there is a bg distribution)
	if (result->nodesetval)
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
