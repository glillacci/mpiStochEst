/*
 *  mpifind0.c
 *  mpiStochEst
 *
 *  Initial condition finder for
 *  Approximate Bayesian Computation with the DCMS metric
 *
 *  Created by Gabriele Lillacci in February 2012.
 *	Latest revision: May 2012.
 *
 *
 *	This free software is available under the Creative Commons Attribution Share Alike License.
 *	You are permitted to use, redistribute and adapt this software as long as appropriate credit
 *	is given to the original author, and all derivative works are distributed under the same
 *	license or a compatible one.
 *	For more information, visit http://creativecommons.org/licenses/by-sa/3.0/ or send a letter to
 *	Creative Commons, 171 2nd Street, Suite 300, San Francisco, California, 94105, USA.
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

// Custom includes
#include <parestlib.h>
#include <stochmod.h>
#include <auxfun.h>

// Cvxgen includes
#include <solver.h>


// Global variable declarations
stochmod * model;			// The model to estimate
gsl_matrix * data;			// The flow cytometry data
gsl_matrix * out; 			// The output matrix of the model
gsl_vector * times;			// The time points
gsl_vector * init;			// The initial guesses
gsl_matrix * bg;			// The first time point
double * lb, * ub;			// Box constraints on the parameters
size_t M;					// Number of experimental samples to use for identification
size_t P;					// Number of measured species
size_t K;					// Number of time points
size_t R;					// Number of parameters
size_t S; 					// Number of simulations
double epsilon;				// Distance tolerance
double mufp = 220.0;		// Mean fluorescence level
double sigfp = 390.0;		// Sd of fluorescence level

// Tags for message passing
#define WORKTAG 1
#define DIETAG 0

// Functions for parallel simulation
static void master (int argc, char * argv[]);
static void slave ();
static double parallel_sim (double theta[], const gsl_rng * r);

// Functions for problem definition processing
static int processPriorData (xmlDocPtr doc, double lb[], double ub[]);
static int processDataSet (xmlDocPtr doc, gsl_vector * times, gsl_matrix * data);
static int processBackground (xmlDocPtr doc, gsl_matrix * bg);

// cvxgen data structures (global)
Vars vars;
Params params;
Workspace work;
Settings settings;


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

	// Initialize cvxgen
	set_defaults ();
	setup_indexing ();


	// =============================== PROCESS COMMAND LINE ARGUMENTS ========================================

	// Find out my identity in the default communicator
	int myrank;
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

	// Check correctness of command line arguments
	if (argc != 2)
	{
		if (myrank == 0)
			printf ("Usage: %s flowdata\n\n", argv[0]);
		return EXIT_FAILURE;
	}

	// Put command line arguments into correct variables
	char * problemfile = argv[1];

	// Inform the user that we are ready to start
	if (myrank == 0)
	{
		printf ("\nThis is mpiFind0 v1.0\nI am the master process\n\n");
		printf ("Problem file path:\n%s\n", problemfile);
		printf ("I will try to find an initial condition for: %s.\n\n", model->name);
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
		printf ("\n\nInitial condition FOUND.\nElapsed time: %g sec.\n\n", (((double) (toc - tic)) / CLOCKS_PER_SEC));


	// ===================================== CLEAN UP =============================================
	// Free objects
	free (model);
	free (lb);
	free (ub);
	gsl_matrix_free (out);
	gsl_matrix_free (data);
	gsl_vector_free (times);
	gsl_vector_free (init);
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
	int ntasks, rank;

	// Find out how many processes there are in the default communicator
	MPI_Comm_size (MPI_COMM_WORLD, &ntasks);

	// Set up the random number generator
	gsl_rng_env_setup ();
	gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (r, (unsigned long int) time (NULL));

	// Set metric parameters
	double beta = 0.05; 						// "External" confidence level
	double alpha = 1 - sqrt (1 - beta);			// "Internal" condifence level
	S = 1*47;									// Number of simulations

	// Compute critical epsilon
	epsilon = kolmogorov_cdf_inverse (S, 1 - alpha) + 0.0049250; // kolmogorov_cdf_inverse (90000, 1-alpha);

	// Visualize metric parameters
	printf("\nThe metric parameters are as follows\n");
	printf("Confidence level beta = %f\n", beta);
	printf("Confidence level alpha = %f\n", alpha);
	printf("Number of simulations S = %d\n", (int) S);
	printf("Number of experimental samples M = %d\n", (int) M);
	printf("Distance tolerance epsilon = %f\n\n", epsilon);

	// Set initial trust region bound
	double D = 1;

	// Initialize iteration counter
	size_t iter = 1;

	// Calculate number of interpolation points needed for quadratic model
	size_t N = 0.5*(R+1)*(R+2);

	// Setup necessary matrices and vectors
	gsl_matrix * Pt = gsl_matrix_calloc (R, N); 	// Interpolation points
	gsl_matrix * AA = gsl_matrix_calloc (N, N);		// Model matrix
	gsl_vector * BB = gsl_vector_alloc (N);			// Function values
	gsl_matrix * H = gsl_matrix_alloc (R, R);		// Hessian
	gsl_matrix * Hpsd = gsl_matrix_alloc (R, R);	// Positive semidefinite part of Hessian
	gsl_vector * g = gsl_vector_alloc (R);			// Gradient
	gsl_vector * gp = gsl_vector_alloc (R);			// Gradient of convex part
	double c = 0.0;									// The constant term
	gsl_matrix * I = gsl_matrix_alloc (R, R);		// Identity
	gsl_matrix_set_identity (I);

	// Initialize vectors to hold current and previous values of optimization variables and generate a random initial condition
	double xc[R], xp[R];
	for (size_t i = 0; i < R; i++)
	{
		xc[i] = (lb[i] + ub[i])/2;
		xp[i] = xc[i];
	}

	// Initalize other useful variable
	double fxc = K*P;

	// Main loop
	do
	{
		printf ("\n\nIteration %d\n", (int) iter);

		// Randomly sample the interpolation points
		while (interpolation_space (Pt, xc, lb, ub, D, r) != GSL_SUCCESS)
		{
		}

		printf ("Interpolation space sampled\n");

		// Evaluate the cost on the points of the interpolations space
		for (size_t j = 0; j < N; j++)
		{
			// Form current point by adding the step to the current iterate xc
			double pc[R + model->nin];
			pc[R+1]=0.0;
			for (size_t i = 0; i < R; i++)
			{
				pc[i] = xc[i] + gsl_matrix_get (Pt, i, j);
			}

			// Evaluate the metric at the current point
			double cumD = parallel_sim (pc, r);

			// Put the cost value in the appropriate spot in the vector BB
			gsl_vector_set (BB, j, cumD);

			// Check if we randomly found a good point
			if (cumD == 0.0)
			{
				printf ("The point\n");
				print_matrix (pc, R, 1);
				printf ("satisfies the metric condition.\n");
				break;
			}

			if (j == 0)
			{
				fxc = cumD;
			}

			printf ("%d ", (int) (j+1));
			fflush (stdout);
		}

		// Create the model matrix
		model_matrix (AA, Pt);

		// Fit the quadratic model
		quadratic_model (H, g, &c, AA, BB);

		/*
		// Check if the model and the cost match on the interp space
		printf ("\n");
		for (size_t l = 0; l < N; l++)
		{
			gsl_vector_view point = gsl_matrix_column (Pt, l);
			double mm = quad (H, g, c, &point.vector);
			double cc = gsl_vector_get (BB, l);
			printf ("Model: %f - Cost: %f\n", mm, cc);
		)
		break;
		*/

		// Extract PSD part of H
		psd_part (Hpsd, H);

		// Compute the new gradient
		gsl_vector_view xcv = gsl_vector_view_array (xc, R);
		gsl_blas_dgemv (CblasNoTrans, 1.0, H, &xcv.vector, 0.0, gp);
		gsl_vector_add (gp, g);

		// Set parameters for cvxgen
		for (size_t i = 0; i < H->size1; i++)
		{
			for (size_t j = 0; j < H->size2; j++)
			{
				// Hessian
				params.H[i+j*(H->size1)] = gsl_matrix_get (Hpsd, i, j);
			}
			// Gradient
			params.g[i] = gsl_vector_get (g, i);
			// Constant term
			params.c[0] = c;
			// Box constraints
			params.lb[i] = lb[i] - xc[i];
			params.ub[i] = ub[i] - xc[i];
			// Trust region bound
			params.d[0] = D;
		}

		// Compute solution with cvxgen
		printf ("\n");
		size_t niter = solve ();

		// Check that the solution was computed correctly
		if (work.converged != 1)
		{
			printf ("\n\nQuadratic subproblem solution not valid...\n\n");
			continue;
		}

		// Form the next candidate iterate
		double xnext[R];
		for (size_t k = 0; k < R; k++)
		{
			xnext[k] = xc[k] + vars.p[k];
		}

		// Evaluate effect of step
		double fxn = parallel_sim (xnext, r);
		double fqxn = quad (H, g, c, vars.p);
		double rho = (fxc - fxn)/(c - fqxn);

		printf ("\n                         FUNCTION\t\tMODEL\n");
		printf ("Cost at current point:   %f\t\t%f\n", fxc, c);
		printf ("Cost at curr pt + step:  %f\t\t%f\n", fxn, fqxn);
		printf ("Reduction ratio: %f\n\n", rho);

		// Decide what to do
		if (rho <= 0.0)
		{
			// Reject the step and decrese the trust region size
			D = (D/2 > 0.01) ? D/2 : 0.01;
			printf ("Model is unreliable, shriking trust region to D = %f\n\n", D);
		}
		else if ((rho > 0.0) && (rho <= 0.75))
		{
			// Accept the step and leave the trust region alone
			// Take the step
			for (size_t k = 0; k < R; k++)
			{
				xp[k] = xc[k];
				xc[k] = xc[k] + vars.p[k];
			}

			printf ("Step accepted. New iterate:\n");
			pm (xc, R, 1);
		}
		else
		{
			// Accept the step and increase the trust region size
			// Take the step
			for (size_t k = 0; k < R; k++)
			{
				xp[k] = xc[k];
				xc[k] = xc[k] + vars.p[k];
			}

			D = D*2;
			printf ("Expanding trust region to D = %f\n", D);
			printf ("Step accepted. New iterate:\n");
			pm (xc, R, 1);
		}

		// Increase iteration counter
		iter++;
	}
	while (1);

	// The loop is done

	// Tell all the slaves to exit by sending an empty message with the DIETAG
	for (rank = 1; rank < ntasks; rank++)
	{
		MPI_Send (0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
	}

	// Save the results

	/*
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
	*/

	// Free manually allocated resources

	gsl_matrix_free (Pt);
	gsl_matrix_free (AA);
	gsl_vector_free (BB);
	gsl_matrix_free (H);
	gsl_matrix_free (Hpsd);
	gsl_vector_free (g);
	gsl_vector_free (gp);
	gsl_matrix_free (I);
	gsl_rng_free (r);

	return;
}


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


// This function evaluates the metric at the point in parameter space theta.
static double parallel_sim (double theta[], const gsl_rng * r)
{
	// Declare some necessary variables
	int ntasks, rank, count = 0;
	MPI_Status status;
	size_t res_no = (model->nspecies)*(times->size);
	double results[res_no];

	// Find out how many processes there are in the default communicator
	MPI_Comm_size (MPI_COMM_WORLD, &ntasks);

	// Allocate rxn_ensemble object to store SCRN simulations
	rxn_ensemble * res = rxn_ensemble_alloc (S, model->nspecies, times->size);
	// Allocate matrix to store SCRN output
	gsl_matrix * counts = gsl_matrix_calloc (res->nreplic, model->nout*res->ntimes);
	// Declare necessary views
	gsl_vector_view drow, crow;

	// Request simulations from the slaves by sending the current parameters to simulate
	for (rank = 1; rank < ntasks; rank++)
	{
		MPI_Send (theta, 							// current parameters
				  model->nparams + model->nin,		// # of data items
				  MPI_DOUBLE,			    		// data items are of type double
				  rank,								// destination process rank
				  WORKTAG + count,					// user chosen message tag
				  MPI_COMM_WORLD);					// default communicator

		count++;
	}

	// Keep requesting simulations untils we have enough
	while (count < S)
	{
		// Receive results from a slave
		MPI_Recv (results,      	   // message buffer
				  res_no,              // # of data items
				  MPI_DOUBLE,          // of type double real
				  MPI_ANY_SOURCE,      // receive from any sender
				  MPI_ANY_TAG,     	   // any type of message
				  MPI_COMM_WORLD,      // default communicator
				  &status);            // info about the received message

		// Process the results
		for (size_t i = 0; i < res_no; i++)
		{
			// The TAG from the sender is simulation #
			res->data[status.MPI_TAG]->counts->data[i] = results[i];
		}

		// Send the slave a new work unit
		MPI_Send (theta,							// message buffer
				  model->nparams + model->nin,		// # of data items
				  MPI_DOUBLE,						// data items are of type double
				  status.MPI_SOURCE,				// destination process rank (the slave we just received from)
				  WORKTAG + count,					// user chosen message tag
				  MPI_COMM_WORLD);					// default communicator

		// Get the next unit of work to be done
		count++;
	}

	// We now have enough simulations, so receive all the outstanding results from the slaves
	for (rank = 1; rank < ntasks; rank++)
	{
		MPI_Recv (results, res_no, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		// Process the results
		for (size_t i = 0; i < res_no; i++)
		{
			res->data[status.MPI_TAG]->counts->data[i] = results[i];
		}
	}

	// We are now ready to evaluate the metric

	// Extract the output
	rxn_ensemble_counts (res, counts, out);

	// Compute the flourescence levels
	for (size_t i = 0; i < counts->size1; i++)
	{
		for (size_t j = 0; j < counts->size2; j++)
		{
			// Resample a background level
			double sambg = gsl_matrix_get (bg, gsl_rng_uniform_int (r, M), 0);

			// Current count
			double z = gsl_matrix_get (counts, i, j);

			// Sample a fluorescence level
			double samfl = mufp*z + gsl_ran_gaussian (r, sqrt (sigfp*sigfp*z));

			gsl_matrix_set (counts, i, j, samfl + sambg);
		}
	}

	// Compute the Kolmogorov distances
	double cumD = 0, dc = 0;
	for (size_t k = 1; k < data->size2; k++)
	{
		drow = gsl_matrix_column (data, k);
		crow = gsl_matrix_column (counts, k);
		dc = ksdist_two_sample_2 (&crow.vector, &drow.vector);
		cumD += (dc > epsilon) ? dc : 0.0;
	}

	// Free manually allocated resources
	rxn_ensemble_free (res);
	gsl_matrix_free (counts);

	return cumD;
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
				// Keep reading until the desired parameter is > 1.0
				do
				{
					fread_floats_swap (fbuf, par, fcs);
				}
				while (fbuf[p-1]<=1.1);

				// Place the data in the proper spot
				gsl_matrix_set (data, m, P*(j-1)+l-1, log10 (fbuf[p-1]));
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
