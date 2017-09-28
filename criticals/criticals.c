/*
 *  criticals.c
 *  mpiStochEst
 *
 *	Evaluation of critical number of simulations using exact Kolmogorov distribution
 *
 *  Created by Gabriele Lillacci in March 2012.
 *	Latest revision: March 2012.
 *
 *
 *	This free software is available under the Creative Commons Attribution Share Alike License.
 *	You are permitted to use, redistribute and adapt this software as long as appropriate credit
 *	is given to the original author, and all derivative works are distributed under the same
 *	license or a compatible one.
 *	For more information, visit http://creativecommons.org/licenses/by-sa/3.0/ or send a letter to
 *	Creative Commons, 171 2nd Street, Suite 300, San Francisco, California, 94105, USA.
 */

#include <stdio.h>
#include <string.h>

#include <parestlib.h>

// Macros
#define EPS 2.2204e-16	// Machine precision


int main (int argc, char * argv[])
{
	// Specify some parameters
	size_t npoints = 601;
	int M = 80000;
	double beta = 0.05;
	double start = 0.04;
	double end = 0.1;

	// Find alpha
	double alpha = 1 - sqrt (1 - beta);

	// Perform the first inversion of K_M
	double inv1 = kolmogorov_cdf_inverse (M, 1-alpha);

	printf ("\n\nThe result of the first inversion is %f\n\n", inv1);

	return EXIT_SUCCESS;

	gsl_vector * S = gsl_vector_alloc (601);

	for (size_t i = 0; i < npoints; i++)
	{
		double epsilon = start + ((double) i)*(end-start)/((double) npoints-1);
		printf ("epsilon = %f\n", epsilon);
		int Sc = kolmogorov_cdf_inverse_M (epsilon - inv1, 1 - alpha);
		gsl_vector_set (S, i, Sc);
	}

	printf ("\n\nThe final result is:\n");

	for (size_t i = 0; i < npoints; i++)
	{
		printf ("%d, ", (int) gsl_vector_get (S, i));
	}

	printf("\n\n");

	gsl_vector_free (S);

	return EXIT_SUCCESS;
}
