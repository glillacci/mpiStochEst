/*
 *  criticals.c
 *  mpiStochEst
 *
 *	Evaluation of critical number of simulations using exact Kolmogorov distribution
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
