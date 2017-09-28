/*
 *  grditest.c
 *  mpiStochEst
 *
 *  Test program for GRID
 *
 *  Created by Gabriele Lillacci in May 2012.
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

// GSL includes
#include <gsl/gsl_matrix.h>

// Mpi include
#include <mpi.h>


int main (int argc, char * argv[])
{
	// Test standard output
	printf ("Hello, standard output!\n");

	// Test standard error
	fprintf (stderr, "Hello, standard error!\n");

	// Test the GSL library
	gsl_matrix * m = gsl_matrix_alloc (3, 3);
	gsl_matrix_set_identity (m);
	printf ("The element (3,3) of the matrix m is: %f.\n", gsl_matrix_get (m, 2, 2));
	gsl_matrix_free (m);

	// Test OpenMPI
	MPI_Init (&argc, &argv);
	int ntasks;
	MPI_Comm_size (MPI_COMM_WORLD, &ntasks);
	printf ("There are %d processes in this parallel run.\n", ntasks);
	MPI_Finalize ();

	return EXIT_SUCCESS;
}
