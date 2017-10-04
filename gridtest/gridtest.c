/*
 *  grditest.c
 *  mpiStochEst
 *
 *  Test program for GRID
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
