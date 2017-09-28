/*
 *  auxfun.c
 *
 *  Auxiliary functions for mpiFind0
 *
 *  Created by Gabriele Lillacci in February 2012.
 *	Latest revision: February 2012.
 *
 *
 *	This free software is available under the Creative Commons Attribution Non-Commercial Share Alike License.
 *	You are permitted to use, share, copy, redistribute and adapt this software as long as appropriate credit
 *	is given to the original author, all derivative works are distributed under the same license or a compatible one,
 *	and this software and its derivatives are not used for commercial purposes.
 *	For more information, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or contact
 *	Creative Commons, 171 2nd Street, Suite 300, San Francisco, California, 94105, USA.
 */


#include <auxfun.h>


// This function computes the model matrix corresponding to a set of interpolation points
int model_matrix (gsl_matrix * AA, const gsl_matrix * P)
{
    // The matrix P contains R points, each of dimension N.
    size_t N = P->size1;
    size_t R = P->size2;

    // Scan the rows of A
    for (size_t r = 0; r < R; r++)
    {
        // Counter for quadratic terms
        size_t cnt = 0;

        // Quadratic terms
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = i; j < N; j++)
            {
                double cv = (i==j ? 0.5 : 1.0) * gsl_matrix_get (P, i, r) * gsl_matrix_get (P, j, r);
                gsl_matrix_set (AA, r, cnt, cv);
                cnt++;
            }
        }

        // Linear terms
        for (size_t i = 0; i < N; i++)
        {
            double cv = gsl_matrix_get (P, i, r);
            gsl_matrix_set (AA, r, i+R-N-1, cv);
        }

        // Constant term
        gsl_matrix_set (AA, r, R-1, 1.0);
    }

    // Signal that the computation was completed correctly
    return GSL_SUCCESS;
}


// This function fits the quadratic model m(x) = 1/2*x'*H*x + g'*x + c
int quadratic_model (gsl_matrix * H, gsl_vector * g, double * c, const gsl_matrix * AA, const gsl_vector * BB)
{
    // Detect the number of optimization variables
    size_t N = H->size1;
    // Detect the number of interpolations points
    size_t R = AA->size1;

    // Make a copy of the matrix AA
    gsl_matrix * AAc = gsl_matrix_alloc (R, R);
    gsl_matrix_memcpy (AAc, AA);

    // Create a vector to store the solution, i.e. the model coefficients
    gsl_vector * q = gsl_vector_alloc (R);

    // Solve the linear system in the model coefficients AA*q=BB
    gsl_linalg_HH_solve (AAc, BB, q);

    // Initialize a counter for the elements of q
    size_t cnt = 0;

    // Extract the quadratic coefficients and put them into the matrix H
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = i; j < N; j++)
        {
            gsl_matrix_set (H, i, j, gsl_vector_get (q, cnt));
            gsl_matrix_set (H, j, i, gsl_vector_get (q, cnt));
            cnt++;
        }
    }

    // Extract the linear coefficients and put them into the vector g
    for (size_t i = 0; i < N; i++)
    {
        gsl_vector_set (g, i, gsl_vector_get (q, cnt));
        cnt++;
    }

    // Extract the constant term and put into c
    * c = gsl_vector_get (q, R-1);

    // Free all the manually allocated resources
    gsl_matrix_free (AAc);
    gsl_vector_free (q);

    // Signal that the computation was completed correctly
    return GSL_SUCCESS;
}

// This function evaulates the quadratic form 1/2*x'*H*x + g'*x + c
double quad (const gsl_matrix * H, const gsl_vector * g, double c, const double x[])
{
	// Allocate necessary objects
	gsl_vector_const_view xview = gsl_vector_const_view_array (x, g->size);
	gsl_vector * y = gsl_vector_alloc (g->size);

	// Compute the product y = H^T*x
	gsl_blas_dgemv (CblasTrans, 1.0, H, &xview.vector, 0.0, y);

	// Compute the dot product r1 = y^T*x = (H^T*x)^T*x = x^T*H*x
	double r1 = 0.0;
	gsl_blas_ddot (y, &xview.vector, &r1);

	// Compute the dot product r2 = g^T*x
	double r2 = 0.0;
	gsl_blas_ddot (g, &xview.vector, &r2);

	// Free all the manually allocated resources
	gsl_vector_free (y);

	// Return
	return 0.5*r1 + r2 + c;
}


// This function computes an estimate of the smallest eigenvalue of the matrix A (based on the bound of Zhan 2006)
// ATTENTION!! A must be REAL and SYMMETRIC.
double lambda_min (const gsl_matrix * A)
{
    // Detect the dimension of A
    size_t N = A->size1;

    // Find the interval that contains the elements of the matrix A
    double a = gsl_matrix_min (A);
    double b = gsl_matrix_max (A);

    double lambda1 = 0.0;

    if (fabs(a)>=b)
        lambda1 = ((double) N)*a;
    else
    {
        if (N%2 == 0)
        {
            // N is even
            lambda1 = ((double) N)*(a-b)/2;
        }
        else
        {
            // N is odd
            lambda1 = (((double) N)*a - sqrt (a*a + (((double) N*N) - 1))*b*b)/2;
        }
    }

    // Return
    return lambda1;
}


// This function computes an estimate of the largest eigenvalue of the matrix A (based on the bound of Zhan 2006)
// ATTENTION!! A must be REAL and SYMMETRIC.
double lambda_max (const gsl_matrix * A)
{
    // Detect the dimension of A
    size_t N = A->size1;

    // Find the interval that contains the elements of the matrix A
    double a = gsl_matrix_min (A);
    double b = gsl_matrix_max (A);

    double lambdaN = 0.0;

    if (a >= -fabs(b))
        lambdaN = ((double) N)*b;
    else
    {
        if (N%2 == 0)
        {
            // N is even
            lambdaN = ((double) N)*(b-a)/2;
        }
        else
        {
            // N is odd
            lambdaN = (((double) N)*b - sqrt (b*b + (((double) N*N) - 1))*a*a)/2;
        }
    }

    // Return
    return lambdaN;
}


// This function computes the positive semidefinite part of a real symmetric matrix A.
int psd_part (gsl_matrix * Apsd, const gsl_matrix * A)
{
	// Detect the size of the square matrices
	size_t N = A->size1;

	// Make a copy of the matrix A
	gsl_matrix * Ac = gsl_matrix_alloc (N, N);
	gsl_matrix_memcpy (Ac, A);

	// Allocate eigenvalue workspace
	gsl_eigen_symmv_workspace * wsp = gsl_eigen_symmv_alloc (N);

	// Allocate objects to store eigenvalues and eigenvector
	gsl_vector * eval = gsl_vector_alloc (N);
	gsl_matrix * V = gsl_matrix_alloc (N, N);
	gsl_matrix * Dp = gsl_matrix_alloc (N, N);

	// Compute the eigenvalues and eigenvectors
	gsl_eigen_symmv (Ac, eval, V, wsp);

	// Reset the matrix Dp
	gsl_matrix_set_all (Dp, 0.0);

	// Scan for the non-zero eigenvalues and put them in the matrix Dp
	for (size_t i = 0; i < N; i++)
	{
		double evc = gsl_vector_get (eval, i);
		gsl_matrix_set (Dp, i, i, (evc > 1.0e-10) ? evc : 0.0);
	}

	// Compute the product C = V*Dp.
	gsl_matrix * C = gsl_matrix_alloc (N, N);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, V, Dp, 0.0, C);

	// Compute the product Apsd = V*Dp*V' = C*V'
	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, C, V, 0.0, Apsd);

	// Free the allocated resources
	gsl_matrix_free (C);
	gsl_matrix_free (Dp);
	gsl_matrix_free (V);
	gsl_vector_free (eval);
	gsl_eigen_symmv_free (wsp);
	gsl_matrix_free (Ac);

	// Return
	return GSL_SUCCESS;
}


// This function samples an interpolation space, subject to the box constraints lb and ub
// and to a trust region bound D on the step. The trust region is given in infinity norm.
int interpolation_space (gsl_matrix * Pt, double xc[], double lb[], double ub[], double D, const gsl_rng * r)
{
	// Detect number of points
	size_t N = Pt->size2;
	// Detect number of components in each point
	size_t R = Pt->size1;

	// Reset the matrix Pt
	gsl_matrix_set_all (Pt, 0.0);

	// Trust region scale factor
	int s = 10;

	// The first point in the space is the current iterate (step = 0)
	for (size_t i = 0; i < R; i++)
	{
		for (size_t j = 1; j < N; j++)
		{
			// Find the largest lower bound between the box contraints and the trust region
			double lbc = (lb[i] > (xc[i]-s*D)) ? lb[i] : (xc[i]-s*D);

			// Find the smallest upper bound between the box contraints and the trust region
			double ubc = (ub[i] < (xc[i]+s*D)) ? ub[i] : (xc[i]+s*D);

			// Generate a uniform random number with these bounds
			double rnd = ran_uniform (r, lbc, ubc);

			// The step is this random number minus the current value of the iterate
			gsl_matrix_set (Pt, i, j, rnd - xc[i]);
		}
	}

	// Allocate objects for the computation of the SVD
	gsl_matrix * U = gsl_matrix_alloc (Pt->size2, Pt->size1);
	gsl_matrix * V = gsl_matrix_alloc (Pt->size1, Pt->size1);
	gsl_vector * S = gsl_vector_alloc (Pt->size1);

	// Make a transpose copy of Pt
	gsl_matrix_transpose_memcpy (U, Pt);

	// Compute the SVD
	gsl_linalg_SV_decomp_jacobi (U, V, S);

	// Scan for singular values that are close to zero
	int success = GSL_SUCCESS;
	for (size_t k = 0; k < S->size; k++)
	{
		if (gsl_vector_get (S, k) <= 1.0e-8)
			success = GSL_EFAILED;
	}

	// Free the manually allocated resources
	gsl_vector_free (S);
	gsl_matrix_free (V);
	gsl_matrix_free (U);

	// Return
	return success;
}


// This function computes the reciprocal condition number of the matrix A using the singular value decomposition.
// This function is more precise and works for general matrices, but it's much more expensive to evaluate.
double rcond_svd (const gsl_matrix * A)
{
    // Allocate the necessary objects for the SVD
    gsl_matrix * U = gsl_matrix_alloc (A->size1, A->size2);
    gsl_matrix * V = gsl_matrix_alloc (A->size2, A->size2);
    gsl_vector * S = gsl_vector_alloc (A->size2);

    // Copy the A matrix into U
    gsl_matrix_memcpy (U, A);

    // Compute the SVD
    gsl_linalg_SV_decomp_jacobi (U, V, S);

    // Compute the reciprocal condition number
    double k = gsl_vector_get (S, (S->size)-1) / gsl_vector_get (S, 0);

    // Free the manually allocated resources
    gsl_matrix_free (U);
    gsl_matrix_free (V);
    gsl_vector_free (S);

    // Return
    return k;
}


// Evaluate the uniform density
double unif_dist (double x, double a, double b)
{
	return (x >= a && x <=b) ? 1/(b-a) : 0.0;
}


// Generate a normal random number with mean m and ***variance*** v
double ran_gaussian (const gsl_rng * r, double m, double v)
{
	return m + gsl_ran_gaussian (r, sqrt (v));
}


// Generate a uniform random number between a and b
double ran_uniform (const gsl_rng * r, double a, double b)
{
	return a + (b-a) * gsl_rng_uniform (r);
}
