/*
 *  auxfun.h
 *
 *  Auxiliary functions for mpiFind0 (header file)
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

/*
 Includes
 */

#include <stdio.h>
#include <parestlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>


/*
 Exported functions prototype declarations
 */

int model_matrix (gsl_matrix * AA, const gsl_matrix * P);
int quadratic_model (gsl_matrix * H, gsl_vector * g, double * c, const gsl_matrix * AA, const gsl_vector * BB);
double quad (const gsl_matrix * H, const gsl_vector * g, double c, const double x[]);
double lambda_min (const gsl_matrix * A);
double lambda_max (const gsl_matrix * A);
double rcond_svd (const gsl_matrix * A);
double unif_dist (double x, double a, double b);
double ran_gaussian (const gsl_rng * r, double m, double v);
double ran_uniform (const gsl_rng * r, double a, double b);
int psd_part (gsl_matrix * Apsd, const gsl_matrix * A);
int interpolation_space (gsl_matrix * Pt, double xc[], double lb[], double ub[], double D, const gsl_rng * r);
