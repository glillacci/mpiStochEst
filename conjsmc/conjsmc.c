/*
 *  conjsmc.c
 *  mpiStochEst
 *
 *  Conjugate test for INSIGHT SMC sampler
 *  This program propagates a population of particles from one tolerance to the next.
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
#include <math.h>

// GSL library includes
#include <gsl/gsl_sf_gamma.h>

// Custom library includes
#include <parestlib.h>

// Auxiliary function declarations
double rinvgamma (const gsl_rng * r, double a, double b);
double rnorm (const gsl_rng * r, double m, double v);
double runif (const gsl_rng * r, double a, double b);
double dinvgamma (double x, double a, double b);
double dnorm (double x, double m, double v);
double dunif (double x, double a, double b);
double rperturb (const gsl_rng * r, double theta, double wd);
double dperturb (double x, double theta, double wd);
double rperturb2 (const gsl_rng * r, double theta, double wd);
double dperturb2 (double x, double theta, double wd);


int main (int argc, char * argv[])
{
	// ********************************************** USER-SPECIFIED STUFF *****************************************************
	// Number of populations
	size_t T = 14;
	// Number of particles per population
	size_t U = 10000;

	// Set tolerance schedule
	double epsilon[] = {0.5, 0.4, 0.3, 0.2, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01};

	// Set "external" confidence level
	double b = 0.01;

	// Set prior parameters
	double nu_0 = 30;
	double sig_0 = 200;
	double mu_0 = 200;
	double k_0 = 0.1;

	// Set number of experimental samples
	size_t M = 5000;


	// *************************************************** INITIALIZATION *******************************************************
	// Set up the random number generator
	gsl_rng_env_setup ();
	gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (r, (unsigned long int) time (NULL));

	// Calculate "internal" confidence level
	double a = 1 - sqrt (1-b);

	// Invert the Kolmogorov distribution
	double inv = kolmogorov_cdf_inverse (M, 1-a);

	// Generate the "experimental samples"
	double Y[M];
	for (size_t i = 0; i < M; i++)
	{
		Y[i] = rnorm (r, 100, 100);
		printf("%f\n", Y[i]);
	}
	gsl_vector_view Yv = gsl_vector_view_array (Y, M);
	gsl_sort_vector (&Yv.vector);

	// Declare arrays for mean and variance
	double muc[U], mup[U];
	double sigc[U], sigp[U];

	// Declare arrays for the weights
	double wc[U], wp[U];

	// Initialize the arrays
	for (size_t u = 0; u < U; u++)
	{
		sigc[u] = 0.0;
		sigp[u] = 0.0;
		muc[u] = 0.0;
		mup[u] = 0.0;
		wc[u]= 0.0;
		wp[u]= 0.0;
	}


	// *************************************************** MAIN SMC LOOP ********************************************************
	// Loop for tolerances
	for (size_t t = 0; t < T; t++)
	{
		// Reset the particle counter
		size_t u = 0;

		// Compute number of simulations for current stage
		//size_t S = (size_t) ceil (-log(a/2)/(2*pow(epsilon[t] - sqrt (-1/(2*(double)M) * log (a/2)), 2)));
		size_t S = (size_t) ceil (-log(a/2)/(2*pow(epsilon[t] - inv, 2)));
		S = (S>10000) ? 10000 : S;
		//if (epsilon[t]<=0.05 && S<10000)
			//S = 10000;

		// Print stage info
		printf ("\n\nStage %d: epsilon = %f, S = %d\n", (int) t+1, epsilon[t], (int) S);

		// Declare array for the simulations
		double X[S];
		gsl_vector_view Xv;

		// Declare a variable for the sum of the weights
		double wsum = 0.0;

		// Calculate the variances of the previous populations
		double vm = gsl_stats_variance (mup, 1, U);
		double vs = gsl_stats_variance (sigp, 1, U);

		printf ("vm = %f, vs = %f\n", vm, vs);

		// Loop for particles
		do
		{
			// Propose a new particle
			if (t == 0)
			{
				// Sample a particle from the priors
				sigc[u] = rinvgamma (r, nu_0/2, nu_0*sig_0/2);
				muc[u] = rnorm (r, mu_0, sigc[u]/k_0);
			}
			else
			{
				// Sample from the previous population
				do
				{
					// Generate a uniform random number
					double uu = gsl_rng_uniform (r);

					// Accumulate the weights until they are larger than uu
					double mu = 0.0;
					size_t ix = 0;
					while ((mu < uu) && (ix < (U-1)))
					{
						mu += wp[ix];
						ix++;
					}

					// Take the particle corresponding to ix
					muc[u] = mup[ix];
					sigc[u] = sigp[ix];

					// Introduce perturbation
					//sigc[u] = rperturb (r, sigc[u], epsilon[t]);
					//muc[u] = rperturb (r, muc[u], epsilon[t]);
					sigc[u] = rperturb2 (r, sigc[u], 2*vs);
					muc[u] = rperturb2 (r, muc[u], 2*vm);

					// Visualize resampled particle
					// printf ("Resampled particle: [%f, %f]\n", muc[u], sigc[u]);

					// We must ensure that the perturbation does not violate the prior!
				}
				while (sigc[u] <= 0.0);
			}

			// Generate data using the proposed particle
			for (size_t s = 0; s < S; s++)
			{
				X[s] = rnorm (r, muc[u], sigc[u]);
			}

			// Calculate the Kolmogorov distance
			Xv = gsl_vector_view_array (X, S);
			gsl_sort_vector (&Xv.vector);
			double dSM = ksdist_two_sample_nr_presort (&Xv.vector, &Yv.vector);

			// Decide on acceptance/rejection
			if (dSM <= epsilon[t])
			{
				// Particle is accepted

				// Compute the weight
				if (t == 0)
				{
					wc[u] = 1.0;
					wsum += wc[u];
				}
				else
				{
					double num = dnorm (muc[u], mu_0, sigc[u]/k_0) * dinvgamma (sigc[u], nu_0/2, nu_0*sig_0/2);
					double den = 0.0;
					for (size_t l = 0; l < U; l++)
					{
						// den += wp[l] * dperturb (muc[u], mup[l], epsilon[t]) * dperturb (sigc[u], sigp[l], epsilon[t]);
						den += wp[l] * dperturb2 (muc[u], mup[l], 2*vm) * dperturb2 (sigc[u], sigp[l], 2*vs);
					}
					wc[u] = num/den;
					wsum += wc[u];
				}

				// Debug info
				// printf ("Accepted particle: [%f, %f] - Weight: %f\n", muc[u], sigc[u], wc[u]);

				if ((u+1) % 100 == 0)
				{
					printf (" %d ", (int) u+1);
					fflush (stdout);
				}

				// Increase accepted particle counter
				u++;
			}
		}
		while (u < U);

		printf ("\n wsum=%f\n", wsum);

		for (size_t u = 0; u < U; u++)
		{
			// Normalize the weights
			wc[u] = wc[u]/wsum;

			// Copy "current" arrays into "previous" arrays
			sigp[u] = sigc[u];
			mup[u] = muc[u];
			wp[u] = wc[u];
		}
	}

	// ***************************************************** OUTPUT *************************************************************
	// Write to file
	FILE * Yfile = fopen ("./conjsmc/conjsmc_Y.txt", "w");
	if (Yfile == NULL || gsl_vector_fprintf (Yfile, &Yv.vector, "%f"))
	{
		printf (">> error writing Y output.. program exiting\n");
		return EXIT_FAILURE;
	}
	fclose (Yfile);

	// Write to file
	FILE * mfile = fopen ("./conjsmc/conjsmc_mu.txt", "w");
	gsl_vector_view muv = gsl_vector_view_array (muc, U);
	if (mfile == NULL || gsl_vector_fprintf (mfile, &muv.vector, "%f"))
	{
		printf (">> error writing mu output.. program exiting\n");
		return EXIT_FAILURE;
	}
	fclose (mfile);

	// Write to file
	FILE * sfile = fopen ("./conjsmc/conjsmc_sig.txt", "w");
	gsl_vector_view sigv = gsl_vector_view_array (sigc, U);
	if (sfile == NULL || gsl_vector_fprintf (sfile, &sigv.vector, "%f"))
	{
		printf (">> error writing sigma output.. program exiting\n");
		return EXIT_FAILURE;
	}
	fclose (sfile);

	// Write to file
	FILE * wfile = fopen ("./conjsmc/conjsmc_w.txt", "w");
	gsl_vector_view wv = gsl_vector_view_array (wc, U);
	if (wfile == NULL || gsl_vector_fprintf (wfile, &wv.vector, "%f"))
	{
		printf (">> error writing weights output.. program exiting\n");
		return EXIT_FAILURE;
	}
	fclose (wfile);


	// ***************************************************** CLEAN UP ***********************************************************
	// Free the random number generator
	gsl_rng_free (r);

	// Return
	printf("\n\n");
	return EXIT_SUCCESS;
}


// Generate a random number from the inverse Gamma distribution with shape parameter a and scale parameter b
double rinvgamma (const gsl_rng * r, double a, double b)
{
	return 1/gsl_ran_gamma (r, a, 1/b);
}


// Generate a normal random number with mean m and ***variance*** v
double rnorm (const gsl_rng * r, double m, double v)
{
	return m + gsl_ran_gaussian (r, sqrt (v));
}


// Generate a uniform random number between a and b
double runif (const gsl_rng * r, double a, double b)
{
	return a + (b-a) * gsl_rng_uniform (r);
}


// Evaluate the inverse gamma pdf with parameters a and b at the point x
double dinvgamma (double x, double a, double b)
{
	return gsl_sf_gammainv (a) * pow (b, a) * pow (x, -a-1) * exp (-b/x);
}


// Evaluate the normal pdf with mean m and ***variance*** v at the point x
double dnorm (double x, double m, double v)
{
	return gsl_ran_gaussian_pdf (x - m, sqrt (v));
}


// Evaluate the uniform density
double dunif (double x, double a, double b)
{
	return (x >= a && x <=b) ? 1/(b-a) : 0;
}


// Generate a randomly perturbed resampled particle
double rperturb (const gsl_rng * r, double theta, double wd)
{
	return runif (r, theta * (1 - wd), theta * (1 + wd));
}


// Evaluate the perturbation density
double dperturb (double x, double theta, double wd)
{
	return dunif (x, theta * (1 - wd), theta * (1 + wd));
}


// Generate a randomly perturbed resampled particle
double rperturb2 (const gsl_rng * r, double theta, double wd)
{
	return rnorm (r, theta, wd);
}


// Evaluate the perturbation density
double dperturb2 (double x, double theta, double wd)
{
	return dnorm (x, theta, wd);
}
