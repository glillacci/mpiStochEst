#	Autoconf script for
#
#	MPISTOCHEST
#
# 
#	Created by Gabriele Lillacci in January 2012.
#	Latest revision: January 2012.
#
#
#	This free software is available under the Creative Commons Attribution Non-Commercial Share Alike License.
#	You are permitted to use, share, copy, redistribute and adapt this software as long as appropriate credit
#	is given to the original author, all derivative works are distributed under the same license or a compatible one,
#	and this software and its derivatives are not used for commercial purposes.
#	For more information, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or contact
#	Creative Commons, 171 2nd Street, Suite 300, San Francisco, California, 94105, USA.


AC_PREREQ(2.59)
AC_INIT(mpistochest, 1.0)

AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE(foreign -Wall -Werror)

AC_CONFIG_MACRO_DIR([m4])

# Look for mpicc and set it to C99 mode
AC_PROG_CC(mpicc)
AC_PROG_CC_C99

# Check for required libraries
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
AC_CHECK_LIB([parestlib],[rxn_sample_path_alloc])
AC_CHECK_LIB([stochmod],[syncirc_mod_setup])
AC_CHECK_LIB([xml2],[xmlReadFile])

AC_CONFIG_HEADER(config.h)
AC_CONFIG_FILES(Makefile mpiabck/Makefile mpismc/Makefile mpifind0/Makefile criticals/Makefile fcsread/Makefile gridtest/Makefile conjsmc/Makefile)

AC_OUTPUT