# Autoconf script for mpiStochEst
#
# This file is part of mpiStochEst.
# Copyright 2011-2017 Gabriele Lillacci.
#
# mpiStochEst is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# mpiStochEst is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with mpiStochEst.  If not, see <http://www.gnu.org/licenses/>.


AC_PREREQ(2.59)
AC_INIT(mpistochest, 1.0)

AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE(foreign -Wall -Werror)

AC_CONFIG_MACRO_DIR([m4])

# Look for mpicc and set it to C99 mode
AC_PROG_CC(mpicc)
AC_PROG_CC_C99

# Check for required libraries
AC_SEARCH_LIBS([cos], [m], [], [
  AC_MSG_ERROR([unable to find the m library])
])
AC_SEARCH_LIBS([cblas_dgemm], [gslcblas], [], [
  AC_MSG_ERROR([unable to find the gslcblas library])
])
AC_SEARCH_LIBS([gsl_blas_dgemm], [gsl], [], [
  AC_MSG_ERROR([unable to find the gsl library])
])
AC_SEARCH_LIBS([rxn_sample_path_alloc], [parestlib], [], [
  AC_MSG_ERROR([unable to find the parestlib library])
])
AC_SEARCH_LIBS([syncirc_mod_setup], [stochmod], [], [
  AC_MSG_ERROR([unable to find the stochmod library])
])
AC_SEARCH_LIBS([xmlReadFile], [xml2], [], [
  AC_MSG_ERROR([unable to find the xml2 library])
])

AC_CONFIG_HEADER(config.h)
AC_CONFIG_FILES(Makefile mpiabck/Makefile mpismc/Makefile mpifind0/Makefile criticals/Makefile fcsread/Makefile gridtest/Makefile conjsmc/Makefile)

AC_OUTPUT
