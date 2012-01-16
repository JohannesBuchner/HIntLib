## Process this file with automake to produce Makefile.in

COMMON_TESTS = \
   abs.test \
   partition.test \
   bitop.test \
   polynomial.test \
   eratosthenes.test \
   prime.test \
   arithmetic.test \
   linearalgebra.test \
   genmatrix.test \
   shiftnet.test \
   tparameter.test \
   prng.test \
   rule.test
FORTRAN_TESTS = \
   compare_with_halton_in_toms659.test \
   compare_with_faure_in_toms659.test \
   compare_with_sobol_in_toms659.test \
   compare_with_nied_2_in_toms738.test \
   compare_with_nied_prime_in_toms738.test \
   compare_with_nied_pp_in_toms738.test
TAIL_SORT_TESTS = qrnsequence.test

EXTRA_DIST = $(TESTS) $(TAIL_SORT_TESTS) $(FORTRAN_TESTS)

if HAVE_FORTRAN
USED_FORTRAN_TESTS = $(FORTRAN_TESTS)
SUBDIRS = toms659 toms738
else
USED_FORTRAN_TESTS =
endif

if HAVE_TAIL_SORT
USED_TAIL_SORT_TESTS = $(TAIL_SORT_TESTS)
else
USED_TAIL_SORT_TESTS =
endif

TESTS = $(COMMON_TESTS) $(USED_FORTRAN_TESTS) $(USED_TAIL_SORT_TESTS)

check_PROGRAMS = \
   numericdiff \
   test_abs \
   test_arithmetic \
   test_bitop \
   test_eratosthenes \
   test_genmatrix \
   test_linearalgebra \
   test_partition \
   test_polynomial \
   test_prng \
   test_prime \
   test_qrnsequence \
   test_rule \
   test_shiftnet \
   test_tparameter

com = test.cpp

numericdiff_SOURCES        = numericdiff.cpp
test_abs_SOURCES           = test_abs.cpp $(com)
test_arithmetic_SOURCES    = test_arithmetic.cpp $(com)
test_bitop_SOURCES         = test_bitop.cpp $(com)
test_eratosthenes_SOURCES  = test_eratosthenes.cpp $(com)
test_genmatrix_SOURCES     = test_genmatrix.cpp $(com)
test_linearalgebra_SOURCES = test_linearalgebra.cpp $(com)
test_partition_SOURCES     = test_partition.cpp $(com)
test_polynomial_SOURCES    = test_polynomial.cpp $(com)
test_prime_SOURCES         = test_prime.cpp $(com)
test_prng_SOURCES          = test_prng.cpp $(com)
test_qrnsequence_SOURCES   = test_qrnsequence.cpp $(com)
test_rule_SOURCES          = test_rule.cpp $(com)
test_shiftnet_SOURCES      = test_shiftnet.cpp $(com)
test_tparameter_SOURCES    = test_tparameter.cpp $(com)

noinst_HEADERS = test.h

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_builddir)/src
AM_CXXFLAGS = @COMPILERWARNINGS@
LDADD = ../src/libhintlib.la

clean-local:
	-rm -rf ii_files

../src/libhintlib.la:
	cd ../src && $(MAKE) libhintlib.la

