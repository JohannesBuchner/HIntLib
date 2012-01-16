## Process this file with automake to produce Makefile.in

SUBDIRS = @fortran_check_dirs@

TESTS = bitop polynomial eratosthenes prime arithmetic rule genmatrix @diff_check_progs@ @fortran_check_progs@

EXTRA_DIST = $(TESTS)

check_PROGRAMS = test_bitop test_polynomial test_eratosthenes test_prime test_arithmetic test_rule test_genmatrix test_qrnsequence numericdiff

com = test.cpp

test_bitop_SOURCES = test_bitop.cpp $(com)
test_rule_SOURCES  = test_rule.cpp $(com)
test_genmatrix_SOURCES = test_genmatrix.cpp $(com)
test_polynomial_SOURCES = test_polynomial.cpp $(com)
test_eratosthenes_SOURCES = test_eratosthenes.cpp $(com)
test_prime_SOURCES = test_prime.cpp $(com)
test_arithmetic_SOURCES = test_arithmetic.cpp $(com)
test_qrnsequence_SOURCES = test_qrnsequence.cpp $(com)
numericdiff_SOURCES = numericdiff.cpp

noinst_HEADERS = test.h

INCLUDES = -I$(top_srcdir)/src -I$(top_builddir)/src @GNUFLAGS@
# LDFLAGS= @SGIFLAGS@
LDADD = ../src/libhintlib.la -lm

clean-local:
	rm -rf ii_files


