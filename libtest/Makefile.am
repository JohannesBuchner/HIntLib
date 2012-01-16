## Process this file with automake to produce Makefile.in

EXTRA_DIST = create

noinst_HEADERS = include.h

lib_LTLIBRARIES = libhintlibtest.la
libhintlibtest_la_SOURCES = func.cpp lib.cpp
libhintlibtest_la_LDFLAGS = @NOUNDEFINED@

noinst_PROGRAMS = main
main_SOURCES = main.cpp

AM_CXXFLAGS = @COMPILERWARNINGS@
LDADD = libhintlibtest.la

result.h: main@EXEEXT@ create
	$(srcdir)/create

clean-local:
	rm -rf ii_files result.h

