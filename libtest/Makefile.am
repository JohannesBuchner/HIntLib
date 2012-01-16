## Process this file with automake to produce Makefile.in

EXTRA_DIST = create

noinst_HEADERS = include.h

lib_LTLIBRARIES = libhintlibtest.la
libhintlibtest_la_SOURCES = func.cpp lib.cpp

noinst_PROGRAMS = main
main_SOURCES = main.cpp

INCLUDES =  @GNUFLAGS@
# LDFLAGS= @SGIFLAGS@
LDADD = libhintlibtest.la

result.h: main create
	$(srcdir)/create

clean-local:
	rm -rf ii_files .libs .deps result.h

