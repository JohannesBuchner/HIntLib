## Process this file with automake to produce Makefile.in

files = \
   fkmat04.cmx fkmat05.cmx fkmat06.cmx fkmat07.cmx fkmat08.cmx fkmat09.cmx \
   fkmat10.cmx fkmat11.cmx fkmat12.cmx fkmat13.cmx fkmat14.cmx fkmat15.cmx \
   fkmat16.cmx fkmat17.cmx fkmat18.cmx fkmat19.cmx fkmat20.cmx fkmat21.cmx \
   fkmat22.cmx fkmat23.cmx fkmat24.cmx fkmat25.cmx fkmat26.cmx fkmat27.cmx \
   fkmat28.cmx fkmat29.cmx fkmat30.cmx fkmat31.cmx fkmat32.cmx

# fkmat511b64.cmx fkmat63b16.cmx

filesbin = \
   fkmat04.bin fkmat05.bin fkmat06.bin fkmat07.bin fkmat08.bin fkmat09.bin \
   fkmat10.bin fkmat11.bin fkmat12.bin fkmat13.bin fkmat14.bin fkmat15.bin \
   fkmat16.bin fkmat17.bin fkmat18.bin fkmat19.bin fkmat20.bin fkmat21.bin \
   fkmat22.bin fkmat23.bin fkmat24.bin fkmat25.bin fkmat26.bin fkmat27.bin \
   fkmat28.bin fkmat29.bin fkmat30.bin fkmat31.bin fkmat32.bin

#   fkmat511b64.bin fkmat63b16.bin

pkgdata_DATA = $(filesbin)
EXTRA_DIST = $(files) convert_matrix.cpp README
CLEANFILES = $(filesbin)

noinst_PROGRAMS = convert_matrix
convert_matrix_SOURCES = convert_matrix.cpp

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_builddir)/src
AM_CXXFLAGS = @COMPILERWARNINGS@
LDADD = ../src/libhintlib.la

clean-local:
	-rm -rf ii_files

SUFFIXES = .cmx .bin

$(filesbin): convert_matrix@EXEEXT@

.cmx.bin:
	@echo Converting $< to $@ ...
	test -x ./convert_matrix@EXEEXT@ || sleep 1
	./convert_matrix@EXEEXT@ $< $@  ||  { rm -f $@ ; exit 1 ; }

../src/libhintlib.la:
	cd ../src && $(MAKE) libhintlib.la

