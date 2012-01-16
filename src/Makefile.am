## Process this file with automake to produce Makefile.in

SUBDIRS = HIntLib

ALLSRC= \
   cubepartitioner.cpp \
   kahanadd.cpp \
   testfunction.cpp \
   qrnsequencebase.cpp \
   shiftscale.cpp \
   digitalnet2.cpp \
   digitalnetgen.cpp \
   sobolmatrix.cpp \
   faure.cpp \
   niederreitermatrix.cpp \
   halton.cpp \
   orbitrule.cpp \
   pseudoembeddedrule.cpp \
   rule1midpoint.cpp \
   rule1trapezoidal.cpp \
   rule2simplex.cpp \
   rule2thacher.cpp \
   rule3octahedron.cpp \
   rule3cross.cpp \
   rule3tyler.cpp \
   rule3gauss.cpp \
   rule3ewing.cpp \
   rule5hammer.cpp \
   rule5stroud.cpp \
   rule7phillips.cpp \
   rule9stenger.cpp \
   rule75genzmalik.cpp \
   integrator.cpp \
   rulebasedintegrator.cpp \
   ruleintegrator.cpp \
   compruleintegrator.cpp \
   adaptintegrator.cpp \
   adapt.cpp \
   regioncollection.cpp \
   qmcintegrator.cpp \
   mcintegrator.cpp \
   vegas.cpp \
   miser.cpp \
   mersennetwister.cpp \
   builtinprng.cpp \
   virtualprng.cpp \
   pointset.cpp \
   mcpointset.cpp \
   loadlibseq.cpp \
   tparameter.cpp \
   make.cpp makecubaturerulefactory.cpp makeqrnsequence.cpp

ALLSRCSEQ=$(ALLSRC) \
   region.cpp

ALLSRCMPI=$(ALLSRC) \
   mympi.cpp \
   region_MPI.cpp \
   hypercube_MPI.cpp \
   exception_MPI.cpp \
   staticloadbalancer.cpp \
   buffer.cpp \
   adaptintegratorms.cpp adaptintegratorlocallist.cpp \
   mcintegrator_MPI.cpp qmcintegrator_MPI.cpp

noinst_LTLIBRARIES = libcreatematrix.la
lib_LTLIBRARIES = libhintlib.la @LIB_HINTLIB_MPI@
EXTRA_LTLIBRARIES = libhintlibmpi.la

libcreatematrix_la_SOURCES = \
   bitop.cpp \
   exception.cpp \
   hypercube.cpp \
   generatormatrix.cpp \
   prime.cpp \
   prime_generated.cpp \
   mymath.cpp \
   integerring.cpp \
   niederreitermatrixgen.cpp \
   modulararithmetic.cpp \
   polynomial.cpp \
   polynomial2.cpp \
   lookupfield.cpp \
   lookupgaloisfield.cpp \
   galoisfield.cpp \
   gcd.cpp \
   output.cpp \
   stl.cpp

libcreatematrix_la_LIBADD = -lm

libhintlib_la_SOURCES = $(ALLSRCSEQ)
libhintlib_la_LDFLAGS = @NOUNDEFINED@
libhintlib_la_LIBADD = libcreatematrix.la

libhintlibmpi_la_SOURCES = $(ALLSRCMPI)
libhintlibmpi_la_LDFLAGS = @NOUNDEFINED@ -rpath $(libdir)
libhintlibmpi_la_LIBADD = libcreatematrix.la

noinst_PROGRAMS = create_sobol create_niederreiter create_bitop
create_sobol_SOURCES = create_sobol.cpp
create_sobol_LDADD = libcreatematrix.la
create_niederreiter_SOURCES = create_niederreiter.cpp
create_niederreiter_LDADD = libcreatematrix.la
create_bitop_SOURCES = create_bitop.cpp

INCLUDES = @GNUFLAGS@

CLEANFILES = prime_generated.cpp sobolmatrix.cpp niederreitermatrix.cpp bitop.cpp

clean-local:
	rm -rf ii_files

dist-hook:
	cd $(distdir) && rm -f $(CLEANFILES)


# Produce program-generated Source files

bitop.cpp: create_bitop@EXEEXT@
	@echo Creating source file \"bitop.cpp\"
	./create_bitop@EXEEXT@ > bitop.cpp

niederreitermatrix.cpp: create_niederreiter@EXEEXT@
	@echo Creating source file \"niederreitermatrix.cpp\"
	./create_niederreiter@EXEEXT@ > niederreitermatrix.cpp

sobolmatrix.cpp: create_sobol@EXEEXT@
	@echo Creating source file \"sobolmatrix.cpp\"
	./create_sobol@EXEEXT@ > sobolmatrix.cpp

prime_generated.cpp: HIntLib/create_prime@EXEEXT@
	@echo Creating source file \"prime_generated.cpp\"
	HIntLib/create_prime@EXEEXT@ cpp > prime_generated.cpp

HIntLib/create_prime@EXEEXT@:
	cd HIntLib && $(MAKE) create_prime@EXEEXT@

