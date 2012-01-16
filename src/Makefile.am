## Process this file with automake to produce Makefile.in

SUBDIRS = HIntLib

ALLSRC= \
   bitop.cpp \
   old_bitop.cpp \
   prime.cpp \
   prime_generated.cpp \
   modulararithmetic.cpp \
   polynomial.cpp \
   polynomial2.cpp \
   integerring.cpp \
   galoisfield.cpp \
   precalculatedgaloisfield.cpp \
   precalculatedfield.cpp \
   gcd.cpp \
   mymath.cpp \
   cubepartitioner.cpp \
   kahanadd.cpp \
   testfunction.cpp \
   qrnsequencebase.cpp \
   shiftscale.cpp \
   digitalnet2.cpp \
   digitalnetgen.cpp \
   sobolmatrix.cpp \
   faure.cpp \
   niederreitermatrixgen.cpp \
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
   vegasintegrator.cpp \
   miserintegrator.cpp \
   mersennetwister.cpp \
   builtinprng.cpp \
   virtualprng.cpp \
   pointset.cpp \
   mcpointset.cpp \
   generatormatrix.cpp \
   loadlibseq.cpp \
   tparameter.cpp \
   stl.cpp \
   exception.cpp \
   make.cpp makecubaturerulefactory.cpp makeqrnsequence.cpp \
   output.cpp

ALLSRCSEQ=$(ALLSRC) \
   hypercube.cpp region.cpp

ALLSRCMPI=$(ALLSRC) \
   mympi.cpp \
   region_MPI.cpp \
   hypercube_MPI.cpp \
   exception_MPI.cpp \
   staticloadbalancer.cpp \
   buffer.cpp \
   adaptintegratorms.cpp adaptintegratorlocallist.cpp \
   mcintegrator_MPI.cpp qmcintegrator_MPI.cpp

lib_LTLIBRARIES = libhintlib.la @LIB_HINTLIB_MPI@
EXTRA_LTLIBRARIES = libhintlibmpi.la
libhintlib_la_SOURCES = $(ALLSRCSEQ)
libhintlibmpi_la_SOURCES = $(ALLSRCMPI)
libhintlibmpi_la_LDFLAGS = -rpath $(libdir)

noinst_PROGRAMS = create_sobol create_niederreiter create_bitop
GNS= bitop.cpp exception.cpp hypercube.cpp generatormatrix.cpp prime.cpp prime_generated.cpp mymath.cpp polynomial2.cpp
create_sobol_SOURCES = create_sobol.cpp $(GNS)
create_niederreiter_SOURCES = create_niederreiter.cpp $(GNS) integerring.cpp niederreitermatrixgen.cpp modulararithmetic.cpp polynomial.cpp precalculatedfield.cpp precalculatedgaloisfield.cpp galoisfield.cpp gcd.cpp output.cpp
create_bitop_SOURCES = create_bitop.cpp old_bitop.cpp

INCLUDES = @GNUFLAGS@
# LDFLAGS=
LDADD = -lm

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
	cd HIntLib && $(MAKE) create_prime

