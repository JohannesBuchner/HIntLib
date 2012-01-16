## Process this file with automake to produce Makefile.in

SUBDIRS = HIntLib

noinst_LTLIBRARIES = libcreatelookupfield.la libcreatematrix.la
lib_LTLIBRARIES = libhintlib.la @LIB_HINTLIB_MPI@
EXTRA_LTLIBRARIES = libhintlibmpi.la

libcreatelookupfield_la_LIBADD = @STL_OBJECT@
libcreatelookupfield_la_DEPENDENCIES = @STL_OBJECT@
libcreatelookupfield_la_SOURCES = \
   bitop.cpp \
   exception.cpp \
   factorring_lf.cpp \
   galoisfield.cpp \
   linearalgebra_lf.cpp \
   modulararithmetic.cpp \
   hlmath.cpp \
   output.cpp \
   polynomial.cpp \
   polynomial_lf.cpp \
   prime.cpp \
   prime_generated.cpp
EXTRA_libcreatelookupfield_la_SOURCES = stl.cpp

libcreatematrix_la_LDFLAGS = -lm
libcreatematrix_la_LIBADD = libcreatelookupfield.la
libcreatematrix_la_SOURCES = \
   generatormatrix.cpp \
   generatormatrix2.cpp \
   generatormatrixgen.cpp \
   linearalgebra.cpp \
   lookupfield.cpp \
   lookupfield_generated.cpp \
   lookupgaloisfield.cpp \
   niederreitermatrixgen.cpp \
   polynomial_gm.cpp \
   polynomial2.cpp

ALLSRC= \
   adaptintegrator.cpp \
   builtinprng.cpp \
   compruleintegrator.cpp \
   counter.cpp \
   cubepartitioner.cpp \
   digitalnet2.cpp \
   digitalnetgen.cpp \
   factorring_x.cpp \
   faure.cpp \
   generatormatrix2row.cpp \
   generatormatrixvec.cpp \
   halton.cpp \
   hypercube.cpp \
   integerring.cpp \
   integrator.cpp \
   kahanadd.cpp \
   loadlibseq.cpp \
   make.cpp \
   makecubaturerulefactory.cpp \
   makeqrnsequence.cpp \
   mcintegrator.cpp \
   mcpointset.cpp \
   mersennetwister.cpp \
   miser.cpp \
   niederreitermatrix.cpp \
   orbitrule.cpp \
   pointset.cpp \
   polynomial_field.cpp \
   polynomial_ring.cpp \
   pseudoembeddedrule.cpp \
   qmcintegrator.cpp \
   qrnsequencebase.cpp \
   quotientfield.cpp \
   realfield.cpp \
   regioncollection.cpp \
   rule1midpoint.cpp \
   rule1trapezoidal.cpp \
   rule2ionescu.cpp \
   rule2simplex.cpp \
   rule2thacher.cpp \
   rule3cross.cpp \
   rule3ewing.cpp \
   rule3gauss.cpp \
   rule3octahedron.cpp \
   rule3simpson.cpp \
   rule3tyler.cpp \
   rule5gauss.cpp \
   rule5hammer.cpp \
   rule5stroud.cpp \
   rule75genzmalik.cpp \
   rule7phillips.cpp \
   rule9stenger.cpp \
   rulebasedintegrator.cpp \
   ruleintegrator.cpp \
   shiftnet.cpp \
   shiftscale.cpp \
   sobolmatrix.cpp \
   testintegrand.cpp \
   tparameter.cpp \
   vegas.cpp \
   virtualprng.cpp

ALLSRCSEQ=$(ALLSRC) \
   region.cpp \
   staticloadbalancer.cpp

ALLSRCMPI=$(ALLSRC) \
   adaptintegratorlocallist.cpp \
   adaptintegratorms.cpp \
   buffer.cpp \
   exception_MPI.cpp \
   hypercube_MPI.cpp \
   mcintegrator_MPI.cpp \
   hlmpi.cpp \
   qmcintegrator_MPI.cpp \
   region_MPI.cpp \
   staticloadbalancer_MPI.cpp

libhintlib_la_SOURCES = $(ALLSRCSEQ)
libhintlib_la_LDFLAGS = @NOUNDEFINED@
libhintlib_la_LIBADD = libcreatematrix.la

libhintlibmpi_la_SOURCES = $(ALLSRCMPI)
libhintlibmpi_la_LDFLAGS = @NOUNDEFINED@ -rpath $(libdir)
libhintlibmpi_la_LIBADD = libcreatematrix.la

noinst_PROGRAMS = create_sobol create_niederreiter create_bitop create_lookupfield
create_bitop_SOURCES = create_bitop.cpp
create_lookupfield_SOURCES = create_lookupfield.cpp
create_lookupfield_LDADD = libcreatelookupfield.la
create_sobol_SOURCES = create_sobol.cpp
create_sobol_LDADD = libcreatematrix.la libcreatelookupfield.la
create_niederreiter_SOURCES = create_niederreiter.cpp
create_niederreiter_LDADD = libcreatematrix.la libcreatelookupfield.la

AM_CXXFLAGS = @COMPILERWARNINGS@

CLEANFILES = prime_generated.cpp sobolmatrix.cpp niederreitermatrix.cpp bitop.cpp lookupfield_generated.cpp

clean-local:
	rm -rf ii_files

dist-hook:
	cd $(distdir) && rm -f $(CLEANFILES)


# Produce program-generated Source files

bitop.cpp: create_bitop@EXEEXT@
	@echo Creating source file \"bitop.cpp\"
	./create_bitop@EXEEXT@ > bitop.cpp \
           || { rm -f bitop.cpp ; exit 1 ; }

niederreitermatrix.cpp: create_niederreiter@EXEEXT@
	@echo Creating source file \"niederreitermatrix.cpp\"
	./create_niederreiter@EXEEXT@ > niederreitermatrix.cpp \
            || { rm -f niederreitermatrix.cpp ; exit 1 ; }

sobolmatrix.cpp: create_sobol@EXEEXT@
	@echo Creating source file \"sobolmatrix.cpp\"
	./create_sobol@EXEEXT@ > sobolmatrix.cpp \
           || { rm -f sobolmatrix.cpp ; exit 1 ; }

prime_generated.cpp: HIntLib/create_prime@EXEEXT@
	@echo Creating source file \"prime_generated.cpp\"
	HIntLib/create_prime@EXEEXT@ cpp > prime_generated.cpp \
           || { rm -f prime_generated.cpp ; exit 1 ; }
HIntLib/prime_generated.h:
	cd HIntLib && $(MAKE) prime_generated.h

lookupfield_generated.cpp: create_lookupfield@EXEEXT@
	@echo Creating source file \"lookupfield_generated.cpp\"
	./create_lookupfield@EXEEXT@ > lookupfield_generated.cpp \
           || { rm lookupfield_generated.cpp ; exit 1 ; }

HIntLib/create_prime@EXEEXT@:
	cd HIntLib && $(MAKE) create_prime@EXEEXT@

