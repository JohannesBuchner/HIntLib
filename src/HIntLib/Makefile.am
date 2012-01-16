## Process this file with automake to produce Makefile.in

pkginclude_HEADERS = \
   adaptintegrator.h \
   adaptintegratorlocallist.h \
   adaptintegratorms.h \
   algebra.h \
   array.h \
   bitop.h \
   buffer.h \
   builtinprng.h \
   compruleintegrator.h \
   counter.h \
   cubaturerule.h \
   cubepartitioner.h \
   defaultcubaturerulefactory.h \
   defaults.h \
   digitalnet2.h \
   digitalnetgen.h digitalnetgen.tcc \
   distribution.h \
   eratosthenes.h \
   esterr.h \
   exception.h \
   exception_MPI.h \
   factorring.h \
   fallback_limits.h \
   fallback_sstream.h \
   faure.h \
   fourthdiff.h \
   galoisfield.h \
   gcd.h gcd.tcc \
   generatormatrix.h \
   generatormatrix2.h \
   generatormatrix2row.h \
   generatormatrixgen.h \
   generatormatrixvec.h \
   halton.h \
   hypercube.h \
   integerring.h \
   integrand.h \
   integrator.h \
   kahanadd.h \
   lcg_pow2.h \
   lcg_prime.h \
   linearalgebra.h linearalgebra.tcc \
   lookupfield.h \
   make.h \
   mcintegration.h \
   mcintegrator.h \
   mcpointset.h \
   mersennetwister.h \
   minmaxfinder.h \
   miser.h \
   myalgorithm.h \
   mymath.h \
   mympi.h \
   niederreitermatrix.h \
   niederreitermatrixgen.h \
   old_bitop.h \
   onedimvectorspace.h \
   orbitrule.h \
   pointset.h \
   polynomial.h polynomial.tcc \
   polynomial2.h \
   prime.h \
   prime_generated.h \
   pseudoembeddedrule.h \
   qmcintegration.h \
   qmcintegrator.h \
   qmcroutinesimp.h \
   qrnsequence.h \
   qrnsequencebase.h \
   quotientfield.h \
   region.h \
   regioncollection.h \
   rule1midpoint.h \
   rule1trapezoidal.h \
   rule2ionescu.h \
   rule2simplex.h \
   rule2thacher.h \
   rule3cross.h \
   rule3ewing.h \
   rule3gauss.h \
   rule3octahedron.h \
   rule3simpson.h \
   rule3tyler.h \
   rule5gauss.h \
   rule5hammer.h \
   rule5stroud.h \
   rule75genzmalik.h \
   rule7phillips.h \
   rule9stenger.h \
   rulebasedintegrator.h \
   ruleintegrator.h \
   shiftnet.h \
   shiftscale.h \
   sobolmatrix.h \
   staticloadbalancer.h \
   statistic.h \
   testintegrand.h \
   tparameter.h \
   vegas.h \
   virtualprng.h

CLEANFILES = defaults.h prime_generated.h
AM_CXXFLAGS = @COMPILERWARNINGS@

noinst_HEADERS = _defaults1.h _defaults2.h

noinst_PROGRAMS = create_prime
create_prime_SOURCES = create_prime.cpp
create_prime.obj create_prime.o create_prime.lo: defaults.h

dist-hook:
	cd $(distdir) && rm -f $(CLEANFILES)

# defaults.h

# Is created by concatenating _*.h files, ../config.h, and create.h produced
# by libtest

defaults.h: $(srcdir)/_defaults1.h $(srcdir)/_defaults2.h \
            $(top_builddir)/libtest/result.h $(top_builddir)/config.h
	@echo Creating header file \"defaults.h\"
	cat $(srcdir)/_defaults1.h > defaults.h
	sed -e 's/#define [ ]*\([A-Z_]\)/#define HINTLIB_\1/' \
            -e 's/#undef [ ]*\([A-Z_]\)/#undef HINTLIB_\1/' \
            $(top_builddir)/config.h >> defaults.h
	cat $(top_builddir)/libtest/result.h \
            $(srcdir)/_defaults2.h >> defaults.h

$(top_builddir)/libtest/result.h:
	cd $(top_builddir)/libtest && $(MAKE) result.h

# prime_generated.h

prime_generated.h: create_prime@EXEEXT@
	@echo Creating header file \"prime_generated.h\"
	./create_prime@EXEEXT@ h > prime_generated.h \
           || { rm -f prime_generated.h ; exit 1; }


