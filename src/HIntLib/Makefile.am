## Process this file with automake to produce Makefile.in

pkginclude_HEADERS = \
   miser.h \
   adaptintegrator.h             modulararithmetic.h \
   adaptintegratorlocallist.h    myalgorithm.h \
   adaptintegratorms.h           mymath.h gcd.h \
   array.h                       mympi.h \
   bitop.h                       niederreitermatrix.h niederreitermatrixgen.h \
   faure.h \
   buffer.h                      old_bitop.h \
   builtinprng.h                 orbitrule.h \
   compruleintegrator.h          polynomial2.h polynomial.h \
   galoisfield.h \
   lookupfield.h \
                                 prime.h integerring.h \
   cubaturerule.h                pseudoembeddedrule.h \
   cubepartitioner.h             qmcintegration.h \
   defaultcubaturerulefactory.h  qmcintegrator.h \
   qmcroutinesimp.h \
   digitalnet2.h                 qrnsequencebase.h \
   digitalnetgen.h digitalnetgen.tcc  qrnsequence.h \
   distribution.h                regioncollection.h \
   embeddedrule.h                region.h \
   eratosthenes.h                rule1midpoint.h \
   esterr.h                      rule1trapezoidal.h \
   exception.h                   rule2simplex.h rule2thacher.h rule2ionescu.h \
   exception_MPI.h               rule3gauss.h rule3cross.h \
   fourthdiff.h                  rule3octahedron.h rule3tyler.h rule3ewing.h \
   integrand.h                    rule5hammer.h \
   generatormatrix.h generatormatrix2.h \
   generatormatrixgen.h generatormatrixvec.h \
   rule5stroud.h \
   halton.h                      rule75genzmalik.h \
         rule7phillips.h \
                rule9stenger.h \
   hypercube.h                   rulebasedintegrator.h \
   integrator.h                  ruleintegrator.h \
   kahanadd.h                    shiftscale.h \
   lcg_pow2.h                    sobolmatrix.h \
   lcg_prime.h                   staticloadbalancer.h \
   make.h                        statistic.h \
   mcintegration.h               \
   mcintegrator.h                \
   testintegrand.h \
   mcpointset.h               vegas.h \
   mersennetwister.h             virtualprng.h \
   minmaxfinder.h \
   prime_generated.h \
   defaults.h  pointset.h algebra.h \
   fallback_limits.h fallback_sstream.h \
   onedimvectorspace.h tparameter.h

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


