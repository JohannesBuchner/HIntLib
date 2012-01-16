## Process this file with automake to produce Makefile.in

pkginclude_HEADERS = \
   adaptintegrator.h \
   adaptintegratorlocallist.h \
   adaptintegratorms.h \
   algebra.h \
   array.h \
   bitop.h bitop.tcc \
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
   discrepancy.h \
   distribution.h \
   eratosthenes.h \
   esterr.h \
   exception.h \
   exception_MPI.h \
   factorring.h factorring.tcc \
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
   generatormatrixgenrow.h \
   generatormatrixvec.h \
   generatormatrixvirtual.h \
   gf2.h \
   gf2vectorspace.h gf2vectorspace.tcc \
   halton.h \
   hlalgorithm.h \
   hlmath.h \
   hlmpi.h \
   hypercube.h \
   integerring.h \
   integrand.h \
   integrator.h \
   kahanadd.h \
   lcg_pow2.h \
   lcg_prime.h \
   linearalgebra.h \
   linearalgebragen.h linearalgebragen.tcc \
   linearalgebra2.h linearalgebra2.tcc \
   linereader.h \
   lookupfield.h \
   make.h \
   mcintegration.h \
   mcintegrator.h \
   mcpointset.h \
   mersennetwister.h \
   minmaxfinder.h \
   miser.h \
   modulararithmetic.h \
   niederreitermatrix.h \
   niederreitermatrixgen.h \
   old_bitop.h \
   onedimvectorspace.h \
   orbitrule.h \
   output.h \
   pointset.h \
   polynomialbase.h \
   polynomial.h polynomial.tcc polynomial_ring.tcc polynomial_field.tcc \
      polynomial_gf.tcc polynomial_rational.tcc polynomial_real.tcc \
   polynomial2base.h \
   polynomial2.h polynomial2.tcc \
   prime.h \
   prime_generated.h \
   productrule.h \
   pseudoembeddedrule.h \
   qmcintegration.h \
   qmcintegrator.h \
   qmcroutinesimp.h \
   qrnsequence.h \
   qrnsequencebase.h \
   quotientfield.h quotientfield.tcc \
   realfield.h \
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
   rule5mustardlynessblatt.h \
   rule5stroud.h \
   rule5stroud2.h \
   rule75genzmalik.h \
   rule7phillips.h \
   rule9stenger.h \
   rulebasedintegrator.h \
   rulegauss.h \
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
AM_CXXFLAGS = @COMPILERWARNINGS@ -I..

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
	test -x ./create_prime@EXEEXT@ || sleep 1
	./create_prime@EXEEXT@ h > prime_generated.h \
           || { rm -f prime_generated.h ; exit 1; }


