/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002  Rudolf Schürer <rudolf.schuerer@sbg.ac.at>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA
 */

/**
 *  QMCIntegrator
 *  QMCIntegratorNoLB    (when compiled with PARALLEL defined)
 *
 *  The basic Quasi Monte Carlo integration routine
 *
 *  In Parllel mode, static load balancing ist provided
 */

// Check different flags, depending on parallel mode

#if defined(PARALLEL) && !defined(QMCINTEGRATOR_MPI_H) || !defined(PARALLEL) && !defined(QMCNINTEGRATOR_H)

#include <HIntLib/integrator.h>

// Define Name macro and set flag according to parallel mode

#ifdef PARALLEL
   #define QMCINTEGRATOR_MPI_H 1
   #define NAME(x) x##StaticLB
#else
   #define QMCINTEGRATOR_H 1
   #define NAME(x) x
#endif


namespace HIntLib
{

class PointSet;
class PartitionablePointSet;

class NAME(QMCIntegrator) : public Integrator
{
private:

#if PARALLEL
   typedef PartitionablePointSet PS;
#else
   typedef PointSet PS;
#endif

public:
   NAME(QMCIntegrator) (PS* _ps) : ps(_ps) {}

   virtual
   Status integrate (
         Function &, const Hypercube &, Index maxEval,
         real reqAbsError, real reqRelError, EstErr &ee);

private:
   PS* ps;
};



#if 0
/**
 *  QMCIntegratorBase(NoLB)
 *
 *  Abstract base class for all QMCIntegrators.
 *
 *  Contains the non-timecritical part of the integration algorithm.
 *
 *  For the actual sample-loop, a virtual function integr() is called, which
 *  has to be supplied by the derived class.
 */
 
class NAME(QMCIntegratorBase) : public Integrator
{
public:

   NAME(QMCIntegratorBase) (Index s) : skip (s) {}

   virtual
   Status integrate (Function &,
                     const Hypercube &,
                     Index maxEval,
                     real reqRelError, real reqAbsError,
                     EstErr &ee);
 
protected:
 
   virtual void integr (
      Function&, const Hypercube&, Statistic<> &,
      Index n, Index begin, Index end) = 0;

   virtual Index getOptimalNumber (unsigned dim, Index) = 0;

   Index skip;
};
 
 
/**
 *  QMCIntegrator(NoLB)<>
 *
 *  Actual impelmentation for a given QRNGenerator
 *
 *  Generator as well as Summation type are template parameter, allowing
 *  inlining of the corresponding function calls in the main loop.
 */
 
template<class G, class Sum = real>
class NAME(QMCIntegrator) : public NAME(QMCIntegratorBase)
{
public:
   NAME(QMCIntegrator) (Index s = 1) : NAME(QMCIntegratorBase)(s) {}

protected:
 
   virtual void integr (
      Function &, const Hypercube&, Statistic<> &,
      Index n, Index begin, Index end);

   virtual Index getOptimalNumber (unsigned dim, Index);
};

template<class G, class Sum>
inline
Index NAME(QMCIntegrator)<G,Sum>::getOptimalNumber(unsigned dim, Index n)
{
   G g (dim);

   return g.getOptimalNumber (n + skip) - skip;
}
 
template<class G, class Sum>
inline
void NAME(QMCIntegrator)<G,Sum>::integr (
   Function &f, const Hypercube &h, Statistic<> &s,
   Index n, Index begin, Index end)
{
   G g (h);
   Statistic<real,Sum> stat;
 
   qmcIntegration (g, f, begin + skip, end + skip, stat);
 
   s = stat;
}
#endif

#if 0

/**
 *  QMCNetIntegrator(NoLB)<>
 *
 *  Actual impelmentation of a Net for a given QRNGenerator
 *
 *  Generator as well as Summation type are template parameter, allowing
 *  inlining of the corresponding function calls in the main loop.
 */
 
template<class G, class Sum = real>
class NAME(QMCNetIntegrator) : public NAME(QMCIntegratorBase)
{
public:
   NAME(QMCNetIntegrator) (Index s = 1) : NAME(QMCIntegratorBase)(s) {}

protected:
 
   virtual void integr (
      Function &, const Hypercube&, Statistic<> &,
      Index n, Index begin, Index end);

   virtual Index getOptimalNumber (unsigned dim, Index);
};
 
template<class G, class Sum>
inline
void NAME(QMCNetIntegrator)<G,Sum>::integr (
   Function &f, const Hypercube &h, Statistic<> &s,
   Index n, Index begin, Index end)
{
   QRNNet<G> g (h.getDimension(), n + skip);
   g.setIndex(begin + skip);
 
   Statistic<real,Sum> stat;
 
   qmcIntegration (g, h, f, end - begin, stat);
 
   s = stat;
}

template<class G, class Sum>
inline
Index NAME(QMCNetIntegrator)<G,Sum>::getOptimalNumber(unsigned dim, Index n)
{
   G g (dim);

   return g.getOptimalNumber (n + skip) - skip;
}

#endif

}  // namespace HIntLib

#undef NAME

#endif

