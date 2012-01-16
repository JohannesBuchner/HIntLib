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
 *  Digital Net for arbitrary bases and its corresponding QMCRoutines
 */

#ifndef HINTLIB_DIGITAL_NET_GEN_H
#define HINTLIB_DIGITAL_NET_GEN_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/generatormatrixvec.h>
#include <HIntLib/qrnsequencebase.h>
#include <HIntLib/shiftscale.h>
#include <HIntLib/pointset.h>
#include <HIntLib/qmcintegration.h>

namespace HIntLib
{

/**
 *  Digital Net Gen
 *
 *  A Digital Net with arbitrary
 *
 *  A .. the Arithmetic used for all calculations (a vector space)
 *  S .. unsigned integer type that can hold base^prec-1
 */

template<class A, class S>
class DigitalNetGen: public QRNSequenceBase, public DigitalNet
{
public:
   typedef typename A::type T;

   const GeneratorMatrixVec<T> & getGeneratorMatrix() const  { return c; }
   Index getOptimalNumber(Index max) const;

   void setCube (const Hypercube &);
   void randomize (PRNG &g);

protected:

   const A arith;
   const typename A::scalar_algebra scalArith;
   const T base;
   const unsigned totalPrec;
   GeneratorMatrixVecCopy<T> c;
   const unsigned prec;
   const S vecBase;
   Array<T> x;      // current vector (size dim)
   Array<T> xStart; // Inital values for x (size dim)
   const Truncation trunc;
   ShiftScale ss;
   const real trivialScale;

   DigitalNetGen
      (const A &, const GeneratorMatrix &, const Hypercube &,
       unsigned m, Index index, bool equi, Truncation);

   void resetX          (Index n, real* p) { resetX(n); copyXtoP(p); }
   void resetXDontScale (Index n, real* p) { resetX(n); copyXtoPDontScale(p); }

   void copyXtoP          (real*);
   void copyXtoPDontScale (real*);

private:
   void resetX (Index);
};


/**
 *  Digital Net Gen Naive
 *
 *  Traverses the net in normal order - naive implementation
 */

template<class A, class S>
class DigitalNetGenNaive : public DigitalNetGen<A,S>
{
public:
   DigitalNetGenNaive
      (const A& _arith, const GeneratorMatrix &gm, const Hypercube &_h,
       unsigned m, Index i, bool equi, DigitalNet::Truncation t)
      : DigitalNetGen<A,S> (_arith, gm, _h, m, i, equi, t) {}

   DigitalNetGenNaive
      (const A &_arith, const GeneratorMatrix &gm, const Hypercube &_h)
      : DigitalNetGen<A,S> (_arith, gm, _h, gm.getM(), 0, false, FULL) {}

   void first          (real *p, Index new_n){ resetX          (n = new_n, p); }
   void firstDontScale (real *p, Index new_n){ resetXDontScale (n = new_n, p); }
   void next           (real *p)             { resetX          (++n, p); }
   void nextDontScale  (real *p)             { resetXDontScale (++n, p); }
};


/**
 *  Digital Net Gen Normal
 *
 *  Traverses the net in normal order
 */

template<class A, class S>
class DigitalNetGenNormal : public DigitalNetGen<A,S>
{
public:
   DigitalNetGenNormal
      (const A& _arith, const GeneratorMatrix &gm, const Hypercube &_h,
       unsigned m, Index i, bool equi, DigitalNet::Truncation t)
      : DigitalNetGen<A,S> (_arith, gm, _h, m, i, equi, t) {}

   DigitalNetGenNormal
      (const A& _arith, const GeneratorMatrix &gm, const Hypercube &_h)
      : DigitalNetGen<A,S> (_arith, gm, _h, gm.getM(), 0, false, FULL) {}

   void first          (real *p, Index new_n){ resetX (n = new_n, p); }
   void firstDontScale (real *p, Index new_n){ resetXDontScale (n = new_n, p); }

   void next           (real *p)  { updateX(); copyXtoP (p); }
   void nextDontScale  (real *p)  { updateX(); copyXtoPDontScale (p); }

private:
   void updateX ();
};


/**
 *  Digital Net Gen Gray
 *
 *  Traverses the net in gray-code order
 */

template<class A, class S>
class DigitalNetGenGray : public DigitalNetGen<A,S>
{
public:
   DigitalNetGenGray
      (const A& _arith, const GeneratorMatrix &gm, const Hypercube &_h,
       unsigned m, Index i, bool equi, DigitalNet::Truncation t)
      : DigitalNetGen<A,S> (_arith, gm, _h, m, i, equi, t) {}

   DigitalNetGenGray
      (const A& _arith, const GeneratorMatrix &gm, const Hypercube &_h)
      : DigitalNetGen<A,S> (_arith, gm, _h, gm.getM(), 0, false, FULL) {}

   void first          (real *p, Index new_n){ resetX (n = new_n, p); }
   void firstDontScale (real *p, Index new_n){ resetXDontScale (n = new_n, p); }

   void next           (real *p)  { updateX(); copyXtoP (p); }
   void nextDontScale  (real *p)  { updateX(); copyXtoPDontScale (p); }

private:
   void updateX ();
};


/**
 *  Digital Net Gen Cyclic Gray
 *
 *  Traverses the net in gray-code order
 *
 *  This optimized routine assumes that A is a filed with cyclic additive
 *  group (i.e.  base is prime).
 */

template<class A, class S>
class DigitalNetGenCyclicGray : public DigitalNetGen<A,S>
{
public:
   DigitalNetGenCyclicGray
      (const A& _arith, const GeneratorMatrix &gm, const Hypercube &_h,
       unsigned m, Index i, bool equi, DigitalNet::Truncation t)
      : DigitalNetGen<A,S> (_arith, gm, _h, m, i, equi, t) {}

   DigitalNetGenCyclicGray
      (const A& _arith, const GeneratorMatrix &gm, const Hypercube &_h)
      : DigitalNetGen<A,S> (_arith, gm, _h, gm.getM(), 0, false, FULL) {}

   void first          (real *p, Index new_n){ resetX (n = new_n, p); }
   void firstDontScale (real *p, Index new_n){ resetXDontScale (n = new_n, p); }

   void next           (real *p)  { updateX(); copyXtoP (p); }
   void nextDontScale  (real *p)  { updateX(); copyXtoPDontScale (p); }

private:
   void updateX ();
};


#if 0

/**
 *  QMC Routines for Digital Net
 */

// Non-template base

class DigitalNet2QMCRoutinesBase : public QMCRoutines
{
protected:
   const GeneratorMatrix2<Index>* dm;
   int m;
   bool equi;
   DigitalNet::Truncation trunc;
   Index index;
   bool fast;

   void checkSize (const DigitalNet2<Index> &, Index);

public:
   DigitalNet2QMCRoutinesBase (
      const GeneratorMatrix2<Index>* s, bool equi, DigitalNet::Truncation,
      Index index, bool fast);

   Index getOptimalNumber (Index max, const Hypercube &);
};

// Template

template<class Sum>
class DigitalNet2QMCRoutines: public DigitalNet2QMCRoutinesBase
{
private:
   typedef StatisticVar<real, Sum, Index> S;

public:
   DigitalNet2QMCRoutines (
      const GeneratorMatrix2<Index>* _dm,
      bool e = true, DigitalNet::Truncation t = DigitalNet2Base::TRUNCATE,
      Index i = 0, bool fast = true)
   : DigitalNet2QMCRoutinesBase(_dm, e, t, i, fast) {}

   virtual void integrate (
      Integrand &f, const Hypercube &h, Index begin, Index end,
      Statistic<real,real,Index>& stat);
};

// Specializatioin for Sum=real

template<>
class DigitalNet2QMCRoutines<real>: public DigitalNet2QMCRoutinesBase
{
public:
   DigitalNet2QMCRoutines (
      const GeneratorMatrix2<Index>* _dm,
      bool e = true, DigitalNet::Truncation t = DigitalNet2Base::TRUNCATE,
      Index i = 0, bool fast = true)
   : DigitalNet2QMCRoutinesBase(_dm, e, t, i, fast) {}

   virtual void integrate (
      Integrand &f, const Hypercube &h, Index begin, Index end,
      Statistic<real,real,Index>& stat);
};

#endif

} // namespace HIntLib




/*****  Implementation  **************/


#if 0

/**
 *  Gray :: next()
 *
 *  The next method updates the array x depending on n and the vectors in v
 *  Once x is updated, these values are converted to real numbers and returned.
 */

template<class T>
inline
void HIntLib::DigitalNet2Gray<T>::next (real* point)
{
   // Determine which bit changed in gray-code notation of n
 
   unsigned r = ls0 (n++);
 
   // Update x and points

   const ShiftScale1* ssp = &ss[0];
   T *xp = x;
   const T *vp = c(r);
 
   for (unsigned i = getDimension(); i; --i)
   {
      *point++ = (*ssp++) (*xp++ ^= *vp++);
   }
} 
#endif


#if 0

/**
 *  integrate()
 */

template<class Sum>
void HIntLib::DigitalNet2QMCRoutines<Sum>::integrate (
   Integrand &f, const Hypercube &h, Index begin, Index end,
   Statistic<real,real,Index>& stat)
{
   if (m < 0)  return;

   if (fast)
   {
      DigitalNet2Gray<Index> net (*dm, h, m, index, equi, trunc);

      checkSize(net, end);

      S s;
      qmcIntegration (net, f, begin, end, s);
      stat = s;
   }
   else
   {
      DigitalNet2Normal<Index> net (*dm, h, m, index, equi, trunc);

      checkSize(net, end);

      S s;
      qmcIntegration (net, f, begin, end, s);
      stat = s;
   }
}
#endif

#endif

