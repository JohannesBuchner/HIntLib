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

#include <algorithm>

#include <HIntLib/digitalnetgen.h>

#include <HIntLib/exception.h>
#include <HIntLib/hlalgorithm.h>

namespace HIntLib
{

/**
 *  The constructor initializes all member variables
 */

template<class A, typename S>
DigitalNetGen<A,S>::DigitalNetGen (
   const A &_arith, const GeneratorMatrix &_c, const Hypercube &_h,
   unsigned _m, Index index, bool equi, Truncation _trunc) 
: QRNSequenceBase (_h),
  DigitalNet (_c.getBase(), _m),
  arith (_arith),
  scalArith (arith.getScalarAlgebra()),
  base (_c.getBase()),
  totalPrec (min3 (
     (_trunc == FULL) ? _c.getTotalPrec() : m,  // what we want
     digitsRepresentable (S(scalArith.size())), // what we can get for this S
     unsigned (HINTLIB_MN ceil (
           HINTLIB_MN log(2.0) / HINTLIB_MN log(double(scalArith.size()))
               * double(std::numeric_limits<real>::digits - 1)))
                                                // what can be stored in real
     )),
  c (_c, GMCopy().dim(getDimension()).m(m).totalPrec(totalPrec)
                 .equi(equi).vec(arith.dimension())),
  prec   (c.getPrec()),
  vecBase(c.getVecBase()),
  x      (getDimension() * prec),
  xStart (getDimension() * prec, 0),
  trunc (_trunc),
  ss (h),
  trivialScale (1.0 / HINTLIB_MN pow(real(base), int (totalPrec)))
{
   if (arith.size() != c.getVecBase())  throw FIXME (__FILE__, __LINE__);
   if (arith.dimension() != c.getVectorization())
         throw FIXME (__FILE__, __LINE__);

   setCube (h);

#if 0
   c.dump(cerr);
   c.vectorDump(cerr);
   cerr << " totalPrec=" << totalPrec
        << " prec=" << prec
        << " vecBase=" << vecBase
        << " base=" << int (base)
        << " h=" << h
        << endl
        << " arith.size()=" << arith.size()
        << " arith.dimension()=" << arith.dimension()
        << endl
        << " scalArith.size()=" << scalArith.size()
        << endl
        << " trivalScale=" << trivialScale << " " << 1/trivialScale
        << endl;
#endif

   // Initialize xStart vector

   if (index > 0)  throw InternalError (__FILE__, __LINE__);   // FIXME

#if 0
   const int dd = equi;

   unsigned i = m;
   while (index)
   {
      if (index & 1)
      {
         for (unsigned d = dd; d < c.getDimension(); ++d)
         {
            xStart [d] ^= c(d-dd,i);
         }
      }

      index /= 2;
      ++i;
   }
#endif
}


/**
 *  setCube()
 *
 *  Reinitializes the ShiftScale for a new destination cube
 */

template<class A, typename S>
void DigitalNetGen<A,S>::setCube (const Hypercube &h)
{
   real shift = trunc==CENTER ? -1./base : .0;
   ss.set (h, shift, HINTLIB_MN pow(real(base), int (totalPrec)) + shift);
}


/**
 *  randomize()
 */

template<class A, typename S>
void DigitalNetGen<A,S>::randomize (PRNG &g)
{
   for (unsigned i = 0; i < prec * getDimension(); ++i)
   {
      xStart [i] = g.equidist (int (vecBase));
   }
}


/**
 *  copyXtoP()
 *
 *  Copies the content of the array _x_ to the real array _point_
 */

template<class A, typename S>
void DigitalNetGen<A,S>::copyXtoP (real* point)
{
   for (unsigned d = 0; d < getDimension(); ++d)
   {
      S sum = 0;

      for (unsigned b = 0; b < prec; ++b)
      {
         sum = vecBase * sum + x[d * prec + b];
      }

      point[d] = ss[d] (sum);
   }
}

template<class A, typename S>
void DigitalNetGen<A,S>::copyXtoPDontScale (real* point)
{
   for (unsigned d = 0; d < getDimension(); ++d)
   {
      S sum = 0;

      for (unsigned b = 0; b < prec; ++b)
      {
         sum = vecBase * sum + x[d * prec + b];
      }

      point[d] = sum * trivialScale;
   }
}


/**
 *  resetX()
 */

template<class A, typename S>
void DigitalNetGen<A,S>::resetX (Index nn)
{
   std::copy (&xStart[0], &xStart[getDimension() * prec], &x[0]);
   const Index bas = base;
   const unsigned p = prec;

   for (unsigned r = 0; nn != 0; ++r)
   {
      typename A::scalar_type digit = nn % bas; nn /= bas;

      if (! scalArith.is0(digit))
      {
         for (unsigned d = 0; d < getDimension(); ++d)
         {
            // x [d] ^= c(d,r);

            for (unsigned b = 0; b < p; ++b)
            {
               arith.addTo (x [d * p + b], arith.mul(c(d,r,b), digit));
            }
         }
      }
   }
}


/**
 *  getOptimalNumber()
 */

template<class A, typename S>
Index DigitalNetGen<A,S>::getOptimalNumber(Index max) const
{
   if (max > getSize())  return getSize();

   return (max == 0) ? 0 : roundDownToPower (max, Index(base));
}


/********************  Normal  ***********************************************/

/**
 *  Normal :: updateX()
 *
 *  The next method updates the array x depending on n and the vectors in c
 *  Once x is updated, these values are converted to real numbers and returned.
 */

template<class A, typename S>
void HIntLib::DigitalNetGenNormal<A,S>::updateX ()
{
   Index n1 = n;
   Index n2 = ++n;
   const Index bas = base;
   const unsigned p = prec;

   for (unsigned r = 0; n1 != n2; ++r)
   {
      typename A::scalar_type digit = scalArith.sub (n2 % bas, n1 % bas);

      n1 /= bas; n2 /= bas;

      for (unsigned d = 0; d != getDimension(); ++d)
      {
         // x [d] ^= c(d,r);

         for (unsigned b = 0; b != p; ++b)
         {
            arith.addTo (x [d * p + b], arith.mul(c(d,r,b), digit));
         }
      }
   }
} 


/*******************  Gray  **************************************************/


/**
 *  Gray :: updateX()
 *
 *  The next method updates the array x depending on n and the vectors in c
 *  Once x is updated, these values are converted to real numbers and returned.
 */

template<class A, typename S>
void HIntLib::DigitalNetGenGray<A,S>::updateX ()
{
   const Index bas = base;
   const unsigned p = prec;
   Index n1 = n;
   Index n2 = ++n;

   int r = -1;
   typename A::scalar_type rem1, rem2;

   do
   {
      rem1 = n1 % bas; n1 /= bas;
      rem2 = n2 % bas; n2 /= bas;
      ++r;
   }
   while (n1 != n2);

   const typename A::scalar_type digit = scalArith.sub (rem2, rem1);
   
   for (unsigned d = 0; d != getDimension(); ++d)
   {
      // x [d] ^= c(d,r);

      for (unsigned b = 0; b != p; ++b)
      {
         arith.addTo (x [d * p + b], arith.mul(c(d,r,b), digit));
      }
   }
} 


/********************  Cyclic Gray  ******************************************/


/**
 *  Cyclic Gray :: updateX()
 *
 *  The next method updates the array x depending on n and the vectors in c
 *  Once x is updated, these values are converted to real numbers and returned.
 */

template<class A, typename S>
void HIntLib::DigitalNetGenCyclicGray<A,S>::updateX ()
{
   const Index bas = base;
   const unsigned p = prec;
   Index n1 = n;
   Index n2 = ++n;

   unsigned r = 0;

   // Determine digit that has to be incremented 

   while ((n1 /= bas) != (n2 /= bas))  ++r;

   for (unsigned d = 0; d != getDimension(); ++d)
   {
      // x [d] ^= c(d,r);

      for (unsigned b = 0; b != p; ++b)
      {
         arith.addTo (x [d * p + b], c(d,r,b));
      }
   }
} 

#define HINTLIB_INSTANTIATE_DIGITALNETGEN(X) \
   template class DigitalNetGen<X,Index>; \
   template void DigitalNetGenNormal<X,Index>::updateX (); \
   template void DigitalNetGenGray<X,Index>::updateX (); \
   template void DigitalNetGenCyclicGray<X,Index>::updateX ();


#if 0

/**
 * Digital Net QMC Routines
 */

//  Contructor

DigitalNet2QMCRoutinesBase::DigitalNet2QMCRoutinesBase (
   const GeneratorMatrix2<Index>* _dm, bool e, DigitalNet::Truncation t,
   Index i, bool f)
: dm(_dm), m(-1), equi(e), trunc(t), index (i), fast(f) {}

//  get Optimal Number()

Index DigitalNet2QMCRoutinesBase::getOptimalNumber
   (Index max, const Hypercube &)
{
   m = ms1(max);

   if (m == -1)  return 0;
   else return Index(1) << m;
}

// checkSize()

void DigitalNet2QMCRoutinesBase::checkSize (
   const DigitalNet2<Index> &net, Index end)
{
   if (end > net.getSize()) throw NumbersExhausted();
}

// integrate()

void DigitalNet2QMCRoutines<HIntLib::real>::integrate (
   Function &f, const Hypercube &h, Index begin, Index end,
   Statistic<real,real,Index>& stat)
{
   if (m < 0)  return;

   if (fast)
   {
      DigitalNet2Gray<Index> net (*dm, h, m, index, equi, trunc);
      checkSize(net, end);
      qmcIntegration (net, f, begin, end, stat);
   }
   else
   {
      DigitalNet2Normal<Index> net (*dm, h, m, index, equi, trunc);
      checkSize(net, end);
      qmcIntegration (net, f, begin, end, stat);
   }
}

#endif
}
