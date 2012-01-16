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

#ifdef __GNUG__
#pragma implementation
#endif

#include <algorithm>
#include <iostream>

#include <HIntLib/digitalnetgen.h>

#include <HIntLib/exception.h>
#include <HIntLib/modulararithmetic.h>
#include <HIntLib/precalculatedfield.h>
#include <HIntLib/bitvectorring.h>

namespace L = HIntLib;

using L::Index;

/**
 *  Digital Net
 */

L::DigitalNet::DigitalNet (unsigned base, unsigned _m)
   : m (std::min (_m,
                  unsigned(log(2.0) * double(std::numeric_limits<Index>::digits)
                                    / log(double(base)) - 1e-7))),
     size (powInt (base, m))
{}


/**
 *  The constructor initializes all member variables
 */

template<class A, typename S>
L::DigitalNetGen<A,S>::DigitalNetGen (
   const GeneratorMatrixGen<T> &_c, const A &_arith, const Hypercube &_h,
   unsigned _m, Index index, bool equi, Truncation _trunc) 
: QRNSequenceBase (_h),
  DigitalNet (_c.getBase(), _m),
  arith (_arith),
  prec (std::min ((_trunc == FULL) ? _c.getPrecision() : m,
                  digitsRepresentable (S(base()))  )),
  c (_c, getDimension(), m, prec, equi),
  x      (getDimension() * prec),
  xStart (getDimension() * prec, 0),
  trunc (_trunc),
  ss (h),
  trivialScale (1.0 / powInt(real(base()), prec))
{
   if (arith.size() != c.getBase())  throw FIXME (__FILE__, __LINE__);

   setCube (h);

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
void L::DigitalNetGen<A,S>::setCube (const Hypercube &h)
{
  ss.set (h, trunc==CENTER ? -1./base() : .0,
             powInt(real(base()), prec) - (trunc==CENTER ? 1./base() : .0));
}


/**
 *  randomize()
 */

template<class A, typename S>
void L::DigitalNetGen<A,S>::randomize (PRNG &g)
{
   for (unsigned i = 0; i < prec * getDimension(); ++i)
   {
      xStart [i] = g.equidist (int (base()));
   }
}


/**
 *  copyXtoP()
 *
 *  Copies the content of the array _x_ to the real array _point_
 */

template<class A, typename S>
void L::DigitalNetGen<A,S>::copyXtoP (real* point)
{
   for (unsigned d = 0; d < getDimension(); ++d)
   {
      S sum = 0;

      for (unsigned b = 0; b < prec; ++b)
      {
         sum = base() * sum + x[d * prec + b];
      }

      point[d] = ss[d] (sum);
   }
}

template<class A, typename S>
void L::DigitalNetGen<A,S>::copyXtoPDontScale (real* point)
{
   for (unsigned d = 0; d < getDimension(); ++d)
   {
      S sum = 0;

      for (unsigned b = 0; b < prec; ++b)
      {
         sum = base() * sum + x[d * prec + b];
      }

      point[d] = sum * trivialScale;
   }
}


/**
 *  resetX()
 */

template<class A, typename S>
void L::DigitalNetGen<A,S>::resetX (Index nn)
{
   std::copy (&xStart[0], &xStart[getDimension() * c.getPrecision()], &x[0]);
   const Index bas = base();

   for (unsigned r = 0; nn != 0; ++r)
   {
      T digit = nn % bas; nn /= bas;

      if (! arith.is0(digit))
      {
         for (unsigned d = 0; d < getDimension(); ++d)
         {
            // x [d] ^= c(d,r);

            for (unsigned b = 0; b < prec; ++b)
            {
               arith.addTo (x [d * prec + b], arith.mul(c(d,r,b), digit));
            }
         }
      }
   }
}


/**
 *  getOptimalNumber()
 */

template<class A, typename S>
Index L::DigitalNetGen<A,S>::getOptimalNumber(Index max) const
{
   if (max > getSize())  return getSize();

   return (max == 0) ? 0 : roundDownToPower (max, Index(base()));
}


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
   const Index bas = base();

   for (unsigned r = 0; n1 != n2; ++r)
   {
      typename DigitalNetGenNormal<A,S>::T digit
         = arith.sub (n2 % bas, n1 % bas); // determine change

      n1 /= bas; n2 /= bas;

      for (unsigned d = 0; d != getDimension(); ++d)
      {
         // x [d] ^= c(d,r);

         for (unsigned b = 0; b != prec; ++b)
         {
            arith.addTo (x [d * prec + b], arith.mul(c(d,r,b), digit));
         }
      }
   }
} 


/**
 *  Gray :: updateX()
 *
 *  The next method updates the array x depending on n and the vectors in c
 *  Once x is updated, these values are converted to real numbers and returned.
 */

template<class A, typename S>
void HIntLib::DigitalNetGenGray<A,S>::updateX ()
{
   const Index bas = base();
   Index n1 = n;
   Index n2 = ++n;

   int r = -1;
   Index rem1, rem2;

   do
   {
      rem1 = n1 % bas; n1 /= bas;
      rem2 = n2 % bas; n2 /= bas;
      ++r;
   }
   while (n1 != n2);

   const typename DigitalNetGenGray<A,S>::T digit = arith.sub (rem2, rem1);
   
   for (unsigned d = 0; d != getDimension(); ++d)
   {
      // x [d] ^= c(d,r);

      for (unsigned b = 0; b != prec; ++b)
      {
         arith.addTo (x [d * prec + b], arith.mul(c(d,r,b), digit));
      }
   }
} 


/**
 *  Cyclic Gray :: updateX()
 *
 *  The next method updates the array x depending on n and the vectors in c
 *  Once x is updated, these values are converted to real numbers and returned.
 */

template<class A, typename S>
void HIntLib::DigitalNetGenCyclicGray<A,S>::updateX ()
{
   const Index bas = base();
   Index n1 = n;
   Index n2 = ++n;

   unsigned r = 0;

   // Determine digit that has to be incremented 

   while ((n1 /= bas) != (n2 /= bas))  ++r;

   for (unsigned d = 0; d != getDimension(); ++d)
   {
      // x [d] ^= c(d,r);

      for (unsigned b = 0; b != prec; ++b)
      {
         arith.addTo (x [d * prec + b], c(d,r,b));
      }
   }
} 


#if 0

/**
 * Digital Net QMC Routines
 */

//  Contructor

L::DigitalNet2QMCRoutinesBase::DigitalNet2QMCRoutinesBase (
   const GeneratorMatrix2<Index>* _dm, bool e, DigitalNet::Truncation t,
   Index i, bool f)
: dm(_dm), m(-1), equi(e), trunc(t), index (i), fast(f) {}

//  get Optimal Number()

Index L::DigitalNet2QMCRoutinesBase::getOptimalNumber
   (Index max, const Hypercube &)
{
   m = ms1(max);

   if (m == -1)  return 0;
   else return Index(1) << m;
}

// checkSize()

void L::DigitalNet2QMCRoutinesBase::checkSize (
   const DigitalNet2<Index> &net, Index end)
{
   if (end > net.getSize()) throw NumbersExhausted();
}

// integrate()

void L::DigitalNet2QMCRoutines<HIntLib::real>::integrate (
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

#define HINTLIB_INSTANTIATE(X) \
   template class L::DigitalNetGen<X, Index>; \
   template void L::DigitalNetGenNormal<X, Index>::updateX (); \
   template void L::DigitalNetGenGray<X, Index>::updateX (); \
   template void L::DigitalNetGenCyclicGray<X, Index>::updateX ();

// HINTLIB_INSTANTIATE (L::ModularIntegerField<unsigned char>);
HINTLIB_INSTANTIATE (L::PrecalculatedField<unsigned char>);
#undef HINTLIB_INSTANTIATE

#if 0
template class L::DigitalNetGen<L::BitVectorRing<Index>, Index>;
template class L::DigitalNetGenNormal<L::BitVectorRing<Index>, Index>;
template class L::DigitalNetGenGray<L::BitVectorRing<Index>, Index>;

template class L::DigitalNetGen<L::BitVectorRing<unsigned char>, Index>;
template class L::DigitalNetGenNormal<L::BitVectorRing<unsigned char>, Index>;
template class L::DigitalNetGenGray<L::BitVectorRing<unsigned char>, Index>;
#endif

