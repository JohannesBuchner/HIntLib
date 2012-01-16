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

#include <HIntLib/digitalnet2.h>

#include <HIntLib/exception.h>

namespace L = HIntLib;

using L::Index;

/*****************************************************************************/
/***         Digital Net                                                   ***/
/*****************************************************************************/

L::DigitalNet::DigitalNet (unsigned _m)
   : m (std::min (_m, unsigned (std::numeric_limits<Index>::digits) - 1)),
     size (Index(1) << m)
{}


/*****************************************************************************/
/***         Digital Net 2 <T>                                             ***/
/*****************************************************************************/

/**
 *  Constructor
 *
 *  Initialize all member variables
 */

template<class T>
L::DigitalNet2<T>::DigitalNet2 (
   const GeneratorMatrix2<T> &_c, const Hypercube &_h,
   unsigned _m, Index index, bool equi, Truncation _trunc) 
: QRNSequenceBase(_h),
  DigitalNet (_m),
  prec (std::min((_trunc == FULL) ? _c.getPrecision() : m,
                 unsigned (std::numeric_limits<real>::digits - 1))),
  c (_c, getDimension(), m, prec, equi),
  x      (getDimension()),
  xStart (getDimension(), 0),
  trunc (_trunc),
  trivialScale (powInt(.5, c.getPrecision())),
  ss (h.getDimension())
#ifdef HINTLIB_IEEE_MAGIC_WORKS
  , ssMagic (h.getDimension()),
  mode (STANDARD)
#endif
{
   // Initialize Shift Scale

   setCube (h);
   
   // Initialize xStart vector

   const int dd = equi;
   unsigned i = m;
   Index indexCopy = index;
   int shift = _c.getPrecision() - c.getPrecision();

   while (index)
   {
      if (i >= _c.getM())
      {
         throw NetIndexToHigh (indexCopy, _c.getM(), m);
      }
      
      if (index & 1)
      {
         for (unsigned d = dd; d < c.getDimension(); ++d)
         {
            xStart [d] ^= (_c(d-dd,i) >> shift);
         }
      }

      index /= 2;
      ++i;
   }

#ifdef HINTLIB_IEEE_MAGIC_WORKS
   // Prepare for direct updates, if we can

   if (c.getPrecision() <= floatBits)
   {
      mode = DIRECT;

      const int shift = floatBits - c.getPrecision();
      union { floatType f; T i; } x;
      x.f = 1.;

      for (unsigned d = 0; d < getDimension(); ++d)
      {
         xStart [d] = (xStart[d] << shift) | x.i;
      }

      c.adjustPrecision (floatBits);
   }
#endif
}


/**
 *  setCube()
 *
 *  Reinitializes the ShiftScale(s) for a new destination cube
 */

template<class T>
void L::DigitalNet2<T>::setCube (const Hypercube &h)
{
  ss.set (h, trunc == CENTER ? -.5 : .0, 
             powInt(2.0, c.getPrecision()) - (trunc == CENTER ? .5 : .0));
#ifdef HINTLIB_IEEE_MAGIC_WORKS
  ssMagic.set (h, trunc == CENTER ? -powInt(.5, c.getPrecision()+1) : .0, 
                1.0 - (trunc == CENTER ?  powInt(.5, c.getPrecision()+1) : .0));
#endif
}


/**
 *  randomize()
 */

template<class T>
void L::DigitalNet2<T>::randomize (PRNG &g)
{
#ifdef HINTLIB_IEEE_MAGIC_WORKS
   const int shift = floatBits - prec;
   union { floatType f; T i; } xx;
   xx.f = 1.;
#endif

   for (unsigned d = 0; d < getDimension(); ++d)
   {
      T x = 0;
      for (unsigned b = 0; b < prec; ++b)
      {
         x = 2 * x + g.equidist(2);
      }

#ifdef HINTLIB_IEEE_MAGIC_WORKS
      if (      std::numeric_limits<floatType>::digits
             >= std::numeric_limits<real>::digits
          || mode == DIRECT)
      {
         xStart [d] = (x << shift) | xx.i;
      }
      else
#endif
      xStart [d] = x;
   }
}


/**
 *  resetX()
 */

template<class T>
void L::DigitalNet2<T>::resetX (Index nn)
{
   std::copy (&xStart[0], &xStart[getDimension()], &x[0]);

   // Recalculate x

   for (unsigned r = 0; nn != 0; ++r, nn >>= 1)
   {
      if (nn & 1)
      {
         for (unsigned d = 0; d < getDimension(); ++d)  x [d] ^= c(d,r);
      }
   }
}


/**
 *  getOptimalNumber()
 */

template<class T>
Index L::DigitalNet2<T>::getOptimalNumber(Index max) const
{
   if (max >= getSize())  return getSize();
   if (trunc == FULL)  return max;

   return (max == 0) ? 0 : Index(1) << ms1(max);
}


/*****************************************************************************/
/***         Digital Net 2 Gray <T>                                        ***/
/*****************************************************************************/

/**
 *  Constructor
 */

template<class T>
L::DigitalNet2Gray<T>::DigitalNet2Gray
   (const GeneratorMatrix2<T> &gm, const Hypercube &_h,
    unsigned m, Index i, bool equi, DigitalNet::Truncation t, bool correct)
: DigitalNet2<T> (gm, _h, m, i, equi, t)
{
   if (correct)  c.prepareForGrayCode();
}

template<class T>
L::DigitalNet2Gray<T>::DigitalNet2Gray
   (const GeneratorMatrix2<T> &gm, const Hypercube &_h, bool correct)
: DigitalNet2<T> (gm, _h, gm.getM(), 0, false, FULL)
{
   if (correct)  c.prepareForGrayCode();
}


/*****************************************************************************/
/***         Digital Net 2 Point Set Base                                  ***/
/*****************************************************************************/

// setCube()

void L::DigitalNet2PointSetBase::setCube (const Hypercube *_h)
{
   h = _h;
}

//  get Optimal Number()

Index L::DigitalNet2PointSetBase::getOptimalNumber
   (Index max, const Hypercube &)
{
   if (trunc == DigitalNet::FULL)  return max;
   if (max == 0)  return 0;

   int m = std::min (ms1(max), int (gm.getM()));

   return Index(1) << m;
}

// calculate M()

int L::DigitalNet2PointSetBase::calculateM (Index num) const
{
   int m = ms1 (num);

   if (m >= 0)
   {
      Index size = Index(1) << m;

      if (num > size)  ++m;
   }

   return m;
}

// doJobPartition()

void L::DigitalNet2PointSetBase::doJobPartition
   (real *point, Job &job, Index num, Index begin, Index end)
{
   DigitalNet2Gray<BaseType> net (gm, *h, calculateM (num), index, equi, trunc);
   qmcDoJob (point, net, job, begin, end);
}

// doJobRep()

bool L::DigitalNet2PointSetBase::doJobRep
   (real *point, ReportingJob &job, Index num)
{
   DigitalNet2Gray<BaseType> net (gm, *h, calculateM (num), index, equi, trunc);
   return qmcDoJob (point, net, job, 0, num);
}

/*****************************************************************************/
/***         Digital Net 2 Point Set <real>                                ***/
/*****************************************************************************/

// integrate()

void L::DigitalNet2PointSet<L::real>::integratePartition
   (real* point, Function &f, Index num, Index begin, Index end, Stat& stat)
{
   DigitalNet2Gray<BaseType> net (gm, *h, calculateM (num), index, equi, trunc);
   qmcIntegration (point, net, f, begin, end, stat);
}

void L::DigitalNet2PointSet<L::real>::integratePartition
   (real* point, Function &f, Index num, Index begin, Index end, StatVar& stat)
{
   DigitalNet2Gray<BaseType> net (gm, *h, calculateM (num), index, equi, trunc);
   qmcIntegration (point, net, f, begin, end, stat);
}


/*****************************************************************************/
/***         Digital Seq 2 Point Set Base                                  ***/
/*****************************************************************************/

// Destructor

L::DigitalSeq2PointSetBase::~DigitalSeq2PointSetBase ()
{
   delete net;
}

//  get Optimal Number()

Index L::DigitalSeq2PointSetBase::getOptimalNumber
   (Index max, const Hypercube &)
{
   return std::min (max, Index(1) << gm.getM());
}

// checkSize()

void L::DigitalSeq2PointSetBase::checkSize (
   const DigitalNet2<BaseType> &net, Index end) const
{
   if (end > net.getSize())  throw NumbersExhausted();
}

// setCube()

void L::DigitalSeq2PointSetBase::setCube (const Hypercube *h)
{
   if (net)  net->setCube (*h);
   else      net = new DigitalNet2Gray<BaseType> (gm, *h);
}

// doJobPartition()

void L::DigitalSeq2PointSetBase::doJobPartition
   (real *point, Job &job, Index num, Index begin, Index end)
{
   if (! reset)
   {
      begin  += offset;
      end    += offset;
      offset += num;
   }

   checkSize(*net, num);
   qmcDoJob (point, *net, job, begin, end);
}

// doJobRep()

bool L::DigitalSeq2PointSetBase::doJobRep
   (real *point, ReportingJob &job, Index num)
{
   if (! reset)  offset += num;

   checkSize (*net, num);
   return qmcDoJob (point, *net, job, 0, num);
}

/*****************************************************************************/
/***         Digital Seq 2 Point Set <real>                                ***/
/*****************************************************************************/

// integrate()

void L::DigitalSeq2PointSet<L::real>::integratePartition
   (real* point, Function &f, Index num, Index begin, Index end, Stat& stat)
{
   if (! reset)
   {
      begin  += offset;
      end    += offset;
      offset += num;
   }

   checkSize (*net, num);
   qmcIntegration (point, *net, f, begin, end, stat);
}

void L::DigitalSeq2PointSet<L::real>::integratePartition
   (real* point, Function &f, Index num, Index begin, Index end, StatVar& stat)
{
   if (! reset)
   {
      begin  += offset;
      end    += offset;
      offset += num;
   }

   checkSize (*net, num);
   qmcIntegration (point, *net, f, begin, end, stat);
}

#ifdef HINTLIB_32BIT
template class L::DigitalNet2<L::u32>;
template class L::DigitalNet2Gray<L::u32>;
#endif

template class L::DigitalNet2<L::u64>;
template class L::DigitalNet2Gray<L::u64>;


