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

#include <HIntLib/mersennetwister.h>
#include <HIntLib/mcpointset.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;

using L::Index;


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
  DigitalNet (_c.getBase(), _m),
  alg(),
  scalAlg (alg.getScalarAlgebra()),
  totalPrec (std::min(
     (_trunc == FULL) ? _c.getTotalPrec() : m,   // what we want
     unsigned (std::numeric_limits<real>::digits - 1)
  )), // what makes sense

  c (_c, GMCopy().dim(getDimension()).m(m).totalPrec(totalPrec).equi(equi)),
  x      (getDimension()),
  xStart (getDimension(), 0),
  trunc (_trunc),
  trivialScale (HINTLIB_MN pow(real(.5), int(c.getTotalPrec()))),
  ss (h.getDimension())
#ifdef HINTLIB_IEEE_MAGIC_WORKS
  , ssMagic (h.getDimension()),
  mode (STANDARD)
#endif
{
   if (c.getBase() != 2)  throw FIXME (__FILE__, __LINE__);

   // Initialize Shift Scale

   setCube (h);
   
   // Initialize xStart vector

   const int dd = equi;
   unsigned i = m;
   Index indexCopy = index;
   int shift = _c.getTotalPrec() - c.getTotalPrec();

   while (index)
   {
      if (i >= _c.getM())
      {
         throw NetIndexTooHigh (indexCopy, _c.getM(), m);
      }
      
      if (index & 1)
      {
         for (unsigned d = dd; d < c.getDimension(); ++d)
         {
            alg.addTo (xStart [d], (_c(d-dd, i) >> shift));
         }
      }

      index >>= 1;
      ++i;
   }

#ifdef HINTLIB_IEEE_MAGIC_WORKS
   // Prepare for direct updates, if we can

   if (c.getTotalPrec() <= floatBits)
   {
      mode = DIRECT;

      const int shift = floatBits - c.getTotalPrec();
      union { floatType f; T i; } x;
      x.f = 1.;

      for (unsigned d = 0; d < getDimension(); ++d)
      {
         xStart [d] = (xStart[d] << shift) | x.i;
      }

      c.adjustTotalPrec (floatBits);
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
  double centerShift = trunc == CENTER ? -.5 : .0;

  ss.set (h, centerShift, 
             centerShift + HINTLIB_MN pow(real(2.0), int (c.getTotalPrec())));

#ifdef HINTLIB_IEEE_MAGIC_WORKS
  double centerShiftMagic = trunc == CENTER
     ? - HINTLIB_MN pow (real(.5), int (c.getTotalPrec()) + 1) : real(0);

  ssMagic.set (h, 1.0 + centerShiftMagic, 2.0 + centerShiftMagic);
#endif
}


/**
 *  randomize()
 */

template<class T>
void L::DigitalNet2<T>::randomize (PRNG &g)
{
#ifdef HINTLIB_IEEE_MAGIC_WORKS
   const int shift = floatBits - totalPrec;
   union { floatType f; T i; } xx;
   xx.f = 1.;
#endif

   for (unsigned d = 0; d < getDimension(); ++d)
   {
      T x = 0;
      for (unsigned b = 0; b < totalPrec; ++b)
      {
         x = (x << 1) | g.equidist(scalAlg.size());
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
      const unsigned char digit = nn & 1;

      if (! scalAlg.is0 (digit))
      {
         for (unsigned d = 0; d < getDimension(); ++d)
         {
            alg.addTo (x [d], alg.mul (c(d,r), digit));
         }
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
/***         Digital Net 2 Point Set Base                                  ***/
/*****************************************************************************/

// performRandomization

void
L::Digital2PointSet::performRandomization (DigitalNet2Gray<BaseType> &net)
{
   if (isRandomized)
   {
      PRNGImp<MersenneTwister> mt (seed);
      net.randomize (mt);
   }
}

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
   DigitalNet2Gray<BaseType> net
      (gm, *h, calculateM (num), index, equi, trunc);
   performRandomization (net);
   qmcDoJob (point, net, job, begin, end);
}

// doJobRep()

bool L::DigitalNet2PointSetBase::doJobRep
   (real *point, ReportingJob &job, Index num)
{
   DigitalNet2Gray<BaseType> net (gm, *h, calculateM (num), index, equi, trunc);
   performRandomization (net);
   return qmcDoJob (point, net, job, 0, num);
}

/*****************************************************************************/
/***         Digital Net 2 Point Set <real>                                ***/
/*****************************************************************************/

// integrate()

void L::DigitalNet2PointSet<L::real>::integratePartition
   (real* point, Integrand &f, Index num, Index begin, Index end, Stat& stat)
{
   DigitalNet2Gray<BaseType> net (gm, *h, calculateM (num), index, equi, trunc);
   performRandomization (net);
   qmcIntegration (point, net, f, begin, end, stat);
}

void L::DigitalNet2PointSet<L::real>::integratePartition
   (real* point, Integrand &f, Index num, Index begin, Index end, StatVar& stat)
{
   DigitalNet2Gray<BaseType> net (gm, *h, calculateM (num), index, equi, trunc);
   performRandomization (net);
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
   else
   {
      net = new DigitalNet2Gray<BaseType> (gm, *h);
      performRandomization (*net);
   }
}

// randomzie()

void
L::DigitalSeq2PointSetBase::randomize (unsigned _seed)
{
   Digital2PointSet::randomize (_seed);
   if (net)  performRandomization (*net);
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
   (real* point, Integrand &f, Index num, Index begin, Index end, Stat& stat)
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
   (real* point, Integrand &f, Index num, Index begin, Index end, StatVar& stat)
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

namespace HIntLib
{
   template class DigitalNet2 <u32>;
#ifdef HINTLIB_U32_NOT_EQUAL_U64
   template class DigitalNet2 <u64>;
#endif
}


