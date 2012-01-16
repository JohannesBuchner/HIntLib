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

#include <iostream>
#include <iomanip>
#include <algorithm>

#include <HIntLib/generatormatrix.h>
#include <HIntLib/mymath.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;

using std::numeric_limits;


/*****************************************************************************/
/***     GM Copy                                                           ***/
/*****************************************************************************/

L::GMCopy::GMCopy()
   : dimValue          (-1),
     maxDimValue       (-1),
     mValue            (-1),
     maxMValue         (-1),
     totalPrecValue    (-1),
     maxTotalPrecValue (-1),
     vecValue          (-1),
     equiValue      (false)
{}

unsigned L::GMCopy::getDimension (const GeneratorMatrix& m) const
{
   if (dimValue >= 0)  return dimValue;

   int def = m.getDimension() + equiValue;

   return (maxDimValue >= 0 && maxDimValue < def)  ? maxDimValue : def;
}

unsigned L::GMCopy::getM (const GeneratorMatrix& m) const
{
   if (mValue >= 0)  return mValue;
   if (maxMValue >= 0 && unsigned(maxMValue) < m.getM())  return maxMValue;
   return m.getM();
}

unsigned L::GMCopy::getTotalPrec (const GeneratorMatrix& m, unsigned max) const
{
   if (totalPrecValue >= 0)  return totalPrecValue;
   if (maxTotalPrecValue >= 0)
   {
      return std::min (unsigned (maxTotalPrecValue), m.getTotalPrec());
   }

   if (max > 0)
      return std::min (m.getTotalPrec(), max / ms1(m.getBase()));
   else
      return m.getTotalPrec();
}

unsigned L::GMCopy::getVectorization 
   (const GeneratorMatrix& m, unsigned /* bits */) const
{
   switch (vecValue)
   {
   case -1: return 1;
   case -2: return m.getVectorization(); 
   default: return vecValue;
   };
}

void L::GMCopy::checkNoVec () const
{
   if (vecValue >= 0)  throw FIXME (__FILE__, __LINE__);
}


/*****************************************************************************/
/***     Generator Matrix                                                  ***/
/*****************************************************************************/

/**
 *  Copy Constructor
 */

L::GeneratorMatrix::GeneratorMatrix (const GeneratorMatrix &gm)
: base      (gm.base),
  vec       (gm.vec),
  prec      (gm.prec),
  dim       (gm.dim),
  m         (gm.m),
  totalPrec (gm.totalPrec)
{}


/**
 *  setDigit()
 *  setVector()
 *
 *  By default, a Generator Matrix can not be written to.
 *
 *  Non-dummy implementations are available in
 *     MutableGeneratorMatrixGen<T> and
 *     MutalbeGeneratorMatrix2<T>.
 */

void L::GeneratorMatrix::setDigit (unsigned, unsigned, unsigned, unsigned)
{
   throw InternalError (__FILE__, __LINE__);
}

void L::GeneratorMatrix::setVector (unsigned, unsigned, unsigned, u64)
{
   throw InternalError (__FILE__, __LINE__);
}


/**
 *  dump()
 */

void L::GeneratorMatrix::dump (std::ostream &o) const
{
   o << "Base=" << getBase()
     << " Dim=" << getDimension()
     << " m=" << getM()
     << " TotalPrec=" << getTotalPrec()
     << " (stored in " << getPrec() << " blocks of " << getVectorization()
     << " digits)\n";

   unsigned size = logInt (base, 10u) + 1;
   if (base >= 10)  ++size;
   
   for (unsigned d = 0; d < getDimension(); ++d)
   {
      o << "Dimension " << d << ":\n";

      for (unsigned p = 0; p < getTotalPrec(); ++p)
      {
         for (unsigned r = 0; r < getM(); ++r)
         {
            o << std::setw (size) << getDigit (d,r,p);
         }

         o << '\n';
      }
   }
}


/**
 *  vectorDump()
 */

void L::GeneratorMatrix::vectorDump (std::ostream &o) const
{
   for (unsigned d = 0; d < dim; ++d)
   {
      o << "Dimension " << d << ":\n";

      for (unsigned r = 0; r < m; ++r)
      {
         for (unsigned p = 0; p < getPrec(); ++p)
         {
            if (p > 0)  o << ' ';
            o << getVector (d,r,p);
         }
         o << '\n';
      }
   }
}


/**
 *  operator==
 */

bool L::operator== (const L::GeneratorMatrix &gm1,
                    const L::GeneratorMatrix &gm2)
{
   if (gm1.getDimension() != gm2.getDimension()
    || gm1.getM() != gm2.getM()
    || gm1.getBase() != gm2.getBase()
    || gm1.getTotalPrec() != gm2.getTotalPrec())
   {
      return false;
   }

   if (gm1.getVectorization() == gm2.getVectorization())
   {
      for (unsigned d = 0; d < gm1.getDimension(); ++d)
      {
         for (unsigned r = 0; r < gm1.getM(); ++r)
         {
            for (unsigned b = 0; b < gm1.getPrec(); ++b)
            {
               if (gm1.getVector(d, r, b) != gm2.getVector(d, r, b))
                  return false;
            }
         }
      }
   }
   else
   {
      for (unsigned d = 0; d < gm1.getDimension(); ++d)
      {
         for (unsigned r = 0; r < gm1.getM(); ++r)
         {
            for (unsigned b = 0; b < gm1.getTotalPrec(); ++b)
            {
               if (gm1.getDigit(d, r, b) != gm2.getDigit(d, r, b)) return false;
            }
         }
      }
   }

   return true;
}


/**
 *  checkCopy()
 */

void L::GeneratorMatrix::checkCopy
   (const GeneratorMatrix &o, const GeneratorMatrix &n)
{
   if (o.getBase() != n.getBase())
      throw GM_CopyBase (n.getBase(), o.getBase());

   if (o.getTotalPrec() < n.getTotalPrec())
      throw GM_CopyPrec (n.getTotalPrec(), o.getTotalPrec());

   if (o.getM() < n.getM())
      throw GM_CopyM (n.getM(), o.getM());
}

void L::GeneratorMatrix::checkCopyDim (
   const GeneratorMatrix &o, const GeneratorMatrix &n, int offset)
{
   checkCopy (o, n);

   if (int (o.getDimension()) + offset < int (n.getDimension()))
      throw GM_CopyDim (n.getDimension(), o.getDimension());
}


/**
 *  assign Matrix
 *
 *  Assigns a submatrix of a Generator Matrix to another Generator Matrix
 */

void L::assign (
   const GeneratorMatrix & src, unsigned srcDim,
         GeneratorMatrix & dst, unsigned dstDim)
{
   GeneratorMatrix::checkCopy (src, dst);

   if (dstDim >= dst.getDimension())
   {
      throw DimensionTooHigh (dstDim, dst.getDimension() - 1);
   }

   if (srcDim >= src.getDimension())
   {
      throw DimensionTooHigh (dstDim, dst.getDimension() - 1);
   }

   for (unsigned r = 0; r < dst.getM(); ++r)
   {
      for (unsigned b = 0; b < dst.getTotalPrec(); ++b)
      {
         dst.setDigit (dstDim, r, b, src.getDigit (srcDim, r, b));
      }
   }
}


void L::assign (const GeneratorMatrix &src, GeneratorMatrix &dst)
{
   for (unsigned d = 0; d < dst.getDimension(); ++d)  assign (src, d, dst, d);
}


