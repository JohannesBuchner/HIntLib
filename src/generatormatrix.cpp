/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_OSTREAM
#include <ostream>
#else
#include <iostream>
#endif

#include <iomanip>

#include <HIntLib/generatormatrix.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/hlmath.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;


/*****************************************************************************/
/***     Generator Matrix                                                  ***/
/*****************************************************************************/

/**
 *  setParameters()
 */

void L::GeneratorMatrix::setParameters (const GeneratorMatrix& gm)
{
   base = gm.base;
   dim  = gm.dim;
   m    = gm.m;
   prec = gm.prec;
}


/**
 *  getDefaultM ()
 */

int L::GeneratorMatrix::getDefaultM (int base)
{
   const Index maxNetSize
      = std::numeric_limits<Index>::digits <= 49
      ? std::numeric_limits<Index>::max()
      : (Index(1) << 48) - 1;

   return logInt (maxNetSize, Index (base));
}

/**
 *  getDefaultPrec ()
 */

int L::GeneratorMatrix::getDefaultPrec (int base)
{
   return int (HINTLIB_MN ceil(
      HINTLIB_MN log(2.0) / HINTLIB_MN log(double(base))
                          * double(std::numeric_limits<real>::digits - 1)));
}


/**
 *  setDigit()
 *  vSetPackedRowVector()
 *
 *  By default, a Generator Matrix cannot be written to.
 */

void L::GeneratorMatrix::setDigit (int, int, int, int)
{
   throw InternalError (__FILE__, __LINE__);
}

void L::GeneratorMatrix::vSetPackedRowVector (int, int, u64)
{
   throw InternalError (__FILE__, __LINE__);
}


/**
 *  print()
 */

void L::GeneratorMatrix::print (std::ostream &o) const
{
   o << "Base=" << getBase()
     << " Dim=" << getDimension()
     << " m=" << getM()
     << " prec=" << getPrec() << '\n';

   for (int d = 0; d < getDimension(); ++d)
   {
      o << "Dimension " << d << ":\n";

      printDimension (o, d);
   }
}

#ifdef HINTLIB_BUILD_WCHAR
void L::GeneratorMatrix::print (std::wostream &o) const
{
   o << L"Base=" << getBase()
     << L" Dim=" << getDimension()
     << L" m=" << getM()
     << L" prec=" << getPrec() << L'\n';

   for (int d = 0; d < getDimension(); ++d)
   {
      o << L"Dimension " << d << L":\n";

      printDimension (o, d);
   }
}
#endif


/**
 *  printDimension()
 */

void L::GeneratorMatrix::printDimension (std::ostream &o, int d) const
{
   for (int b = 0; b < getPrec(); ++b)
   {
      printRowVector (o, d, b);
      o << '\n';
   }
}

#ifdef HINTLIB_BUILD_WCHAR
void L::GeneratorMatrix::printDimension (std::wostream &o, int d) const
{
   for (int b = 0; b < getPrec(); ++b)
   {
      printRowVector (o, d, b);
      o << L'\n';
   }
}
#endif


/**
 *  printRowVector ()
 */

void L::GeneratorMatrix::printRowVector (std::ostream &o, int d, int b) const
{
   int size = (base < 10) ? 0 : logInt (base, 10) + 2;

   for (int r = 0; r < getM(); ++r)
   {
      if (size)  o << std::setw (size - (r == 0));
      o << getDigit (d,r,b);
   }
}

#ifdef HINTLIB_BUILD_WCHAR
void L::GeneratorMatrix::printRowVector (std::wostream &o, int d, int b) const
{
   int size = (base < 10) ? 0 : logInt (base, 10) + 2;

   for (int r = 0; r < getM(); ++r)
   {
      if (size)  o << std::setw (size - (r == 0));
      o << getDigit (d,r,b);
   }
}
#endif


/**
 *  printColumnVector ()
 */

void L::GeneratorMatrix::printColumnVector (std::ostream &o, int d, int r) const
{
   int size = (base < 10) ? 1 : logInt (base, 10) + 2;

   for (int b = 0; b < getPrec(); ++b)
   {
      o << std::setw (size) << getDigit (d,r,b);
   }
}

#ifdef HINTLIB_BUILD_WCHAR
void L::GeneratorMatrix::printColumnVector (
      std::wostream &o, int d, int r) const
{
   int size = (base < 10) ? 1 : logInt (base, 10) + 2;

   for (int b = 0; b < getPrec(); ++b)
   {
      o << std::setw (size) << getDigit (d,r,b);
   }
}
#endif


/**
 *  operator==
 */

bool L::operator== (const L::GeneratorMatrix &gm1,
                    const L::GeneratorMatrix &gm2)
{
   if (gm1.getDimension() != gm2.getDimension()
    || gm1.getM()         != gm2.getM()
    || gm1.getBase()      != gm2.getBase()
    || gm1.getPrec()      != gm2.getPrec())
   {
      return false;
   }

   for (int d = 0; d < gm1.getDimension(); ++d)
   {
      for (int r = 0; r < gm1.getM(); ++r)
      {
         for (int b = 0; b < gm1.getPrec(); ++b)
         {
            if (gm1.getDigit(d, r, b) != gm2.getDigit(d, r, b)) return false;
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

   if (o.getPrec() < n.getPrec())
      throw GM_CopyPrec (n.getPrec(), o.getPrec());

   if (o.getM() < n.getM())
      throw GM_CopyM (n.getM(), o.getM());
}

void L::GeneratorMatrix::checkCopyDim (
   const GeneratorMatrix &o, const GeneratorMatrix &n, int offset)
{
   checkCopy (o, n);

   if (o.getDimension() + offset < n.getDimension())
      throw GM_CopyDim (n.getDimension(), o.getDimension());
}


/**
 *  assign Matrix
 *
 *  Assigns a submatrix of a Generator Matrix to another Generator Matrix
 */

void L::assign (
   const GeneratorMatrix & src, int srcDim,
         GeneratorMatrix & dst, int dstDim)
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

   for (int r = 0; r < dst.getM(); ++r)
   {
      for (int b = 0; b < dst.getPrec(); ++b)
      {
         dst.setDigit (dstDim, r, b, src.getDigit (srcDim, r, b));
      }
   }
}


void L::assign (const GeneratorMatrix &src, GeneratorMatrix &dst)
{
   for (int d = 0; d < dst.getDimension(); ++d)  assign (src, d, dst, d);
}

