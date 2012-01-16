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

#include <HIntLib/generatormatrixgen.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;

using std::numeric_limits;


/**
 *  assign Matrix
 *
 *  Assigns a submatrix of a Generator Matrix to another Generator Matrix
 */

template<class T>
void L::assign (
   const GeneratorMatrix & src, unsigned srcDim,
         MutableGeneratorMatrixGen<T> & dst, unsigned dstDim)
{
#if 0
   cerr <<"A3" <<endl;
   cerr << "dim = " << srcM.getDimension() << " "<<dstM.getDimension() << endl;
   cerr << "m = " << srcM.getM() << " "<<dstM.getM() << endl;
   cerr << "totalPrec = " << srcM.getTotalPrec() << " "<<dstM.getTotalPrec() << endl;
   cerr << "base = " << srcM.getBase() << " "<<dstM.getBase() << endl;
   cerr << "prec = " << srcM.getPrec() << " "<<dstM.getPrec() << endl;
   cerr << "vec = " << srcM.getVectorization() << " "<<dstM.getVectorization() << endl;
#endif

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
      for (unsigned b = 0; b < dst.getPrec(); ++b)
      {
         dst.setd (dstDim, r, b, src.getDigit (srcDim, r, b));
      }
   }

   // FIXME  Optimization for src.getVectorization() != 1 is missing.
}

template<class T>
void L::assign (const GeneratorMatrix &src, MutableGeneratorMatrixGen<T> &dst)
{
   for (unsigned d = 0; d < dst.getDimension(); ++d)  assign (src, d, dst, d);
}

/*****************************************************************************/
/***     Generator Matrix Gen <T>                                          ***/
/*****************************************************************************/

/**
 *  getDigit ()
 */

template<class T>
unsigned
L::GeneratorMatrixGen<T>::getDigit (unsigned d, unsigned r, unsigned b) const
{
   if (numeric_limits<T>::digits > numeric_limits<unsigned>::digits)
   {
      throw InternalError (__FILE__, __LINE__);
   }

   return unsigned(operator() (d, r, b));
}


/**
 *  getVector ()
 */

template<class T>
L::u64
L::GeneratorMatrixGen<T>::getVector (unsigned d, unsigned r, unsigned b) const
{
   if (numeric_limits<T>::digits > numeric_limits<u64>::digits)
   {
      throw InternalError (__FILE__, __LINE__);
   }

   return u64(operator() (d, r, b));
}


/**
 *  operator==
 */

template<class T>
bool L::GeneratorMatrixGen<T>::operator==
   (const L::GeneratorMatrixGen<T> &gm) const
{
   if (dim != gm.dim || m != gm.m || totalPrec != gm.totalPrec)  return false;

   return std::equal(c, c + dimPrec*m , gm.c);
}


/**
 *  checkBase()
 */

template<class T>
void L::GeneratorMatrixGen<T>::checkBase() const
{
   if (getBase()-1 > numeric_limits<T>::max())
   {
      throw GM_BaseTooLarge (getBase(), unsigned (numeric_limits<T>::max()));
   }
}


/*****************************************************************************/
/***     Mutable Generator Matrix Gen <T>                                  ***/
/*****************************************************************************/


/**
 *  setDigit ()
 */

template<class T>
void L::MutableGeneratorMatrixGen<T>::setDigit
   (unsigned d, unsigned r, unsigned b, unsigned x)
{
   if (   numeric_limits<T>::digits < numeric_limits<unsigned>::digits
       && x >= getBase())
   {
      throw InternalError (__FILE__, __LINE__);
   }

   setd (d, r, b, T(x));
}


/**
 *  setVector ()
 */

template<class T>
void L::MutableGeneratorMatrixGen<T>::setVector
   (unsigned d, unsigned r, unsigned b, u64 x)
{
   if (   numeric_limits<T>::digits < numeric_limits<u64>::digits
       && x >= getBase())
   {
      throw InternalError (__FILE__, __LINE__);
   }

   setd (d, r, b, T(x));
}


/**
 *  make Equidistributed Coordinate()
 */

template<class T>
void L::MutableGeneratorMatrixGen<T>::makeEquidistributedCoordinate (unsigned d)
{
   makeZeroMatrix (d);
   
   for (unsigned r = 0; r < m; ++r)
   {
      if (m-(r+1) < totalPrec)  setd (d, r, m-(r+1), 1);
   }
}


/**
 *  make Identity Matrix ()
 */

template<class T>
void L::MutableGeneratorMatrixGen<T>::makeIdentityMatrix  (unsigned d)
{
   makeZeroMatrix (d);
   
   const unsigned ub = std::min (m, totalPrec);
   for (unsigned r = 0; r < ub; ++r)  setd (d, r, r, 1);
}


/**
 *  make Zero Matrix ()
 */

template<class T>
void L::MutableGeneratorMatrixGen<T>::makeZeroMatrix (unsigned d)
{
   for (unsigned r = 0; r < m; ++r)
   {
      T* base = c + r*dimPrec + d*prec;
      std::fill (base, base + prec, 0);
   }
}

template<class T>
void L::MutableGeneratorMatrixGen<T>::makeZeroMatrix ()
{
   std::fill (c, c + dimPrec * m, 0);
}


/*****************************************************************************/
/***     Heap Allocated  Matrix Gen <T>                                    ***/
/*****************************************************************************/

/**
 *  Constructor
 */

template<class T>
L::HeapAllocatedGeneratorMatrixGen<T>::HeapAllocatedGeneratorMatrixGen
   (unsigned _base, unsigned _dim, unsigned _m)
: MutableGeneratorMatrixGen<T> (
      _base, _dim, _m,
      std::min (_m, unsigned (ceil(log(2.0) / log(double(_base))
                              * double(numeric_limits<real>::digits - 1)))))
{
   allocate();
   makeZeroMatrix();
}

template<class T>
L::HeapAllocatedGeneratorMatrixGen<T>::HeapAllocatedGeneratorMatrixGen
   (unsigned _base, unsigned _dim)
: MutableGeneratorMatrixGen<T> (
   _base, _dim,
   logInt (std::numeric_limits<Index>::max(), Index (_base)),
   unsigned (ceil(log(2.0) / log(double(_base))
               * double(numeric_limits<real>::digits - 1))))
{
   allocate();
   makeZeroMatrix();
}


/**
 *  allocate()
 */

template<class T>
void L::HeapAllocatedGeneratorMatrixGen<T>::allocate()
{
   if (c)  throw InternalError (__FILE__, __LINE__);
   c = new T [dimPrec * m];
}


/*****************************************************************************/
/***     Generator Matrix Gen Copy <T>                                     ***/
/*****************************************************************************/

/**
 *  Copy constructor
 */

template<class T>
L::GeneratorMatrixGenCopy<T>::GeneratorMatrixGenCopy
   (const GeneratorMatrixGen<T> &gm)
: HeapAllocatedGeneratorMatrixGen<T> (gm, true)
{
   std::copy (gm.getMatrix(), gm.getMatrix() + m*dimPrec, c);
}


/**
 *  Copy from arbitrary Generator Matrix
 */

template<class T>
L::GeneratorMatrixGenCopy<T>::GeneratorMatrixGenCopy (const GeneratorMatrix &gm)
   : HeapAllocatedGeneratorMatrixGen<T> (
      gm.getBase(),
      gm.getDimension(),
      gm.getM(),
      gm.getTotalPrec(),
      true)
{
   for (unsigned d = 0; d < dim; ++d)  assign (gm, d, *this, d);
}


/**
 *  Constructor  (using GMCopy)
 */

template<class T>
L::GeneratorMatrixGenCopy<T>::GeneratorMatrixGenCopy (
   const GeneratorMatrixGen<T> &gm, const GMCopy &c)
: HeapAllocatedGeneratorMatrixGen<T>
   (gm.getBase(),
    c.getDimension (gm),
    c.getM (gm),
    c.getTotalPrec (gm),
    false)
{
   const int dd = c.getEqui();
   checkCopyDim (gm, *this, dd);

   allocate();

   for (unsigned d = 0; d < dim - dd; ++d)  assign (gm, d, *this, d + dd);

   if (dd && dim)  makeEquidistributedCoordinate (0);
}

template<class T>
L::GeneratorMatrixGenCopy<T>::GeneratorMatrixGenCopy (
   const GeneratorMatrix &gm, const GMCopy& c)
: HeapAllocatedGeneratorMatrixGen<T>
   (gm.getBase(),
    c.getDimension (gm),
    c.getM (gm),
    c.getTotalPrec (gm),
    false)
{
   const int dd = c.getEqui();
   checkCopyDim (gm, *this, dd);

   allocate();

   for (unsigned d = 0; d < dim - dd; ++d)  assign (gm, d, *this, d + dd);

   if (dd && dim)  makeEquidistributedCoordinate (0);
}


namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template class GeneratorMatrixGen<X>; \
   template class MutableGeneratorMatrixGen<X>; \
   template class HeapAllocatedGeneratorMatrixGen<X>; \
   template class GeneratorMatrixGenCopy<X>; \
   template void assign ( \
      const GeneratorMatrix &, unsigned, \
      MutableGeneratorMatrixGen<X> &, unsigned); \
   template void assign ( \
      const GeneratorMatrix &, MutableGeneratorMatrixGen<X> &);

   HINTLIB_INSTANTIATE (unsigned char)
   HINTLIB_INSTANTIATE (unsigned short)
   HINTLIB_INSTANTIATE (u32)
#undef HINTLIB_INSTANTIATE
}

