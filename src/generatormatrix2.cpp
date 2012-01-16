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

#include <HIntLib/generatormatrix2.h>
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
         MutableGeneratorMatrix2<T> & dst, unsigned dstDim)
{
#if 0
   cerr <<"A2" <<endl;
   cerr << "dim = " << src.getDimension() << " "<<dst.getDimension() << endl;
   cerr << "m = " << src.getM() << " "<<dst.getM() << endl;
   cerr << "totalPrec = " << src.getTotalPrec() << " "<<dst.getTotalPrec() << endl;
   cerr << "base = " << src.getBase() << " "<<dst.getBase() << endl;
   cerr << "prec = " << src.getPrec() << " "<<dst.getPrec() << endl;
   cerr << "vec = " << src.getVectorization() << " "<<dst.getVectorization() << endl;
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

   const unsigned srcVectorization = src.getVectorization();
   const unsigned leadingDigits = src.getNumOfLeadingDigits();

   if (leadingDigits >= dst.getTotalPrec())
   {
      const unsigned shift = leadingDigits - dst.getTotalPrec();

      for (unsigned r = 0; r < dst.getM(); ++r)
      {
         dst.setv (dstDim, r, src.getVector (srcDim, r, 0) >> shift);
      }
   }
   else   // leadingDigits < dst.getTotalPrec()
   {
      for (unsigned r = 0; r < dst.getM(); ++r)
      {
         T x = src.getVector (srcDim, r, 0);

         unsigned missingDigits = dst.getTotalPrec() - leadingDigits;
         unsigned next = 1;

         while (missingDigits >= srcVectorization)
         {
            x = (x << srcVectorization) | src.getVector (srcDim, r, next++);
            missingDigits -= srcVectorization;
         }

         if (missingDigits)
         {
            x = (x << missingDigits)
              | (src.getVector (srcDim, r, next)
                   >> (srcVectorization - missingDigits));
         }

         dst.setv (dstDim, r, x);
      }
   }
}

template<class T>
void L::assign (const GeneratorMatrix &src, MutableGeneratorMatrix2<T> &dst)
{
   for (unsigned d = 0; d < dst.getDimension(); ++d)  assign (src, d, dst, d);
}


/*****************************************************************************/
/***     Generator Matrix 2 Base                                           ***/
/*****************************************************************************/

/**
 *  CArrayDump()
 */

void L::GeneratorMatrix2Base::CArrayDump (std::ostream &o) const
{
   const char* suffix;
   if (vec > 32)       suffix = "ull";
   else if (vec > 16)  suffix = "ul";
   else                suffix = "u";
   
   o << "{\n";

   for (unsigned r = 0; r < m; ++r)
   {
      o << "   {  // bit " << r << '\n';

      for (unsigned d = 0; d < dim; ++d)
      {
         o << std::setw(22) << getVector(d,r,0) << suffix << ", // dim ="
           << std::setw(3) << d+1 << '\n';
      }
      o << "   },\n";
   }

   o << "};\n";
}


/*****************************************************************************/
/***     Generator Matrix 2 <T>                                            ***/
/*****************************************************************************/

/**
 *  checkVec()
 */

template<class T>
void L::GeneratorMatrix2<T>::checkTotalPrec() const
{
   if (totalPrec > unsigned (std::numeric_limits<T>::digits))
   {
      throw GM_PrecTooHigh (base, numeric_limits<T>::digits, totalPrec);
   }
}

/**
 *  operator==
 */

template<class T>
bool L::GeneratorMatrix2<T>::operator== (const L::GeneratorMatrix2<T> &gm) const
{
   return dim       == gm.dim
       && m         == gm.m
       && totalPrec == gm.totalPrec
       && std::equal(c, c + dim*m , gm.c);
}


/**
 *  getDigit()
 */

template<class T>
unsigned L::GeneratorMatrix2<T>::getDigit
   (unsigned d, unsigned r, unsigned b) const
{
   return unsigned (operator() (d, r, b));
}


/**
 *  getVector()
 */

template<class T>
L::u64 L::GeneratorMatrix2<T>::getVector(unsigned d, unsigned r, unsigned) const
{
   if (numeric_limits<T>::digits > numeric_limits<u64>::digits)
   {
      throw InternalError (__FILE__, __LINE__);
   }

   return u64 (operator() (d, r));
}


/*****************************************************************************/
/***     Mutable Generator Matrix 2 <T>                                    ***/
/*****************************************************************************/

/**
 *  setd ()
 */

template<class T>
void L::MutableGeneratorMatrix2<T>::setd
   (unsigned d, unsigned r, unsigned b, unsigned char x)
{
   T mask = T(1) << (totalPrec - b - 1);
   if (x)  c[r*dim + d] |=  mask;
   else    c[r*dim + d] &= ~mask;
}


/**
 *  setDigit ()
 */

template<class T>
void L::MutableGeneratorMatrix2<T>::setDigit
   (unsigned d, unsigned r, unsigned b, unsigned x)
{
   setd (d, r, b, x);
}


/**
 *  setVector()
 */

template<class T>
void L::MutableGeneratorMatrix2<T>::setVector
   (unsigned d, unsigned r, unsigned, L::u64 x)
{
   if (numeric_limits<T>::digits > numeric_limits<u64>::digits)
   {
      throw InternalError (__FILE__, __LINE__);
   }

   setv (d, r, x);
}


/**
 *  adjust Total Prec()
 */

template<class T>
void L::MutableGeneratorMatrix2<T>::adjustTotalPrec (unsigned n)
{
   if (n != totalPrec)
   {
      if (n < totalPrec)
      {
         unsigned shift = totalPrec - n;
         for (T* p = c; p < c + m*dim; ++p)  *p >>= shift;
         totalPrec = n;
      }
      else  // _vec > vec
      {
         unsigned shift = n - totalPrec;
         totalPrec = n;
         checkTotalPrec();
         for (T* p = c; p < c + m*dim; ++p)  *p <<= shift;
      }

   }
}


/**
 *  make Equidistributed Coordinate()
 *
 *  Sets the matrix for dimension  d  to such that an equidistributed
 *  coordinate is produiced.
 *  The Matrix is equal to a mirrored identity matrix.
 */

template<class T>
void L::MutableGeneratorMatrix2<T>::makeEquidistributedCoordinate (unsigned d)
{
   if (totalPrec >= m)
   {
      T a = T(1) << (totalPrec - m);
      for (unsigned r = 0; r < m; ++r)  setv (d, r, a), a <<= 1;
   }
   else
   {
      for (unsigned r = 0; r < m - totalPrec; ++r)  c[r*dim + d] = 0;
      T a (1);
      for (unsigned r = m - totalPrec; r < m; ++r) setv (d, r, a), a <<= 1;
   }
}


/**
 *  make Identity Matrix ()
 *
 *  Sets the matrix   d  to the identity matrix.
 *  Sets the matrix for dimensin  d  to the Halton sequence for 2.
 */

template<class T>
void L::MutableGeneratorMatrix2<T>::makeIdentityMatrix (unsigned d)
{
   T a = T(1) << (totalPrec - 1);

   for (unsigned r = 0; r < m; ++r)  c[r*dim + d] = a, a >>= 1;
}


/**
 *  make Zero Matrix ()
 *
 *  Fills the matrix for dimension  d  with 0.
 */

template<class T>
void L::MutableGeneratorMatrix2<T>::makeZeroMatrix (unsigned d)
{
   for (unsigned r = 0; r < m; ++r)  c[r*dim + d] = 0;
}

template<class T>
void L::MutableGeneratorMatrix2<T>::makeZeroMatrix ()
{
   std::fill (c, c + dim * m, 0);
}


/**
 *  prepare For Gray Code()
 *
 *  Changes the Matrix such that applying the Gray-Code algorithm leads to the
 *  same sequence as the normal algorithm applied to the original matrix.
 *
 *  This operation is equivalent with multiplying each matrix with
 *
 *       1 0 ..... 0
 *       1 1 0 ... 0
 *       1 1 1 0 . 0
 *       : : :     :
 *       : : :     0
 *       1 1 1 1 1 1
 *
 *  The operation can be reversed using restoreFromGrayCode().
 */

template<class T>
void L::MutableGeneratorMatrix2<T>::prepareForGrayCode ()
{
   for (unsigned d = 0; d < dim; ++d)
   {
      T a = c[0*dim + d];

      for (unsigned r = 1; r < m; ++r)  a = (c[r*dim + d] ^= a);
   }
}


/**
 *  restore From Gray Code()
 *
 *  This is the reverse operation to prepareForGrayCode().
 *
 *  It is equivalent with multiplying each matrix with
 *
 *       1 0 ..... 0
 *       1 1 0 ... 0
 *       0 1 1 0 . 0
 *       : : :     :
 *       : : :     0
 *       0 0 0 0 1 1
 */

template<class T>
void L::MutableGeneratorMatrix2<T>::restoreFromGrayCode ()
{
   for (unsigned d = 0; d < dim; ++d)
   {
      for (unsigned r = m - 1; r > 0; --r)
      {
         c[r*dim + d] ^= c[(r-1)*dim + d];
      }
   }
}


/*****************************************************************************/
/***     Heap Allocated Generator Matrix 2 <T>                             ***/
/*****************************************************************************/

/**
 *  Constructor
 */

template<class T>
L::HeapAllocatedGeneratorMatrix2<T>::HeapAllocatedGeneratorMatrix2
   (unsigned _dim)
: MutableGeneratorMatrix2<T> (
   _dim,
   std::numeric_limits<Index>::digits - 1,
   std::numeric_limits<real>::digits - 1)
{
   allocate();
   makeZeroMatrix();
}

template<class T>
L::HeapAllocatedGeneratorMatrix2<T>::HeapAllocatedGeneratorMatrix2
   (unsigned _dim, unsigned _m)
: MutableGeneratorMatrix2<T> (
   _dim, _m, std::min (_m, unsigned (std::numeric_limits<real>::digits - 1)))
{
   allocate();
   makeZeroMatrix();
}


/**
 *  allocate()
 */

template<class T>
void L::HeapAllocatedGeneratorMatrix2<T>::allocate()
{
   if (c)  throw InternalError (__FILE__, __LINE__);
   c = new T [dim * m];
}


/*****************************************************************************/
/***     Generator Matrix 2 Copy <T>                                       ***/
/*****************************************************************************/

/**
 *  Copy constructor
 */

template<class T>
L::GeneratorMatrix2Copy<T>::GeneratorMatrix2Copy
   (const GeneratorMatrix2Copy<T> &gm)
: HeapAllocatedGeneratorMatrix2<T> (gm, true)
{
   std::copy (gm.getMatrix(), gm.getMatrix() + m*dim, c);
}

template<class T>
L::GeneratorMatrix2Copy<T>::GeneratorMatrix2Copy (const GeneratorMatrix2<T> &gm)
: HeapAllocatedGeneratorMatrix2<T> (gm, true)
{
   std::copy (gm.getMatrix(), gm.getMatrix() + m*dim, c);
}


/**
 *  Copy from Generator Matrix
 */

template<class T>
L::GeneratorMatrix2Copy<T>::GeneratorMatrix2Copy (const GeneratorMatrix &gm)
   : HeapAllocatedGeneratorMatrix2<T> (
       gm.getDimension(),
       gm.getM(),
       std::min (gm.getTotalPrec(), unsigned (std::numeric_limits<T>::digits)),
       false)
{
   if (gm.getBase() != 2)  throw GM_CopyBase (2, gm.getBase());
   allocate();
   for (unsigned d = 0; d < dim; ++d)  assign (gm, d, *this, d);
}


/**
 *  Constructor  (using GMCopy)
 */

template<class T>
L::GeneratorMatrix2Copy<T>::GeneratorMatrix2Copy (
   const GeneratorMatrix &gm, const GMCopy &c)
: HeapAllocatedGeneratorMatrix2<T>
   (c.getDimension (gm),
    c.getM         (gm),
    c.getTotalPrec (gm, std::numeric_limits<T>::digits),
    false)
{
   const int dd = c.getEqui();
   checkCopyDim (gm, *this, dd);   // dim is checked during assignment

   allocate();

   for (unsigned d = 0; d < dim - dd; ++d)  assign (gm, d, *this, d + dd);

   if (dim && dd)  makeEquidistributedCoordinate(0);
}


namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template class GeneratorMatrix2<X>; \
   template class MutableGeneratorMatrix2<X>; \
   template class HeapAllocatedGeneratorMatrix2<X>; \
   template class GeneratorMatrix2Copy<X>; \
   template void assign ( \
      const GeneratorMatrix &, unsigned, \
      MutableGeneratorMatrix2<X> &, unsigned); \
   template void assign ( \
      const GeneratorMatrix &, MutableGeneratorMatrix2<X> &);

   HINTLIB_INSTANTIATE (u32)
   #ifdef HINTLIB_U32_NOT_EQUAL_U64
   HINTLIB_INSTANTIATE (u64)
   #endif
#undef HINTLIB_INSTANTIATE
}


