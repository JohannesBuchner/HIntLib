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

#include <algorithm>

#include <HIntLib/generatormatrix2row.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/exception.h>


namespace L = HIntLib;


/**
 *  assign()
 *
 *  Assigns a submatrix of a Generator Matrix to another Generator Matrix
 */

template<typename T>
void L::assign (
   const GeneratorMatrix & src, int srcDim,
         GeneratorMatrix2Row<T> & dst, int dstDim)
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

   if (src.getM() <= std::numeric_limits<u64>::digits)
   {
      const T mask = ~T(0) >> (dst.MAX_M - dst.getM());

      for (int b = 0; b < dst.getPrec(); ++b)
      {
         dst.setPackedRowVector (srcDim, b,
               T(src.vGetPackedRowVector (dstDim, b)) & mask);
      }
   }
   else
   {
      for (int b = 0; b < dst.getPrec(); ++b)
      {
         T x = 0;
         for (int r = dst.getM() - 1; r >= 0; --r)
         {
            x = (x << 1) | src.getDigit (srcDim, r, b);
         }
         dst.setPackedRowVector (dstDim, b, x);
      }
   }
}

template<typename T>
void L::assign (const GeneratorMatrix &src, GeneratorMatrix2Row<T> &dst)
{
   for (int d = 0; d < dst.getDimension(); ++d)  assign (src, d, dst, d);
}

/**
 *  Constructor
 */

template<typename T>
L::GeneratorMatrix2Row<T>::GeneratorMatrix2Row (int _dim)
   : GeneratorMatrix2RowBase (
      _dim, CORR_DEFAULT_M_BASE2, DEFAULT_TOTALPREC_BASE2)
{
   allocate();
}

template<typename T>
L::GeneratorMatrix2Row<T>::GeneratorMatrix2Row (int _dim, int _m)
   : GeneratorMatrix2RowBase (_dim, _m,
         std::min (_m, int(DEFAULT_TOTALPREC_BASE2)))
{
   checkM();
   allocate();
}

template<typename T>
L::GeneratorMatrix2Row<T>::GeneratorMatrix2Row (int _dim, int _m, int _prec)
   : GeneratorMatrix2RowBase (_dim, _m, _prec)
{
   checkM();
   allocate();
}


/**
 *  Copy constructor
 */

template<typename T>
L::GeneratorMatrix2Row<T>::GeneratorMatrix2Row (
      const GeneratorMatrix2Row<T> &gm)
   : GeneratorMatrix2RowBase (gm)
{
   allocate();
   setMatrix (gm.getMatrix());
}

template<typename T>
L::GeneratorMatrix2Row<T>::GeneratorMatrix2Row (const GeneratorMatrix &gm)
   : GeneratorMatrix2RowBase (
       gm.getDimension(),
       gm.getM(),
       gm.getPrec())
{
   checkM();
   if (gm.getBase() != 2)  throw GM_CopyBase (2, gm.getBase());
   allocate();
   for (int d = 0; d < dim; ++d)  assign (gm, d, *this, d);
}


/**
 *  allocate()
 */

template<typename T>
void L::GeneratorMatrix2Row<T>::allocate()
{
   c = new T [dim * prec];
   makeZeroMatrix();
}


/**
 *  checkM()
 */

template<typename T>
void L::GeneratorMatrix2Row<T>::checkM() const
{
   if (m > MAX_M)
   {
      // throw GM_PrecTooHigh (base, numeric_limits<T>::digits, m);
      throw FIXME (__FILE__, __LINE__);
   }
}


/**
 *  set Matrix ()
 */

template<typename T>
void L::GeneratorMatrix2Row<T>::setMatrix (const T* source)
{
   std::copy (source, source + prec * dim, c);
}


/**
 *  operator==
 */

template<typename T>
bool L::GeneratorMatrix2Row<T>::operator==
      (const L::GeneratorMatrix2Row<T> &gm) const
{
   return dim  == gm.dim
       && m    == gm.m
       && prec == gm.prec
       && std::equal(c, c + dim*prec , gm.c);
}


/**
 *  getDigit()
 */

template<typename T>
int
L::GeneratorMatrix2Row<T>::getDigit (int d, int r, int b) const
{
   return (*this) (d, r, b);
}


/**
 *  setDigit ()
 */

template<typename T>
void L::GeneratorMatrix2Row<T>::setDigit (int d, int r, int b, int x)
{
   setd (d, r, b, x);
}


/**
 *  vSet/GetPackedRowVector ()
 */

template<typename T>
L::u64
L::GeneratorMatrix2Row<T>::vGetPackedRowVector (int d, int b) const
{
   return getPackedRowVector (d, b);
}

template<typename T>
void
L::GeneratorMatrix2Row<T>::vSetPackedRowVector (int d, int b, u64 x)
{
   setPackedRowVector (d, b, T(x));
}


/**
 *  make Hammersley()
 *
 *  Sets the matrix for dimension  d  to such that an equidistributed
 *  coordinate is produiced.
 *  The Matrix is equal to a mirrored identity matrix.
 */

template<typename T>
void
L::GeneratorMatrix2Row<T>::makeHammersley (int d)
{
   if (prec >= m)
   {
      T a = mask (m - 1);
      for (int b = 0; b < m; ++b)
      {
         setPackedRowVector (d, b, a);
         a <<= 1;
      }
      for (int b = m; b < prec; ++b)  makeZeroRowVector (d, b);
   }
   else  // prec < m
   {
      T a = mask (m - 1);
      for (int b = 0; b < prec; ++b)
      {
         setPackedRowVector (d, b, a);
         a <<= 1;
      }
   }
}


/**
 *  make Identity Matrix ()
 *
 *  Sets the matrix   d  to the identity matrix.
 *  Sets the matrix for dimensin  d  to the Halton sequence for 2.
 */

template<typename T>
void L::GeneratorMatrix2Row<T>::makeIdentityMatrix (int d)
{
   T a = mask (0);
   int ub = std::min (m, prec);

   for (int b = 0; b < ub; ++b)
   {
      setPackedRowVector (d, b, a);
      a >>= 1;
   }
   for (int b = ub; b < prec; ++b)  makeZeroRowVector (d, b);
}


/**
 *  make Zero Matrix ()
 *
 *  Fills the matrix for dimension  d  with 0.
 */

template<typename T>
void L::GeneratorMatrix2Row<T>::makeZeroMatrix (int d)
{
   std::fill (c + d * prec, c + (d+1) * prec, 0);
}

template<typename T>
void L::GeneratorMatrix2Row<T>::makeZeroMatrix ()
{
   std::fill (c, c + dim * prec, 0);
}


#if 0
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

template<typename T>
void L::GeneratorMatrix2Row<T>::prepareForGrayCode ()
{
   for (int d = 0; d < dim; ++d)
   {
      T a = (*this)(d,0);

      for (int r = 1; r < m; ++r)  a = ((*this)(d,r) ^= a);
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

template<typename T>
void L::GeneratorMatrix2Row<T>::restoreFromGrayCode ()
{
   for (int d = 0; d < dim; ++d)
   {
      for (int r = m - 1; r > 0; --r)
      {
         (*this)(d,r) ^= (*this)(d,r-1);
      }
   }
}
#endif


/**
 *  make Shift Net ()
 *
 *  Creates a shift net based on the matrix for dimension 0
 */

template<typename T>
void L::GeneratorMatrix2Row<T>::makeShiftNet(int b)
{
   const T leadingDigit = mask (m - 1);
   const T allDigits = ~T(0) >> (MAX_M - m);
   T x = getPackedRowVector (0, b);

   for (int d = 1; d < dim; ++d)
   {
      setPackedRowVector (d, b,
            x = ((x << 1) & allDigits) | T((x & leadingDigit) != 0));
   }
}

template<typename T>
void L::GeneratorMatrix2Row<T>::makeShiftNet()
{
   for (int b = 0; b < prec; ++b)  makeShiftNet (b);
}


namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template class GeneratorMatrix2Row<X>; \
   template void assign ( \
      const GeneratorMatrix &, int, GeneratorMatrix2Row<X> &, int); \
   template void assign (const GeneratorMatrix &, GeneratorMatrix2Row<X> &);

   HINTLIB_INSTANTIATE (u32)
#  ifdef HINTLIB_U32_NOT_EQUAL_U64
   HINTLIB_INSTANTIATE (u64)
#  endif
#undef HINTLIB_INSTANTIATE
}


