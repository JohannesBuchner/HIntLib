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

#include <iomanip>
#include <algorithm>

#include <HIntLib/generatormatrix2.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/exception.h>

#ifdef HINTLIB_HAVE_OSTREAM
#include <ostream>
#else
#include <iostream>
#endif

namespace L = HIntLib;


/**
 *  assign Matrix
 *
 *  Assigns a submatrix of a Generator Matrix to another Generator Matrix
 */

template<typename T>
void L::assign (
   const GeneratorMatrix & src, int srcDim,
         GeneratorMatrix2<T> & dst, int dstDim)
{
#if 0
   cerr <<"A2" <<endl;
   cerr << "dim = " << src.getDimension() << " "<<dst.getDimension() << endl;
   cerr << "m = " << src.getM() << " "<<dst.getM() << endl;
   cerr << "prec = " << src.getPrec() << " "<<dst.getPrec() << endl;
   cerr << "base = " << src.getBase() << " "<<dst.getBase() << endl;
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

   for (int r = 0; r < dst.getM(); ++r)
   {
      T x = 0;

      for (int b = 0; b < dst.getPrec(); ++b)
      {
         x = (x << 1) | src.getDigit (srcDim, r, b);
      }

      dst.setv (dstDim, r, x);
   }
}

template<typename T>
void L::assign (const GeneratorMatrix &src, GeneratorMatrix2<T> &dst)
{
   for (int d = 0; d < dst.getDimension(); ++d)  assign (src, d, dst, d);
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
   if (prec > 32)      suffix = "ull";
   else if (prec > 16) suffix = "ul";
   else                     suffix = "u";

   o << "{\n";

   for (int r = 0; r < m; ++r)
   {
      o << "   {  // bit " << r << '\n';

      for (int d = 0; d < dim; ++d)
      {
         o << std::setw(22) << getVector(d,r) << suffix << ", // dim ="
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
 *  Constructor
 */

template<typename T>
L::GeneratorMatrix2<T>::GeneratorMatrix2 (int _dim)
   : GeneratorMatrix2Base (
         MAX_TOTALPREC, _dim, DEFAULT_M_BASE2, CORR_DEFAULT_TOTALPREC_BASE2)
{
   allocate();
}

template<typename T>
L::GeneratorMatrix2<T>::GeneratorMatrix2 (int _dim, int _m)
   : GeneratorMatrix2Base (
         MAX_TOTALPREC, _dim, _m,
            std::min (_m, int(CORR_DEFAULT_TOTALPREC_BASE2)))
{
   allocate();
}

template<typename T>
L::GeneratorMatrix2<T>::GeneratorMatrix2 (int _dim, int _m, int _prec)
   : GeneratorMatrix2Base (MAX_TOTALPREC, _dim, _m, _prec)
{
   checkPrec();
   allocate();
}


/**
 *  Copy constructor
 */

template<typename T>
L::GeneratorMatrix2<T>::GeneratorMatrix2 (const GeneratorMatrix2<T> &gm)
   : GeneratorMatrix2Base (gm)
{
   allocate();
   setMatrix (gm.getMatrix());
}

template<typename T>
L::GeneratorMatrix2<T>::GeneratorMatrix2 (const GeneratorMatrix &gm)
   : GeneratorMatrix2Base (
       MAX_TOTALPREC,
       gm.getDimension(),
       gm.getM(),
       std::min (gm.getPrec(), int(MAX_TOTALPREC)))
{
   if (gm.getBase() != 2)  throw GM_CopyBase (2, gm.getBase());
   allocate();
   for (int d = 0; d < dim; ++d)  assign (gm, d, *this, d);
}


/**
 *  allocate()
 */

template<typename T>
void L::GeneratorMatrix2<T>::allocate()
{
   c = new T [dim * m];
   makeZeroMatrix();
}


/**
 *  check Total Prec()
 */

template<typename T>
void L::GeneratorMatrix2<T>::checkPrec() const
{
   if (prec > MAX_TOTALPREC)
   {
      throw GM_PrecTooHigh (base, MAX_TOTALPREC, prec);
   }
}


/**
 *  set Matrix ()
 */

template<typename T>
void L::GeneratorMatrix2<T>::setMatrix (const T* source)
{
   std::copy (source, source + m * dim, c);
}


/**
 *  getPackedRowVector ()
 */

template<typename T>
L::u64
L::GeneratorMatrix2<T>::getPackedRowVector (int d, int b) const
{
   const T ma = mask (b);
   u64 result (0);

   for (int r = m - 1; r >= 0; --r)
   {
      result = (result << 1) | getdMask (d, r, ma);
   }

   return result;
}


/**
 *  setPackedRowVector ()
 */

template<typename T>
void
L::GeneratorMatrix2<T>::setPackedRowVector (int d, int b, u64 x)
{
   const T ma = mask (b);

   for (int r = 0; r < m; ++r)
   {
      setdMask (d, r, ma, unsigned(x) & 1u);
      x >>= 1;
   }
}


/**
 *  makeZeroRowVector (()
 */

template<typename T>
void L::GeneratorMatrix2<T>::makeZeroRowVector (int d, int b)
{
   const T ma = mask (b);

   for (int r = 0; r < m; ++r)  setd0Mask (d, r, ma);
}


/**
 *  getDigit()
 */

template<typename T>
int
L::GeneratorMatrix2<T>::getDigit (int d, int r, int b) const
{
   return (*this)(d, r, b);
}


/**
 *  setDigit ()
 */

template<typename T>
void L::GeneratorMatrix2<T>::setDigit (int d, int r, int b, int x)
{
   setd (d, r, b, x);
}


/**
 *  vSet/GetPackedRowVector ()
 */

template<typename T>
L::u64
L::GeneratorMatrix2<T>::vGetPackedRowVector (int d, int b) const
{
   return getPackedRowVector (d, b);
}

template<typename T>
void
L::GeneratorMatrix2<T>::vSetPackedRowVector (int d, int b, u64 x)
{
   setPackedRowVector (d, b, x);
}

template<typename T>
L::u64
L::GeneratorMatrix2<T>::getVector (int d, int b) const
{
   return u64 ((*this)(d,b));
}


/**
 *  adjust Total Prec()
 */

template<typename T>
void L::GeneratorMatrix2<T>::adjustPrec (int n)
{
   if (n != prec)
   {
      if (n < prec)
      {
         int shift = prec - n;
         for (T* p = c; p < c + m*dim; ++p)  *p >>= shift;
         prec = n;
      }
      else  // _vec > vec
      {
         int shift = n - prec;
         prec = n;
         checkPrec();
         for (T* p = c; p < c + m*dim; ++p)  *p <<= shift;
      }

   }
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
L::GeneratorMatrix2<T>::makeHammersley (int d)
{
   if (prec >= m)
   {
      T a = mask (m - 1);
      for (int r = 0; r < m; ++r)
      {
         setv (d, r, a);
         a <<= 1;
      }
   }
   else  // prec < m
   {
      for (int r = 0; r < m - prec; ++r)  makeZeroColumnVector (d,r);
      T a (1);
      for (int r = m - prec; r < m; ++r)  setv (d, r, a), a <<= 1;
   }
}


/**
 *  make Identity Matrix ()
 *
 *  Sets the matrix   d  to the identity matrix.
 *  Sets the matrix for dimensin  d  to the Halton sequence for 2.
 */

template<typename T>
void L::GeneratorMatrix2<T>::makeIdentityMatrix (int d)
{
   T a = mask (0);

   for (int r = 0; r < m; ++r)
   {
      setv (d,r, a);
      a >>= 1;
   }
}


/**
 *  make Zero Matrix ()
 *
 *  Fills the matrix for dimension  d  with 0.
 */

template<typename T>
void L::GeneratorMatrix2<T>::makeZeroMatrix (int d)
{
   for (int r = 0; r < m; ++r)  setv (d,r, 0);
}

template<typename T>
void L::GeneratorMatrix2<T>::makeZeroMatrix ()
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

template<typename T>
void L::GeneratorMatrix2<T>::prepareForGrayCode ()
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
void L::GeneratorMatrix2<T>::restoreFromGrayCode ()
{
   for (int d = 0; d < dim; ++d)
   {
      for (int r = m - 1; r > 0; --r)
      {
         (*this)(d,r) ^= (*this)(d,r-1);
      }
   }
}


/**
 *  make Shift Net ()
 *
 *  Creates a shift net based on the matrix for dimension 0
 */

template<typename T>
void L::GeneratorMatrix2<T>::makeShiftNet(int b)
{
   const T ma = mask (b);

   for (int r = 0; r < m; ++r)
   {
      if (getdMask (0, r, ma))
      {
         for (int d = 1;       d < dim - r; ++d) setd1Mask (d, d+r,    ma);
         for (int d = dim - r; d < dim;     ++d) setd1Mask (d, d+r-dim,ma);
      }
      else
      {
         for (int d = 1;       d < dim - r; ++d) setd0Mask (d, d+r,    ma);
         for (int d = dim - r; d < dim;     ++d) setd0Mask (d, d+r-dim,ma);
      }
   }
}

template<typename T>
void L::GeneratorMatrix2<T>::makeShiftNet()
{
   for (int r = 0; r < m; ++r)
   {
      T x = (*this)(0, r);

      for (int d = 1;       d < dim - r; ++d)  setv (d, d+r,     x);
      for (int d = dim - r; d < dim;     ++d)  setv (d, d+r-dim, x);
   }
}


/**
 *  operator==
 */

template<typename T>
bool L::operator== (
      const L::GeneratorMatrix2<T> &m1, const L::GeneratorMatrix2<T> &m2)
{
   return m1.getDimension() == m2.getDimension()
       && m1.getM()         == m2.getM()
       && m1.getPrec()      == m2.getPrec()
       && std::equal (m1.getMatrix(),
                      m1.getMatrix() + m1.getDimension()*m1.getM(),
                      m2.getMatrix());
}


namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template class GeneratorMatrix2<X>; \
   template void assign ( \
      const GeneratorMatrix &, int, GeneratorMatrix2<X> &, int); \
   template void assign (const GeneratorMatrix &, GeneratorMatrix2<X> &); \
   template bool operator== ( \
      const GeneratorMatrix2<X> &, const GeneratorMatrix2<X> &);

   HINTLIB_INSTANTIATE (u32)
#  ifdef HINTLIB_U32_NOT_EQUAL_U64
   HINTLIB_INSTANTIATE (u64)
#  endif
#undef HINTLIB_INSTANTIATE
}


