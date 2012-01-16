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

#include <HIntLib/generatormatrixgen.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

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
   const GeneratorMatrix & src, int srcDim,
         GeneratorMatrixGen<T> & dst, int dstDim)
{
#if 0
   cerr <<"A3" <<endl;
   cerr << "dim = " << srcM.getDimension() << " "<<dstM.getDimension() << endl;
   cerr << "m = " << srcM.getM() << " "<<dstM.getM() << endl;
   cerr << "prec = " << srcM.getPrec() << " "<<dstM.getPrec() << endl;
   cerr << "base = " << srcM.getBase() << " "<<dstM.getBase() << endl;
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
      for (int b = 0; b < dst.getPrec(); ++b)
      {
         dst.setd (dstDim, r, b, src.getDigit (srcDim, r, b));
      }
   }
}

template<class T>
void L::assign (const GeneratorMatrix &src, GeneratorMatrixGen<T> &dst)
{
   for (int d = 0; d < dst.getDimension(); ++d)  assign (src, d, dst, d);
}

/*****************************************************************************/
/***     Generator Matrix Gen Base                                         ***/
/*****************************************************************************/

L::GeneratorMatrixGenBase::GeneratorMatrixGenBase (const GeneratorMatrix &gm)
   : GeneratorMatrix (
        gm.getBase(),
        gm.getDimension(),
        gm.getM(),
        gm.getPrec()),
     dimPrec (dim * prec)
{}


/*****************************************************************************/
/***     Generator Matrix Gen <T>                                          ***/
/*****************************************************************************/

/**
 *  Constructor
 */

template<class T>
L::GeneratorMatrixGen<T>::GeneratorMatrixGen (int _base, int _dim)
   : GeneratorMatrixGenBase (
      _base, _dim, getDefaultM (_base), getDefaultPrec (_base))
{
   checkBase();
   allocate();
   makeZeroMatrix();
}

template<class T>
L::GeneratorMatrixGen<T>::GeneratorMatrixGen (int _base, int _dim, int _m)
   : GeneratorMatrixGenBase (
      _base, _dim, _m, std::min (_m, getDefaultPrec (_base)))
{
   checkBase();
   allocate();
   makeZeroMatrix();
}


template<class T>
L::GeneratorMatrixGen<T>::GeneratorMatrixGen
      (int _base, int _dim, int _m, int _prec)
   : GeneratorMatrixGenBase (_base, _dim, _m, _prec)
{
   checkBase();
   allocate();
   makeZeroMatrix();
}


/**
 *  Copy constructor
 */

template<class T>
L::GeneratorMatrixGen<T>::GeneratorMatrixGen (const GeneratorMatrixGen<T> &gm)
   : GeneratorMatrixGenBase (gm)
{
   allocate();
   std::copy (gm.getMatrix(), gm.getMatrix() + m*dimPrec, c);
}


/**
 *  Copy from arbitrary Generator Matrix
 */

template<class T>
L::GeneratorMatrixGen<T>::GeneratorMatrixGen (const GeneratorMatrix &gm)
   : GeneratorMatrixGenBase (gm)
{
   checkBase();
   allocate();
   for (int d = 0; d < dim; ++d)  assign (gm, d, *this, d);
}


/**
 *  checkBase()
 */

template<class T>
void L::GeneratorMatrixGen<T>::checkBase() const
{
   if (getBase() < 2 || getBase() - 1 > numeric_limits<T>::max())
   {
      throw GM_BaseTooLarge (getBase(), numeric_limits<T>::max());
   }
}


/**
 *  allocate()
 */

template<class T>
void L::GeneratorMatrixGen<T>::allocate()
{
   c = new T [dimPrec * m];
}


/**
 *  makeZeroColumnVector ()
 */

template<class T>
void L::GeneratorMatrixGen<T>::makeZeroColumnVector (int d, int r)
{
   T* base = (*this)(d,r);
   std::fill (base, base + prec, 0);
}


/**
 *  getPackedRowVector ()
 */

template<class T>
L::u64
L::GeneratorMatrixGen<T>::getPackedRowVector (int d, int b) const
{
   u64 result (0);
   for (int r = m - 1; r >= 0; --r)  result = result*base + (*this)(d,r,b);
   return result;
}


/**
 *  setPackedRowVector ()
 */

template<class T>
void
L::GeneratorMatrixGen<T>::setPackedRowVector (int d, int b, u64 x)
{
   for (int r = 0; r < m; ++r)
   {
      setd (d, r, b, T(x % base));

      if (! (x /= base))
      {
         for (++r; r < m; ++r)  setd (d, r, b, 0);
         break;
      }
   }
}


/**
 *  makeZeroRowVector ()
 */

template<class T>
void L::GeneratorMatrixGen<T>::makeZeroRowVector (int d, int b)
{
   for (int r = 0; r < m; ++r)  setd (d, r, b, 0);
}


/**
 *  getDigit ()
 */

template<class T>
int
L::GeneratorMatrixGen<T>::getDigit (int d, int r, int b) const
{
   if (numeric_limits<T>::digits > numeric_limits<int>::digits)
   {
      throw InternalError (__FILE__, __LINE__);
   }

   return int((*this) (d, r, b));
}


/**
 *  setDigit ()
 */

template<class T>
void L::GeneratorMatrixGen<T>::setDigit (int d, int r, int b, int x)
{
   if (   numeric_limits<T>::digits < numeric_limits<int>::digits
       && x >= getBase())
   {
      throw InternalError (__FILE__, __LINE__);
   }

   setd (d, r, b, T(x));
}


/**
 *  vSet/GetPackedRowVector ()
 */

template<typename T>
L::u64
L::GeneratorMatrixGen<T>::vGetPackedRowVector (int d, int b) const
{
   return getPackedRowVector (d, b);
}

template<typename T>
void
L::GeneratorMatrixGen<T>::vSetPackedRowVector (int d, int b, u64 x)
{
   setPackedRowVector (d, b, x);
}


/**
 *  make Hammersley()
 */

template<class T>
void L::GeneratorMatrixGen<T>::makeHammersley (int d)
{
   makeZeroMatrix (d);

   for (int r = 0; r < m; ++r)
   {
      if (m-(r+1) < prec)  setd (d, r, m-(r+1), 1);
   }
}


/**
 *  make Identity Matrix ()
 */

template<class T>
void L::GeneratorMatrixGen<T>::makeIdentityMatrix  (int d)
{
   makeZeroMatrix (d);

   const int ub = std::min (m, prec);
   for (int r = 0; r < ub; ++r)  setd (d, r, r, 1);
}


/**
 *  make Zero Matrix ()
 */

template<class T>
void L::GeneratorMatrixGen<T>::makeZeroMatrix (int d)
{
   for (int r = 0; r < m; ++r)  makeZeroColumnVector (d, r);
}

template<class T>
void L::GeneratorMatrixGen<T>::makeZeroMatrix ()
{
   std::fill (c, c + dimPrec * m, 0);
}


/**
 *  make Shift Net ()
 *
 *  Creates a shift net based on the matrix for dimension 0
 */

template<typename T>
void L::GeneratorMatrixGen<T>::makeShiftNet (int b)
{
   for (int r = 0; r < m; ++r)
   {
      T x = (*this)(0, r, b);

      for (int d = 1;       d < dim - r; ++d)  setd (d, d+r,     b, x);
      for (int d = dim - r; d < dim;     ++d)  setd (d, d+r-dim, b, x);
   }
}

template<typename T>
void L::GeneratorMatrixGen<T>::makeShiftNet()
{
   for (int b = 0; b < prec; ++b)  makeShiftNet (b);
}


/**
 *  operator==
 */

template<class T>
bool L::operator==
   (const L::GeneratorMatrixGen<T> &m1, const L::GeneratorMatrixGen<T> &m2)
{
   return m1.getDimension() == m2.getDimension()
       && m1.getM()         == m2.getM()
       && m1.getPrec()      == m2.getPrec()
       && std::equal(
            m1.getMatrix(),
            m1.getMatrix() + m1.getDimension() * m1.getPrec() * m1.getM(),
            m2.getMatrix());
}


namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template class GeneratorMatrixGen<X>; \
   template void assign ( \
      const GeneratorMatrix &, int, GeneratorMatrixGen<X> &, int); \
   template void assign (const GeneratorMatrix &, GeneratorMatrixGen<X> &); \
   template bool operator== ( \
      const GeneratorMatrixGen<X> &, const GeneratorMatrixGen<X> &);

   HINTLIB_INSTANTIATE (unsigned char)
   HINTLIB_INSTANTIATE (unsigned short)
#undef HINTLIB_INSTANTIATE
}

