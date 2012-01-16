/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schürer <rudolf.schuerer@sbg.ac.at>
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

#define HINTLIB_LIBRARY_OBJECT

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
         GeneratorMatrixGen<T> & dst, unsigned dstDim)
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

   for (unsigned r = 0; r < dst.getM(); ++r)
   {
      for (unsigned b = 0; b < dst.getPrec(); ++b)
      {
         dst.setd (dstDim, r, b, src.getDigit (srcDim, r, b));
      }
   }
}

template<class T>
void L::assign (const GeneratorMatrix &src, GeneratorMatrixGen<T> &dst)
{
   for (unsigned d = 0; d < dst.getDimension(); ++d)  assign (src, d, dst, d);
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
L::GeneratorMatrixGen<T>::GeneratorMatrixGen (unsigned _base, unsigned _dim)
   : GeneratorMatrixGenBase (
      _base, _dim, getDefaultM (_base), getDefaultPrec (_base))
{
   checkBase();
   allocate();
   makeZeroMatrix();
}

template<class T>
L::GeneratorMatrixGen<T>::GeneratorMatrixGen
      (unsigned _base, unsigned _dim, unsigned _m)
   : GeneratorMatrixGenBase (
      _base, _dim, _m, std::min (_m, getDefaultPrec (_base)))
{
   checkBase();
   allocate();
   makeZeroMatrix();
}


template<class T>
L::GeneratorMatrixGen<T>::GeneratorMatrixGen
      (unsigned _base, unsigned _dim, unsigned _m, unsigned _prec)
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
   for (unsigned d = 0; d < dim; ++d)  assign (gm, d, *this, d);
}


/**
 *  checkBase()
 */

template<class T>
void L::GeneratorMatrixGen<T>::checkBase() const
{
   if (getBase() < 2 || getBase() - 1 > numeric_limits<T>::max())
   {
      throw GM_BaseTooLarge (getBase(), unsigned (numeric_limits<T>::max()));
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
void L::GeneratorMatrixGen<T>::makeZeroColumnVector (unsigned d, unsigned r)
{
   T* base = (*this)(d,r);
   std::fill (base, base + prec, 0);
}


/**
 *  getPackedRowVector ()
 */

template<class T>
L::u64
L::GeneratorMatrixGen<T>::getPackedRowVector (unsigned d, unsigned b) const
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
L::GeneratorMatrixGen<T>::setPackedRowVector (unsigned d, unsigned b, u64 x)
{
   for (unsigned r = 0; r < m; ++r)
   {
      setd (d, r, b, x % base);

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
void L::GeneratorMatrixGen<T>::makeZeroRowVector (unsigned d, unsigned b)
{
   for (unsigned r = 0; r < m; ++r)  setd (d, r, b, 0);
}


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

   return unsigned((*this) (d, r, b));
}


/**
 *  setDigit ()
 */

template<class T>
void L::GeneratorMatrixGen<T>::setDigit
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
 *  vSet/GetPackedRowVector ()
 */

template<typename T>
L::u64
L::GeneratorMatrixGen<T>::vGetPackedRowVector (unsigned d, unsigned b) const
{
   return getPackedRowVector (d, b);
}

template<typename T>
void
L::GeneratorMatrixGen<T>::vSetPackedRowVector (unsigned d, unsigned b, u64 x)
{
   setPackedRowVector (d, b, x);
}


/**
 *  make Hammersley()
 */

template<class T>
void L::GeneratorMatrixGen<T>::makeHammersley (unsigned d)
{
   makeZeroMatrix (d);

   for (unsigned r = 0; r < m; ++r)
   {
      if (m-(r+1) < prec)  setd (d, r, m-(r+1), 1);
   }
}


/**
 *  make Identity Matrix ()
 */

template<class T>
void L::GeneratorMatrixGen<T>::makeIdentityMatrix  (unsigned d)
{
   makeZeroMatrix (d);

   const unsigned ub = std::min (m, prec);
   for (unsigned r = 0; r < ub; ++r)  setd (d, r, r, 1);
}


/**
 *  make Zero Matrix ()
 */

template<class T>
void L::GeneratorMatrixGen<T>::makeZeroMatrix (unsigned d)
{
   for (unsigned r = 0; r < m; ++r)  makeZeroColumnVector (d, r);
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
void L::GeneratorMatrixGen<T>::makeShiftNet (unsigned b)
{
   for (unsigned r = 0; r < m; ++r)
   {
      T x = (*this)(0, r, b);

      for (unsigned d = 1;       d < dim - r; ++d)  setd (d, d+r,     b, x);
      for (unsigned d = dim - r; d < dim;     ++d)  setd (d, d+r-dim, b, x);
   }
}

template<typename T>
void L::GeneratorMatrixGen<T>::makeShiftNet()
{
   for (unsigned b = 0; b < prec; ++b)  makeShiftNet (b);
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
      const GeneratorMatrix &, unsigned, GeneratorMatrixGen<X> &, unsigned); \
   template void assign (const GeneratorMatrix &, GeneratorMatrixGen<X> &); \
   template bool operator== ( \
      const GeneratorMatrixGen<X> &, const GeneratorMatrixGen<X> &);

   HINTLIB_INSTANTIATE (unsigned char)
   HINTLIB_INSTANTIATE (unsigned short)
#undef HINTLIB_INSTANTIATE
}

