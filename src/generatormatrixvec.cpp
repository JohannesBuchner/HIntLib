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

#include <HIntLib/generatormatrixvec.h>

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

template<typename T>
void L::assign (
   const GeneratorMatrix & src, int srcDim,
         GeneratorMatrixVec<T> & dst, int dstDim)
{
#if 0
   cerr <<"A3" <<endl;
   cerr << "dim = " << srcM.getDimension() << " "<<dstM.getDimension() << endl;
   cerr << "m = " << srcM.getM() << " "<<dstM.getM() << endl;
   cerr << "prec = " << srcM.getPrec() << " "<<dstM.getPrec() << endl;
   cerr << "base = " << srcM.getBase() << " "<<dstM.getBase() << endl;
   cerr << "vecPrec = " << srcM.getVecPrec() << " "<<dstM.getVecPrec() << endl;
   cerr << "vec = " << srcM.getVec() << " "<<dstM.getVec() << endl;
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

   if (dst.getVec() == 1)
   {
      for (int r = 0; r < dst.getM(); ++r)
      {
         for (int b = 0; b < dst.getPrec(); ++b)
         {
            dst.setv (dstDim, r, b, src.getDigit (srcDim, r, b));
         }
      }
   }
   else  // dst use non-trivial vectorization
   {
      const int base = dst.getBase();
      const int vectorization = dst.getVec();
      const int leadingDigits = dst.getNumOfLeadingDigits();

      for (int r = 0; r < dst.getM(); ++r)
      {
         int srcB = 0;

         // first vector (b == 0)
         {
            T x (0);
            for (int v = 0; v < leadingDigits; ++v)
            {
               x = x * base + src.getDigit (srcDim, r, srcB++);
            }
            dst.setv (dstDim, r, 0, x);
         }

         for (int b = 1; b < dst.getVecPrec(); ++b)
         {
            T x (0);
            for (int v = 0; v < vectorization; ++v)
            {
               x = x * base + src.getDigit (srcDim, r, srcB++);
            }
            dst.setv (dstDim, r, b, x);
         }
      }
   }
}

template<typename T>
void L::assign (
      const GeneratorMatrix &src, GeneratorMatrixVec<T> &dst)
{
   for (int d = 0; d < dst.getDimension(); ++d)  assign (src, d, dst, d);
}


/*****************************************************************************/
/***     Generator Matrix Base                                             ***/
/*****************************************************************************/

inline
L::GeneratorMatrixVecBase::GeneratorMatrixVecBase
      (int _base, int _dim, int _m, int _prec, int _vec)
   : GeneratorMatrix (_base, _dim, _m, _prec),
     vec (_vec),
     vecBase (powInt (base, _vec)),
     vecPrec (div_ru (_prec, _vec)),
     dimPrec (dim*vecPrec)
{}

inline
L::GeneratorMatrixVecBase::GeneratorMatrixVecBase
      (const GeneratorMatrix &gm, int _vec)
   : GeneratorMatrix (gm),
     vec (_vec),
     vecBase (powInt (base, _vec)),
     vecPrec (div_ru (prec, _vec)),
     dimPrec (dim * vecPrec)
{}

inline
L::GeneratorMatrixVecBase::GeneratorMatrixVecBase
      (const GeneratorMatrixVecBase& gm)
   : GeneratorMatrix (gm),
     vec     (gm.vec),
     vecBase (gm.vecBase),
     vecPrec (gm.vecPrec),
     dimPrec (gm.dimPrec)
{}


/*****************************************************************************/
/***     Generator Matrix Vec <T>                                          ***/
/*****************************************************************************/

/**
 *  Constructor
 */

template<typename T>
L::GeneratorMatrixVec<T>::GeneratorMatrixVec (int _base, int _dim, int _vec)
   : GeneratorMatrixVecBase (
         _base, _dim, getDefaultM (_base), getDefaultPrec (_base), _vec)
{
   checkVecBase();
   allocate();
}

template<typename T>
L::GeneratorMatrixVec<T>::GeneratorMatrixVec
      (int _base, int _dim, int _m, int _vec)
   : GeneratorMatrixVecBase (
         _base, _dim, _m, std::min (_m, getDefaultPrec (_base)), _vec)
{
   checkVecBase();
   allocate();
}


/**
 *  Copy constructor
 */

template<typename T>
L::GeneratorMatrixVec<T>::GeneratorMatrixVec (const GeneratorMatrixVec<T> &gm)
   : GeneratorMatrixVecBase (gm)
{
   allocate();
   std::copy (gm.getMatrix(), gm.getMatrix() + m*dimPrec, c);
}


/**
 *  Copy from arbitrary Generator Matrix
 */

template<typename T>
L::GeneratorMatrixVec<T>::GeneratorMatrixVec
      (const GeneratorMatrix &gm)
   : GeneratorMatrixVecBase (
      gm.getBase(),
      gm.getDimension(),
      gm.getM(),
      gm.getPrec(),
      digitsRepresentable (T(gm.getBase())))
{
   checkVecBase();
   allocate();
   for (int d = 0; d < dim; ++d)  assign (gm, d, *this, d);
}


template<typename T>
L::GeneratorMatrixVec<T>::GeneratorMatrixVec
      (const GeneratorMatrix &gm, int _vec)
   : GeneratorMatrixVecBase (
      gm.getBase(),
      gm.getDimension(),
      gm.getM(),
      gm.getPrec(),
      _vec)
{
   checkVecBase();
   if (_vec < 1 || _vec > int(digitsRepresentable (T(gm.getBase()))))
   {
      throw FIXME (__FILE__, __LINE__);
   }

   allocate();
   for (int d = 0; d < dim; ++d)  assign (gm, d, *this, d);
}


/**
 *  allocate()
 */

template<typename T>
void L::GeneratorMatrixVec<T>::allocate()
{
   c = new T [dimPrec * m];
   makeZeroMatrix();
}


/**
 *  checkVecBase()
 */

template<typename T>
void L::GeneratorMatrixVec<T>::checkVecBase() const
{
   if (unsigned(getVecBase() - 1) > numeric_limits<T>::max())
   {
      throw GM_BaseTooLarge (getVecBase(), numeric_limits<T>::max());
   }
}


/**
 *  getd()
 */

template<typename T>
typename L::GeneratorMatrixVec<T>::D
L::GeneratorMatrixVec<T>::getd (int d, int r, int b) const
{
   if (vec == 1) return (*this)(d,r,b);

   b += getNumOfMissingDigits();
   return ((*this)(d, r, b / vec) / powInt(base, vec - b % vec - 1)) % base;
}


/**
 *  setd ()
 */

template<typename T>
void L::GeneratorMatrixVec<T>::setd
   (int d, int r, int b, typename GeneratorMatrixVec<T>::D x)
{
   if (vec == 1)  setv (d,r,b,x);
   else
   {
      b += getNumOfMissingDigits();

      T& vector = c[r*dimPrec + d*vecPrec + b / vec];

      T shift = powInt(base, vec - b % vec - 1);

      T oldDigit = T((vector / shift) % base) * shift;
      vector = vector - oldDigit + x * shift;
      if (vector >= unsigned(vecBase))
      {
         throw InternalError (__FILE__, __LINE__);
      }
   }
}


/**
 *  getPackedRowVector ()
 */

template<typename T>
L::u64
L::GeneratorMatrixVec<T>::getPackedRowVector (int d, int b) const
{
   u64 result (0);
   for (int r = m - 1; r >= 0; --r)  result = result*base + getd(d,r,b);
   return result;
}


/**
 *  setPackedRowVector ()
 */

template<typename T>
void
L::GeneratorMatrixVec<T>::setPackedRowVector (int d, int b, u64 x)
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
 *  getDigit ()
 */

template<typename T>
int
L::GeneratorMatrixVec<T>::getDigit (int d, int r, int b) const
{
   if (numeric_limits<D>::digits > numeric_limits<int>::digits)
   {
      throw InternalError (__FILE__, __LINE__);
   }

   return int(getd(d, r, b));
}


/**
 *  setDigit ()
 */

template<typename T>
void L::GeneratorMatrixVec<T>::setDigit (int d, int r, int b, int x)
{
   if (   numeric_limits<typename GeneratorMatrixVec<T>::D>::digits
        < numeric_limits<int>::digits
       && x >= getBase())
   {
      throw InternalError (__FILE__, __LINE__);
   }

   setd (d, r, b, typename GeneratorMatrixVec<T>::D(x));
}


/**
 *  vSet/GetPackedRowVector ()
 */

template<typename T>
L::u64
L::GeneratorMatrixVec<T>::vGetPackedRowVector (int d, int b) const
{
   return getPackedRowVector (d, b);
}

template<typename T>
void
L::GeneratorMatrixVec<T>::vSetPackedRowVector (int d, int b, u64 x)
{
   setPackedRowVector (d, b, x);
}


/**
 *  make Hammersley()
 */

template<typename T>
void L::GeneratorMatrixVec<T>::makeHammersley (int d)
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

template<typename T>
void L::GeneratorMatrixVec<T>::makeIdentityMatrix  (int d)
{
   makeZeroMatrix (d);

   const int ub = std::min (m, prec);
   for (int r = 0; r < ub; ++r)  setd (d, r, r, 1);
}


/**
 *  make Zero Matrix ()
 */

template<typename T>
void L::GeneratorMatrixVec<T>::makeZeroMatrix (int d)
{
   for (int r = 0; r < m; ++r)
   {
      T* base = c + r*dimPrec + d*vecPrec;
      std::fill (base, base + vecPrec, 0);
   }
}

template<typename T>
void L::GeneratorMatrixVec<T>::makeZeroMatrix ()
{
   std::fill (c, c + dimPrec * m, T(0));
}


/**
 *  operator==
 */

template<typename T>
bool L::operator==
   (const L::GeneratorMatrixVec<T> &m1, const L::GeneratorMatrixVec<T> &m2)
{
   if (m1.getVec() == m2.getVec())
   {
      return m1.getDimension() == m2.getDimension()
          && m1.getM()         == m2.getM()
          && m1.getPrec()      == m2.getPrec()
          && std::equal(
               m1.getMatrix(),
               m1.getMatrix()
                   + m1.getDimension() * m1.getPrec() * m1.getM(),
               m2.getMatrix());
   }

   return static_cast<const GeneratorMatrix&>(m1)
       == static_cast<const GeneratorMatrix&>(m2);
}


namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template class GeneratorMatrixVec<X>; \
   template void assign ( \
      const GeneratorMatrix &, int, GeneratorMatrixVec<X> &, int); \
   template void assign (const GeneratorMatrix &, GeneratorMatrixVec<X> &); \
   template bool operator== ( \
      const GeneratorMatrixVec<X> &, const GeneratorMatrixVec<X> &);

   HINTLIB_INSTANTIATE (unsigned char)
   HINTLIB_INSTANTIATE (unsigned short)
   HINTLIB_INSTANTIATE (u32)
#undef HINTLIB_INSTANTIATE
}

