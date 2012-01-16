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

#include <HIntLib/generatormatrixvec.h>
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
   const GeneratorMatrix & src, unsigned srcDim,
         GeneratorMatrixVec<T> & dst, unsigned dstDim)
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

   if (dst.getVectorization() == 1 && src.getVectorization() == 1)
   {
      for (unsigned r = 0; r < dst.getM(); ++r)
      {
         for (unsigned b = 0; b < dst.getTotalPrec(); ++b)
         {
            dst.setv (dstDim, r, b, T(src.getVector (srcDim, r, b)));
         }
      }
   }
   else if (dst.getVectorization() == 1)
   {
      // FIXME Optimization for src.getVectorization() > 1 missing
      // cf. GeneratorMatrixGen::assign()

      for (unsigned r = 0; r < dst.getM(); ++r)
      {
         for (unsigned b = 0; b < dst.getTotalPrec(); ++b)
         {
            dst.setv (dstDim, r, b, src.getDigit (srcDim, r, b));
         }
      }
   }
   else  // dst use non-trivial vectorization
   {
      const unsigned base = dst.getBase();
      const unsigned vectorization = dst.getVectorization();
      const unsigned leadingDigits = dst.getNumOfLeadingDigits();

      for (unsigned r = 0; r < dst.getM(); ++r)
      {
         unsigned srcB = 0;

         // first vector (b == 0)
         {
            T x (0);
            for (unsigned v = 0; v < leadingDigits; ++v)
            {
               x = x * base + src.getDigit (srcDim, r, srcB++);
            }
            dst.setv (dstDim, r, 0, x);
         }

         for (unsigned b = 1; b < dst.getPrec(); ++b)
         {
            T x (0);
            for (unsigned v = 0; v < vectorization; ++v)
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
   for (unsigned d = 0; d < dst.getDimension(); ++d)  assign (src, d, dst, d);
}


/*****************************************************************************/
/***     Generator Matrix Vec <T>                                          ***/
/*****************************************************************************/

/**
 *  Constructor
 */

template<typename T>
L::GeneratorMatrixVec<T>::GeneratorMatrixVec
      (unsigned _base, unsigned _vec, unsigned _dim)
   : GeneratorMatrixVecBase (
         _base, _vec, _dim, getDefaultM (_base), getDefaultTotalPrec (_base))
{
   checkVecBase();
   allocate();
}

template<typename T>
L::GeneratorMatrixVec<T>::GeneratorMatrixVec
      (unsigned _base, unsigned _vec, unsigned _dim, unsigned _m)
   : GeneratorMatrixVecBase (
         _base, _vec, _dim, _m, std::min (_m, getDefaultTotalPrec (_base)))
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
L::GeneratorMatrixVec<T>::GeneratorMatrixVec (const GeneratorMatrix &gm)
   : GeneratorMatrixVecBase (
      gm.getBase(),
      digitsRepresentable (T(gm.getBase())),
      gm.getDimension(),
      gm.getM(),
      gm.getTotalPrec())
{
   checkVecBase();
   allocate();
   for (unsigned d = 0; d < dim; ++d)  assign (gm, d, *this, d);
}


/**
 *  Constructor  (truncating a given Generator Matrix)
 */

template<typename T>
L::GeneratorMatrixVec<T>::GeneratorMatrixVec (
   const GeneratorMatrixVec<T> &gm, const GMCopy &c)
: GeneratorMatrixVecBase
   (gm.getBase(),
    c.getVectorization (gm, std::numeric_limits<T>::digits),
    c.getDimension (gm),
    c.getM (gm),
    c.getTotalPrec (gm))
{
   checkVecBase();

   const int dd = c.getEqui();
   checkCopyDim (gm, *this, dd);

   allocate();

   for (unsigned d = 0; d < dim - dd; ++d)  assign (gm, d, *this, d + dd);

   if (dd && dim)  makeEquidistributedCoordinate (0);
}

template<typename T>
L::GeneratorMatrixVec<T>::GeneratorMatrixVec (
   const GeneratorMatrix &gm, const GMCopy& c)
: GeneratorMatrixVecBase
   (gm.getBase(),
    c.getVectorization(gm, std::numeric_limits<T>::digits),
    c.getDimension (gm),
    c.getM (gm),
    c.getTotalPrec (gm))
{
   checkVecBase();

   const int dd = c.getEqui();
   checkCopyDim (gm, *this, dd);

   allocate();

   for (unsigned d = 0; d < dim - dd; ++d)  assign (gm, d, *this, d + dd);

   if (dd && dim)  makeEquidistributedCoordinate (0);
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
   if (getVecBase()-1 > numeric_limits<T>::max())
   {
      throw GM_BaseTooLarge (getVecBase(), unsigned (numeric_limits<T>::max()));
   }
}


/**
 *  getd()
 */

template<typename T>
typename L::GeneratorMatrixVec<T>::D
L::GeneratorMatrixVec<T>::getd (unsigned d, unsigned r, unsigned b) const
{
   if (vec == 1) return operator()(d,r,b);

   b += getNumOfMissingDigits();
   return (operator()(d, r, b / vec) / powInt(base, vec - b % vec - 1)) % base;
}


/**
 *  setd ()
 */

template<typename T>
void L::GeneratorMatrixVec<T>::setd
   (unsigned d, unsigned r, unsigned b, typename GeneratorMatrixVec<T>::D x)
{
   if (vec == 1)  setv (d,r,b,x);
   else
   {
      b += getNumOfMissingDigits();

      T& vector = c[r*dimPrec + d*prec + b / vec];

      T shift = powInt(base, vec - b % vec - 1);

      T oldDigit = T((vector / shift) % base) * shift;
      vector = vector - oldDigit + x * shift;
      if (vector >= vecBase)  throw 1;
   }
}


/**
 *  getPackedRowVector ()
 */

template<typename T>
L::u64
L::GeneratorMatrixVec<T>::getPackedRowVector (unsigned d, unsigned b) const
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
L::GeneratorMatrixVec<T>::setPackedRowVector (unsigned d, unsigned b, u64 x)
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
 *  getDigit ()
 */

template<typename T>
unsigned
L::GeneratorMatrixVec<T>::getDigit (unsigned d, unsigned r, unsigned b) const
{
   if (numeric_limits<D>::digits > numeric_limits<unsigned>::digits)
   {
      throw InternalError (__FILE__, __LINE__);
   }

   return unsigned(getd(d, r, b));
}


/**
 *  getVector ()
 */

template<typename T>
L::u64
L::GeneratorMatrixVec<T>::getVector (unsigned d, unsigned r, unsigned b) const
{
   if (numeric_limits<T>::digits > numeric_limits<u64>::digits)
   {
      throw InternalError (__FILE__, __LINE__);
   }

   return u64(operator() (d, r, b));
}


/**
 *  setDigit ()
 */

template<typename T>
void L::GeneratorMatrixVec<T>::setDigit
   (unsigned d, unsigned r, unsigned b, unsigned x)
{
   if (   numeric_limits<typename GeneratorMatrixVec<T>::D>::digits
        < numeric_limits<unsigned>::digits
       && x >= getBase())
   {
      throw InternalError (__FILE__, __LINE__);
   }

   setd (d, r, b, typename GeneratorMatrixVec<T>::D(x));
}


/**
 *  setVector ()
 */

template<typename T>
void L::GeneratorMatrixVec<T>::setVector
   (unsigned d, unsigned r, unsigned b, u64 x)
{
   if (   numeric_limits<T>::digits < numeric_limits<u64>::digits
       && x >= getVecBase())
   {
      throw InternalError (__FILE__, __LINE__);
   }

   setv (d, r, b, T(x));
}


/**
 *  vSet/GetPackedRowVector ()
 */

template<typename T>
L::u64
L::GeneratorMatrixVec<T>::vGetPackedRowVector (unsigned d, unsigned b) const
{
   return getPackedRowVector (d, b);
}

template<typename T>
void
L::GeneratorMatrixVec<T>::vSetPackedRowVector (unsigned d, unsigned b, u64 x)
{
   setPackedRowVector (d, b, x);
}


/**
 *  make Equidistributed Coordinate()
 */

template<typename T>
void L::GeneratorMatrixVec<T>::makeEquidistributedCoordinate (unsigned d)
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

template<typename T>
void L::GeneratorMatrixVec<T>::makeIdentityMatrix  (unsigned d)
{
   makeZeroMatrix (d);
   
   const unsigned ub = std::min (m, totalPrec);
   for (unsigned r = 0; r < ub; ++r)  setd (d, r, r, 1);
}


/**
 *  make Zero Matrix ()
 */

template<typename T>
void L::GeneratorMatrixVec<T>::makeZeroMatrix (unsigned d)
{
   for (unsigned r = 0; r < m; ++r)
   {
      T* base = c + r*dimPrec + d*prec;
      std::fill (base, base + prec, 0);
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
   if (m1.getVectorization() == m2.getVectorization())
   {
      return m1.getDimension() == m2.getDimension()
          && m1.getM()         == m2.getM()
          && m1.getTotalPrec() == m2.getTotalPrec()
          && std::equal(
               m1.getMatrix(),
               m1.getMatrix()
                   + m1.getDimension() * m1.getTotalPrec() * m1.getM(),
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
      const GeneratorMatrix &, unsigned, GeneratorMatrixVec<X> &, unsigned); \
   template void assign (const GeneratorMatrix &, GeneratorMatrixVec<X> &); \
   template bool operator== ( \
      const GeneratorMatrixVec<X> &, const GeneratorMatrixVec<X> &);

   HINTLIB_INSTANTIATE (unsigned char)
   HINTLIB_INSTANTIATE (unsigned short)
   HINTLIB_INSTANTIATE (u32)
#undef HINTLIB_INSTANTIATE
}

