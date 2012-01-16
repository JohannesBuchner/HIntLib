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
   (const GeneratorMatrix& m, unsigned bits) const
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
 *  Constructor
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
   o << "Base: " << getBase()
     << " Dim: " << getDimension()
     << " M: " << getM()
     << " TotalPrec: " << getTotalPrec()
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
 *  operator=
 */

L::GeneratorMatrix& L::GeneratorMatrix::operator= (const GeneratorMatrix &m)
{
   for (unsigned d = 0; d < getDimension(); ++d)  assign (m, d, *this, d);
   
   return *this;
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


namespace
{
   void checkCopy (const L::GeneratorMatrix &o, const L::GeneratorMatrix &n)
   {
      if (o.getBase() != n.getBase())
         throw L::GM_CopyBase (n.getBase(), o.getBase());

      if (o.getTotalPrec() < n.getTotalPrec())
         throw L::GM_CopyPrec (n.getTotalPrec(), o.getTotalPrec());

      if (o.getM() < n.getM())
         throw L::GM_CopyM (n.getM(), o.getM());
   }

   void checkCopyDim (
      const L::GeneratorMatrix &o,
      const L::GeneratorMatrix &n, int offset = 0)
   {
      checkCopy (o, n);

      if (int (o.getDimension()) + offset < int (n.getDimension()))
         throw L::GM_CopyDim (n.getDimension(), o.getDimension());
   }
}


/**
 *  assign Matrix
 *
 *  Assigns a submatrix of a Generator Matrix to another Generator Matrix
 */

void L::assign (
   const GeneratorMatrix & srcM, unsigned srcDim,
         GeneratorMatrix & dstM, unsigned dstDim)
{
   checkCopy (srcM, dstM);

   if (dstDim >= dstM.getDimension())
   {
      throw DimensionTooHigh (dstDim, dstM.getDimension() - 1);
   }

   if (srcDim >= srcM.getDimension())
   {
      throw DimensionTooHigh (dstDim, dstM.getDimension() - 1);
   }

   for (unsigned r = 0; r < dstM.getM(); ++r)
   {
      for (unsigned b = 0; b < dstM.getTotalPrec(); ++b)
      {
         dstM.setDigit (dstDim, r, b, srcM.getDigit (srcDim, r, b));
      }
   }
}

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

   checkCopy (src, dst);

   if (dstDim >= dst.getDimension())
   {
      throw DimensionTooHigh (dstDim, dst.getDimension() - 1);
   }

   if (srcDim >= src.getDimension())
   {
      throw DimensionTooHigh (dstDim, dst.getDimension() - 1);
   }

   const unsigned srcVectorization = src.getVectorization();
   unsigned leadingDigits = src.getTotalPrec() % srcVectorization;
   if (! leadingDigits)  leadingDigits = srcVectorization;

   if (leadingDigits >= dst.getTotalPrec())
   {
      const unsigned shift = leadingDigits - dst.getTotalPrec();

      for (unsigned r = 0; r < dst.getM(); ++r)
      {
         dst.setv (dstDim, r, src.getVector (srcDim, r, 0) >> shift);
      }
   }
   else
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

   checkCopy (src, dst);

   if (dstDim >= dst.getDimension())
   {
      throw DimensionTooHigh (dstDim, dst.getDimension() - 1);
   }

   if (srcDim >= src.getDimension())
   {
      throw DimensionTooHigh (dstDim, dst.getDimension() - 1);
   }

   if (dst.getVectorization() == 1)
   {
      for (unsigned r = 0; r < dst.getM(); ++r)
      {
         for (unsigned b = 0; b < dst.getPrec(); ++b)
         {
            dst.setv (dstDim, r, b, src.getDigit (srcDim, r, b));
         }
      }
   }
   else
   {
      for (unsigned r = 0; r < dst.getM(); ++r)
      {
         for (unsigned b = 0; b < dst.getTotalPrec(); ++b)
         {
            dst.setd (dstDim, r, b, src.getDigit (srcDim, r, b));
         }
      }
   }
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
 *  set ()
 */

template<class T>
void L::MutableGeneratorMatrix2<T>::set
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
   set (d, r, b, x);
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
   (unsigned _dim, unsigned _m, unsigned _totalPrec)
   : MutableGeneratorMatrix2<T> (_dim, _m, _totalPrec)
{
   allocate();
   std::fill (c, c + dim * m, 0);
}


/**
 *  allocate()
 */

template<class T>
void L::HeapAllocatedGeneratorMatrix2<T>::allocate()
{
   delete c;
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
   if (gm.getBase() != 2)  throw FIXME (__FILE__, __LINE__);
   allocate();

   for (unsigned d = 0; d < dim; ++d)  assign (gm, d, *this, d);
}


/**
 *  Constructor  (truncating a given Generator Matrix)
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
   checkCopyDim (gm, *this, c.getEqui());   // dim is checked during assignment

   allocate();

   const int dd = c.getEqui();
   for (unsigned d = 0; d < dim - dd; ++d)  assign (gm, d, *this, d + dd);

   if (c.getEqui())  makeEquidistributedCoordinate(0);
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
      MutableGeneratorMatrix2<X> &, unsigned dstDim);

   HINTLIB_INSTANTIATE (u32)
   #ifdef HINTLIB_U32_NOT_EQUAL_U64
   HINTLIB_INSTANTIATE (u64)
   #endif
#undef HINTLIB_INSTANTIATE
}


/*****************  General Base  ********************************************/


/*****************************************************************************/
/***     Generator Matrix Gen <T>                                          ***/
/*****************************************************************************/

/**
 *  getd()
 */

template<class T>
typename L::GeneratorMatrixGen<T>::D
L::GeneratorMatrixGen<T>::getd (unsigned d, unsigned r, unsigned b) const
{
   if (vec == 1) return operator()(d,r,b);

   b += prec * vec - totalPrec;
   return (operator()(d, r, b / vec) / powInt(base, vec - b % vec - 1)) % base;
}


/**
 *  getDigit ()
 */

template<class T>
unsigned
L::GeneratorMatrixGen<T>::getDigit (unsigned d, unsigned r, unsigned b) const
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

   if (vec == gm.vec)  return std::equal(c, c + dimPrec*m , gm.c);

   return static_cast<const GeneratorMatrix&>(*this)
       == static_cast<const GeneratorMatrix&>(gm);
}


/**
 *  checkVecBase()
 */

template<class T>
void L::GeneratorMatrixGen<T>::checkVecBase() const
{
   if (getVecBase()-1 > numeric_limits<T>::max())
   {
      throw GM_BaseTooLarge (getVecBase(), unsigned (numeric_limits<T>::max()));
   }
}


/*****************************************************************************/
/***     Mutable Generator Matrix Gen <T>                                  ***/
/*****************************************************************************/


/**
 *  setd ()
 */

template<class T>
void L::MutableGeneratorMatrixGen<T>::setd
   (unsigned d, unsigned r, unsigned b, typename GeneratorMatrixGen<T>::D x)
{
   if (vec == 1)
   {
      setv (d,r,b,x);
   }
   else
   {
      b += prec * vec - totalPrec;

      T& vector = c[r*dimPrec + d*prec + b / vec];

      T shift = powInt(base, vec - b % vec - 1);

      T oldDigit = T((vector / shift) % base) * shift;
      vector = vector - oldDigit + x * shift;
   }
}


/**
 *  setDigit ()
 */

template<class T>
void L::MutableGeneratorMatrixGen<T>::setDigit
   (unsigned d, unsigned r, unsigned b, unsigned x)
{
   if (   numeric_limits<typename GeneratorMatrixGen<T>::D>::digits
        < numeric_limits<unsigned>::digits)
   {
      throw InternalError (__FILE__, __LINE__);
   }

   setd (d, r, b, typename GeneratorMatrixGen<T>::D(x));
}


/**
 *  setVector ()
 */

template<class T>
void L::MutableGeneratorMatrixGen<T>::setVector
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
   
   for (unsigned r = 0; r < m; ++r)  if (r < totalPrec)  setd (d, r, r, 1);
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


/*****************************************************************************/
/***     Heap Allocated  Matrix Gen <T>                                    ***/
/*****************************************************************************/

/**
 *  Constructor
 */

template<class T>
L::HeapAllocatedGeneratorMatrixGen<T>::HeapAllocatedGeneratorMatrixGen
   (unsigned _base, unsigned _vec,
    unsigned _dim, unsigned _m, unsigned _totalPrec)
   : MutableGeneratorMatrixGen<T> (_base, _vec, _dim, _m, _totalPrec)
{
   allocate();
   std::fill (c, c + dimPrec * m, 0);
}

template<class T>
L::HeapAllocatedGeneratorMatrixGen<T>::HeapAllocatedGeneratorMatrixGen
   (unsigned _base, unsigned _vec, unsigned _dim)
   : MutableGeneratorMatrixGen<T> (
         _base, _vec, _dim,
         logInt (std::numeric_limits<Index>::max(), Index (_base)),
         unsigned (ceil(log(2.0) / log(double(_base))
                        * double(numeric_limits<real>::digits - 1))))
{
   allocate();
   std::fill (c, c + dimPrec * m, 0);
}


/**
 *  allocate()
 */

template<class T>
void L::HeapAllocatedGeneratorMatrixGen<T>::allocate()
{
   delete c;
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
      1,
      gm.getDimension(),
      gm.getM(),
      gm.getTotalPrec(),
      true)
{
   for (unsigned d = 0; d < dim; ++d)  assign (gm, d, *this, d);
}


/**
 *  Constructor  (truncating a given Generator Matrix)
 */

template<class T>
L::GeneratorMatrixGenCopy<T>::GeneratorMatrixGenCopy (
   const GeneratorMatrixGen<T> &gm, const GMCopy &c)
: HeapAllocatedGeneratorMatrixGen<T>
   (gm.getBase(),
    c.getVectorization (gm, std::numeric_limits<T>::digits),
    c.getDimension (gm),
    c.getM (gm),
    c.getTotalPrec (gm),
    false)
{
   checkCopyDim (gm, *this);

   allocate();

   const int dd = c.getEqui();
   for (unsigned d = 0; d < dim - dd; ++d)  assign (gm, d, *this, d + dd);

   if (c.getEqui())  makeEquidistributedCoordinate (0);
}

template<class T>
L::GeneratorMatrixGenCopy<T>::GeneratorMatrixGenCopy (
   const GeneratorMatrix &gm, const GMCopy& c)
: HeapAllocatedGeneratorMatrixGen<T>
   (gm.getBase(),
    c.getVectorization(gm, std::numeric_limits<T>::digits),
    c.getDimension (gm),
    c.getM (gm),
    c.getTotalPrec (gm),
    false)
{
   checkCopyDim (gm, *this);

   allocate();

   const int dd = c.getEqui();
   for (unsigned d = 0; d < dim - dd; ++d)  assign (gm, d, *this, d + dd);

   if (c.getEqui())  makeEquidistributedCoordinate (0);
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
      MutableGeneratorMatrixGen<X> &, unsigned);

   HINTLIB_INSTANTIATE (unsigned char)
   HINTLIB_INSTANTIATE (unsigned short)
   HINTLIB_INSTANTIATE (u32)
#undef HINTLIB_INSTANTIATE
}

