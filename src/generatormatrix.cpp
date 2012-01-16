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
/***     Generator Matrix                                                  ***/
/*****************************************************************************/

/**
 *  Constructor
 */

L::GeneratorMatrix::GeneratorMatrix (const GeneratorMatrix &gm)
: base(gm.base),
  dim (gm.dim),
  m   (gm.m),
  prec(gm.prec)
{}


/**
 *  setDigit()
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


/**
 *  dump()
 */

void L::GeneratorMatrix::dump (std::ostream &o) const
{
   unsigned size = logInt (base, 10u) + 1;
   if (base >= 10)  ++size;
   
   for (unsigned d = 0; d < dim; ++d)
   {
      o << "Dimension " << d << ":\n";

      for (unsigned r = 0; r < m; ++r)
      {
         for (unsigned p = 0; p < prec; ++p)
         {
            o << std::setw (size) << getDigit (d,r,p);
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
   checkDimensionEqual (m.getDimension(), getDimension());

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
       || gm1.getPrecision() != gm2.getPrecision()
       || gm1.getBase() != gm2.getBase())
   {
      return false;
   }

   for (unsigned d = 0; d < gm1.getDimension(); ++d)
   {
      for (unsigned r = 0; r < gm1.getM(); ++r)
      {
         for (unsigned b = 0; b < gm1.getPrecision(); ++b)
         {
            if (gm1.getDigit (d, r, b) != gm2.getDigit (d, r, b))  return false;
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

      if (o.getPrecision() < n.getPrecision())
         throw L::GM_CopyPrec (n.getPrecision(), o.getPrecision());

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
      for (unsigned b = 0; b < dstM.getPrecision(); ++b)
      {
         dstM.setDigit (dstDim, r, b, srcM.getDigit (srcDim, r, b));
      }
   }
}

template<class T>
void L::assign (
   const GeneratorMatrix & srcM, unsigned srcDim,
         MutableGeneratorMatrix2<T> & dstM, unsigned dstDim)
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
      T x (0);

      for (unsigned b = 0; b < dstM.getPrecision(); ++b)
      {
         x = 2 * x + srcM.getDigit (srcDim, r, b);
      }

      dstM.set (dstDim, r, x);
   }
}

template<class T>
void L::assign (
   const GeneratorMatrix & srcM, unsigned srcDim,
         MutableGeneratorMatrixGen<T> & dstM, unsigned dstDim)
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
      for (unsigned b = 0; b < dstM.getPrecision(); ++b)
      {
         dstM.set (dstDim, r, b, srcM.getDigit (srcDim, r, b));
      }
   }
}


/*****************************************************************************/
/***     Generator Matrix 2 Base                                           ***/
/*****************************************************************************/

/**
 *  binaryDump()
 */

void L::GeneratorMatrix2Base::binaryDump (std::ostream &o) const
{
   for (unsigned d = 0; d < dim; ++d)
   {
      o << "Dimension " << d << ":\n";

      for (unsigned r = 0; r < m; ++r)  o << getVector (d,r) << '\n';
   }
}


/**
 *  dumpAsCArray()
 */

void L::GeneratorMatrix2Base::dumpAsCArray (std::ostream &o) const
{
   const char* suffix;
   if (prec > 32)       suffix = "ull";
   else if (prec > 16)  suffix = "ul";
   else                 suffix = "u";
   
   o << "{\n";

   for (unsigned r = 0; r < m; ++r)
   {
      o << "   {  // bit " << r << '\n';

      for (unsigned d = 0; d < dim; ++d)
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
 *  checkPrec()
 */

template<class T>
void L::GeneratorMatrix2<T>::checkPrec() const
{
   if (prec > unsigned (numeric_limits<T>::digits))
      throw GM_PrecTooHigh (2, numeric_limits<T>::digits, prec);
}

/**
 *  operator==
 */

template<class T>
bool L::GeneratorMatrix2<T>::operator== (const L::GeneratorMatrix2<T> &gm) const
{
   return dim  == gm.dim
       && m    == gm.m
       && prec == gm.prec
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
L::u64 L::GeneratorMatrix2<T>::getVector (unsigned d, unsigned r) const
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
 *  setDigit ()
 */

template<class T>
void L::MutableGeneratorMatrix2<T>::setDigit
   (unsigned d, unsigned r, unsigned b, unsigned x)
{
   if (x > 1)  throw InternalError (__FILE__, __LINE__);

   set (d, r, b, x);
}


/**
 *  adjustPrecision()
 */

template<class T>
void L::MutableGeneratorMatrix2<T>::adjustPrecision (unsigned _prec)
{
   if (_prec == prec)  return;
   else
   {
      if (_prec < prec)
      {
         for (T* p = c; p < c + m*dim; ++p)  *p >>= (prec - _prec);
         prec = _prec;
      }
      else  // _prec > prec
      {
         unsigned shift = _prec - prec;
         prec = _prec;
         checkPrec();
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
   if (prec >= m)
   {
      T a = T(1) << (prec - m);
      for (unsigned r = 0; r < m; ++r)  c[r*dim + d] = a, a <<= 1;
   }
   else
   {
      for (unsigned r = 0; r < m - prec; ++r)  c[r*dim + d] = 0;
      T a (1);
      for (unsigned r = m - prec; r < m; ++r)  c[r*dim + d] = a, a <<= 1;
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
   T a = T(1) << (prec - 1);

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
   (unsigned _dim, unsigned _m, unsigned _prec)
   : MutableGeneratorMatrix2<T> (_dim,_m,_prec)
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
: HeapAllocatedGeneratorMatrix2<T> (gm.getDimension(), gm.getM(),
                                    gm.getPrecision(), false)
{
   if (gm.getBase() != 2)  throw GM_CopyBase (2, gm.getBase());

   allocate();

   for (unsigned d = 0; d < dim; ++d)  assign (gm, d, *this, d);
}

/**
 *  Constructor  (truncating a given Generator Matrix)
 */

template<class T>
L::GeneratorMatrix2Copy<T>::GeneratorMatrix2Copy (
   const GeneratorMatrix &gm,
   unsigned _dim, unsigned _m, unsigned _prec, bool equi)
: HeapAllocatedGeneratorMatrix2<T>
   (_dim  ? _dim : gm.getDimension() + equi,
    _m    ? _m   : gm.getM(),
    _prec ? _prec : gm.getPrecision(),
    false)
{
   checkCopyDim (gm, *this, equi);   // dim is checked during assignment

   allocate();

   const int dd = equi;
   for (unsigned d = 0; d < dim - dd; ++d)  assign (gm, d, *this, d + dd);

   if (equi)  makeEquidistributedCoordinate(0);
}

template<class T>
L::GeneratorMatrix2Copy<T>::GeneratorMatrix2Copy (
   const GeneratorMatrix2<T> &gm,
   unsigned _dim, unsigned _m, unsigned _prec, bool equi)
: HeapAllocatedGeneratorMatrix2<T>
   (_dim  ? _dim  : gm.getDimension() + equi,
    _m    ? _m    : gm.getM(),
    _prec ? _prec : gm.getPrecision(),
    false)
{
   checkCopyDim (gm, *this, equi);

   allocate();

   const int shift = gm.getPrecision() - prec;
   const int dd = equi;

   for (unsigned d = 0; d < dim - dd; ++d)
   {
      for (unsigned r = 0; r < m; ++r)  c[r*dim + d + dd] = gm(d,r) >> shift; 
   }

   if (equi)  makeEquidistributedCoordinate(0);
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
 *  getDigit ()
 */

template<class T>
unsigned L::GeneratorMatrixGen<T>::getDigit
   (unsigned d, unsigned r, unsigned b) const
{
   if (numeric_limits<T>::digits <= numeric_limits<unsigned>::digits)
   {
      return unsigned (operator() (d, r, b));
   }
   else  throw InternalError (__FILE__, __LINE__);
}


/**
 *  operator==
 */

template<class T>
bool L::GeneratorMatrixGen<T>::operator==
   (const L::GeneratorMatrixGen<T> &gm) const
{
   return dim  == gm.dim
       && m    == gm.m
       && prec == gm.prec
       && base == gm.base
       && std::equal(c, c + dimPrec*m , gm.c);
}


/**
 *  checkBase()
 */

template<class T>
void L::GeneratorMatrixGen<T>::checkBase() const
{
   if (base-1 > numeric_limits<T>::max())
   {
      throw GM_BaseTooLarge (base, unsigned (numeric_limits<T>::max()));
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
       && x >= base)
   {
      throw InternalError (__FILE__, __LINE__);
   }

   set (d, r, b, T(x));
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
      if (m-(r+1) < prec)  set (d, r, m-(r+1), 1);
   }
}


/**
 *  make Identity Matrix ()
 */

template<class T>
void L::MutableGeneratorMatrixGen<T>::makeIdentityMatrix  (unsigned d)
{
   makeZeroMatrix (d);
   
   for (unsigned r = 0; r < m; ++r)  if (r < prec)  set (d, r, r, 1);
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
   (unsigned _base, unsigned _dim, unsigned _m, unsigned _prec)
   : MutableGeneratorMatrixGen<T> (_base, _dim,_m,_prec)
{
   allocate();
   std::fill (c, c + dimPrec * m, 0);
}

template<class T>
L::HeapAllocatedGeneratorMatrixGen<T>::HeapAllocatedGeneratorMatrixGen
   (unsigned _base, unsigned _dim)
   : MutableGeneratorMatrixGen<T> (
         _base, _dim,
         logInt (std::numeric_limits<Index>::max(), Index (_base)),
         unsigned (  log(2.0) * double(numeric_limits<real>::digits - 1)
                   / log(double(_base))))
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
: HeapAllocatedGeneratorMatrixGen<T> (gm, true)
{
   for (unsigned d = 0; d < dim; ++d)  assign (gm, d, *this, d);
}


/**
 *  Constructor  (truncating a given Generator Matrix)
 */

template<class T>
L::GeneratorMatrixGenCopy<T>::GeneratorMatrixGenCopy (
   const GeneratorMatrixGen<T> &gm,
   unsigned _dim, unsigned _m, unsigned _prec, bool equi)
: HeapAllocatedGeneratorMatrixGen<T>
   (gm.getBase(),
    _dim  ? _dim  : gm.getDimension(),
    _m    ? _m    : gm.getM(),
    _prec ? _prec : gm.getPrecision(),
    false)
{
   checkCopyDim (gm, *this);

   allocate();

   const int dd = equi;
   for (unsigned d = 0; d < dim - dd; ++d)  assign (gm, d, *this, d + dd);

   if (equi)  makeEquidistributedCoordinate (0);
}

template<class T>
L::GeneratorMatrixGenCopy<T>::GeneratorMatrixGenCopy (
   const GeneratorMatrix &gm,
   unsigned _dim, unsigned _m, unsigned _prec, bool equi)
: HeapAllocatedGeneratorMatrixGen<T>
   (gm.getBase(),
    _dim  ? _dim  : gm.getDimension(),
    _m    ? _m    : gm.getM(),
    _prec ? _prec : gm.getPrecision(),
    false)
{
   checkCopyDim (gm, *this);

   allocate();

   const int dd = equi;
   for (unsigned d = 0; d < dim - dd; ++d)  assign (gm, d, *this, d + dd);

   if (equi)  makeEquidistributedCoordinate (0);
}


/*****************************************************************************/
/***     Generator Matrix Gen Vectorize <T>                                ***/
/*****************************************************************************/

template<class T>
L::GeneratorMatrixGenVectorize<T>::GeneratorMatrixGenVectorize
   (const GeneratorMatrixGen<T> &gm, unsigned k)
   : HeapAllocatedGeneratorMatrixGen<T> (
       powInt (gm.getBase(), k),  // base
       gm.getDimension(),
       gm.getM(),
       (gm.getPrecision() - 1) / k + 1, // prec
       true)        // performe memory allocation
{
   for (unsigned d = 0; d < dim; ++d)
   {
      for (unsigned r = 0; r < m; ++r)
      {
         for (unsigned b = 0; b < prec; ++b)
         {
            T x = 0;
            for (unsigned i = b*k; i < gm.getPrecision() && i < (b+1) * k; ++i)
            {
               x = x * k + gm (d, r, i);
            }
            set (d, r, b, x);
         }
      }
   }
}


namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template class GeneratorMatrixGen<X>; \
   template class MutableGeneratorMatrixGen<X>; \
   template class HeapAllocatedGeneratorMatrixGen<X>; \
   template class GeneratorMatrixGenCopy<X>; \
   template class GeneratorMatrixGenVectorize<X>; \
   template void assign ( \
      const GeneratorMatrix &, unsigned, \
      MutableGeneratorMatrixGen<X> &, unsigned);

   HINTLIB_INSTANTIATE (unsigned char)
   HINTLIB_INSTANTIATE (unsigned short)
   HINTLIB_INSTANTIATE (u32)
#undef HINTLIB_INSTANTIATE
}

