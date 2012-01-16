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


#ifndef HINTLIB_GENERATOR_MATRIX_2_H
#define HINTLIB_GENERATOR_MATRIX_2_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/generatormatrix.h>

#ifdef HINTLIB_HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/fallback_limits.h>
#endif


namespace HIntLib
{

/**
 *  Generator Matrix 2 Base
 *
 *  The non-template base class for GeneratorMatrix2
 *
 *  A special implementation for base-2 is provided that uses bits to store
 *  the matrix elements and packs a full matrix-colum into a single work.
 */

class GeneratorMatrix2Base : public GeneratorMatrix
{
public:
   void CArrayDump (std::ostream &) const;

   unsigned getPrec() const  { return 1; }
   unsigned getBase() const  { return 2; }
   unsigned getNumOfLeadingDigits() const { return totalPrec; }
   unsigned getNumOfMissingDigits() const { return vec - totalPrec; }

protected:
   GeneratorMatrix2Base
      (unsigned availableBits, unsigned _dim, unsigned _m, unsigned _totalPrec)
      : GeneratorMatrix (2, availableBits, _dim, _m, _totalPrec) {} 

   GeneratorMatrix2Base (const GeneratorMatrix2Base &gm)
      : GeneratorMatrix (gm)  {}
};


/**
 *  Generator Matrix 2
 *
 *  A base-2 Generator Matrix.  It proviedes all access functions to the table.
 *
 *  However, there is no constructor for actually creating the matrices.
 */

template<class T>
class GeneratorMatrix2 : public GeneratorMatrix2Base
{
public:
   typedef T BaseType;

   // no public constructor!

   // geting matrix entries
   
   const T* getMatrix() const  { return c; }
   const T* operator() (unsigned r) const  { return &c[r*dim]; }
   T operator() (unsigned d, unsigned r) const  { return c[r*dim + d]; }
   unsigned char operator() (unsigned d, unsigned r, unsigned b) const
      { return (operator()(d,r) >> (totalPrec - b - 1)) & T(1); }

   virtual unsigned getDigit (unsigned d, unsigned r, unsigned b) const;
   virtual u64 getVector (unsigned d, unsigned r, unsigned) const;

   // Comparison
   
   bool operator== (const GeneratorMatrix2<T> &) const;
   bool operator!= (const GeneratorMatrix2<T> &gm) const
      { return ! operator==(gm); }

protected:

   GeneratorMatrix2 (
          unsigned _dim, unsigned _m, unsigned _totalPrec,
          const T* _c = 0)
      : GeneratorMatrix2Base
           (std::numeric_limits<T>::digits, _dim, _m, _totalPrec),
        c (const_cast<T*>(_c))
      { checkTotalPrec(); }
   GeneratorMatrix2 (const GeneratorMatrix2<T> &gm, const T* p)
      : GeneratorMatrix2Base (gm),
        c(const_cast<T*>(p)) {}

   void checkTotalPrec() const;

   T* c; 

private:
   GeneratorMatrix2 (const GeneratorMatrix2<T> &);   // no copy
};


/**
 *  Mutable Matrix 2
 */

template<class T>
class MutableGeneratorMatrix2 : public GeneratorMatrix2<T>
{
public:
   ~MutableGeneratorMatrix2 () {}   // make GCC 3.2 happy

   void setv (unsigned d, unsigned r, T x)  { c[r*dim + d] = x; }
   void setv (unsigned d, unsigned r, unsigned, T x)  { setv (d,r,x); }
   void setd (unsigned d, unsigned r, unsigned b, unsigned char x);

   virtual void setDigit  (unsigned d, unsigned r, unsigned b, unsigned x);
   virtual void setVector (unsigned d, unsigned r, unsigned b, u64 x);

   void adjustTotalPrec (unsigned);
   void makeEquidistributedCoordinate (unsigned);
   void makeIdentityMatrix (unsigned);
   void makeZeroMatrix (unsigned);
   void makeZeroMatrix ();
   void prepareForGrayCode ();
   void restoreFromGrayCode ();

protected:
   MutableGeneratorMatrix2 (unsigned _dim, unsigned _m, unsigned _totalPrec)
      : GeneratorMatrix2<T> (_dim, _m, _totalPrec, 0)  {}
   MutableGeneratorMatrix2 (const GeneratorMatrix2<T> &gm, const T* p)
      : GeneratorMatrix2<T> (gm, p)  {}
};


/**
 *  Heap Allocated Generator Matrix 2
 *
 *  A Mutable Generator Matrix with *c allocated on the heap
 */

template<class T>
class HeapAllocatedGeneratorMatrix2: public MutableGeneratorMatrix2<T>
{
public:
   HeapAllocatedGeneratorMatrix2 (unsigned _dim);
   HeapAllocatedGeneratorMatrix2 (unsigned _dim, unsigned _m);
   HeapAllocatedGeneratorMatrix2 (
         unsigned _dim, unsigned _m, unsigned _totalPrec)
   : MutableGeneratorMatrix2<T> (_dim, _m, _totalPrec)
   { allocate(); makeZeroMatrix(); }

   ~HeapAllocatedGeneratorMatrix2 () { delete[] c; }

protected:
   HeapAllocatedGeneratorMatrix2 (
      unsigned _dim, unsigned _m, unsigned _totalPrec, bool alloc)
   : MutableGeneratorMatrix2<T> (_dim, _m, _totalPrec)
      { if (alloc)  allocate(); }

   HeapAllocatedGeneratorMatrix2 (const GeneratorMatrix2<T> &gm,  bool alloc)
      : MutableGeneratorMatrix2<T> (gm, 0)  { if (alloc)  allocate(); }

   void allocate();
};


/**
 *  assign()
 */

template<class T>
void assign (const GeneratorMatrix &, unsigned,
                   MutableGeneratorMatrix2<T> &, unsigned);
template<class T>
void assign (const GeneratorMatrix &, MutableGeneratorMatrix2<T> &);


/**
 *  Generator Matrix 2 Copy
 *
 *  Creates a Generator Matrix from a given Generator Matrix.
 *
 *  Memory for the new matrix is automatically allocated from free store
 */

template<class T>
class GeneratorMatrix2Copy : public HeapAllocatedGeneratorMatrix2<T>
{
public:

   // Normal copy
   GeneratorMatrix2Copy (const GeneratorMatrix2<T> &);
   GeneratorMatrix2Copy (const GeneratorMatrix2Copy<T> &);

   // Copy from arbitrary GeneratorMatrix
   GeneratorMatrix2Copy (const GeneratorMatrix &);

   // Truncate a given matrix
   GeneratorMatrix2Copy (const GeneratorMatrix &, const GMCopy &);
};

}  // namespace HIntLib


#endif

