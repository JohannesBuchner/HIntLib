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


#ifndef HINTLIB_GENERATOR_MATRIX_GEN_H
#define HINTLIB_GENERATOR_MATRIX_GEN_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/generatormatrix.h>

#ifdef HINTLIB_HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/fallback_limits.h>
#endif

#include <HIntLib/mymath.h>


namespace HIntLib
{

/**
 *  Generator Matrix Gen Base
 *
 *  The non-template base class for all Generator Matrixes with general base
 */

class GeneratorMatrixGenBase : public GeneratorMatrix
{
public:
   // no public constructor!

   unsigned getVectorization() const       { return 1; }
   unsigned getNumOfLeadingDigits() const  { return 1; }
   unsigned getNumOfMissingDigits() const  { return 0; }

protected:
   // Constructor with direct initialization
   GeneratorMatrixGenBase
      (unsigned _base, unsigned _dim, unsigned _m, unsigned _totalPrec)
   : GeneratorMatrix (_base, 1, _dim, _m, _totalPrec),
     dimPrec (dim*prec) {}

   // protected copy constructor

   GeneratorMatrixGenBase (const GeneratorMatrix &gm)
      : GeneratorMatrix (gm), dimPrec (dim*prec)  {}
   GeneratorMatrixGenBase (const GeneratorMatrixGenBase &gm)
      : GeneratorMatrix (gm), dimPrec (gm.dimPrec)  {}

   const unsigned dimPrec;
};


/**
 *  Generator Matrix Gen
 *
 *  Generator Matrix with general base
 *
 *  However, there is no constructor for actually creating the matrices.
 */

template<class T>
class GeneratorMatrixGen : public GeneratorMatrixGenBase
{
public:
   // no public constructor!

   const T* getMatrix() const  { return c; }

   const T* operator() (unsigned r) const  { return &c[r*dimPrec]; }
   const T* operator() (unsigned d, unsigned r) const
      { return &c[r*dimPrec + d*prec]; }
   T operator() (unsigned d, unsigned r, unsigned b) const
      { return c[r*dimPrec + d*prec + b]; }

   virtual unsigned getDigit  (unsigned d, unsigned r, unsigned b) const;
   virtual u64      getVector (unsigned d, unsigned r, unsigned b) const;

   // Comparision
   
   bool operator== (const GeneratorMatrixGen<T> &) const;
   bool operator!= (const GeneratorMatrixGen<T> &gm) const
      { return ! this->operator==(gm); }

   ~GeneratorMatrixGen () {}   // Make GCC 3.2 happy

protected:

   GeneratorMatrixGen (
      unsigned _base, unsigned _dim, unsigned _m, unsigned _totalPrec)
      : GeneratorMatrixGenBase (_base, _dim, _m, _totalPrec), c(0)
   { checkBase(); }

   GeneratorMatrixGen (const GeneratorMatrixGen<T> &gm, const T* _c)
      : GeneratorMatrixGenBase (gm), c(const_cast<T*>(_c))  {}
   GeneratorMatrixGen (const GeneratorMatrix &gm, const T* _c)
      : GeneratorMatrixGenBase (gm), c(const_cast<T*>(_c))  { checkBase(); }

   T* c; 

private:
   void checkBase() const;
   GeneratorMatrixGen (const GeneratorMatrixGen<T> &);  // don't copy
};


/**
 *  Mutable Generator Matrix Gen
 */

template<class T>
class MutableGeneratorMatrixGen : public GeneratorMatrixGen<T>
{
public:
   void setd (unsigned d, unsigned r, unsigned b, T x)
      { c[r*dimPrec + d*prec + b] = x; }

   virtual void setDigit  (unsigned d, unsigned r, unsigned b, unsigned x);
   virtual void setVector (unsigned d, unsigned r, unsigned b, u64 x);

   void makeEquidistributedCoordinate (unsigned);
   void makeIdentityMatrix (unsigned);
   void makeZeroMatrix ();
   void makeZeroMatrix (unsigned);

protected:
   MutableGeneratorMatrixGen (
      unsigned _base, unsigned _dim, unsigned _m, unsigned _totalPrec)
      : GeneratorMatrixGen<T> (_base, _dim, _m, _totalPrec)  {}

   MutableGeneratorMatrixGen
      (const MutableGeneratorMatrixGen<T> &gm, const T* _c)
      : GeneratorMatrixGen<T> (gm, _c)  {}
   MutableGeneratorMatrixGen (const GeneratorMatrix &gm, const T* _c)
      : GeneratorMatrixGen<T> (gm, _c)  {}
};


/**
 *  Heap Allocated Generator Matrix Gen
 *
 *  A Mutable Generator Matrix with *c allocated on the heap
 */

template<class T>
class HeapAllocatedGeneratorMatrixGen: public MutableGeneratorMatrixGen<T>
{
public:
   HeapAllocatedGeneratorMatrixGen (unsigned _base, unsigned _dim);
   HeapAllocatedGeneratorMatrixGen (unsigned _base, unsigned _dim, unsigned _m);
   HeapAllocatedGeneratorMatrixGen
      (unsigned _base, unsigned _dim, unsigned _m, unsigned _totalPrec)
   : MutableGeneratorMatrixGen<T> (_base, _dim, _m, _totalPrec)
      { allocate(); makeZeroMatrix(); }

   ~HeapAllocatedGeneratorMatrixGen () { delete[] c; }

protected:
   HeapAllocatedGeneratorMatrixGen (
      unsigned _base,
      unsigned _dim, unsigned _m, unsigned _totalPrec, bool alloc)
      : MutableGeneratorMatrixGen<T> (_base, _dim, _m, _totalPrec)
   { if (alloc)  allocate(); }

   HeapAllocatedGeneratorMatrixGen
      (const HeapAllocatedGeneratorMatrixGen<T> &gm, bool alloc)
      : MutableGeneratorMatrixGen<T> (gm, 0)  { if (alloc)  allocate(); }
   HeapAllocatedGeneratorMatrixGen (const GeneratorMatrix &gm,bool alloc)
      : MutableGeneratorMatrixGen<T> (gm, 0)  { if (alloc)  allocate(); }

   void allocate();
};

MutableGeneratorMatrixGen<unsigned char>* loadLibSeq (std::istream &);
MutableGeneratorMatrixGen<unsigned char>* loadLibSeq (const char*);
MutableGeneratorMatrixGen<unsigned char>* loadBinary (std::istream &);
MutableGeneratorMatrixGen<unsigned char>* loadBinary (const char*);
MutableGeneratorMatrixGen<unsigned char>* loadNiederreiterXing (unsigned dim);


/**
 *  Generator Matrix Gen Copy
 *
 *  Creates a Generator Matrix from a given Generator Matrix.
 *
 *  Memory for the new matrix is automatically allocated from free store
 */

template<class T>
class GeneratorMatrixGenCopy : public HeapAllocatedGeneratorMatrixGen<T>
{
public:

   // Copy constructor
   GeneratorMatrixGenCopy (const GeneratorMatrixGen<T> &);

   // Copy constructor for arbitrary Generator Matrix
   GeneratorMatrixGenCopy (const GeneratorMatrix &);

   // Truncate a given Matrix
   GeneratorMatrixGenCopy (const GeneratorMatrixGen<T> &, const GMCopy &);
   GeneratorMatrixGenCopy (const GeneratorMatrix &, const GMCopy &);
};


template<class T>
void assign (const GeneratorMatrix &, unsigned,
                   MutableGeneratorMatrixGen<T> &, unsigned);
template<class T>
void assign (const GeneratorMatrix &, MutableGeneratorMatrixGen<T> &);

typedef GeneratorMatrixGenCopy<unsigned char> GeneratorMatrixGenCopy8;
typedef GeneratorMatrixGenCopy<unsigned short> GeneratorMatrixGenCopy16;
typedef GeneratorMatrixGenCopy<u32> GeneratorMatrixGenCopy32;
typedef GeneratorMatrixGenCopy<u64> GeneratorMatrixGenCopy64;

}  // namespace HIntLib

#endif

