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


#ifndef HINTLIB_GENERATOR_MATRIX_VEC_H
#define HINTLIB_GENERATOR_MATRIX_VEC_H 1

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
 *  Generator Matrix Vec Base
 *
 *  The non-template base class for all Generator Matrixes with general base
 */

class GeneratorMatrixVecBase : public GeneratorMatrix
{
public:
   unsigned getVecBase() const  { return vecBase; }

   // no public constructor!

protected:
   // Constructor with direct initialization
   GeneratorMatrixVecBase
      (unsigned _base, unsigned _vec,
       unsigned _dim, unsigned _m, unsigned _totalPrec)
   : GeneratorMatrix (_base, _vec, _dim, _m, _totalPrec),
     dimPrec (dim*prec),
     vecBase (powInt (base, vec)) {}

   // protected copy constructor

   GeneratorMatrixVecBase (const GeneratorMatrix &gm)
      : GeneratorMatrix (gm), dimPrec (dim*prec), vecBase (powInt (base, vec))
      {}
   GeneratorMatrixVecBase (const GeneratorMatrixVecBase &gm)
      : GeneratorMatrix (gm), dimPrec (gm.dimPrec), vecBase (gm.vecBase) {}

   const unsigned dimPrec;
   const unsigned vecBase;
};


/**
 *  Generator Matrix Vec
 *
 *  Generator Matrix with general base
 *
 *  However, there is no constructor for actually creating the matrices.
 */

template<class T>
class GeneratorMatrixVec : public GeneratorMatrixVecBase
{
public:
   typedef T D;

   // no public constructor!

   const T* getMatrix() const  { return c; }

   const T* operator() (unsigned r) const  { return &c[r*dimPrec]; }
   const T* operator() (unsigned d, unsigned r) const
      { return &c[r*dimPrec + d*prec]; }
   T operator() (unsigned d, unsigned r, unsigned b) const
      { return c[r*dimPrec + d*prec + b]; }
   D getd (unsigned d, unsigned r, unsigned b) const;
   

   virtual unsigned getDigit  (unsigned d, unsigned r, unsigned b) const;
   virtual u64      getVector (unsigned d, unsigned r, unsigned b) const;

   // Comparision
   
   bool operator== (const GeneratorMatrixVec<T> &) const;
   bool operator!= (const GeneratorMatrixVec<T> &gm) const
      { return ! this->operator==(gm); }

   ~GeneratorMatrixVec () {}   // Make GCC 3.2 happy

protected:

   GeneratorMatrixVec (
      unsigned _base, unsigned _vec,
      unsigned _dim, unsigned _m, unsigned _totalPrec)
      : GeneratorMatrixVecBase (_base, _vec, _dim, _m, _totalPrec), c(0)
   { checkVecBase(); }

   GeneratorMatrixVec (const GeneratorMatrixVec<T> &gm, const T* _c)
      : GeneratorMatrixVecBase (gm), c(const_cast<T*>(_c))  {}
   GeneratorMatrixVec (const GeneratorMatrix &gm, const T* _c)
      : GeneratorMatrixVecBase (gm), c(const_cast<T*>(_c))  { checkVecBase(); }

   T* c; 

private:
   void checkVecBase() const;
   GeneratorMatrixVec (const GeneratorMatrixVec<T> &);  // don't copy
};


/**
 *  Mutable Generator Matrix Vec
 */

template<class T>
class MutableGeneratorMatrixVec : public GeneratorMatrixVec<T>
{
public:
   void setv (unsigned d, unsigned r, unsigned b, T x)
      { c[r*dimPrec + d*prec + b] = x; }
   void setd (unsigned d, unsigned r, unsigned b,
              typename GeneratorMatrixVec<T>::D x);

   virtual void setDigit  (unsigned d, unsigned r, unsigned b, unsigned x);
   virtual void setVector (unsigned d, unsigned r, unsigned b, u64 x);

   void makeEquidistributedCoordinate (unsigned);
   void makeIdentityMatrix (unsigned);
   void makeZeroMatrix (unsigned);
   void makeZeroMatrix ();

protected:
   MutableGeneratorMatrixVec (
      unsigned _base, unsigned _vec,
      unsigned _dim, unsigned _m, unsigned _totalPrec)
      : GeneratorMatrixVec<T> (_base, _vec, _dim, _m, _totalPrec)  {}

   MutableGeneratorMatrixVec
      (const MutableGeneratorMatrixVec<T> &gm, const T* _c)
      : GeneratorMatrixVec<T> (gm, _c)  {}
   MutableGeneratorMatrixVec (const GeneratorMatrix &gm, const T* _c)
      : GeneratorMatrixVec<T> (gm, _c)  {}
};


/**
 *  Heap Allocated Generator Matrix Vec
 *
 *  A Mutable Generator Matrix with *c allocated on the heap
 */

template<class T>
class HeapAllocatedGeneratorMatrixVec: public MutableGeneratorMatrixVec<T>
{
public:
   HeapAllocatedGeneratorMatrixVec
      (unsigned _base, unsigned _vec, unsigned _dim);
   HeapAllocatedGeneratorMatrixVec
      (unsigned _base, unsigned _vec, unsigned _dim, unsigned _m);
   HeapAllocatedGeneratorMatrixVec
      (unsigned _base, unsigned _vec,
       unsigned _dim, unsigned _m, unsigned _totalPrec)
    : MutableGeneratorMatrixVec<T> (_base, _vec, _dim, _m, _totalPrec)
      { allocate(); }

   ~HeapAllocatedGeneratorMatrixVec () { delete[] c; }

protected:
   HeapAllocatedGeneratorMatrixVec (
      unsigned _base, unsigned _vec,
      unsigned _dim, unsigned _m, unsigned _totalPrec, bool alloc)
      : MutableGeneratorMatrixVec<T> (_base, _vec, _dim, _m, _totalPrec)
   { if (alloc)  allocate(); }

   HeapAllocatedGeneratorMatrixVec
      (const HeapAllocatedGeneratorMatrixVec<T> &gm, bool alloc)
      : MutableGeneratorMatrixVec<T> (gm, 0)  { if (alloc)  allocate(); }
   HeapAllocatedGeneratorMatrixVec (const GeneratorMatrix &gm,bool alloc)
      : MutableGeneratorMatrixVec<T> (gm, 0)  { if (alloc)  allocate(); }

   void allocate();
};

#if 0
MutableGeneratorMatrixVec<unsigned char>* loadLibSeq (std::istream &);
MutableGeneratorMatrixVec<unsigned char>* loadLibSeq (const char*);
MutableGeneratorMatrixVec<unsigned char>* loadBinary (std::istream &);
MutableGeneratorMatrixVec<unsigned char>* loadBinary (const char*);
MutableGeneratorMatrixVec<unsigned char>* loadNiederreiterXing (unsigned dim);
#endif


/**
 *  Generator Matrix Vec Copy
 *
 *  Creates a Generator Matrix from a given Generator Matrix.
 *
 *  Memory for the new matrix is automatically allocated from free store
 */

template<class T>
class GeneratorMatrixVecCopy : public HeapAllocatedGeneratorMatrixVec<T>
{
public:

   // Copy constructor
   GeneratorMatrixVecCopy (const GeneratorMatrixVec<T> &);

   // Copy constructor for arbitrary Generator Matrix
   GeneratorMatrixVecCopy (const GeneratorMatrix &);

   // Truncate a given Matrix
   GeneratorMatrixVecCopy (const GeneratorMatrixVec<T> &, const GMCopy &);
   GeneratorMatrixVecCopy (const GeneratorMatrix &, const GMCopy &);
};


template<class T>
void assign (const GeneratorMatrix &, unsigned,
                   MutableGeneratorMatrixVec<T> &, unsigned);
template<class T>
void assign (const GeneratorMatrix &, MutableGeneratorMatrixVec<T> &);

typedef GeneratorMatrixVecCopy<unsigned char> GeneratorMatrixVecCopy8;
typedef GeneratorMatrixVecCopy<unsigned short> GeneratorMatrixVecCopy16;
typedef GeneratorMatrixVecCopy<u32> GeneratorMatrixVecCopy32;
typedef GeneratorMatrixVecCopy<u64> GeneratorMatrixVecCopy64;

}  // namespace HIntLib

#endif

