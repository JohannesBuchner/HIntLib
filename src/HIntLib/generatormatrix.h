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


#ifndef GENERATOR_MATRIX_H
#define GENERATOR_MATRIX_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <algorithm>
#include <iosfwd>

#include <HIntLib/defaults.h>

#ifdef HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/hintlib_limits.h>
#endif

#include <HIntLib/bitop.h>


namespace HIntLib
{

/**
 *  Generator Matrix
 *
 *  Base class for all Generator Matrices
 *
 *  A Generator Matrix contains the set of matrices used for the generation of
 *  a digital net.
 *
 *  A Generator Matrix contains _dim_ matrices (for dimension 0,..,dim-1).
 *  Each matrix hold _m_*_prec_ entries from a ring with _base_ elements.
 *  The maximum number of points that can be generated from the matrix is given
 *  by  _base_ ^ _m_.
 *  Each point has _prec_ significant base-_base_ digits.
 */

class GeneratorMatrix
{
public:

   // no public constructor!
   virtual ~GeneratorMatrix() {};

   unsigned getBase() const      { return base; }
   unsigned getDimension() const { return dim; }
   unsigned getM() const         { return m; }
   unsigned getPrecision() const { return prec; }

   virtual void     setDigit (unsigned d, unsigned r, unsigned b, unsigned x);
   virtual unsigned getDigit (unsigned d, unsigned r, unsigned b) const = 0;

   void dump       (std::ostream &) const;
   void dumpLibSeq (std::ostream &) const;
   void dumpLibSeq (const char *) const;
   void dumpBinary (std::ostream &) const;
   void dumpBinary (const char *) const;

   GeneratorMatrix& operator= (const GeneratorMatrix &);

protected:
   GeneratorMatrix (
      unsigned _base, unsigned _dim, unsigned _m, unsigned _prec)
   : base (_base),
     dim (_dim),
     m   (_m),
     prec(_prec) {}
   GeneratorMatrix (const GeneratorMatrix &);

   const unsigned base;
   const unsigned dim;
   const unsigned m;
         unsigned prec;
};

bool operator== (const GeneratorMatrix &, const GeneratorMatrix &);

/**
 *  Generator Matrix 2 Base
 *
 *  The non-template base class for GeneratorMatrix2
 *
 *  A special implementation for base-2 is provided that uses bits to store
 *  the matrix elements and packs a full matrix-colum into a single work.
 *
 *  MAX_M defines the maximal allowed value of _m_.
 */

class GeneratorMatrix2Base : public GeneratorMatrix
{
public:
   virtual u64 getVector (unsigned d, unsigned r) const = 0;

   void binaryDump   (std::ostream &) const;
   void dumpAsCArray (std::ostream &) const;

protected:
   GeneratorMatrix2Base (unsigned _dim, unsigned _m, unsigned _prec)
      : GeneratorMatrix (2, _dim, _m, _prec) {}
   GeneratorMatrix2Base (const GeneratorMatrix2Base &gm)
      : GeneratorMatrix (gm) {}
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

   const T* getMatrix() const    { return c; }

   // geting matrix entries
   
   const T* operator() (unsigned r) const  { return &c[r*dim]; }
   T operator() (unsigned d, unsigned r) const  { return c[r*dim + d]; }
   char operator() (unsigned d, unsigned r, unsigned b) const
      { return bit (operator()(d,r), prec-(b+1)) != 0; }

   unsigned getDigit (unsigned d, unsigned r, unsigned b) const;
   u64 getVector (unsigned d, unsigned r) const;

   // Comparison
   
   bool operator== (const GeneratorMatrix2<T> &) const;
   bool operator!= (const GeneratorMatrix2<T> &gm) const
      { return ! this->operator==(gm); }

protected:

   GeneratorMatrix2 (
      unsigned _dim, unsigned _m, unsigned _prec, const T* _c = 0)
      : GeneratorMatrix2Base (_dim,_m,_prec), c(const_cast<T*>(_c))
      { checkPrec(); }
   GeneratorMatrix2 (const GeneratorMatrix2<T> &gm, const T* p)
      : GeneratorMatrix2Base (gm), c(const_cast<T*>(p))  {}
   GeneratorMatrix2 (const GeneratorMatrix2Base &gm, const T* p)
      : GeneratorMatrix2Base (gm), c(const_cast<T*>(p))  { checkPrec(); }

   void checkPrec() const;

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

   void set      (unsigned d, unsigned r, T x)  { c[r*dim + d] = x; }
   void setDigit (unsigned d, unsigned r, unsigned b, unsigned x);
   void set      (unsigned d, unsigned r, unsigned b, unsigned x)
   {
      T mask = T(1) << (prec-(b+1));
      if (x == 0)  c[r*dim + d] &= ~mask;
      else         c[r*dim + d] |= mask;
   }

   void adjustPrecision (unsigned);
   void makeEquidistributedCoordinate (unsigned);
   void makeIdentityMatrix (unsigned);
   void makeZeroMatrix (unsigned);
   void prepareForGrayCode ();
   void restoreFromGrayCode ();

protected:
   MutableGeneratorMatrix2 (
      unsigned _dim, unsigned _m, unsigned _prec)
      : GeneratorMatrix2<T> (_dim, _m, _prec, 0)  {}
   MutableGeneratorMatrix2 (const GeneratorMatrix2<T> &gm, const T* p)
      : GeneratorMatrix2<T> (gm, p)  {}
   MutableGeneratorMatrix2 (const GeneratorMatrix2Base &gm, const T* p)
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
   HeapAllocatedGeneratorMatrix2 (
         unsigned _dim,
         unsigned _m    = std::numeric_limits<Index>::digits - 1,
         unsigned _prec = std::numeric_limits<real>::digits - 1);
   ~HeapAllocatedGeneratorMatrix2 () { delete[] c; }

protected:
   HeapAllocatedGeneratorMatrix2 (
      unsigned _dim, unsigned _m, unsigned _prec, bool alloc)
      : MutableGeneratorMatrix2<T> (_dim,_m,_prec)
   { if (alloc)  allocate(); }

   HeapAllocatedGeneratorMatrix2 (const GeneratorMatrix2<T> &gm, bool alloc)
      : MutableGeneratorMatrix2<T> (gm, 0)  { if (alloc)  allocate(); }
   HeapAllocatedGeneratorMatrix2 (const GeneratorMatrix2Base &gm,bool alloc)
      : MutableGeneratorMatrix2<T> (gm, 0)  { if (alloc)  allocate(); }

   void allocate();
};


/**
 *  Generator Matrix 2 Copy
 *
 *  Creates a Generator Matrix from a given Generator Matrix.
 *
 *  Memory for the new matrix is automatically allocated from free store
 */

template<class T>
class GeneratorMatrixGen;

template<class T>
class GeneratorMatrix2Copy : public HeapAllocatedGeneratorMatrix2<T>
{
public:

   // Normal copy
   GeneratorMatrix2Copy (const GeneratorMatrix2<T> &);
   GeneratorMatrix2Copy (const GeneratorMatrix2Copy<T> &);

   // Normal copy for different base type
   template<class TT>
   GeneratorMatrix2Copy (const GeneratorMatrix2<TT> &);

   // Copy from arbitrary GeneratorMatrix
   GeneratorMatrix2Copy (const GeneratorMatrix &);

   // Truncate a given matrix
   GeneratorMatrix2Copy (
      const GeneratorMatrix2<T> &, unsigned, unsigned, unsigned, bool = false);
   GeneratorMatrix2Copy (
      const GeneratorMatrix &, unsigned, unsigned, unsigned, bool = false);
};

}  // namespace HIntLib


// *** Implementation ***


template<class T>
template<class TT>
inline
HIntLib::GeneratorMatrix2Copy<T>::GeneratorMatrix2Copy (
   const GeneratorMatrix2<TT> &gm)
: HeapAllocatedGeneratorMatrix2<T> (
     gm.getDimension(),
     gm.getM(),
     std::min(gm.getPrecision(), unsigned(std::numeric_limits<T>::digits)),
     true)
{
   int shift = gm.getPrecision() > prec ? gm.getPrecision() - prec : 0;

   for (unsigned r = 0; r < m; ++r)
   {
      for (unsigned d = 0; d < dim; ++d)  c[r*dim + d] = T(gm(d,r) >> shift);
   }
}


/*****************************************************************************/

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

protected:
   // Constructor with direct initialization
   GeneratorMatrixGenBase
      (unsigned _base, unsigned _dim, unsigned _m, unsigned _prec)
   : GeneratorMatrix (_base, _dim, _m, _prec), dimPrec (dim*prec) {}

   // protected copy constructor
   GeneratorMatrixGenBase (const GeneratorMatrix &gm)
      : GeneratorMatrix (gm), dimPrec (dim*prec) {}
   GeneratorMatrixGenBase (const GeneratorMatrixGenBase &gm)
      : GeneratorMatrix (gm), dimPrec (dim*prec) {}

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

   unsigned getDigit (unsigned d, unsigned r, unsigned b) const;

   // Comparision
   
   bool operator== (const GeneratorMatrixGen<T> &) const;
   bool operator!= (const GeneratorMatrixGen<T> &gm) const
      { return ! this->operator==(gm); }

   ~GeneratorMatrixGen () {}   // Make GCC 3.2 happy

protected:

   GeneratorMatrixGen (
      unsigned _base, unsigned _dim, unsigned _m, unsigned _prec)
      : GeneratorMatrixGenBase (_base,_dim,_m,_prec), c(0)
   { checkBase(); }

   GeneratorMatrixGen (const GeneratorMatrixGen<T> &gm, const T* _c)
      : GeneratorMatrixGenBase (gm), c(const_cast<T*>(_c))  {}
   GeneratorMatrixGen (const GeneratorMatrix &gm, const T* _c)
      : GeneratorMatrixGenBase (gm), c(const_cast<T*>(_c))  { checkBase(); }

   T* c; 

private:
   void checkBase() const;
   GeneratorMatrixGen (const GeneratorMatrixGen<T> &);  // dont copy
};


template<class T> int t_parameter (
   const GeneratorMatrix2<T> &, int lb, int ub, bool dimOpt = false);
template<class T> int t_parameter (
   const GeneratorMatrixGen<T> &, int lb, int ub, bool dimOpt = false);

template<class T> int t_parameter (const GeneratorMatrix2<T> &gm)
{
   return t_parameter (gm, 0, gm.getM());
}
template<class T> int t_parameter (const GeneratorMatrixGen<T> &gm)
{
   return t_parameter (gm, 0, gm.getM());
}


/**
 *  Mutable Generator Matrix Gen
 */

template<class T>
class MutableGeneratorMatrixGen : public GeneratorMatrixGen<T>
{
public:
   void setDigit (unsigned d, unsigned r, unsigned b, unsigned x);
   void set      (unsigned d, unsigned r, unsigned b, unsigned x)
      { c[r*dimPrec + d*prec + b] = x; }

   void makeEquidistributedCoordinate (unsigned);
   void makeIdentityMatrix (unsigned);
   void makeZeroMatrix (unsigned);

protected:
   MutableGeneratorMatrixGen (
      unsigned _base, unsigned _dim, unsigned _m, unsigned _prec)
      : GeneratorMatrixGen<T> (_base,_dim,_m,_prec)  {}

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
   HeapAllocatedGeneratorMatrixGen
      (unsigned _base, unsigned _dim, unsigned _m, unsigned _prec);
   HeapAllocatedGeneratorMatrixGen (unsigned _base, unsigned _dim);
   ~HeapAllocatedGeneratorMatrixGen () { delete[] c; }

protected:
   HeapAllocatedGeneratorMatrixGen (
      unsigned _base, unsigned _dim, unsigned _m, unsigned _prec, bool alloc)
      : MutableGeneratorMatrixGen<T> (_base,_dim,_m,_prec)
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
   GeneratorMatrixGenCopy (
      const GeneratorMatrixGen<T> &,
      unsigned, unsigned, unsigned, bool = false);

   GeneratorMatrixGenCopy (
      const GeneratorMatrix &, unsigned, unsigned, unsigned, bool = false);
};


void assign (const GeneratorMatrix &, unsigned,
                   GeneratorMatrix &, unsigned);
template<class T>
void assign (const GeneratorMatrix &, unsigned,
                   MutableGeneratorMatrix2<T> &, unsigned);
template<class T>
void assign (const GeneratorMatrix &, unsigned,
                   MutableGeneratorMatrixGen<T> &, unsigned);


/**
 *  Generator Matrix Gen Vectorize
 *
 *  Combines a number of digits into one digit
 */

template<class T>
class GeneratorMatrixGenVectorize : public HeapAllocatedGeneratorMatrixGen<T>
{
public:
   GeneratorMatrixGenVectorize (const GeneratorMatrixGen<T> &, unsigned);
};

typedef GeneratorMatrixGenCopy<unsigned char> GeneratorMatrixGenCopy8;
typedef GeneratorMatrixGenCopy<unsigned short> GeneratorMatrixGenCopy16;
typedef GeneratorMatrixGenCopy<u32> GeneratorMatrixGenCopy32;
typedef GeneratorMatrixGenCopy<u64> GeneratorMatrixGenCopy64;

}  // namespace HIntLib

#endif

