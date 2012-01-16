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


#ifndef HINTLIB_GENERATOR_MATRIX_H
#define HINTLIB_GENERATOR_MATRIX_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <algorithm>
#include <iosfwd>

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/fallback_limits.h>
#endif

#include <HIntLib/bitop.h>
#include <HIntLib/mymath.h>


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

   unsigned getBase() const          { return base; }
   unsigned getVectorization() const { return vec; }
   unsigned getPrec() const          { return prec; }      // # vecBase digits
   unsigned getDimension() const     { return dim; }
   unsigned getM() const             { return m; }
   unsigned getTotalPrec() const     { return totalPrec; } // # base digits

   virtual void     setDigit  (unsigned d, unsigned r, unsigned b, unsigned x);
   virtual void     setVector (unsigned d, unsigned r, unsigned b, u64 x);
   virtual unsigned getDigit  (unsigned d, unsigned r, unsigned b) const = 0;
   virtual u64      getVector (unsigned d, unsigned r, unsigned b) const = 0;

   void dump       (std::ostream &) const;
   void vectorDump (std::ostream &) const;
   void libSeqDump (std::ostream &) const;
   void libSeqDump (const char *) const;
   void binaryDump (std::ostream &) const;
   void binaryDump (const char *) const;

   GeneratorMatrix& operator= (const GeneratorMatrix &);

protected:
   GeneratorMatrix (unsigned _base, unsigned _vec,
                    unsigned _dim, unsigned _m, unsigned _totalPrec)
      : base (_base),
        vec (_vec),
        prec ((_totalPrec-1) / _vec + 1),
        dim (_dim),
        m   (_m),
        totalPrec (_totalPrec) {}
   GeneratorMatrix (const GeneratorMatrix &);

   const unsigned base;
   const unsigned vec;
   const unsigned prec;
   const unsigned dim;
   const unsigned m;
         unsigned totalPrec;
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
   void CArrayDump (std::ostream &) const;

   unsigned getPrec() const  { return 1; }

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
   void set  (unsigned d, unsigned r, unsigned b, unsigned char x);

   virtual void setDigit  (unsigned d, unsigned r, unsigned b, unsigned x);
   virtual void setVector (unsigned d, unsigned r, unsigned b, u64 x);

   void adjustTotalPrec (unsigned);
   void makeEquidistributedCoordinate (unsigned);
   void makeIdentityMatrix (unsigned);
   void makeZeroMatrix (unsigned);
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
   HeapAllocatedGeneratorMatrix2 (
         unsigned _dim,
         unsigned _m         = std::numeric_limits<Index>::digits - 1,
         unsigned _totalPrec = std::numeric_limits<real> ::digits - 1);
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
 *  GM Copy
 */

class GMCopy
{
public:
   GMCopy ();

   GMCopy& dim (int x)          { dimValue = x; return *this; } 
   GMCopy& maxDim (int x)       { maxDimValue = x; return *this; }
   GMCopy& m (int x)            { mValue = x; return *this; }
   GMCopy& maxM (int x)         { maxMValue = x; return *this; }
   GMCopy& totalPrec (int x)    { totalPrecValue = x; return *this; }
   GMCopy& maxTotalPrec (int x) { maxTotalPrecValue = x; return *this; }
   GMCopy& vec (int x)   { vecValue = x; return *this; }
   GMCopy& noVec ()      { vecValue = 1; return *this; }
   GMCopy& keepVec()     { vecValue = -2; return *this; }
   // GMCopy& optimalVec()  { vecValue = -3; return *this; }
   GMCopy& equi()        { equiValue = true; return *this; }
   GMCopy& equi(bool e)  { equiValue = e; return *this; }

   unsigned getDimension (const GeneratorMatrix &) const;
   unsigned getM         (const GeneratorMatrix &) const;
   unsigned getTotalPrec (const GeneratorMatrix &, unsigned max = 0) const;
   unsigned getVectorization (const GeneratorMatrix &, unsigned bits) const;
   bool getEqui () const  { return equiValue; }
   void checkNoVec () const;

private:
   int dimValue;
   int maxDimValue;
   int mValue;
   int maxMValue;
   int totalPrecValue;
   int maxTotalPrecValue;
   int vecValue;
   bool equiValue;
};

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
   unsigned getVecBase() const  { return vecBase; }

   // no public constructor!

protected:
   // Constructor with direct initialization
   GeneratorMatrixGenBase
      (unsigned _base, unsigned _vec,
       unsigned _dim, unsigned _m, unsigned _totalPrec)
   : GeneratorMatrix (_base, _vec, _dim, _m, _totalPrec),
     dimPrec (dim*prec),
     vecBase (powInt (base, vec)) {}

   // protected copy constructor

   GeneratorMatrixGenBase (const GeneratorMatrix &gm)
      : GeneratorMatrix (gm), dimPrec (dim*prec), vecBase (powInt (base, vec))
      {}
   GeneratorMatrixGenBase (const GeneratorMatrixGenBase &gm)
      : GeneratorMatrix (gm), dimPrec (gm.dimPrec), vecBase (gm.vecBase) {}

   const unsigned dimPrec;
   const unsigned vecBase;
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
   
   bool operator== (const GeneratorMatrixGen<T> &) const;
   bool operator!= (const GeneratorMatrixGen<T> &gm) const
      { return ! this->operator==(gm); }

   ~GeneratorMatrixGen () {}   // Make GCC 3.2 happy

protected:

   GeneratorMatrixGen (
      unsigned _base, unsigned _vec,
      unsigned _dim, unsigned _m, unsigned _totalPrec)
      : GeneratorMatrixGenBase (_base, _vec, _dim, _m, _totalPrec), c(0)
   { checkVecBase(); }

   GeneratorMatrixGen (const GeneratorMatrixGen<T> &gm, const T* _c)
      : GeneratorMatrixGenBase (gm), c(const_cast<T*>(_c))  {}
   GeneratorMatrixGen (const GeneratorMatrix &gm, const T* _c)
      : GeneratorMatrixGenBase (gm), c(const_cast<T*>(_c))  { checkVecBase(); }

   T* c; 

private:
   void checkVecBase() const;
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

template<class T> int t_parameter2 (
   const GeneratorMatrix2<T> &, int lb, int ub, int maxRows);
template<class T> int t_parameter2 (
   const GeneratorMatrixGen<T> &, int lb, int ub, int maxRows);


/**
 *  Mutable Generator Matrix Gen
 */

template<class T>
class MutableGeneratorMatrixGen : public GeneratorMatrixGen<T>
{
public:
   void setv (unsigned d, unsigned r, unsigned b, T x)
      { c[r*dimPrec + d*prec + b] = x; }
   void setd (unsigned d, unsigned r, unsigned b,
              typename GeneratorMatrixGen<T>::D x);

   virtual void setDigit  (unsigned d, unsigned r, unsigned b, unsigned x);
   virtual void setVector (unsigned d, unsigned r, unsigned b, u64 x);

   void makeEquidistributedCoordinate (unsigned);
   void makeIdentityMatrix (unsigned);
   void makeZeroMatrix (unsigned);

protected:
   MutableGeneratorMatrixGen (
      unsigned _base, unsigned _vec,
      unsigned _dim, unsigned _m, unsigned _totalPrec)
      : GeneratorMatrixGen<T> (_base, _vec, _dim, _m, _totalPrec)  {}

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
      (unsigned _base, unsigned _vec,
       unsigned _dim, unsigned _m, unsigned _totalPrec);
   HeapAllocatedGeneratorMatrixGen
      (unsigned _base, unsigned _vec, unsigned _dim);
   ~HeapAllocatedGeneratorMatrixGen () { delete[] c; }

protected:
   HeapAllocatedGeneratorMatrixGen (
      unsigned _base, unsigned _vec,
      unsigned _dim, unsigned _m, unsigned _totalPrec, bool alloc)
      : MutableGeneratorMatrixGen<T> (_base, _vec, _dim, _m, _totalPrec)
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


void assign (const GeneratorMatrix &, unsigned,
                   GeneratorMatrix &, unsigned);
template<class T>
void assign (const GeneratorMatrix &, unsigned,
                   MutableGeneratorMatrix2<T> &, unsigned);
template<class T>
void assign (const GeneratorMatrix &, unsigned,
                   MutableGeneratorMatrixGen<T> &, unsigned);


typedef GeneratorMatrixGenCopy<unsigned char> GeneratorMatrixGenCopy8;
typedef GeneratorMatrixGenCopy<unsigned short> GeneratorMatrixGenCopy16;
typedef GeneratorMatrixGenCopy<u32> GeneratorMatrixGenCopy32;
typedef GeneratorMatrixGenCopy<u64> GeneratorMatrixGenCopy64;

}  // namespace HIntLib

#endif

