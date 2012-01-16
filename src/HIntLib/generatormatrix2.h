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


#ifndef HINTLIB_GENERATOR_MATRIX_2_H
#define HINTLIB_GENERATOR_MATRIX_2_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/generatormatrix.h>


namespace HIntLib
{

/**
 *  Generator Matrix 2 Base
 *
 *  The non-template base class for GeneratorMatrix2
 *
 *  A special implementation for base-2 is provided that uses bits to store
 *  the matrix elements and packs a full matrix-colum into a single word.
 */

class GeneratorMatrix2Base : public GeneratorMatrix
{
public:
   ~GeneratorMatrix2Base () {}

   void CArrayDump (std::ostream &) const;
#ifdef HINTLIB_BUILD_WCHAR
   void CArrayDump (std::wostream &) const;
#endif

   int getBase() const  { return 2; }

   virtual u64 getVector (int d, int r) const = 0;

protected:
   GeneratorMatrix2Base
      (int availableBits, int _dim, int _m, int _prec)
      : GeneratorMatrix (2, _dim, _m, _prec) {} 

   GeneratorMatrix2Base (const GeneratorMatrix2Base &gm)
      : GeneratorMatrix (gm)  {}
};


/**
 *  Generator Matrix 2
 *
 *  A base-2 Generator Matrix.  It proviedes all access functions to the table.
 */

template<typename T>
class GeneratorMatrix2 : public GeneratorMatrix2Base
{
public:
   typedef T BaseType;

   enum { MAX_TOTALPREC = std::numeric_limits<T>::digits };
   enum { CORR_DEFAULT_TOTALPREC_BASE2 =
        unsigned(MAX_TOTALPREC) < unsigned(DEFAULT_TOTALPREC_BASE2) 
      ? unsigned(MAX_TOTALPREC) : unsigned(DEFAULT_TOTALPREC_BASE2) };
   
   // Constructors

   explicit GeneratorMatrix2 (int _dim);
   GeneratorMatrix2 (int _dim, int _m);
   GeneratorMatrix2 (int _dim, int _m, int _prec);

   // Copy constructors

   GeneratorMatrix2 (const GeneratorMatrix2<T> &);
   GeneratorMatrix2 (const GeneratorMatrix &);

   ~GeneratorMatrix2 () { delete[] c; }

   // get (parts of) Matrix

   const T* getMatrix() const  { return c; }
   const T* operator() (int r) const  { return &c[r*dim]; }
         T* operator() (int r)        { return &c[r*dim]; }

   // get/set column vectors

   T  operator() (int d, int r) const  { return c[r*dim + d]; }
   T& operator() (int d, int r)        { return c[r*dim + d]; }
   void setv (int d, int r,      T x)  { c[r*dim + d] = x; }
   void setv (int d, int r, int, T x)  { setv (d,r,x); }
   void makeZeroColumnVector (int d, int r) { setv (d,r,0); }

   // get/set row vectors

   u64  getPackedRowVector (int d, int b) const;
   void setPackedRowVector (int d, int b, u64 x);
   void makeZeroRowVector (int d, int b);

   // get/set digits

private:
   T mask (int b) const  { return T(1) << (prec - b - 1); }

   unsigned char getdMask (int d, int r, T ma) const
      { return ((*this)(d,r) & ma) != 0; }

   void setd0Mask (int d, int r, T ma)
      { c[r*dim + d] &= ~ma; }
   void setd1Mask (int d, int r, T ma)
      { c[r*dim + d] |=  ma; }
   void setdiMask (int d, int r, T ma)
      { c[r*dim + d] ^=  ma; }
   void setdMask (int d, int r, T ma, int x)
      { if (x) setd1Mask (d, r, ma); else setd0Mask (d, r, ma); }

public:
   unsigned char operator() (int d, int r, int b) const
      { return getdMask (d, r, mask (b)); }
   void setd0 (int d, int r, int b) { setd0Mask (d,r,mask(b)); }
   void setd1 (int d, int r, int b) { setd1Mask (d,r,mask(b)); }
   void setdi (int d, int r, int b) { setdiMask (d,r,mask(b)); }
   void setd  (int d, int r, int b, int x)
      { setdMask (d, r, mask(b), x); }

   // virtual get/set

   virtual int getDigit (int d, int r, int b) const;
   virtual void setDigit  (int d, int r, int b, int x);

   virtual u64  vGetPackedRowVector (int d, int b) const;
   virtual void vSetPackedRowVector (int d, int b, u64 x);

   virtual u64 getVector (int d, int r) const;

   // Manipulations

   void adjustPrec (int);
   void makeHammersley (int d);
   void makeIdentityMatrix (int d);
   void makeZeroMatrix (int d);
   void makeZeroMatrix ();
   void prepareForGrayCode ();
   void restoreFromGrayCode ();
   void makeShiftNet ();
   void makeShiftNet (int b);

protected:
   void setMatrix (const T*);

private:
   void checkPrec() const;
   void allocate();

   T* c; 
};


// Comparison

template<typename T>
bool operator== (const GeneratorMatrix2<T> &, const GeneratorMatrix2<T> &);

/**
 *  assign()
 */

template<typename T>
void assign (
      const GeneratorMatrix &, int, GeneratorMatrix2<T> &, int);
template<typename T>
void assign (const GeneratorMatrix &, GeneratorMatrix2<T> &);

}  // namespace HIntLib

#endif

