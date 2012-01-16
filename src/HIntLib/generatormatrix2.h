/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schürer <rudolf.schuerer@sbg.ac.at>
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

   unsigned getBase() const  { return 2; }

   virtual u64 getVector (unsigned d, unsigned r) const = 0;

protected:
   GeneratorMatrix2Base
      (unsigned availableBits, unsigned _dim, unsigned _m, unsigned _prec)
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

   static const unsigned MAX_TOTALPREC = std::numeric_limits<T>::digits;
   static const unsigned CORR_DEFAULT_TOTALPREC_BASE2 =
        MAX_TOTALPREC < DEFAULT_TOTALPREC_BASE2 
      ? MAX_TOTALPREC : DEFAULT_TOTALPREC_BASE2;
   
   // Constructors

   explicit GeneratorMatrix2 (unsigned _dim);
   GeneratorMatrix2 (unsigned _dim, unsigned _m);
   GeneratorMatrix2 (unsigned _dim, unsigned _m, unsigned _prec);

   // Copy constructors

   GeneratorMatrix2 (const GeneratorMatrix2<T> &);
   GeneratorMatrix2 (const GeneratorMatrix &);

   ~GeneratorMatrix2 () { delete[] c; }

   // get (parts of) Matrix

   const T* getMatrix() const  { return c; }
   const T* operator() (unsigned r) const  { return &c[r*dim]; }
         T* operator() (unsigned r)        { return &c[r*dim]; }

   // get/set column vectors

   T  operator() (unsigned d, unsigned r) const  { return c[r*dim + d]; }
   T& operator() (unsigned d, unsigned r)        { return c[r*dim + d]; }
   void setv (unsigned d, unsigned r,           T x)  { c[r*dim + d] = x; }
   void setv (unsigned d, unsigned r, unsigned, T x)  { setv (d,r,x); }
   void makeZeroColumnVector (unsigned d, unsigned r) { setv (d,r,0); }

   // get/set row vectors

   u64  getPackedRowVector (unsigned d, unsigned b) const;
   void setPackedRowVector (unsigned d, unsigned b, u64 x);
   void makeZeroRowVector (unsigned d, unsigned b);

   // get/set digits

private:
   T mask (unsigned b) const  { return T(1) << (prec - b - 1); }

   unsigned char getdMask (unsigned d, unsigned r, T ma) const
      { return ((*this)(d,r) & ma) != 0; }

   void setd0Mask (unsigned d, unsigned r, T ma)
      { c[r*dim + d] &= ~ma; }
   void setd1Mask (unsigned d, unsigned r, T ma)
      { c[r*dim + d] |=  ma; }
   void setdiMask (unsigned d, unsigned r, T ma)
      { c[r*dim + d] ^=  ma; }
   void setdMask (unsigned d, unsigned r, T ma, unsigned char x)
      { if (x) setd1Mask (d, r, ma); else setd0Mask (d, r, ma); }

public:
   unsigned char operator() (unsigned d, unsigned r, unsigned b) const
      { return getdMask (d, r, mask (b)); }
   void setd0 (unsigned d, unsigned r, unsigned b) { setd0Mask (d,r,mask(b)); }
   void setd1 (unsigned d, unsigned r, unsigned b) { setd1Mask (d,r,mask(b)); }
   void setdi (unsigned d, unsigned r, unsigned b) { setdiMask (d,r,mask(b)); }
   void setd (unsigned d, unsigned r, unsigned b, unsigned char x)
      { setdMask (d, r, mask(b), x); }

   // virtual get/set

   virtual unsigned getDigit (unsigned d, unsigned r, unsigned b) const;
   virtual void setDigit  (unsigned d, unsigned r, unsigned b, unsigned x);

   virtual u64  vGetPackedRowVector (unsigned d, unsigned b) const;
   virtual void vSetPackedRowVector (unsigned d, unsigned b, u64 x);

   virtual u64 getVector (unsigned d, unsigned r) const;

   // Manipulations

   void adjustPrec (unsigned);
   void makeHammersley (unsigned d);
   void makeIdentityMatrix (unsigned d);
   void makeZeroMatrix (unsigned d);
   void makeZeroMatrix ();
   void prepareForGrayCode ();
   void restoreFromGrayCode ();
   void makeShiftNet ();
   void makeShiftNet (unsigned b);

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
      const GeneratorMatrix &, unsigned, GeneratorMatrix2<T> &, unsigned);
template<typename T>
void assign (const GeneratorMatrix &, GeneratorMatrix2<T> &);

}  // namespace HIntLib

#endif

