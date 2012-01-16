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


#ifndef HINTLIB_GENERATOR_MATRIX_2_ROW_H
#define HINTLIB_GENERATOR_MATRIX_2_ROW_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/generatormatrix.h>


namespace HIntLib
{

/**
 *  Generator Matrix 2 Row Base
 *
 *  The non-template base class for GeneratorMatrix2Row
 *
 *  A special implementation for base-2 is provided that uses bits to store
 *  the matrix elements and packs a full matrix-row into a single word.
 */

class GeneratorMatrix2RowBase : public GeneratorMatrix
{
public:
   unsigned getBase() const  { return 2; }
   unsigned getVectorization() const      { return 1; }
   unsigned getNumOfLeadingDigits() const { return 1; }
   unsigned getNumOfMissingDigits() const { return 0; }

protected:
   GeneratorMatrix2RowBase (unsigned _dim, unsigned _m, unsigned _totalPrec)
      : GeneratorMatrix (2, 1, _dim, _m, _totalPrec) {} 

   GeneratorMatrix2RowBase (const GeneratorMatrix2RowBase &gm)
      : GeneratorMatrix (gm)  {}
};


/**
 *  Generator Matrix 2 Row
 *
 *  A base-2 Generator Matrix.  It proviedes all access functions to the table.
 */

template<class T>
class GeneratorMatrix2Row : public GeneratorMatrix2RowBase
{
public:
   static const unsigned MAX_M = std::numeric_limits<T>::digits;
   static const unsigned CORR_DEFAULT_M_BASE2
      = MAX_M < DEFAULT_M_BASE2
      ? MAX_M : DEFAULT_M_BASE2;
   
   typedef T BaseType;

   // Constructors

   explicit GeneratorMatrix2Row (unsigned _dim);
   GeneratorMatrix2Row (unsigned _dim, unsigned _m);
   GeneratorMatrix2Row (unsigned _dim, unsigned _m, unsigned _totalPrec);

   // Copy constructors

   GeneratorMatrix2Row (const GeneratorMatrix2Row<T> &);
   GeneratorMatrix2Row (const GeneratorMatrix &);

   // Truncate a given matrix

   GeneratorMatrix2Row (const GeneratorMatrix &, const GMCopy &);

   ~GeneratorMatrix2Row () { delete[] c; }

   // get (parts of) Matrix

   const T* getMatrix() const  { return c; }
   const T* operator() (unsigned d) const  { return &c[d*totalPrec]; }
         T* operator() (unsigned d)        { return &c[d*totalPrec]; }

   // get/set column vectors

   T  operator() (unsigned d, unsigned b) const  { return c[d*totalPrec + b]; }
   T& operator() (unsigned d, unsigned b)        { return c[d*totalPrec + b]; }
#if 0
   void setv (unsigned d, unsigned r,           T x)  { c[r*dim + d] = x; }
   void setv (unsigned d, unsigned r, unsigned, T x)  { setv (d,r,x); }
   void makeZeroColumnVector (unsigned d, unsigned r) { setv (d,r,0); }
#endif

   // get/set row vectors

   T getPackedRowVector (unsigned d, unsigned b) const
      { return operator()(d,b); }
   void setPackedRowVector (unsigned d, unsigned b, T x)
      { c[d*totalPrec + b] = x; }
   void makeZeroRowVector (unsigned d, unsigned b)
      { setPackedRowVector (d,b,0); }

   // get/set digits

private:
   T mask (unsigned r) const  { return T(1) << r; }

   unsigned char getdMask (unsigned d, T ma, unsigned b) const
      { return (operator()(d,b) & ma) != 0; }

   void setd0Mask (unsigned d, T ma, unsigned b)
      { c[d*totalPrec + b] &= ~ma; }
   void setd1Mask (unsigned d, T ma, unsigned b)
      { c[d*totalPrec + b] |=  ma; }
   void setdiMask (unsigned d, T ma, unsigned b)
      { c[d*totalPrec + b] ^=  ma; }
   void setdMask (unsigned d, T ma, unsigned b, unsigned char x)
      { if (x) setd1Mask (d, ma, b); else setd0Mask (d, ma, b); }

public:
   unsigned char operator() (unsigned d, unsigned r, unsigned b) const
      { return getdMask (d, mask (r), b); }
   void setd0 (unsigned d, unsigned r, unsigned b) { setd0Mask (d,mask(r),b); }
   void setd1 (unsigned d, unsigned r, unsigned b) { setd1Mask (d,mask(r),b); }
   void setdi (unsigned d, unsigned r, unsigned b) { setdiMask (d,mask(r),b); }
   void setd (unsigned d, unsigned r, unsigned b, unsigned char x)
   {
      setdMask (d, mask(r), b, x);
   }

   // virtual get

   virtual unsigned getDigit (unsigned d, unsigned r, unsigned b) const;
   virtual u64 getVector (unsigned d, unsigned r, unsigned) const;
   virtual void setDigit  (unsigned d, unsigned r, unsigned b, unsigned x);
   virtual void setVector (unsigned d, unsigned r, unsigned b, u64 x);

   virtual u64  vGetPackedRowVector (unsigned d, unsigned b) const;
   virtual void vSetPackedRowVector (unsigned d, unsigned b, u64 x);

   // Comparison
   
   bool operator== (const GeneratorMatrix2Row<T> &) const;

   // Manipulations

   void makeEquidistributedCoordinate (unsigned d);
   void makeIdentityMatrix (unsigned d);
   void makeZeroMatrix (unsigned d);
   void makeZeroMatrix ();
   // void prepareForGrayCode ();
   // void restoreFromGrayCode ();
   void makeShiftNet ();
   void makeShiftNet (unsigned b);

protected:
   void setMatrix (const T*);

private:
   void checkM() const;
   void allocate();

   T* c; 
};

/**
 *  assign()
 */

template<typename T>
void assign (
      const GeneratorMatrix &, unsigned, GeneratorMatrix2Row<T> &, unsigned);
template<typename T>
void assign (const GeneratorMatrix &, GeneratorMatrix2Row<T> &);

}  // namespace HIntLib

#endif

