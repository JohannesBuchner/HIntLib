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


#ifndef HINTLIB_GENERATOR_MATRIX_2_ROW_H
#define HINTLIB_GENERATOR_MATRIX_2_ROW_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
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
   int getBase() const  { return 2; }

protected:
   GeneratorMatrix2RowBase (int _dim, int _m, int _prec)
      : GeneratorMatrix (2, _dim, _m, _prec) {} 

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
   enum { MAX_M = std::numeric_limits<T>::digits };
   enum { CORR_DEFAULT_M_BASE2
      = unsigned(MAX_M) < unsigned(DEFAULT_M_BASE2)
      ? unsigned(MAX_M) : unsigned(DEFAULT_M_BASE2) };
   
   typedef T BaseType;

   // Constructors

   explicit GeneratorMatrix2Row (int _dim);
   GeneratorMatrix2Row (int _dim, int _m);
   GeneratorMatrix2Row (int _dim, int _m, int _prec);

   // Copy constructors

   GeneratorMatrix2Row (const GeneratorMatrix2Row<T> &);
   GeneratorMatrix2Row (const GeneratorMatrix &);

   ~GeneratorMatrix2Row () { delete[] c; }

   // get (parts of) Matrix

   const T* getMatrix() const  { return c; }
   const T* operator() (int d) const  { return &c[d*prec]; }
         T* operator() (int d)        { return &c[d*prec]; }

   // get/set column vectors

   T  operator() (int d, int b) const  { return c[d*prec + b]; }
   T& operator() (int d, int b)        { return c[d*prec + b]; }

   // get/set row vectors

   T getPackedRowVector (int d, int b) const
      { return (*this)(d,b); }
   void setPackedRowVector (int d, int b, T x)
      { c[d*prec + b] = x; }
   void makeZeroRowVector (int d, int b)
      { setPackedRowVector (d,b,0); }

   // get/set digits

private:
   T mask (int r) const  { return T(1) << r; }

   unsigned char getdMask (int d, T ma, int b) const
      { return ((*this)(d,b) & ma) != 0; }

   void setd0Mask (int d, T ma, int b)
      { c[d*prec + b] &= ~ma; }
   void setd1Mask (int d, T ma, int b)
      { c[d*prec + b] |=  ma; }
   void setdiMask (int d, T ma, int b)
      { c[d*prec + b] ^=  ma; }
   void setdMask (int d, T ma, int b, int x)
      { if (x) setd1Mask (d, ma, b); else setd0Mask (d, ma, b); }

public:
   unsigned char operator() (int d, int r, int b) const
      { return getdMask (d, mask (r), b); }
   void setd0 (int d, int r, int b) { setd0Mask (d,mask(r),b); }
   void setd1 (int d, int r, int b) { setd1Mask (d,mask(r),b); }
   void setdi (int d, int r, int b) { setdiMask (d,mask(r),b); }
   void setd  (int d, int r, int b, int x)
   {
      setdMask (d, mask(r), b, x);
   }

   // virtual get

   virtual int  getDigit (int d, int r, int b) const;
   virtual void setDigit (int d, int r, int b, int x);

   virtual u64  vGetPackedRowVector (int d, int b) const;
   virtual void vSetPackedRowVector (int d, int b, u64 x);

   // Comparison
   
   bool operator== (const GeneratorMatrix2Row<T> &) const;

   // Manipulations

   void makeHammersley (int d);
   void makeIdentityMatrix (int d);
   void makeZeroMatrix (int d);
   void makeZeroMatrix ();
   // void prepareForGrayCode ();
   // void restoreFromGrayCode ();
   void makeShiftNet ();
   void makeShiftNet (int b);

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
      const GeneratorMatrix &, int, GeneratorMatrix2Row<T> &, int);
template<typename T>
void assign (const GeneratorMatrix &, GeneratorMatrix2Row<T> &);

}  // namespace HIntLib

#endif

