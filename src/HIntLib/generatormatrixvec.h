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


#ifndef HINTLIB_GENERATOR_MATRIX_VEC_H
#define HINTLIB_GENERATOR_MATRIX_VEC_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/generatormatrix.h>

#include <HIntLib/hlmath.h>


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
   int getVec()     const  { return vec; }
   int getVecBase() const  { return vecBase; }
   int getVecPrec() const  { return vecPrec; }

   int getNumOfLeadingDigits() const { return prec - (vec * (vecPrec - 1)); }
   int getNumOfMissingDigits() const { return vecPrec * vec - prec; }

   // no public constructor!

protected:
   // Constructor with direct initialization

   GeneratorMatrixVecBase (int _base, int _dim, int _m, int _prec, int _vec);

   // protected copy constructor

   GeneratorMatrixVecBase (const GeneratorMatrix&, int _vec);
   GeneratorMatrixVecBase (const GeneratorMatrixVecBase&);

   const int vec;
   const int vecBase;
   const int vecPrec;

   const int dimPrec;
};


/**
 *  Generator Matrix Vec
 *
 *  Generator Matrix with general base, vectorized
 */

template<typename T>
class GeneratorMatrixVec : public GeneratorMatrixVecBase
{
public:
   typedef T D;

   // Constructors

   GeneratorMatrixVec (int _base, int _dim, int _vec);
   GeneratorMatrixVec (int _base, int _dim, int _m, int _vec);
   GeneratorMatrixVec (int _base, int _dim, int _m, int _prec, int _vec)
    : GeneratorMatrixVecBase (_base, _dim, _m, _prec, _vec)
      { checkVecBase(); allocate(); }

   // Copy constructors

   GeneratorMatrixVec (const GeneratorMatrixVec<T>&);
   GeneratorMatrixVec (const GeneratorMatrix&);
   GeneratorMatrixVec (const GeneratorMatrix&, int _vec);

   ~GeneratorMatrixVec ()  { delete[] c; }

   // get and set elements

   const T* getMatrix() const  { return c; }

   const T* operator() (int r) const  { return &c[r*dimPrec]; }
   const T* operator() (int d, int r) const
      { return &c[r*dimPrec + d*vecPrec]; }
   T operator() (int d, int r, int b) const
      { return c[r*dimPrec + d*vecPrec + b]; }
   D getd (int d, int r, int b) const;
   void setv (int d, int r, int b, T x) { c[r*dimPrec + d*vecPrec + b] = x; }
   void setd (int d, int r, int b, typename GeneratorMatrixVec<T>::D x);
   
   u64  getPackedRowVector (int d, int b) const;
   void setPackedRowVector (int d, int b, u64 x);

   // Virtual set/get

   virtual int  getDigit (int d, int r, int b) const;
   virtual void setDigit (int d, int r, int b, int x);

   virtual u64  vGetPackedRowVector (int d, int b) const;
   virtual void vSetPackedRowVector (int d, int b, u64 x);

   // Manipulation

   void makeHammersley (int);
   void makeIdentityMatrix (int);
   void makeZeroMatrix (int);
   void makeZeroMatrix ();

private:
   void allocate();
   void checkVecBase() const;

   T* c; 
};


// Comparision

template<typename T>
bool operator== (const GeneratorMatrixVec<T> &, const GeneratorMatrixVec<T> &);

/**
 *  assign ()
 */

template<typename T>
void assign (const GeneratorMatrix &, int, GeneratorMatrixVec<T> &, int);

template<typename T>
void assign (const GeneratorMatrix &, GeneratorMatrixVec<T> &);

}  // namespace HIntLib

#endif

