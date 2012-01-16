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
 *  Generator Matrix with general base, vectorized
 */

template<typename T>
class GeneratorMatrixVec : public GeneratorMatrixVecBase
{
public:
   typedef T D;

   // Constructors

   GeneratorMatrixVec (unsigned _base, unsigned _vec, unsigned _dim);
   GeneratorMatrixVec
      (unsigned _base, unsigned _vec, unsigned _dim, unsigned _m);
   GeneratorMatrixVec
      (unsigned _base, unsigned _vec,
       unsigned _dim, unsigned _m, unsigned _totalPrec)
    : GeneratorMatrixVecBase (_base, _vec, _dim, _m, _totalPrec)
      { checkVecBase(); allocate(); }

   // Copy constructors

   GeneratorMatrixVec (const GeneratorMatrixVec<T> &);
   GeneratorMatrixVec (const GeneratorMatrix &);

   // Truncate a given Matrix

   GeneratorMatrixVec (const GeneratorMatrixVec<T> &, const GMCopy &);
   GeneratorMatrixVec (const GeneratorMatrix &, const GMCopy &);

   ~GeneratorMatrixVec ()  { delete[] c; }

   // get and set elements

   const T* getMatrix() const  { return c; }

   const T* operator() (unsigned r) const  { return &c[r*dimPrec]; }
   const T* operator() (unsigned d, unsigned r) const
      { return &c[r*dimPrec + d*prec]; }
   T operator() (unsigned d, unsigned r, unsigned b) const
      { return c[r*dimPrec + d*prec + b]; }
   D getd (unsigned d, unsigned r, unsigned b) const;
   void setv (unsigned d, unsigned r, unsigned b, T x)
      { c[r*dimPrec + d*prec + b] = x; }
   void setd (unsigned d, unsigned r, unsigned b,
              typename GeneratorMatrixVec<T>::D x);
   
   u64  getPackedRowVector (unsigned d, unsigned b) const;
   void setPackedRowVector (unsigned d, unsigned b, u64 x);

   // Virtual set/get

   virtual unsigned getDigit  (unsigned d, unsigned r, unsigned b) const;
   virtual u64      getVector (unsigned d, unsigned r, unsigned b) const;
   virtual void setDigit  (unsigned d, unsigned r, unsigned b, unsigned x);
   virtual void setVector (unsigned d, unsigned r, unsigned b, u64 x);

   virtual u64  vGetPackedRowVector (unsigned d, unsigned b) const;
   virtual void vSetPackedRowVector (unsigned d, unsigned b, u64 x);

   // Manipulation

   void makeEquidistributedCoordinate (unsigned);
   void makeIdentityMatrix (unsigned);
   void makeZeroMatrix (unsigned);
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
void assign (
   const GeneratorMatrix &, unsigned, GeneratorMatrixVec<T> &, unsigned);

template<typename T>
void assign (const GeneratorMatrix &, GeneratorMatrixVec<T> &);

}  // namespace HIntLib

#endif

