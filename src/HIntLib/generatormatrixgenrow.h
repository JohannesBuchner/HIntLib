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


#ifndef HINTLIB_GENERATOR_MATRIX_GEN_ROW_H
#define HINTLIB_GENERATOR_MATRIX_GEN_ROW_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <memory>

#include <HIntLib/generatormatrix.h>

#include <HIntLib/hlmath.h>
#include <HIntLib/linearalgebra.h>


namespace HIntLib
{

/**
 *  Generator Matrix Gen Row Base
 *
 *  The non-template base class for Generator Matrixes Gen Row
 */

class GeneratorMatrixGenRowBase : public GeneratorMatrix
{
public:
   LinearAlgebra& la() const  { return *linAlg; }

protected:
   // Constructor with direct initialization

   GeneratorMatrixGenRowBase
      (unsigned _base, unsigned _dim, unsigned _m, unsigned _prec)
      : GeneratorMatrix (_base, _dim, _m, _prec),
        mPrec (m * _prec),
        linAlg (LinearAlgebra::make (_base))
    {}

   // protected copy constructor

   GeneratorMatrixGenRowBase (const GeneratorMatrix &gm);
   GeneratorMatrixGenRowBase (const GeneratorMatrixGenRowBase &gm)
      : GeneratorMatrix (gm),
        mPrec (gm.mPrec),
        linAlg (LinearAlgebra::make(base))
   {}

   const unsigned mPrec;
   std::auto_ptr<LinearAlgebra> linAlg;
};


/**
 *  Generator Matrix Gen Row
 *
 *  Generator Matrix with general base
 */

template<typename T>
class GeneratorMatrixGenRow : public GeneratorMatrixGenRowBase
{
public:
   // Constructing a new matrix

   GeneratorMatrixGenRow (unsigned _base, unsigned _dim);
   GeneratorMatrixGenRow (unsigned _base, unsigned _dim, unsigned _m);
   GeneratorMatrixGenRow
      (unsigned _base, unsigned _dim, unsigned _m, unsigned _prec);

   // Copying a given matrix

   GeneratorMatrixGenRow (const GeneratorMatrixGenRow<T> &);
   GeneratorMatrixGenRow (const GeneratorMatrix &);

   ~GeneratorMatrixGenRow()  { delete[] c; }

   // get (parts of) the matrix

   const T* getMatrix() const  { return c; }
         T* getMatrix()        { return c; }
   const T* operator() (unsigned d) const  { return c + d * mPrec; }
         T* operator() (unsigned d)        { return c + d * mPrec; }
   const T* operator() (unsigned d, unsigned b) const
      { return c + d * mPrec + b * m; }
         T* operator() (unsigned d, unsigned b)
      { return c + d * mPrec + b * m; }

   // get/set column vectors

   void makeZeroColumnVector (unsigned d, unsigned r);

   // get/set row vectors

   u64  getPackedRowVector (unsigned d, unsigned b) const;
   void setPackedRowVector (unsigned d, unsigned b, u64 x);
   void makeZeroRowVector (unsigned d, unsigned b);

   // get/set digits

   T operator() (unsigned d, unsigned r, unsigned b) const
      { return c [d*mPrec + b*m + r]; }
   void setd (unsigned d, unsigned r, unsigned b, T x)
      { c [d*mPrec + b*m + r] = x; }

   // Virtual set/get

   virtual unsigned getDigit  (unsigned d, unsigned r, unsigned b) const;
   virtual void setDigit  (unsigned d, unsigned r, unsigned b, unsigned x);

   virtual u64  vGetPackedRowVector (unsigned d, unsigned b) const;
   virtual void vSetPackedRowVector (unsigned d, unsigned b, u64 x);

   // Manipulation

   void makeHammersley (unsigned d);
   void makeIdentityMatrix (unsigned d);
   void makeZeroMatrix ();
   void makeZeroMatrix (unsigned d);
   void makeShiftNet (unsigned b);
   void makeShiftNet ();

private:
   void allocate();
   void checkBase() const;

   T* c; 
};

template<typename T>
bool operator== (const GeneratorMatrixGenRow<T> &,
                 const GeneratorMatrixGenRow<T> &);


template<typename T>
void
assign (const GeneratorMatrix &, unsigned, GeneratorMatrixGenRow<T> &, unsigned);

template<typename T>
void assign (const GeneratorMatrix &, GeneratorMatrixGenRow<T> &);

}  // namespace HIntLib

#endif

