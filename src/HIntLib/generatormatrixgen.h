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


#ifndef HINTLIB_GENERATOR_MATRIX_GEN_H
#define HINTLIB_GENERATOR_MATRIX_GEN_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/generatormatrix.h>

#include <HIntLib/hlmath.h>


namespace HIntLib
{

/**
 *  Generator Matrix Gen Base
 *
 *  The non-template base class for Generator Matrixes Gen
 */

class GeneratorMatrixGenBase : public GeneratorMatrix
{
protected:
   // Constructor with direct initialization

   GeneratorMatrixGenBase
      (int _base, int _dim, int _m, int _prec)
   : GeneratorMatrix (_base, _dim, _m, _prec), dimPrec (dim * prec) {}

   // protected copy constructor

   GeneratorMatrixGenBase (const GeneratorMatrix &gm);
   GeneratorMatrixGenBase (const GeneratorMatrixGenBase &gm)
      : GeneratorMatrix (gm), dimPrec (gm.dimPrec)  {}

   const int dimPrec;
};


/**
 *  Generator Matrix Gen
 *
 *  Generator Matrix with general base
 */

template<typename T>
class GeneratorMatrixGen : public GeneratorMatrixGenBase
{
public:
   // Constructing a new matrix

   GeneratorMatrixGen (int _base, int _dim);
   GeneratorMatrixGen (int _base, int _dim, int _m);
   GeneratorMatrixGen (int _base, int _dim, int _m, int _prec);

   // Copying a given matrix

   GeneratorMatrixGen (const GeneratorMatrixGen<T> &);
   GeneratorMatrixGen (const GeneratorMatrix &);

   ~GeneratorMatrixGen()  { delete[] c; }

   // get (parts of) the matrix

   const T* getMatrix() const  { return c; }
   const T* operator() (int r) const  { return c + r*dimPrec; }

   // get/set column vectors

   const T* operator() (int d, int r) const { return c + r*dimPrec + d*prec; }
         T* operator() (int d, int r)       { return c + r*dimPrec + d*prec; }
   void makeZeroColumnVector (int d, int r);

   // get/set row vectors

   u64  getPackedRowVector (int d, int b) const;
   void setPackedRowVector (int d, int b, u64 x);
   void makeZeroRowVector  (int d, int b);

   // get/set digits

   T operator() (int d, int r, int b) const
      { return c[r*dimPrec + d*prec + b]; }
   void setd (int d, int r, int b, T x)  { c[r*dimPrec + d*prec + b] = x; }

   // Virtual set/get

   virtual int  getDigit (int d, int r, int b) const;
   virtual void setDigit (int d, int r, int b, int x);

   virtual u64  vGetPackedRowVector (int d, int b) const;
   virtual void vSetPackedRowVector (int d, int b, u64 x);

   // Manipulation

   void makeHammersley (int d);
   void makeIdentityMatrix (int d);
   void makeZeroMatrix ();
   void makeZeroMatrix (int d);
   void makeShiftNet (int b);
   void makeShiftNet ();

private:
   void allocate();
   void checkBase() const;

   T* c; 
};

template<typename T>
bool operator== (const GeneratorMatrixGen<T> &, const GeneratorMatrixGen<T> &);

GeneratorMatrixGen<unsigned char>* loadLibSeq (std::istream &);
GeneratorMatrixGen<unsigned char>* loadEdel   (std::istream &, int);
GeneratorMatrixGen<unsigned char>* loadBinary (std::istream &);
GeneratorMatrixGen<unsigned char>* loadNiederreiterXing (int dim);


template<typename T>
void
assign (const GeneratorMatrix &, int, GeneratorMatrixGen<T> &, int);

template<typename T>
void assign (const GeneratorMatrix &, GeneratorMatrixGen<T> &);

}  // namespace HIntLib

#endif

