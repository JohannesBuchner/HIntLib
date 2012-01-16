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


#ifndef HINTLIB_GENERATOR_MATRIX_GEN_H
#define HINTLIB_GENERATOR_MATRIX_GEN_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/generatormatrix.h>

#include <HIntLib/mymath.h>


namespace HIntLib
{

/**
 *  Generator Matrix Gen Base
 *
 *  The non-template base class for Generator Matrixes Gen
 */

class GeneratorMatrixGenBase : public GeneratorMatrix
{
public:
   // no public constructor!

   unsigned getVectorization()      const  { return 1; }
   unsigned getNumOfLeadingDigits() const  { return 1; }
   unsigned getNumOfMissingDigits() const  { return 0; }

protected:
   // Constructor with direct initialization

   GeneratorMatrixGenBase
      (unsigned _base, unsigned _dim, unsigned _m, unsigned _totalPrec)
   : GeneratorMatrix (_base, 1, _dim, _m, _totalPrec), dimPrec (dim * prec) {}

   // protected copy constructor

   GeneratorMatrixGenBase (const GeneratorMatrix &gm);
   GeneratorMatrixGenBase (const GeneratorMatrixGenBase &gm)
      : GeneratorMatrix (gm), dimPrec (gm.dimPrec)  {}

   const unsigned dimPrec;
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

   GeneratorMatrixGen (unsigned _base, unsigned _dim);
   GeneratorMatrixGen (unsigned _base, unsigned _dim, unsigned _m);
   GeneratorMatrixGen
      (unsigned _base, unsigned _dim, unsigned _m, unsigned _totalPrec);

   // Copying a given matrix

   GeneratorMatrixGen (const GeneratorMatrixGen<T> &);
   GeneratorMatrixGen (const GeneratorMatrix &);

   // Truncate a given matrix

   GeneratorMatrixGen (const GeneratorMatrixGen<T> &, const GMCopy &);
   GeneratorMatrixGen (const GeneratorMatrix &, const GMCopy &);

   ~GeneratorMatrixGen()  { delete[] c; }

   // get (parts of ) the matrix

   const T* getMatrix() const  { return c; }
   const T* operator() (unsigned r) const  { return &c[r*dimPrec]; }

   // get/set column vectors

   const T* operator() (unsigned d, unsigned r) const
      { return &c[r*dimPrec + d*prec]; }
         T* operator() (unsigned d, unsigned r)
      { return &c[r*dimPrec + d*prec]; }
   void makeZeroColumnVector (unsigned d, unsigned r);

   // get/set row vectors

   u64  getPackedRowVector (unsigned d, unsigned b) const;
   void setPackedRowVector (unsigned d, unsigned b, u64 x);
   void makeZeroRowVector (unsigned d, unsigned b);

   // get/set digits

   T operator() (unsigned d, unsigned r, unsigned b) const
      { return c[r*dimPrec + d*prec + b]; }
   void setd (unsigned d, unsigned r, unsigned b, T x)
      { c[r*dimPrec + d*prec + b] = x; }

   // Virtual set/get

   virtual unsigned getDigit  (unsigned d, unsigned r, unsigned b) const;
   virtual u64      getVector (unsigned d, unsigned r, unsigned b) const;
   virtual void setDigit  (unsigned d, unsigned r, unsigned b, unsigned x);
   virtual void setVector (unsigned d, unsigned r, unsigned b, u64 x);

   virtual u64  vGetPackedRowVector (unsigned d, unsigned b) const;
   virtual void vSetPackedRowVector (unsigned d, unsigned b, u64 x);

   // Manipulation

   void makeEquidistributedCoordinate (unsigned d);
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
bool operator== (const GeneratorMatrixGen<T> &, const GeneratorMatrixGen<T> &);

GeneratorMatrixGen<unsigned char>* loadLibSeq (std::istream &);
GeneratorMatrixGen<unsigned char>* loadLibSeq (const char*);
GeneratorMatrixGen<unsigned char>* loadBinary (std::istream &);
GeneratorMatrixGen<unsigned char>* loadBinary (const char*);
GeneratorMatrixGen<unsigned char>* loadNiederreiterXing (unsigned dim);


template<typename T>
void
assign (const GeneratorMatrix &, unsigned, GeneratorMatrixGen<T> &, unsigned);

template<typename T>
void assign (const GeneratorMatrix &, GeneratorMatrixGen<T> &);

}  // namespace HIntLib

#endif

