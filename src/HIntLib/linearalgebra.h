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


#ifndef HINTLIB_LINEAR_ALGEBRA_H
#define HINTLIB_LINEAR_ALGEBRA_H 1

#ifdef __GNUG__
#pragma interface
#endif

namespace HIntLib
{

/**
 *  Linear Algebra
 *
 *  Linear Algebra objects provide virtual methods implementing the optimal
 *  linear algebra code for various finite fields
 *  (p = 2, p = prime, p = 2^k, p = general prime power).
 */

class LinearAlgebra
{
public:
   virtual ~LinearAlgebra() {}

   virtual bool isLinearlyIndependent (
      unsigned char*, unsigned numRows, unsigned numCols) = 0;
   virtual unsigned matrixRank (
      unsigned char*, unsigned numRows, unsigned numCols) = 0;
   virtual unsigned numLinearlyIndependentVectors (
      unsigned char*, unsigned numRows, unsigned numCols) = 0;
   virtual unsigned nullSpace (
      unsigned char*, unsigned numRows, unsigned numCols, unsigned char*) = 0;
   virtual unsigned basisSupplement (
      unsigned char*, unsigned numRows, unsigned numCols, unsigned*) = 0;
   virtual unsigned basisSupplement (
      unsigned char*, unsigned numRows, unsigned numCols) = 0;
   virtual bool matrixInverse  (unsigned char*, unsigned num) = 0;
   virtual void matrixMul (
      const unsigned char*, const unsigned char*,
      unsigned numRows1, unsigned numRowsCols, unsigned numCols2,
      unsigned char*) = 0;
   virtual void matrixVectorMul (
      const unsigned char*, const unsigned char*,
      unsigned numRows, unsigned numCols, unsigned char* res) = 0;
   virtual void vectorMatrixMul (
      const unsigned char*, const unsigned char*,
      unsigned numRows, unsigned numCols, unsigned char* res) = 0;

   bool isIdentityMatrix (const unsigned char*, unsigned num);
   bool isZeroMatrix (const unsigned char*, unsigned numCols, unsigned numRows);

   void matrixMul (
      const unsigned char* m1, const unsigned char* m2,
      unsigned num, unsigned char* res)
   {
      matrixMul (m1, m2, num, num, num, res);
   }

   static LinearAlgebra* make (unsigned base);

protected:
   LinearAlgebra () {}
};

}  // namespace HIntLib

#endif

