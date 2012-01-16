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

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/linearalgebra.h>
#include <HIntLib/linearalgebragen.tcc>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/lookupfield.h>
#include <HIntLib/gf2.h>
#include <HIntLib/prime.h>

namespace L = HIntLib;


/***********  Linear Algebra  ************************************************/


namespace HIntLib
{

/**
 *  LinearAlgebra
 *
 *  The non-virtual members of LinearAlgebra
 */

bool
L::LinearAlgebra::isZeroMatrix (
      const unsigned char* m, unsigned numCols, unsigned numRows)
{
   const unsigned char* end = m + numCols * numRows;
   while (m != end)  if (*m++)  return false;
   return true;
}

bool
L::LinearAlgebra::isIdentityMatrix (const unsigned char* m, unsigned size)
{
   const unsigned char* end = m + size * size;
   unsigned pos = 0;

   while (m != end)
   {
      if (pos == 0)
      {
         if (*m != 1)  return false;
         pos = size;
      }
      else
      {
         if (*m != 0)  return false;
         --pos;
      }
      ++m;
   }

   return true;
}

/**
 *  LinearAlgebraImp
 *
 *  Template implementations of all the virtual members of LinearAlgebra
 */

template<class A>
class LinearAlgebraImp : public LinearAlgebra
{
public:
   LinearAlgebraImp (unsigned base) : field (base) {}
   bool isLinearlyIndependent (
      unsigned char *m, unsigned numRows, unsigned numCols);
   unsigned matrixRank (unsigned char*, unsigned numRows, unsigned numCols);
   unsigned numLinearlyIndependentVectors (
      unsigned char*, unsigned numRows, unsigned numCols);
   unsigned nullSpace (
      unsigned char*, unsigned numRows, unsigned numCols, unsigned char*);
   unsigned basisSupplement (
      unsigned char*, unsigned numRows, unsigned numCols, unsigned*);
   unsigned basisSupplement (
      unsigned char*, unsigned numRows, unsigned numCols);
   bool matrixInverse  (unsigned char*, unsigned num);
   void matrixMul (
      const unsigned char*, const unsigned char*,
      unsigned numRows1, unsigned numRowsCols, unsigned numCols2,
      unsigned char*);
   void matrixVectorMul (
      const unsigned char*, const unsigned char*,
      unsigned numRows, unsigned numCols, unsigned char*);
   void vectorMatrixMul (
      const unsigned char*, const unsigned char*,
      unsigned numRows, unsigned numCols, unsigned char*);

private:
   A field;
};

}  // namespace HIntLib

template<class A>
bool L::LinearAlgebraImp<A>::isLinearlyIndependent (
   unsigned char* m, unsigned numRows, unsigned numCols)
{
   return L::isLinearlyIndependent (field, m, numRows, numCols);
}

template<class A>
unsigned L::LinearAlgebraImp<A>::matrixRank (
      unsigned char* m, unsigned numRows, unsigned numCols)
{
   return L::matrixRank (field, m, numRows, numCols);
}

template<class A>
unsigned L::LinearAlgebraImp<A>::numLinearlyIndependentVectors (
      unsigned char* m, unsigned numRows, unsigned numCols)
{
   return L::numLinearlyIndependentVectors (field, m, numRows, numCols);
}

template<class A>
unsigned L::LinearAlgebraImp<A>::nullSpace (
      unsigned char* m, unsigned numRows, unsigned numCols, unsigned char* res)
{
   return L::nullSpace (field, m, numRows, numCols, res);
}

template<class A>
unsigned L::LinearAlgebraImp<A>::basisSupplement (
      unsigned char* m, unsigned numRows, unsigned numCols, unsigned* res)
{
   return L::basisSupplement (field, m, numRows, numCols, res);
}

template<class A>
unsigned L::LinearAlgebraImp<A>::basisSupplement (
      unsigned char* m, unsigned numRows, unsigned numCols)
{
   return L::basisSupplement (field, m, numRows, numCols);
}

template<class A>
bool L::LinearAlgebraImp<A>::matrixInverse (unsigned char* m, unsigned num)
{
   return L::matrixInverse (field, m, num);
}

template<class A>
void L::LinearAlgebraImp<A>::matrixMul (
      const unsigned char* m1, const unsigned char* m2,
      unsigned numRows1, unsigned numRowsCols, unsigned numCols2,
      unsigned char* res)
{
   return L::matrixMul (field, m1, m2, numRows1, numRowsCols, numCols2, res);
}

template<class A>
void L::LinearAlgebraImp<A>::matrixVectorMul (
      const unsigned char* m, const unsigned char* v,
      unsigned numRows, unsigned numCols, unsigned char* res)
{
   return L::matrixVectorMul (field, m, v, numRows, numCols, res);
}

template<class A>
void L::LinearAlgebraImp<A>::vectorMatrixMul (
      const unsigned char* v, const unsigned char* m,
      unsigned numRows, unsigned numCols, unsigned char* res)
{
   return L::vectorMatrixMul (field, v, m, numRows, numCols, res);
}


L::LinearAlgebra* L::LinearAlgebra::make (unsigned base)
{
   unsigned prime;
   unsigned power;
   Prime::factorPrimePower (base, prime, power);

   if (power == 1)
   {
      if (prime == 2) return new LinearAlgebraImp<GF2> (2);
      return new LinearAlgebraImp<LookupGaloisFieldPrime<unsigned char> >(base);
   }
   else
   {
      if (prime == 2)
         return new LinearAlgebraImp<LookupGaloisFieldPow2<unsigned char> >
            (base);
      return new LinearAlgebraImp<LookupGaloisField<unsigned char> > (base);
   }
}

namespace HIntLib
{
   HINTLIB_INSTANTIATE_LINEARALGEBRAGEN (GF2)
   HINTLIB_INSTANTIATE_LINEARALGEBRAGEN (LookupField<unsigned char>)
   HINTLIB_INSTANTIATE_LINEARALGEBRAGEN (LookupFieldPow2<unsigned char>)
   HINTLIB_INSTANTIATE_LINEARALGEBRAGEN (LookupFieldPrime<unsigned char>)

   HINTLIB_INSTANTIATE_LINEARALGEBRAGEN_T (unsigned char)
}

