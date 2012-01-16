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

#ifdef __GNUG__
#pragma implementation
#endif

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/linearalgebra.tcc>

#include <HIntLib/lookupfield.h>
#include <HIntLib/gf2.h>
#include <HIntLib/prime.h>
#include <HIntLib/bitop.h>

namespace L = HIntLib;


#if 0
/**
 *  rank()
 *
 *  Determines the rank of a packed matrix over GF(2).
 *
 *  Determines the rank of the matrix composed of the vectors
 *  *first,..,*(last-1).
 *  Each vector hold _columns_ elements from {0,1}, which have to be stored in
 *  the lower order bits.
 */

template<class Bi>
unsigned rank (Bi first, Bi last, unsigned columns)
{
   typedef typename std::iterator_traits<Bi>::value_type T;

   T mask (1);          // mask for current column
   unsigned rank = 0;   // number of l.i. vectors we have found so far

   while (columns != 0 && first != last)  // Any rows and columns left?
   {
      // cout <<  "   A, mask=" << mask << "  column=" << columns << endl;

      // find pivot row

      Bi i = first;

      while ((*i & mask) == 0)  // if the coefficient == 0, we have to go on
      {
         if (++i == last)    // if we run out of rows, try a differnt column
         {
            // try next column

            if (--columns == 0)  return rank;
            mask <<= 1;
            i = first;
            // cout << "   B, mask=" << mask << endl;
         }
      }

      ++ rank;  // we have found a new l.i. vector!

      // exchange it with first row

      T pivot = *i;
      *i = *first;
      // *first = pivot; // we are not using *first again, so no need to update

      // subtract pivot from all remaining rows

      while (++i != last)
      {
         if (*i & mask)  *i ^= pivot;
      }

      // continue with smaller matrix

      mask <<= 1;
      --columns;
      ++first;
   }

   return rank;
}
#endif




/**
 *  linearlyIndependent()
 *
 *  Checks if the given set of vectors over GF2 (packed into the bits of a
 *  word) are linearly independent.
 *
 *  The empty set is independent.
 */

template<class Bi>
bool L::isLinearlyIndependent (Bi first, Bi last)
{
   typedef typename std::iterator_traits<Bi>::value_type T;

   switch (last - first)
   {
   case 0: return true;
   case 1: return first[0];
   case 2: return first[0] && first[1] && (first[0] ^ first[1]);
   case 3: return first[0] && first[1] && first[2] &&
                  (first[0] ^ first[1]) &&
                  (first[1] ^ first[2]) &&
                  (first[0] ^ first[2]) &&
                  (first[0] ^ first[1] ^ first[2]);
   }

   while (last - first > 4)  // More than 4 rows left?
   {
      // Get the current row

      T pivot = *first;

      // find a non-zero column in the current row

      if (!pivot)  return false;
      T mask = T(1) << ls1 (pivot);

      // subtract pivot from all remaining rows to clear a column

      for (Bi i = ++first; i < last; ++i)  if (*i & mask)  *i ^= pivot;
   }

   // exactly four rows are left

   return first[0] && first[1] && first[2] && first[3] &&
         (first[0] ^ first[1]) &&
         (first[0] ^ first[2]) &&
         (first[0] ^ first[3]) &&
         (first[1] ^ first[2]) &&
         (first[1] ^ first[3]) &&
         (first[2] ^ first[3]) &&
         (first[1] ^ first[2] ^ first[3]) &&
         (first[0] ^ first[2] ^ first[3]) &&
         (first[0] ^ first[1] ^ first[3]) &&
         (first[0] ^ first[1] ^ first[2]) &&
         (first[0] ^ first[1] ^ first[2] ^ first[3]);
}


/***********  Linear Algebra  ************************************************/


namespace HIntLib
{

template<class A>
class LinearAlgebraImp : public LinearAlgebra
{
public:
   LinearAlgebraImp () {}
   LinearAlgebraImp (unsigned base) : field (base) {}
   bool isLinearlyIndependent (
      unsigned char *m, unsigned numCols, unsigned numRows);
   unsigned matrixRank (unsigned char*, unsigned numCols, unsigned numRows);
   bool matrixInverse  (unsigned char*, unsigned num);
   void matrixMul (
      unsigned char*, const unsigned char*, const unsigned char*, unsigned num);

private:
   A field;
};

}  // namespace HIntLib

template<class A>
bool L::LinearAlgebraImp<A>::isLinearlyIndependent (
   unsigned char* m, unsigned numCols, unsigned numRows)
{
   return L::isLinearlyIndependent (field, m, numCols, numRows);
}

template<class A>
unsigned L::LinearAlgebraImp<A>::matrixRank (
      unsigned char* m, unsigned numCols, unsigned numRows)
{
   return L::matrixRank (field, m, numCols, numRows);
}

template<class A>
bool L::LinearAlgebraImp<A>::matrixInverse (unsigned char* m, unsigned num)
{
   return L::matrixInverse (field, m, num);
}

template<class A>
void L::LinearAlgebraImp<A>::matrixMul (
      unsigned char* m, const unsigned char* m1, const unsigned char* m2,
      unsigned num)
{
   return L::matrixMul (field, m, m1, m2, num);
}


L::LinearAlgebra* L::makeLinearAlgebra (unsigned base)
{
   unsigned prime;
   unsigned power;
   Prime::factorPrimePower (base, prime, power);

   if (power == 1)
   {
      if (prime == 2) return new LinearAlgebraImp<GF2> ();
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
#define HINTLIB_INSTANTIATE(X) \
   template bool LinearAlgebraImp<X >::isLinearlyIndependent ( \
      unsigned char*, unsigned, unsigned); \
   template unsigned LinearAlgebraImp<X >::matrixRank ( \
      unsigned char*, unsigned, unsigned); \
   template bool LinearAlgebraImp<X >::matrixInverse  ( \
      unsigned char*, unsigned); \
   template void LinearAlgebraImp<X >::matrixMul ( \
      unsigned char*, const unsigned char*, const unsigned char*, unsigned);

   HINTLIB_INSTANTIATE (LookupField<unsigned char>)
   HINTLIB_INSTANTIATE (LookupFieldPow2<unsigned char>)
   HINTLIB_INSTANTIATE (LookupFieldPrime<unsigned char>)
#undef HINTLIB_INSTANTIATE

   HINTLIB_INSTANTIATE_LINEARALGEBRA (LookupField<unsigned char>)
   HINTLIB_INSTANTIATE_LINEARALGEBRA (LookupFieldPow2<unsigned char>)
   HINTLIB_INSTANTIATE_LINEARALGEBRA (LookupFieldPrime<unsigned char>)

   // base 2

#define HINTLIB_INSTANTIATE(X) \
   template bool isLinearlyIndependent(X, X);

   HINTLIB_INSTANTIATE (u32*)
#ifdef HINTLIB_U32_NOT_EQUAL_U64
   HINTLIB_INSTANTIATE (u64*)
#endif

#undef HINTLIB_INSTANTIATE
}

