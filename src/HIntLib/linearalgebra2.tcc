/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration
 *
 *  Copyright (C) 2002  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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

#ifndef HINTLIB_LINEAR_ALGEBRA_2_TCC
#define HINTLIB_LINEAR_ALGEBRA_2_TCC 1

#include <HIntLib/linearalgebra2.h>

#include <HIntLib/bitop.h>
#include <HIntLib/array.h>

/*
 *  linearalgebra2
 *
 *  All functions defined in this file work on matrices over GF(2).
 *  The rows of the matrix are stored at *first,...,*(last-1).
 *  The coefficients of each row are stored in the lower-order bits of *i.
 */


/**
 *  packMatrix()
 */

template<typename Bi>
void
HIntLib::packMatrix (
      const unsigned char* m, int numCols, Bi first, Bi last)
{
   typedef typename std::iterator_traits<Bi>::value_type T;

   while (first != last)
   {
      T mask = 1;
      *first = 0;

      for (int col = 0; col < numCols; ++col)
      {
         if (*m++)  *first |= mask;
         mask <<= 1;
      }

      ++first;
   }
}


/**
 *  unpackMatrix()
 */

template<typename Bi>
void
HIntLib::unpackMatrix (Bi first, Bi last, int numCols, unsigned char* res)
{
   typedef typename std::iterator_traits<Bi>::value_type T;

   while (first != last)
   {
      T row = *first++;

      for (int col = 0; col < numCols; ++col)
      {
         *res++ = row & 1;
         row >>= 1;
      }
   }
}


/**
 *  isZeroMatrix()
 */

template<typename Bi>
bool
HIntLib::isZeroMatrix (Bi first, Bi last)
{
   while (first != last)  if (*first++ != 0)  return false;
   return true;
}


/**
 *  isIdentityMatrix()
 */

template<typename Bi>
bool
HIntLib::isIdentityMatrix (Bi first, Bi last)
{
   typedef typename std::iterator_traits<Bi>::value_type T;
   T x = 1;
   
   while (first != last)
   {
      if (*first++ != x)  return false;
      x <<= 1;
   }

   return true;
}


/**
 *  matrixVectorMul()
 *
 *  Returns the product of a matrix A and a column vector v.
 *
 *  This requires  O(nm)  steps.
 */

template<typename Bi>
typename std::iterator_traits<Bi>::value_type
HIntLib::matrixVectorMul (
      Bi first, Bi last, typename std::iterator_traits<Bi>::value_type v)
{
   typedef typename std::iterator_traits<Bi>::value_type T;

   T result = 0;

   --first;
   while (--last != first)
   {
      result <<= 1;
      T x = *last & v;   // coefficient wise multiplication

      // count number of 1s

      while (x)
      { 
         result ^= (x & 1);
         x >>= 1;
      }
   }

   return result;
}


/**
 *  vectorMatrixMul()
 *
 *  Returns the product of a row vector v and a matrix A.
 *
 *  This requires only  O(n)  steps.
 */

template<typename Bi>
typename std::iterator_traits<Bi>::value_type
HIntLib::vectorMatrixMul (
      typename std::iterator_traits<Bi>::value_type v, Bi first, Bi last)
{
   typedef typename std::iterator_traits<Bi>::value_type T;

   T result = 0;

   for (; first != last; ++first)
   {
      if (v & 1)  result ^= *first;
      v >>= 1;
   }

   return result;
}


/**
 *  matrixRank()
 *
 *  Determines the rank of a packed matrix over GF(2), i.e. the largest number
 *  of row vectors that are linearly independent.
 */

template<typename Bi>
int
HIntLib::matrixRank (Bi first, Bi last)
{
   typedef typename std::iterator_traits<Bi>::value_type T;

   int rank = 0;

   for (; first != last; ++first)
   {
      // Get the current row

      const T pivot = *first;

      if (pivot)
      {
         // find the non-zero column in the current row

         const T mask = pivot & ~(pivot - 1);

         // subtract pivot from all remaining rows to clear a column

         for (Bi i = first + 1; i < last; ++i)  if (*i & mask)  *i ^= pivot;

         ++ rank;
      }
   }

   return rank;
}


/**
 *  linearlyIndependent()
 *
 *  Checks if the given packed row vectors over GF2 are linearly independent.
 *
 *  The empty set is linearly independent.
 */

template<typename Bi>
bool
HIntLib::isLinearlyIndependent (Bi first, Bi last)
{
   typedef typename std::iterator_traits<Bi>::value_type T;

   switch (last - first)
   {
   case 0: return true;
   case 1: return first[0];
   case 2: return first[0] && first[1] && (first[0] ^ first[1]);
   case 3:
      {
         T x 
            = first[0]; if (! x)  return false;
         x ^= first[1]; if (! x)  return false;
         x ^= first[0]; if (! x)  return false;
         x ^= first[2]; if (! x)  return false;
         x ^= first[0]; if (! x)  return false;
         x ^= first[1]; if (! x)  return false;
         x ^= first[0]; if (! x)  return false;
         return true;
      }
   }

   while (last - first > 4)  // More than 4 rows left?
   {
      // If the next vector is zero, the vectors cannot be linearly independent

      if (! *first)  return false;

      // find a non-zero column in the current row

      const T* pivot = first;
      const T mask = *pivot & ~(*pivot - 1);

      // subtract pivot from all remaining rows to clear a column

      for (Bi i = ++first; i < last; ++i)  if (*i & mask)  *i ^= *pivot;
   }

   // exactly four rows are left

   T x
      = first[0]; if (! x)  return false;
   x ^= first[1]; if (! x)  return false;
   x ^= first[0]; if (! x)  return false;
   x ^= first[2]; if (! x)  return false;
   x ^= first[0]; if (! x)  return false;
   x ^= first[1]; if (! x)  return false;
   x ^= first[0]; if (! x)  return false;
   x ^= first[3]; if (! x)  return false;
   x ^= first[0]; if (! x)  return false;
   x ^= first[1]; if (! x)  return false;
   x ^= first[0]; if (! x)  return false;
   x ^= first[2]; if (! x)  return false;
   x ^= first[0]; if (! x)  return false;
   x ^= first[1]; if (! x)  return false;
   x ^= first[0]; if (! x)  return false;
   return true;
}


/**
 *  nullSpace()
 *
 *  Calculates a base of the null space of a matrix  M, i.e. the set of all
 *  vectors  x  with  M x = 0.
 *
 *  The base vectors of the null space are stored in the  numCols x numCols
 *  (sic!) matrix result as the first  dim   row vectors.
 *
 *  The dimension of the null space is returned.
 *
 *  See CANT, Alg. 2.3.1, which is based on TACP, 4.6.2, Alg. N
 */

template<typename Bi>
int
HIntLib::nullSpace (Bi first, Bi last, int numCols, Bi result)
{
   typedef typename std::iterator_traits<Bi>::value_type T;

   Array<bool> rowSelected (last - first, false);
   int col2row [std::numeric_limits<T>::digits];
   int* col2rowP = col2row;

   // Loop over all columns

   T lastColMask = T(1) << (numCols - 1);
   for (T colMask = 1; colMask <= lastColMask; colMask <<= 1)
   {
      // find pivot row in current column

      Bi rowI = first;
      int row = 0;

      while (rowI != last && (((*rowI & colMask) == 0) || rowSelected[row]))
      {
         ++row;
         ++rowI;
      }

      // Did we find a valid pivot?

      if (rowI != last)  // pivot found  ->  Gauss elimination
      {
         // Subtract pivot row from other rows to clear out 1s in current column

         for (Bi r = first; r != last; ++r)
         {
            if ((r != rowI) && (*r & colMask))  *r ^= *rowI;
         }

         rowSelected [row] = true;
         *col2rowP++ = row;
      }
      else // sorry, no pivot -> new base vector for nullspace
      {
         *col2rowP++ = -1;
      }
   }

   // output base of null space

   int dimNullSpace = 0;

   T kMask = 1;
   for (int k = 0; k < numCols; ++k)
   {
      if (col2row[k] == -1)
      {
         ++ dimNullSpace;

         *result = kMask;
         T mask = 1;

         for (int i = 0; i < numCols; ++i)
         {
            if (col2row[i] >= 0)
            {
               if (first [col2row[i]] & kMask)  *result |= mask;
            }

            mask <<= 1;
         }

         ++ result;
      }
      kMask <<= 1;
   }

   return dimNullSpace;
}

template<typename Bi>
int
HIntLib::nullSpaceT (Bi first, Bi last, Bi result)
{
   typedef typename std::iterator_traits<Bi>::value_type T;

   T rowsAllowed = ~T(0);

   int col2row [std::numeric_limits<T>::digits];
   int* col2rowP = col2row;

   // Loop over all columns

   for (Bi col = first; col != last; ++col)
   {
      // cout << "column: " << *col;
      // Try to find pivot element (row) in the current column

      if (T allowedCoefficients = *col & rowsAllowed)
      {
         int row = ls1 (allowedCoefficients);
         T rowMask = T(1) << row;

         // cout << "   row " << row << " (" << rowMask << ")" << endl;

         T x = *col ^ rowMask;  // these rows need to be updated

         // Subtract pivot row from other rows to clear out 1s in current column

         for (Bi c = col + 1; c != last; ++c)  if (*c & rowMask)  *c ^= x;

         *col = rowMask;
         rowsAllowed ^= rowMask;
         *col2rowP++ = row;
      }
      else // sorry, no pivot -> new base vector for nullspace
      {
         // cout << "    no" << endl;
         *col2rowP++ = -1;
      }
   }

   // output base of null space

   int dimNullSpace = 0;

   col2rowP = col2row;
   for (Bi col = first; col != last; ++col)
   {
      if (*col2rowP++ == -1)
      {
         ++ dimNullSpace;

         *result = T(1) << (col - first);
         T mask = 1;

         for (int i = 0; i < (last - first); ++i)
         {
            if (col2row[i] >= 0)
            {
               if (*col & (T(1) << col2row[i]))  *result |= mask;
            }

            mask <<= 1;
         }

         ++ result;
      }
   }

   return dimNullSpace;
}

#define HINTLIB_INSTANTIATE_LINEARALGEBRA2(X) \
   template void packMatrix (const unsigned char*,int,X,X);\
   template void unpackMatrix (X,X,int,unsigned char*);\
   template bool isZeroMatrix    (X,X);\
   template bool isIdentityMatrix(X,X);\
   template std::iterator_traits<X >::value_type \
      matrixVectorMul (X,X,std::iterator_traits<X >::value_type); \
   template std::iterator_traits<X >::value_type \
      vectorMatrixMul (std::iterator_traits<X >::value_type,X,X); \
   template bool isLinearlyIndependent(X,X);\
   template int matrixRank (X,X);\
   template int nullSpace  (X,X,int,X); \
   template int nullSpaceT (X,X,X);

#endif

