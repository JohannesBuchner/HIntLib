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

#include <algorithm>

#include <HIntLib/linearalgebra.h>

#include <HIntLib/array.h>


/**
 *  linearlyIndependent()
 *
 *  Determins of the set of _numRows_ vectors with _numCols_ entries from _A_
 *  is linearly independent.
 *
 *  The vectors are accessed using _m_ [vector * _numcols + coord]
 */

template<class A>
bool HIntLib::isLinearlyIndependent (
      const A &a, typename A::type *m, unsigned numCols, unsigned numRows)
{
   typedef typename A::type T;

   unsigned col = 0;
   unsigned row = 0;
   
   while (row < numRows && col < numCols)  // Any rows/columns left?
   {
      // find pivot row
      
      unsigned i = row;
          
      while (a.is0 (m [i * numCols + col]))  // if the coefficient is 0, go on
      {
         if (++i == numRows)    // if we run out of rows, try a differnt column
         {
            i = row;

            if (++col == numCols)  return false;
         }
      }

      // subtract pivot from all remaining rows
      
      T recip = a.recip (m [i * numCols + col]);
      
      for (unsigned j = i + 1; j < numRows; ++j)
      {
         T lc = m [j * numCols + col];
         if (! a.is0 (lc))
         {
            T x = a.neg (a.mul (lc, recip));
            for (unsigned k = col + 1; k < numCols; ++k)
            {
               a.addTo (m [j * numCols + k], a.mul (x, m [i * numCols + k]));
            }
         }
      }

      // exchange it with first row

      if (i != row)
      {
         for (unsigned k = col + 1; k < numCols; ++k)
         {
            m [i * numCols + k] = m [row * numCols + k];
         }
      }

      // continue with smaller matrix
      
      ++col;
      ++row;
   }

   return row == numRows;
}


/**
 *  matrix Rank ()
 */

template<class A>
unsigned HIntLib::matrixRank (
   const A& a, typename A::type* m, unsigned numCols, unsigned numRows)
{
   typedef typename A::type T;

   unsigned col = 0;
   unsigned row = 0;
   
   while (row < numRows && col < numCols)  // Any rows/columns left?
   {
      // find pivot row
      
      unsigned i = row;
          
      while (a.is0 (m [i * numCols + col]))  // if the coefficient is 0, go on
      {
         if (++i == numRows)    // if we run out of rows, try a differnt column
         {
            i = row;

            if (++col == numCols)  return row;
         }
      }

      // subtract pivot from all remaining rows
      
      T recip = a.recip (m [i * numCols + col]);
      
      for (unsigned j = i + 1; j < numRows; ++j)
      {
         T lc = m [j * numCols + col];
         if (! a.is0 (lc))
         {
            T x = a.neg (a.mul (lc, recip));
            for (unsigned k = col + 1; k < numCols; ++k)
            {
               a.addTo (m [j * numCols + k], a.mul (x, m [i * numCols + k]));
            }
         }
      }

      // exchange it with first row

      if (i != row)
      {
         for (unsigned k = col + 1; k < numCols; ++k)
         {
            m [i * numCols + k] = m [row * numCols + k];
         }
      }

      // continue with smaller matrix
      
      ++col;
      ++row;
   }

   return row;
}


/**
 *  matrix Inverse ()
 *
 *  Calculates the inverse of a square matrix.
 *
 *  The inverse matrix is constructed in-place using Gauss-Jordan elimination.
 *
 *  The algorithm is based on NR, 2.1, with the difference that we do not use
 *  full pivoting, since we have exact arithmetic. Searching for the maximal
 *  element to be used as pivot is replaced by searching for the first
 *  non-zero element.
 */

template<class A>
bool HIntLib::matrixInverse (const A &a, typename A::type *m, unsigned num)
{
   // This array records where each column can be found in the output matrix

   Array<unsigned> permutation (num);

   typedef typename A::type T;

   for (unsigned pos = 0; pos < num; ++pos)
   {
      {
         // search for pivot row

         unsigned row = pos;

         while (a.is0 (m[row * num + pos]))  if (++row == num)  return false;

         permutation [pos] = row;

         // swap current row with pivot row

         if (row != pos)
         {
            for (unsigned col = 0; col < num; ++col)
            {
               std::swap (m [pos * num + col], m [row * num + col]);
            }
         }
      }

      // update current row
      
      T pivot = a.recip (m[pos * num + pos]);
      m[pos * num + pos] = a.one();

      for (unsigned col = 0; col < num; ++col)
      {
         a.mulBy (m[pos * num + col], pivot);
      }

      // update remaining rows

      for (unsigned row = 0; row < num; ++row)
      {
         if (row == pos)  continue;

         T lc = a.neg (m [row * num + pos]);
         m[row * num + pos] = T();

         if (! a.is0 (lc))
         {
            for (unsigned col = 0; col < num; ++col)
            {
               a.addTo (m [row * num + col], a.mul (lc, m [pos * num + col]));
            }
         }
      }
   }

   // unscramble

   for (int col = num - 1; col >= 0; --col)
   {
      if (int (permutation[col]) == col)  continue;

      for (unsigned row = 0; row < num; ++row)
      {
         std::swap (m[row * num + col], m[row * num + permutation[col]]);
      }
   }
   
   return true;
}


/**
 *  matrix Mul ()
 */

template<class A>
void HIntLib::matrixMul (
      const A &a,
      typename A::type *m,
      const typename A::type *m1, const typename A::type *m2,
      unsigned size)
{
   for (unsigned row = 0; row < size; ++row)
   {
      for (unsigned col = 0; col < size; ++col)
      {
         typename A::type x = 0;

         for (unsigned i = 0; i < size; ++i)
         {
            a.addTo (x, a.mul (m1[row * size + i], m2[i * size + col]));
         }

         m [row * size + col] = x;
      }
   }
}


#define HINTLIB_INSTANTIATE_LINEARALGEBRA(X) \
   template bool matrixInverse        (const X&, X::type*, unsigned); \
   template bool isLinearlyIndependent(const X&, X::type*, unsigned, unsigned);\
   template unsigned matrixRank       (const X&, X::type*, unsigned, unsigned);\
   template void matrixMul ( \
         const X&, X::type *, const X::type*, const X::type*, unsigned size);

