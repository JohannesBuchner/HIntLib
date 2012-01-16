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

#ifndef HINTLIB_LINEAR_ALGEBRA_GEN_TCC
#define HINTLIB_LINEAR_ALGEBRA_GEN_TCC 1

#include <algorithm>

#include <HIntLib/linearalgebragen.h>

#include <HIntLib/array.h>

/**
 *  linearalgebragen.tcc
 *
 *  A number of templates for doing linear algebra on C-style matrices.
 *
 *  The first argument is allways a const reference to the algebra to be used
 *  (in most cases, this has to be a field).
 *
 *  Square matrices are defined by a (const) pointer to the first element and
 *  the number of rows/columns.
 *
 *  Rectangular matrices are defined by a (const) pointer to the first element
 *  followed by the number of rows, followed by the number of columns.
 *
 *  In both cases, elements are accessed using
 *
 *           _base_ [ row * numColums + col ]
 *
 *  i.e. the matrices are stored as an array of row vectors.
 */

/**
 *  isZeroMatrix ()
 *
 *  Returns true if the given matrix is the zero matrix.
 *
 *  Used operations: is0()
 */

template<class A>
bool
HIntLib::isZeroMatrix (
      const A& a, const typename A::type* m, int numRows, int numCols)
{
   const typename A::type* end = m + numRows * numCols;
   while (m != end)  if (! a.is0 (*m++))  return false;
   return true;
}


/**
 *  isIdentityMatrix ()
 *
 *  Returns true if the given square matrix is the identity matrix.
 *
 *  Used operations: is0(), is1()
 */

template<class A>
bool
HIntLib::isIdentityMatrix (const A& a, const typename A::type* m, int size)
{
   const typename A::type* end = m + size * size;
   int pos = 0;

   while (m != end)
   {
      if (pos == 0)
      {
         if (! a.is1 (*m))  return false;
         pos = size;
      }
      else
      {
         if (! a.is0 (*m))  return false;
         --pos;
      }
      ++m;
   }

   return true;
}


/**
 *  matrixTranspose()
 */

template<typename T>
void
HIntLib::matrixTranspose (const T* in, int numRows, int numCols, T* out)
{
   for (int r = 0; r < numRows; ++r)
   {
      T* p = out++;
      for (int c = 0; c < numCols; ++c)
      {
         *p = *in++;
         p += numRows;
      }
   }
}


/**
 *  matrixMul ()
 *
 *  Multiplies two matrices.
 *
 *  m1 is a numRows1 x numRowsCols matrix.
 *  m2 is a numRowsCols x numCols2 matrix.
 *  The result is a numRows1 x numCols2 matrix.
 *
 *  Since the matrices we are dealing with are usually not too large, we use
 *  the naive algorithm requiering  n^3  multiplications and  n^3 - n^2
 *  additions.
 *
 *  Used operations:  addTo(), mul()    (therefore it also works for rings)
 */

template<class A>
void
HIntLib::matrixMul (
   const A &a,
   const typename A::type *m1, const typename A::type *m2,
   int numRows1, int numRowsCols, int numCols2,
   typename A::type *dest)
{
   typedef typename A::type T;

   for (int row = 0; row < numRows1; ++row)
   {
      const T* const m1end = &m1 [row * numRowsCols + numRowsCols];

      for (int col = 0; col < numCols2; ++col)
      {
         const T* m1i   = &m1 [row * numRowsCols + 0];
         const T* m2i   = &m2 [0 * numCols2 + col];

         *dest = a.mul (*m1i, *m2i);

         while (++m1i != m1end)
         {
            m2i += numCols2;
            a.addTo (*dest, a.mul (*m1i, *m2i));
         }

         ++dest;
      }
   }
}


/**
 *  matrixVectorMul ()
 *
 *  Multiplies a matrices with a vector.
 *
 *  m    is a numRows x numCols matrix.
 *  v    is a vector with numCols entries.
 *  dest is a vector with numRows entries.
 *
 *  Used operations:  addTo(), mul()    (therefore it also works for rings)
 */

template<class A>
void
HIntLib::matrixVectorMul (
   const A &a,
   const typename A::type* m, const typename A::type* v,
   int numRows, int numCols,
   typename A::type *dest)
{
   typedef typename A::type T;

   const T* const matrixEnd = m + numRows * numCols;
   const T* const vectorEnd = v + numCols;

   while (m != matrixEnd)
   {
      const T* vi = v;

      *dest = a.mul (*m++, *vi++);
      while (vi != vectorEnd)  a.addTo (*dest, a.mul (*m++, *vi++));

      ++dest;
   }
}


/**
 *  vectorMatrixMul ()
 *
 *  Multiplies a matrices with a vector.
 *
 *  m    is a numRows x numCols matrix.
 *  v    is a vector with numCols entries.
 *  dest is a vector with numRows entries.
 *
 *  Used operations:  addTo(), mul()    (therefore it also works for rings)
 */

template<class A>
void
HIntLib::vectorMatrixMul (
   const A &a,
   const typename A::type* v, const typename A::type* m,
   int numRows, int numCols,
   typename A::type *dest)
{
   typedef typename A::type T;

   const T* const matrixEnd = m + numCols;
   const T* const vectorEnd = v + numRows;

   while (m != matrixEnd)
   {
      const T* vi = v;
      const T* mi = m++;

      *dest = a.mul (*mi, *vi++);
      mi += numCols;

      while (vi != vectorEnd)
      {
         a.addTo (*dest, a.mul (*mi, *vi++));
         mi += numCols;
      }

      ++dest;
   }
}


/**
 *  islinearlyIndependent()
 *
 *  Determines if the set of row vectors is linearly independent.
 *
 *  Used operations: is0(), addTo(), neg(), mul(), recip()
 */

template<class A>
bool
HIntLib::isLinearlyIndependent (
   const A &a, typename A::type *m, int numRows, int numCols)
{
   typedef typename A::type T;

   if (numRows > numCols)  return false;

   // Loop over all rows

   for (int col = 0, row = 0; row != numRows; ++col, ++row)
   {
      // I have tried two different implementations. The second one seems to be
      // slightly faster.

#if 0
      // find a pivot element in the current row

      int pivotCol = col;

      while (a.is0 (m [row * numCols + pivotCol]))
      {
         if (++pivotCol == numCols)  return false;
      }

      // subtract pivot from all remaining rows
      
      const T recip = a.recip (m [row * numCols + pivotCol]);
      
      for (int r = row + 1; r < numRows; ++r)
      {
         const T lc = m [r * numCols + pivotCol];

         if (! a.is0 (lc))
         {
            const T x = a.neg (a.mul (lc, recip));

            for (int c = pivotCol + 1; c < numCols; ++c)
            {
               a.addTo (m [r * numCols + c], a.mul (x, m [row * numCols + c]));
            }
         }
      }

      // exchange pivot column with column col

      if (pivotCol != col)
      {
         for (int r = row + 1; r < numRows; ++r)
         {
            m [r * numCols + pivotCol] = m [r * numCols + col];
         }
      }
#else
      // find pivot element in current column (or the next...)
      
      int pivotRow = row;

      while (a.is0 (m [pivotRow * numCols + col]))
      {
         if (++pivotRow == numRows)  // if we run out of rows, try a next column
         {
            pivotRow = row;

            if (numCols - ++col < numRows - row)  return false;
         }
      }

      // subtract pivot from all remaining rows
      
      const T recip = a.recip (m [pivotRow * numCols + col]);
      
      for (int r = pivotRow + 1; r < numRows; ++r)
      {
         const T lc = m [r * numCols + col];

         if (! a.is0 (lc))
         {
            const T x = a.neg (a.mul (lc, recip));

            for (int c = col + 1; c < numCols; ++c)
            {
               a.addTo (m [r * numCols + c],
                        a.mul (x, m [pivotRow * numCols + c]));
            }
         }
      }

      // exchange it with first row

      if (pivotRow != row)
      {
         for (int c = col + 1; c < numCols; ++c)
         {
            m [pivotRow * numCols + c] = m [row * numCols + c];
         }
      }
#endif
   }

   return true;
}


/**
 *  matrixRank ()
 *
 *  Determines the rank of a rectangular matrix.
 *
 *  Used operations: is0(), addTo(), neg(), mul(), recip()
 */

template<class A>
int
HIntLib::matrixRank (const A& a, typename A::type* m, int numRows, int numCols)
{
   typedef typename A::type T;

   int col = 0;
   int row = 0;
   
   while (row < numRows && col < numCols)  // Any rows/columns left?
   {
      // find pivot row
      
      int i = row;
          
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
      
      for (int j = i + 1; j < numRows; ++j)
      {
         T lc = m [j * numCols + col];
         if (! a.is0 (lc))
         {
            T x = a.neg (a.mul (lc, recip));
            for (int k = col + 1; k < numCols; ++k)
            {
               a.addTo (m [j * numCols + k], a.mul (x, m [i * numCols + k]));
            }
         }
      }

      // exchange it with first row

      if (i != row)
      {
         for (int k = col + 1; k < numCols; ++k)
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
 *  numLinearlyIndependentVectors ()
 *
 *  Determines the largest  n  such that the first  n   row vectors of M are
 *  linearly independent.
 *
 *  Used operations: is0(), addTo(), neg(), mul(), recip()
 */

template<class A>
int
HIntLib::numLinearlyIndependentVectors (
   const A& a, typename A::type* m, int numRows, int numCols)
{
   typedef typename A::type T;

   if (numRows > numCols)  numRows = numCols;

   // Loop over all rows

   for (int col = 0, row = 0; ; ++col, ++row)
   {
      // find a pivot element in the current row

      int pivotCol = col;

      while (a.is0 (m [row * numCols + pivotCol]))
      {
         if (++pivotCol == numCols)  return row;
      }

      // is this the last row?

      int r = row + 1;
      if (r == numRows)  return numRows;

      // subtract pivot from all remaining rows
      
      const T recip = a.recip (m [row * numCols + pivotCol]);
      
      for ( ; r < numRows; ++r)
      {
         const T lc = m [r * numCols + pivotCol];

         if (! a.is0 (lc))
         {
            const T x = a.neg (a.mul (lc, recip));

            for (int c = pivotCol + 1; c < numCols; ++c)
            {
               a.addTo (m [r * numCols + c], a.mul (x, m [row * numCols + c]));
            }
         }
      }

      // exchange pivot column with column col

      if (pivotCol != col)
      {
         for (r = row + 1; r < numRows; ++r)
         {
            m [r * numCols + pivotCol] = m [r * numCols + col];
         }
      }
   }
}


/**
 *  nullSpace ()
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
 *
 *  Used operations: one(), is0(), neg(), addTo(), mul(), mulBy(), recip()
 */

template<class A>
int
HIntLib::nullSpace (
   const A& a, typename A::type* m, int numRows, int numCols,
   typename A::type* result)
{
   typedef typename A::type T;

   const T minusOne = a.neg (a.one());

   Array<bool> rowSelected (numRows, false);
   Array<int> col2row (numCols);
   int* col2rowP = col2row.begin();

   for (int col = 0; col < numCols; ++col)
   {
      // find pivot row in current column
      
      int row = 0;

      while (   row < numRows
             && (a.is0 (m [row * numCols + col]) || rowSelected[row]))
      {
         ++ row;
      }

      // Did we find a valid pivot?

      if (row < numRows)  // pivot found -> Gauss elimination
      {
         // multiply pivot row to make pivot = -1
      
         T recip = a.neg (a.recip (m [row * numCols + col]));
         m [row * numCols + col] = minusOne;
      
         for (int c = col + 1; c < numCols; ++c) 
         {
            a.mulBy (m [row * numCols + c], recip);
         }

         // subtract pivot row from all remaining rows to clear out 1s in
         // current column

         for (int r = 0; r < numRows; ++r) 
         {
            if (r == row)  continue;
   
            const T t = m [r * numCols + col];

            if (! a.is0 (t))
            {
               m [r * numCols + col] = T();

               for (int c = col + 1; c < numCols; ++c) 
               {
                  a.addTo (m [r * numCols + c],
                           a.mul (t, m [row * numCols + c]));
               }
            }
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

   for (int k = 0; k < numCols; ++k)
   {
      if (col2row[k] == -1)
      {
         ++ dimNullSpace;

         for (int i = 0; i < numCols; ++i)
         {
            if (col2row[i] >= 0) *result = m [col2row[i] * numCols + k];
            else if (k == i)     *result = a.one();
            else                 *result = T();

            ++ result;
         }
      }
   }

   return dimNullSpace;
}


/**
 *  basisSupplement ()
 *
 *  Given a  k x n  matrix  M, this function
 *
 *  1) determines the largest number k' such that the first k' row vectors of M
 *     are linearly independent (same as numLinearlyIndependentVectors())
 *  2) Determines indices result[k'],...,result[n-1] in {0,...,n} such that 
 *     the corresponding canonical base vectors supplement M to a complete
 *     basis of K^n.
 *
 *  See CANT, Alg. 2.3.6.
 *
 *  Used operations: one(), is0(), subFrom(), mul(), recip()
 */

template<class A>
int
HIntLib::basisSupplement (
   const A& a, typename A::type* m, int numRows, int n, int* result)
{
   typedef typename A::type T;

   // set result to identity matrix

   for (int i = 0; i < n; ++i)  result [i] = i;

   // process all rows

   T* mrow = m;
   const T* ub = m + n * numRows;

   for (int row = 0; row < numRows; mrow += n, ++row)
   {
      // find non-zero element in current row

      int col = row;

      while (a.is0 (mrow [col]))  if (++col == n)  return row;

      const T d = a.recip (mrow [col]);

      // update result

      if (col != row)  result [col] = result [row];
      
      // update m

      for (T* mr = mrow + n; mr != ub; mr += n)
      {
         if (row != col)  std::swap (mr [col], mr [row]); 
         a.mulBy (mr [row], d);

         for (int c = 0; c < n; ++c)
         {
            if (c != row && c != col)
            {
               a.subFrom (mr [c], a.mul (mrow [c], mr [row]));
            }
         }
      }
   }

   return numRows;
}


/**
 *  basisSupplement ()
 *
 *  Given a  k x n  matrix  M, this function
 *
 *  1) determines the largest number k' such that the first k' row vectors of M
 *     are linearly independent (same as numLinearlyIndependentVectors())
 *  2) supplements M to a regular n x n - Matrix, with the first k' row vectors
 *     identical to the ones in the original matrix.
 *
 *  See CANT, Alg. 2.3.6.
 *
 *  Used operations: one(), is0(), subFrom(), mul(), recip()
 */

template<class A>
int
HIntLib::basisSupplement (const A& a, typename A::type* m, int numRows, int n)
{
   typedef typename A::type T;

   // Copy matrix

   Array<T> work (numRows * n);
   std::copy (m, m + n * numRows, work.begin());

   Array<int> result (n);

   // Determine the vectors that have to be replaced by canonical basis vectors

   int validRows =
      basisSupplement (a, work.begin(), numRows, n, result.begin());

   // fill remaining vectors with zero

   m += n * validRows;
   std::fill (m, m + n * (n - validRows), T());

   // set ones at appropriate places

   const int* ub = &result[n];
   for (int* i = &result[validRows]; i != ub; ++i)
   {
      m [*i] = a.one();
      m += n;
   }

   return validRows;
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
 *
 *  Used operations:  is0(), one(), neg(), addTo(), mul(), mulBy(), recip()
 */

template<class A>
bool
HIntLib::matrixInverse (const A &a, typename A::type *m, int num)
{
   // This array records where each column can be found in the output matrix

   Array<int> permutation (num);

   typedef typename A::type T;

   for (int pos = 0; pos < num; ++pos)
   {
      {
         // search for pivot row

         int row = pos;

         while (a.is0 (m[row * num + pos]))  if (++row == num)  return false;

         permutation [pos] = row;

         // swap current row with pivot row

         if (row != pos)
         {
            for (int col = 0; col < num; ++col)
            {
               std::swap (m [pos * num + col], m [row * num + col]);
            }
         }
      }

      // update current row
      
      T pivot = a.recip (m[pos * num + pos]);
      m[pos * num + pos] = a.one();

      for (int col = 0; col < num; ++col)
      {
         a.mulBy (m[pos * num + col], pivot);
      }

      // update remaining rows

      for (int row = 0; row < num; ++row)
      {
         if (row == pos)  continue;

         T lc = a.neg (m [row * num + pos]);
         m[row * num + pos] = T();

         if (! a.is0 (lc))
         {
            for (int col = 0; col < num; ++col)
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

      for (int row = 0; row < num; ++row)
      {
         std::swap (m[row * num + col], m[row * num + permutation[col]]);
      }
   }
   
   return true;
}


#define HINTLIB_INSTANTIATE_LINEARALGEBRAGEN(X) \
   template bool isIdentityMatrix (const X&, const X::type*, int); \
   template bool isZeroMatrix (const X&, const X::type*, int, int); \
   template bool matrixInverse(const X&, X::type*, int); \
   template bool isLinearlyIndependent(const X&, X::type*, int, int);\
   template int matrixRank(const X&, X::type*, int, int);\
   template int numLinearlyIndependentVectors (const X&, X::type*, int, int); \
   template int nullSpace (const X&, X::type*, int, int, X::type*); \
   template int basisSupplement (const X&, X::type*, int, int);\
   template int basisSupplement (const X&, X::type*, int, int, int*);\
   template void matrixMul ( \
      const X&, const X::type*, const X::type*, int, int, int, X::type*); \
   template void matrixVectorMul ( \
      const X&, const X::type*, const X::type*, int, int, X::type*);\
   template void vectorMatrixMul ( \
      const X&, const X::type*, const X::type*, int, int, X::type*);

#define HINTLIB_INSTANTIATE_LINEARALGEBRAGEN_T(X) \
   template void matrixTranspose (const X*, int, int, X*);

#endif

