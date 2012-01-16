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


#include <HIntLib/generatormatrix.h>
#include <HIntLib/array.h>
#include <HIntLib/myalgorithm.h>
#include <HIntLib/exception.h>
#include <HIntLib/mymath.h>
#include <HIntLib/lookupfield.h>

#include <iomanip>

namespace L = HIntLib;

using std::min;
using std::max;

namespace
{

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
      
      while (++i != last)  if (*i & mask)  *i ^= pivot;

      // continue with smaller matrix
      
      mask <<= 1;
      --columns;
      ++first;
   }
   
   return rank;
}


/**
 *  linearlyIndependent()
 *
 *  Checks if a packed matrix over GF(2) has full row-rank.
 *
 *  rank must be larger that (last-first) -3 !!!
 *
 *  Each vector hold _columns_ elements from {0,1}, which have to be stored in
 *  the lower order bits.
 */

template<class Bi>
bool isLinearlyIndependent (Bi first, Bi last)
{
   typedef typename std::iterator_traits<Bi>::value_type T;

   switch (last - first)
   {
   case 0: return true;
   case 1: return *first != 0;
   case 2: return first[0] && first[1] && (first[0] ^ first[1]);
   case 3: return first[0] && first[1] && first[2] &&
                  (first[0] ^ first[1]) &&
                  (first[1] ^ first[2]) &&
                  (first[0] ^ first[2]) &&
                  (first[0] ^ first[1] ^ first[2]);
   }

   T mask (1);          // mask for current column
   
   while (last - first > 4)  // Any rows left?
   {
      // find pivot row
      
      Bi i = first;
          
      while ((*i & mask) == 0)  // if the coefficient == 0, we have to go on
      {
         if (++i == last)    // if we run out of rows, try a differnt column
         {
            mask <<= 1;  // try next column
            i = first;
         }
      }

      // exchange it with first row

      T pivot = *i;
      *i = *first;
      // *first = pivot; // we are not using *first again, so no need to update

      // subtract pivot from all remaining rows
      
      while (++i != last)  if (*i & mask)  *i ^= pivot;

      // continue with smaller matrix
      
      mask <<= 1;
      ++first;
   }

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


/**
 *  linearlyIndependent()
 *
 *  Checks if a packed matrix over GF(2) has full row-rank.
 *
 *  rank must be larger that (last-first) -3 !!!
 *
 *  Each vector hold _columns_ elements from {0,1}, which have to be stored in
 *  the lower order bits.
 */

template<class A>
bool isLinearlyIndependent (
      const A &a, typename A::type *m, unsigned numCols, unsigned numRows)
{
   typedef typename A::type T;

   unsigned col = 0;
   unsigned row = 0;
   
   while (row < numRows && col < numCols)  // Any rows left?
   {
#if 0
      // print

      cout << endl << row << "/" << col << endl;
      for (unsigned j = 0; j < numRows; ++j)
      {
         for (unsigned k = 0; k < numCols; ++k)
         {
            cout << setw(3) << unsigned (m [j *numCols + k]);
         }
         cout << endl;
      }
#endif

      // find pivot row
      
      unsigned i = row;
          
      while (a.is0 (m [i * numCols + col]))  // if the coefficient == 0, go on
      {
         if (++i == numRows)    // if we run out of rows, try a differnt column
         {
            i = row;

            if (++col == numCols)
            {
               // cout << "x"<< endl;
               return false;  // try next column
            }
         }
      }

      // subtract pivot from all remaining rows
      
      T recip = a.recip (m [i * numCols + col]);
      
      for (unsigned j = i + 1; j < numRows; ++j)
      {
         T lc = m [j * numCols + col];
         if (! a.is0 (lc))
         {
            T x = a.mul (lc, recip);
            for (unsigned k = col + 1; k < numCols; ++k)
            {
               a.subFrom (m [j * numCols + k], a.mul (x, m [i * numCols + k]));
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
   // cout << ((row==numRows) ? "ok" : "sorry") << endl;

   return row == numRows;
}


#if 0
/**
 *  isLinearlyIndependent()
 *
 *  determines if vectors over GF(2) are linearly independent by building all
 *  possible linear combinations.
 *
 *  This is O(2^n), whit n denoting the number of vectors!
 *
 *  Use rank() instead!
 */

template<class Bi>
bool isLinearlyIndependent (Bi first, Bi last)
{
   typedef typename std::iterator_traits<Bi>::value_type T;

   T x (0);

   for (T i = 0; i < (T(1) << (last-first)) - 1; ++i)
   {
      x ^= first[L::ls0(i)];

      if (! x)  return false;
   }

   return true;
}
#endif

}  // namespace


/***********  Base 2  ********************************************************/


template<class T>
int L::t_parameter (
   const L::GeneratorMatrix2<T> &gm, int lb, int ub, bool dimOpt)
{
   const unsigned DIM = gm.getDimension();
   const int M   = gm.getM();

   // theoretical lower bound:
   //
   // ceil (log_b (b*s - s + 1)) - 2 =
   //   ceil (log_b (s + 1)) - 2 =
   //   floor (log_b (s)) - 1 =
   //   ms1 (s) - 1

   lb = max (lb,
        max (0,
        max (M - int (gm.getTotalPrec()),
             min (M - 1, ms1 (DIM) - 1))));
   ub = min (ub, M);

   if (lb > ub)  throw InternalError (__FILE__, __LINE__);
   if (lb == ub)  return lb;

   Array<T> c (M * DIM);
   Array<T> selection (M-lb);
   Array<int> partition (DIM);

   // Copy Generator Matrix into c-array (containing tranposed matrices)

   for (unsigned d = 0; d < DIM; ++d)
   {
      for (int b = 0; b < M; ++b)
      {
         T x = 0;

         for (int r = 0; r < M; ++r)  x = (x << 1) + gm(d,r,b);

         c [d * M + b] = x;
      }
   }

   // try to select 1, 2, 3, ... M vectors (modulo ub/lb)
   
   for (int numVectors = M - ub + 1; numVectors <= M - lb; ++numVectors)
   {
      // choose these _numVectors_ vectors from c
      
      initial_partition (&partition[0], &partition[DIM], numVectors);

      // if we have checkt a low-dimensional submatrix already, at least one
      // vector has to come from the new dimension

      if (dimOpt)
      {
         ++partition [DIM-1];
         --partition [0];
      }

      do
      {
#if 0
         cout << "Partition: ";
         for (unsigned i = 0; i < DIM; ++i)  cout << partition[i] << " ";
         cout << endl;
#endif

         // copy the selected vectors
         // the last vectors for each dimension MUST be part of the linear
         // combination. We add them up and put them into selection[0].

         selection [0] = 0;
         unsigned l = 1;

         for (unsigned d = 0; d < DIM; ++d)
         {
            if (partition[d] > 0)
            {
               selection [0] ^= c [d * M + partition[d]-1];

               for (int i = 0; i < partition[d]-1; ++i)
               {
                  selection[l++] = c [d * M + i];
               }
            }
         }

#if 0
         cout << "Selected vectors: ";
         for (unsigned i = 0; i < l; ++i) cout << selection [i] << " ";
         cout << endl;
#endif
#if 0
         unsigned r = rank(&selection[0], &selection[l], M);
         if (r != l)  return M - numVectors + 1;
#else

         if (! isLinearlyIndependent (&selection[0], &selection[l]))
         {
            return M - numVectors + 1;
         }
#endif
      }
      while (next_partition (&partition[0], &partition[DIM]));
   }

   return lb;
}


template<class T>
int L::t_parameter2 (
   const L::GeneratorMatrix2<T> &gm, int lb, int ub, int maxRows)
{
   const int DIM = gm.getDimension();
   const int M   = gm.getM();

   if (maxRows > int (gm.getTotalPrec()))
   {
      throw InternalError (__FILE__, __LINE__);
   }

   lb = max (lb, 0);
   ub = min (ub, M);

   if (lb > ub)  throw InternalError (__FILE__, __LINE__);
   if (lb == ub)  return lb;

   Array<T> c (M * DIM);
   Array<T> selection (M-lb);
   Array<int> partition (DIM);

   // Copy Generator Matrix into c-array (containing tranposed matrices)

   for (int d = 0; d < DIM; ++d)
   {
      for (int b = 0; b < M; ++b)
      {
         T x = 0;

         for (int r = 0; r < M; ++r)  x = (x << 1) + gm(d,r,b);

         c [d * M + b] = x;
      }
   }

   // try to select 1, 2, 3, ... M vectors (modulo ub/lb)
   
   for (int numVectors = M - ub + 1; numVectors <= M - lb; ++numVectors)
   {
      if (numVectors > maxRows * DIM)  break;

      // choose these _numVectors_ vectors from c
      
      initial_partition (&partition[0], &partition[DIM], maxRows, numVectors);

#if 0
      // if we have checkt a low-dimensional submatrix already, at least one
      // vector has to come from the new dimension

      if (dimOpt)
      {
         ++partition [DIM-1];
         --partition [0];
      }
#endif

      do
      {
         // copy the selected vectors
         // the last vectors for each dimension MUST be part of the linear
         // combination. We add them up and put them into selection[0].

         selection [0] = 0;
         unsigned l = 1;

         for (int d = 0; d < DIM; ++d)
         {
            if (partition[d] > 0)
            {
               selection [0] ^= c [d * M + partition[d]-1];

               for (int i = 0; i < partition[d]-1; ++i)
               {
                  selection[l++] = c [d * M + i];
               }
            }
         }

         if (! isLinearlyIndependent (&selection[0], &selection[l]))
         {
            return M - numVectors + 1;
         }
      }
      while (next_partition (&partition[0], &partition[DIM], maxRows));
   }

   return lb;
}


/***********  General Case  **************************************************/


template<class T>
int L::t_parameter (
   const L::GeneratorMatrixGen<T> &gm, int lb, int ub, bool dimOpt)
{
   if (gm.getVectorization() != 1)  throw FIXME (__FILE__, __LINE__);

   const int DIM = gm.getDimension();
   const int M   = gm.getM();

   // theoretical lower bound:
   //
   // ceil (log_b (b*s - s + 1)) - 2 =
   //   ceil (log_b (s + 1)) - 2 =
   //   floor (log_b (s)) - 1

   lb = max (lb,
        max (0,
        max (M - int (gm.getTotalPrec()),
             min (M - 1,
                  int (logInt ((gm.getBase()-1) * DIM, gm.getBase()) - 1)))));
   ub = min (ub, M);

   if (lb > ub)  throw InternalError (__FILE__, __LINE__);
   if (lb == ub)  return lb;

   Array<T> selection ((M-lb) * M);
   Array<int> partition (DIM);
   LookupGaloisField<unsigned char> field (gm.getBase());

   // try to select 1, 2, 3, ... M vectors (modulo ub/lb)
   
   for (int numVectors = M - ub + 1; numVectors <= M - lb; ++numVectors)
   {
      // choose these _numVectors_ vectors from c
      
      initial_partition (&partition[0], &partition[DIM], numVectors);

      // if we have checked a low-dimensional submatrix already, at least one
      // vector must come from the new dimension

      if (dimOpt)
      {
         ++partition [DIM-1];
         --partition [0];
      }

      do
      {
         // copy the selected vectors
         // the last vectors for each dimension MUST be part of the linear
         // combination. We add them up and put them into selection[0].

         unsigned l = 0;

         for (int d = 0; d < DIM; ++d)
         {
            for (int i = 0; i < partition[d]; ++i)
            {
               for (int j = 0; j < M; ++j)
               {
                  selection [l * M + j] = gm (d, j, i);
               }
               ++l;
            }
         }

         if (! isLinearlyIndependent (field, &selection[0], M, l))
         {
            return M - numVectors + 1;
         }
      }
      while (next_partition (&partition[0], &partition[DIM]));
   }

   return lb;
}

template<class T>
int L::t_parameter2 (
   const L::GeneratorMatrixGen<T> &gm, int lb, int ub, int maxRows)
{
   if (gm.getVectorization() != 1)  throw FIXME (__FILE__, __LINE__);

   const int DIM = gm.getDimension();
   const int M   = gm.getM();

   if (maxRows > int (gm.getTotalPrec()))
   {
      throw InternalError (__FILE__, __LINE__);
   }

   lb = max (lb, 0);
   ub = min (ub, M);

   if (lb > ub)  throw InternalError (__FILE__, __LINE__);
   if (lb == ub)  return lb;

   Array<T> selection ((M-lb) * M);
   Array<int> partition (DIM);
   LookupGaloisField<unsigned char> field (gm.getBase());

   // try to select 1, 2, 3, ... M vectors (modulo ub/lb)
   
   for (int numVectors = M - ub + 1; numVectors <= M - lb; ++numVectors)
   {
      // choose these _numVectors_ vectors from c
      
      if (numVectors > maxRows * DIM)  break;

      initial_partition (&partition[0], &partition[DIM], maxRows, numVectors);

#if 0
      // if we have checked a low-dimensional submatrix already, at least one
      // vector must come from the new dimension

      if (dimOpt)
      {
         ++partition [DIM-1];
         --partition [0];
      }
#endif

      do
      {
         // copy the selected vectors
         // the last vectors for each dimension MUST be part of the linear
         // combination. We add them up and put them into selection[0].

         unsigned l = 0;

         for (int d = 0; d < DIM; ++d)
         {
            for (int i = 0; i < partition[d]; ++i)
            {
               for (int j = 0; j < M; ++j)
               {
                  selection [l * M + j] = gm (d, j, i);
               }
               ++l;
            }
         }

         if (! isLinearlyIndependent (field, &selection[0], M, l))
         {
            return M - numVectors + 1;
         }
      }
      while (next_partition (&partition[0], &partition[DIM], maxRows));
   }

   return lb;
}

namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template int t_parameter<> (const GeneratorMatrix2<X> &, int, int, bool); \
   template int t_parameter2<>(const GeneratorMatrix2<X> &, int, int, int);

   HINTLIB_INSTANTIATE (u32)
   HINTLIB_INSTANTIATE (u64)
#undef HINTLIB_INSTANTIATE

#define HINTLIB_INSTANTIATE(X) \
   template int t_parameter<> (const GeneratorMatrixGen<X> &, int, int, bool); \
   template int t_parameter2<>(const GeneratorMatrixGen<X> &, int, int, int);

   HINTLIB_INSTANTIATE (unsigned char)
#undef HINTLIB_INSTANTIATE
}

