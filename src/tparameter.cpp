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

#include <HIntLib/tparameter.h>
#include <HIntLib/generatormatrix2.h>
#include <HIntLib/array.h>
#include <HIntLib/myalgorithm.h>
#include <HIntLib/mymath.h>
#include <HIntLib/exception.h>
#include <HIntLib/lookupfield.h>

namespace L = HIntLib;

using std::min;
using std::max;

namespace
{

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


/*********************  linearly Independent ()  *****************************/


/**
 *  linearlyIndependent()
 *
 *  Checks if the given set of vectors over GF2 (packed into the bits of a
 *  word) are linearly independent.
 *
 *  The empty set is independent.
 */

template<class Bi>
bool isLinearlyIndependent (Bi first, Bi last)
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

   T mask (1);          // mask for current column
   
   while (last - first > 4)  // Any rows left?
   {
      // find pivot row
      
      Bi i = first;
          
      if (!*i)  return false;
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
      
      while (++i != last)
      {
         if (*i & mask)  if (! (*i ^= pivot)) return false;
      }

      // continue with smaller matrix
      
      mask <<= 1;
      ++first;
   }

   // exactly four vectors are left

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
 *  Determins of the set of _numRows_ vectors with _numCols_ entries from _A_
 *  is linearly independent.
 *
 *  The vectors are accessed useing _m_ [vector * _numcols + coord]
 */

template<class A>
bool isLinearlyIndependent (
      const A &a, typename A::type *m, unsigned numCols, unsigned numRows)
{
   typedef typename A::type T;

   unsigned col = 0;
   unsigned row = 0;
   
   while (row < numRows && col < numCols)  // Any rows/columns left?
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

   return row == numRows;
}


/***********  TCalc2  ********************************************************/


/**
 * TCalc2
 *
 * Contains the infrastructure for calculating the t-parameter in base b=2.
 *
 * This class is not multithreading-save.
 */

template<typename T>
class TCalc2
{
public:
   TCalc2 (const L::GeneratorMatrixGen<unsigned char> &, unsigned prec);
   // TCalc2 (const L::GeneratorMatrix2<T> &, unsigned prec);
   ~TCalc2 ();

   bool check (int thickness, bool dimOpt);
   bool checkTO (int thickness, bool dimOpt);
   bool checkRestricted (int thickness, int maxRows);
   bool checkRestrictedRO (int thickness, int maxRows);
   bool checkRestrictedTO (int thickness, int maxRows);

private:
   void init (const L::GeneratorMatrix &);

   const unsigned totalPrec;
   const int dim;

   T* c;
   int* partition;

   static const unsigned MAX_C_SIZE = 5000;
   static const unsigned MAX_PARTITION_SIZE = 200;

   static T selection [std::numeric_limits<T>::digits];
   static T staticC [MAX_C_SIZE];
   static int staticPartition [MAX_PARTITION_SIZE];
};


/**
 *  static data
 */

template<typename T> T   TCalc2<T>::selection [];
template<typename T> T   TCalc2<T>::staticC [];
template<typename T> int TCalc2<T>::staticPartition [];

/**
 *  init ()
 */

template<typename T>
void TCalc2<T>::init (const L::GeneratorMatrix &gm)
{
   if (gm.getTotalPrec() < totalPrec)  throw L::FIXME (__FILE__, __LINE__);

   c = (totalPrec * dim > MAX_C_SIZE) ? new T [totalPrec * dim] : staticC;

   partition = (dim > int (MAX_PARTITION_SIZE))
        ?  new int [dim] : staticPartition;
}

template<typename T>
TCalc2<T>::~TCalc2 ()
{
   if (c != staticC)  delete[] c;
   if (partition != staticPartition)  delete[] partition;
}

/**
 *  Constructor
 */

#if 0
template<typename T>
TCalc2<T>::TCalc2 (const L::GeneratorMatrix2<T> &gm, unsigned prec)
   : totalPrec (prec),
     dim (gm.getDimension())
{
   init (gm);

   // Copy Generator Matrix into c-array (containing tranposed matrices)

   for (int d = 0; d < dim; ++d)
   {
      for (unsigned b = 0; b < totalPrec; ++b)
      {
         T x = 0;

         for (unsigned r = 0; r < gm.getM(); ++r)  x = (x << 1) | gm(d,r,b);

         c [d * totalPrec + b] = x;
      }
   }
}
#endif

template<typename T>
TCalc2<T>::TCalc2 (
   const L::GeneratorMatrixGen<unsigned char> &gm, unsigned prec)
   : totalPrec (prec),
     dim (gm.getDimension())
{
   init (gm);

   // Copy Generator Matrix into c-array (containing tranposed matrices)

   for (int d = 0; d < dim; ++d)
   {
      for (unsigned b = 0; b < totalPrec; ++b)
      {
         T x = 0;

         for (unsigned r = 0; r < gm.getM(); ++r)  x = (x << 1) | gm(d,r,b);

         c [d * totalPrec + b] = x;
      }
   }
}

/**
 *  check()
 *
 *  Checks if (at least) a certain thickness is present.
 *
 *  If _dimOpt_ is true, only partitions containing a vector from the last
 *  matrix are checked.
 */

template<typename T>
bool TCalc2<T>::check (int thickness, bool dimOpt)
{
   // choose these _thickness_ vectors from c
      
   L::initial_partition (&partition[0], &partition[dim], thickness);

   // if we have checkt a low-dimensional submatrix already, at least one
   // vector has to come from the new dimension

   if (dimOpt)
   {
      partition [dim-1] = 1;
      --partition [0];
   }

   do
   {
      // copy the selected vectors

      unsigned l = 0;

      for (int d = 0; d < dim; ++d)
      {
         for (int i = 0; i < partition[d]; ++i)
         {
            selection[l++] = c [d * totalPrec + i];
         }
      }

      if (! isLinearlyIndependent (&selection[0], &selection[l]))
      {
         return false;
      }
   }
   while (L::next_partition (&partition[0], &partition[dim]));

   return true;
}


/**
 *  checkTO ()
 *
 *  Checks if (at least) a certain thickness is present.
 *
 *  The check is performed under the assumption that a thickness of
 *  _thickness_-1 is present.
 *
 *  If _dimOpt_ is true, only partitions containing a vector from the last
 *  matrix are checked.
 */

template<typename T>
bool TCalc2<T>::checkTO (int thickness, bool dimOpt)
{
   // choose these _thickness_ vectors from c
      
   L::initial_partition (&partition[0], &partition[dim], thickness);

   // if we have checkt a low-dimensional submatrix already, at least one
   // vector has to come from the new dimension

   if (dimOpt)
   {
      partition [dim-1] = 1;
      --partition [0];
   }

   do
   {
      // copy the selected vectors
      // the last vectors for each dimension MUST be part of the linear
      // combination. We add them up and put them into selection[0].

      selection [0] = 0;
      unsigned l = 1;

      for (int d = 0; d < dim; ++d)
      {
         if (partition[d] > 0)
         {
            selection [0] ^= c [d * totalPrec + partition[d]-1];

            for (int i = 0; i < partition[d]-1; ++i)
            {
               selection[l++] = c [d * totalPrec + i];
            }
         }
      }

      if (! isLinearlyIndependent (&selection[0], &selection[l]))
      {
         return false;
      }
   }
   while (L::next_partition (&partition[0], &partition[dim]));

   return true;
}


/**
 *  checkRestrictedT0
 *
 *  Checks if (at least) a certain thickness is present.  Only partitions with
 *  at most _maxRows_ per matrix are considered.
 *
 *  The check is performed under the assumption that a thickness of
 *  _thickness_-1 is present.
 *
 *  If _dimOpt_ is true, only partitions containing a vector from the last
 *  matrix are checked.
 */

template<typename T>
bool TCalc2<T>::checkRestrictedTO (int thickness, int maxRows)
{
   if (thickness > maxRows * dim)  return true;

   // choose these _numVectors_ vectors from c
   
   L::initial_partition (&partition[0], &partition[dim], maxRows, thickness);

   do
   {
      // copy the selected vectors
      // the last vectors for each dimension MUST be part of the linear
      // combination. We add them up and put them into selection[0].

      selection [0] = 0;
      unsigned l = 1;

      for (int d = 0; d < dim; ++d)
      {
         if (partition[d] > 0)
         {
            selection [0] ^= c [d * totalPrec + partition[d]-1];

            for (int i = 0; i < partition[d]-1; ++i)
            {
               selection[l++] = c [d * totalPrec + i];
            }
         }
      }

      if (! isLinearlyIndependent (&selection[0], &selection[l]))
      {
         return false;
      }
   }
   while (L::next_partition (&partition[0], &partition[dim], maxRows));

   return true;
}

/**
 *  check()
 *
 *  Checks if (at least) a certain thickness is present. Only partitions with
 *  at least _maxRows_ per matrix considered.
 *
 *  If _dimOpt_ is true, only partitions containing a vector from the last
 *  matrix are checked.
 */

template<typename T>
bool TCalc2<T>::checkRestricted (int thickness, int maxRows)
{
   if (thickness > maxRows * dim)  return true;;

   // choose these _numVectors_ vectors from c
   
   L::initial_partition (&partition[0], &partition[dim], maxRows, thickness);

   do
   {
      // copy the selected vectors

      unsigned l = 0;

      for (int d = 0; d < dim; ++d)
      {
         for (int i = 0; i < partition[d]; ++i)
         {
            selection [l++] = c [d * totalPrec + i];
         }
      }

      if (! isLinearlyIndependent (&selection[0], &selection[l]))
      {
         return false;
      }
   }
   while (L::next_partition (&partition[0], &partition[dim], maxRows));

   return true;
}

/**
 *  check()
 *
 *  Checks if (at least) a certain thickness is present. Only partitions with
 *  at most _maxRows_ per matrix considered.
 *
 *  The check is performed under the assumption that the thickness is present
 *  considerung at most _maxRows_-1 per matrix. 
 *
 *  If _dimOpt_ is true, only partitions containing a vector from the last
 *  matrix are checked.
 */

template<typename T>
bool TCalc2<T>::checkRestrictedRO (int thickness, int maxRows)
{
   if (thickness > maxRows * dim)  return true;;

   // choose these _numVectors_ vectors from c
   
   L::initial_partition (&partition[0], &partition[dim], maxRows, thickness);

   do
   {
      // copy the selected vectors

      unsigned l = 0;
      bool fullColumn = false;

      for (int d = 0; d < dim; ++d)
      {
         if (partition[d] == maxRows)  fullColumn = true;

         for (int i = 0; i < partition[d]; ++i)
         {
            selection [l++] = c [d * totalPrec + i];
         }
      }

      if (! fullColumn)  continue;

      if (! isLinearlyIndependent (&selection[0], &selection[l]))
      {
         return false;
      }
   }
   while (L::next_partition (&partition[0], &partition[dim], maxRows));

   return true;
}


/***********  T Calc Gen  ****************************************************/


class TCalcGen
{
public:
   TCalcGen (const L::GeneratorMatrixGen<unsigned char> &);

   bool check (int thickness, bool dimOpt);
   bool checkRestricted (int thickness, int maxRows);
   bool checkRestrictedRO (int thickness, int maxRows);

private:
   const L::GeneratorMatrixGen<unsigned char> & gm;
   L::Array<unsigned char> selection;
   L::Array<int> partition;
   L::LookupGaloisField<unsigned char> field;
};

TCalcGen::TCalcGen (const L::GeneratorMatrixGen<unsigned char> &_gm)
   : gm (_gm),
     selection (L::sqr (_gm.getM())),
     partition (_gm.getDimension()),
     field (_gm.getBase())
{}

bool TCalcGen::check (int thickness, bool dimOpt)
{
   const int M   = gm.getM();
   const int DIM = gm.getDimension();
      
   L::initial_partition (&partition[0], &partition[DIM], thickness);

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

      if (! isLinearlyIndependent (field, &selection[0], M, l))  return false;
   }
   while (L::next_partition (&partition[0], &partition[DIM]));

   return true;
}

bool TCalcGen::checkRestricted (int thickness, int maxRows)
{
   const int M   = gm.getM();
   const int DIM = gm.getDimension();
   
   if (thickness > maxRows * DIM)  return true;

   L::initial_partition (&partition[0], &partition[DIM], maxRows, thickness);

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

      if (! isLinearlyIndependent (field, &selection[0], M, l))  return false;
   }
   while (L::next_partition (&partition[0], &partition[DIM], maxRows));

   return true;
}

bool TCalcGen::checkRestrictedRO (int thickness, int maxRows)
{
   const int M   = gm.getM();
   const int DIM = gm.getDimension();
   
   if (thickness > maxRows * DIM)  return true;

   L::initial_partition (&partition[0], &partition[DIM], maxRows, thickness);

   do
   {
      // copy the selected vectors
      // the last vectors for each dimension MUST be part of the linear
      // combination. We add them up and put them into selection[0].

      unsigned l = 0;
      bool fullColumn = false;

      for (int d = 0; d < DIM; ++d)
      {
         if (partition[d] == maxRows)  fullColumn = true;

         for (int i = 0; i < partition[d]; ++i)
         {
            for (int j = 0; j < M; ++j)
            {
               selection [l * M + j] = gm (d, j, i);
            }
            ++l;
         }
      }

      if (! fullColumn)  continue;

      if (! isLinearlyIndependent (field, &selection[0], M, l))  return false;
   }
   while (L::next_partition (&partition[0], &partition[DIM], maxRows));

   return true;
}

   inline
   bool use32bit (const L::GeneratorMatrixGen<unsigned char> &gm)
   {
      return gm.getBase() == 2
          && int(gm.getM()) <= std::numeric_limits<L::u32>::digits;
   }

}  // anonymous namespace


/*****************  Public routines  *****************************************/


bool L::confirmT (
   const L::GeneratorMatrixGen<unsigned char> &gm, int t, TOption opts)
{
   if (opts & LOWER_RESTRICTION_OK)  throw FIXME (__FILE__, __LINE__);

   if (t < 0 || t > int (gm.getM()))  throw FIXME (__FILE__, __LINE__);
   if (gm.getM() - t > gm.getTotalPrec())  return false;

   if (use32bit (gm))
   {
      TCalc2<u32> calc (gm, gm.getM() - t);

      if (opts & LARGER_T_OK)
         return calc.check (gm.getM() - t, opts & LOWER_DIM_OK);
      else
         return calc.checkTO (gm.getM() - t, opts & LOWER_DIM_OK);
   }
   else
   {
      TCalcGen calc (gm);
      return calc.check (gm.getM() - t, opts & LOWER_DIM_OK);
   }
}

bool L::confirmTRestricted (
   const L::GeneratorMatrixGen<unsigned char> &gm, int t, int maxRows,
   TOption opts)
{
   // missing: LARGER_T_OK, LOWER_DIM_OK

   if (t < 0 || t > int (gm.getM()))  throw FIXME (__FILE__, __LINE__);
   if (maxRows > int (gm.getTotalPrec()))  throw FIXME (__FILE__, __LINE__);

   if (use32bit (gm))
   {
      TCalc2<u32> calc (gm, std::min (unsigned (maxRows), gm.getM() - t));

      if (opts & LOWER_RESTRICTION_OK)
         return calc.checkRestrictedRO (gm.getM() - t, maxRows);
      else
         return calc.checkRestricted   (gm.getM() - t, maxRows);
   }
   else
   {
      TCalcGen calc (gm);

      if (opts & LOWER_RESTRICTION_OK)
         return calc.checkRestrictedRO (gm.getM() - t, maxRows);
      else
         return calc.checkRestricted   (gm.getM() - t, maxRows);
   }
}

int L::tParameter (
   const L::GeneratorMatrixGen<unsigned char> &gm, int lb, int ub, TOption opts)
{
   if (opts & (LOWER_RESTRICTION_OK | LARGER_T_OK))
   {
      throw FIXME (__FILE__, __LINE__);
   }

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
                  int (L::logInt ((gm.getBase()-1) * gm.getDimension(),
                                  gm.getBase()) - 1)))));
   ub = min (ub, M);

   if (lb > ub)  throw InternalError (__FILE__, __LINE__);
   if (lb == ub)  return lb;

   if (use32bit (gm))
   {
      TCalc2<u32> calc (gm, M - lb);
   
      for (int numVectors = M - ub + 1; numVectors <= M - lb; ++numVectors)
      {
         if (! calc.checkTO (numVectors, opts & LOWER_DIM_OK))
         {
            return M - numVectors + 1;
         }
      }
   }
   else
   {
      TCalcGen calc (gm);
   
      for (int numVectors = M - ub + 1; numVectors <= M - lb; ++numVectors)
      {
         if (! calc.check (numVectors, opts & LOWER_DIM_OK))
         {
            return M - numVectors + 1;
         }
      }
   }

   return lb;
}

int L::tParameterRestricted (
   const L::GeneratorMatrixGen<unsigned char> &gm,
   int lb, int ub, int maxRows, TOption opts)
{
   if (opts)  throw FIXME (__FILE__, __LINE__);

   const int M = gm.getM();

   if (maxRows > int (gm.getTotalPrec()))
   {
      throw InternalError (__FILE__, __LINE__);
   }

   lb = max (lb, 0);
   ub = min (ub, M);

   if (lb > ub)  throw InternalError (__FILE__, __LINE__);
   if (lb == ub)  return lb;

   if (use32bit (gm))
   {
      TCalc2<u32> calc (gm, std::min (maxRows, M - lb));
   
      for (int numVectors = M - ub + 1; numVectors <= M - lb; ++numVectors)
      {
         if (! calc.checkRestrictedTO (numVectors, maxRows))
         {
            return M - numVectors + 1;
         }
      }
   }
   else
   {
      TCalcGen calc (gm);
   
      for (int numVectors = M - ub + 1; numVectors <= M - lb; ++numVectors)
      {
         if (! calc.checkRestricted (numVectors, maxRows))
         {
            return M - numVectors + 1;
         }
      }
   }

   return lb;
}


