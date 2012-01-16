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

#include <memory>

#include <HIntLib/tparameter.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/generatormatrix2row.h>
#include <HIntLib/generatormatrixgen.h>
#include <HIntLib/generatormatrixvirtual.h>
#include <HIntLib/minmaxfinder.h>
#include <HIntLib/hlalgorithm.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;

using std::min;
using std::max;

namespace
{
   typedef L::GeneratorMatrixGen<unsigned char> GMGen;

   inline
   bool use32bit (const GMGen &gm)
   {
      return gm.getBase() == 2
          && int(gm.getM()) <= std::numeric_limits<L::u32>::digits;
   }

   inline
   bool use64bit (const GMGen &gm)
   {
      return gm.getBase() == 2
          && int(gm.getM()) <= std::numeric_limits<L::u64>::digits;
   }

}  // anonymous namespace


/***********  TCalc  *********************************************************/


/**
 *  Constructor and destructor
 */

int L::TCalc::staticPartition [];

L::TCalc::TCalc (int _dim)
   : dim (_dim),
     partition ((dim > int (MAX_PARTITION_SIZE))
             ?  new int [dim] : staticPartition)
{}

L::TCalc::~TCalc ()
{
   if (partition != staticPartition)  delete[] partition;
}


/***********  TCalc2  ********************************************************/


/**
 *  static data
 */

template<typename T> T L::TCalc2<T>::selection [];

/**
 *  check()
 *
 *  Checks if (at least) a certain strength is present.
 *
 *  If _dimOpt_ is true, only partitions containing a vector from the last
 *  matrix are checked.
 */

template<typename T>
bool L::TCalc2<T>::check (int strength, bool dimOpt)
{
   // choose these _strength vectors from c

   initial_partition (&partition[0], &partition[dim], strength);

   // if we have checkt a low-dimensional submatrix already, at least one
   // vector has to come from the new dimension

   if (dimOpt)
   {
      partition [dim-1] = 1;
      --partition [0];
   }

   do
   {
      copyToSelection ();
      if (singular (strength))  return false;
   }
   while (next_partition (&partition[0], &partition[dim]));

   return true;
}


/**
 *  checkTO ()
 *
 *  Checks if (at least) a certain strength is present.
 *
 *  The check is performed under the assumption that a strength of
 *  _strength-1 is present.
 *
 *  If _dimOpt_ is true, only partitions containing a vector from the last
 *  matrix are checked.
 */

template<typename T>
bool L::TCalc2<T>::checkTO (int strength, bool dimOpt)
{
   // choose these _strength vectors from c

   initial_partition (&partition[0], &partition[dim], strength);

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
      T* l = &selection[1];

      for (int d = 0; d < dim; ++d)
      {
         if (partition[d] > 0)
         {
            selection [0] ^= gm (d, partition[d]-1);

            copyToSelection (d, partition[d] - 1, l);
         }
      }

      if (singular (l))  return false;
   }
   while (next_partition (&partition[0], &partition[dim]));

   return true;
}


/**
 *  checkRestricted()
 *
 *  Checks if (at least) a certain strength is present. Only partitions with
 *  at most _maxRows_ per matrix considered.
 */

// #include <iostream>
template<typename T>
bool L::TCalc2<T>::checkRestricted (int strength, int maxRows)
{
   if (strength > maxRows * dim)  return true;;

   // choose these _numVectors_ vectors from c

   initial_partition (&partition[0], &partition[dim], maxRows, strength);

   do
   {
      copyToSelection ();

      if (singular (strength))
      {
#if 0
         for (int i = 0; i < dim; ++i)  std::cout << partition[i];
         std::cout <<'\n';
#endif
         
         return false;
      }
   }
   while (next_partition (&partition[0], &partition[dim], maxRows));

   return true;
}


/**
 *  checkRestrictedRO()
 *
 *  Checks if (at least) a certain strength is present. Only partitions with
 *  at most _maxRows_ per matrix considered.
 *
 *  The check is performed under the assumption that the strength is present
 *  considerung at most _maxRows_-1 per matrix.
 */

template<typename T>
bool L::TCalc2<T>::checkRestrictedRO (int strength, int maxRows)
{
   if (strength > maxRows * dim)  return true;;

   // choose these _numVectors_ vectors from c

   initial_partition (&partition[0], &partition[dim], maxRows, strength);

   do
   {
      // Make sure we have a matrix with maxRows rows

      if (! partitionContains (maxRows))  continue;

      copyToSelection ();

      if (singular (strength))  return false;
   }
   while (next_partition (&partition[0], &partition[dim], maxRows));

   return true;
}


/**
 *  checkRestrictedT0
 *
 *  Checks if (at least) a certain strength is present.  Only partitions with
 *  at most _maxRows_ per matrix are considered.
 *
 *  The check is performed under the assumption that a strength of
 *  _strength-1 is present.
 */

template<typename T>
bool L::TCalc2<T>::checkRestrictedTO (int strength, int maxRows)
{
   if (strength > maxRows * dim)  return true;

   // choose these _numVectors_ vectors from c

   initial_partition (&partition[0], &partition[dim], maxRows, strength);

   do
   {
      // copy the selected vectors
      // the last vectors for each dimension MUST be part of the linear
      // combination. We add them up and put them into selection[0].

      selection [0] = 0;
      T* l = &selection[1];

      for (int d = 0; d < dim; ++d)
      {
         if (partition[d] > 0)
         {
            selection [0] ^= gm (d, partition[d]-1);
            copyToSelection (d, partition[d] - 1, l);
         }
      }

      if (singular (l))  return false;
   }
   while (next_partition (&partition[0], &partition[dim], maxRows));

   return true;
}


/**
 *  checkRestrictedT0RO
 *
 *  Checks if (at least) a certain strength is present.  Only partitions with
 *  at most _maxRows_ per matrix are considered.
 *
 *  The check is performed under the assumption that a strength of
 *  _strength-1 is present.
 *
 *  The check is performed under the assumption that the strength is present
 *  considerung at most _maxRows_-1 per matrix.
 */

template<typename T>
bool L::TCalc2<T>::checkRestrictedTORO (int strength, int maxRows)
{
   if (strength > maxRows * dim)  return true;

   // choose these _numVectors_ vectors from c

   initial_partition (&partition[0], &partition[dim], maxRows, strength);

   do
   {
      // Make sure we have a matrix with maxRows rows

      if (! partitionContains (maxRows))  continue;

      // copy the selected vectors
      // the last vectors for each dimension MUST be part of the linear
      // combination. We add them up and put them into selection[0].

      selection [0] = 0;
      T* l = &selection[1];

      for (int d = 0; d < dim; ++d)
      {
         if (partition[d] > 0)
         {
            selection [0] ^= gm (d, partition[d]-1);
            copyToSelection (d, partition[d] - 1, l);
         }
      }

      if (singular (l))  return false;
   }
   while (next_partition (&partition[0], &partition[dim], maxRows));

   return true;
}


/***********  T Calc Gen  ****************************************************/


/**
 *  static data
 */

unsigned char L::TCalcGen::staticSelection [];


/**
 *  Constructor
 */

L::TCalcGen::TCalcGen (const GMGen &_gm)
   : TCalc(_gm.getDimension()), gm (_gm), M (_gm.getM()),
     la (LinearAlgebra::make (_gm.getBase()))
{
   selection = (sqr(_gm.getM()) > MAX_SELECTION_SIZE)
      ? new unsigned char [sqr(M)] : staticSelection;
}

L::TCalcGen::~TCalcGen ()
{
   if (selection != staticSelection)  delete[] selection;
   delete la;
}


/**
 *  copy To Selection ()
 */

inline
void
L::TCalcGen::copyToSelection (int d, int num, unsigned& pos)
{
   for (int i = 0; i < num; ++i)
   {
      for (int j = 0; j < M; ++j)  selection [pos * M + j] = gm (d, j, i);
      ++pos;
   }
}

inline
void
L::TCalcGen::copyToSelection ()
{
   unsigned pos = 0;
   for (int d = 0; d < dim; ++d)  copyToSelection (d, partition[d], pos);
}


/**
 *  check()
 */

bool L::TCalcGen::check (int strength, bool dimOpt)
{
   initial_partition (&partition[0], &partition[dim], strength);

   // if we have checked a low-dimensional submatrix already, at least one
   // vector must come from the new dimension

   if (dimOpt)
   {
      ++partition [dim-1];
      --partition [0];
   }

   do
   {
      copyToSelection ();
      if (singular (strength))  return false;
   }
   while (next_partition (&partition[0], &partition[dim]));

   return true;
}


/**
 *  checkRestricted()
 *
 *  Checks if (at least) a certain strength is present. Only partitions with
 *  at most _maxRows_ per matrix considered.
 */

bool L::TCalcGen::checkRestricted (int strength, int maxRows)
{
   if (strength > maxRows * dim)  return true;

   initial_partition (&partition[0], &partition[dim], maxRows, strength);

   do
   {
      copyToSelection ();
      if (singular (strength))  return false;
   }
   while (next_partition (&partition[0], &partition[dim], maxRows));

   return true;
}


/**
 *  checkRestrictedRO()
 *
 *  Checks if (at least) a certain strength is present. Only partitions with
 *  at most _maxRows_ per matrix considered.
 *
 *  The check is performed under the assumption that the strength is present
 *  considerung at most _maxRows_-1 per matrix.
 */

bool L::TCalcGen::checkRestrictedRO (int strength, int maxRows)
{
   if (strength > maxRows * dim)  return true;

   // Partition  strength-maxRows  elements on  dim-1  stacks.
   // Then, add an additional stack with  maxRows  elements.
   // Adding this extra stack ensures that at least one stack with  maxRows
   //   elements is present.
   // Adding the extra stack only left of any other stacks that contain  maxRows
   //   elements ensures that no combination is created twice.

   initial_partition (
         &partition[0], &partition[dim-1], maxRows, strength - maxRows);

   do
   {
      for (int pos = 0; pos < dim; ++pos)
      {
         if (pos > 0 && partition [pos-1] == maxRows)  break;

         // create selection

         unsigned l = 0;

         for (int d = 0; d < pos; ++d)  copyToSelection (d, partition[d], l);
         copyToSelection (pos, maxRows, l);
         for (int d = pos + 1; d < dim; ++d)
            copyToSelection (d, partition[d-1], l);

         if (singular (strength))
         {
#if 0
            for (int d = 0; d < pos; ++d)  std::cout << partition[d];
            std::cout << maxRows;
            for (int d = pos + 1; d < dim; ++d)std::cout << partition[d];
            std::cout << '\n';
#endif
            return false;
         }
      }
   }
   while (next_partition (&partition[0], &partition[dim-1], maxRows));

   return true;
}


/*****************  Public routines  *****************************************/


bool L::confirmT (const GeneratorMatrix& gm, int t, TOption opts)
{
   if (opts & LOWER_RESTRICTION_OK)  throw FIXME (__FILE__, __LINE__);

   if (t < 0 || t > int (gm.getM()))  throw FIXME (__FILE__, __LINE__);
   const unsigned k = gm.getM() - t;
   if (k > gm.getPrec())  return false;

   AdjustPrec gm1 (k, gm);

   if (use32bit (gm))
   {
      GeneratorMatrix2Row<u32> gm2 (gm1);
      TCalc2<u32> calc (gm2);

      if (opts & LARGER_T_OK)
         return calc.checkTO (k, opts & LOWER_DIM_OK);
      else
         return calc.check (k, opts & LOWER_DIM_OK);
   }
   else if (use64bit (gm))
   {
      GeneratorMatrix2Row<u64> gm2 (gm1);
      TCalc2<u64> calc (gm2);

      if (opts & LARGER_T_OK)
         return calc.checkTO (k, opts & LOWER_DIM_OK);
      else
         return calc.check (k, opts & LOWER_DIM_OK);
   }
   else
   {
      GMGen gm2 (gm1);
      TCalcGen calc (gm2);
      return calc.check (k, opts & LOWER_DIM_OK);
   }
}


bool L::confirmTRestricted (
      const GeneratorMatrix& gm, int t, int maxRows, TOption opts)
{
   if (opts & LOWER_DIM_OK)  throw FIXME (__FILE__, __LINE__);

   if (t < 0 || t > int (gm.getM()))  throw FIXME (__FILE__, __LINE__);
   if (maxRows > int (gm.getPrec()))  throw FIXME (__FILE__, __LINE__);

   const int strength = gm.getM() - t;
   const int numRows = std::min (maxRows, strength);

   AdjustPrec gm1 (numRows, gm);

   if (use32bit (gm))
   {
      GeneratorMatrix2Row<u32> gm2 (gm1);
      TCalc2<u32> calc (gm2);

      if ((opts & LOWER_RESTRICTION_OK) && strength > 1)
      {
         if (opts & LARGER_T_OK)
            return calc.checkRestrictedTORO (strength, maxRows);
         else
            return calc.checkRestrictedRO   (strength, maxRows);
      }
      else
      {
         if (opts & LARGER_T_OK)
            return calc.checkRestrictedTO (strength, maxRows);
         else
            return calc.checkRestricted   (strength, maxRows);
      }
   }
   else if (use64bit (gm))
   {
      GeneratorMatrix2Row<u64> gm2 (gm1);
      TCalc2<u64> calc (gm2);

      if ((opts & LOWER_RESTRICTION_OK) && strength > 1)
      {
         if (opts & LARGER_T_OK)
            return calc.checkRestrictedTORO (strength, maxRows);
         else
            return calc.checkRestrictedRO   (strength, maxRows);
      }
      else
      {
         if (opts & LARGER_T_OK)
            return calc.checkRestrictedTO (strength, maxRows);
         else
            return calc.checkRestricted   (strength, maxRows);
      }
   }
   else
   {
      GMGen gm2 (gm1);
      TCalcGen calc (gm2);

      // there is no way to optimize for LARGER_T_OK in base > 2

      if ((opts & LOWER_RESTRICTION_OK) && strength > 1)
      {
         return calc.checkRestrictedRO (strength, maxRows);
      }
      else
         return calc.checkRestricted   (strength, maxRows);
   }
}


int L::tParameter (const GeneratorMatrix& gm, int lb, int ub, TOption opts)
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
        max (M - int (gm.getPrec()),
             min (M - 1,
                  int (L::logInt ((gm.getBase()-1) * gm.getDimension(),
                                  gm.getBase()) - 1)))));
   ub = min (ub, M);

   if (lb > ub)  throw InternalError (__FILE__, __LINE__);
   if (lb == ub)  return lb;

   AdjustPrec gm1 (M - lb, gm);

   if (use32bit (gm1))
   {
      GeneratorMatrix2Row<u32> gm2 (gm1);
      TCalc2<u32> calc (gm2);

      for (int numVectors = M - ub + 1; numVectors <= M - lb; ++numVectors)
      {
         if (! calc.checkTO (numVectors, opts & LOWER_DIM_OK))
         {
            return M - numVectors + 1;
         }
      }
   }
   else if (use64bit (gm1))
   {
      GeneratorMatrix2Row<u64> gm2 (gm1);
      TCalc2<u64> calc (gm2);

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
      GMGen gm2 (gm1);
      TCalcGen calc (gm2);

      while (lb != ub)
      {
         int middle = (ub + lb) / 2;
         if (calc.check (M - middle, opts & LOWER_DIM_OK))
         {
            ub = middle;
         }
         else
         {
            lb = middle + 1;
         }
      }
   }

   return lb;
}


int L::tParameterRestricted (
      const GeneratorMatrix& gm, int lb, int ub, int maxRows, TOption opts)
{
   if (opts)  throw FIXME (__FILE__, __LINE__);

   const int M = gm.getM();

   if (maxRows > int (gm.getPrec()))
   {
      throw InternalError (__FILE__, __LINE__);
   }

   lb = max (lb, 0);
   ub = min (ub, M);

   if (lb > ub)  throw InternalError (__FILE__, __LINE__);
   if (lb == ub)  return lb;

   AdjustPrec gm1 (M - lb, gm);

   if (use32bit (gm1))
   {
      GeneratorMatrix2Row<u32> gm2 (gm1);
      TCalc2<u32> calc (gm2);

      for (int numVectors = M - ub + 1; numVectors <= M - lb; ++numVectors)
      {
         if (! calc.checkRestrictedTO (numVectors, maxRows))
         {
            return M - numVectors + 1;
         }
      }
   }
   else if (use64bit (gm1))
   {
      GeneratorMatrix2Row<u64> gm2 (gm1);
      TCalc2<u64> calc (gm2);

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
      GMGen gm2 (gm1);
      TCalcGen calc (gm2);

      while (lb != ub)
      {
         int middle = (ub + lb) / 2;
         if (calc.checkRestricted (M - middle, maxRows))
         {
            ub = middle;
         }
         else
         {
            lb = middle + 1;
         }
      }
   }

   return lb;
}


/*************** t-Parameter of low-dimensional projections ******************/


/**
 *  tParameter1DimProjection()
 */

int L::tParameter1DimProjection (const GeneratorMatrix& gm, unsigned d)
{
   const unsigned m = gm.getM();
   const unsigned prec = gm.getPrec();

   Array<unsigned char> matrix (m * prec);

   for (unsigned r = 0; r < m; ++r)
   {
      for (unsigned b = 0; b < prec; ++b)
      {
         matrix [r + m * b] = gm.getDigit (d, r, b);
      }
   }
   
   std::auto_ptr<LinearAlgebra> la (LinearAlgebra::make (gm.getBase()));

   return m - la->numLinearlyIndependentVectors (matrix.begin(), prec, m);
}


/**
 *  tParameterMax1DimProjection()
 */

int L::tParameterMax1DimProjection (const GeneratorMatrix& gm)
{
   const unsigned dim = gm.getDimension();
   const unsigned m = gm.getM();
   const unsigned prec = gm.getPrec();

   Array<unsigned char> matrix (m * prec);
   MinFinder<unsigned> mf;

   for (unsigned d = 0; d < dim; ++d)
   {
      for (unsigned r = 0; r < m; ++r)
      {
         for (unsigned b = 0; b < prec; ++b)
         {
            matrix [r + m * b] = gm.getDigit (dim, r, b);
         }
      }
   
      std::auto_ptr<LinearAlgebra> la (LinearAlgebra::make (gm.getBase()));

      mf << la->numLinearlyIndependentVectors (matrix.begin(), prec, m);
   }

   return m - mf.getMinimum();
}


/**
 *  tParameter2DimProjection()
 */

int L::tParameter2DimProjection
      (const GeneratorMatrix& gm, unsigned d1, unsigned d2)
{
   unsigned dimensions [2];
   dimensions [0] = d1;
   dimensions [1] = d2;

   SelectDimensions gm2 (dimensions, dimensions + 2, gm);

   GeneratorMatrixGen<unsigned char> gm3 (gm2);

   return tParameter (gm3);
}


/**
 *  tParameterMax2DimProjection()
 */

int L::tParameterMax2DimProjection (const GeneratorMatrix& gm)
{
   const unsigned dim = gm.getDimension();

   if (dim < 2)  return tParameterMax1DimProjection (gm);

   MaxFinder<unsigned> mf;

   SelectDimensions gm2 (2, gm);

   for (unsigned d1 = 1; d1 < dim; ++d1)
   {
      gm2.selectDimension (0, d1);

      for (unsigned d2 = 0; d2 < d1; ++d2)
      {
         gm2.selectDimension (1, d2);
         mf << tParameter (gm2);
      }
   }

   return mf.getMaximum();
}


/**
 *  tParameterMax3DimProjection()
 */

int L::tParameterMax3DimProjection (const GeneratorMatrix& gm)
{
   const unsigned dim = gm.getDimension();

   if (dim < 3)  return tParameterMax2DimProjection (gm);

   MaxFinder<unsigned> mf;

   SelectDimensions gm2 (3, gm);

   for (unsigned d1 = 2; d1 < dim; ++d1)
   {
      gm2.selectDimension (0, d1);

      for (unsigned d2 = 1; d2 < d1; ++d2)
      {
         gm2.selectDimension (1, d2);

         for (unsigned d3 = 0; d3 < d2; ++d3)
         {
            gm2.selectDimension (2, d3);
            mf << tParameter (gm2);
         }
      }
   }

   return mf.getMaximum();
}

