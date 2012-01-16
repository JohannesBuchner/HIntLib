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

#include <HIntLib/tparameter.h>

#include <HIntLib/myalgorithm.h>
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
 *  Checks if (at least) a certain thickness is present.
 *
 *  If _dimOpt_ is true, only partitions containing a vector from the last
 *  matrix are checked.
 */

template<typename T>
bool L::TCalc2<T>::check (int thickness, bool dimOpt)
{
   // choose these _thickness_ vectors from c
      
   initial_partition (&partition[0], &partition[dim], thickness);

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
      if (singular (thickness))  return false;
   }
   while (next_partition (&partition[0], &partition[dim]));

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
bool L::TCalc2<T>::checkTO (int thickness, bool dimOpt)
{
   // choose these _thickness_ vectors from c
      
   initial_partition (&partition[0], &partition[dim], thickness);

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
 *  Checks if (at least) a certain thickness is present. Only partitions with
 *  at most _maxRows_ per matrix considered.
 */

template<typename T>
bool L::TCalc2<T>::checkRestricted (int thickness, int maxRows)
{
   if (thickness > maxRows * dim)  return true;;

   // choose these _numVectors_ vectors from c
   
   initial_partition (&partition[0], &partition[dim], maxRows, thickness);

   do
   {
      copyToSelection ();

      if (singular (thickness))  return false;
   }
   while (next_partition (&partition[0], &partition[dim], maxRows));

   return true;
}

/**
 *  checkRestrictedRO()
 *
 *  Checks if (at least) a certain thickness is present. Only partitions with
 *  at most _maxRows_ per matrix considered.
 *
 *  The check is performed under the assumption that the thickness is present
 *  considerung at most _maxRows_-1 per matrix. 
 */

template<typename T>
bool L::TCalc2<T>::checkRestrictedRO (int thickness, int maxRows)
{
   if (thickness > maxRows * dim)  return true;;

   // choose these _numVectors_ vectors from c
   
   initial_partition (&partition[0], &partition[dim], maxRows, thickness);

   do
   {
      // Make sure we have a matrix with maxRows rows

      if (! partitionContains (maxRows))  continue;

      copyToSelection ();

      if (singular (thickness))  return false;
   }
   while (next_partition (&partition[0], &partition[dim], maxRows));

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
 */

template<typename T>
bool L::TCalc2<T>::checkRestrictedTO (int thickness, int maxRows)
{
   if (thickness > maxRows * dim)  return true;

   // choose these _numVectors_ vectors from c
   
   initial_partition (&partition[0], &partition[dim], maxRows, thickness);

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
 *  Checks if (at least) a certain thickness is present.  Only partitions with
 *  at most _maxRows_ per matrix are considered.
 *
 *  The check is performed under the assumption that a thickness of
 *  _thickness_-1 is present.
 *
 *  The check is performed under the assumption that the thickness is present
 *  considerung at most _maxRows_-1 per matrix. 
 */

template<typename T>
bool L::TCalc2<T>::checkRestrictedTORO (int thickness, int maxRows)
{
   if (thickness > maxRows * dim)  return true;

   // choose these _numVectors_ vectors from c
   
   initial_partition (&partition[0], &partition[dim], maxRows, thickness);

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
     la (makeLinearAlgebra (_gm.getBase()))
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
 *  check()
 */

bool L::TCalcGen::check (int thickness, bool dimOpt)
{
   initial_partition (&partition[0], &partition[dim], thickness);

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
      if (singular (thickness))  return false;
   }
   while (next_partition (&partition[0], &partition[dim]));

   return true;
}

/**
 *  checkRestricted()
 *
 *  Checks if (at least) a certain thickness is present. Only partitions with
 *  at most _maxRows_ per matrix considered.
 */

bool L::TCalcGen::checkRestricted (int thickness, int maxRows)
{
   if (thickness > maxRows * dim)  return true;

   initial_partition (&partition[0], &partition[dim], maxRows, thickness);

   do
   {
      copyToSelection ();
      if (singular (thickness))  return false;
   }
   while (next_partition (&partition[0], &partition[dim], maxRows));

   return true;
}

/**
 *  checkRestrictedRO()
 *
 *  Checks if (at least) a certain thickness is present. Only partitions with
 *  at most _maxRows_ per matrix considered.
 *
 *  The check is performed under the assumption that the thickness is present
 *  considerung at most _maxRows_-1 per matrix. 
 */

bool L::TCalcGen::checkRestrictedRO (int thickness, int maxRows)
{
   if (thickness > maxRows * dim)  return true;

   // Partition  thickness-maxRows  elements on  dim-1  stacks.
   // Then, add an additional stack with  maxRows  elements.
   // Adding this extra stack ensures that at least one stack with  maxRows
   //   elements is present.
   // Adding the extra stack only left of any other stacks that contain  maxRows
   //   elements ensures that no combination is created twice.

   initial_partition (
         &partition[0], &partition[dim-1], maxRows, thickness - maxRows);

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

         if (singular (thickness))  return false;
      }
   }
   while (next_partition (&partition[0], &partition[dim-1], maxRows));

   return true;
}


/*****************  Public routines  *****************************************/


bool L::confirmT (const GMGen &gm, int t, TOption opts)
{
   if (opts & LOWER_RESTRICTION_OK)  throw FIXME (__FILE__, __LINE__);

   if (t < 0 || t > int (gm.getM()))  throw FIXME (__FILE__, __LINE__);
   if (gm.getM() - t > gm.getTotalPrec())  return false;

   if (use32bit (gm))
   {
      GMCopy copy; copy.totalPrec (gm.getM() - t);
      GeneratorMatrix2Row<u32> gm2 (gm, copy);
      TCalc2<u32> calc (gm2);

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

bool L::confirmTRestricted (const GMGen &gm, int t, int maxRows, TOption opts)
{
   if (opts & LOWER_DIM_OK)  throw FIXME (__FILE__, __LINE__);

   if (t < 0 || t > int (gm.getM()))  throw FIXME (__FILE__, __LINE__);
   if (maxRows > int (gm.getTotalPrec()))  throw FIXME (__FILE__, __LINE__);

   const int thickness = gm.getM() - t;

   if (use32bit (gm))
   {
      GMCopy copy; copy.totalPrec (std::min (maxRows, thickness));
      GeneratorMatrix2Row<u32> gm2 (gm, copy);
      TCalc2<u32> calc (gm2);

      if ((opts & LOWER_RESTRICTION_OK) && thickness > 1)
      {
         if (opts & LARGER_T_OK)
            return calc.checkRestrictedTORO (thickness, maxRows);
         else
            return calc.checkRestrictedRO   (thickness, maxRows);
      }
      else
      {
         if (opts & LARGER_T_OK)
            return calc.checkRestrictedTO (thickness, maxRows);
         else
            return calc.checkRestricted   (thickness, maxRows);
      }
   }
   else
   {
      TCalcGen calc (gm);

      // there is no way to optimize for LARGER_T_OK in base > 2

      if ((opts & LOWER_RESTRICTION_OK) && thickness > 1)
      {
         return calc.checkRestrictedRO (thickness, maxRows);
      }
      else
         return calc.checkRestricted   (thickness, maxRows);
   }
}

int L::tParameter (const GMGen &gm, int lb, int ub, TOption opts)
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
      GMCopy copy; copy.totalPrec (M - lb);
      GeneratorMatrix2Row<u32> gm2 (gm, copy);
      TCalc2<u32> calc (gm2);
   
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
      const GMGen &gm, int lb, int ub, int maxRows, TOption opts)
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
      GMCopy copy; copy.totalPrec (M - lb);
      GeneratorMatrix2Row<u32> gm2 (gm, copy);
      TCalc2<u32> calc (gm2);
   
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

