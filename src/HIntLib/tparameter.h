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


#ifndef HINTLIB_TPARAMETER_H
#define HINTLIB_TPARAMETER_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/generatormatrixgen.h>
#include <HIntLib/generatormatrix2row.h>
#include <HIntLib/lookupfield.h>
#include <HIntLib/linearalgebra.h>

namespace HIntLib
{


/***********  TCalc  *********************************************************/


/**
 *  TCalc
 *
 *  Base class for TCalc2 and TCalcGen
 */

class TCalc
{
protected:
   explicit TCalc (int _dim);
   ~TCalc ();

   bool partitionContains (int) const;

   const int dim;
   int *const partition;

private:
   static const unsigned MAX_PARTITION_SIZE = 200;
   static int staticPartition [MAX_PARTITION_SIZE];
};

inline
bool TCalc::partitionContains (int value) const
{
   for (int d = 0; d < dim; ++d)  if (partition[d] == value)  return true;
   return false;
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
class TCalc2 : public TCalc
{
public:
   TCalc2 (const HIntLib::GeneratorMatrix2Row<T> &_gm)
      : TCalc (_gm.getDimension()), gm(_gm) {}

   bool check (int thickness, bool dimOpt);
   bool checkTO (int thickness, bool dimOpt);
   bool checkRestricted (int thickness, int maxRows);
   bool checkRestrictedRO (int thickness, int maxRows);
   bool checkRestrictedTO (int thickness, int maxRows);
   bool checkRestrictedTORO (int thickness, int maxRows);

private:
   void init (const GeneratorMatrix &);

protected:
   void copyToSelection (int d, int num, T*& pos);
   void copyToSelection ();
   bool singular (T* l)
      { return ! isLinearlyIndependent (&selection[0], l); }
   bool singular (int thickness)
      { return ! isLinearlyIndependent (&selection[0], &selection[thickness]); }

   const GeneratorMatrix2Row<T> & gm;
   static T selection [std::numeric_limits<T>::digits];
};


/**
 *  copyToSeclection ()
 */

template<typename T>
inline
void TCalc2<T>::copyToSelection (int d, int num, T*& pos)
{
   for (int i = 0; i < num; ++i)  *(pos++) = gm(d,i);
}

template<typename T>
inline
void TCalc2<T>::copyToSelection ()
{
   T* pos = &selection[0];
   for (int d = 0; d < dim; ++d)  copyToSelection (d, partition[d], pos);
}


/***********  T Calc Gen  ****************************************************/


class TCalcGen : public TCalc
{
public:
   TCalcGen (const GeneratorMatrixGen<unsigned char> &);
   ~TCalcGen ();

   bool check (int thickness, bool dimOpt);
   bool checkRestricted (int thickness, int maxRows);
   bool checkRestrictedRO (int thickness, int maxRows);

protected:
   void copyToSelection (int d, int num, unsigned& pos);
   void copyToSelection ();
   bool singular (int thickness, int newM);
   bool singular (int thickness);

   void dumpSelection (std::ostream &) const;

   unsigned char* selection;

   const GeneratorMatrixGen<unsigned char> & gm;
   const int M;
   
private:
   LinearAlgebra* la;

   static const unsigned MAX_SELECTION_SIZE = 1000;
   static unsigned char staticSelection [MAX_SELECTION_SIZE];

   TCalcGen ();
};


/**
 * singular()
 */

inline
bool TCalcGen::singular (int thickness, int newM)
{
   return ! la->isLinearlyIndependent (&selection[0], newM, thickness);
}

inline
bool TCalcGen::singular (int thickness)
{
   return ! la->isLinearlyIndependent (&selection[0], M, thickness);
}


/**
 *  copy To Selection ()
 */

inline
void TCalcGen::copyToSelection (int d, int num, unsigned& pos)
{
   for (int i = 0; i < num; ++i)
   {
      for (int j = 0; j < M; ++j)  selection [pos * M + j] = gm (d, j, i);
      ++pos;
   }
}

inline
void TCalcGen::copyToSelection ()
{
   unsigned pos = 0;
   for (int d = 0; d < dim; ++d)  copyToSelection (d, partition[d], pos);
}


/*****************  Interface to public routines  ****************************/

enum TOption
{
   DEFAULT = 0,
   LOWER_RESTRICTION_OK = 1, // Asuume that any smaller restriction is ok
   LARGER_T_OK = 2,          // Assume that any larger t (lower thickness) is ok
   LOWER_DIM_OK = 4          // Assume that its ok without the last matrix
};


/**
 *  tParameter()
 */

int tParameter (
      const GeneratorMatrixGen<unsigned char> &, int lb, int ub,
      TOption = DEFAULT);

inline
int tParameter (
      const GeneratorMatrixGen<unsigned char> &gm, TOption opts = DEFAULT)
{
   return tParameter (gm, 0, gm.getM(), opts);
}

/**
 *  tParameterRestricted()
 */

int tParameterRestricted (
      const GeneratorMatrixGen<unsigned char> &,
      int lb, int ub, int maxRows, TOption opts = DEFAULT);

/**
 *  confirmT()
 */

bool confirmT (
      const GeneratorMatrixGen<unsigned char> &gm, int t,
      TOption opts = DEFAULT);

/**
 *  confirmTRestricted ()
 */

bool confirmTRestricted (
      const GeneratorMatrixGen<unsigned char> &gm, int t, int maxRows,
      TOption opts = DEFAULT);


}  // namespace HIntLib

#endif

