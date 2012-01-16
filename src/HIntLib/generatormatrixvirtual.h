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


#ifndef HINTLIB_GENERATOR_MATRIX_VIRTUAL_H
#define HINTLIB_GENERATOR_MATRIX_VIRTUAL_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/generatormatrix.h>

#include <HIntLib/array.h>


namespace HIntLib
{

/**
 *  ZeroMatrices
 */

class ZeroMatrices : public GeneratorMatrix
{
public:
   ZeroMatrices (int _base, int _dim)
      : GeneratorMatrix (_base, _dim, getDefaultM (_base),
                                         getDefaultPrec (_base)) {}
   ZeroMatrices (int _base, int _dim, int _m)
      : GeneratorMatrix (_base, _dim, _m,
                         std::min (_m, getDefaultPrec (_base))) {}
   ZeroMatrices (int _base, int _dim, int _m, int _prec)
      : GeneratorMatrix (_base, _dim, _m, _prec) {}

   virtual int getDigit  (int d, int r, int b) const;
   virtual u64 vGetPackedRowVector (int d, int b) const;
};


/**
 *  Identity Matrices
 */

class IdentityMatrices : public GeneratorMatrix
{
public:
   IdentityMatrices (int _base, int _dim)
      : GeneratorMatrix (_base, _dim, getDefaultM (_base),
                                         getDefaultPrec (_base)) {}
   IdentityMatrices (int _base, int _dim, int _m)
      : GeneratorMatrix (_base, _dim, _m,
                         std::min (_m, getDefaultPrec (_base))) {}
   IdentityMatrices (int _base, int _dim, int _m, int _prec)
      : GeneratorMatrix (_base, _dim, _m, _prec) {}

   virtual int getDigit  (int d, int r, int b) const;
   virtual u64 vGetPackedRowVector (int d, int b) const;
};


/**
 *  Adjust Prec
 *
 *  Adjusts the total precision.
 *
 *  If the original precision is larger, rows are discarded.
 *  If the original precision is lower, additional 0 rows are added.
 *
 *  Vectorization is set to 1.
 */

class AdjustPrec : public GeneratorMatrix
{
public:
   AdjustPrec (int p, const GeneratorMatrix& _gm)
      : GeneratorMatrix (_gm), gm (&_gm), oldPrec (_gm.getPrec())
   {
      if (p >= 0)  prec = p;
   }

   virtual int getDigit  (int d, int r, int b) const;
   virtual u64 vGetPackedRowVector (int d, int b) const;

private:
   const GeneratorMatrix* gm;
   const int oldPrec;
};

class AdjustPrecSquare : public AdjustPrec
{
public:
   AdjustPrecSquare (const GeneratorMatrix& _gm)
      : AdjustPrec (_gm.getM(), _gm) {}
};


/**
 *  Add Last Row
 */

class AddLastRow : public GeneratorMatrix
{
public:
   AddLastRow (const GeneratorMatrix& _gm)
      : GeneratorMatrix (_gm), gm (&_gm)
   {
      prec = _gm.getPrec() + 1;
   }

   virtual int getDigit  (int d, int r, int b) const;
   virtual u64 vGetPackedRowVector (int d, int b) const;

private:
   const GeneratorMatrix* gm;
};


/**
 *  Adjust m
 *
 *  Increasing m leads to net duplication
 *  Decreasing m iterates only through the first base^m points. Only for
 *     sequences a good idea.
 */

class AdjustM : public GeneratorMatrix
{
public:
   AdjustM (int _m, const GeneratorMatrix&);

   virtual int getDigit  (int d, int r, int b) const;
   virtual u64 vGetPackedRowVector (int d, int b) const;

private:
   const GeneratorMatrix* gm;
   const int oldM;
   const u64 bToTheM;
};


/**
 *  MReduction
 *
 *  Produces a submatrix with the same t-Parameter.
 *  Matrix number _identityDim_ (or default 0) must be the identity matrix.
 *
 *  The new matrix is obtained by 
 *  a) removing the first _diff_ columns from each matrix (which amounts to
 *     only selecting points with x_identityDim < 1/b^_diff_), and
 *  b) removing the first _diff_ rows from matrix number _identityDim_
 *     (which amounts to scaling each coordinate in direction _identityDim_ by
 *     a factor of b^_diff_).
 */

class MReduction : public GeneratorMatrix
{
public:
   MReduction (int _m, const GeneratorMatrix&);
   MReduction (int _m, int _identityDim, const GeneratorMatrix&);

   virtual int getDigit  (int d, int r, int b) const;
   virtual u64 vGetPackedRowVector (int d, int b) const;

private:
   const GeneratorMatrix* gm;
   const int identityDim;
   const int diff;
   const u64 bToTheDiff;
   const u64 bToTheM;
};


/**
 *  Net From Sequence
 */

class NetFromSequence : public GeneratorMatrix
{
public:
   NetFromSequence (int _m, bool _equi, const GeneratorMatrix& _gm)
      : GeneratorMatrix (_gm), gm (&_gm), equi(_equi)
   {
      m = _m;
      if (equi)  ++dim;
   }

   NetFromSequence (int _m, const GeneratorMatrix& _gm)
      : GeneratorMatrix (_gm), gm (&_gm), equi(true)
   {
      m = _m;
      ++dim;
   }

   virtual int getDigit  (int d, int r, int b) const;
   virtual u64 vGetPackedRowVector (int d, int b) const;

private:
   const GeneratorMatrix* gm;
   const bool equi;
};


/**
 *  Discard Dimensions
 */

class DiscardDimensions : public GeneratorMatrix
{
public:
   DiscardDimensions (int _dim, const GeneratorMatrix& _gm)
      : GeneratorMatrix (_gm), gm (&_gm)
   {
      if (_dim >= 0)  dim = _dim;
   }

   virtual int getDigit (int d, int r, int b) const;
   virtual u64 vGetPackedRowVector (int d, int b) const;

private:
   const GeneratorMatrix* gm;
};


/**
 *  Select Dimensions
 *
 *  Selects certain dimensions from a Generator Matrix
 */

class SelectDimensions : public GeneratorMatrix
{
public:
   SelectDimensions (int _dim, const GeneratorMatrix&);
   SelectDimensions (int dim1, int dim2, const GeneratorMatrix&);
   SelectDimensions (const int* begin, const int* end, const GeneratorMatrix&);

   void selectDimension (int dResult, int dOri)
   {
      dimensions [dResult] = dOri;
   }

   virtual int getDigit  (int d, int r, int b) const;
   virtual u64 vGetPackedRowVector (int d, int b) const;

private:
   const GeneratorMatrix* gm;
   Array<int> dimensions;
};


/**
 *  VirtualMatrixBase
 */

template<class T> class GeneratorMatrixGenRow;

class VirtualMatrixBase : public GeneratorMatrix
{
public:
   ~VirtualMatrixBase();

   virtual int getDigit  (int d, int r, int b) const;
   virtual u64 vGetPackedRowVector (int d, int b) const;

protected:
   VirtualMatrixBase () : gm(0) {}
   void setMatrix (GeneratorMatrixGenRow<unsigned char>*);

private:
   GeneratorMatrixGenRow<unsigned char>* gm;
};


/**
 *  Base Reduction
 *
 *  Given a Matrix in base p^r, a new matrix in base p yielding the same point
 *  set is generated.
 */

class BaseReduction : public VirtualMatrixBase
{
public:
   BaseReduction (const GeneratorMatrix&);
};


/**
 *  With Identity Matrix
 */

class WithIdentityMatrix : public VirtualMatrixBase
{
public:
   WithIdentityMatrix (const GeneratorMatrix& _gm)  { init (0, _gm); }
   WithIdentityMatrix (int d, const GeneratorMatrix& _gm)
      { init (d, _gm); }

private:
   void init (int d, const GeneratorMatrix&);
};

}  // namespace HIntLib

#endif

