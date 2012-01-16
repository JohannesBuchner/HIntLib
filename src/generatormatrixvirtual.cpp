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

#include <HIntLib/generatormatrixvirtual.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/hlmath.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;


/**
 * Zero Matrices
 */

int
L::ZeroMatrices::getDigit (int /* d */, int /* r */, int /* b */) const
{
   return 0;
}

L::u64
L::ZeroMatrices::vGetPackedRowVector (int /* d */, int /* b */) const
{
   return 0;
}


/**
 * Identity Matrices
 */

int
L::IdentityMatrices::getDigit (int /* d */, int r, int b) const
{
   return int (r == b);
}

L::u64
L::IdentityMatrices::vGetPackedRowVector (int /* d */, int b) const
{
   return powInt (u64 (base), b);
}


/**
 * Adjust Total Precision
 */

int
L::AdjustPrec::getDigit  (int d, int r, int b) const
{
   return (b < oldPrec) ? gm->getDigit (d, r, b) : 0;
}

L::u64
L::AdjustPrec::vGetPackedRowVector (int d, int b) const
{
   return (b < oldPrec) ? gm->vGetPackedRowVector (d,b) : 0;
}


/**
 *  Add Last Row
 */

int
L::AddLastRow::getDigit  (int d, int r, int b) const
{
   if (b == prec - 1)
   {
      b = 0;
      if (++d == dim)  d = 0;
   }
      
   return gm->getDigit (d, r, b);
}

L::u64
L::AddLastRow::vGetPackedRowVector (int d, int b) const
{
   if (b == prec - 1)
   {
      b = 0;
      if (++d == dim)  d = 0;
   }

   return gm->vGetPackedRowVector (d,b);
}


/**
 *  Adjust M
 */

L::AdjustM::AdjustM (int _m, const GeneratorMatrix& _gm)
   : GeneratorMatrix (_gm),
     gm (&_gm),
     oldM (_gm.getM()),
     bToTheM (powInt (_gm.getBase(), _m))
{
   if (_m >= 0)  m = _m;
}

int
L::AdjustM::getDigit  (int d, int r, int b) const
{
   return (r < oldM) ? gm->getDigit (d, r, b) : 0;
}

L::u64
L::AdjustM::vGetPackedRowVector (int d, int b) const
{
   return gm->vGetPackedRowVector (d, b) % bToTheM;
}


/**
 *  MReduction
 */

L::MReduction::MReduction (int _m, const GeneratorMatrix& _gm)
   : GeneratorMatrix (_gm),
     gm (&_gm),
     identityDim (0),
     diff (_gm.getM() - _m),
     bToTheDiff (powInt (_gm.getBase(), diff)),
     bToTheM    (powInt (_gm.getBase(), _m))
{
   if (_gm.getM() < _m)  throw FIXME(__FILE__, __LINE__);
   m = _m;
}

L::MReduction::MReduction (int _m, int _identityDim, const GeneratorMatrix& _gm)
   : GeneratorMatrix (_gm),
     gm (&_gm),
     identityDim (_identityDim),
     diff (_gm.getM() - _m),
     bToTheDiff (powInt (_gm.getBase(), diff)),
     bToTheM    (powInt (_gm.getBase(), _m))
     
{
   if (_gm.getM() < _m)  throw FIXME(__FILE__, __LINE__);
   m = _m;
}

int
L::MReduction::getDigit (int d, int r, int b) const
{
   if (d == identityDim)
   {
      const int bb = b + diff;
      return (bb >= prec) ? 0 : gm->getDigit (d, r + diff, bb);
   }
   else
   {
      return gm->getDigit (d, r + diff, b);
   }
}

L::u64
L::MReduction::vGetPackedRowVector (int d, int b) const
{
   if (d == identityDim)
   {
      const int bb = b + diff;
      return (bb >= prec) ? 0 :
                 (gm->vGetPackedRowVector (d, bb) / bToTheDiff);
   }
   else
   {
      return gm->vGetPackedRowVector (d, b) % bToTheM;
   }
}


/**
 *  Net From Sequence
 */

int
L::NetFromSequence::getDigit (int d, int r, int b) const
{
   if (equi && d == 0)  return int (r + b + 1 == m);
   else return gm->getDigit (d - equi, r, b);
}

L::u64
L::NetFromSequence::vGetPackedRowVector (int d, int b) const
{
   if (equi && d == 0)
   {
      return (b <= m)  ?  powInt (u64(base), m - 1 - b)  :  0;
   }
   else return gm->vGetPackedRowVector (d - equi, b);
}


/**
 *  Discard Dimensions
 */

int
L::DiscardDimensions::getDigit  (int d, int r, int b) const
{
   return gm->getDigit (d, r, b);
}

L::u64
L::DiscardDimensions::vGetPackedRowVector (int d, int b) const
{
   return gm->vGetPackedRowVector (d, b);
}


/**
 *  Select Dimensions
 */

L::SelectDimensions::SelectDimensions (
      int _dim, const GeneratorMatrix& _gm)
   : GeneratorMatrix (_gm), gm (&_gm), dimensions (_dim)
{
   dim = _dim;
   for (int d = 0; d < _dim; ++d)  dimensions [d] = d;
}

L::SelectDimensions::SelectDimensions (
      int dim1, int dim2, const GeneratorMatrix& _gm)
   : GeneratorMatrix (_gm), gm (&_gm), dimensions (dim2 - dim1)
{
   dim = dim2 - dim1;
   for (int d = 0; d < dim; ++d)  dimensions [d] = dim1++;
}

L::SelectDimensions::SelectDimensions (
   const int* begin, const int* end, const GeneratorMatrix& _gm)
   : GeneratorMatrix (_gm), gm (&_gm), dimensions (end - begin)
{
   dim = end - begin;
   std::copy (begin, end, dimensions.begin());
}

int
L::SelectDimensions::getDigit  (int d, int r, int b) const
{
   return gm->getDigit (dimensions [d], r, b);
}

L::u64
L::SelectDimensions::vGetPackedRowVector (int d, int b) const
{
   return gm->vGetPackedRowVector (dimensions[d], b);
}


