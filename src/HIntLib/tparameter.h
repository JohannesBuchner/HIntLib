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

namespace HIntLib
{

enum TOption
{
   DEFAULT = 0,
   LOWER_RESTRICTION_OK = 1, // Asuume that a smaller restriction is ok
   LARGER_T_OK = 2,          // Assume that a larger t (lower thickness) is ok
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
 *  confirmT()
 */

bool confirmT (
      const GeneratorMatrixGen<unsigned char> &gm, int t, TOption opts);
bool confirmTRestricted (
   const GeneratorMatrixGen<unsigned char> &gm, int t, int maxRows,
   TOption opts);


/**
 *  tParameterRestricted()
 */

int tParameterRestricted (
   const GeneratorMatrixGen<unsigned char> &,
   int lb, int ub, int maxRows, TOption opts);

}  // namespace HIntLib

#endif

