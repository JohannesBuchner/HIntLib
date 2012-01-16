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


#ifndef HINTLIB_LINEAR_ALGEBRA_GEN_H
#define HINTLIB_LINEAR_ALGEBRA_GEN_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

namespace HIntLib
{

/*
 * Function for C-style matrices over arbitrary fields
 *
 * Definitions are given in linearalgebragen.tcc
 */

template<class A>
bool
isZeroMatrix (const A&, const typename A::type*, int numRows, int numCols);

template<class A>
bool
isIdentityMatrix (const A&, const typename A::type*, int num);

template<typename T>
void
matrixTranspose (const T*, int numRows, int numCols, T*);

template<class A>
void
matrixMul (const A&,
      const typename A::type*, const typename A::type*,
      int numRows, int numRowsCols, int numCols2, typename A::type*);

template<class A>
inline
void
matrixMul (const A& a, const typename A::type* m1, const typename A::type* m2,
      int size, typename A::type* res)
{
   matrixMul (a, m1, m2, size, size, size, res);
}

template<class A>
void
matrixVectorMul (const A& a,
      const typename A::type* m, const typename A::type* v,
      int numRows, int numCols, typename A::type* res);

template<class A>
void
vectorMatrixMul (const A& a,
      const typename A::type* v, const typename A::type* m,
      int numRows, int numCols, typename A::type* res);

template<class A>
bool
isLinearlyIndependent (const A&, typename A::type*, int numRows, int numCols);

template<class A>
int
matrixRank (const A&, typename A::type*, int numRows, int numCols);

template<class A>
int
numLinearlyIndependentVectors (
      const A&, typename A::type*, int numRows, int numCols);

template<class A>
int
nullSpace (const A&, typename A::type*, int numRows, int numCols,
           typename A::type*);

template<class A>
int
basisSupplement (const A&, typename A::type*, int numRows, int numCols);

template<class A>
int
basisSupplement (const A&, typename A::type*, int numRows, int numCols, int*);

template<class A>
bool
matrixInverse (const A&, typename A::type*, int num);

}  // namespace HIntLib

#endif

