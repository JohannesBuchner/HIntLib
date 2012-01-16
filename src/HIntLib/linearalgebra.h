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


#ifndef HINTLIB_LINEAR_ALGEBRA_H
#define HINTLIB_LINEAR_ALGEBRA_H 1

#ifdef __GNUG__
#pragma interface
#endif

namespace HIntLib
{

template<class Bi> bool isLinearlyIndependent (Bi first, Bi last);

template<class A>
bool isLinearlyIndependent (
      const A&, typename A::type*, unsigned numCols, unsigned numRows);

template<class A>
unsigned matrixRank (
      const A&, typename A::type*, unsigned numCols, unsigned numRows);

template<class A>
void matrixMul (const A&, typename A::type*,
      const typename A::type*, const typename A::type*, unsigned num);

template<class A>
bool matrixInverse (const A&, typename A::type*, unsigned num);

class LinearAlgebra
{
public:
   virtual ~LinearAlgebra() {}
   virtual bool isLinearlyIndependent (
      unsigned char* m, unsigned numCols, unsigned numRows) = 0;
   virtual unsigned matrixRank (
      unsigned char*, unsigned numCols, unsigned numRows) = 0;
   virtual bool matrixInverse  (unsigned char*, unsigned num) = 0;
   virtual void matrixMul (
      unsigned char*, const unsigned char*, const unsigned char*,
      unsigned num) = 0;

protected:
   LinearAlgebra () {}
};

LinearAlgebra* makeLinearAlgebra (unsigned base);



}  // namespace HIntLib

#endif

