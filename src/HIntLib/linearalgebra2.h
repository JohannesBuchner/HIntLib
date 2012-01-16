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


#ifndef HINTLIB_LINEAR_ALGEBRA_2_H
#define HINTLIB_LINEAR_ALGEBRA_2_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <iterator>

namespace HIntLib
{
   /*
    * Functions for packed matrices over GF(2)
    */

   template<typename Bi>
      void packMatrix (const unsigned char*, int numCols, Bi,Bi);
   template<typename Bi>
      void unpackMatrix (Bi,Bi, int numCols, unsigned char*);

   template<typename Bi> bool isZeroMatrix          (Bi,Bi);
   template<typename Bi> bool isIdentityMatrix      (Bi,Bi);

   template<typename Bi>
      typename std::iterator_traits<Bi>::value_type
      matrixVectorMul (Bi, Bi, typename std::iterator_traits<Bi>::value_type);
   template<typename Bi>
      typename std::iterator_traits<Bi>::value_type
      vectorMatrixMul (typename std::iterator_traits<Bi>::value_type, Bi, Bi);

   template<typename Bi> bool isLinearlyIndependent (Bi,Bi);
   template<typename Bi> int matrixRank        (Bi,Bi);
   template<typename Bi> int nullSpace  (Bi, Bi, int, Bi);
   template<typename Bi> int nullSpaceT (Bi, Bi, Bi);

}  // namespace HIntLib

#endif

