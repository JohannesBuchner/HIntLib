/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schürer <rudolf.schuerer@sbg.ac.at>
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

#include <HIntLib/polynomial_ring.tcc>

#include <HIntLib/integerring.h>
#include <HIntLib/lookupfield.h>
#include <HIntLib/modulararithmetic.h>

namespace HIntLib
{
   // Polynomial<>

   HINTLIB_INSTANTIATE_POLYNOMIAL (int)

   HINTLIB_INSTANTIATE_POLYNOMIAL (Polynomial<int>)
   HINTLIB_INSTANTIATE_POLYNOMIAL (Polynomial<unsigned char>)

   // PolynomialRing<>

   HINTLIB_INSTANTIATE_POLYNOMIALRING_UFD (PolynomialRing<IntegerRing<int> >)
   HINTLIB_INSTANTIATE_POLYNOMIALRING_UFD
      (PolynomialRing<LookupField<unsigned char> > )

   HINTLIB_INSTANTIATE_POLYNOMIALRING_UFD (IntegerRing<int>)

   HINTLIB_INSTANTIATE_POLYNOMIALRING_RING
      (ModularArithmeticRing<unsigned char>)
   HINTLIB_INSTANTIATE_POLYNOMIALRING_RING
      (ModularArithmeticRing<unsigned short>)
}

