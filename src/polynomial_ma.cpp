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

#include <HIntLib/polynomial.tcc>

#include <HIntLib/factorring.h>

namespace HIntLib
{
   HINTLIB_INSTANTIATE_POLYNOMIAL (int)
   HINTLIB_INSTANTIATE_POLYNOMIAL (unsigned char)
   HINTLIB_INSTANTIATE_POLYNOMIAL (unsigned short)

   HINTLIB_INSTANTIATE_POLYNOMIALRING_RING (FactorRing<unsigned char>)
   HINTLIB_INSTANTIATE_POLYNOMIALRING_RING (FactorRing<unsigned short>)

   HINTLIB_INSTANTIATE_POLYNOMIALRING_DOMAIN (IntegerRing<int>)

   HINTLIB_INSTANTIATE_POLYNOMIALRING_GF (FactorField<unsigned char>)
   HINTLIB_INSTANTIATE_POLYNOMIALRING_GF (FactorField<unsigned short>)
}


/**
 *  funny Sum ()
 *
 *  Return the first  n+1  terms of the sum
 *
 *      1 + 1 + base + base^2 + base^3 + ...
 */

unsigned HIntLib::Private::funnySum (int n, unsigned base)
{
   if (n < 0)  return 0;

   unsigned res = 1;
   unsigned k = 1;

   for (int i = 0; i < n; ++i)
   {
      res += k;
      k *= base;
   }
   return res;
}

