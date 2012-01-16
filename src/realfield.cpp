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

#include <HIntLib/integerring.h>
#include <HIntLib/mymath.h>

namespace L = HIntLib;

/**
 *  RealField
 */

template<typename T>
typename L::RealField<T>::type
L::RealField<T>::element (unsigned i) const
{
   if (! i)  return type();
   --i;

   T sign = (i & 1) ? -1 : 1;

   unsigned i1, i2;
   unthread (i / 2 + 1, i1, i2);

   return sign * (T(i1) + radicalInverseFunction2 (i2));
}

template<typename T>
unsigned L::RealField<T>::index (type r) const
{
   if (is0 (r)) return 0;

#if 0
   real integerPart;
   real rest = modf (abs (r), &integerPart);
              // we need an extra real() because abs() is broken in old glibc
#else
   const real absr = abs (r);
   const real integerPart = floor (absr);
   real rest = absr - integerPart;
#endif

   unsigned i2 = 0;
   unsigned mask = 1;

   while (rest > 1e-7 && mask < (1 << 16))
   {
      if (rest >= .5 - 1e-7)
      {
         i2 |= mask;
         rest -= .5;
      }
      mask <<= 1;
      rest *= 2.0;
   }

   return (thread (unsigned (integerPart), i2) - 1) * 2 + ((r < 0) ? 2 : 1);
}

namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template RealField<X>::type RealField<X>::element(unsigned) const; \
   template unsigned RealField<X>::index(type) const; \

   HINTLIB_INSTANTIATE(real)
#undef HINTLIB_INSTANTIATE
}

