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

#ifdef __GNUG__
#pragma implementation
#endif

#include <iostream>

#include <HIntLib/integerring.h>

#include <HIntLib/prime.h>


namespace L = HIntLib;

/**
 *  IntegerRing
 */

template<typename T>
T L::IntegerRing<T>::element (unsigned i) const
{
   return odd(i)  ?  T(i/2 + 1)  :  -T(i/2);
}

template<typename T>
unsigned L::IntegerRing<T>::index (T x) const
{
   return (x > 0)  ?  x*2 - 1  :  x * -2;
}

template<typename T>
bool L::IntegerRing<T>::isPrime (T x) const
{
   return Prime::test (unsigned (abs(x)));
}

template<typename T>
bool L::IntegerRing<T>::isComposit (T x) const
{
   unsigned xx = abs (x);
   return xx > 3 && ! Prime::test (xx);
}


/**
 *  RealField
 */

template<typename T>
T L::RealField<T>::element (unsigned i) const
{
   if (i == 0)  return T(0);
   else return odd(i)  ?  T(i/2 + 1)  :  -T(i/2);
}

template<typename T>
unsigned L::RealField<T>::index (T r) const
{
   int x = int (r); return x > 0  ?  x*2 - 1 :  x * -2;
}

namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template X IntegerRing<X>::element(unsigned) const; \
   template unsigned IntegerRing<X>::index(X) const; \
   template bool IntegerRing<X>::isPrime(X) const; \
   template bool IntegerRing<X>::isComposit(X) const;

   HINTLIB_INSTANTIATE(int)
#undef HINTLIB_INSTANTIATE

#define HINTLIB_INSTANTIATE(X) \
   template X RealField<X>::element(unsigned) const; \
   template unsigned RealField<X>::index(X) const; \

   HINTLIB_INSTANTIATE(real)
#undef HINTLIB_INSTANTIATE
}

