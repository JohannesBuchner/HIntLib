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

#ifndef HINTLIB_GCD_H
#define HINTLIB_GCD_H 1

#include <HIntLib/integerring.h>

namespace HIntLib
{

/**
 *  genGcd()
 *
 *  Calculates the greatest common divisor
 *
 *  A has to be an Euclidean Ring
 */

template<class A>
typename A::type genGcd (
   const A&,
   typename A::type    u, typename A::type    v,
   typename A::type & mu, typename A::type & mv);
 
template<class A>
typename A::type genGcd (
   const A&,
   typename A::type    u, typename A::type    v,
   typename A::type & mu);

template<class A>
typename A::type genGcd (
   const A&,
   typename A::type    u, typename A::type    v);


/**
 *  gen Lcm ()
 */

template<class A>
typename A::type genLcm (const A&a, typename A::type u, typename A::type v)
{
   return a.quot (a.mul (u, v), genGcd (a, u, v));
}


/**
 *  gcd()
 *
 *  Calculates the greatest common divisor
 *
 *  Works for
 *     - integers
 *     - GenPolynomial2
 */
 
template<typename T> T gcd(T u, T v)  HINTLIB_GNU_CONST;
template<typename T> T gcd(T u, T v)
{
   while (v)
   {
      T r = u % v;
      u = v;
      v = r;
   }

   return u;
}

template<typename T> inline T gcd (T a, T b, T &ma)
{
   IntegerRing<T> r;
   return genGcd (r, a, b, ma);
}

template<typename T> inline T gcd (T a, T b, T &ma, T &mb)
{
   IntegerRing<T> r;
   return genGcd (r, a, b, ma, mb);
}


/**
 *  lcm ()
 */

template<typename T> T lcm(T u, T v)  HINTLIB_GNU_CONST;
template<typename T> inline T lcm (T u, T v)
{
   return u * v / gcd (u,v);
}


}  // namespace HIntLib

#endif
