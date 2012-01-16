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

#include <HIntLib/gcd.h>

#include <HIntLib/polynomial2.h>
#include <HIntLib/modulararithmetic.h>
#include <HIntLib/polynomial.h>

namespace L = HIntLib;

template<class A>
typename A::type HIntLib::genGcd (
   const A &a,
   typename A::type   u3, typename A::type   v3,
   typename A::type & mu, typename A::type & mv)
{
   typedef typename A::type T;
   
   T u1 = a.one();
   T u2 = a.zero();

   T v1 = a.zero();
   T v2 = a.one();

   while (! a.is0 (v3))
   {
      T q, re;
      a.div (u3, v3, q, re);

      T
      t = a.sub (u1, a.mul (q, v1)); u1 = v1; v1 = t;
      t = a.sub (u2, a.mul (q, v2)); u2 = v2; v2 = t;
      u3 = v3; v3 = re;
   }

   mu = u1; mv = u2;

   return u3;
}

template<class A>
typename A::type HIntLib::genGcd (
   const A &a,
   typename A::type u3, typename A::type v3, typename A::type& mu)
{
   typedef typename A::type T;
   
   T u1 = a.one();
   T v1 = a.zero();

   while (! a.is0 (v3))
   {
      T q, re;
      a.div (u3, v3, q, re);

      T t = a.sub (u1, a.mul (q, v1)); u1 = v1; v1 = t;
      u3 = v3; v3 = re;
   }

   mu = u1;

   return u3;
}

template<class A>
typename A::type HIntLib::genGcd (
   const A &a,
   typename A::type u, typename A::type v)
{
   while (! a.is0 (v))
   {
      typename A::type q, re;
      a.div (u, v, q, re);
      u = v; v = re;
   }

   return u;
}

namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template X::type genGcd (X const &, X::type, X::type); \
   template X::type genGcd (X const &, X::type, X::type, X::type &); \
   template X::type genGcd (X const &, X::type, X::type, X::type &, X::type &);

   HINTLIB_INSTANTIATE (HIntLib::IntegerRing<int>);
   HINTLIB_INSTANTIATE (HIntLib::PolynomialRingField<HIntLib::ModularIntegerField<unsigned char> >)
   HINTLIB_INSTANTIATE (HIntLib::PolynomialRingField<HIntLib::ModularIntegerField<unsigned short> >)
#undef HINTLIB_INSTANTIATE
}

