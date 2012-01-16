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

#include <HIntLib/defaults.h>

#ifdef HAVE_OSTREM
  #include <ostream>
#else
  #include <iostream>
#endif

#include <HIntLib/modulararithmetic.h>

#include <HIntLib/prime.h>
#include <HIntLib/polynomial.h>
#include <HIntLib/gcd.h>

namespace L = HIntLib;


/**
 *  Modular Integer Ring
 */

std::ostream& L::operator<< (std::ostream &o, const ModularIntegerRingBase &a)
{
   return o << "Z_" << a.modulus();
}


/**
 *  Modular Arithmetic Ring
 */

template<class A>
typename L::ModularArithmeticRing<A>::T
L::ModularArithmeticRing<A>::rem (const T& x) const
{
   T q, r;
   a.div (x, m, q, r);
   return r;
}

template<class A>
std::ostream&
L::operator<< (std::ostream &o, const ModularArithmeticRing<A> &a)
{
   return o << a.arithmetic() << " mod " << a.modulus();
}


/**
 *  Modular Integer Field
 */

template<class T>
L::ModularIntegerField<T>::ModularIntegerField (const T _m)
   : ModularIntegerRing<T> (_m)
{
   if (! Prime::test (unsigned (_m)))  throw InvalidModularFieldSize (_m);
}

template<class T>
T L::ModularIntegerField<T>::recip (const T& x) const
{
   int a;
   gcd (int (x), int (m), a);
   if (a < 0)  a+= m;
   return a;
}


/**
 *  Modular Arithmetic Field
 */

template<class A>
typename L::ModularArithmeticField<A>::T
L::ModularArithmeticField<A>::recip (const T& x) const
{
   T b, bb;
   T res = genGcd (a, x, m, b, bb);
   if (! a.isUnit (res))
   {
      throw InternalError (__FILE__, __LINE__);
   }
   return mul (a.unitRecip(res), b);
}

namespace HIntLib
{
   // Instantiate ModularIntegerField<>

#define HINTLIB_INSTANTIATE(X) \
   template ModularIntegerField<X>::ModularIntegerField (const X); \
   template X ModularIntegerField<X>::recip (const X&) const; 

   HINTLIB_INSTANTIATE (unsigned char);
   HINTLIB_INSTANTIATE (unsigned short);
#undef HINTLIB_INSTANTIATE

   // Instantiate ModularArithmeticRing<>
   //             ModularArithmeticField<>

#define HINTLIB_INSTANTIATE(X) \
   template ModularArithmeticRing<X>::type \
            ModularArithmeticRing<X>::rem (const T&) const; \
   template ModularArithmeticField<X>::type \
            ModularArithmeticField<X>::recip (const T&) const; \
   template std::ostream& \
      operator<< (std::ostream &, const ModularArithmeticRing<X> &);

   HINTLIB_INSTANTIATE (PolynomialRingField<ModularIntegerField<unsigned char> >);
   HINTLIB_INSTANTIATE (PolynomialRingField<ModularIntegerField<unsigned short> >);
#undef HINTLIB_INSTANTIATE

}

