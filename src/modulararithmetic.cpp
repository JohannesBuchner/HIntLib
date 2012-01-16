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

#ifdef HINTLIB_HAVE_OSTREM
  #include <ostream>
#else
  #include <iostream>
#endif

#ifdef HINTLIB_HAVE_SSTREAM
  #include <sstream>
#else
  #include <HIntLib/fallback_sstream.h>
#endif

#include <HIntLib/modulararithmetic.h>

#include <HIntLib/prime.h>
#include <HIntLib/polynomial.h>
#include <HIntLib/gcd.h>

namespace L = HIntLib;


/**
 *  Modular Integer Ring Base
 */

std::ostream&
L::Priv::operator<< (std::ostream &o, const ModularIntegerRingBase &a)
{
   return o << "Z_" << a.modulus();
}

void L::Priv::ModularIntegerRingBase::checkField() const
{
   if (! Prime::test (m))  throw InvalidModularFieldSize (m);
}

unsigned L::Priv::ModularIntegerRingBase::fieldRecip (unsigned x) const
{
   if (x == 0)  throw DivisionByZero();
   if (x == 1)  return 1;

   if (m < 8)
   {
      for (unsigned i = 2; ; ++i)
      {
         if ((i * x) % m == 1)  return i;
      }
   }
   else
   {
      int a;
      gcd (int (x), int (m), a);
      if (a < 0)  a+= m;
      return a;
   }
}


/**
 *  Modular Arith Ring
 */

template<typename A>
std::ostream&
L::operator<< (std::ostream &o, const ModularArith<A> &a)
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   ss << a.arithmetic() << " mod ";
   a.arithmetic().printShort (ss, a.modulus());
   return o << ss.str().c_str();
}

template<typename A>
void
L::ModularArith<A>::print (std::ostream &o, const type& x) const
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   a.printShort (ss, x);
   ss << " (";
   a.printShort (ss, m);
   ss << ')';
   a.printSuffix (ss);

   o << ss.str().c_str();
}

template<typename A>
void
L::ModularArith<A>::printShort (std::ostream &o, const type& x) const
{
   a.printShort (o, x);
}

template<typename A>
void
L::ModularArith<A>::printSuffix (std::ostream &o) const
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   ss << '(';
   a.printShort (ss, m);
   ss << ')';
   a.printSuffix (ss);

   o << ss.str().c_str();
}


/**
 *  Modular Arith Field
 */

template<typename A>
typename L::ModularArithField<A>::type
L::ModularArithField<A>::recip (const type& x) const
{
   type b;
   type res = genGcd (a, x, m, b);
   if (! a.isUnit (res))
   {
      throw InternalError (__FILE__, __LINE__);
   }
   return mul (a.unitRecip(res), b);
}


template<typename A>
unsigned
L::ModularArithField<A>::characteristic () const
{
   if (size() == 0)  return 0;

   unsigned prime;
   unsigned power;

   Prime::factorPrimePower (size(), prime, power);

   return prime;
}

namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template ModularArithField<X>::type \
            ModularArithField<X>::recip (const type&) const; \
   template unsigned ModularArithField<X>::characteristic() const; \
   template std::ostream& operator<< (std::ostream&, const ModularArith<X>&); \
   template void ModularArith<X>::print (std::ostream&, const type&) const; \
   template void ModularArith<X>::printShort(std::ostream&, const type&) const;\
   template void ModularArith<X>::printSuffix (std::ostream&) const;

   HINTLIB_INSTANTIATE (PolynomialRing<ModularArithField<unsigned char> >);
   HINTLIB_INSTANTIATE (PolynomialRing<ModularArithField<unsigned short> >);
#undef HINTLIB_INSTANTIATE

}

