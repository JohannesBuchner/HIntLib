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

#include <HIntLib/factorring.h>

#include <HIntLib/gcd.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;


/**
 *  Modular Integer Ring Base
 */

std::ostream&
L::Priv::operator<< (std::ostream &o, const ModularIntegerRingBase &a)
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   ss << "Z/(" << a.modulus() << ")";

   return o << ss.str().c_str();
}

void L::Priv::ModularIntegerRingBase::checkField() const
{
   if (! Prime::test (m))  throw InvalidModularFieldSize (m);
}

unsigned L::Priv::ModularIntegerRingBase::fieldRecip (unsigned x) const
{
   throwDivisionByZero (x);

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

unsigned L::Priv::ModularIntegerRingBase::ringAdditiveOrder (unsigned a) const
{
   return m / gcd (a, m);
}

unsigned L::Priv::ModularIntegerRingBase::fieldOrder (unsigned a) const
{
   throwDivisionByZero (a);

   unsigned x = a;
   unsigned n = 1;

   while (x != 1)
   {
      x = (x * a) % m;
      ++n;
   }
   return n;
}



/**
 *  Factor Ring
 */

template<typename A>
std::ostream&
L::operator<< (std::ostream &o, const FactorRing<A> &a)
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   ss << a.arithmetic() << "/(";
   a.arithmetic().printShort (ss, a.modulus());
   ss << ")";

   return o << ss.str().c_str();
}

template<typename A>
void
L::FactorRing<A>::print (std::ostream &o, const type& x) const
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
L::FactorRing<A>::printShort (std::ostream &o, const type& x) const
{
   a.printShort (o, x);
}

template<typename A>
void
L::FactorRing<A>::printSuffix (std::ostream &o) const
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

template<typename A>
unsigned
L::FactorRing<A>::additiveOrder (const type& u) const
{
   type x (u);
   unsigned n = 1;

   while (! is0 (x))
   {
      addTo (x, u);
      ++n;
   }
   return n;
}


/**
 *  Modular Arith Field
 */

template<class A>
void
L::FactorField<A>::throwIfZero (const type& x) const
{
   if (is0 (x)) throw DivisionByZero();
}

template<class A>
typename L::FactorField<A>::type
L::FactorField<A>::recip (const type& x) const
{
   throwIfZero (x);

   type b;
   type res = genGcd (a, x, m, b);
   if (! a.isUnit (res))  throw InternalError (__FILE__, __LINE__);
   return a.mulByUnit (b, a.unitRecip(a.toUnit (res)));
}


template<class A>
unsigned L::FactorField<A>::order (const type& u) const
{
   throwIfZero (u);

   type x (u);
   unsigned n = 1;

   while (! is1 (x))
   {
      mulBy (x, u);
      ++n;
   }
   return n;
}

template<class A>
bool L::FactorField<A>::isPrimitiveElement (const type& u) const
{
   throwIfZero (u);

   const unsigned q = size() - 1;
   if (Prime::test(q))  return ! is1(u);

   PrimeDivisors pd (q);

   while (unsigned prime = pd.next())
   {
      if (is1 (power(u, q / prime)))  return false;
   }
   return true;
}


#include <HIntLib/polynomial.h>

namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template FactorField<X >::type \
            FactorField<X >::recip (const type&) const; \
   template unsigned FactorField<X >::order(const type&) const; \
   template bool FactorField<X >::isPrimitiveElement(const type&) const; \
   template void FactorField<X >::throwIfZero(const type&) const; \
   template std::ostream& operator<< (std::ostream&, const FactorRing<X >&); \
   template void FactorRing<X >::print (std::ostream&, const type&) const; \
   template void FactorRing<X >::printShort(std::ostream&,const type&) const;\
   template void FactorRing<X >::printSuffix (std::ostream&) const; \
   template unsigned FactorRing<X >::additiveOrder (const type&) const;

   HINTLIB_INSTANTIATE (PolynomialRing<FactorField<unsigned char> >)
   HINTLIB_INSTANTIATE (PolynomialRing<FactorField<unsigned short> >)
#undef HINTLIB_INSTANTIATE

}

