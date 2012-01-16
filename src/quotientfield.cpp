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

#include <HIntLib/quotientfield.h>

#include <HIntLib/gcd.h>

namespace L = HIntLib;


/**
 *  index()
 */

template<class A>
unsigned
L::QuotientField<A>::index (const type& u) const
{
   if (is0 (u)) return 0;

   // examine unit multiplier

   unsigned numUnits = a.numUnits() ? a.numUnits() : 6;

   base_type num = u.num;
   base_type den = u.den;
   unsigned unit = a.unitIndex (a.makeCanonical (num));

   // examine irreducible factors

   unsigned primes = 0;

   typename A::PrimeGenerator gen (a);
   base_type irred;
   unsigned power = 1;

   while (! a.is1 (num) || ! a.is1 (den))
   {
      irred = gen.next();

      int exponent = 0;

      while (a.is0 (a.rem (num, irred)))
      {
         num = a.quot (num, irred);
         ++exponent;
      }

      while (a.is0 (a.rem (den, irred)))
      {
         den = a.quot (den, irred);
         --exponent;
      }

      if (exponent < 0)  exponent = 3 - exponent;

      primes += power * exponent;
      power *= 7;
   }

   return unit + numUnits * primes + 1;
}


/**
 *  order()
 */

template<class A>
unsigned
L::QuotientField<A>::order (const type& u) const
{
   if (is0 (u))  throwDivisionByZero();

   unsigned o1 = a.order (u.num);
   if (o1 == 0)  return 0;

   unsigned o2 = a.order (u.den);
   if (o2 == 0)  return 0;

   return lcm (o1, o2);
}


/**
 *  element()
 */

template<class A>
typename L::QuotientField<A>::type
L::QuotientField<A>::element(unsigned i) const
{
   // check for 0 (resulting in 0)
   
   if (!i)  return type();
   --i;

   // determine unit multiplier

   type x = one();

   unsigned numUnits = a.numUnits() ? a.numUnits() : 6;
   a.mulByUnit (x.num, a.unitElement (i % numUnits));
   i /= numUnits;

   // determine irreducible multipliers 

   typename A::PrimeGenerator gen (a);
   base_type irred;

   while (i)
   {
      int exp = i % 7;
      i /= 7;

      irred = gen.next();

      if (exp >= 1 && exp <= 3)  a.mulBy (x.num, a.power (irred, exp));
      else if (exp >= 4)         a.mulBy (x.den, a.power (irred, exp - 3));
   }

   return x;
}

template<class A>
typename L::QuotientField<A>::type
L::QuotientField<A>::makeElement(const base_type& num) const
{
   if (a.is0 (num))
      return type();
   else
      return type (num, a.one());
}

/**
 *  toLowestTerms()
 */

template<class A>
typename L::QuotientField<A>::type
L::QuotientField<A>::toLowestTerms (const base_type& u, const base_type& v) const
{
   if (a.is0 (u))  return type();

   const base_type& g = genGcd (a, u, v);

   if (a.isUnit (g))  return type (u, v);

   base_type den = a.quot (v, g);
   const typename A::unit_type& unit = a.unitRecip (a.makeCanonical (den));
   return type (a.mulUnit (a.quot (u, g), unit), den);
}


/**
 *  toLowestTermsAndNormalize()
 */

template<class A>
typename L::QuotientField<A>::type
L::QuotientField<A>::toLowestTermsAndNormalize (const base_type& u, const base_type& v) const
{
   if (a.is0 (u))  return type();

   const base_type& g = genGcd (a, u, v);

   base_type den = a.quot (v, g);
   const typename A::unit_type& unit = a.unitRecip (a.makeCanonical (den));
   return type (a.mulUnit (a.quot (u, g), unit), den);
}


/**
 *  add()
 */

template<class A>
typename L::QuotientField<A>::type
L::QuotientField<A>::add (const type& u, const type& v) const
{
   // Is one of the summands 0?

   if (is0 (u))  return v;
   if (is0 (v))  return u;

   // Are the denominators equal?

   if (u.den == v.den)  return toLowestTerms (a.add (u.num, v.num), u.den);

   return toLowestTerms (
      a.add (a.mul (u.num, v.den), a.mul (v.num, u.den)),
      a.mul (u.den, v.den));
}


/**
 *  sub()
 */

template<class A>
typename L::QuotientField<A>::type
L::QuotientField<A>::sub (const type& u, const type& v) const
{
   // Is one of the summands 0?

   if (is0 (u))  return neg (v);
   if (is0 (v))  return u;

   // Are the denominators equal?

   if (u.den == v.den)  return toLowestTerms (a.sub (u.num, v.num), u.den);

   return toLowestTerms (
      a.sub (a.mul (u.num, v.den), a.mul (v.num, u.den)),
      a.mul (u.den, v.den));
}


/**
 *  recip()
 */

template<class A>
typename L::QuotientField<A>::type
L::QuotientField<A>::recip (const type& u) const
{
   if (is0 (u))  throwDivisionByZero();

   base_type den = u.num; 
   typename A::unit_type unit = a.unitRecip (a.makeCanonical (den));
   return type (a.mulUnit (u.den, unit), den);
}


/**
 *  div()
 */

template<class A>
typename L::QuotientField<A>::type
L::QuotientField<A>::div (const type& u, const type& v) const
{
   if (is0 (u))  return type();
   if (is0 (v))  throwDivisionByZero();

   return toLowestTermsAndNormalize (a.mul(u.num, v.den), a.mul(u.den, v.num));
}


/**
 *  print Short()
 */

template<typename A>
void
L::QuotientField<A>::printShort (std::ostream& o, const type& u) const
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   a.printShort (ss, u.num);

   if (! is0 (u) && ! a.is1 (u.den))
   {
      ss << "/";
      a.printShort (ss, u.den);
   }
   
   o << ss.str().c_str();
}

template<typename A>
void
L::QuotientField<A>::print (std::ostream& o, const type& u) const
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   printShort (ss, u);
   ss << ' ';
   printSuffix (ss);

   o << ss.str().c_str();
}

/**
 *  operator<<
 */

template<typename A>
std::ostream&
L::operator<< (std::ostream &o, const QuotientField<A> &field)
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   const A alg = field.getBaseAlgebra();

   ss << alg << '/' << alg;

   return o << ss.str().c_str();
}

#include <HIntLib/integerring.h>
#include <HIntLib/polynomial.h>
#include <HIntLib/lookupfield.h>

namespace HIntLib
{

#define HINTLIB_INSTANTIATE(X) \
template QuotientField<X >::type \
         QuotientField<X >::toLowestTerms \
            (const base_type&, const base_type&) const; \
template QuotientField<X >::type \
         QuotientField<X >::toLowestTermsAndNormalize \
            (const base_type&, const base_type&) const; \
template QuotientField<X >::type \
         QuotientField<X >::element(unsigned) const; \
template QuotientField<X >::type \
         QuotientField<X >::makeElement(const base_type&) const; \
template unsigned QuotientField<X >::index (const type&) const; \
template QuotientField<X >::type \
         QuotientField<X >::add (const type&, const type&) const; \
template QuotientField<X >::type \
         QuotientField<X >::sub (const type&, const type&) const; \
template QuotientField<X >::type \
         QuotientField<X >::div (const type&, const type&) const; \
template QuotientField<X >::type QuotientField<X >::recip (const type&) const; \
template unsigned QuotientField<X >::order (const type&) const; \
template void QuotientField<X >::printShort (std::ostream&, const type&) const; \
template void QuotientField<X >::print      (std::ostream&, const type&) const; \
template std::ostream& operator<< (std::ostream &, const QuotientField<X > &);

   HINTLIB_INSTANTIATE (IntegerRing<int>)
   HINTLIB_INSTANTIATE (PolynomialRing<RealField<real> >)
   HINTLIB_INSTANTIATE (PolynomialRing<LookupField<unsigned char> >)
#undef HINTLIB_INSTANTIATE
}


