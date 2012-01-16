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

#include <HIntLib/quotientfield.h>

#include <HIntLib/gcd.h>
#include <HIntLib/output.h>
#include <HIntLib/bitop.h>
#include <HIntLib/exception.h>


/**
 *  element()
 */

template<class A>
typename HIntLib::Private::QFB<A>::type
HIntLib::Private::QFB<A>::element(unsigned ind) const
{
   // check for 0 (resulting in 0)
   
   if (!ind)  return type();
   --ind;

   // determine unit multiplier

   unsigned unitIndex;

   if (const unsigned numUnits = a.numUnits())
   {
      unitIndex = ind % numUnits;
      ind /= numUnits;
   }
   else
   {
      unthread (ind, ind, unitIndex, 2, 1);
   }

   type x (a.fromUnit (a.unitElement (unitIndex)), a.one());

   // determine irreducible multipliers 

   typename A::PrimeGenerator gen (a);
   unsigned indices [10];
   unsigned num = unthreadinf (ind, indices);

   for (unsigned i = 0; i < num; ++i)
   {
      base_type irred = gen.next();

      // produce proper power in the order 1, x, 1/x, x^2, 1/x^2,...

      if (const unsigned index = indices[i])
      {
         unsigned exponent = (index + 1) / 2;

         if (index & 1)  a.mulBy (x.num, a.power (irred, exponent));
         else            a.mulBy (x.den, a.power (irred, exponent));
      }
   }

   return x;
}

/**
 *  index()
 */

template<class A>
unsigned
HIntLib::Private::QFB<A>::index (const type& u) const
{
   if (is0 (u)) return 0;

   // Split into numerator, denominator und unit

   base_type num = u.num;
   base_type den = u.den;
   unsigned unit = a.unitIndex (a.makeCanonical (num));

   // examine irreducible factors

   const unsigned MAX_NUM_PRIMES = 10;
   unsigned indices [MAX_NUM_PRIMES];  // sqrt(digits) is sufficient
   unsigned n = 0;
   typename A::PrimeGenerator gen (a);

   while (! (a.is1 (num) && a.is1 (den)) && n < MAX_NUM_PRIMES)
   {
      base_type irred = gen.next();

      int exponent = 0;
      int sign = 0;

      base_type divisor;

      while (a.isDivisor (num, irred, divisor))
      {
         sign = 1;
         destructiveAssign (num, divisor);
         ++exponent;
      }

      if (sign == 0)
      {
         while (a.isDivisor (den, irred, divisor))
         {
            destructiveAssign (den, divisor);
            ++exponent;
         }
      }

      indices[n++] = 2 * exponent - sign;
   }

   const unsigned primes = threadinf (indices, n);

   // Combine unit-index and prime-index

   const unsigned numUnits = a.numUnits();

   return numUnits ? (unit + numUnits * primes + 1)
                   : (thread (primes, unit, 2, 1) + 1);
}


/**
 *  makeElement()
 */

template<class A>
typename HIntLib::Private::QFB<A>::type
HIntLib::Private::QFB<A>::makeElement(const base_type& num) const
{
   if (a.is0 (num))
      return type();
   else
      return type (num, a.one());
}


/**
 *  order()
 */

template<class A>
unsigned
HIntLib::Private::QFB<A>::order (const type& u) const
{
   if (is0 (u))  throwDivisionByZero();
   return (a.is1(u.den) && a.isUnit(u.num)) ? a.order (u.num) : 0;
}


/**
 *  toLowestTerms()
 */

template<class A>
typename HIntLib::Private::QFB<A>::type
HIntLib::Private::QFB<A>::
toLowestTerms (const base_type& u, const base_type& v) const
{
   if (a.is0 (u))  return type();

   const base_type& g = genGcd (a, v, u);

   if (a.isUnit (g))  return type (u, v);

   base_type den = a.div (v, g);
   const typename A::unit_type& unit = a.unitRecip (a.makeCanonical (den));
   return type (a.mulUnit (a.div (u, g), unit), den);
}


/**
 *  toLowestTermsAndNormalize()
 */

template<class A>
typename HIntLib::Private::QFB<A>::type
HIntLib::Private::QFB<A>::
toLowestTermsAndNormalize (const base_type& u, const base_type& v) const
{
   if (a.is0 (u))  return type();

   const base_type& g = genGcd (a, v, u);

   base_type den = a.div (v, g);
   const typename A::unit_type& unit = a.unitRecip (a.makeCanonical (den));
   return type (a.mulUnit (a.div (u, g), unit), den);
}


/**
 *  add()
 *
 *  See TACP, vol 2, 4.5.1
 */

template<class A>
typename HIntLib::Private::QFB<A>::type
HIntLib::Private::QFB<A>::add (const type& u, const type& v) const
{
   // Is one of the summands 0?

   if (is0 (u))  return v;
   if (is0 (v))  return u;

   // Are the denominators equal?

   if (u.den == v.den)  return toLowestTerms (a.add (u.num, v.num), u.den);

   // Calculate gcd of the denomiators

   const base_type& d1 = genGcd (a, u.den, v.den);

   if (a.is1 (d1))
   {
      return type (
         a.add (a.mul (u.num, v.den), a.mul (v.num, u.den)),
         a.mul (u.den, v.den));
   }

   const base_type& t = a.add (a.mul (u.num, a.div (v.den, d1)),
                               a.mul (v.num, a.div (u.den, d1)));
   const base_type& d2 = genGcd (a, d1, t);

   base_type den = a.mul (a.div (u.den, d1), a.div (v.den, d2));

   const typename A::unit_type& unit = a.unitRecip (a.makeCanonical (den));
   return type (a.mulUnit (a.div (t, d2), unit), den);
}


/**
 *  sub()
 *
 *  See TACP, vol 2, 4.5.1
 */

template<class A>
typename HIntLib::Private::QFB<A>::type
HIntLib::Private::QFB<A>::sub (const type& u, const type& v) const
{
   // Is one of the summands 0?

   if (is0 (u))  return neg (v);
   if (is0 (v))  return u;

   // Are the denominators equal?

   if (u.den == v.den)  return toLowestTerms (a.sub (u.num, v.num), u.den);

   // Calculate gcd of the denomiators

   const base_type& d1 = genGcd (a, u.den, v.den);

   if (a.is1 (d1))
   {
      return type (
         a.sub (a.mul (u.num, v.den), a.mul (v.num, u.den)),
         a.mul (u.den, v.den));
   }

   const base_type& t = a.sub (a.mul (u.num, a.div (v.den, d1)),
                               a.mul (v.num, a.div (u.den, d1)));
   const base_type& d2 = genGcd (a, d1, t);

   base_type den = a.mul (a.div (u.den, d1), a.div (v.den, d2));

   const typename A::unit_type& unit = a.unitRecip (a.makeCanonical (den));
   return type (a.mulUnit (a.div (t, d2), unit), den);
}


/**
 *  mul()
 *
 *  See TACP, vol 2, 4.5.1
 */

template<class A>
typename HIntLib::Private::QFB<A>::type
HIntLib::Private::QFB<A>::mul (const type& u, const type& v) const
{
   if (is0 (u) || is0 (v))  return type();

   const base_type& d1 = genGcd (a, u.num, v.den);
   const base_type& d2 = genGcd (a, v.num, u.den);

   const bool d1Unit = a.isUnit (d1);
   const bool d2Unit = a.isUnit (d2);

   if (d1Unit && d2Unit)
   {
      return type (a.mul (u.num, v.num), a.mul (u.den, v.den));
   }
   else if (d1Unit)
   {
      base_type den = a.mul (a.div (u.den, d2), v.den);

      const typename A::unit_type& unit = a.unitRecip (a.makeCanonical (den));

      return type (a.mulUnit (a.mul (u.num, a.div (v.num, d2)), unit), den);
   }
   else if (d2Unit)
   {
      base_type den = a.mul (u.den, a.div (v.den, d1));

      const typename A::unit_type& unit = a.unitRecip (a.makeCanonical (den));

      return type (a.mulUnit (a.mul (a.div (u.num, d1), v.num), unit), den);
   }
   else
   {
      base_type den = a.mul (a.div (u.den, d2), a.div (v.den, d1));

      const typename A::unit_type& unit = a.unitRecip (a.makeCanonical (den));

      return type (a.mulUnit (
               a.mul (a.div (u.num, d1), a.div (v.num, d2)),
            unit), den);
   }
}


/**
 *  div()
 *
 *  See TACP, vol 2, 4.5.1
 */

template<class A>
typename HIntLib::Private::QFB<A>::type
HIntLib::Private::QFB<A>::div (const type& u, const type& v) const
{
   if (is0 (u))  return type();
   if (is0 (v))  throwDivisionByZero();

   const base_type& d1 = genGcd (a, u.num, v.num);
   const base_type& d2 = genGcd (a, v.den, u.den);  // d2 is always canonical

   const bool d1Unit = a.isUnit (d1);
   const bool d2One  = a.isUnit (d2);

   if (d1Unit && d2One)
   {
      base_type den = a.mul (u.den, v.num);

      const typename A::unit_type& unit = a.unitRecip (a.makeCanonical (den));

      return type (a.mulUnit (a.mul (u.num, v.den), unit), den);
   }
   else if (d1Unit)
   {
      base_type den = a.mul (a.div (u.den, d2), v.num);

      const typename A::unit_type& unit = a.unitRecip (a.makeCanonical (den));

      return type (a.mulUnit (a.mul (u.num, a.div (v.den, d2)), unit), den);
   }
   else if (d2One)
   {
      base_type den = a.mul (u.den, a.div (v.num, d1));

      const typename A::unit_type& unit = a.unitRecip (a.makeCanonical (den));

      return type (a.mulUnit (a.mul (a.div (u.num, d1), v.den), unit), den);
   }
   else
   {
      base_type den = a.mul (a.div (u.den, d2), a.div (v.num, d1));

      const typename A::unit_type& unit = a.unitRecip (a.makeCanonical (den));

      return type (a.mulUnit (
               a.mul (a.div (u.num, d1), a.div (v.den, d2)),
            unit), den);
   }
}


/**
 *  recip()
 */

template<class A>
typename HIntLib::Private::QFB<A>::type
HIntLib::Private::QFB<A>::recip (const type& u) const
{
   if (is0 (u))  throwDivisionByZero();

   type result (u.den, u.num);
   a.mulByUnit (result.num, a.unitRecip (a.makeCanonical (result.den)));
   return result;
}


/**
 *  reciprocal()
 */

template<class A>
void
HIntLib::Private::QFB<A>::reciprocal (type& u) const
{
   if (is0 (u))  throwDivisionByZero();

   std::swap (u.num, u.den);
   a.mulByUnit (u.num, a.unitRecip (a.makeCanonical (u.den)));
}


/**
 *  print Short()
 */

template<typename A>
void
HIntLib::Private::QFB<A>::
printShort (std::ostream& o, const type& u, PrintShortFlag f) const
{
   if (! is0 (u) && ! a.is1 (u.den))
   {
      Printer ss (o);

      if (f & FIT_FOR_MUL)
      {
         if (o.flags() & o.showpos)
         {
            ss << '+';
            ss.unsetf (ss.showpos);
         }
         ss << '(';
      }
      // if showpos is set, keep it set
      a.printShort (ss, u.num, PrintShortFlag (f | FIT_FOR_MUL));
      ss << "/";
      ss.unsetf (ss.showpos);
      a.printShort (ss, u.den, PrintShortFlag (f | FIT_FOR_MUL));
      if (f & FIT_FOR_MUL)  ss << ')';
   }
   else
   {
      a.printShort (o, u.num, f);
   }
}

template<typename A>
void
HIntLib::Private::QFB<A>::print (std::ostream& o, const type& u) const
{
   Printer ss (o);

   printShort (ss, u);
   ss << ' ';
   printSuffix (ss);
}


/************  Quotient Field Base 2  -  integers  ***************************/


/**
 *  dbl()
 */

template<class A>
typename HIntLib::Private::QFB<A>::type
HIntLib::Private::QFB2<A,HIntLib::integer_tag>::dbl (const type& u) const
{
   return this->a.isEven (u.den) ?
      type (u.num, this->a.div (u.den, this->a.dbl (this->a.one()))) :
      type (this->a.dbl (u.num), u.den);
}


/**
 *  times2()
 */

template<class A>
void
HIntLib::Private::QFB2<A,HIntLib::integer_tag>::times2 (type& u) const
{
   if (this->a.isEven (u.den))
   {
      this->a.divBy (u.den, this->a.dbl (this->a.one()));
   }
   else
   {
      this->a.times2 (u.num);
   }
}


/**
 *  times()
 */

template<class A>
typename HIntLib::Private::QFB<A>::type
HIntLib::Private::QFB2<A,HIntLib::integer_tag>::
times (const type& u, unsigned k) const
{
   typedef typename A::type base_type;

   if (k == 0 || is0 (u))  return type();

   const A& aa (this->a);
   const base_type d = genGcd (aa, base_type (k), u.den);

   if (aa.isUnit (d) || aa.is1 (u.den))
   {
      return type (aa.times (u.num, k), u.den);
   }
   else
   {
      return type (aa.mul (u.num, aa.div (base_type(k), d)), aa.div (u.den, d));
   }
}


/**
 *  operator<<
 */

template<typename A>
std::ostream&
HIntLib::Private::
operator<< (std::ostream& o,
            const HIntLib::Private::QFB2<A,HIntLib::integer_tag>&)
{
   return o << "Q";
}

/************  Quotient Field Base 2  -  polynomials  ************************/


/**
 *  times()
 */

// for polynomials over fields, characteristic finite

template<class A>
typename HIntLib::Private::QFB<A>::type
HIntLib::Private::QFB2<A,HIntLib::polyoverfield_tag>::
timesImp (const type& u, unsigned k, char_prime) const
{
   if (k)
   {
      const typename A::type num = this->a.times (u.num, k);
      if (! this->a.is0 (num)) return type (num, u.den);
   }
   return type();
}


/**
 *  operator<<
 */

template<typename A>
std::ostream&
HIntLib::Private::
operator<< (std::ostream& o,
            const HIntLib::Private::QFB2<A,HIntLib::polyoverfield_tag>& a)
{
   Printer ss (o);

   const A& alg = a.getBaseAlgebra(); 
   ss << alg.getCoeffAlgebra() << '(';
   alg.printVariable (ss);
   ss << ')';

   return o;
}


#define HINTLIB_INSTANTIATE_QUOTIENTFIELD(X) \
namespace Private { \
template QFB<X >::type QFB<X >::toLowestTerms \
            (const base_type&, const base_type&) const; \
template QFB<X >::type QFB<X >::toLowestTermsAndNormalize \
            (const base_type&, const base_type&) const; \
template unsigned QFB<X >::index (const type&) const; \
template QFB<X >::type QFB<X >::element(unsigned) const; \
template QFB<X >::type QFB<X >::makeElement(const base_type&) const; \
template QFB<X >::type QFB<X >::add (const type&, const type&) const; \
template QFB<X >::type QFB<X >::sub (const type&, const type&) const; \
template QFB<X >::type QFB<X >::mul (const type&, const type&) const; \
template QFB<X >::type QFB<X >::div (const type&, const type&) const; \
template QFB<X >::type QFB<X >::recip (const type&) const; \
template void QFB<X >::reciprocal (type&) const; \
template unsigned QFB<X >::order (const type&) const; \
template void QFB<X >::printShort \
            (std::ostream&, const type&, PrintShortFlag) const; \
template void QFB<X >::print (std::ostream&, const type&) const; \
}

#define HINTLIB_INSTANTIATE_QUOTIENTFIELD_INT(X) \
   HINTLIB_INSTANTIATE_QUOTIENTFIELD(X) \
   namespace Private { \
   template QFB<X >::type QFB2<X,integer_tag>::times \
               (const type&, unsigned k) const; \
   template QFB<X >::type QFB2<X,integer_tag>::dbl (const type&) const; \
   template void QFB2<X,integer_tag>::times2 (type&) const; \
   template std::ostream& operator<< \
           (std::ostream&, const QFB2<X,integer_tag> &); \
   }

#define HINTLIB_INSTANTIATE_QUOTIENTFIELD_POL(X) \
   HINTLIB_INSTANTIATE_QUOTIENTFIELD(X) \
   namespace Private { \
   template QFB<X >::type QFB2<X,polyoverfield_tag>::timesImp \
               (const type&, unsigned k, X::char_category) const; \
   template std::ostream& operator<< \
               (std::ostream&, const QFB2<X,polyoverfield_tag> &); \
   }


