/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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
 *  printShort()
 */

template<typename A>
void
HIntLib::Private::QFB2<A,HIntLib::integer_tag>::
printShort (std::ostream& o, const type& u, PrintShortFlag f) const
{
   const A& a (this->a);

   if (is0 (u) || a.is1 (u.den))
   {
      a.printShort (o, u.num, f);
      return;
   }

   Printer ss (o);
#ifdef HINTLIB_ENCODING_LOCALE
   const bool utf8 = ss.utf8();
#endif

#ifdef HINTLIB_UTF8_SELECT
   typedef typename A::type base_type;
   base_type posNum (u.num);
   int sign = a.makeCanonical (posNum);
   const base_type max = base_type(9);

   if (u.den < max && posNum < u.den)
   {
      int num = int(posNum);
      int den = int(u.den);
      const char* x = 0;

#ifdef HINTLIB_ENCODING_LOCALE
      if (utf8)
#endif
#if defined(HINTLIB_ENCODING_LOCALE) || defined(HINTLIB_ENCODING_UTF8)
      {
         switch (num * 10 + den)
         {
         case 12: x = "\xc2\xbd";     break;

#if HINTLIB_CHARACTER_SET >= 4
         case 13: x = "\xe2\x85\x93"; break;
         case 23: x = "\xe2\x85\x94"; break;
#endif

         case 14: x = "\xc2\xbc";     break;
         case 34: x = "\xc2\xbe";     break;

#if HINTLIB_CHARACTER_SET >= 4
         case 15: x = "\xe2\x85\x95"; break;
         case 25: x = "\xe2\x85\x96"; break;
         case 35: x = "\xe2\x85\x97"; break;
         case 45: x = "\xe2\x85\x98"; break;

         case 16: x = "\xe2\x85\x99"; break;
         case 56: x = "\xe2\x85\x9a"; break;
#endif

#if HINTLIB_CHARACTER_SET >= 3
         case 18: x = "\xe2\x85\x9b"; break;
         case 38: x = "\xe2\x85\x9c"; break;
         case 58: x = "\xe2\x85\x9d"; break;
         case 78: x = "\xe2\x85\x9e"; break;
#endif
         }
      }
#endif
#ifdef HINTLIB_ENCODING_LOCALE
      else
#endif
#if defined(HINTLIB_ENCODING_LOCALE) || defined(HINTLIB_ENCODING_LATIN1)
      {
         switch (num * 10 + den)
         {
         case 12: x = "\xbd"; break;
         case 14: x = "\xbc"; break;
         case 34: x = "\xbe"; break;
         }
      }
#endif

      if (x)
      {
         if (sign == -1)
         {
            ss.minusSign();
         }
         else if (ss.flags() & ss.showpos)  ss << '+';

         ss << x;
         return;
      }
   }
#endif

   if (f & FIT_FOR_MUL)
   {
      if (ss.flags() & ss.showpos)
      {
         ss << '+';
         ss.unsetf (ss.showpos);
      }
      ss << '(';
   }
   // if showpos is set, keep it set
   a.printShort (ss, u.num, PrintShortFlag (f | FIT_FOR_MUL));
#if HINTLIB_CHARACTER_SET >= 3 && defined HINTLIB_UTF8_SELECT
   // DIVISION SLASH
   HINTLIB_UTF8_SELECT(ss.utf8(), ss << "\xe2\x88\x95", ss << '/' )
#else
   ss << '/';
#endif
   ss.unsetf (ss.showpos);
   a.printShort (ss, u.den, PrintShortFlag (f | FIT_FOR_MUL));
   if (f & FIT_FOR_MUL)  ss << ')';
}

#ifdef HINTLIB_BUILD_WCHAR
template<typename A>
void
HIntLib::Private::QFB2<A,HIntLib::integer_tag>::
printShort (std::wostream& o, const type& u, PrintShortFlag f) const
{
   const A& a (this->a);

   if (is0 (u) || a.is1 (u.den))
   {
      a.printShort (o, u.num, f);
      return;
   }

   WPrinter ss (o);

#if HINTLIB_CHARACTER_SET >= 2
   typedef typename A::type base_type;
   base_type posNum (u.num);
   int sign = a.makeCanonical (posNum);
   const base_type max = base_type(9);

   if (u.den < max && posNum < u.den)
   {
      int num = int(posNum);
      int den = int(u.den);
      wchar_t x = 0;

      switch (num * 10 + den)
      {
      case 12: x = L'\x00bd'; break;

#if HINTLIB_CHARACTER_SET >= 4
      case 13: x = L'\x2153'; break;
      case 23: x = L'\x2154'; break;
#endif

      case 14: x = L'\x00bc'; break;
      case 34: x = L'\x00be'; break;

#if HINTLIB_CHARACTER_SET >= 4
      case 15: x = L'\x2155'; break;
      case 25: x = L'\x2156'; break;
      case 35: x = L'\x2157'; break;
      case 45: x = L'\x2158'; break;

      case 16: x = L'\x2159'; break;
      case 56: x = L'\x215A'; break;
#endif

#if HINTLIB_CHARACTER_SET >= 3
      case 18: x = L'\x215B'; break;
      case 38: x = L'\x215C'; break;
      case 58: x = L'\x215D'; break;
      case 78: x = L'\x215E'; break;
#endif
      }

      if (x)
      {
         if (sign == -1)
         {
            ss.minusSign();
         }
         else if (ss.flags() & ss.showpos)  ss << L'+';

         ss << x;
         return;
      }
   }
#endif

   if (f & FIT_FOR_MUL)
   {
      if (ss.flags() & ss.showpos)
      {
         ss << L'+';
         ss.unsetf (ss.showpos);
      }
      ss << L'(';
   }
   // if showpos is set, keep it set
   a.printShort (ss, u.num, PrintShortFlag (f | FIT_FOR_MUL));
#if HINTLIB_CHARACTER_SET >= 3
   ss << L'\x2215';  // DIVISION SLASH
#else
   ss << L'/';
#endif
   ss.unsetf (ss.showpos);
   a.printShort (ss, u.den, PrintShortFlag (f | FIT_FOR_MUL));
   if (f & FIT_FOR_MUL)  ss << L')';
}
#endif


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
 *  printShort()
 */

template<typename A>
void
HIntLib::Private::QFB2<A,HIntLib::polyoverfield_tag>::
printShort (std::ostream& o, const type& u, PrintShortFlag f) const
{
   const A& a (this->a);

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
#if HINTLIB_CHARACTER_SET >= 3 && defined HINTLIB_UTF8_SELECT
      // DIVISION SLASH
      HINTLIB_UTF8_SELECT(ss.utf8(), ss << "\xe2\x88\x95", ss << '/' )
#else
      ss << '/';
#endif
      ss.unsetf (ss.showpos);
      a.printShort (ss, u.den, PrintShortFlag (f | FIT_FOR_MUL));
      if (f & FIT_FOR_MUL)  ss << ')';
   }
   else
   {
      a.printShort (o, u.num, f);
   }
}

#ifdef HINTLIB_BUILD_WCHAR
template<typename A>
void
HIntLib::Private::QFB2<A,HIntLib::polyoverfield_tag>::
printShort (std::wostream& o, const type& u, PrintShortFlag f) const
{
   const A& a (this->a);

   if (! is0 (u) && ! a.is1 (u.den))
   {
      WPrinter ss (o);

      if (f & FIT_FOR_MUL)
      {
         if (o.flags() & o.showpos)
         {
            ss << L'+';
            ss.unsetf (ss.showpos);
         }
         ss << L'(';
      }
      // if showpos is set, keep it set
      a.printShort (ss, u.num, PrintShortFlag (f | FIT_FOR_MUL));
#if HINTLIB_CHARACTER_SET >= 3
      ss << L'\x2215';  // DIVISION SLASH
#else
      ss << L'/';
#endif
      ss.unsetf (ss.showpos);
      a.printShort (ss, u.den, PrintShortFlag (f | FIT_FOR_MUL));
      if (f & FIT_FOR_MUL)  ss << L')';
   }
   else
   {
      a.printShort (o, u.num, f);
   }
}
#endif


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

#ifdef HINTLIB_BUILD_WCHAR
template<typename A>
std::wostream&
HIntLib::Private::
operator<< (std::wostream& o,
            const HIntLib::Private::QFB2<A,HIntLib::polyoverfield_tag>& a)
{
   WPrinter ss (o);

   const A& alg = a.getBaseAlgebra(); 
   ss << alg.getCoeffAlgebra() << L'(';
   alg.printVariable (ss);
   ss << L')';

   return o;
}
#endif


/*****************  Quotient Field *******************************************/

/**
 *  print()
 */

template<typename A>
void
HIntLib::QuotientField<A>::print (std::ostream& o, const type& u) const
{
   Private::Printer ss (o);

   printShort (ss, u);
   ss << ' ';
   this->printSuffix (ss);
}

#ifdef HINTLIB_BUILD_WCHAR
template<typename A>
void
HIntLib::QuotientField<A>::print (std::wostream& o, const type& u) const
{
   Private::WPrinter ss (o);

   printShort (ss, u);
   ss << L' ';
   this->printSuffix (ss);
}
#endif


#ifdef HINTLIB_BUILD_WCHAR
#define HINTLIB_INSTANTIATE_QUOTIENTFIELD_W(X) \
   template void QuotientField<X >::print (std::wostream&, const type&) const;
#else
#define HINTLIB_INSTANTIATE_QUOTIENTFIELD_W(X)
#endif

#define HINTLIB_INSTANTIATE_QUOTIENTFIELD(X) \
   HINTLIB_INSTANTIATE_QUOTIENTFIELD_W(X) \
   template void QuotientField<X >::print (std::ostream&, const type&) const; \
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
   }

#ifdef HINTLIB_BUILD_WCHAR
#define HINTLIB_INSTANTIATE_QUOTIENTFIELD_INT_W(X) \
   template void QFB2<X,integer_tag>::printShort \
               (std::wostream&, const type&, PrintShortFlag) const;
#else
#define HINTLIB_INSTANTIATE_QUOTIENTFIELD_INT_W(X)
#endif

#define HINTLIB_INSTANTIATE_QUOTIENTFIELD_INT(X) \
   HINTLIB_INSTANTIATE_QUOTIENTFIELD(X) \
   namespace Private { \
   HINTLIB_INSTANTIATE_QUOTIENTFIELD_INT_W(X) \
   template void QFB2<X,integer_tag>::printShort \
               (std::ostream&, const type&, PrintShortFlag) const; \
   template QFB<X >::type QFB2<X,integer_tag>::times \
               (const type&, unsigned k) const; \
   template QFB<X >::type QFB2<X,integer_tag>::dbl (const type&) const; \
   template void QFB2<X,integer_tag>::times2 (type&) const; \
   }

#ifdef HINTLIB_BUILD_WCHAR
#define HINTLIB_INSTANTIATE_QUOTIENTFIELD_POL_W(X) \
   template std::wostream& operator<< \
               (std::wostream&, const QFB2<X,polyoverfield_tag> &); \
   template void QFB2<X,polyoverfield_tag>::printShort \
               (std::wostream&, const type&, PrintShortFlag) const;
#else
#define HINTLIB_INSTANTIATE_QUOTIENTFIELD_POL_W(X)
#endif

#define HINTLIB_INSTANTIATE_QUOTIENTFIELD_POL(X) \
   HINTLIB_INSTANTIATE_QUOTIENTFIELD(X) \
   namespace Private { \
   HINTLIB_INSTANTIATE_QUOTIENTFIELD_POL_W(X) \
   template QFB<X >::type QFB2<X,polyoverfield_tag>::timesImp \
               (const type&, unsigned k, X::char_category) const; \
   template std::ostream& operator<< \
               (std::ostream&, const QFB2<X,polyoverfield_tag> &); \
   template void QFB2<X,polyoverfield_tag>::printShort \
               (std::ostream&, const type&, PrintShortFlag) const; \
   }


