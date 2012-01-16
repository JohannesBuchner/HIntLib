/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schürer <rudolf.schuerer@sbg.ac.at>
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

#ifndef HINTLIB_INTEGER_RING_H
#define HINTLIB_INTEGER_RING_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <vector>
#include <utility>
#include <iosfwd>

#include <HIntLib/hlmath.h>
#include <HIntLib/algebra.h>
#include <HIntLib/prime.h>
#include <HIntLib/gcd.h>

namespace HIntLib
{

/**
 * ZRing
 *
 * Common base class for all integer rings
 */

namespace Private
{
   class ZRing { protected: ZRing() {} };

   std::ostream& operator<< (std::ostream&, const ZRing&);

}  // namespace Private


/**
 *  Integer Ring
 */

template <typename T = int>
class IntegerRing : public Private::ZRing
{
public:
   // types

   typedef T type;
   typedef int unit_type;
   typedef std::vector<std::pair<type,unsigned> > Factorization;

   // Algebra traits

   typedef integer_tag algebra_category; 
   typedef factor_tag primedetection_category;
   typedef infinite_tag size_category;
   typedef nopolynomial_tag polynomial_category;
   typedef char_zero char_category;
   typedef nozerodivisor_tag zerodivisor_category;

   // Properties of Z

   static unsigned size()  { return 0; }
   static unsigned characteristic()  { return 0; }
   static unsigned numUnits()  { return 2; }

   // Creation of elements

   static T one()  { return T(1); }
   static T element(unsigned);
   static unsigned index (const T&);

   // Properties of integers

   static bool is0 (T a)  { return a == T(0); }
   static bool is1 (T a)  { return a == T(1); }
   static bool isUnit (T a)  { return a == 1 || a == -1; }
   static bool isEven (T a)  { return ! (a & 1); }
   static bool isPrime (T a)  { return Prime::test (unsigned (abs(a))); }
   static bool isComposit (T a)
   {
      unsigned pos = abs (a);
      return pos > 3 && ! Prime::test (pos);
   }

   static unsigned numOfRemainders (T x)  { return abs (x); }
   static unsigned norm (T x)  { return abs (x); }

   // Additive arithmetic

   static T    add (const T& a, const T& b)  { return a +  b; }
   static void addTo (T& a,     const T& b)  { a += b; }

   static T sub (const T& a, const T& b)  { return a -  b; }
   static void subFrom (T& a,  const T& b)  { a -= b; }

   static T neg (const T& a)  { return     -a; }
   static void negate (T& a)    { a = -a; }

   static T  dbl (const T& a)  { return a << 1; }
   static void times2 (T& a)     { a <<= 1; }
   
   static T times (const T& a, unsigned k)  { return a * T(k); }

   // Multiplicative arithmetic

   static T mul (const T& a, const T& b)  { return a *  b; }
   static void mulBy (T& a,    const T& b)  { a *= b; }

   static T  sqr (const T& a)  { return a *  a; }
   static void square (T& a)     { a *= a; }
   
   static T power (const T& a, unsigned k)  { return powInt (a, k); }

   static void div (const T& a, const T& b, T& q, T& r)
      { q = a / b; r = a % b; }
   static T rem  (const T& a, const T& b)  { return a % b; }
   static T quot (const T& a, const T& b)  { return a / b; }
   static T div  (const T& a, const T& b)  { return a / b; }
   static void reduce   (T& a, const T& b)  { a %= b; }
   static void quotient (T& a, const T& b)  { a /= b; }
   static void divBy    (T& a, const T& b)  { a /= b; }

   static bool isAssociate (T a, T b)  { return a == b || a == -b; }
   static bool isAssociate (T a, T b, unit_type& u)
      {      if (a ==  b)  { u = 1;  return true; }
        else if (a == -b)  { u = -1; return true; }
        else return false;
      }
   static bool isDivisor (const T& a, const T& b)  { return a % b == 0; }
   static bool isDivisor (const T& a, const T& b, T& c)
      { if (a % b == 0) { c = a / b; return true; } else { return false; } }

   static unsigned order (T a)  { return is1(a) ? 1 : 0; }

   // I/O

   static void print (std::ostream &o, const T& a)  {  o << a; }
   static void printShort (
         std::ostream &o, const T& a, PrintShortFlag = PrintShortFlag())
      { print (o, a); }
   static void printSuffix (std::ostream &)  {}

   // units

   static unit_type unitElement (unsigned i)  { return i ? -1 : 1; }
   static unsigned unitIndex (unit_type a)  { return a > 0 ? 0 : 1; }
   
   static bool isCanonical (T a)  { return a >= 0; }
   static unit_type makeCanonical (T& a)
      { if (a >= 0) return 1; else { a = -a; return -1; } }

   static T  mulUnit   (const T& a, unit_type b)  { return a *  b; }
   static void mulByUnit (      T& a, unit_type b)  { a *= b; }

   static T       fromUnit (const unit_type& u)  { return u; }
   static unit_type toUnit (const T& a)  { return a; }

   static unit_type unitRecip (const unit_type& a)  { return a; }

   // primes

   class PrimeGenerator
   {
   public:
      PrimeGenerator (const IntegerRing<T> &) : prime (1) {}
      T next()  { return prime = Prime::next (prime + 1); }
   private:
      PrimeGenerator (const PrimeGenerator&);
      PrimeGenerator& operator= (const PrimeGenerator&);
      unsigned prime;
   };

   static unit_type factor (Factorization &, type);
   
   HINTLIB_TRIVIAL_DOMAIN_MEMBERS
};

// Specialization for powerMod()

namespace Private
{
   template<typename T> struct type_traits;
   template<>
   struct type_traits<int>   { typedef unsigned unsigned_type; };
   template<>
   struct type_traits<short> { typedef unsigned short unsigned_type; };
   template<>
   struct type_traits<long>  { typedef unsigned long unsigned_type; };
#ifdef HINTLIB_HAVE_LONG_LONG_INT
#ifdef HINTLIB_HAVE_UNSIGNED_LONG_LONG_INT
   template<>
   struct type_traits<long long> { typedef unsigned long long unsigned_type; };
#endif
#endif
}  // namespace Private

template<typename T>
inline
typename IntegerRing<T>::type
powerMod (
   const IntegerRing<T> &,
   const typename IntegerRing<T>::type& x, unsigned exponent,
   const typename IntegerRing<T>::type& m)
{
   typedef typename Private::type_traits<T>::unsigned_type TT;
   return powerMod (TT(x % m), exponent, TT(m));
}

// Specialization for genGcd()

template<typename T>
inline
bool
genIsCoprime (const IntegerRing<T>&,
     const typename IntegerRing<T>::type& u,
     const typename IntegerRing<T>::type& v)
{
   return isCoprime (u, v);
}

template<typename T>
inline
typename IntegerRing<T>::type
genGcd (const IntegerRing<T>&,
        const typename IntegerRing<T>::type& u,
        const typename IntegerRing<T>::type& v)
{
   return gcd (u, v);
}

template<typename T>
inline
typename IntegerRing<T>::type
genGcd (const IntegerRing<T>&,
        const typename IntegerRing<T>::type& u,
        const typename IntegerRing<T>::type& v,
              typename IntegerRing<T>::type& mu)
{
   return gcd (u, v, mu);
}

template<typename T>
inline
typename IntegerRing<T>::type
genGcd (const IntegerRing<T>&,
        const typename IntegerRing<T>::type& u,
        const typename IntegerRing<T>::type& v,
              typename IntegerRing<T>::type& mu,
              typename IntegerRing<T>::type& mv)
{
   return gcd (u, v, mu, mv);
}

}  // namespace HIntLib

#endif

