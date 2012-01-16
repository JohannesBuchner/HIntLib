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

#ifndef HINTLIB_INTEGER_RING_H
#define HINTLIB_INTEGER_RING_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <iosfwd>

#include <HIntLib/mymath.h>
#include <HIntLib/prime.h>
#include <HIntLib/algebra.h>
#include <HIntLib/exception.h>


namespace HIntLib
{

template <typename T>
class NumberRing
{
protected:
   NumberRing() {}
public:
   typedef T type;
   typedef nopolynomial_tag polynomial_category;

   static unsigned size()  { return 0; }
   static unsigned characteristic()  { return 0; }

   static bool is1 (T a)  { return a == T(1); }

   static T one()  { return T(1); }

   static T add (const T& a, const T& b)  { return a +  b; }
   static T& addTo (T& a,    const T& b)  { return a += b; }

   static T neg (const T& a)  { return     -a; }
   static T& negate (T& a)    { return a = -a; }

   static T sub (const T& a, const T& b)  { return a -  b; }
   static T& subFrom (T& a,  const T& b)  { return a -= b; }

   static T mul (const T& a, const T& b)  { return a *  b; }
   static T& mulBy (T& a,    const T& b)  { return a *= b; }

   static T times (const T& a, unsigned k)  { return a * T(k); }
   static T power (const T& a, unsigned k)  { return powInt (a, k); }

   static void print (std::ostream &o, const T& a)  {  o << a; }
   static void printShort (std::ostream &o, const T& a)  { print (o, a); }
   static void printSuffix (std::ostream &)  {}
};

class ZRing
{
protected:
   ZRing () {}
};

class RRing
{
protected:
   RRing () {}
};


/**
 *  Integer Ring
 */

template <typename T = int>
class IntegerRing : public NumberRing<T>, public ZRing
{
public:
   typedef integer_tag algebra_category; 
   typedef T type;
   typedef type unit_type;

   class PrimeGenerator
   {
   public:
      PrimeGenerator (const IntegerRing<T> &) : prime (1) {}
      T next()  { return prime = Prime::next (prime + 1); }
   private:
      PrimeGenerator ();
      PrimeGenerator& operator= (const PrimeGenerator&);
      unsigned prime;
   };
   
   IntegerRing()
   {
      if (   ! std::numeric_limits<T>::is_signed
          || ! std::numeric_limits<T>::is_integer)
      {
         throw InvalidType ("IntegerRing");
      }
   } 

   static bool is0 (T a)  { return a == T(0); }

   static T element(unsigned);
   static unsigned index (T);
   static bool isCanonical (T a)  { return a >= 0; }
   static unit_type makeCanonical (T& a)
      { if (a >= 0) return 1; else { a = -a; return -1; } }

   static unsigned numUnits()  { return 2; }
   static unit_type unitElement (unsigned i)  { return i ? -1 : 1; }
   static unsigned unitIndex (unit_type a)  { return a > 0 ? 0 : 1; }
   
   static void div (const T& a, const T& b, T& q, T& r)
      { q = a / b; r = a % b; }
   static T quot (const T& a, const T& b)  { return a / b; }
   static T rem  (const T& a, const T& b)  { return a % b; }

   static unsigned numOfRemainders (T x)  { return abs (x); }
   static unsigned norm (T x)  { return abs (x); }

   static bool isUnit  (T a)  { return a == 1 || a == -1; }
   static bool isPrime (T a);
   static bool isComposit (T a);
   static unit_type unitRecip (unit_type a)  { return a; }

   static T  mulUnit   (const T& a, const unit_type& b)  { return a *  b; }
   static T& mulByUnit (      T& a, const unit_type& b)  { return a *= b; }

   static T       fromUnit (const unit_type& u)  { return u; }
   static unit_type toUnit (const T& a)  { return a; }

   static unsigned order (T a)  { return is1(a) ? 1 : 0; }

   HINTLIB_TRIVIAL_DOMAIN_MEMBERS
};


/**
 *  Real
 */

template<typename T>
class Real
{
public:
   Real () : x(0) {}
   Real (T a) : x(a) {}

   operator T () const  { return x; }

   Real operator- () const  { return Real (-x); }
   
   Real operator+ (const Real& a) const  { return Real (x + a.x); }
   Real operator- (const Real& a) const  { return Real (x - a.x); }
   Real operator* (const Real& a) const  { return Real (x * a.x); }
   Real operator/ (const Real& a) const  { return Real (x / a.x); }

   Real& operator+= (const Real& a)  { x += a.x; return *this; }
   Real& operator-= (const Real& a)  { x -= a.x; return *this; }
   Real& operator*= (const Real& a)  { x *= a.x; return *this; }
   Real& operator/= (const Real& a)  { x /= a.x; return *this; }

private:
   T x;
};

template<typename T>
inline
bool operator== (const Real<T>& a, const Real<T>& b)
{
   return approx (T(a), T(b));
}

template<typename T>
inline
Real<T> powInt (const Real<T>& a, unsigned k)
{
   return Real<T> (powInt(T(a), k));
}


/**
 *  Real Field
 */

template<typename T = real>
class RealField : public NumberRing<Real<T> >, public RRing
{
public:
   typedef real_tag algebra_category; 
   typedef Real<T> type;

   RealField()
   {
      if (std::numeric_limits<T>::is_integer)  throw InvalidType ("RealField");
   } 

   bool is0 (const type& a) const
      { return abs(a) < std::numeric_limits<T>::epsilon() * 100; }

   type element(unsigned) const;
   unsigned index (type r) const;

   type recip (const type& a)     const  { return type(1) / a; }

   type div  (const type& a, const type& b) const  { return a / b; }
   type& divBy (type& a, const type& b) const  { return a /= b; }

   unsigned order (const type& a) const  { return is1(a) ? 1 : 0; }

   HINTLIB_TRIVIAL_DOMAIN_MEMBERS
};

std::ostream & operator<< (std::ostream &, const ZRing &);
std::ostream & operator<< (std::ostream &, const RRing &);


}  // namespace HIntLib

#endif

