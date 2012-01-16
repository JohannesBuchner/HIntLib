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
#include <HIntLib/algebra.h>


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

   unsigned size() const  { return 0; }

   T zero() const  { return T(0); }
   T one() const  { return T(1); }

   bool is0 (T a) const  { return a == T(0); }
   bool is1 (T a) const  { return a == T(1); }

   T add (const T& a, const T& b) const  { return a +  b; }
   T& addTo (T& a,    const T& b) const  { return a += b; }

   T neg (const T& a) const  { return     -a; }
   T& negate (T& a) const    { return a = -a; }

   T sub (const T& a, const T& b) const  { return a -  b; }
   T& subFrom (T& a,  const T& b) const  { return a -= b; }

   T mul (const T& a, const T& b) const  { return a *  b; }
   T& mulBy (T& a,    const T& b) const  { return a *= b; }

   T times (const T& a, unsigned k) const  { return a * T(k); }
   T power (const T& a, unsigned k) const  { return powInt (a, k); }

   void print (std::ostream &o, const T& a) const  {  o << a; }
   void printShort (std::ostream &o, const T& a) const  { print (o, a); }
   void printSuffix (std::ostream &) const  {}
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

template <typename T = int>
class IntegerRing : public NumberRing<T>, public ZRing
{
public:
   typedef euclidean_tag algebra_category; 
   
   IntegerRing()
   {
      if (   ! std::numeric_limits<T>::is_signed
          || ! std::numeric_limits<T>::is_integer)
      {
         throw InvalidType ("IntegerRing");
      }
   } 

   T element(unsigned) const;
   unsigned index (T) const;

   void div (const T& a, const T& b, T& q, T& r) const { q = a / b; r = a % b; }
   T quot (const T& a, const T& b) const  { return a / b; }
   T rem  (const T& a, const T& b) const  { return a % b; }

   unsigned numOfRemainders (T x) const  { return abs (x); }
   unsigned norm (T x) const  { return abs (x); }

   bool isUnit  (T a) const  { return a == 1 || a == -1; }
   bool isPrime (T a) const;
   bool isIrreducible (T a) const  { return isPrime (a); }
   bool isComposit (T a) const;
   T unitRecip (T a) const  { return a; }
};


template <typename T = real>
class RealField : public NumberRing<T>, public RRing,
                  public TrivialFieldMembers<T>
{
public:
   typedef field_tag algebra_category; 

   RealField()
   {
      if (std::numeric_limits<T>::is_integer)  throw InvalidType ("RealField");
   } 

   T element(unsigned) const;
   unsigned index (T r) const;

   T recip (T a) const  { return T(1) / a; }

   void div (T a, T b, T& q, T& r) const  { q = a / b; r = T(0); }
   T div  (T a, T b) const  { return a / b; }
   T quot (T a, T b) const  { return a / b; }
   T& divBy (T& a, T b) const  { return a /= b; }

   bool isUnit  (T a) const  { return ! is0 (a); }
};

std::ostream & operator<< (std::ostream &, const ZRing &);
std::ostream & operator<< (std::ostream &, const RRing &);


}  // namespace HIntLib

#endif

