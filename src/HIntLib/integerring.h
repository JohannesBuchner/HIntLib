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


namespace HIntLib
{

template <class T>
class NumberRing
{
public:
   typedef T type;

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

   T times (const T& a, int k) const  { return a * T(k); }
   T power (const T& a, int k) const  { return powInt (a, k); }

   unsigned size() const  { return 0; }
};

class ZRing {};
class RRing {};

template <class T>
class IntegerRing : public NumberRing<T>, public ZRing
{
public:
   IntegerRing()
   {
      if (   ! std::numeric_limits<T>::is_signed
          || ! std::numeric_limits<T>::is_integer)
      {
         throw 1;
      }
   } 

   T element(unsigned i) const
      { return odd(i)  ?  T(i/2 + 1)  :  -T(i/2); }
   unsigned index (T x) const
      { return x > 0  ?  x*2 - 1 :  x * -2; }

   void div (const T& a, const T& b, T& q, T& r) const { q = a / b; r = a % b; }

   unsigned numOfRemainders (T x) const  { return abs (x); }

   bool isUnit  (T a) const  { return a == 1 || a == -1; }
   bool isPrime (T a) const;
   bool isIrreducible (T a) const  { return isPrime (a); }
   T unitRecip (T a) const  { return a; }
};


template <class T = real>
class RealField : public NumberRing<T>, public RRing
{
public:
   RealField()
   {
      if (std::numeric_limits<T>::is_integer)  throw 1;
   } 

   T element(unsigned i) const
      { if (i == 0)  return T(0);
        else return odd(i)  ?  T(i/2 + 1)  :  -T(i/2); }
   unsigned index (T r) const
      { int x = int (r); return x > 0  ?  x*2 - 1 :  x * -2; }

   T recip (T a) const  { return T(1) / a; }

   void div (T a, T b, T& q, T& r) const  { q = a / b; r = T(0); }
   T div (T a, T b) const  { return a /  b; }
   T& divBy (T& a, T b) const  { return a /= b; }

   bool isUnit  (T a) const  { return ! is0 (a); }
   bool isPrime (T) const  { return false; }
   bool isIrreducible (T) const { return false; }
};

std::ostream & operator<< (std::ostream &, const ZRing &);
std::ostream & operator<< (std::ostream &, const RRing &);


}  // namespace HIntLib

#endif

