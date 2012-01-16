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

#ifndef HINTLIB_MODULAR_ARITHMETIC_H
#define HINTLIB_MODULAR_ARITHMETIC_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <iosfwd>

#include <HIntLib/mymath.h>


namespace HIntLib
{

/**
 * Modular Integer Ring
 */

class ModularIntegerRingBase
{
public:
   unsigned size() const  { return m; }
   unsigned modulus() const  { return m; }

protected:
   ModularIntegerRingBase (unsigned _m) : m(_m) {}
   const unsigned m;
};

std::ostream& operator<< (std::ostream &, const ModularIntegerRingBase &);

template <class T>
class ModularIntegerRing : public ModularIntegerRingBase
{
public:

   typedef T type;

   ModularIntegerRing(const T _m)
      : ModularIntegerRingBase (unsigned (_m))
   {
      if (   std::numeric_limits<T>::is_signed
          || ! std::numeric_limits<T>::is_integer
          || 2 * std::numeric_limits<T>::digits
               > std::numeric_limits<unsigned>::digits)
      {
         throw 1;
      }
   } 

   T zero() const  { return T(0); }
   T one()  const  { return T(1); }

   bool is0 (T a) const  { return a == 0; }
   bool is1 (T a) const  { return a == 1; }

   T element(unsigned i) const  { return T(i); }
   unsigned index (T x) const   { return x; }

   T add (const T& a, const T& b) const
      { return     (unsigned(a) + unsigned(b)) % m; }
   T& addTo (T& a,    const T& b) const
      { return a = (unsigned(a) + unsigned(b)) % m; }

   T neg (const T& a) const  { return a > 0  ?  m - a  :  0; }
   T& negate (T& a) const    { if (a == 0)  return a; return a = m - a; }

   T sub (const T& a, const T& b) const { return   (m + int(a) - int(b)) % m; }
   T& subFrom (T& a,  const T& b) const { return a=(m + int(a) - int(b)) % m; }

   T mul (const T& a, const T& b) const
      { return (unsigned(a) * unsigned(b)) % m; }
   T& mulBy (T& a,    const T& b) const
      { return a = (unsigned(a) * unsigned(b)) % m; }

   T times (const T& a, unsigned k) const  { return (a * (k % m)) % m; }
   T power (const T& a, unsigned k) const
      { return powerMod (unsigned(a), k, unsigned(m)); }
};


/**
 *  Modular Arithmetic Ring
 */

template <class A>
class ModularArithmeticRing
{
private:
   typedef typename A::type T;

protected:
   const A a;
   const T m;

   T rem (const T& x) const;

public:
   typedef typename A::type type;

   ModularArithmeticRing(const A _a, const T _m) : a(_a), m(_m) {}

   T modulus() const  { return m; }
   A arithmetic() const  { return a; }

   T zero() const  { return a.zero(); }
   T one()  const  { return a.one(); }

   bool is0 (const T& u) const  { return a.is0(u); }
   bool is1 (const T& u) const  { return a.is1(u); }

   T element(unsigned i) const  { return rem (a.element(i)); }
   unsigned index (const T& u) const   { return a.index (u); }

   T  add (const T& u, const T& v) const  { return     rem (a.add (u, v)); }
   T& addTo (    T& u, const T& v) const  { return u = rem (a.add (u, v)); }

   T  neg (const T& u) const  { return     rem (a.neg (u)); }
   T& negate    (T& u) const  { return u = rem (a.neg (u)); }

   T  sub (const T& u, const T& v) const { return     rem (a.sub (u, v)); }
   T& subFrom   (T& u, const T& v) const { return u = rem (a.sub (u, v)); }

   T  mul (const T& u, const T& v) const { return     rem (a.mul (u, v)); }
   T& mulBy     (T& u, const T& v) const { return u = rem (a.mul (u, v)); }

   T times (const T& u, unsigned k) const  { return rem (a.times (u, k)); }
   T power (const T& u, unsigned k) const  { return powerMod (a, u, k, m); }
   unsigned size() const  { return a.numOfRemainders(m); }
};

template<class A>
std::ostream& operator<< (std::ostream &, const ModularArithmeticRing<A> &);


/**
 *  Modular Integer Field
 */

template <class T>
class ModularIntegerField : public ModularIntegerRing<T>
{
public:
   ModularIntegerField (const T);

   T recip (const T&) const;
   T unitRecip (const T& u) const  { return recip (u); }

   void div (const T& a, const T& b, T& q, T& r) const
      { q = mul (a, recip (b)); r = 0; }
   T div (const T& a, const T& b) const  { return mul (a, recip(b)); }
   T& divBy (T& a,  const T& b) const  { return mulBy (a, recip(b)); }

   bool isUnit  (T a) const  { return a != 0; }
   bool isPrime (T)   const  { return false; }
   bool isIrreducible (T) const { return false; }

   T power (const T& a, unsigned k) const
      { return a ? powerMod (unsigned(a), k % (m-1), unsigned(m)) : 0; }
};


/**
 *  Modular Arithmetic Field
 */

template <class A>
class ModularArithmeticField : public ModularArithmeticRing<A>
{
private:
   typedef typename A::type T;

public:
   ModularArithmeticField (const A &_a, const T &_m)
      : ModularArithmeticRing<A> (_a, _m) {}
   // ModularArithmeticField (unsigned deg, unsigned n = 0);

   T recip (const T&) const;
   T unitRecip (const T& u) const  { return recip (u); }

   void div (const T& u, const T& v, T& q, T& r) const
      { q = mul (u, recip (v)); r = a.zero(); }
   T div (const T& u, const T& v) const  { return mul (u, recip(v)); }
   T& divBy (T& u,  const T& v) const  { return mulBy (u, recip(v)); }

   bool isUnit  (const T &u) const  { return ! is0 (u); }
   bool isPrime (const T &)  const  { return false; }
   bool isIrreducible (const T &) const { return false; }
};

}  // namespace HIntLib

#endif

