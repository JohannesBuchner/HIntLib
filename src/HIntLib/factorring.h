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

#ifndef HINTLIB_FACTOR_RING_H
#define HINTLIB_FACTOR_RING_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <iosfwd>

#include <HIntLib/defaults.h>
#include <HIntLib/algebra.h>
#include <HIntLib/prime.h>

#ifdef HINTLIB_HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/fallback_limits.h>
#endif

namespace HIntLib
{
namespace Priv
{

/**
 * Modular Integer Ring Base
 */

class ModularIntegerRingBase
{
public:
   unsigned size() const     { return m; }
   unsigned modulus() const  { return m; }


protected:
   ModularIntegerRingBase (unsigned _m) : m(_m) {}
   void checkField() const;

   unsigned fieldRecip (unsigned) const;
   unsigned fieldOrder (unsigned) const;
   unsigned ringAdditiveOrder (unsigned) const;

   const unsigned m;

          void prnSuffix (std::ostream &) const;
   static void prnShort (std::ostream &, unsigned);
          void prn (std::ostream &, unsigned) const;
};

std::ostream& operator<< (std::ostream &, const ModularIntegerRingBase &);


/**
 * Modular Integer Ring
 */

template <class T>
class ModularIntegerRing : public ModularIntegerRingBase
{
protected:

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

public:

   typedef T type;
   typedef ring_tag algebra_category;
   typedef nopolynomial_tag polynomial_category;

   T one()  const  { return T(1); }

   bool is0 (T a) const  { return a == 0; }
   bool is1 (T a) const  { return a == 1; }

   T element(unsigned i) const  { return T(i); }
   unsigned index (T x) const   { return x; }

   T add (const T& a, const T& b) const
      { unsigned t = unsigned(a) + unsigned(b); return     (t>=m) ? t - m : t; }
   T& addTo (T& a,    const T& b) const
      { unsigned t = unsigned(a) + unsigned(b); return a = (t>=m) ? t - m : t; }

   T neg (const T& a) const  { return a != 0  ?  m - a  :  0; }
   T& negate (T& a) const    { if (a != 0) a = m - a;  return a; }

   T sub (const T& a, const T& b) const
      { int t = int(a) - int(b); return     (t<0) ? t+m : t; }
   T& subFrom (T& a,  const T& b) const
      { int t = int(a) - int(b); return a = (t<0) ? t+m : t; }

   T mul (const T& a, const T& b) const
      { return (unsigned(a) * unsigned(b)) % m; }
   T& mulBy (T& a,    const T& b) const
      { return a = (unsigned(a) * unsigned(b)) % m; }

   T times (const T& a, unsigned k) const
      { if (k >= m) k %= m;  return (a * k) % m; }
   T power (const T& a, unsigned k) const
      { return powerMod (unsigned(a), k, unsigned(m)); }

   unsigned additiveOrder (T a) const  { return ringAdditiveOrder (a); }

   void printSuffix (std::ostream &o) const  {  prnSuffix (o); }
   void printShort (std::ostream &o, const T& a) const  { prnShort (o, a); }
   void print (std::ostream &o, const T& a) const  { prn (o, a); }
};

}  // namespace Priv


/**
 *  Modular Arith
 */

template <class A>
class FactorRing
{
public:
   typedef typename A::type type;
   typedef ring_tag algebra_category;
   typedef nopolynomial_tag polynomial_category;

protected:
   const A a;
   const type m;

public:
   FactorRing (const A& _a, const type& _m) : a(_a), m(_m) {}

   type modulus() const  { return m; }
   A arithmetic() const  { return a; }
   unsigned size() const  { return a.numOfRemainders(m); }

   type one()  const  { return a.one(); }

   bool is0 (const type& u) const  { return a.is0(u); }
   bool is1 (const type& u) const  { return a.is1(u); }

   type element(unsigned i) const  { return a.rem (a.element(i), m); }
   unsigned index (const type& u) const   { return a.index (u); }

protected:
   typedef typename A::polynomial_category PC;

   // non-polynomials

   type  add (const type& u, const type& v, nopolynomial_tag) const
      { return     a.rem (a.add (u, v), m); }
   type& addTo (    type& u, const type& v, nopolynomial_tag) const
      { return u = a.rem (a.add (u, v), m); }

   type  neg (const type& u, nopolynomial_tag) const
      { return     a.rem (a.neg (u), m); }
   type& negate    (type& u, nopolynomial_tag) const
      { return u = a.rem (a.neg (u), m); }

   type  sub (const type& u, const type& v, nopolynomial_tag) const
      { return     a.rem (a.sub (u, v), m); }
   type& subFrom   (type& u, const type& v, nopolynomial_tag) const
      { return u = a.rem (a.sub (u, v), m); }

   type times (const type& u, unsigned k, nopolynomial_tag) const
      { return a.rem (a.times (u, k), m); }

   // polynomials

   type  add (const type& u, const type& v, polynomial_tag) const
      { return a.add   (u, v); }
   type& addTo (    type& u, const type& v, polynomial_tag) const
      { return a.addTo (u, v); }

   type  neg (const type& u, polynomial_tag) const  { return a.neg (u); }
   type& negate    (type& u, polynomial_tag) const  { return a.negate (u); }

   type  sub (const type& u, const type& v, polynomial_tag) const
      { return a.sub     (u, v); }
   type& subFrom   (type& u, const type& v, polynomial_tag) const
      { return a.subFrom (u, v); }

   type times (const type& u, unsigned k, polynomial_tag) const
      { return a.times (u, k); }

public:
   type  add (const type& u, const type& v) const { return add   (u,v, PC()); }
   type& addTo (    type& u, const type& v) const { return addTo (u,v, PC()); }

   type  neg (const type& u) const  { return neg    (u, PC()); }
   type& negate    (type& u) const  { return negate (u, PC()); }

   type  sub (const type& u, const type& v) const { return sub    (u,v, PC()); }
   type& subFrom   (type& u, const type& v) const { return subFrom(u,v, PC()); }

   type  mul (const type& u, const type& v) const
      { return     a.rem (a.mul (u, v), m); }
   type& mulBy     (type& u, const type& v) const
      { return u = a.rem (a.mul (u, v), m); }

   type times (const type& u, unsigned k) const  { return times (u,k, PC()); }
   type power (const type& u, unsigned k) const
      { return powerMod (a, u, k, m); }

   unsigned additiveOrder (const type& u) const;

   void print (std::ostream &, const type&) const;
   void printShort (std::ostream &, const type&) const;
   void printSuffix (std::ostream &) const;
};

template <>
class FactorRing<unsigned char>
   : public Priv::ModularIntegerRing<unsigned char>
{
public:
   FactorRing(unsigned char m)
      : Priv::ModularIntegerRing<unsigned char> (m) {}
};

template <>
class FactorRing<unsigned short>
   : public Priv::ModularIntegerRing<unsigned short>
{
public:
   FactorRing(unsigned short m)
      : Priv::ModularIntegerRing<unsigned short> (m) {}
};


/**
 *  operator<<
 */

template<class A>
std::ostream& operator<< (std::ostream &, const FactorRing<A> &);

template<>
inline
std::ostream& operator<< (
      std::ostream &o, const FactorRing<unsigned char> &a)
{
   return o << static_cast<const Priv::ModularIntegerRingBase&> (a);
}

template<>
inline
std::ostream& operator<< (
      std::ostream &o, const FactorRing<unsigned short> &a)
{
   return o << static_cast<const Priv::ModularIntegerRingBase&> (a);
}


/**
 *  Modular Integer Field
 */

namespace Priv
{
template<typename T>
class ModularIntegerField : public ModularIntegerRing<T>
{
protected:

   ModularIntegerField (const T _m) : ModularIntegerRing<T> (_m)
      { checkField(); }

public:
   typedef cyclic_tag algebra_category;
   typedef nopolynomial_tag polynomial_category;
   typedef T type;

   unsigned characteristic() const  { return m; }

   T recip (const T& x) const { return T (fieldRecip (x)); }

   T div  (const T& a, const T& b) const  { return mul (a, recip(b)); }
   T& divBy (T& a,  const T& b) const  { return mulBy (a, recip(b)); }

   unsigned order (T a) const    { return fieldOrder (a); }
   bool isPrimitiveElement (T a) const  { return isPrimitiveRoot (a, m); }
   T power (const T& a, unsigned k) const
      { return powerModReduce (unsigned(a), k, unsigned(m)); }

   HINTLIB_TRIVIAL_CYCLIC_MEMBERS
};

}  // namespace Priv


/**
 *  Modular Arith Field
 */

template <class A>
class FactorField : public FactorRing<A>
{
public:
   typedef gf_tag algebra_category;
   typedef nopolynomial_tag polynomial_category;
   typedef typename A::type type;

   FactorField (const A &_a, const type &_m) : FactorRing<A> (_a, _m) {}

   unsigned characteristic() const  { return a.characteristic(); }

   type recip (const type&) const;

   type div (const type& u, const type& v) const  { return mul (u, recip(v)); }
   type& divBy (type& u,  const type& v) const  { return mulBy (u, recip(v)); }

   unsigned order (const type& u) const;
   bool isPrimitiveElement (const type& u) const;
   unsigned extensionDegree() const
      { return Prime::factorPrimePowerPower (size()); }

   HINTLIB_TRIVIAL_DOMAIN_MEMBERS

private:
   void throwIfZero (const type& u) const;
};

template<>
class FactorField<unsigned char>
      : public Priv::ModularIntegerField<unsigned char>
{
public:
   FactorField (unsigned char m)
      : Priv::ModularIntegerField<unsigned char> (m) {}
};

template<>
class FactorField<unsigned short>
      : public Priv::ModularIntegerField<unsigned short>
{
public:
   FactorField (unsigned short m)
      : Priv::ModularIntegerField<unsigned short> (m) {}
};

}  // namespace HIntLib

#endif

