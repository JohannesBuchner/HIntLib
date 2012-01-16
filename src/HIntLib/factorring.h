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
#include <vector>
#include <utility>

#include <HIntLib/algebra.h>
#include <HIntLib/gcd.h>

namespace HIntLib
{
namespace Private
{

/*************  Factor Base  ***************************************************/


/**
 *  Factor BB
 */

template <class A>
class FactorBB : protected A
{
public:
   typedef typename A::type type;
   typedef nopolynomial_tag polynomial_category;

protected:
   FactorBB (const A& a, const type& modulus) : A(a), m(modulus) {}
   ~FactorBB();
   void throwIfZero (const type& u) const;

   const type m;

public:
   const type& modulus() const  { return m; }
   const A& arithmetic() const  { return *this; }
   unsigned size() const  { return numOfRemainders (m); }

   using A::one;
   using A::is0;
   using A::is1;

   type  mul (const type& u, const type& v) const
      { return     A::rem (A::mul (u, v), m); }
   void mulBy     (type& u, const type& v) const
      { A::mulBy (u, v); A::reduce (u, m); }

   type sqr (const type& u) const  { return A::rem (A::sqr (u), m); }
   void square (type& u) const    { A::square (u); A::reduce(u, m); }
   type power (const type& u, unsigned k) const
      { return powerMod (arithmetic(), u, k, m); }

   using A::printShort;
   void print (std::ostream &, const type&) const;
   void printSuffix (std::ostream &) const;
};


/**
 *  operator<<
 */

template<typename A>
std::ostream& operator<< (std::ostream &, const FactorBB<A> &);


/**
 *  Factor Poly B
 */

template<class A, typename SC> class FactorPolyB;

// finite base field

template<class A>
class FactorPolyB<A,finite_tag> : public FactorBB<A>
{
public:
   typedef finite_tag size_category;

   using A::element;
   using A::index;

protected:
   FactorPolyB (const A& a, const typename A::type& modulus)
      : FactorBB<A> (a, modulus) {}
};

// infinite base field

template<class A>
class FactorPolyB<A,infinite_tag> : public FactorBB<A>
{
public:
   typedef typename FactorBB<A>::type type;
   typedef infinite_tag size_category;

   type element (unsigned i) const;
   unsigned index (const type& u) const;

protected:
   FactorPolyB (const A& a, const type& modulus) : FactorBB<A> (a, modulus) {}
};


/**
 *  Factor Poly <A>
 */

template<class A>
class FactorPoly
   : public FactorPolyB<A,typename A::coeff_algebra::size_category>
{
public:
   typedef typename A::type type;
   typedef typename A::char_category char_category;

   using A::add;
   using A::addTo;
   using A::neg;
   using A::negate;
   using A::sub;
   using A::subFrom;
   using A::dbl;
   using A::times2;
   using A::times;
   using A::additiveOrder;
   using A::characteristic;

protected:
   FactorPoly (const A& a, const type& modulus)
      : FactorPolyB<A,typename A::coeff_algebra::size_category> (a, modulus) {}

   void recipImp (const type& x, type& r) const;
};


/**
 *  Factor Integer <A>
 */

template<class A>
class FactorInteger : public FactorBB<A>
{
public:
   typedef typename A::type type;
   typedef finite_tag size_category;

   type element (unsigned i) const  { return int (i); }
   unsigned index (const type& u) const  { return int (u); }

   type add (const type& u, const type& v) const
      { type t = A::add(u,v); return (t >= this->m) ? A::sub(t, this->m) : t; }
   void addTo (    type& u, const type& v) const
      { A::addTo (u,v); if (u >= this->m) A::subFrom (u, this->m); }

   type neg (const type& u) const
      { return ! is0(u) ? A::sub(this->m, u) : type(); }
   void negate    (type& u) const  { if (! is0(u)) u = A::sub(this->m, u); }

   type  sub (const type& u, const type& v) const
      { type t = A::sub(u,v); return (t<0) ? A::add(t, this->m) : t; }
   void subFrom   (type& u, const type& v) const
      { A::subFrom (u,v); if (u < 0)  A::addTo (u, this->m); }

   type dbl (const type& u) const
      { type t = A::dbl(u); return t >= this->m ? A::sub (t, this->m) : t; } 
   void times2 (type& u) const
      { A::times2 (u); if (u >= this->m) A::subFrom (u, this->m); }
   type times (const type& u, unsigned k) const
      { return A::rem (A::times (u, k), this->m); }

protected:
   FactorInteger (const A& a, const type& modulus) : FactorBB<A> (a, modulus) {}

   void recipImp (const type& x, type& r) const;
};


/*************  Factor Ring  *************************************************/


/**
 *  Factor Ring B <A,AC,PC>
 */

template<class A, typename AC> class FactorRingB;

// polynomial

template<class A>
class FactorRingB<A,polyoverfield_tag> : public FactorPoly<A>
{
public:
   unsigned numNilpotents() const;

protected:
   FactorRingB (const A& a, const typename A::type& modulus)
      : FactorPoly<A> (a, modulus)
   {
      init (typename A::coeff_algebra::size_category());
   }

   typename A::type nilradical;
   unsigned numberOfUnits;

private:
   void init (infinite_tag);
   void init (finite_tag);
};

// integers

template<class A>
class FactorRingB<A,integer_tag> : public FactorInteger<A>
{
public:
   typedef char_non char_category;

   unsigned additiveOrder (const typename A::type& u) const;
   unsigned numNilpotents () const
      { return int (A::div (this->m, nilradical)); }

protected:
   FactorRingB (const A&, const typename A::type&);

   typename A::type nilradical;
   unsigned numberOfUnits;
};

}  // namespace Private


/**
 *  Factor Ring
 */

template<class A>
class FactorRing
   : public Private::FactorRingB<A,typename A::algebra_category>
{
public:
   typedef typename A::type type;
   typedef type unit_type;
   typedef ring_tag algebra_category;
   typedef zerodivisor_tag zerodivisor_category;
   
   FactorRing (const A& a, const type& modulus)
      : Private::FactorRingB<A,typename A::algebra_category> (a, modulus) {}

   unsigned numUnits() const  { return this->numberOfUnits; }

   bool isUnit (const type& u) const
      { return genIsCoprime (this->arithmetic(), u, this->m); }
   bool isZerodivisor (const type& u) const
      { return ! genIsCoprime (this->arithmetic(), u, this->m); }
   bool isNilpotent (const type& u) const
      { return A::isDivisor (u, this->nilradical); }

   type unitRecip (const type& u) const  { type r; recipImp (u,r); return r; }

   type mulUnit (const type& v, const unit_type& u) const { return mul  (v,u); }
   void mulByUnit     (type& v, const unit_type& u) const { mulBy(v,u); }

   const      type& fromUnit (const unit_type& unit) const  { return unit; }
   const unit_type&   toUnit (const      type& unit) const  { return unit; }
};


/*************  Factor Field  ************************************************/


namespace Private
{

/**
 *  Factor Field Poly B
 */

template<class A, typename SC> class FactorFieldPolyB;

// finite base field

template<class A>
class FactorFieldPolyB<A,finite_tag> : public FactorPoly<A>
{
public:
   typedef typename A::type type;

   unsigned order (const type&) const;
   bool isPrimitiveElement (const type& u) const;

protected:
   FactorFieldPolyB (const A& a, const type& modulus)
      : FactorPoly<A> (a, modulus) {}

private:
   typedef std::vector<std::pair<unsigned, unsigned> > Fac;
   typedef Fac::const_iterator FacI;

   mutable Fac facNMinus1;
};

// infinite base field

template<class A>
class FactorFieldPolyB<A,infinite_tag> : public FactorPoly<A>
{
public:
   typedef typename A::type type;

   unsigned order (const type&) const;

protected:
   FactorFieldPolyB (const A& a, const type& modulus)
      : FactorPoly<A> (a, modulus) {}
};


/**
 *  Factor Field B <A,AC,PC>
 */

template<class A, typename AC> class FactorFieldB;

// polynomial

template<class A>
class FactorFieldB<A,polyoverfield_tag>
   : public FactorFieldPolyB<A,typename A::coeff_algebra::size_category>
{
public:
   typedef typename A::type type;
   typedef typename A::coeff_algebra::algebra_category::F algebra_category;

   FactorFieldB (const A& a, const type& modulus)
      : FactorFieldPolyB<A,typename A::coeff_algebra::size_category>
           (a, modulus)
   {}

   unsigned extensionDegree() const
      { return this->getCoeffAlgebra().extensionDegree() * this->m.degree(); }

   type frobenius (const type& u) const
      { return power (u, this->characteristic()); }
   type invFrobenius (const type& u) const
   {
      return power (u, powInt (this->characteristic(), extensionDegree() - 1));
   }
};

// integers

template<class A>
class FactorFieldB<A,integer_tag> : public FactorInteger<A>
{
public:
   typedef typename A::type type;
   typedef cyclic_tag algebra_category;
   typedef char_prime char_category;

   FactorFieldB (const A& a, const type& modulus)
      : FactorInteger<A> (a, modulus) {}

   unsigned order (const type&) const;
   bool isPrimitiveElement (const type& u) const;

   unsigned characteristic() const  { return int (this->m); }

   HINTLIB_TRIVIAL_CYCLIC_MEMBERS

private:
   typedef typename A::Factorization::const_iterator FacI;
   mutable typename A::Factorization facNMinus1;
};

}  // namespace Private


/**
 *  Factor Field
 *
 *  Based on A::algebra_categoyr and A::polynomial_category, the proper base
 *  class is included.
 */

template<class A>
class FactorField
   : public Private::FactorFieldB<A,typename A::algebra_category>
{
public:
   typedef typename Private::FactorFieldB<A, typename A::algebra_category>
                           ::algebra_category algebra_category;
   typedef typename A::type type;
   typedef nozerodivisor_tag zerodivisor_category;

   FactorField (const A& a, const type& modulus)
      : Private::FactorFieldB<A,typename A::algebra_category> (a, modulus) {}

   type recip (const type& u) const  { type r; recipImp (u,r); return r; }
   void reciprocal (type& u) const  { recipImp (u, u); }

   type div (const type& u, const type& v) const
      { type r; recipImp (v, r); return mul (u, r); }
   void divBy (type& u,  const type& v) const
      { type r; recipImp (v, r); mulBy (u, r); }
};

}  // namespace HIntLib

#endif

