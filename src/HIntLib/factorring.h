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
 *  Factor B <A,AC,PC>
 */

template<class A, class AC, class PC> class FactorB;

// for polynomials

template<class A>
class FactorB<A,euclidean_tag,polynomial_tag>
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
   FactorB (const A& a, const type& modulus)
      : FactorPolyB<A,typename A::coeff_algebra::size_category> (a, modulus) {}

   void recipImp (const type& x, type& r) const;
};

// for integers

template<class A>
class FactorB<A,integer_tag,nopolynomial_tag> : public FactorBB<A>
{
public:
   typedef typename A::type type;
   typedef finite_tag size_category;

   type element (unsigned i) const  { return int (i); }
   unsigned index (const type& u) const  { return int (u); }

   type  add (const type& u, const type& v) const
      { type t = A::add(u,v); return (t>=m) ? A::sub(t,m) : t; }
   void addTo (    type& u, const type& v) const
      { A::addTo (u,v); if (u >= m) A::subFrom (u,m); }

   type  neg (const type& u) const  { return ! is0(u) ? A::sub(m,u) : type(); }
   void negate    (type& u) const  { if (! is0(u)) u = A::sub(m,u); }

   type  sub (const type& u, const type& v) const
      { type t = A::sub(u,v); return (t<0) ? A::add(t,m) : t; }
   void subFrom   (type& u, const type& v) const
      { A::subFrom (u,v); if (u < 0)  A::addTo (u,m); }

   type dbl (const type& u) const
      { type t = A::dbl(u); return t >= m ? A::sub (t,m) : t; } 
   void times2 (type& u) const
      { A::times2 (u); if (u >= m) A::subFrom (u,m); }
   type times (const type& u, unsigned k) const
      { return A::rem (A::times (u, k), m); }

protected:
   FactorB (const A& a, const type& modulus) : FactorBB<A> (a, modulus) {}

   void recipImp (const type& x, type& r) const;
};


/*************  Factor Ring  *************************************************/


/**
 *  Factor Ring B <A,AC,PC>
 */

template<class A, typename AC, typename PC> class FactorRingB;

// polynomial

template<class A>
class FactorRingB<A,euclidean_tag,polynomial_tag>
 : public FactorB<A,euclidean_tag,polynomial_tag>
{
public:
   unsigned numNilpotents() const;

protected:
   FactorRingB (const A&, const typename A::type&);

   typename A::type nilradical;
};

// integers

template<class A>
class FactorRingB<A,integer_tag,nopolynomial_tag>
 : public FactorB<A,integer_tag,nopolynomial_tag>
{
public:
   typedef char_non char_category;

   unsigned additiveOrder (const typename A::type& u) const;
   unsigned numNilpotents () const
      { return int (A::div (m, nilradical)); }


protected:
   FactorRingB (const A&, const typename A::type&);

   typename A::type nilradical;
};

}  // namespace Private


/**
 *  Factor Ring
 */

template<class A>
class FactorRing
   : public Private::FactorRingB<A,typename A::algebra_category,
                                   typename A::polynomial_category>
{
public:
   typedef typename A::type type;
   typedef type unit_type;
   typedef ring_tag algebra_category;
   typedef zerodivisor_tag zerodivisor_category;
   
   FactorRing (const A& a, const type& modulus)
      : Private::FactorRingB<A,typename A::algebra_category,
                               typename A::polynomial_category> (a, modulus) {}

   bool isUnit (const type& u) const
      { return genIsCoprime (arithmetic(), u, m); }
   bool isZerodivisor (const type& u) const
      { return ! genIsCoprime (arithmetic(), u, m); }
   bool isNilpotent (const type& u) const
      { return A::isDivisor (u, nilradical); }

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
 *  Factor Field B <A,AC,PC>
 */

template<class A, typename AC, typename PC> class FactorFieldB;

// polynomial

template<class A>
class FactorFieldB<A,euclidean_tag,polynomial_tag>
  : public FactorB<A,euclidean_tag,polynomial_tag>
{
public:
   typedef typename A::type type;
   typedef typename A::coeff_algebra::algebra_category::F algebra_category;

   FactorFieldB (const A& a, const type& modulus)
      : FactorB<A,euclidean_tag,polynomial_tag> (a, modulus) {}

   unsigned extensionDegree() const
      { return getCoeffAlgebra().extensionDegree() * m.degree(); }

   type frobenius (const type& u) const  { return power (u, characteristic()); }
   type invFrobenius (const type& u) const
      { return power (u, powInt (characteristic(), extensionDegree() - 1)); }
};

// integers

template<class A>
class FactorFieldB<A,integer_tag,nopolynomial_tag>
  : public FactorB<A,integer_tag,nopolynomial_tag>
{
public:
   typedef typename A::type type;
   typedef cyclic_tag algebra_category;
   typedef char_prime char_category;

   FactorFieldB (const A& a, const type& modulus)
      : FactorB<A,integer_tag,nopolynomial_tag> (a, modulus) {}

   unsigned characteristic() const  { return int (m); }

   HINTLIB_TRIVIAL_CYCLIC_MEMBERS
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
   : public Private::FactorFieldB<A,typename A::algebra_category,
                                    typename A::polynomial_category>
{
public:
   typedef typename Private::FactorFieldB<A,
      typename A::algebra_category,
      typename A::polynomial_category>::algebra_category algebra_category;
   typedef typename A::type type;
   typedef nozerodivisor_tag zerodivisor_category;

   FactorField (const A& a, const type& modulus)
      : Private::FactorFieldB<A,typename A::algebra_category,
                                typename A::polynomial_category> (a, modulus) {}

   type recip (const type& u) const  { type r; recipImp (u,r); return r; }
   void reciprocal (type& u) const  { recipImp (u, u); }

   type div (const type& u, const type& v) const
      { type r; recipImp (v, r); return mul (u, r); }
   void divBy (type& u,  const type& v) const
      { type r; recipImp (v, r); mulBy (u, r); }

   bool isPrimitiveElement (const type& u) const;

private:
   unsigned orderImp (const type&, field_tag) const;
   unsigned orderImp (const type&, gf_tag) const;
public:
   unsigned order (const type& u) const
      { return orderImp (u, algebra_category()); }
};

}  // namespace HIntLib

#endif

