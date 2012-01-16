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

/**
 *  A Polynomial over an arbitrary ring or field.
 */

#ifndef HINTLIB_POLYNOMIAL_H
#define HINTLIB_POLYNOMIAL_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <vector>
#include <iosfwd>

#include <HIntLib/algebra.h>

namespace HIntLib
{

/**
 *  Polynomial
 *
 *  Polynomial with an arbitrary number of terms, coefficients of type T
 */

template <class T>
class Polynomial
{
private:
   typedef Polynomial<T> P;
   typedef std::vector<T> V;

   V c;

public:
   typedef T coeff_type;
   typedef typename V::reference coeff_reference;
   typedef typename V::const_reference coeff_const_reference;

   enum Consts { ZERO, ONE };

   explicit Polynomial (unsigned size = 10) { c.reserve (size); }
   template<typename I> Polynomial (I a0, I an)
      : c (std::reverse_iterator<I> (an), std::reverse_iterator<I> (a0)) {}
   Polynomial (const P& p) : c (p.c) {}

   P& operator= (const P& p);

   int degree() const  { return int (c.size()) - 1; }

   coeff_const_reference operator[] (unsigned i) const { return c[degree()-i]; }
   coeff_reference       operator[] (unsigned i)       { return c[degree()-i]; }
       
   coeff_const_reference lc() const  { return c[0]; }
   coeff_reference       lc()        { return c[0]; }

   // Comparation

   bool operator== (const P& p) const { return c == p.c; }
   bool operator!= (const P& p) const { return ! operator==(p); }

   P& divByX () { if (degree() >= 0)  c.pop_back(); return *this; }
   P& mulByX () { c.push_back(0); return *this; }
   P& divByXPow (unsigned);
   P& mulByXPow (unsigned);
   P& mulAndAdd (const T& x)  { c.push_back (x); return *this; }
};

template<class A>
std::ostream& operator<< (std::ostream &, const Polynomial<A> &);


/**
 *  Polynomial Ring Base Base
 */

template <typename A>
class PolynomialRingBB
{
public:
   typedef A coeff_algebra;
   typedef typename A::type coeff_type;
   typedef Polynomial<coeff_type> type;
   typedef typename type::coeff_reference coeff_reference;

   PolynomialRingBB (const A& _a) : a(_a) {}

   unsigned size() const  { return 0; }
   coeff_algebra getCoeffAlgebra() const  { return a; }
   
   type zero() const  { return type(0);}
   type one()  const;

   bool is0 (const type &p) const  { return p.degree() == -1; }
   bool is1 (const type &)  const;

   type element(unsigned i) const;
   unsigned index (const type &p) const;
   
   type& negate (type&)     const;
   type neg (const type& p) const;

   type  add (const type&, const type&) const;
   type& addTo (type& p1,  const type& p2) const
      { p1 = add (p1, p2); return p1; }

   type  sub (const type& p1, const type& p2) const
      { return add (p1, neg (p2)); }
   type& subFrom (  type& p1, const type& p2) const
      { return addTo (p1, neg (p2)); }

   type mul (const type&, const type&) const;
   type& mulBy (type& p1, const type& p2) const  { return p1 = mul (p1, p2); }

   type times (const type& p, unsigned k) const;
   type power (const type& p, unsigned k) const { return powInt (*this, p, k); }

   bool isMonic (const type&) const;
   coeff_type evaluate (const type&, const coeff_type&) const;

   void print (std::ostream &, const type&) const;
   void printShort (std::ostream &, const type&) const;
   void printSuffix (std::ostream &o) const  { a.printSuffix(o); }

protected:
   const A a;
};

template<class A>
std::ostream& operator<< (std::ostream &, const PolynomialRingBB<A> &);


/**
 *  Polynomial Ring Base
 */

template <typename A, typename TAG> class PolynomialRingBase {};

template <typename A>
class PolynomialRingBase<A,ring_tag> : public PolynomialRingBB<A>
{
public:
   typedef ring_tag algebra_category;
   typedef polynomial_tag polynomial_category;
   PolynomialRingBase (const A& _a) : PolynomialRingBB<A> (_a) {}
};

template <typename A>
class PolynomialRingBase<A,domain_tag> : public PolynomialRingBB<A>
{
public:
   typedef domain_tag algebra_category;
   typedef polynomial_tag polynomial_category;
   PolynomialRingBase (const A& _a) : PolynomialRingBB<A> (_a) {}
};

template <typename A>
class PolynomialRingBase<A,ufd_tag> : public PolynomialRingBB<A>
{
public:
   //typedef ufd_tag algebra_category;
   typedef domain_tag algebra_category;
   typedef polynomial_tag polynomial_category;

   typedef typename A::type coeff_type;
   typedef Polynomial<coeff_type> type;

   PolynomialRingBase (const A& _a) : PolynomialRingBB<A> (_a) {}

   bool isUnit  (const type&p) const  { return p.degree() == 0;}
   type unitRecip (const type& p) const
      { type q (1); q.mulAndAdd (a.recip (p[0])); return q; }

#if 0
   // This could be implemented (sse TACP, 4.6.1)

   bool isPrime (const type&)  const;
   bool isIrreducible (const type &p) const  { return isPrime (p); }
   bool isComposit (const type& p) const
      { return p.degree () > 1 && ! isPrime(p); }
#endif
};

template <typename A>
class PolynomialRingBase<A,field_tag> : public  PolynomialRingBase<A,ufd_tag>
{
public:
   typedef euclidean_tag algebra_category;
   typedef polynomial_tag polynomial_category;

   typedef typename A::type coeff_type;
   typedef Polynomial<coeff_type> type;

   PolynomialRingBase (const A& _a) :  PolynomialRingBase<A,ufd_tag> (_a) {}
   void div (const type&, const type&, type&, type&) const;
   type quot (const type&, const type&) const;
   type rem  (const type&, const type&) const;
   unsigned numOfRemainders (const type&) const;
   unsigned norm (const type& p) const  { return p.degree() + 1; }

   // The following should be in the UFD case...
   bool isPrime (const type&)  const;
   bool isIrreducible (const type &p) const  { return isPrime (p); }
   bool isComposit (const type& p) const
      { return p.degree () > 1 && ! isPrime(p); }
};


/**
 *  Polynomial Ring
 */

template<typename T> struct RingId {};
template<typename T> struct RingId<T*>
   { typedef ring_tag algebra_category; };
template<typename T> struct RingId<T**>
   { typedef domain_tag algebra_category; };
template<typename T> struct RingId<T***>
   { typedef ufd_tag algebra_category; };
template<typename T> struct RingId<T*****>
   { typedef field_tag algebra_category; };

template <typename A>
class PolynomialRing : public PolynomialRingBase<A,
        typename RingId<typename A::algebra_category::MAGIC>::algebra_category>
{
public:
   PolynomialRing (const A& _a)
      : PolynomialRingBase<A,
        typename RingId<typename A::algebra_category::MAGIC>::algebra_category>
        (_a) {}
};

} // namespace HIntLib

#endif

