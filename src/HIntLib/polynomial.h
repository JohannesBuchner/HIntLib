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

#include <HIntLib/defaults.h>
#include <HIntLib/algebra.h>
#include <HIntLib/integerring.h>

namespace HIntLib
{

/**
 *  Polynomial
 *
 *  Polynomial with an arbitrary number of terms, coefficients of type T
 */

template <typename T>
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

   template<typename I>
   void toCoeff (I a0) const  { std::copy (c.rbegin(), c.rend(), a0); }

   // Comparation

   template<typename TT>
   friend
   bool operator== (const Polynomial<TT>&, const Polynomial<TT>&);

   P& divByX () { if (degree() >= 0)  c.pop_back(); return *this; }
   P& mulByX () { if (degree() >= 0)  c.push_back(coeff_type()); return *this; }
   P& divByX (unsigned);
   P& mulByX (unsigned);

   P& mulAndAdd (const T& x)  { c.push_back (x); return *this; }
};

// Comparision

template<typename T>
inline
bool operator== (const Polynomial<T>& p1, const Polynomial<T>& p2)
{
   return p1.c == p2.c;
}

template<typename T>
bool operator== (const Polynomial<Real<T> >&, const Polynomial<Real<T> >&);


/**
 *  Polynomial Ring Base Base
 */

template <class A>
class PolynomialRingBB
{
public:
   typedef A coeff_algebra;
   typedef typename A::type coeff_type;
   typedef Polynomial<coeff_type> type;
   typedef typename type::coeff_reference coeff_reference;

   PolynomialRingBB (const A& _a, char _var) : a(_a), var(_var) {}

   coeff_algebra getCoeffAlgebra() const  { return a; }
   
   unsigned size() const  { return 0; }

   type x (unsigned k = 1) const;
   type one() const  { return x (0); }

   bool is0 (const type &p) const  { return p.degree() == -1; }
   bool is1 (const type &)  const;
   bool isCanonical (const type& p) const;

   type element(unsigned) const;
   type elementMonic (unsigned) const;
   unsigned index (const type&) const;
   
   type& negate (type&)     const;
   type neg (const type& p) const;

   type  add (const type&, const type&) const;
   type& addTo (type& p1,  const type& p2) const  { return p1 = add (p1, p2); }

   type  sub (const type& p1, const type& p2) const
      { return add (p1, neg (p2)); }
   type& subFrom (  type& p1, const type& p2) const
      { return addTo (p1, neg (p2)); }

   type times (const type& p, unsigned k) const;
   type derivative (const type& p) const;

   coeff_type evaluate (const type&, const coeff_type&) const;

   char getVar() const  { return var; }
   void print (std::ostream &, const type&) const;
   void printShort (std::ostream &, const type&) const;
   void printSuffix (std::ostream &o) const  { a.printSuffix(o); }

protected:
   const A a;
   const char var;
};

template<class A>
std::ostream& operator<< (std::ostream &, const PolynomialRingBB<A> &);


/**
 *  Polynomial Ring Base
 */

// Basis case

template <class A, typename TAG> class PolynomialRingBase {};


// ring_tag

template <class A>
class PolynomialRingBase<A,ring_tag> : public PolynomialRingBB<A>
{
public:
   typedef ring_tag algebra_category;
   typedef polynomial_tag polynomial_category;
   typedef typename A::type coeff_type;
   typedef Polynomial<coeff_type> type;

   PolynomialRingBase (const A& _a, char _var)
      : PolynomialRingBB<A> (_a, _var) {}

   type mul (const type&, const type&) const;
   type& mulBy (type& p1, const type& p2) const  { return p1 = mul (p1, p2); }
   type power (const type& p, unsigned k) const { return powInt (*this, p, k); }

   unsigned additiveOrder (const type& a) const;
};


// domain_tag

template <class A>
class PolynomialRingBase<A,domain_tag> : public PolynomialRingBB<A>
{
public:
   typedef domain_tag algebra_category;
   typedef polynomial_tag polynomial_category;
   typedef typename A::type coeff_type;
   typedef Polynomial<coeff_type> type;

   PolynomialRingBase (const A& _a, char _var)
      : PolynomialRingBB<A> (_a, _var) {}

   type mul (const type&, const type&) const;
   type& mulBy (type& p1, const type& p2) const  { return p1 = mul (p1, p2); }
   type power (const type& p, unsigned k) const { return powInt (*this, p, k); }

   unsigned order(const type&) const;
   unsigned characteristic() const  { return a.characteristic(); }

   HINTLIB_TRIVIAL_DOMAIN_MEMBERS
};


// field_tag

template<class A>
class PolynomialRingBase<A,field_tag> : public  PolynomialRingBase<A,domain_tag>
{
public:
   typedef domain_tag algebra_category;
   typedef typename A::type coeff_type;
   typedef coeff_type unit_type;
   typedef Polynomial<coeff_type> type;

   PolynomialRingBase (const A& _a, char _var)
      : PolynomialRingBase<A,domain_tag> (_a, _var) {}

   unit_type unitRecip (const unit_type& u) const  { return a.recip (u); }
   void div (const type&, const type&, type&, type&) const;
   type quot (const type&, const type&) const;
   type rem  (const type&, const type&) const;
   unsigned numOfRemainders (const type&) const;
   unsigned norm (const type& p) const  { return p.degree() + 1; }

   unit_type makeCanonical (type &) const;

   // units

   bool isUnit  (const type&p) const  { return p.degree() == 0;}
   unsigned numUnits() const  { return a.size() ? a.size() - 1 : 0; }
   unit_type unitElement (unsigned i) const  { return a.element(i + 1); }
   unsigned unitIndex (const unit_type& u) const  { return a.index (u) - 1; }

   type    fromUnit (const unit_type& u) const { return type().mulAndAdd(u); }
   unit_type toUnit (const type& p)      const { return p[0]; }

   type  mulUnit   (const type&, const unit_type&) const;
   type& mulByUnit (      type&, const unit_type&) const;

   unit_type  mulUnit   (const unit_type& u1, const unit_type& u2) const
      { return a.mul (u1, u2); }
   unit_type& mulByUnit (      unit_type& u1, const unit_type& u2) const
      { return a.mulBy (u1, u2); }
};


// rational_tag

template<class A>
class PolynomialRingBase<A,rational_tag> : public PolynomialRingBase<A,field_tag>
{
public:
   typedef typename A::type coeff_type;
   typedef Polynomial<coeff_type> type;

   PolynomialRingBase (const A& _a, char _var)
      : PolynomialRingBase<A,field_tag> (_a, _var) {}

   bool isPrime    (const type& p) const;
   bool isComposit (const type& p) const
      { return p.degree() > 1 && ! isPrime (p); }

   class PrimeGenerator;
};

template<class A>
class PolynomialRingBase<A,rational_tag>::PrimeGenerator
{
public:
   PrimeGenerator(const PolynomialRingBase<A,rational_tag> &_alg)
      : alg (_alg), n(2) {}
   type next();
private:
   PrimeGenerator ();
   PrimeGenerator& operator= (const PrimeGenerator&);
   const PolynomialRingBase<A,rational_tag> alg;
   unsigned n;
};


// real_tag

template<class A>
class PolynomialRingBase<A,real_tag> : public PolynomialRingBase<A,field_tag>
{
public:
   typedef typename A::type coeff_type;
   typedef Polynomial<coeff_type> type;

   PolynomialRingBase (const A& _a, char _var)
      : PolynomialRingBase<A,field_tag> (_a, _var) {}

   bool isPrime    (const type& p) const;
   bool isComposit (const type& p) const
      { return p.degree() > 1 && ! isPrime (p); }

   class PrimeGenerator;
};

template<class A>
class PolynomialRingBase<A,real_tag>::PrimeGenerator
{
public:
   PrimeGenerator(const PolynomialRingBase<A,real_tag> &_alg)
      : alg (_alg.getCoeffAlgebra()), n(0) {}
   type next();
private:
   PrimeGenerator ();
   PrimeGenerator& operator= (const PrimeGenerator&);
   const A alg;
   unsigned n;
};


// complex_tag

template<class A>
class PolynomialRingBase<A,complex_tag> : public PolynomialRingBase<A,field_tag>
{
public:
   typedef typename A::type coeff_type;
   typedef Polynomial<coeff_type> type;

   PolynomialRingBase (const A& _a, char _var)
      : PolynomialRingBase<A,field_tag> (_a, _var) {}

   bool isPrime    (const type& p) const  { return p.degree() == 1; };
   bool isComposit (const type& p) const  { return p.degree() >  1; };

   class PrimeGenerator;
};

template<class A>
class PolynomialRingBase<A,complex_tag>::PrimeGenerator
{
public:
   PrimeGenerator(const PolynomialRingBase<A,complex_tag> &_alg)
      : alg (_alg.getCoeffAlgebra()), n(0) {}
   type next();
private:
   PrimeGenerator ();
   PrimeGenerator& operator= (const PrimeGenerator&);
   const A alg;
   unsigned n;
};


// gf_tag

namespace Private
{
   unsigned funnySum (int, unsigned);
}

template<class A>
class PolynomialRingBase<A,gf_tag> : public  PolynomialRingBase<A,field_tag>
{
public:
   typedef typename A::type coeff_type;
   typedef Polynomial<coeff_type> type;

   PolynomialRingBase (const A& _a, char _var)
      : PolynomialRingBase<A,field_tag> (_a, _var) {}

   bool isPrimitive (const type&) const;
   bool isPrime (const type&) const;
   bool isComposit (const type& p) const
      { return p.degree() > 1 && ! isPrime (p); }

   class PrimeGenerator;
};

template<class A>
class PolynomialRingBase<A,gf_tag>::PrimeGenerator
{
public:
   PrimeGenerator(const PolynomialRingBase<A,gf_tag> &_alg)
      : alg (_alg), n(2) {}
   PrimeGenerator(const PolynomialRingBase<A,gf_tag> &_alg, int deg)
      : alg (_alg), n(Private::funnySum (deg, alg.getCoeffAlgebra().size()))
      {}

   type next();
private:
   const PolynomialRingBase<A,gf_tag> alg;
   unsigned n;
   PrimeGenerator ();
   PrimeGenerator& operator= (const PrimeGenerator&);
};


/**
 *  Polynomial Ring
 */

template<typename T> struct RingId {};
template<typename T> struct RingId<T*>    // Ring 
   { typedef ring_tag algebra_category; };
template<typename T> struct RingId<T**>   // Domain
   { typedef domain_tag algebra_category; };
template<typename T> struct RingId<T***>  // Field
   { typedef field_tag algebra_category; };
template<> struct RingId<unsigned***>     // Rational Field
   { typedef rational_tag algebra_category; };
template<> struct RingId<int***>          // Real Field
   { typedef real_tag algebra_category; };
template<> struct RingId<char***>         // Complex Field
   { typedef complex_tag algebra_category; };
template<typename T> struct RingId<T****> // Finite field
   { typedef gf_tag algebra_category; };

template<class A>
class PolynomialRing : public PolynomialRingBase<A,
        typename RingId<typename A::algebra_category::MAGIC>::algebra_category>
{
public:
   PolynomialRing (const A& _a, char _var = 'x')
      : PolynomialRingBase<A,
        typename RingId<typename A::algebra_category::MAGIC>::algebra_category>
        (_a, _var) {}
};

} // namespace HIntLib

#endif

