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

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <vector>
#include <iosfwd>

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

   std::vector<T> c;

public:
   enum Consts { ZERO, ONE };

   Polynomial (unsigned size, const T*);
   Polynomial (unsigned size = 10) { c.reserve (size); }
   Polynomial (const P& p) : c (p.c) {}

   P& operator= (const P& p);

   int degree() const { return int (c.size()) - 1; }

   T  operator[] (unsigned i) const  { return c[degree() - i]; }
   T& operator[] (unsigned i)        { return c[degree() - i]; }
       
         T& lc()        { return c[0]; }
   const T& lc() const  { return c[0]; }

   // Comparation

   bool operator== (const P& p) const { return c == p.c; }
   bool operator!= (const P& p) const { return ! operator==(p); }

   bool is0() const  { return degree() == -1; }

   // Compare degree

   bool operator<  (const P& p) const { return degree() <  p.degree(); }
   bool operator<= (const P& p) const { return degree() <= p.degree(); }
   bool operator>  (const P& p) const { return degree() >  p.degree(); }
   bool operator>= (const P& p) const { return degree() >= p.degree(); }

   P& divByX () { if (degree() >= 0)  c.pop_back(); return *this; }
   P& mulByX () { c.push_back(0); return *this; }
   P& divByXPow (unsigned);
   P& mulByXPow (unsigned);
   P& mulAndAdd (const T& x)  { c.push_back (x); return *this; }
};

template<class A>
std::ostream& operator<< (std::ostream &, const Polynomial<A> &);


/**
 *  Polynomial Ring
 */

template <class A>
class PolynomialRing
{
private:
   typedef typename A::type T;
   typedef Polynomial<T> P;

protected:
   const A a;

public:
   typedef Polynomial<T> type;
   
   PolynomialRing (A _a) : a(_a) {}

   A arithmetic() const  { return a; }
   
   P zero() const  { return P(0);}
   P one()  const;

   bool is0 (const P &p) const  { return p.degree() == -1; }
   bool is1 (const P &)  const;

   P element(unsigned i) const;
   unsigned index (const P &p) const;
   
   P& negate (P&)     const;
   P neg (const P& p) const;

   P  add (const P&, const P&) const;
   P& addTo (P& p1,  const P& p2) const  { p1 = add (p1, p2); return p1; }

   P  sub (const P& p1, const P& p2) const  { return add (p1, neg (p2)); }
   P& subFrom (  P& p1, const P& p2) const  { return addTo (p1, neg (p2)); }

   P mul (const P&, const P&) const;
   P& mulBy (P& p1, const P& p2) const  { return p1 = mul (p1, p2); }

   P times (const P& p, unsigned k) const;
   P power (const P& p, unsigned k) const { return powInt (*this, p, k); }

   bool isMonic (const P&) const;

   unsigned size() const  { return 0; }
};

template<class A>
std::ostream& operator<< (std::ostream &, const PolynomialRing<A> &);


/**
 *  Polynomials over a Field
 */

template <class A>
class PolynomialRingField : public PolynomialRing<A>
{
private:
   typedef typename A::type T;
   typedef Polynomial<T> P;
public:
   PolynomialRingField (A _a) : PolynomialRing<A> (_a) {}
   void div (const P&, const P&, P&, P&) const;
   unsigned numOfRemainders (const P&) const;

   bool isUnit  (const P&p) const  { return p.degree() == 0;}
   bool isPrime (const P&)  const;
   bool isIrreducible (const P &p) const  { return isPrime (p); }

   P unitRecip (const P& p) const
      { P q (1); q.mulAndAdd (a.recip (p[0])); return q; }
};

} // namespace HIntLib

#endif

