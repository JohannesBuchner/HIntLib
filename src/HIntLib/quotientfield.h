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

#ifndef HINTLIB_RATIONAL_FUNCTION_FIELD_H
#define HINTLIB_RATIONAL_FUNCTION_FIELD_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <iosfwd>

#include <HIntLib/algebra.h>

namespace HIntLib
{

/**
 *  Quotient<>
 *
 *  A concrete object storing numerator and denominator of a certain type T.
 */

template<typename T>
struct Quotient
{
   Quotient () : num(), den() {}
   Quotient (const T& a, const T& b) : num(a), den(b) {}

   T numerator () const    { return num; }
   T denominator () const  { return den; }

   // Default copy constructor and assignment is alright

   T num;
   T den;
};

template<typename T>
bool operator== (const Quotient<T>& u, const Quotient<T>& v)
{
   return (u.num == v.num && u.den == v.den);
}


/**
 *  QuotientField<>
 *
 *  The field of rational numbers over a given Euclidean ring
 */

namespace Private
{
   template<typename T>
   struct QuotientFieldFieldId { typedef field_tag algebra_category; };

   template<>
   struct QuotientFieldFieldId<integer_tag>
   { typedef rational_tag algebra_category; };
}

template<class A>
class QuotientField
{
public:
   typedef typename Private::QuotientFieldFieldId<typename A::algebra_category>
      ::algebra_category algebra_category;

   typedef nopolynomial_tag polynomial_category;
   typedef A base_algebra;
   typedef typename A::type base_type;
   typedef Quotient<base_type> type; 

   QuotientField (const A& _a) : a(_a) {}

   base_algebra getBaseAlgebra () const  { return a; }

   static unsigned size()  { return 0; }
   unsigned characteristic() const  { return a.characteristic(); }

   type element (unsigned) const;
   unsigned index (const type &) const;

   type one () const  { return type (a.one(), a.one()); }

   type makeElement (const base_type&) const;
   type makeElement (const base_type& num, const base_type& den) const
      { return toLowestTermsAndNormalize (num, den); }

   bool is0 (const type& u) const  { return a.is0 (u.num); }
   bool is1 (const type& u) const  { return a.is1 (u.num) && a.is1 (u.den); }

   unsigned order (const type&) const;

   type& negate (type& u) const  { a.negate (u.num); return u; }
   type  neg (const type& u) const  { return type (a.neg (u.num), u.den); }

   type add (const type&, const type&) const;
   type sub (const type&, const type&) const;

   type mul (const type& u, const type& v) const
      { return toLowestTerms (a.mul(u.num, v.num), a.mul(u.den, v.den)); }
   type div   (const type&, const type&) const;
   type recip (const type&) const;

   type& addTo   (type& u, const type& v) const { return u = add (u, v); }
   type& subFrom (type& u, const type& v) const { return u = sub (u, v); }
   type& divBy   (type& u, const type& v) const { return u = div (u, v); }
   type& mulBy   (type& u, const type& v) const { return u = mul (u, v); }

   type times (const type& u, unsigned k) const
      { return toLowestTerms (a.times (u.num, k), u.den); }
   type power (const type& u, unsigned k) const
      { return type (a.power (u.num, k), a.power (u.den,k)); }

   void print (std::ostream &, const type&) const;
   void printShort (std::ostream &, const type&) const;
   void printSuffix (std::ostream &o) const  { a.printSuffix(o); }

   HINTLIB_TRIVIAL_DOMAIN_MEMBERS

private:
   type toLowestTerms (const base_type&, const base_type&) const;
   type toLowestTermsAndNormalize (const base_type&, const base_type&) const;

   A a;
};

template<class A>
std::ostream& operator<< (std::ostream &, const QuotientField<A> &);

} // namespace HIntLib

#endif

