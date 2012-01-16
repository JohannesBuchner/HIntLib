/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schürer <rudolf.schuerer@sbg.ac.at>
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
#include <algorithm>

#include <HIntLib/algebra.h>

namespace HIntLib
{


/************  Quotient <>  ****************************************************/


/**
 *  Quotient<>
 *
 *  A concrete object storing numerator and denominator of a certain type T.
 *
 *  A value of 0 is always denoted by 0/0 to allow default construction.
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

} // namespace HIntLib

namespace std
{
   template<typename T>
   inline
   void swap (::HIntLib::Quotient<T>& q1, ::HIntLib::Quotient<T> q2)
   {
      swap (q1.num, q2.num);
      swap (q1.den, q2.den);
   }
}  // namespace std

namespace HIntLib
{

/**
 *  destructiveAssign()
 *
 *  The most efficient way to destructively assign a quotient is probably a
 *  destructive assignment of numerator and denominator.
 */

template<typename T>
inline
void destructiveAssign (Quotient<T>& dst, Quotient<T>& src)
{
   destructiveAssign (dst.num, src.num);
   destructiveAssign (dst.den, src.den);
}


/**
 *  operator== ()
 */

template<typename T>
bool operator== (const Quotient<T>& u, const Quotient<T>& v)
{
   return (u.num == v.num && u.den == v.den);
}


/*********************  Quotient Field <>  *************************************/


namespace Private
{

/**
 * Quotient Field Base
 */

template<typename A>
class QFB
{
public:
   typedef typename A::char_category char_category;
   typedef A base_algebra;
   typedef typename A::type base_type;
   typedef Quotient<base_type> type; 

   const base_algebra& getBaseAlgebra () const  { return a; }

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

   void negate (type& u) const  { a.negate (u.num); }
   type neg (const type& u) const  { return type (a.neg (u.num), u.den); }

   type add (const type&, const type&) const;
   type sub (const type&, const type&) const;

   type mul (const type&, const type&) const;
   type div (const type&, const type&) const;

   type recip (const type&) const;
   void reciprocal (type&) const;

   type sqr (const type& u) const { return type (a.sqr(u.num), a.sqr(u.den)); }
   void square (type& u) const    { a.square (u.num); a.square (u.den); }
   type power (const type& u, unsigned k) const
      { return type (a.power (u.num, k), a.power (u.den,k)); }

   void addTo   (type& u, const type& v) const  { u = add (u, v); }
   void subFrom (type& u, const type& v) const  { u = sub (u, v); }
   void divBy   (type& u, const type& v) const  { u = div (u, v); }
   void mulBy   (type& u, const type& v) const  { u = mul (u, v); }

   void print (std::ostream &, const type&) const;
   void printShort (
         std::ostream &, const type&, PrintShortFlag = PrintShortFlag()) const;
   void printSuffix (std::ostream &o) const  { a.printSuffix(o); }

   HINTLIB_TRIVIAL_DOMAIN_MEMBERS

protected:
   QFB (const A& _a) : a (_a) {}

   type toLowestTerms (const base_type&, const base_type&) const;
   type toLowestTermsAndNormalize (const base_type&, const base_type&) const;

   A a;
};


/**
 * Quotient Field Base 2
 *
 * Some methods (doubeling and times) are can be implemented more efficiently
 * if we know if we have quotients of integers or quotients of polynomials.
 */

template<typename A, typename AC> class QFB2;

template<class A, typename AC>
std::ostream& operator<< (std::ostream &, const QFB2<A,AC> &);

// integers

template<typename A>
class QFB2<A,integer_tag> : public QFB<A>
{
protected:
   QFB2 (const A& _a) : QFB<A> (_a) {}

public:
   typedef typename QFB<A>::type type;

   type dbl (const type&) const;
   void times2 (type&) const;
   type times (const type&, unsigned) const;
};

template<class A>
std::ostream&
operator<< (std::ostream&, const QFB2<A,integer_tag>&);

// polynomials over a field

template<typename A>
class QFB2<A,polyoverfield_tag> : public QFB<A>
{
protected:
   QFB2 (const A& _a) : QFB<A> (_a) {}

public:
   typedef typename QFB<A>::type type;

private:
   type timesImp (const type& u, unsigned, char_prime) const;
   type timesImp (const type& u, unsigned k, char_two) const
      { return (k & 1) ? u : type(); }
   type timesImp (const type& u, unsigned k, char_zero) const
      { return k ? type (this->a.times (u.num, k), u.den) : type(); }

public:
   type times (const type& u, unsigned k) const
      { return timesImp(u,k, typename QFB<A>::char_category()); }

   // for char_two and char_zero, we assume that the compiler is smart enough
   //    to get rid of the condition

   type dbl (const type& u) const
   {
      return this->characteristic() == 2 ?
         type() : type (this->a.dbl (u.num), u.den);
   }
   void times2 (type& u)   const
   {
      if (this->characteristic() == 2) u = type();
      else this->a.times2 (u.num);
   }
};

template<class A>
std::ostream&
operator<< (std::ostream&, const QFB2<A,polyoverfield_tag>&);


/**
 *  Quotient Field Field Id
 *
 *  A template determining the algebra category of the quotient field based on
 *  the algebra category of the base field.
 */

template<typename AC>
struct QuotientFieldFieldId { typedef field_tag algebra_category; };

template<>
struct QuotientFieldFieldId<integer_tag>
{ typedef rational_tag algebra_category; };

template<>
struct QuotientFieldFieldId<polyoverfield_tag>
{ typedef ratfunfield_tag algebra_category; };

} // namespace Private


/**
 *  QuotientField<>
 *
 *  The field of rational numbers over a given Euclidean ring
 */

template<class A>
class QuotientField 
   : public Private::QFB2<A,typename A::algebra_category>
{
public:
   typedef typename Private::QuotientFieldFieldId<typename A::algebra_category>
      ::algebra_category algebra_category;

   typedef nopolynomial_tag polynomial_category;
   typedef infinite_tag size_category;
   typedef nozerodivisor_tag zerodivisor_category;

   QuotientField (const A& _a)
      : Private::QFB2<A,typename A::algebra_category> (_a) {}
};

} // namespace HIntLib

#endif

