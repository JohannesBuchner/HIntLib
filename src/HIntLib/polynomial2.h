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
 *  This Class models the polynomiols over the field {0,1} = GF_2
 *
 *  An integer type has to be supplied. The bits of this type are used to 
 *  store and manage the coefficients of the polynom, i.e.
 *  unsigned long int can store 32 coeffizients.
 *
 */

#ifndef POLYNOMIAL2_H
#define POLYNOMIAL2_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <iosfwd>

#include <HIntLib/defaults.h>

#ifdef HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/hintlib_limits.h>
#endif

#include <HIntLib/bitop.h>
#include <HIntLib/prime.h>
#include <HIntLib/exception.h>


namespace HIntLib
{

using ::abs;

/**
 *  Polynomial2 Class
 */

template <class T>
class GenPolynomial2
{
private:
   typedef GenPolynomial2<T> P;

   bool msb() const  { return d & (T(1) << std::numeric_limits<T>::digits-1); }
public:

   // Constructors

   GenPolynomial2 ()             : d (T()) {}   // Default is f(x) = 0
   GenPolynomial2 (const P &p)   : d (p.d) {}   // Copy constructor
   explicit GenPolynomial2 (T d) : d (d)   {}   // Conversion from int

   // You can convert a polynomial to an int to get the coefficient bitmap

   operator T () { return d; }

   // Returns a certain coefficient (0 or 1).

   int operator[] (unsigned i) const { return (d >> i) & 1; }

   // Degree of a polynomial. degree(f(x)=0) = -1

   int degree() const { return ms1 (d); }

   // Addition and Substraction

   P operator+ (const P p) const  { return P (d ^ p.d); }
   P operator- (const P p) const  { return P (d ^ p.d); }
   P operator+ () const   { return *this; }   // Unary plus
   P operator- () const   { return *this; }   // Unary minus
   P& operator+= (const P p)  { d ^= p.d; return *this; }
   P& operator-= (const P p)  { d ^= p.d; return *this; }

   // Multiplikation and Division

   static void div (const P u, const P v, P &q, P &r);

   P  operator* (const P p) const;
   P  operator/ (const P p) const { P q, r; div (*this, p, q, r); return q; }
   P  operator% (const P p) const { P q, r; div (*this, p, q, r); return r; }
   P& operator*= (const P p) { return *this = *this * p; }
   P& operator/= (const P p) { P r; div (*this, p, *this, r); return *this; }
   P& operator%= (const P p) { P q; div (*this, p, q, *this); return *this; }

   P& divByX () { d >>= 1; return *this; }
   P& mulByX () { if (msb()) throw Overflow(); d <<= 1; return *this; }
   P& divByXPow (unsigned k)  { d >>= k; return *this; }
   P& mulByXPow (unsigned k)  { d <<= k; return *this; }

   // Comparation

   bool operator== (const P p) const { return p.d == d; }
   bool operator!= (const P p) const { return p.d != d; }

   bool is0() const  { return d == 0; }
   bool is1() const  { return d == 1; }
   bool isX() const  { return d == 2; }

   // Compare degree

   bool operator<  (const P p) const { return degree () <  p.degree (); }
   bool operator<= (const P p) const { return degree () <= p.degree (); }
   bool operator>  (const P p) const { return degree () >  p.degree (); }
   bool operator>= (const P p) const { return degree () >= p.degree (); }

   // Constants
   // FIX ME!!!  For some reason this does not work with GCC :-((

   static P zero () { return P (0); }
   static P one ()  { return P (1); }
   static P x ()    { return P (2); }
   static P xPow (unsigned k) { return P (T(1) << k); }

   bool isPrimitive() const;
   bool isIrreducible() const;

private:
   T d;      // coefficient bitmap
};

template<class T>
std::ostream& operator<< (std::ostream &, const GenPolynomial2<T>);

// Define a type for the standard implementation (using unsigned int)

typedef GenPolynomial2<u32> Polynomial2;


// Define other operations


/**
 *  Multiplication with an integer
 */

template<class T>
inline
GenPolynomial2<T> operator* (const GenPolynomial2<T> p, int k)
{
   return (k & 1) ? p : zero();
}

template<class T>
inline
GenPolynomial2<T> operator* (int k, const GenPolynomial2<T> p)
{
   return (k & 1) ? p : zero();
}


class Polynomial2RingBase
{
public:
   unsigned size() const  { return 0; }
};

template <class T>
class Polynomial2Ring : public Polynomial2RingBase
{
private:
   typedef GenPolynomial2<T> P;

public:
   typedef P type;

   P zero() const  { return P(0); }
   P one() const  { return P(1); }

   bool is0 (P p) const  { return p.is0(); }
   bool is1 (P p) const  { return p.is1(); }

   P element(unsigned i) const { return P(i); }
   unsigned index (P x) const { return T(x); }

   P add (const P& a, const P& b) const  { return a +  b; }
   P& addTo (P& a,    const P& b) const  { return a += b; }

   P neg (const P& a) const  { return     -a; }
   P& negate (P& a) const    { return a = -a; }

   P sub (const P& a, const P& b) const  { return a -  b; }
   P& subFrom (P& a,  const P& b) const  { return a -= b; }

   P mul (const P& a, const P& b) const  { return a *  b; }
   P& mulBy (P& a,    const P& b) const  { return a *= b; }

   void div (const P& a, const P& b, P& q, P& r) const { P::div (a, b, q, r); }

   P times (const P& p, unsigned k) const  { return odd (k) ? p : P(0); }
   P power (const P& p, unsigned k) const  { return powInt (p, k); }

   bool isUnit (const P& p)  const  { return p.is1(); }
   bool isPrime (const P& p) const  { return p.isIrreducible(); }
   P unitRecip (const P& p)  const  { return p; }
};

std::ostream& operator<< (std::ostream &, const Polynomial2RingBase &);

} // namespace HIntLib

#endif

