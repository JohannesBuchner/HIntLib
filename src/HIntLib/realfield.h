/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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

#ifndef HINTLIB_REAL_FIELD_H
#define HINTLIB_REAL_FIELD_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <complex>
#include <iosfwd>

#include <HIntLib/hlmath.h>
#include <HIntLib/algebra.h>

namespace HIntLib
{
namespace Private
{
   /**
    *  Non-template base classes for printing
    */

   class RRing { protected: RRing() {} };
   class CRing { protected: CRing() {} };

   std:: ostream& operator<< (std:: ostream&, const RRing&);
   std:: ostream& operator<< (std:: ostream&, const CRing&);
#ifdef HINTLIB_BUILD_WCHAR
   std::wostream& operator<< (std::wostream&, const RRing&);
   std::wostream& operator<< (std::wostream&, const CRing&);
#endif

   /**
    *  Number Field
    */

   template<typename T>
   class NumberField
   {
   public:
      // types

      typedef T type;

      // Algebra traits

      typedef infinite_tag size_category;
      typedef nopolynomial_tag polynomial_category;
      typedef char_zero char_category;
      typedef nozerodivisor_tag zerodivisor_category;

      // Properties of Z

      static unsigned size()  { return 0; }
      static unsigned characteristic()  { return 0; }

      // Additive arithmetic

      static T add (const T& a, const T& b)  { return T (a.x +  b.x); }
      static void addTo (T& a,  const T& b)  { a.x += b.x; }

      static T sub (const T& a, const T& b)  { return T (a.x - b.x); }
      static void subFrom (T& a,const T& b)  { a.x -= b.x; }

      static T neg (const T& a)  { return T(-a.x); }
      static void negate (T& a)  { a.x = - a.x; }

      static T  dbl (const T& a)  { return T (a.x + a.x); }
      static void times2 (T& a)   { a.x += a.x; }

      // Multiplicative arithmetic

      static T mul (const T& a, const T& b)  { return T (a.x * b.x); }
      static void mulBy (T& a,  const T& b)  { a.x *= b.x; }

      static T div (const T& a, const T& b)  { return T (a.x / b.x); }
      static void divBy (T& a,  const T& b)  { a.x /= b.x; }

      static T  sqr (const T& a)  { return T (a.x * a.x); }
      static void square (T& a)   { a.x *= a.x; }

      // I/O

      static void printSuffix (std::ostream&)  {}
#ifdef HINTLIB_BUILD_WCHAR
      static void printSuffix (std::wostream&)  {}
#endif
   };

}  // namespace Private

// Forward declarations

template<typename T> class RealField;
template<typename T> class ComplexField;

/**
 *  Real
 */

template<typename T>
class Real
{
public:
   Real () : x(0) {}
   Real (T a) : x(a) {}
   // Real (const Real& a) : x(a.x) {}   // Default copy constructor

   operator T () const  { return x; }

   friend class Private::NumberField<Real<T> >;
   friend class RealField<T>;

   const T& data() const  { return x; }

private:
   T x;
};

template<typename T>
inline
bool operator== (const Real<T>& a, const Real<T>& b)
{
   return approx (a.data(), b.data());
}

template<typename T>
inline
bool operator< (const Real<T>& a, const Real<T>& b)
{
   return a.data() < b.data() && ! (a == b);
}

template<typename T>
inline
Real<T> abs (const Real<T>& a)
{
   return Real<T> (abs (a.data()));
}


/**
 *  Complex
 */

template<typename T>
class Complex
{
public:
   Complex () : x(0) {}
   Complex (const T& a) : x(a) {}
   Complex (const T& a, const T& b) : x(a, b) {}
   Complex (const Real<T>& a) : x (a.data()) {}
   Complex (const Real<T>& a, const Real<T>& b) : x (a.data(), b.data()) {}
   Complex (const std::complex<T>& a) : x(a) {}
   // Complex (const Complex<T>& a) : x(a.x) {}  // Default copy constructor

   operator std::complex<T> () const  { return x; }

   friend class Private::NumberField<Complex<T> >;
   friend class ComplexField<T>;

   const std::complex<T>& data() const  { return x; }
private:
   std::complex<T> x;
};


template<class T>
bool
approxc (const std::complex<T>&, const std::complex<T>&, T factor = 10.0);

template<typename T>
inline
bool operator== (const Complex<T>& a, const Complex<T>& b)
{
   return approxc (a.data(), b.data());
}

template<typename T>
inline
Real<T> abs (const Complex<T>& a)
{
   return Real<T> (abs (a.data()));
}


/**
 *  Specialization of   operator==(Polynomial<>, Polynomial<>)
 */

template<typename T> class Polynomial;

template<typename T>
bool operator==(const Polynomial<Real<T> >&, const Polynomial<Real<T> >&);

template<typename T>
bool operator==(const Polynomial<Complex<T> >&, const Polynomial<Complex<T> >&);


/**
 *  Real Field
 */

template<typename T = real>
class RealField : public Private::NumberField<Real<T> >, public Private::RRing
{
public:
   typedef Real<T> type;
   typedef real_tag algebra_category;

   typedef ComplexField<T> complex_field;
   typedef Complex<T> complex_type;

   static complex_field getComplexField();

   static type one()  { return T(1); }
   static type element(unsigned);
   static unsigned index (type r);

   static bool is0 (const type& a)
      { return abs(a.x) < std::numeric_limits<T>::epsilon() * 100; }
   static bool is1 (const type& a)  { return approx (a.x, T(1)); }

   static type times (const type& a, unsigned k)  { return type (a.x * T(k)); }

   static type recip (const type& a)  { return type (T(1) / a.x); }
   static void reciprocal (type& a)  { a.x = T(1) / a.x; }

   static type power (const type& a, unsigned k)
      { return HINTLIB_MN pow (a.x, int(k)); }

   static unsigned order (const type&);

   // I/O

   static void print (std::ostream &o, const type& a)  { o << a.x; }
   static void printShort (std::ostream &o, const type& a)  { print (o, a); }
   static void printShort (std::ostream &o, const type& a, PrintShortFlag)
      { print (o, a); }
#ifdef HINTLIB_BUILD_WCHAR
   static void print (std::wostream &o, const type& a)  { o << a.x; }
   static void printShort (std::wostream &o, const type& a)  { print (o, a); }
   static void printShort (std::wostream &o, const type& a, PrintShortFlag)
      { print (o, a); }
#endif

   HINTLIB_TRIVIAL_DOMAIN_MEMBERS
};


/**
 *  Complex Field
 */

template<typename T = real>
class ComplexField
   : public Private::NumberField<Complex<T> >, public Private::CRing
{
public:
   typedef Complex<T> type; // Complex<T>
   typedef complex_tag algebra_category;

private:
   typedef std::complex<T> TT;

public:
   typedef RealField<T> real_field;
   typedef Real<T> real_type;

   static real_field getRealField();

   static type one()  { return TT(T(1)); }
   static type element(unsigned);
   static unsigned index (const type& r);

   static bool is0 (const type& a)
      { return abs(a.x) < std::numeric_limits<T>::epsilon() * 100; }
   static bool is1 (const type& a)  { return approxc (a.x, TT(T(1))); }

   static type times (const type& a, unsigned k)
      { return type (a.x * TT(T(k))); }

   static type recip (const type& a)  { return type (TT(T(1)) / a.x); }
   static void reciprocal (type& a)  { a.x = TT(T(1)) / a.x; }

   static type power (const type&, unsigned);

   static unsigned order (const type&);

   static real_type re  (const type& a)  { return a.x.real(); }
   static real_type im  (const type& a)  { return a.x.imag(); }

   // I/O

   static void print (std::ostream &, const type&);
   static void printShort (std::ostream &o, const type& a)  { print (o, a); }
   static void printShort (std::ostream &o, const type& a, PrintShortFlag)
      { print (o, a); }
#ifdef HINTLIB_BUILD_WCHAR
   static void print (std::wostream &, const type&);
   static void printShort (std::wostream &o, const type& a)  { print (o, a); }
   static void printShort (std::wostream &o, const type& a, PrintShortFlag)
      { print (o, a); }
#endif

   HINTLIB_TRIVIAL_DOMAIN_MEMBERS
};

template<typename T>
inline
typename ComplexField<T>::type
ComplexField<T>::power (const type& a, unsigned k)
{
   return std::pow (a.x, int(k));
}

#ifdef HINTLIB_HAVE_LONG_DOUBLE
#ifdef HINTLIB_COMPLEX_POW_BUG

/**
 *  If pow(complex<long double>,int)  is broken, we need to provide a
 *  work around.
 */

template<>
ComplexField<long double>::type
ComplexField<long double>::power (const type&, unsigned);

#endif
#endif

/**
 *  getComplexField()
 *  getRealField()
 *
 *  These inline members are defined down here because we need both class
 *  definitions.
 */

template<typename T>
inline
typename ComplexField<T>::real_field
ComplexField<T>::getRealField()
{
   return real_field();
}

template<typename T>
inline
typename RealField<T>::complex_field
RealField<T>::getComplexField()
{
   return complex_field();
}


}  // namespace HIntLib

#endif

