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

#ifdef __GNUG__
#pragma implementation
#endif

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/realfield.h>

#ifdef HINTLIB_HAVE_OSTREAM
  #include <ostream>
#else
  #include <iostream>
#endif

#include <HIntLib/polynomial.h>
#include <HIntLib/output.h>

namespace L = HIntLib;
namespace P = HIntLib::Private;


/*****************   Real Field   ********************************************/


/**
 *  operator<<
 */

std::ostream&
P::operator<< (std::ostream &o, const RRing &)
{
   return o << "R";
}


/**
 *  element()
 */

template<typename T>
typename L::RealField<T>::type
L::RealField<T>::element (unsigned i)
{
   if (! i)  return type();
   --i;

   T sign = (i & 1) ? -1 : 1;

   unsigned i1, i2;
   unthread (i / 2 + 1, i1, i2);

   return sign * (T(i1) + radicalInverseFunction2 (i2));
}


/**
 *  index()
 */

template<typename T>
unsigned
L::RealField<T>::index (type r)
{
   if (is0 (r)) return 0;

#if 0
   real integerPart;
   real rest = HINTLIB_MN modf (abs (r), &integerPart);
              // we need an extra real() because abs() is broken in old glibc
#else
   const real absr = abs (r.x);
   const real integerPart = HINTLIB_MN floor (absr);
   real rest = absr - integerPart;
#endif

   unsigned i2 = 0;
   unsigned mask = 1;

   while (rest > 1e-7 && mask < (1 << 16))
   {
      if (rest >= .5 - 1e-7)
      {
         i2 |= mask;
         rest -= .5;
      }
      mask <<= 1;
      rest *= 2.0;
   }

   return (thread (unsigned (integerPart), i2) - 1) * 2 + ((r.x < 0) ? 2 : 1);
}


/**
 *  order()
 */

template<typename T>
unsigned
L::RealField<T>::order (const type& r)
{
   if (is1(r))  return 1;
   if (is1(neg(r)))  return 2;
   return 0;
}


/**
 *  operator==(Polynomial<>, Polynomial<>)
 */

template<typename T>
bool
L::operator== (const Polynomial<Real<T> >& p1, const Polynomial<Real<T> >& p2)
{
   // determine magnitude of coefficients

   int deg1 = p1.degree();
   int deg2 = p2.degree();
   int deg = std::max (deg1, deg2);

   T mag = 0.0;

   for (int i = 0; i <= deg; ++i)
   {
      if (i <= deg1)  mag = std::max (mag, abs (T(p1[i])));
      if (i <= deg2)  mag = std::max (mag, abs (T(p2[i])));
   }

   mag *= T(deg + 1) * std::numeric_limits<T>::epsilon() * 100.0;

   for (int i = 0; i <= deg; ++i)
   {
      T x1 = (i <= deg1) ? T(p1[i]) : T(0);
      T x2 = (i <= deg2) ? T(p2[i]) : T(0);
      if (abs (x1 - x2) > mag)  return false;
   }

   return true;
}


/****************   Complex Field   ******************************************/


/**
 *  operator<<
 */

std::ostream&
P::operator<< (std::ostream &o, const CRing&)
{
   return o << "C";
}


/**
 *  element()
 */

template<typename T>
typename L::ComplexField<T>::type
L::ComplexField<T>::element (unsigned i)
{
   unsigned i1, i2;
   unthread (i, i1, i2);

   RealField<T> f;

   return type (T(f.element(i1)), T(f.element(i2)));
}


/**
 *  index()
 */

template<typename T>
unsigned
L::ComplexField<T>::index (const type& r)
{
   RealField<T> f;

   const unsigned i1 = f.index (r.x.real());
   const unsigned i2 = f.index (r.x.imag());

   return thread (i1, i2);
}


/**
 *  approxc()
 */

template<class T>
bool
L::approxc (const std::complex<T>& a, const std::complex<T>& b, T factor)
{
   return std::abs(a - b)
       <= factor * std::numeric_limits<T>::epsilon()
                 * (std::abs(a) + std::abs(b));
}


#ifdef HINTLIB_HAVE_LONG_DOUBLE
#ifdef HINTLIB_COMPLEX_POW_BUG
/**
 *  power()
 */

template<>
L::ComplexField<long double>::type
L::ComplexField<long double>::power (const type& a, unsigned k)
{
   std::complex<long double> x = a.x;
   std::complex<long double> y = std::complex<long double>(1);

   while (k)
   {
      if (k & 1)  y *= x;
      x *= x;
      k >>= 1;
   }

   return y;
}
#endif
#endif


/**
 *  order()
 */

template<typename T>
unsigned
L::ComplexField<T>::order (const type& r)
{
   RealField<T> f;

   if (! f.is1 (std::abs(r.x)))  return 0;

   T phi = std::arg(r.x) / (T(2) * Constants<T>::pi());
   phi = phi - HINTLIB_MN floor (phi);

   // XXX no we should use Euclides Algorithm and determine the continued
   // fraction expansion of phi. If it is a rational number we have a finite
   // order.

   for (int i = 1; i < 100; ++i)
   {
      T diff = phi * T(i) - HINTLIB_MN floor (phi * T(i));
      if (diff > T(.5)) diff = 1 - diff;
      if (f.is0 (diff))  return i;
   }

   return 0;
}


/**
 *  print()
 */

template<typename T>
void
L::ComplexField<T>::print (std::ostream& o, const type& x)
{
   RealField<T> f;

   if (f.is0 (x.x.imag()))
   {
      o << x.x.real();
   }
   else if (f.is0 (x.x.real()))
   {
      if (f.is1 (x.x.imag()))
      {
         if (o.flags() & o.showpos)  o << "+i";
         else                        o << "i";
      }
      else if (f.is1 (f.neg (x.x.imag())))
      {
         o << "-i";
      }
      else
      {
         Private::Printer p (o);
         p << x.x.imag() << 'i';
      }
   }
   else
   {
      if (o.flags() & o.showpos)
      {
         Private::Printer p (o);
         p.unsetf (o.showpos);
         if (x.x.real() < 0 && x.x.imag() < 0)
         {
            p << '-' << - x.x;
         }
         else
         {
            p << '+' << x.x;
         }
      }
      else
      {
         o << x.x;
      }
   }
}


/**
 *  operator==(Polynomial<>, Polynomial<>)
 */

template<typename T>
bool
L::operator==(const Polynomial<Complex<T> >& pp1,
              const Polynomial<Complex<T> >& pp2)
{
   typedef std::complex<T> TT;

   // Make sure that  deg p1 <= deg p2

   bool swap = pp1.degree() > pp2.degree();

   const Polynomial<Complex<T> >& p1 = swap ? pp2 : pp1;
   const Polynomial<Complex<T> >& p2 = swap ? pp1 : pp2;

   // determine magnitude of coefficients

   int deg1 = p1.degree();
   int deg2 = p2.degree();

   T mag = 0.0;

   for (int i = 0; i <= deg1; ++i)
   {
      mag = std::max (mag, std::abs (p1[i].data()));
      mag = std::max (mag, std::abs (p2[i].data()));
   }
   for (int i = deg1 + 1; i <= deg2; ++i)
   {
      mag = std::max (mag, std::abs (p2[i].data()));
   }

   // Check distances

   mag *= T(deg2 + 1) * std::numeric_limits<T>::epsilon() * 100.0;

   for (int i = 0; i <= deg1; ++i)
   {
      if (std::abs (TT(p1[i].data()) - TT(p2[i].data())) > mag)  return false;
   }
   for (int i = deg1 + 1; i <= deg2; ++i)
   {
      if (std::abs (TT(p2[i].data())) > mag)  return false;
   }

   return true;
}


namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template RealField<X >::type RealField<X >::element(unsigned); \
   template unsigned RealField<X >::index(type); \
   template unsigned RealField<X >::order(const type&); \
   template bool operator== \
      (const Polynomial<Real<X > >&, const Polynomial<Real<X > >&);

   HINTLIB_INSTANTIATE(real)
#undef HINTLIB_INSTANTIATE

#define HINTLIB_INSTANTIATE(X) \
   template ComplexField<X >::type ComplexField<X >::element(unsigned); \
   template unsigned ComplexField<X >::index(const type&); \
   template unsigned ComplexField<X >::order(const type&); \
   template void ComplexField<X >::print(std::ostream&, const type&); \
   template bool approxc (const std::complex<X >&, const std::complex<X >&, X);\
   template bool operator== \
      (const Polynomial<Complex<X > >&, const Polynomial<Complex<X > >&);

   HINTLIB_INSTANTIATE(real)
#undef HINTLIB_INSTANTIATE
}

