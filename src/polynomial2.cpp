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

#ifdef __GNUG__
#pragma implementation
#endif

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_OSTREM
  #include <ostream>
#else
  #include <iostream>
#endif

#ifdef HINTLIB_HAVE_SSTREAM
  #include <sstream>
#else
  #include <HIntLib/fallback_sstream.h>
#endif

#include <HIntLib/polynomial2.h>

#include <HIntLib/integerring.h>

namespace L = HIntLib;

/******************************  GF2  ****************************************/


std::ostream &
L::operator<< (std::ostream &o, const GF2 &)
{
   return o << "GF2";
}

void
L::GF2::print (std::ostream& o, type x)
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   ss << unsigned (x) << " (2)";

   o << ss.str().c_str();
}

void
L::GF2::printShort (std::ostream& o, type x)
{
   o << unsigned (x);
}

void
L::GF2::printSuffix (std::ostream& o)
{
   o << "(2)";
}


/******************************  GF2 Vector Space  ***************************/

std::ostream &
L::operator<< (std::ostream &o, const GF2VectorSpaceBase &v)
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   ss << "GF2^" << v.dimension();

   return o << ss.str().c_str();
}

template<typename T>
void
L::GF2VectorSpace<T>::printShort (std::ostream& o, const type& v) const
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   ss << '(';
   for (unsigned i = 0; i < dim; ++i)
   {
      if (i > 0)  ss << ',';
      ss << unsigned (coord(v, i));
   }
   ss << ')';

   o << ss.str().c_str();
}

template<typename T>
void
L::GF2VectorSpace<T>::print (std::ostream& o, const type& v) const
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   printShort (o, v);
   ss << " (2)";

   o << ss.str().c_str();
}


/*************************  Polynomial2  *************************************/


/**
 *  Prints a polynomial to an output stream.
 */

template<typename T>
void L::Polynomial2<T>::printShort (std::ostream& o, char var) const
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   bool output = false;

   for (int i = std::numeric_limits<T>::digits - 1; i >= 0; --i)
   {
      if ((*this)[i])
      {
         if (output) ss << '+';

         switch (i)
         {
         case 0:  ss << '1'; break;
         case 1:  ss << var; break;
         case 2:  ss << var << '\262'; break;
         case 3:  ss << var << '\263'; break;
         default: ss << var << '^' << i;
         }

         output = true;
      }
   }

   if (! output)  ss << '0';

   o << ss.str().c_str();
}

/**
 *  Multiply two polynomials
 */

template<typename T>
L::Polynomial2<T> L::Polynomial2<T>::operator* (const P p) const
{
   // Initialize result to f(x)=0

   P result (0);

   // Assign bitmap of higher degree Polynomial to da, the lower one to db

   P a, b;

   if (d > p.d)
   {
      a = *this;
      b = p;
   } else {
      a = p;
      b = *this;
   }

   // Loop to handle every 1 bit in bd

   for (;;)
   {
      if (b[0])  result += a;   // ADD da to result, if lsb in db is 1

      b.divByX ();

      if (b.is0())  break;   // Any more 1s in db? NO -> done!

      // Check for overflow

      a.mulByX();
   }

   return result;
}


/**
 *  div
 *
 *  Division of polynomials
 *
 *  See D.E.Knuth, Art o Computer Programming, 4.6.1, Algo D
 */

template<typename T>
void L::Polynomial2<T>::div (P u, const P v, P &q, P &r)
{
   const int m = u.degree();
   const int n = v.degree();

   if (n == -1)  throwDivisionByZero();

   int loops = m - n;

   if (loops < 0)
   {
      r = u;
      q = P();

      return;
   }

   // Kill highest order coefficient of divisor

   P divisor (v + xPow (n));

   // align with dividend

   divisor.mulByX (loops);

   T u_mask = T(1) << m;

   // Cray assumes that u is const.  So we have to make another copy

   // SGI and CRAY seem to require a copy...   FIXME
   P uu (u);
   #define u uu

   for ( ; loops >= 0; --loops)
   {
      // If the current highest coefficient of u is set, v can be subtracted
      //   at this position

      if (u.d & u_mask)  u += divisor;

      // Try divisor v at next position

      divisor.divByX();
      u_mask  >>= 1;
   }

   const T mask = ~T() << n;  // Mask for separating quotient from remainder

   q.d = (u.d & mask) >> n;

   r.d = u.d & ~mask;

   #undef u
}


/**
 *  isPrimitive()
 *
 *  Determine if  p  is primitive, i.e.  x^m = 1 mod p  implies m >= 2^deg - 1
 *
 *  See Lidl/Niederreiter, Finite Fields, Theorem 3.18, and
 *      Knuth, TACP, vol 2, 3.2.2, p.30.
 */

template<typename T>
bool L::Polynomial2<T>::isPrimitive () const
{
   const int deg = degree();

   if (deg <= 0)  return false;

   // i)  (-1)^r p[0]  must be a primitive element

   if (! d & 1)  return false;

   const unsigned r = (1 << degree()) - 1;
   PrimeDivisors pd (r);
   const Polynomial2<T> polyx = x();

   // ii)  x^r mod p  must be 1

   if (! powerMod (polyx, r, *this).is1())  return false;

   // iii)  For all divisors rr of r, x^rr mod p  must have positive degree

   while (unsigned prime = pd.next())
   {
      if (powerMod (polyx, r / prime, *this).degree() <= 0)  return false;
   }

   return true;
}

template<typename T>
bool L::Polynomial2<T>::isPrime() const
{
   if (d == 2) return true;   // x
   if (! (d & 1) || d < 2)  return false;  //  x_n+...+0  or  0 or 1

   const unsigned ub = 2 << (degree() / 2);
   for (unsigned i = 3; i < ub; i += 2)
   {
      if ((*this % P(i)).is0())  return false;
   }
   return true;
}

template<typename T>
unsigned char L::Polynomial2<T>::evaluate (unsigned char x) const
{
   if (x == 0)  return d & 1;

   unsigned char res = 0;
   T dd = d;
   while (dd)
   {
      res ^= (dd & 1);
      dd >>= 1;
   }
   return res;
}

/**********************  Polynomial 2 Ring Base  *****************************/

std::ostream & L::operator<< (std::ostream &o, const Polynomial2RingBase &r)
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   ss << "GF2[" << r.getVar() << ']';
   return o << ss.str().c_str();
}

void
L::Polynomial2RingBase::printSuffix (std::ostream& o)
{
   o << "(2)";
}


/**********************  Polynomial 2 Ring  **********************************/

template<typename T>
void
L::Polynomial2Ring<T>::print (std::ostream& o, const type& p) const
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   printShort (ss, p);
   ss << x << " (2)";

   o << ss.str().c_str();
}

/**
 *  next()
 */

template<typename T>
typename L::Polynomial2Ring<T>::type
L::Polynomial2Ring<T>::PrimeGenerator::next()
{
   for (;;)
   {
      type p (n++);
      if (p.isPrime())  return p;
   }
}

/********************  Instantiations  ***************************************/


namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template Polynomial2<X> Polynomial2<X>::operator* \
     (const Polynomial2<X>) const; \
   template void Polynomial2<X>::div \
     (Polynomial2<X>,const Polynomial2<X>,Polynomial2<X> &,Polynomial2<X> &); \
   template bool Polynomial2<X>::isPrimitive () const; \
   template bool Polynomial2<X>::isPrime() const; \
   template unsigned char Polynomial2<X>::evaluate (unsigned char) const; \
   template void Polynomial2<X>::printShort (std::ostream &, char) const; \
   template void Polynomial2Ring<X>::print (std::ostream&, const type&) const; \
   template Polynomial2<X> Polynomial2Ring<X>::PrimeGenerator::next(); \
   template void GF2VectorSpace<X>::print (std::ostream&, const type&) const; \
   template void GF2VectorSpace<X>::printShort (std::ostream&, const type&) const;

   HINTLIB_INSTANTIATE (u32)
#ifdef HINTLIB_U32_NOT_EQUAL_U64
   HINTLIB_INSTANTIATE (u64)
#endif

#undef HINTLIB_INSTANTIATE
}

