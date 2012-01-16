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
#include <HIntLib/modulararithmetic.h>

namespace L = HIntLib;

/******************************  GF2  ****************************************/


std::ostream &
L::operator<< (std::ostream &o, const GF2 &)
{
   return o << "GF2";
}

void
L::GF2::print (std::ostream& o, type x) const
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
L::GF2::printShort (std::ostream& o, type x) const
{
   o << unsigned (x);
}

void
L::GF2::printSuffix (std::ostream& o) const
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
std::ostream& L::operator<< (std::ostream& o, const L::Polynomial2<T> p)
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
      if (p[i])
      {
         if (output) ss << '+';

         switch (i)
         {
         case 0:  ss << '1'; break;
         case 1:  ss << 'x'; break;
         case 2:  ss << "x\262"; break;
         case 3:  ss << "x\263"; break;
         default: ss << "x^" << i;
         }

         output = true;
      }
   }

   if (! output)  ss << '0';

   return o << ss.str().c_str();
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

   if (n == -1)  throw DivisionByZero();

   int loops = m - n;

   if (loops < 0)
   {
      r = u;
      q = zero();

      return;
   }

   // Kill highest order coefficient of divisor

   P divisor (v + xPow (n));

   // align with dividend

   divisor.mulByXPow (loops);

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


template<typename T>
bool L::Polynomial2<T>::isPrimitive () const
{
   const int deg = degree();

   if (deg < 1)  return false;

   const unsigned r = (1 << degree()) - 1;

   const Polynomial2<T> p = x();

   if (! powerMod (p, r, *this).is1())  return false;

   unsigned rr = r;

   for (unsigned prime = 2; ; prime = Prime::next(prime+1))
   {
      if (prime > rr)  break;

      if (rr % prime != 0)  continue;

      do
      {
         rr /= prime;
      }
      while (rr % prime == 0);

      if (powerMod (p, r / prime, *this).degree() <= 0)  return false;
   }

   return true;
}

template<typename T>
bool L::Polynomial2<T>::isIrreducible() const
{
   if (d < 2)  return false;

   const unsigned ub = 2 << (degree() / 2);
   for (unsigned i = 2; i < ub; ++i)
   {
      P r = *this % P(i);
      if (r.is0())  return false;
   }
   return true;
}

template<typename T>
unsigned char L::Polynomial2<T>::evaluate (unsigned char x) const
{
   if (x)
   {
      unsigned char res = 0;
      T dd = d;
      while (dd)
      {
         if (dd & 1)  res ^= 1;
         dd >>= 1;
      }
      return res;
   }
   else
   {
      return d & 1;
   }
}

/**********************  Polynomial 2 Ring Base  *****************************/

std::ostream & L::operator<< (std::ostream &o, const Polynomial2RingBase &)
{
   return o << "GF2[x]";
}

void
L::Polynomial2RingBase::printSuffix (std::ostream& o) const
{
   o << "(2)";
}


/**********************  Polynomial 2 Ring  **********************************/

template<typename T>
void
L::Polynomial2Ring<T>::print (std::ostream& o, const type& x) const
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   ss << x << " (2)";
   o << ss.str().c_str();
}


/********************  Instantiations  ***************************************/


namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template \
   std::ostream& operator<< (std::ostream &, const Polynomial2<X>); \
   template Polynomial2<X> Polynomial2<X>::operator* \
     (const Polynomial2<X>) const; \
   template void Polynomial2<X>::div \
     (Polynomial2<X>,const Polynomial2<X>,Polynomial2<X> &,Polynomial2<X> &); \
   template bool Polynomial2<X>::isPrimitive () const; \
   template bool Polynomial2<X>::isIrreducible() const; \
   template unsigned char Polynomial2<X>::evaluate (unsigned char) const; \
   template void Polynomial2Ring<X>::print (std::ostream&, const type&) const; \
   template void GF2VectorSpace<X>::print (std::ostream&, const type&) const; \
   template void GF2VectorSpace<X>::printShort (std::ostream&, const type&) const;

   HINTLIB_INSTANTIATE (u32)
#ifdef HINTLIB_U32_NOT_EQUAL_U64
   HINTLIB_INSTANTIATE (u64)
#endif

#undef HINTLIB_INSTANTIATE
}

