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
#pragma implementation "gf2.h"
#pragma implementation "gf2vectorspace.h"
#endif

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/polynomial2.h>
#include <HIntLib/gf2vectorspace.h>

#include <HIntLib/integerring.h>
#include <HIntLib/output.h>
#include <HIntLib/linearalgebra.h>

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
   Private::Printer ss (o);
   ss << unsigned (x) << " (2)";
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
   Private::Printer ss (o);
   ss << "GF2^" << v.dimension();
   return o;
}

template<typename T>
void
L::GF2VectorSpace<T>::printShort (std::ostream& o, const type& v) const
{
   Private::Printer ss (o);

   ss << '(';
   for (unsigned i = 0; i < dim; ++i)
   {
      if (i > 0)  ss << ',';
      ss << unsigned (coord(v, i));
   }
   ss << ')';
}

template<typename T>
void
L::GF2VectorSpace<T>::print (std::ostream& o, const type& v) const
{
   Private::Printer ss (o);

   printShort (ss, v);
   ss << " (2)";
}


/*************************  Polynomial2  *************************************/


/**
 *  Prints a polynomial to an output stream.
 */

template<typename T>
void L::Polynomial2<T>::printShort (
      std::ostream& o, char var, PrintShortFlag f) const
{
   // The zero-polynomial is a special case

   if (is0())
   {
      o << "0";
      return;
   }

   // count non-zero terms

   int nonZeroTerms = 0;
   T dd = d;

   while (dd && nonZeroTerms < 2)
   {
      if (dd & 1)  ++nonZeroTerms;
      dd >>= 1;
   }

   Private::Printer ss (o);

   bool needPlus = o.flags() & o.showpos;

   if ((f & FIT_FOR_MUL) && nonZeroTerms >= 2)
   {
      if (needPlus)
      {
         ss << '+';
         needPlus = false;
      }
      ss << '(';
   }

   for (int i = degree(); i >= 0; --i)
   {
      if ((*this)[i])
      {
         if (needPlus) ss << '+';

         switch (i)
         {
         case 0:  ss << '1'; break;
         case 1:  ss << var; break;
         case 2:  ss << var << '\262'; break;
         case 3:  ss << var << '\263'; break;
         default: ss << var << '^' << unsigned(i);  // avoid sign from showpos
         }

         needPlus = true;
      }
   }

   if ((f & FIT_FOR_MUL) && nonZeroTerms >= 2)  ss << ')';
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
 *  sqr()
 */

template<typename T>
L::Polynomial2<T> L::sqr (Polynomial2<T> p)
{
   T d = T(p);

   if (d & (~T(0) << std::numeric_limits<T>::digits / 2))  throwOverflow();

   T result = 0;
   T mask = 1;

   while (d)
   {
      if (d & 1)  result |= mask;
      d >>= 1;
      mask <<= 2;
   }

   return Polynomial2<T>(result);
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

   for ( ; loops >= 0; --loops)
   {
      // If the current highest coefficient of u is set, v can be subtracted
      //   at this position

      if (uu.d & u_mask)  uu += divisor;

      // Try divisor v at next position

      divisor.divByX();
      u_mask  >>= 1;
   }

   const T mask = ~T() << n;  // Mask for separating quotient from remainder

   q.d = (uu.d & mask) >> n;

   r.d = uu.d & ~mask;
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


/**
 *  isPrime()
 */

template<typename T>
bool L::Polynomial2<T>::isPrime() const
{
   // handle the basic cases

   if (d < 4)  return d & 2;  // degree 0 and 1

   if (! (d & 1))  return false;  // check for ...+0

   if (d < 256)
   {
      int deg = this->degree();

      for (unsigned i = 3; ; i += 2)
      {
         if (2 * P(i).degree() > deg)  return true;
         if ((*this % P(i)).is0())  return false;
      }
   }

   // check for divisors up to degree 3

   for (unsigned i = 3; i < 16; i += 2)
   {
      if ((*this % P(i)).is0())  return false;
   }

   // Make sure it is square free

   if (! isSquarefree())  return false;

   // Berlekamp's algorithm

   P q (1);
   const int numRows = degree() - 1;
   const T overflow = T(1) << numRows + 1;
   T matrix [std::numeric_limits<T>::digits];

   for (int row = 0; row < numRows; ++row)
   {
      // multiply  q  by  x  and reduce (2 times)

      q.d <<= 1;
      if (q.d & overflow)  q -= *this;
      q.d <<= 1;
      if (q.d & overflow)  q -= *this;

      // subtract I

      matrix[row] = q ^ (2 << row);
   }

   // degree - rank (of the untruncated matrix) gives the number of factors

   return isLinearlyIndependent (matrix, matrix + numRows);
}


/**
 * evaluate()
 */

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


/**
 * powInt()
 */

template<typename T>
L::Polynomial2<T> L::powInt (Polynomial2<T> x, unsigned exponent)
{
   Polynomial2<T> result (1);

   for(;;)
   {
      if (exponent & 1)  result *= x;
      if ((exponent >>= 1) == 0)  return result;
      x = sqr(x);
   }
}


/**
 * powerMod()
 */

template<class T>
L::Polynomial2<T>
L::powerMod (Polynomial2<T> x, unsigned exponent, Polynomial2<T> m)
{
   if (x.is0())  return Polynomial2<T>();
   Polynomial2<T> result(1);

   for(;;)
   {
      if (exponent & 1)  result = (result * x) % m;
      if ((exponent >>= 1) == 0)  return result;
      x = sqr(x) % m;
   }
}


/**********************  Polynomial 2 Ring Base  *****************************/


std::ostream & L::operator<< (std::ostream &o, const Polynomial2RingBase &r)
{
   Private::Printer ss (o);

   ss << "GF2[";
   r.printVariable (ss);
   ss << ']';

   return o;
}

void
L::Polynomial2RingBase::printSuffix (std::ostream& o)
{
   o << "(2)";
}

void
L::Polynomial2RingBase::printVariable (std::ostream& o) const
{
   o << var;
}


/**********************  Polynomial 2 Ring  **********************************/

template<typename T>
void
L::Polynomial2Ring<T>::print (std::ostream& o, const type& p) const
{
   Private::Printer ss (o);

   printShort (ss, p);
   ss << x << " (2)";
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


/**
 *  squarefreeFactor()
 */

template<typename T>
L::Polynomial2RingBase::unit_type
L::Polynomial2Ring<T>::squarefreeFactor (Factorization& f, type p) const
{
   if (p.is0())  throwDivisionByZero();
   if (p.isConstant())  return unit_type();

   // Step 1  --  Initialize

   unsigned e = 0;
   unsigned k = 0;

   // v ... each factor appears here once (module char problem)
   // t ... the remaining factors
   // v * t ... the polynomial to be factored

   type t = genGcd (*this, derivative (p), p);
   type v = p / t;

   for(;;)
   {
      if (v.isConstant())
      {
         if (t.isConstant())  break;

         T t0d = 0;
         T mask = 1;
         T td = T(t);

         while (td)
         {
            if (td & 1)  t0d |= mask;
            mask <<= 1;
            td >>= 2;
         }
         type t0 (t0d);

         t = genGcd (*this, derivative (t0), t0);
         v = t0 / t;
         ++e;
         k = 0;
      }
      else
      {
         // Step 4  -- Special case

         ++k;
         if (k & 1 == 0)  // Special treatment for  k  a multiple of char
         {
            t /= v;
            ++k;
         }

         // Step 5  --  Compute factor

         // w = product of factors apearing more than once

         type w = genGcd (*this, t, v);

         if (v.degree() > w.degree())
         {
            f.push_back (std::make_pair (v / w, k << e));
         }
         v = w;
         t /= v;
      }
   }

   return unit_type();
}

/********************  Instantiations  ***************************************/

#include <HIntLib/gcd.tcc>

namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template Polynomial2<X> Polynomial2<X>::operator* \
      (const Polynomial2<X>) const; \
   template Polynomial2<X> sqr(Polynomial2<X>); \
   template void Polynomial2<X>::div \
      (Polynomial2<X>,const Polynomial2<X>,Polynomial2<X> &,Polynomial2<X> &); \
   template bool Polynomial2<X>::isPrimitive () const; \
   template bool Polynomial2<X>::isPrime() const; \
   template unsigned char Polynomial2<X>::evaluate (unsigned char) const; \
   template void Polynomial2<X>::printShort \
      (std::ostream &, char, PrintShortFlag) const; \
   HINTLIB_INSTANTIATE_GENGCD(Polynomial2Ring<X>) \
   template Polynomial2<X> powInt (Polynomial2<X>, unsigned); \
   template Polynomial2<X> powerMod (Polynomial2<X>, unsigned, Polynomial2<X>); \
   template Polynomial2RingBase::unit_type \
      Polynomial2Ring<X>::squarefreeFactor (Factorization&, type) const; \
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

