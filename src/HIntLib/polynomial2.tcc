/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration
 *
 *  Copyright (C) 2002  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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

#ifndef HINTLIB_POLYNOMIAL_2_TCC
#define HINTLIB_POLYNOMIAL_2_TCC 1

#include <HIntLib/polynomial2.h>

#include <HIntLib/output.h>
#include <HIntLib/linearalgebra2.h>
#include <HIntLib/prime.h>
#include <HIntLib/polynomialbase.h>
#include <HIntLib/gcd.tcc>


/*************************  Polynomial2  *************************************/


/**
 *  printShort()
 *
 *  Prints a polynomial to an output stream.
 */

template<typename T>
void
HIntLib::Polynomial2<T>::printShort (
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

         if (i == 0)  ss << '1';
         else
         {
            ss <<var;
            if (i >= 2)  ss.power (i);
         }

         needPlus = true;
      }
   }

   if ((f & FIT_FOR_MUL) && nonZeroTerms >= 2)  ss << ')';
}

#ifdef HINTLIB_BUILD_WCHAR
template<typename T>
void
HIntLib::Polynomial2<T>::printShort (
      std::wostream& o, wchar_t wvar, PrintShortFlag f) const
{
   // The zero-polynomial is a special case

   if (is0())
   {
      o << L"0";
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

   Private::WPrinter ss (o);

   bool needPlus = o.flags() & o.showpos;

   if ((f & FIT_FOR_MUL) && nonZeroTerms >= 2)
   {
      if (needPlus)
      {
         ss << L'+';
         needPlus = false;
      }
      ss << L'(';
   }

   for (int i = degree(); i >= 0; --i)
   {
      if ((*this)[i])
      {
         if (needPlus) ss << L'+';

         if (i == 0)  ss << L'1';
         else
         {
            ss << wvar;
            if (i >= 2)  ss.power (i);
         }

         needPlus = true;
      }
   }

   if ((f & FIT_FOR_MUL) && nonZeroTerms >= 2)  ss << L')';
}
#endif


/**
 *  Multiply two polynomials
 */

template<typename T>
HIntLib::Polynomial2<T>
HIntLib::Polynomial2<T>::operator* (const P p) const
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
HIntLib::Polynomial2<T>
HIntLib::sqr (Polynomial2<T> p)
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
void
HIntLib::Polynomial2<T>::div (P u, const P v, P &q, P &r)
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
bool
HIntLib::Polynomial2<T>::isPrimitive () const
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
bool
HIntLib::Polynomial2<T>::isPrime() const
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
unsigned char
HIntLib::Polynomial2<T>::evaluate (unsigned char x) const
{
   return (x == 0) ? (d & 1) : parity (d);
}


/**
 * powInt()
 */

template<typename T>
HIntLib::Polynomial2<T>
HIntLib::powInt (Polynomial2<T> x, unsigned exponent)
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
HIntLib::Polynomial2<T>
HIntLib::powerMod (Polynomial2<T> x, unsigned exponent, Polynomial2<T> m)
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


/**********************  Polynomial 2 Ring  **********************************/


/**
 *  print()
 */

template<typename T>
void
HIntLib::Polynomial2Ring<T>::print (std::ostream& o, const type& p) const
{
   Private::Printer ss (o);

   printShort (ss, p);
   ss << " (2)";
}

#ifdef HINTLIB_BUILD_WCHAR
template<typename T>
void
HIntLib::Polynomial2Ring<T>::print (std::wostream& o, const type& p) const
{
   Private::WPrinter ss (o);

   printShort (ss, p);
   ss << L" (2)";
}
#endif

/**
 *  next()
 */

template<typename T>
typename HIntLib::Polynomial2Ring<T>::type
HIntLib::Polynomial2Ring<T>::PrimeGenerator::next()
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
HIntLib::Private::Polynomial2RingBase::unit_type
HIntLib::Polynomial2Ring<T>::squarefreeFactor (Factorization& f, type p) const
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

/**
 *  factor()
 */

template<typename T>
HIntLib::Private::Polynomial2RingBase::unit_type
HIntLib::Polynomial2Ring<T>::factor (Factorization& f, type p) const
{
   // Performe squarefree factorization

   Factorization sff;
   squarefreeFactor (sff, p);

   // Final splitting of all squarefree factors

   for (typename Factorization::iterator it = sff.begin(); 
        it != sff.end(); ++it)
   {
      type& p = it->first;
      int degree = p.degree();

      // Check for trivial cases

      if (degree == 1)   // degree 1 is always prime
      {
         f.push_back (std::make_pair (p, it->second));
         continue;
      }

      // There is a large chance for a linear factor.

      if ((T(p) & T(1)) == 0)
      {
         p.divByX();
         f.push_back (std::make_pair (x(), it->second));
         --degree;
      }

      {
         type rest;
         if (isDivisor (p, type(3), rest))
         {
            f.push_back (std::make_pair (type(3), it->second));
            p = rest;
            --degree;
         }
      }

      if (degree == 0)  continue;

      // Everything with degree less than 4 must be prime now

      if (degree < 4)
      {
         f.push_back (std::make_pair (p, it->second));
         continue;
      }

      // Calculate Berlekamp's matrix  M - I

      T q (1);
      const T overflow = T(1) << degree;
      T matrix [std::numeric_limits<T>::digits];

      for (int row = 0; row < degree; ++row)
      {
         // Set and subtract I

         matrix[row] = q ^ (1 << row);

         // multiply  q  by  x  and reduce (2 times)

         q <<= 1;
         if (q & overflow)  q ^= T(p);
         q <<= 1;
         if (q & overflow)  q ^= T(p);
      }

      // Calculate null space

      T ns [std::numeric_limits<T>::digits];
      T* nsp = ns;

      unsigned numFactors
         = nullSpaceT (matrix, matrix + degree, ns);

      // If  p  does not factor, we are done

      if (numFactors == 1)
      {
         f.push_back (std::make_pair (p, it->second));
         continue;
      }

      // Split  p  until we have  numFactors  factors

      type polys [std::numeric_limits<T>::digits];
      polys[0] = p;
      unsigned numPolys = 1;

      for (;;)
      {
         ++nsp;  // Get next splitter

         unsigned ub = numPolys;
         for (unsigned i = 0; i < ub; ++i)
         {
            type& currentPoly = polys[i];
            if (currentPoly.degree() < 4)  continue;

            bool success = false;

            for (unsigned element = 0; element < 2; ++element)
            {
               type splitter (*nsp ^ element);

               const type split = genGcd (*this, currentPoly, splitter);

               if (split.degree() >= 1 && split.degree() < currentPoly.degree())
               {
                  success = true;
                  polys [numPolys++] = split;
               }
            }

            if (success)
            {
               currentPoly = polys [--numPolys];
               if (numPolys == numFactors)  goto end;
            }
         }
      }
end:

      // copy irreducible factors to output Factorization

      for (unsigned i = 0; i < numPolys; ++i)
      {
         f.push_back (std::make_pair (polys[i], it->second));
      }
   }

   // return the standard unit

   return unit_type();
}


/********************  Instantiations  ***************************************/


#ifdef HINTLIB_BUILD_WCHAR
#define HINTLIB_INSTANTIATE_POLYNOMIAL2_W(X) \
   template void Polynomial2<X >::printShort \
      (std::wostream &, wchar_t, PrintShortFlag) const; \
   template void Polynomial2Ring<X >::print \
      (std::wostream&, const type&) const;
#else
#define HINTLIB_INSTANTIATE_POLYNOMIAL2_W(X)
#endif

#define HINTLIB_INSTANTIATE_POLYNOMIAL2(X) \
   HINTLIB_INSTANTIATE_POLYNOMIAL2_W(X) \
   template Polynomial2<X > Polynomial2<X >::operator* \
      (const Polynomial2<X >) const; \
   template Polynomial2<X > sqr(Polynomial2<X >); \
   template void Polynomial2<X >::div \
      (Polynomial2<X >,const Polynomial2<X >, \
       Polynomial2<X > &,Polynomial2<X > &); \
   template bool Polynomial2<X >::isPrimitive () const; \
   template bool Polynomial2<X >::isPrime() const; \
   template unsigned char Polynomial2<X >::evaluate (unsigned char) const; \
   template void Polynomial2<X >::printShort \
      (std::ostream &, char, PrintShortFlag) const; \
   HINTLIB_INSTANTIATE_GENGCD(Polynomial2Ring<X >) \
   template Polynomial2<X > powInt (Polynomial2<X >, unsigned); \
   template Polynomial2<X > \
      powerMod (Polynomial2<X >, unsigned, Polynomial2<X >); \
   template Private::Polynomial2RingBase::unit_type \
      Polynomial2Ring<X >::squarefreeFactor (Factorization&, type) const; \
   template Private::Polynomial2RingBase::unit_type \
      Polynomial2Ring<X >::factor (Factorization&, type) const; \
   template void Polynomial2Ring<X >::print \
      (std::ostream&, const type&) const; \
   template Polynomial2<X > Polynomial2Ring<X >::PrimeGenerator::next();

#endif

