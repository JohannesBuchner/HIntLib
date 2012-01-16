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

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/prime.h>

#include <HIntLib/hlmath.h>


namespace L = HIntLib;


/**
 *  doPrimeTest()
 *
 *  Checks if a number is prime
 *
 *     isPrime(16) = false
 *     isPrime(17) = true
 */

template<class T>
bool L::Prime::doPrimeTest (T n)
{
   for (T prime = T(2); ; prime = next(prime+1))
   {
      if (n % prime == 0)    return false;
      if (prime * prime > n) return true;
   }
}

template bool L::Prime::doPrimeTest (unsigned);
#ifdef HINTLIB_UNSIGNED_NOT_EQUAL_U32
template bool L::Prime::doPrimeTest (u32);
#endif
#ifdef HINTLIB_U32_NOT_EQUAL_U64
template bool L::Prime::doPrimeTest (u64);
#endif


/**
 *  searchForNextPrime()
 *
 *  used by next(), if n > MAX_NUM_FOR_NEXT_PRIME
 */

template<class T>
T L::Prime::searchForNextPrime (T n)
{
   for (T i = n; ;++i)  if (doPrimeTest(i))  return i;
}

template unsigned L::Prime::searchForNextPrime (unsigned);
#ifdef HINTLIB_UNSIGNED_NOT_EQUAL_U32
template L::u32   L::Prime::searchForNextPrime (u32);
#endif
#ifdef HINTLIB_U32_NOT_EQUAL_U64
template L::u64   L::Prime::searchForNextPrime (u64);
#endif


/**
 *  eulerPhi()
 *
 *  Euler's Phi function.  Counts the number of m with 1<=m<n and gcd(m,n)=1
 *
 *     eulerPhi (0) = 0
 *     eulerPhi (1) = 1
 *     eulerPhi (2) = 1
 *     eulerPhi (3) = 2
 *     eulerPhi (4) = 2
 *     eulerPhi (5) = 4
 *     eulerPhi (6) = 2
 *     eulerPhi (7) = 6
 *         :
 */

template<class T>
T L::Prime::eulerPhi (T n)
{
   if (n == 0)  return 0;

   unsigned phi = 1;
   unsigned prime = 1;

   // test all primes

   while (n > 1)
   {
      prime = nextPrimeDivisor (n, prime);
      n /= prime;

      phi *= (prime - 1);

      while (n % prime == 0)
      {
         n /= prime;
         phi *= prime;
      }
   }

   return phi;
}

template unsigned L::Prime::eulerPhi (unsigned);
#ifdef HINTLIB_UNSIGNED_NOT_EQUAL_U32
template L::u32   L::Prime::eulerPhi (u32);
#endif
#ifdef HINTLIB_U32_NOT_EQUAL_U64
template L::u64   L::Prime::eulerPhi (u64);
#endif


/**
 *  isPrimitiveRoot()
 *
 *  See TACP, vol 2, 3.2.1.2
 */

bool L::isPrimitiveRoot (unsigned a, unsigned p)
{
   const unsigned q = p - 1;

   PrimeDivisors pd (q);

   while (unsigned prime = pd.next())
   {
      if (powerMod (a, q / prime, p) == 1)  return false;
   }

   return true;
}


/**
 *  isPrimePower()
 *  factorPrimePower ()
 *
 *  Given a number  x = p^n , with  p  prime, factorPrimePower() determines
 *  p and n.
 *
 *  If  x  is not a prime power, NotAPrimePower is thrown.
 */

bool L::Prime::isPrimePower (unsigned n, unsigned &_prime, unsigned &_power)
{
   if (Prime::test (n))
   {
      _prime = n; _power = 1;
   }
   else
   {
      if (n < 2)  return false;

      unsigned prime = 2;

      while (n % prime != 0)
      {
         if (prime * prime > n)  return false;
         prime = Prime::next (prime + 1);
      }

      unsigned power = 0;
      unsigned s = n;

      do
      {
         if (s % prime != 0)  return false;
         s /= prime;
         ++power;
      }
      while (s > 1);

      _prime = prime; _power = power;
   }

   return true;
}

void L::Prime::factorPrimePower (unsigned n, unsigned &prime, unsigned &power)
{
   if (! isPrimePower (n, prime, power))  throw NotAPrimePower (n);
}

/**
 *  throwPrimeNumberNth ()
 */

void L::Prime::throwPrimeNumberNth (unsigned n)
{
  throw PrimeNumberNth (n);
}


/**
 *  nextPrimeDivisor()
 *
 *  Finds the next prime divisor of n, assuming that n does not have any
 *  divisor <= prime.
 */

unsigned L::Prime::nextPrimeDivisor (unsigned n, unsigned prime)
{
   // check all primes from the table

   while (prime < Prime::MAX_N)
   {
      prime = Prime::nextPrimeArray[prime + 1];
      if (n % prime == 0)  return prime;
      if (prime * prime > n)  return n;
   }

   // increment prime, but make sure that prime = 6 k + 1

   if (prime % 6 == 5)  prime += 2;

   // Check further divisors of the form  6 k + 1  and  6 k + 5

   for (;;)
   {
      if (n % prime == 0)     return prime;
      if (prime * prime > n)  return n;
      prime += 4;
      if (n % prime == 0)     return prime;
      if (prime * prime > n)  return n;
      prime += 2;
   }
}


/**
 *  PrimeDivisors
 */

unsigned L::PrimeDivisors::next()
{
   if (n < 2)  return 0;
   prime = Prime::nextPrimeDivisor (n, prime);
   do { n /= prime; } while (n % prime == 0);
   return prime;
}

unsigned L::PrimeDivisors::next (unsigned& e)
{
   if (n < 2)  return 0;
   prime = Prime::nextPrimeDivisor (n, prime);
   e = 0;
   do { n /= prime; ++e; } while (n % prime == 0);

   return prime;
}


/**
 *  factor()
 */

void
L::Prime::factor (Factorization& f, unsigned n)
{
   PrimeDivisors pd (n);
   unsigned e;
   while (unsigned prime = pd.next(e))
   {
      f.push_back (std::make_pair (prime, e));
   }
}


/**
 *  gcd()  --  two multipliers
 */

template<typename T>
T L::gcd (T v3, T u3, T& mv, T& mu)
{
   T u1 = 1;
   T u2 = 0;

   T v1 = 0;
   T v2 = 1;

   for (;;)
   {
      if (! v3)
      {
         mu = u1; mv = u2;
         return u3;
      }

      T q = u3 / v3; u3 %= v3;
      u1 -= q * v1;
      u2 -= q * v2;

      if (! u3)
      {
         mu = v1; mv = v2;
         return v3;
      }

      q = v3 / u3; v3 %= u3;
      v1 -= q * u1;
      v2 -= q * u2;
   }
}


/**
 *  gcd()  --  one multiplier
 */

template<typename T>
T L::gcd (T v3, T u3, T& mv)
{
   T u2 = 0;
   T v2 = 1;

   for (;;)
   {
      if (! v3)
      {
         mv = u2;
         return u3;
      }

      u2 -= u3 / v3 * v2;
      u3 %= v3;

      if (! u3)
      {
         mv = v2;
         return v3;
      }

      v2 -= v3 / u3 * u2;
      v3 %= u3;
   }
}

#define HINTLIB_INSTANTIATE_GCD(X) \
   template X gcd (X, X, X&); \
   template X gcd (X, X, X&, X&);

namespace HIntLib
{
   HINTLIB_INSTANTIATE_GCD(int)
   HINTLIB_INSTANTIATE_GCD(long)
#ifdef HINTLIB_HAVE_LONG_LONG_INT
   HINTLIB_INSTANTIATE_GCD (long long)
#endif
}


