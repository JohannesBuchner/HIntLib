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

#ifndef HINTLIB_HAVE_LIMITS
#pragma implementation "fallback_limits.h"
#endif
#endif

#include <HIntLib/hlmath.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;


/**
 *  powInt()
 *
 *  Arbitrary powers to an integer exponent
 *
 *  A specialization is provided for Polynomial2<>
 */

template<class T>
T L::powInt (T x, unsigned exponent)
{
   T result(1);

   for(;;)
   {
      if (exponent & 1)  result *= x;
      if ((exponent >>= 1) == 0)  return result;
      x *= x;
   }
}


namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) template X powInt (X, unsigned);

   HINTLIB_INSTANTIATE (int)
   HINTLIB_INSTANTIATE (unsigned)
   HINTLIB_INSTANTIATE (unsigned long)
#ifdef HINTLIB_HAVE_UNSIGNED_LONG_LONG_INT
   HINTLIB_INSTANTIATE (unsigned long long)
#endif
#undef HINTLIB_INSTANTIATE
}


/**
 *  powerMod()
 */

template<class T>
T L::powerMod (T x, unsigned exponent, T p)
{
   if (! x)  return T();

   const T mask = ~T() << (std::numeric_limits<T>::digits / 2);
   T result (1);

   for (;;)
   {
      if (exponent & 1)
      {
         result *= x;
         if (result & mask)  result %= p;
      }
      if ((exponent >>= 1) == 0)  return result >= p ? result % p : result;
      x *= x;
      if (x & mask)  x %= p;
   }
}

namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) template X powerMod (X, unsigned, X);

   HINTLIB_INSTANTIATE (unsigned)
   HINTLIB_INSTANTIATE (unsigned long)
#undef HINTLIB_INSTANTIATE
}


/**
 *  powerModReduce()
 */

template<class T>
T L::powerModReduce (T x, unsigned exponent, T p)
{
   if (! x)  return T();
   if (exponent >= p - 1)  exponent %= (p - 1);

   const T mask = ~T() << (std::numeric_limits<T>::digits / 2);
   T result (1);

   for (;;)
   {
      if (exponent & 1)
      {
         result *= x;
         if (result & mask)  x %= p;
      }
      if ((exponent >>= 1) == 0)  return result >= p ? result % p : result;
      x *= x;
      if (x & mask)  x %= p;
   }
}

namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) template X powerModReduce (X, unsigned, X);

   HINTLIB_INSTANTIATE (unsigned)
   HINTLIB_INSTANTIATE (unsigned long)
#undef HINTLIB_INSTANTIATE
}


/**
 *  logInt()
 *
 *  Determines the base-_base_ logarithm of _x_, truncated to an int
 *
 *  logInt(0, base)  := -1
 */

template<class T> int L::logInt (T x, T base)
{
   if (base < 2)  throw InvalidLogBase (base);
   if (x == 0)  return -1;
   x /= base;
   
   int result = 0;
   T test = 1;

   while (test <= x)
   {
      test *= base;
      ++ result;
   }

   return result;
}

namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) template int logInt (X, X);

   HINTLIB_INSTANTIATE (unsigned)
   HINTLIB_INSTANTIATE (unsigned long)
#ifdef HINTLIB_HAVE_UNSIGNED_LONG_LONG_INT
   HINTLIB_INSTANTIATE (unsigned long long)
#endif
#undef HINTLIB_INSTANTIATE
}


/**
 *  radicalInverseFunction  ()
 *  radicalInverseFunction2 ()
 */

L::real L::radicalInverseFunction (Index n, unsigned base)
{
   const real realBase = real (1.0) / real (base);

   real x = 0.0;
   real b = realBase;

   while (n)
   {
      x += b * (n % base);
      n /= base;
      b *= realBase;
   }

   return x;
}

L::real L::radicalInverseFunction2 (Index n)
{
   const real realBase = 0.5;

   real x = 0.0;
   real b = realBase;

   while (n)
   {
      if (n & 1)  x += b;
      n >>= 1;
      b *= realBase;
   }

   return x;
}

#include <HIntLib/bitop.tcc>

namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template X thread (X, X, unsigned, unsigned); \
   template void unthread (X, X&, X&, unsigned, unsigned); \
   template X threadn (X*, unsigned); \
   template void unthreadn (X, X*, unsigned); \
   template X threadinf (X*, unsigned); \
   template unsigned unthreadinf (X, X*);
   
   HINTLIB_INSTANTIATE (unsigned)
   HINTLIB_INSTANTIATE (unsigned long)
#ifdef HINTLIB_HAVE_UNSIGNED_LONG_LONG_INT
   HINTLIB_INSTANTIATE (unsigned long long)
#endif
#undef HINTLIB_INSTANTIATE
}



