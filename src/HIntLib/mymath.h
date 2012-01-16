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
 *  MyMath
 *
 *  Defines a number of common mathematical functions missing in the standard
 *  library.
 */

#ifndef HINTLIB_MYMATH_H
#define HINTLIB_MYMATH_H 1

#include <stdlib.h>

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_CMATH
  #include <cmath>
#else
  #include <math.h>
#endif

#ifdef HINTLIB_HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/fallback_limits.h>
#endif


namespace HIntLib
{

#ifdef HINTLIB_GLOBAL_ABS_AVAILABLE
   using ::abs;
#else
 #ifdef HINTLIB_STD_ABS_AVAILABLE
   using std::abs;
   using ::abs;
 #else
   using ::abs;
   inline double      abs (double x)       { return fabs (x); }
   inline float       abs (float x)        { return fabsf(x); }
   inline long double abs (long double x)  { return fabsl(x); }
 #endif
#endif


/**
 *  sqr()
 *
 *  Square of a number
 *
 *  Works for
 *     - integers
 *     - floats
 *     - GenPolynomial2
 */

template<class T> inline T sqr (T x)
{
   return x * x;
}

/**
 *  cube()
 *
 *  Cube of a number
 *
 *  Works for
 *     - integers
 *     - floats
 *     - GenPolynomial2
 */

template<class T> inline T cube (T x)
{
   return x * x * x;
}


/**
 *  powInt()
 */

template<class T> T powInt (T x, unsigned exponent);
template<class A>
// not inline
typename A::type powInt (const A &a, typename A::type x, unsigned exponent)
{
   typename A::type result (a.one());

   for(;;)
   {
      if (exponent & 1)  a.mulBy (result, x);
      if ((exponent >>= 1) == 0)  return result;
      a.mulBy (x, x);
   }
}


/**
 * powerMod()
 */

template<class T> T powerMod (T x, unsigned e, T m);

template<class A> typename A::type
// not inline
powerMod (const A &a, typename A::type x, unsigned exponent,
                const typename A::type &m)
{
   typename A::type result (a.one());
   typename A::type q;
   for (;;)
   {
      if (exponent & 1)  a.div (a.mul (result, x), m, q, result);
      if ((exponent >>= 1) == 0)  return result;
      a.div (a.mul (x, x), m, q, x);
   }
}


/**
 *  logInt()
 */

template<class T> int logInt (T x, T base);

template<> inline int logInt (unsigned char x, unsigned char b)
   { return logInt (unsigned(x), unsigned (b)); }
template<> inline int logInt (unsigned short x, unsigned short b)
   { return logInt (unsigned(x), unsigned (b)); }


/**
 *  roundUpToPower()
 *
 *  Finds the smallest number >= _x_ that is a power of _base_
 */

template<class T> inline T roundUpToPower (T x, T base)
{
   if (base <= T(1))  throw InvalidLogBase (base);

   T result (1);

   while (result < x)  result *= base;

   return result;
}

/**
 *  roundDownToPower()
 *
 *  Finds the greatest number <= _x_ that is a power of _base?
 *
 *  roundDownToPower(0,base) = 0
 */

template<class T> inline T roundDownToPower (T x, T base)
{
   return powInt (base, logInt(x, base));
}

/**
 *  digitsRepresentable()
 *
 *  Determines the number of base-_base_ digits that can be stored in a T
 */

template<class T> inline unsigned digitsRepresentable(T base)
{
   if (base == 2)  return std::numeric_limits<T>::digits;    // base=2

   // check for base = 2^k

   int k = ms1(base);

   if (T(1) << k == base)  return std::numeric_limits<T>::digits / k;

   return logInt (std::numeric_limits<T>::max(), base);
}


/**
 *  even()   odd()
 *
 *  Determine if a number is even or odd
 */

template<class T> inline bool odd (T x)  HINTLIB_GNU_CONST;
template<class T> inline bool odd (T x)
{
   return x % 2;
}

template<class T> inline bool even (T x)  HINTLIB_GNU_CONST;
template<class T> inline bool even (T x)
{
   return ! odd(x);
}


/**
 *  approx()
 *
 *  Determine if two (floating-point) numbers are approximately equal
 */

template <class T> inline bool approx (T a, double b)  HINTLIB_GNU_CONST;
template <class T> inline bool approx (T a, double b)
{
   return abs(a - b) < std::numeric_limits<T>::epsilon() * (abs(a) + abs(b));
}

template <class T> inline bool approx (T a, double b, double )HINTLIB_GNU_CONST;
template <class T> inline bool approx (T a, double b, double factor)
{
   return abs(a - b)
        < factor * std::numeric_limits<T>::epsilon() * (abs(a) + abs(b));
}

}  // namespace HIntLib

#endif
