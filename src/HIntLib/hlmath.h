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

#ifndef HINTLIB_HLMATH_H
#define HINTLIB_HLMATH_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_CSTDLIB
  #include <cstdlib>
#else
  #include <stdlib.h>
#endif

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

#include <HIntLib/bitop.h>

namespace HIntLib
{

/**
 *  abs()
 */

#if HINTLIB_ABS_STRATEGY == 1
   inline float abs (const float& x)
   {
      return HINTLIB_ABS_FOR_FLOAT (x);
   }

   inline double abs (const double& x)
   {
      return HINTLIB_ABS_FOR_DOUBLE (x);
   }

   inline long double abs (const long double& x)
   {
      return HINTLIB_ABS_FOR_LONG_DOUBLE (x);
   }

   inline int abs (const int& x)
   {
      return HINTLIB_ABS_FOR_INT (x);
   }

   inline long int abs (const long int& x)
   {
      return HINTLIB_ABS_FOR_LONG_INT (x);
   }

 #if 0
 #ifdef HAVE_LONG_LONG_INT
   inline long long int abs (const long long int& x)
   {
      return HINTLIB_ABS_FOR_LONG_LONG_INT (x);
   }
 #endif
 #endif
#else
 #if HINTLIB_ABS_STRATEGY == 2
   using ::std::abs;
 #else
  #if HINTLIB_ABS_STRATEGY == 3
   using ::abs;
  #else
   using ::std::abs;
   using ::abs;
  #endif
 #endif
#endif


/**
 *  Standard math functions
 */

#if HINTLIB_MATH_FUN_STRATEGY == 1

 // float

 #if HINTLIB_SIN_FOR_FLOAT_NUM == 1
   inline float sin  (const float& x) { return ::std::sin  (x); }
   inline float cos  (const float& x) { return ::std::cos  (x); }
   inline float sqrt (const float& x) { return ::std::sqrt (x); }
   inline float log  (const float& x) { return ::std::log  (x); }
   inline float floor(const float& x) { return ::std::floor(x); }
   inline float ceil (const float& x) { return ::std::ceil (x); }
   inline float modf (const float& x, float* y) { return ::std::modf (x, y); }
   inline float pow  (const float& x, const float& y)
      { return ::std::pow (x, y); }
   inline float pow  (const float& x, int k) { return ::std::pow (x, k); }
 #else
 #if HINTLIB_SIN_FOR_FLOAT_NUM == 2
   inline float sin  (const float& x) { return ::sin  (x); }
   inline float cos  (const float& x) { return ::cos  (x); }
   inline float sqrt (const float& x) { return ::sqrt (x); }
   inline float log  (const float& x) { return ::log  (x); }
   inline float floor(const float& x) { return ::floor(x); }
   inline float ceil (const float& x) { return ::ceil (x); }
   inline float modf (const float& x, float* y) { return ::modf (x, y); }
   inline float pow  (const float& x, const float& y) { return ::pow (x, y); }
   inline float pow  (const float& x, int k) { return ::pow (x, k); }
 #else
 #if HINTLIB_SIN_FOR_FLOAT_NUM == 3
   inline float sin  (const float& x) { return ::std::sinf  (x); }
   inline float cos  (const float& x) { return ::std::cosf  (x); }
   inline float sqrt (const float& x) { return ::std::sqrtf (x); }
   inline float log  (const float& x) { return ::std::logf  (x); }
   inline float floor(const float& x) { return ::std::floorf(x); }
   inline float ceil (const float& x) { return ::std::ceilf (x); }
   inline float modf (const float& x, float* y) { return ::std::modff (x, y); }
   inline float pow  (const float& x, const float& y)
      { return ::std::powf (x, y); }
   inline float pow  (const float& x, int k)
      { return ::std::powf (x, k); }
 #else
   inline float sin  (const float& x) { return ::std::sinf  (x); }
   inline float cos  (const float& x) { return ::std::cosf  (x); }
   inline float sqrt (const float& x) { return ::std::sqrtf (x); }
   inline float log  (const float& x) { return ::std::logf  (x); }
   inline float floor(const float& x) { return ::std::floorf(x); }
   inline float ceil (const float& x) { return ::std::ceilf (x); }
   inline float modf (const float& x, float* y) { return ::std::modff (x, y); }
   inline float pow  (const float& x, const float& y)
      { return ::std::powf (x, y); }
   inline float pow  (const float& x, int k)
      { return ::std::powf (x, k); }
 #endif
 #endif
 #endif

 // double

 #if HINTLIB_SIN_FOR_DOUBLE_NUM == 1 || HINTLIB_SIN_FOR_DOUBLE_NUM == 3
   inline double sin  (const double& x) { return ::std::sin  (x); }
   inline double cos  (const double& x) { return ::std::cos  (x); }
   inline double sqrt (const double& x) { return ::std::sqrt (x); }
   inline double log  (const double& x) { return ::std::log  (x); }
   inline double floor(const double& x) { return ::std::floor(x); }
   inline double ceil (const double& x) { return ::std::ceil (x); }
   inline double modf (const double& x, double* y) { return ::std::modf (x,y); }
   inline double pow  (const double& x, const double& y)
      { return ::std::pow (x, y); }
   inline double pow  (const double& x, int k)
      { return ::std::pow (x, k); }
 #else
   inline double sin  (const double& x) { return ::sin  (x); }
   inline double cos  (const double& x) { return ::cos  (x); }
   inline double sqrt (const double& x) { return ::sqrt (x); }
   inline double log  (const double& x) { return ::log  (x); }
   inline double floor(const double& x) { return ::floor(x); }
   inline double ceil (const double& x) { return ::ceil (x); }
   inline double modf (const double& x, double* y) { return ::modf (x, y); }
   inline double pow  (const double& x, const double& y)
      { return ::pow (x, y); }
   inline double pow  (const double& x, int k) { return ::pow (x, k); }
 #endif

 // long double

 #if HINTLIB_SIN_FOR_LONG_DOUBLE_NUM == 1
   inline long double sin  (const long double& x) { return ::std::sin  (x); }
   inline long double cos  (const long double& x) { return ::std::cos  (x); }
   inline long double sqrt (const long double& x) { return ::std::sqrt (x); }
   inline long double log  (const long double& x) { return ::std::log  (x); }
   inline long double floor(const long double& x) { return ::std::floor(x); }
   inline long double ceil (const long double& x) { return ::std::ceil (x); }
   inline long double modf (const long double& x, long double* y)
      { return ::std::modf (x, y); }
   inline long double pow  (const long double& x, const long double& y)
      { return ::std::pow (x, y); }
   inline long double pow  (const long double& x, int k)
      { return ::std::pow (x, k); }
 #else
 #if HINTLIB_SIN_FOR_LONG_DOUBLE_NUM == 2
   inline long double sin  (const long double& x) { return ::sin  (x); }
   inline long double cos  (const long double& x) { return ::cos  (x); }
   inline long double sqrt (const long double& x) { return ::sqrt (x); }
   inline long double log  (const long double& x) { return ::log  (x); }
   inline long double floor(const long double& x) { return ::floor(x); }
   inline long double ceil (const long double& x) { return ::ceil (x); }
   inline long double modf (const long double& x, long double* y)
      { return ::modf (x, y); }
   inline long double pow  (const long double& x, const long double& y)
      { return ::pow (x, y); }
   inline long double pow  (const long double& x, int k)
      { return ::pow (x, k); }
 #else
 #if HINTLIB_SIN_FOR_LONG_DOUBLE_NUM == 3
   inline long double sin  (const long double& x) { return ::std::sinl  (x); }
   inline long double cos  (const long double& x) { return ::std::cosl  (x); }
   inline long double sqrt (const long double& x) { return ::std::sqrtl (x); }
   inline long double log  (const long double& x) { return ::std::logl  (x); }
   inline long double floor(const long double& x) { return ::std::floorl(x); }
   inline long double ceil (const long double& x) { return ::std::ceill (x); }
   inline long double modf (const long double& x, long double* y)
      { return ::std::modfl (x, y); }
   inline long double pow  (const long double& x, const long double& y)
      { return ::std::powl (x, y); }
   inline long double pow  (const long double& x, int k)
      { return ::std::powl (x, k); }
 #else
   inline long double sin  (const long double& x) { return ::std::sinl  (x); }
   inline long double cos  (const long double& x) { return ::std::cosl  (x); }
   inline long double sqrt (const long double& x) { return ::std::sqrtl (x); }
   inline long double log  (const long double& x) { return ::std::logl  (x); }
   inline long double floor(const long double& x) { return ::std::floorl(x); }
   inline long double ceil (const long double& x) { return ::std::ceill (x); }
   inline long double modf (const long double& x, long double* y)
      { return ::std::modfl (x, y); }
   inline long double pow  (const long double& x, const long double& y)
      { return ::std::powl (x, y); }
   inline long double pow  (const long double& x, int k)
      { return ::std::powl (x, k); }
 #endif
 #endif
 #endif

 #define HINTLIB_MN ::HIntLib::

#else
 #if HINTLIB_MATH_FUN_STRATEGY == 2
  #define HINTLIB_MN ::std::
 #else
  #define HINTLIB_MN ::
 #endif
#endif


/**
 *  sqr()
 *
 *  Square of a number
 *
 *  A specialization is provided for Polynomial2<> in polynomial2.h.
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

template<class T> T powInt (T, unsigned)  HINTLIB_GNU_CONST;


/**
 * powerMod()
 * powerModReduce()
 */

template<class T> T powerMod (T, unsigned, T m)  HINTLIB_GNU_CONST;
template<class T> T powerModReduce (T, unsigned, T m)  HINTLIB_GNU_CONST;


/**
 *  logInt()
 */

template<class T> int logInt (T x, T base);   // may throw

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

template<class T> inline
bool approx (T a, T b, T factor = 10.0)
{
   return abs(a - b)
        <= factor * std::numeric_limits<T>::epsilon() * (abs(a) + abs(b));
}


/**
 *  choose()
 *
 *  Calcultes  "a choose b", i.e.  a! / (b! * (a-b)!).
 *
 *  a has to be larger than b.
 *
 *  The routines takes  O(min(b, a-b))  steps and uses no intermediate results
 *  that are larger than the final result.
 */

template<typename T>
inline
T choose (T a, T b)
{
   if (b < 0 || b > a)  return 0;

   if (2 * b > a)  b = a - b;
   T result = 1;
   ++a;
   for (T i = 1; i <= b; ++i)  result = (result * (a-i)) / i;
   return result;
}


/**
 *  radicalInverseFunction ()
 *  radicalInverseFunction2 ()
 */

real radicalInverseFunction  (Index, unsigned base)  HINTLIB_GNU_CONST;
real radicalInverseFunction2 (Index)  HINTLIB_GNU_CONST;

}  // namespace HIntLib

#endif
