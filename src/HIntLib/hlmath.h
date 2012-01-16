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

/**
 *  HLMath
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
 *  Mathematical constants
 *
 *  For some reason MSVC does not define this
 */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef HINTLIB_HAVE_LONG_DOUBLE
#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L
#endif
#endif

template<typename T>
struct Constants
{
  static T pi ()  { return M_PI; }
};

#ifdef HINTLIB_HAVE_LONG_DOUBLE
template<>
struct Constants<long double>
{
  static long double pi ()  { return M_PIl; }
};
#endif


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

 #ifdef HINTLIB_HAVE_LONG_DOUBLE
   inline long double abs (const long double& x)
   {
      return HINTLIB_ABS_FOR_LONG_DOUBLE (x);
   }
 #endif

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
   #undef HINTLIB_INCLUDE_STD_MATH_FUNCTIONS
   #define HINTLIB_INCLUDE_STD_MATH_FUNCTIONS
 #else
 #if HINTLIB_SIN_FOR_FLOAT_NUM == 2
   #undef HINTLIB_INCLUDE_GLOBAL_MATH_FUNCTIONS
   #define HINTLIB_INCLUDE_GLOBAL_MATH_FUNCTIONS
 #else
 #if HINTLIB_SIN_FOR_FLOAT_NUM == 3
   inline float sin   (const float& x) { return ::std::sinf   (x); }
   inline float cos   (const float& x) { return ::std::cosf   (x); }
   inline float sqrt  (const float& x) { return ::std::sqrtf  (x); }
   inline float log   (const float& x) { return ::std::logf   (x); }
   inline float exp   (const float& x) { return ::std::expf   (x); }
   inline float floor (const float& x) { return ::std::floorf (x); }
   inline float ceil  (const float& x) { return ::std::ceilf  (x); }
   inline float modf  (const float& x, float* y) { return ::std::modff (x, y); }
   inline float pow   (const float& x, const float& y)
      { return ::std::powf (x, y); }
 #else
 #if HINTLIB_SIN_FOR_FLOAT_NUM == 4
   inline float sin   (const float& x) { return ::std::sinf   (x); }
   inline float cos   (const float& x) { return ::std::cosf   (x); }
   inline float sqrt  (const float& x) { return ::std::sqrtf  (x); }
   inline float log   (const float& x) { return ::std::logf   (x); }
   inline float exp   (const float& x) { return ::std::expf   (x); }
   inline float floor (const float& x) { return ::std::floorf (x); }
   inline float ceil  (const float& x) { return ::std::ceilf  (x); }
   inline float modf  (const float& x, float* y) { return ::std::modff (x, y); }
   inline float pow   (const float& x, const float& y)
      { return ::std::powf (x, y); }
 #else
   // for HINTLIB_SIN_FOR_FLOAT_NUM == 0, we do nothing
 #endif
 #endif
 #endif
 #endif

 // long double

 #ifdef HINTLIB_HAVE_LONG_DOUBLE
 #if HINTLIB_SIN_FOR_LONG_DOUBLE_NUM == 1
   #undef HINTLIB_INCLUDE_STD_MATH_FUNCTIONS
   #define HINTLIB_INCLUDE_STD_MATH_FUNCTIONS
 #else
 #if HINTLIB_SIN_FOR_LONG_DOUBLE_NUM == 2
   #undef HINTLIB_INCLUDE_GLOBAL_MATH_FUNCTIONS
   #define HINTLIB_INCLUDE_GLOBAL_MATH_FUNCTIONS
 #else
 #if HINTLIB_SIN_FOR_LONG_DOUBLE_NUM == 3
   inline long double sin   (const long double& x) { return ::std::sinl   (x); }
   inline long double cos   (const long double& x) { return ::std::cosl   (x); }
   inline long double sqrt  (const long double& x) { return ::std::sqrtl  (x); }
   inline long double log   (const long double& x) { return ::std::logl   (x); }
   inline long double exp   (const long double& x) { return ::std::expl   (x); }
   inline long double floor (const long double& x) { return ::std::floorl (x); }
   inline long double ceil  (const long double& x) { return ::std::ceill  (x); }
   inline long double modf  (const long double& x, long double* y)
      { return ::std::modfl (x, y); }
   inline long double pow  (const long double& x, const long double& y)
      { return ::std::powl (x, y); }
 #else
 #if HINTLIB_SIN_FOR_LONG_DOUBLE_NUM == 4
   inline long double sin   (const long double& x) { return ::std::sinl   (x); }
   inline long double cos   (const long double& x) { return ::std::cosl   (x); }
   inline long double sqrt  (const long double& x) { return ::std::sqrtl  (x); }
   inline long double log   (const long double& x) { return ::std::logl   (x); }
   inline long double exp   (const long double& x) { return ::std::expl   (x); }
   inline long double floor (const long double& x) { return ::std::floorl (x); }
   inline long double ceil  (const long double& x) { return ::std::ceill  (x); }
   inline long double modf  (const long double& x, long double* y)
      { return ::std::modfl (x, y); }
   inline long double pow  (const long double& x, const long double& y)
      { return ::std::powl (x, y); }
 #else
   // for HINTLIB_SIN_FOR_LONG_DOUBLE_NUM == 0, we do nothing
 #endif
 #endif
 #endif
 #endif
 #endif

 // double

 #if HINTLIB_SIN_FOR_DOUBLE_NUM == 1 || HINTLIB_SIN_FOR_DOUBLE_NUM == 3
   #undef HINTLIB_INCLUDE_STD_MATH_FUNCTIONS
   #define HINTLIB_INCLUDE_STD_MATH_FUNCTIONS
 #else
   #undef HINTLIB_INCLUDE_GLOBAL_MATH_FUNCTIONS
   #define HINTLIB_INCLUDE_GLOBAL_MATH_FUNCTIONS
 #endif

 // Import functions using "using"

 #ifdef HINTLIB_INCLUDE_STD_MATH_FUNCTIONS
   using ::std::sin;
   using ::std::cos;
   using ::std::sqrt;
   using ::std::log;
   using ::std::exp;
   using ::std::floor;
   using ::std::ceil;
   using ::std::modf;
   using ::std::pow;
   #undef HINTLIB_INCLUDE_STD_MATH_FUNCTIONS
 #endif

 #ifdef HINTLIB_INCLUDE_GLOBAL_MATH_FUNCTIONS
   using ::sin;
   using ::cos;
   using ::sqrt;
   using ::log;
   using ::exp;
   using ::floor;
   using ::ceil;
   using ::modf;
   using ::pow;
   #undef HINTLIB_INCLUDE_GLOBAL_MATH_FUNCTIONS
 #endif

 #define HINTLIB_MN ::HIntLib::

#else
 #if HINTLIB_MATH_FUN_STRATEGY == 2
  #define HINTLIB_MN ::std::
 #else
  #define HINTLIB_MN ::
 #endif
#endif


#if 0
/**
 *  Other math functions:  lgamma()
 */
 
inline double lgamma (const double& x)
{
   return HINTLIB_LGAMMA_FOR_DOUBLE (x);
}

inline float lgamma (const float& x)
{
   return
#if HINTLIB_LGAMMA_FOR_FLOAT_NUM > 0
      HINTLIB_LGAMMA_FOR_FLOAT (x);
#else
      float (HINTLIB_LGAMMA_FOR_DOUBLE (double (x)));
#endif
}

#ifdef HINTLIB_HAVE_LONG_DOUBLE
#if HINTLIB_LGAMMA_FOR_LONG_DOUBLE_NUM > 0
inline long double lgamma (const long double& x)
{
   return HINTLIB_LGAMMA_FOR_LONG_DOUBLE (x);
}
#endif
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
 *     - Polynomial2<>
 */

template<class T> inline T cube (T x)
{
   return sqr (x) * x;
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
 *  divAndRoundUp
 */

template<class T> inline T divAndRoundUp (T a, T b)
{
   return (a + b - 1) / b;
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


#if 0

// The following functions have been removed. They are not used in HIntLib
// anymore, but depend on lgamma(), which is not available on all C compilers
// (e.g. MSVC), even though it is demanded by C99.


/**
 *  lfact()
 *
 *  Returns the natural logarithm of  n!
 */

double lfact (unsigned n);


/**
 *  bico()
 *  lbico()
 */

inline
double
lbico (unsigned n, unsigned k)
{
   return lfact (n) - lfact (k) - lfact (n - k);
}

inline
double
bico (unsigned n, unsigned k)
{
   return HINTLIB_MN floor (.5 + HINTLIB_MN exp (lbico (n, k)));
}
#endif


/**
 *  radicalInverseFunction ()
 *  radicalInverseFunction2 ()
 */

real radicalInverseFunction  (Index, unsigned base)  HINTLIB_GNU_CONST;
real radicalInverseFunction2 (Index)  HINTLIB_GNU_CONST;

}  // namespace HIntLib

#endif
