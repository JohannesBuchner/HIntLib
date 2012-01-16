/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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
 *  fallback/limits.h
 *
 *  This is a crude replacement for the <limits> standard library header.
 *  It is used only on systems that do not provide <limits> (most noteworthy
 *  GCC versions < 3.0).
 *  
 *  The values are taken either from <float.h> (for floating point types) or
 *  from <limits.h> (for integral types).
 *  Information not available in these header files is set to (i386 influenced)
 *  default values.
 */

#ifndef HINTLIB_FALLBACK_LIMITS_H
#define HINTLIB_FALLBACK_LIMITS_H 1

#ifndef HINTLIB_DEFAULTS_H
#error "HIntlib/defaults.h has to be included before HIntLib/fallback_limits.h"
#endif

#ifdef HINTLIB_HAVE_LIMITS
#error "HIntLib/fallback_limits.h should not be included if limits is available"
#endif

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
// Implementation in hlmath.cpp
#endif

#include <float.h>
#include <limits.h>

#ifdef HINTLIB_HAVE_CMATH
  #include <cmath>
#else
  #include <math.h>
#endif

namespace std
{

/**
 * Rounding style
 */

enum float_round_style
{
   round_indeterminate       = -1,
   round_towards_zero        =  0,
   round_to_nearest          =  1,
   round_toward_infinity     =  2,
   round_toward_neg_infinity =  3
};

/**
 * Unspecialized numeric_limits<>
 */

template<class T>
class numeric_limits
{
public:
   static const bool is_specialized = false;

   static const int radix = 0;
   static const int digits = 0;
   static const int digits10 = 0;

   static const bool is_signed = false;
   static const bool is_integer = false;
   static const bool is_exact = false;

   inline static T min() throw()  { return T(); }
   inline static T max() throw()  { return T(); }

   inline static T epsilon() throw()  { return T(); }
   inline static T round_error() throw()  { return T(); }

   static const int min_exponent = 0;
   static const int max_exponent = 0;
   static const int min_exponent10 = 0;
   static const int max_exponent10 = 0;

   static const bool has_infinity = false;
   static const bool has_quiet_NaN = false;
   static const bool has_signaling_NaN = false;
   static const bool has_denorm = false;
   static const bool has_denorm_loss = false;
   inline static T infinity() throw()  { return T(); }
   inline static T quiet_NaN() throw()  { return T(); }
   inline static T signaling_NaN() throw()  { return T(); }
   inline static T denorm_min() throw()  { return T(); }

   static const bool is_iec559 = false;
   static const bool is_bounded = false;
   static const bool is_modulo = false;

   static const bool traps = false;
   static const bool tinyness_before = false;
   static const float_round_style round_style = round_towards_zero;
};


/*****************************************************************************
 *               Floating point numbers
 *****************************************************************************/


namespace Private
{

/**
 *  Type-independent data for floating point numbers
 *  (float, double, long double)
 */

struct real
{
   static const bool is_specialized = true;

   static const int radix = FLT_RADIX;
   static const float_round_style round_style =
#ifdef HINTLIB_FLT_ROUNDS_CONST
      float_round_style (FLT_ROUNDS);
#else
      round_indeterminate;
#endif
      

#ifdef HINTLIB_HAVE_REAL_INFINITY
      static const bool has_infinity = false;
#else
     static const bool has_infinity = true;
#endif

   static const bool has_quiet_NaN = true;

   static const bool is_exact  = false;
   static const bool is_signed = true;
   static const bool is_integer = false;
   static const bool is_bounded = true;
   static const bool is_modulo = false;
};

/**
 *  Type-dependent data for floating point numbers
 */

template<class T>
struct realT : public real
{
   #ifdef HINTLIB_HAVE_REAL_INFINITY
      inline static T infinity() throw()  { return 0.0; }
   #endif

   inline static T round_error() throw()  { return 0.5; }
   inline static T quiet_NaN() throw()  { return T(sqrt(T(-1.0))); }
};

}  // namespace Private


/**
 *  numeric_limits<> for floating point numbers
 */

template<>
class numeric_limits<float> : public Private::realT<float>
{
public:
   static const int digits = FLT_MANT_DIG;
   static const int digits10 = FLT_DIG;

   inline static float min() throw()  { return FLT_MIN; }
   inline static float max() throw()  { return FLT_MAX; }

#ifndef HINTLIB_HAVE_REAL_INFINITY
   inline static float infinity() throw()  { return FLT_MAX / FLT_MIN; }
#endif

   inline static float epsilon() throw()  { return FLT_EPSILON; }

   static const int min_exponent = FLT_MIN_EXP;
   static const int min_exponent10 = FLT_MIN_10_EXP;
   static const int max_exponent = FLT_MAX_EXP;
   static const int max_exponent10 = FLT_MAX_10_EXP;
};

template<>
class numeric_limits<double> : public Private::realT<double>
{
public:
   static const int digits = DBL_MANT_DIG;
   static const int digits10 = DBL_DIG;

   inline static double min() throw()  { return DBL_MIN; }
   inline static double max() throw()  { return DBL_MAX; }

#ifndef HINTLIB_HAVE_REAL_INFINITY
   inline static float infinity() throw()  { return DBL_MAX / DBL_MIN; }
#endif
   inline static double epsilon() throw()  { return DBL_EPSILON; }

   static const int min_exponent = DBL_MIN_EXP;
   static const int min_exponent10 = DBL_MIN_10_EXP;
   static const int max_exponent = DBL_MAX_EXP;
   static const int max_exponent10 = DBL_MAX_10_EXP;
};

#ifdef HINTLIB_HAVE_LONG_DOUBLE
template<>
class numeric_limits<long double> : public Private::realT<long double>
{
public:
   static const int digits = LDBL_MANT_DIG;
   static const int digits10 = LDBL_DIG;

   inline static long double min() throw()  { return LDBL_MIN; }
   inline static long double max() throw()  { return LDBL_MAX; }

#ifndef HINTLIB_HAVE_REAL_INFINITY
   inline static float infinity() throw()  { return LDBL_MAX / LDBL_MIN; }
#endif

   inline static long double epsilon() throw()  { return LDBL_EPSILON; }

   static const int min_exponent = LDBL_MIN_EXP;
   static const int min_exponent10 = LDBL_MIN_10_EXP;
   static const int max_exponent = LDBL_MAX_EXP;
   static const int max_exponent10 = LDBL_MAX_10_EXP;
};
#endif


/*****************************************************************************
 *                 Integral data types
 *****************************************************************************/


namespace Private
{

   struct integer
   {
      static const bool is_specialized = true;

      static const int radix = 2;

      static const bool has_infinity = false;
      static const bool has_quiet_NaN = false;

      static const bool is_exact = true;
      static const bool is_integer = true;
      static const bool is_bounded = true;
      static const bool is_modulo = true;

      static const int min_exponent = 0;
      static const int max_exponent = 0;
      static const int min_exponent10 = 0;
      static const int max_exponent10 = 0;

      static const float_round_style round_style = round_towards_zero;
   };

   template<class T>
   struct integerT : public integer
   {
      inline static T infinity() throw()  { return 0; }
      inline static T quiet_NaN() throw()  { return 0; }
      inline static T epsilon() throw()  { return 0; }
   };

   struct u_integer
   {
      static const bool is_signed = false;
   };

   template<class T>
   struct u_integerT : public integerT<T>, public u_integer
   {
      static const int digits = sizeof(T) * CHAR_BIT;
      static const int digits10
         = int (digits * 0.3010299956639811952137388947245);
      inline static T min() throw()  { return T(0); }
   };

} // namespace Private

template<> class numeric_limits<unsigned char>
   : public Private::u_integerT<unsigned char>
{
public:
   inline static unsigned char max() throw()  { return UCHAR_MAX; }
};
#if CHAR_MIN == 0
template<> class numeric_limits<char> : public Private::u_integerT<char>
{
public:
   inline static char max() throw()  { return CHAR_MAX; }
};
#endif
template<> class numeric_limits<unsigned short>
   : public Private::u_integerT<unsigned short>
{
public:
   inline static unsigned short max() throw()  { return USHRT_MAX; }
};
template<> class numeric_limits<unsigned int>
   : public Private::u_integerT<unsigned int>
{
public:
   inline static unsigned int max() throw()  { return UINT_MAX; }
};
template<> class numeric_limits<unsigned long>
   : public Private::u_integerT<unsigned long>
{
public:
   inline static unsigned long max() throw()  { return ULONG_MAX; }
};
#ifdef HINTLIB_HAVE_UNSIGNED_LONG_LONG_INT
template<> class numeric_limits<unsigned long long>
   : public Private::u_integerT<unsigned long long>
{
public:
   inline static unsigned long long max() throw()  { return
   #ifdef ULLONG_MIN
      ULLONG_MIN;           // ISO C 99 standard
   #else
     #ifdef ULONG_LONG_MIN
        ULONG_LONG_MIN;     // GCC
     #else
	static_cast<unsigned long long>(-1);
     #endif
   #endif
   }
};
#endif


namespace Private
{
   struct s_integer
   {
      static const bool is_signed = true;
   };

   template <class T>
   struct s_integerT : public integerT<T>, public s_integer
   {
      static const int digits = sizeof(T) * CHAR_BIT - 1;
      static const int digits10
         = int (digits * 0.3010299956639811952137388947245);
   };
}  // namespace Private


template<> class numeric_limits<signed char>
   : public Private::s_integerT<signed char>
{
public:
   inline static char min() throw()  { return SCHAR_MIN; }
   inline static char max() throw()  { return SCHAR_MAX; }
};
#if CHAR_MIN < 0
template<> class numeric_limits<char> : public Private::s_integerT<char>
{
public:
   inline static char min() throw()  { return CHAR_MIN; }
   inline static char max() throw()  { return CHAR_MAX; }
};
#else
  #error "Invalid macro CHAR_MIN: "##CHAR_MIN
#endif
template<> class numeric_limits<short> : public Private::s_integerT<short>
{
public:
   inline static short min() throw()  { return SHRT_MIN; }
   inline static short max() throw()  { return SHRT_MAX; }
};
template<> class numeric_limits<int> : public Private::s_integerT<int>
{
public:
   inline static int min() throw()  { return INT_MIN; }
   inline static int max() throw()  { return INT_MAX; }
};
template<> class numeric_limits<long> : public Private::s_integerT<long>
{
public:
   inline static long min() throw()  { return LONG_MIN; }
   inline static long max() throw()  { return LONG_MAX; }
};

#ifdef HINTLIB_HAVE_LONG_LONG_INT
template<> class numeric_limits<long long>
   : public Private::s_integerT<long long>
{
public:
   inline static long long min() throw()  { return
   #ifdef LLONG_MIN
      LLONG_MIN;         // ISO C 99 standard
   #else
     #ifdef LONG_LONG_MIN
        LONG_LONG_MIN;   // GCC
     #else
	(static_cast<unsigned long long>(-1) >> 1) + 1;
     #endif
   #endif
   }
   inline static long long max() throw()  { return
   #ifdef LLONG_MAX
      LLONG_MAX;         // ISO C 99 standard
   #else
     #ifdef LONG_LONG_MAX
        LONG_LONG_MAX;   // GCC
     #else
	(static_cast<unsigned long long>(-1) >> 1);
     #endif
   #endif
   }
};
#endif

} // namespace std

#endif

