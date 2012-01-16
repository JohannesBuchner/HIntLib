dnl ###########################################################################
dnl HIntLib autoconf macro definitions
dnl ###########################################################################

# HL_UNREACHABLE_CALLS_REMOVED
# ----------------------------
# Check if unreachable calls are removed
AC_DEFUN([HL_UNREACHABLE_CALLS_REMOVED],
[AC_CACHE_CHECK([whether unreachable calls are removed],
hl_cv_unreachable_calls_removed,
[AC_LINK_IFELSE(AC_LANG_PROGRAM([[
#ifdef HAVE_LIMITS
#include <limits>
#endif

void missing();

template<typename T> class Y {};

template<> class Y<int>
{
public:
   static bool f ()  { return true; }
   static const int anInt = 5;
};

template<typename T,unsigned a> class X
{
   public: void test();
};

template<typename T, unsigned a>
void X<T,a>::test()
{
   if (1==2) missing();
   if (a==3) missing();
   if (! Y<T>::f())  missing();
   if (~T() >> Y<T>::anInt == 0) missing();
   if (a % 4 == 1)  missing();
#ifdef HAVE_LIMITS
   if (! std::numeric_limits<T>::is_signed)  missing();
   if (~T() >> (std::numeric_limits<T>::digits - 5) == 0) missing();
#endif
}
]],[[
X<int, 8> x;
x.test();
exit (0);
]]), hl_cv_unreachable_calls_removed=yes, hl_cv_unreachable_calls_removed=no)])

if test x"$hl_cv_unreachable_calls_removed" = xyes; then
   AC_DEFINE([UNREACHABLE_CALLS_REMOVED],1,
      [define if unreachable calls are removed by the optimizer])
fi
])


# HL_ABS_AVAILABLE(PREFIX)
# ------------------------
# Check if PREFIX::abs (dobule) is available
AC_DEFUN([HL_ABS_AVAILABLE],
[AC_CACHE_CHECK([whether $1::abs(double) is available], hl_cv_$2_abs_available,
[AC_RUN_IFELSE(AC_LANG_PROGRAM([[
#include <stdlib.h>
#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif
]],[[
   if (   $1::abs (static_cast<float> (-3.75)) == static_cast<float> (3.75)
       && $1::abs (static_cast<double> (-3.75)) == static_cast<double> (3.75)
       && $1::abs (static_cast<long double> (-3.75)) == static_cast<long double> (3.75))
   {
      exit (0);
   }

   exit (1);
]]),
hl_cv_$2_abs_available=yes,
hl_cv_$2_abs_available=no,
AC_MSG_ERROR([HL_ ABS_AVAILABLE can not be used when crosscompiling!]))])

if test x"$hl_cv_$2_abs_available" = xyes; then
   AC_DEFINE(translit([$2],[a-z],[A-Z])_ABS_AVAILABLE, 1,
             [define if $1::abs(double) is available])
fi
])


# HL_FLT_ROUNDS_CONST
# -------------------
# Check if FLT_ROUNDS is constant
AC_DEFUN([HL_FLT_ROUNDS_CONST],
[AC_CACHE_CHECK([whether FLT_ROUNDS is constant], hl_cv_flt_rounds_const,
[AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[
#include <float.h>
]],[[
#if FLT_ROUNDS < 7
  ;
#endif
  ;
]]),
hl_cv_flt_rounds_const=yes, hl_cv_flt_rounds_const=no)])

if test x"$hl_cv_flt_rounds_const" = xyes; then
   AC_DEFINE([FLT_ROUNDS_CONST],1,[define if FLT_ROUNDS is constant])
fi
])


# HL_IEEE_MAGIC_WORKS
# -------------------
# Check if we can use IEEE Magic
AC_DEFUN([HL_IEEE_MAGIC_WORKS],
[AC_CACHE_CHECK([whether IEEE magic works], hl_cv_ieee_magic_works,
[AC_RUN_IFELSE(AC_LANG_PROGRAM([[
#include <stdlib.h>
#if SIZEOF_UNSIGNED_INT >= 4
   #undef UNSIGNED_NOT_EQUAL_U32
   typedef unsigned u32;
#else
#if SIZEOF_UNSIGNED_LONG_INT >= 4
   #define UNSIGNED_NOT_EQUAL_U32 1
   typedef unsigned long u32;
#else
   error "unsigned long int does not have 32 bits!"
#endif
#endif

#if SIZEOF_UNSIGNED_INT >= 8
   #undef U32_NOT_EQUAL_U64
   typedef unsigned u64;
#else
#if SIZEOF_UNSIGNED_LONG_INT >= 8
   #ifdef UNSIGNED_NOT_EQUAL_U32
      #undef U32_NOT_EQUAL_U64
   #else
      #define U32_NOT_EQUAL_U64
   #endif
   typedef unsigned long u64;
#else
#ifdef HAVE_UNSIGNED_LONG_LONG_INT
   #if SIZEOF_UNSIGNED_LONG_LONG_INT >= 8
      #define U32_NOT_EQUAL_U64 1
      typedef unsigned long long u64;
   #else
      #error "Can not determine an unsigned integer type with at least 64 bits!"
   #endif
#else
   #error "Can not determine an unsigned integer type with at least 64 bits!"
#endif
#endif
#endif
]],[[
   union { double d; u64 i; } x64;

   x64.d = 1.0; if (x64.d != 1.0)  exit (1);
   x64.i |= (1ull << 51); if (x64.d != 1.5)   exit (1);
   x64.i |= (1ull << 50); if (x64.d != 1.75)  exit (1);
   x64.i ^= (1ull << 51); if (x64.d != 1.25)  exit (1);
   
   union { float d; u32 i; } x32;

   x32.d = 1.0f; if (x32.d != 1.0f)  exit (1);
   x32.i |= (1ul << 22); if (x32.d != 1.5f)   exit (1);
   x32.i |= (1ul << 21); if (x32.d != 1.75f)  exit (1);
   x32.i ^= (1ul << 22); if (x32.d != 1.25f)  exit (1);
   
   exit (0);
]]),
hl_cv_ieee_magic_works=yes, hl_cv_ieee_magic_works=no,
hl_cv_ieee_magic_works=no)])

if test x"$hl_cv_ieee_magic_works" = xyes; then
   AC_DEFINE([IEEE_MAGIC_WORKS],1,[define if IEEE magic works])
fi
])


# HL_STREAMS_SUPPORT_LOCAL
# ------------------------
# Check if streams support local
AC_DEFUN([HL_STREAMS_SUPPORT_LOCAL],
[
AC_CACHE_CHECK([whether streams support local], hl_cv_streams_support_local,
[AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[
#include <iostream>
]],[[
   std::cin.imbue (std::cout.getloc());
]]),
hl_cv_streams_support_local=yes, hl_cv_streams_support_local=no)])

if test x"$hl_cv_streams_support_local" = xyes; then
   AC_DEFINE(STREAMS_SUPPORT_LOCAL,1,[define if streams support local])
fi
])

# HL_PROG_F77_REALLY_WORKS
# ------------------------
#
# Deterimes if the Fortran 77 compile found by AC_PROG_F77 can comile a
# program without errors.
AC_DEFUN([HL_PROG_F77_REALLY_WORKS],
[
  AC_REQUIRE([AC_PROG_F77])dnl
  AC_CACHE_CHECK(
    [whether ${F77} really works], hl_cv_prog_f77_really_works,
    [
      AC_LANG_PUSH([Fortran 77])
      AC_COMPILE_IFELSE(AC_LANG_PROGRAM(,
        [
          INTEGER I
          I = 1
        ]),
        hl_cv_prog_f77_really_works=yes,
        hl_cv_prog_f77_really_works=no)
      AC_LANG_POP(Fortran 77)
    ])
])

dnl End of acinclude.m4

