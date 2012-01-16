###############################################################################
# Makro definitions
###############################################################################


# Checks abs (dobule)

AC_DEFUN([HL_ABS_AVAILABLE],
[AC_CACHE_CHECK([whether $1::abs(double) is available], hl_cv_$2_abs_available,
[AC_TRY_RUN([
#include <stdlib.h>
#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

int main()
{
   if (   $1::abs (static_cast<float> (-3.75)) == static_cast<float> (3.75)
       && $1::abs (static_cast<double> (-3.75)) == static_cast<double> (3.75)
       && $1::abs (static_cast<long double> (-3.75)) == static_cast<long double> (3.75))
   {
      exit (0);
   }

   exit (1);
}
], hl_cv_$2_abs_available=yes, hl_cv_$2_abs_available=no)])

if test x"$hl_cv_$2_abs_available" = xyes; then
   AC_DEFINE(translit([$2],[a-z],[A-Z])_ABS_AVAILABLE, 1,
             [define if $1::abs(double) is available])
fi
])

AC_DEFUN([HL_FLT_ROUNDS_CONST],
[AC_CACHE_CHECK([whether FLT_ROUNDS is constant], hl_cv_flt_rounds_const,
[AC_TRY_COMPILE([
#include <float.h>
],[
#if FLT_ROUNDS < 7
  ;
#endif
],hl_cv_flt_rounds_const=yes, hl_cv_flt_rounds_const=no)])

if test x"$hl_cv_flt_rounds_const" = xyes; then
   AC_DEFINE([FLT_ROUNDS_CONST],1,[define if FLT_ROUNDS is constant])
fi
])


# Check if we can use IEEE Magic

AC_DEFUN([HL_IEEE_MAGIC_WORKS],
[AC_CACHE_CHECK([whether IEEE magic works], hl_cv_ieee_magic_works,
[AC_TRY_RUN([
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
int main()
{
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
}
],hl_cv_ieee_magic_works=yes, hl_cv_ieee_magic_works=no, hl_cv_ieee_magic_works=no)])

if test x"$hl_cv_ieee_magic_works" = xyes; then
   AC_DEFINE([IEEE_MAGIC_WORKS],1,[define if IEEE magic works])
fi
])

AC_DEFUN([HL_STREAMS_SUPPORT_LOCAL],
[
AC_CACHE_CHECK([whether streams support local], hl_cv_streams_support_local,
[AC_TRY_COMPILE([
#include <iostream>
],[
   std::cin.imbue (std::cout.getloc());
],hl_cv_streams_support_local=yes, hl_cv_streams_support_local=no)])

if test x"$hl_cv_streams_support_local" = xyes; then
   AC_DEFINE(STREAMS_SUPPORT_LOCAL,1,[define if streams support local])
fi
])

