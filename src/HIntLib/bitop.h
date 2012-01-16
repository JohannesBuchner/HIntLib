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
 *  bitop.h
 *
 *  Some functions that perform bit operations on integer-like types
 */

#ifndef HINTLIB_BITOP_H
#define HINTLIB_BITOP_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#ifdef HINTLIB_HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/fallback_limits.h>
#endif

namespace HIntLib
{

/**
 *  This data is private
 *
 *  Do not used directly
 */

namespace Private
{
   extern HINTLIB_DLL_IMPORT const int ms1_data [];
   extern HINTLIB_DLL_IMPORT const int ls0_data [];
}


/**
 *  grayCode()
 *
 *  Converts an integer  x  to its GRAY Code  G(x)
 *
 *  G(0) = 0
 *  G(1) = 1
 *  G(2) = 3
 *  G(3) = 2
 *  G(4) = 6
 *  G(5) = 7
 *       :
 *
 *  G(x) = x iff x in {0,1}
 *
 *  G({0,...,2^k - 1}) = {0,...,2^k - 1}  for all k = 0,1,2,...
 *
 *  popCount (G(x) ^ G(x+1)) = 1  for all x
 */

template<typename T>
inline
T grayCode (const T& n)
{
   return n ^ (n >> 1);
}


/**
 *  bit()
 *
 *  returns bit number k from n
 *
 *  k must be from the range 0 .. numeric_limits<T>::digits - 1
 */

template<typename T>
inline
T bit (const T& n, int k)
{
   return n & (T(1) << k);
}


/**
 *  lsBitsMask()
 *
 *  Returns a mask selecting the  k  least significant bits.
 */

template<typename T>
inline
T lsBitsMask (int k)
{
   // the rhs operand to ">>" has to be strictly less the the number of digits

   return (k == std::numeric_limits<T>::digits) ? ~T(0) : (T(1) << k) - 1;
}
   


/**
 *  lsBits()
 *
 *  Selects the k lower order bits of n
 *
 *  leastSignificantBits (0)  returns 0
 */

template<typename T>
inline
T lsBits (const T& n, int k)
{
   return n & lsBitsMask<T>(k);
}


/**
 *  ls0()
 *
 *  Determines the position (starting with zero) of the Least Significant Zero
 *  (LSZ) of an integer number.
 *  This is the same as the number of trailing Ones.
 *
 *  Example:
 *
 *     ls0 (xxxxxxx0) = 0
 *     ls0 (xxxxxx01) = 1
 *                    :
 *     ls0 (01111111) = 7
 *     ls0 (11111111) = 8
 */

template<typename T>
inline
int ls0 (T n)
{
   int result = 0;

   while ((n & T(0xff)) == T(0xff))
   {
      n >>= 8;
      result += 8;
   }

   return result + Private::ls0_data [n & T(0xff)];
}


#ifdef HINTLIB_HAVE_BUILTIN_CTZ
inline int ls0 (unsigned n)
{
   n = ~n;
   return n ? __builtin_ctz(n) : std::numeric_limits<unsigned>::digits;
}
inline int ls0 (         int n)   { return ls0 (static_cast<unsigned>(n)); }
inline int ls0 (         short n) { return ls0 (static_cast<unsigned>(n)); }
inline int ls0 (unsigned short n) { return ls0 (static_cast<unsigned>(n)); }
inline int ls0 (unsigned char  n) { return ls0 (static_cast<unsigned>(n)); }
inline int ls0 (signed   char  n) { return ls0 (static_cast<unsigned>(n)); }
inline int ls0 (         char  n) { return ls0 (static_cast<unsigned>(n)); }
#endif
#ifdef HINTLIB_HAVE_BUILTIN_CTZL
inline int ls0 (unsigned long n)
{
   n = ~n;
   return n ? __builtin_ctzl(n) : std::numeric_limits<unsigned long>::digits;
}
inline int ls0 (long n) { return ls0 (static_cast<unsigned long>(n)); }
#endif
#ifdef HINTLIB_HAVE_BUILTIN_CTZLL
#ifdef HINTLIB_HAVE_LONG_LONG_INT
#ifdef HINTLIB_HAVE_UNSIGNED_LONG_LONG_INT
inline int ls0 (unsigned long long n)
{
   n = ~n;
   return n ? __builtin_ctzll(n)
            : std::numeric_limits<unsigned long long>::digits;
}
inline int ls0 (long long n) { return ls0(static_cast<unsigned long long>(n)); }
#endif
#endif
#endif


/**
 *  ls1()
 *
 *  Determines the position (starting with zero) of the Least Significant One
 *  of an integer number.
 *  This is the same as the number of trailing zeros.
 *
 *  Example:
 *
 *     ls1 (00000000) = -1
 *     ls1 (xxxxxxx1) = 0
 *     ls1 (xxxxxx10) = 1
 *      :
 *     ls1 (x1000000) = 6
 *     ls1 (10000000) = 7
 *
 *  The value of  ls1(n)  is calculated as  ls0(n-1) .
 */

template<typename T>
inline int ls1 (T n)
{
   return n ? ls0 (n-1) : -1;
}

#ifdef HINTLIB_HAVE_BUILTIN_FFS
inline int ls1 (unsigned n)
{
   return __builtin_ffs(n) - 1;
}
inline int ls1 (         int n)   { return ls1 (static_cast<unsigned>(n)); }
inline int ls1 (         short n) { return ls1 (static_cast<unsigned>(n)); }
inline int ls1 (unsigned short n) { return ls1 (static_cast<unsigned>(n)); }
inline int ls1 (unsigned char  n) { return ls1 (static_cast<unsigned>(n)); }
inline int ls1 (signed   char  n) { return ls1 (static_cast<unsigned>(n)); }
inline int ls1 (         char  n) { return ls1 (static_cast<unsigned>(n)); }
#endif
#ifdef HINTLIB_HAVE_BUILTIN_FFSL
inline int ls1 (unsigned long n)
{
   return __builtin_ffsl(n) - 1;
}
inline int ls1 (long n) { return ls1 (static_cast<unsigned long>(n)); }
#endif
#ifdef HINTLIB_HAVE_BUILTIN_FFSLL
#ifdef HINTLIB_HAVE_LONG_LONG_INT
#ifdef HINTLIB_HAVE_UNSIGNED_LONG_LONG_INT
inline int ls1 (unsigned long long n)
{
   return __builtin_ffsll(n) - 1;
}
inline int ls1 (long long n) { return ls1(static_cast<unsigned long long>(n)); }
#endif
#endif
#endif


/**
 *  ms1()
 *
 *  Determines the position of the Most Significant One of an integer number.
 *
 *  Examples:
 *
 *     ms1 (00000000) = -1
 *     ms1 (00000001) = 0
 *     ms1 (0000001x) = 1
 *      :
 *     ms1 (01xxxxxx) = 6
 *     ms1 (1xxxxxxx) = 7
 */

template<typename T>
inline int ms1 (T n)
{
   int result = 0;

   while ((~T() ^ T(0xff)) & n)
   {
      n >>= 8;
      result += 8;
   }

   return result + Private::ms1_data [n];
}


#ifdef HINTLIB_HAVE_BUILTIN_CLZ
inline int ms1 (unsigned n)
{
   return n ? std::numeric_limits<unsigned>::digits - 1 - __builtin_clz(n) : -1;
}
inline int ms1 (         int n)   { return ms1 (static_cast<unsigned>(n)); }
inline int ms1 (         short n) { return ms1 (static_cast<unsigned>(n)); }
inline int ms1 (unsigned short n) { return ms1 (static_cast<unsigned>(n)); }
inline int ms1 (unsigned char  n) { return ms1 (static_cast<unsigned>(n)); }
inline int ms1 (signed   char  n) { return ms1 (static_cast<unsigned>(n)); }
inline int ms1 (         char  n) { return ms1 (static_cast<unsigned>(n)); }
#endif
#ifdef HINTLIB_HAVE_BUILTIN_CLZL
inline int ms1 (unsigned long n)
{
   return n ? std::numeric_limits<unsigned long>::digits - 1
                - __builtin_clzl(n) : -1;
}
inline int ms1 (long n) { return ms1 (static_cast<unsigned long>(n)); }
#endif
#ifdef HINTLIB_HAVE_BUILTIN_CLZLL
#ifdef HINTLIB_HAVE_LONG_LONG_INT
#ifdef HINTLIB_HAVE_UNSIGNED_LONG_LONG_INT
inline int ms1 (unsigned long long n)
{
   return n ? std::numeric_limits<unsigned long long>::digits - 1
                - __builtin_clzll(n) : -1;
}
inline int ms1 (long long n) { return ms1(static_cast<unsigned long long>(n)); }
#endif
#endif
#endif


/**
 *  popCount()
 *
 *  Determines the number of 1s.
 */

template<typename T>
inline int popCount(T x)
{
   int result = 0;

   while (x)
   {
      ++result;
      x &= x - 1;
   }

   return result;
}

#ifdef HINTLIB_HAVE_BUILTIN_POPCOUNT
inline int popCount (unsigned     n)   { return __builtin_popcount(n); }
inline int popCount (         int n)   { return __builtin_popcount(n); }
inline int popCount (         short n) { return __builtin_popcount(n); }
inline int popCount (unsigned short n) { return __builtin_popcount(n); }
inline int popCount (unsigned char  n) { return __builtin_popcount(n); }
inline int popCount (signed   char  n) { return __builtin_popcount(n); }
inline int popCount (         char  n) { return __builtin_popcount(n); }
#endif
#ifdef HINTLIB_HAVE_BUILTIN_POPCOUNTL
inline int popCount (unsigned long n) { return __builtin_popcountl(n); }
inline int popCount (         long n) { return __builtin_popcountl(n); }
#endif
#ifdef HINTLIB_HAVE_BUILTIN_POPCOUNTLL
#ifdef HINTLIB_HAVE_LONG_LONG_INT
#ifdef HINTLIB_HAVE_UNSIGNED_LONG_LONG_INT
inline int popCount (unsigned long long n) { return __builtin_popcountll(n); }
inline int popCount (         long long n) { return __builtin_popcountll(n); }
#endif
#endif
#endif


/**
 *  parity()
 *
 *  Returns 0 if the number of 1s is even
 *  Returns 1 if the number of 1s is odd
 */

inline
int parity (unsigned char x)
{
   x ^= x >> 4;
   x &= 0xf;
   return (0x6996 >> x) & 1;
}

inline int parity (unsigned x)
{
#ifdef HINTLIB_HAVE_BUILTIN_PARITY
   return __builtin_parity(x);
#else
#if HINTLIB_SIZEOF_UNSIGNED_INT > 16
# error "unsigned too large"
#endif
#if HINTLIB_SIZEOF_UNSIGNED_INT > 8
   x ^= x >> 64;
#endif
#if HINTLIB_SIZEOF_UNSIGNED_INT > 4
   x ^= x >> 32;
#endif
#if HINTLIB_SIZEOF_UNSIGNED_INT > 2
   x ^= x >> 16;
#endif
#if HINTLIB_SIZEOF_UNSIGNED_INT > 1
   x ^= x >> 8;
#endif
   x ^= x >> 4;
   x &= 0xf;
   return (0x6996 >> x) & 1;
#endif
}

inline int parity (unsigned short x) { return parity(unsigned (x)); }

inline int parity (unsigned long x)
{
#ifdef HINTLIB_HAVE_BUILTIN_PARITYL
   return __builtin_parityl(x);
#else
#if HINTLIB_SIZEOF_UNSIGNED_LONG_INT > 16
# error "unsigned long too large"
#endif
#if HINTLIB_SIZEOF_UNSIGNED_LONG_INT > 8
   x ^= x >> 64;
#endif
#if HINTLIB_SIZEOF_UNSIGNED_LONG_INT > 4
   x ^= x >> 32;
#endif
#if HINTLIB_SIZEOF_UNSIGNED_LONG_INT > 2
   x ^= x >> 16;
#endif
#if HINTLIB_SIZEOF_UNSIGNED_LONG_INT > 1
   x ^= x >> 8;
#endif
   x ^= x >> 4;
   x &= 0xf;
   return (0x6996 >> x) & 1;
#endif
}

#ifdef HINTLIB_HAVE_UNSIGNED_LONG_LONG_INT
inline int parity (unsigned long long x)
{
#ifdef HINTLIB_HAVE_BUILTIN_PARITYLL
   return __builtin_parityll(x);
#else
#if HINTLIB_SIZEOF_UNSIGNED_LONG_LONG_INT > 16
# error "unsigned long long too large"
#endif
#if HINTLIB_SIZEOF_UNSIGNED_LONG_LONG_INT > 8
   x ^= x >> 64;
#endif
#if HINTLIB_SIZEOF_UNSIGNED_LONG_LONG_INT > 4
   x ^= x >> 32;
#endif
#if HINTLIB_SIZEOF_UNSIGNED_LONG_LONG_INT > 2
   x ^= x >> 16;
#endif
#if HINTLIB_SIZEOF_UNSIGNED_LONG_LONG_INT > 1
   x ^= x >> 8;
#endif
   x ^= x >> 4;
   x &= 0xf;
   return (0x6996 >> x) & 1;
#endif
}
#endif


/**
 *  thread()
 *
 *  Given integers  ...a3a2a1a0 and ...b3b2b1b0, create ...b3a3b2a2b1a1b0a0
 */

template<typename T>
inline T thread (T a, T b)
{
   T mask = 1;
   T result = 0;

   while (a)
   {
      if (a & 1)  result |= mask;
      mask <<= 2;
      a >>= 1;
   }

   mask = 2;

   while (b)
   {
      if (b & 1)  result |= mask;
      mask <<= 2;
      b >>= 1;
   }

   return result;
}


/**
 *  unthread()
 *
 *  Given an integer  ...a5a4a3a2a1a0, create ...a4a2a0 and ...a5a3a1.
 */

template<typename T>
inline void unthread (T both, T& a, T& b)
{
   a = b = 0;
   T mask = 1;

   while (both)
   {
      if (both & 1)  a |= mask;
      if (both & 2)  b |= mask;
      both >>= 2;
      mask <<= 1;
   }
}


/**
 *  thread()
 *  unthread()
 *
 *  threan()
 *  unthreadn()
 */

template<typename T>
T thread (T a, T b, unsigned na, unsigned nb);

template<typename T>
void unthread (T both, T& a, T& b, unsigned na, unsigned nb);

template<typename T>
T threadn (T* indices, unsigned num);

template<typename T>
void unthreadn (T all, T* indices, unsigned num);

template<typename T>
T threadinf (T* indices, unsigned num);

template<typename T>
unsigned unthreadinf (T all, T* indices);


/**
 *  BitRef
 */

template<typename T>
class BitRef
{
public:
   BitRef (T* _ptr, unsigned bit) : ptr (_ptr), mask (T(1) << bit) {}

   operator unsigned char() const { return (*ptr & mask) != 0; }
   BitRef<T>&  operator= (unsigned char x)
      { if (x) *ptr |= mask; else *ptr &= ~mask; return *this; }
private:
   T* const ptr;
   const T  mask;
};

}  // namespace HIntLib

#endif

