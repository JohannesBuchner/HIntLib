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
 *  bitop.h
 *
 *  Some functions that perform bit operations on integer-like types
 */

#ifndef HINTLIB_BITOP_H
#define HINTLIB_BITOP_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/defaults.h>

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
 *  Converts an integer to its GRAY Code
 */

template<class T>
inline T grayCode (T n)
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

template<class T>
inline T bit (T n, unsigned k)
{
   return n & (T(1) << k);
}


/**
 *  Selects the k lower order bits of n
 *
 *  leastSignificantBits (0)  returns 0
 */

template<class T>
inline T leastSignificantBits (T n, unsigned k)
{
   // the rhs operand to "<<" has to be strictly less the the number of digits

   return (k >= unsigned(std::numeric_limits<T>::digits))
      ? n : n & ((T(1) << k) - 1);
}


/**
 *  ls0()
 *
 *  Determines the position of the Least Significatn Zero (LSZ) of an integer
 *  number.
 *
 *  Example:
 *
 *     ls0 (xxxxxxx0) = 0
 *     ls0 (xxxxxx01) = 1
 *      :
 *     ls0 (01111111) = 7
 *     ls0 (11111111) = 8
 */

template<class T>
inline unsigned ls0 (T n)
{
   unsigned result = 0;

   while ((n & T(0xff)) == T(0xff))
   {
      n >>= 8;
      result += 8;
   }

   return result + Private::ls0_data [n & T(0xff)];
}


/**
 *  ls1()
 *
 *  Determines the position of the Least Significant One of an integer number.
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

template<class T>
inline int ls1 (T n)
{
   return n ? ls0 (n-1) : -1;
}


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

template<class T>
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


/**
 *  thread()
 *
 *  Given integers  ...a3a2a1a0 and ...b3b2b1b0, create ...b3a3b2a2b1a1b0a0
 */

template<class T>
inline T thread (T a, T b)
{
   T mask = 1;
   T result = 0;

   while (a | b)
   {
      if (a & 1)  result |= mask;
      mask <<= 1;
      a >>= 1;
      if (b & 1)  result |= mask;
      mask <<= 1;
      b >>= 1;
   }

   return result;
}


/**
 *  unthread()
 *
 *  Given an integer  ...a5a4a3a2a1a0, create ...a4a2a0 and ...a5a3a1.
 */

template<class T>
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

template<class T>
T thread (T a, T b, unsigned na, unsigned nb);

template<class T>
void unthread (T both, T& a, T& b, unsigned na, unsigned nb);

template<class T>
T threadn (T* indices, unsigned num);

template<class T>
void unthreadn (T all, T* indices, unsigned num);

template<class T>
T threadinf (T* indices, unsigned num);

template<class T>
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

