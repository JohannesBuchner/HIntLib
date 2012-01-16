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
 *  bitop.h
 *
 *  Some functions that perform bit operations on integer-like types
 */

#ifndef BITOP_H
#define BITOP_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/defaults.h>

#ifdef HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/hintlib_limits.h>
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
   extern const int ms1_data [];
   extern const int ls0_data [];
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

}  // namespace HIntLib

#endif

