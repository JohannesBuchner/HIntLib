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

#ifndef HINTLIB_OLD_BITOP_H
#define HINTLIB_OLD_BITOP_H 1

#ifdef __GNUG__
#pragma interface
#endif


namespace HIntLib
{
namespace NormalBitOp
{

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

   while (n & 1)
   {
      n >>= 1;
      result++;
   }

   return result;
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
 */

template<class T>
inline int ls1 (T n)
{
   if (! n) return -1;

   int result = 0;
 
   while (! (n & 1))
   {
      n >>= 1;
      result++;
   }
 
   return result; 
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
   int result = -1;
 
   while (n)
   {
      n >>= 1;
      result++;
   }
 
   return result;
}

}  // namespace NormalBitOp
}  // namespace HIntLib

#endif

