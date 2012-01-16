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
 *   My algorithm
 *
 *   A number of general purpose routines.  Some of them are along the lines
 *   of the STL.
 */

#ifndef HINTLIB_MYALGORITHM_H
#define HINTLIB_MYALGORITHM_H 1

#include <algorithm>


namespace HIntLib
{

/**
 *  Function deleting an object
 *
 *  Used by purge()
 */

template<typename T>
// inline
T* delete_ptr (T* p) { delete p; return 0; }


/**
 *  purge()
 *
 *  Calls delete() on each member of the container
 */

template<typename In>
inline
void purge (In first, In last)
{
   while (first != last)  delete *first++;

   // for_each (first, last, delete_ptr);     // does not compile :-((
}


/**
 *  min3 (x1,x2,x3)
 *  max3 (x1,x3,x3)
 */

template<typename T>
inline
const T& max3 (const T& x1, const T& x2, const T& x3)
{
   if (x1 > x2)
   {
      return x1 > x3 ? x1 : x3;
   }
   else
   {
      return x2 > x3 ? x2 : x3;
   }
}

template<typename T>
inline
const T& min3 (const T& x1, const T& x2, const T& x3)
{
   if (x1 < x2)
   {
      return x1 < x3 ? x1 : x3;
   }
   else
   {
      return x2 < x3 ? x2 : x3;
   }
}


/**
 *  initial_partition()
 *  next_partition ()
 *
 *  Generates the next partion.
 *
 *  We search for the first stack that is not empty.
 *  All elements but one return to the first stack. The remaining elements move
 *  on to the next stack.
 *
 *  If no new partition can be created, false is returned.
 */

template<class Bi>
void initial_partition
   (Bi first, Bi last, typename std::iterator_traits<Bi>::value_type num)
{
   *first = num;
   std::fill (first + 1, last, 0);
}

template<class Bi>
bool next_partition (Bi first, Bi last)
{
   // make sure there are at least two stacks

   if (first + 1 >= last)  return false;
   
   // Check if the first stack is not empty
   // This (frequent) case is treated seperately

   if (*first > 0)
   {
      -- *first;      // move one element to the second stack
      ++ *(first+1);
      return true;
   }

   // search for the first non-empty stack

   for (Bi i = first + 1; i < last-1; ++i)
   {
      if (*i > 0)   // is the current stack non-empty?
      {
         *first = *i - 1;   // all but one elment go back to stack # 0
         *i = 0;

         ++ *(i+1);  // the remaining element moves on to the next stack

         return true;
      }
   }

   return false;
}


/**
 *  initial_partition()
 *  next_partition ()
 *
 *  Generates the next partion, restricting the size of the subsets
 *
 *  We search for the first stack that is not empty.
 *  All elements but one return to the first stack. The remaining elements move
 *  on to the next stack.
 *
 *  If no new partition can be created, false is returned.
 */

template<class Bi>
void initial_partition
   (Bi first, Bi last, typename std::iterator_traits<Bi>::value_type max,
    typename std::iterator_traits<Bi>::value_type num)
{
   Bi i = first;

   while (num > max)
   {
      *(i++) = max;
      num -= max;
   }

   *(i++) = num;

   std::fill (i, last, 0);
}

template<class Bi>
bool next_partition
   (Bi first, Bi last, typename std::iterator_traits<Bi>::value_type max)
{
   // make sure there are at least two stacks

   if (first + 1 >= last)  return false;

   // Skip leading full stacks

   if (*first == max)
   {
      for (Bi i = first + 1; i < last - 1; ++i)
      {
         if (*i < max) break;
         ++first;
      }
   }
   
   // search for the first non-empty stack

   for (Bi i = first; i < last-1; ++i)
   {
      if (*i > 0)   // is the current stack non-empty?
      {
         // track end-position for redistribution

         Bi end = first;

         // search for the next stack with room for an element

         for (Bi j = i + 1; j < last; ++j)
         {
            if (*j < max)
            {
               ++ *j;   // one elment goes to this stack

               // all other elements back to the beginning

               *end = *i - 1;
               std::fill (first, end, max);
               std::fill (end + 1, j, 0);

               return true;
            }
            else
            {
               ++end;
            }
         }

         return false;
      }
   }

   return false;
}


/**
 *  Returns the difference between two unsigned integers of arbitrary size
 */

#define HINTLIB_MYALGO_H_XXX(type) \
inline \
type diff (unsigned type a, unsigned type b) \
{ \
   return static_cast<type>(a) - static_cast<type>(b); \
}

HINTLIB_MYALGO_H_XXX(char)
HINTLIB_MYALGO_H_XXX(int)
HINTLIB_MYALGO_H_XXX(long)
#ifdef HINTLIB_HAVE_LONG_LONG_INT
#ifdef HINTLIB_HAVE_UNSIGNED_LONG_LONG_INT
HINTLIB_MYALGO_H_XXX(long long)
#endif
#endif

#undef HINTLIB_MYALGO_H_XXX

}  // namespace HIntLib

#endif

