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
 *   My algorithm
 *
 *   A number of general purpose routines.  Some of them are along the lines
 *   of the STL.
 */

#ifndef HINTLIB_HLALGORITHM_H
#define HINTLIB_HLALGORITHM_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <algorithm>

#include <HIntLib/bitop.h>


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
 *  purgeArray()
 *
 *  Calls delete[]() on each member of the container
 */

template<typename In>
inline
void purgeArray (In first, In last)
{
   while (first != last)  delete[] *first++;
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
 *  min4 (x1,x2,x3,x4)
 *  max4 (x1,x3,x3,x4)
 */

template<typename T>
inline
const T& max4 (const T& x1, const T& x2, const T& x3, const T& x4)
{
   if (x1 > x2)
   {
      return max3 (x1, x3, x4);
   }
   else
   {
      return max3 (x2, x3, x4);
   }
}

template<typename T>
inline
const T& min4 (const T& x1, const T& x2, const T& x3, const T& x4)
{
   if (x1 < x2)
   {
      return min3 (x1, x3, x4);
   }
   else
   {
      return min3 (x2, x3, x4);
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

namespace Private
{
   void throwPartitionNotPossible() HINTLIB_GNU_NORETURN;
}

template<class Bi>
void initial_partition
   (Bi first, Bi last, typename std::iterator_traits<Bi>::value_type num)
{
   if (first == last)
   {
      if (num)  Private::throwPartitionNotPossible();
      return;
   }

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
 *  If no new partition can be created, false is returned.
 */

template<class Bi>
void initial_partition
   (Bi first, Bi last, typename std::iterator_traits<Bi>::value_type max,
    typename std::iterator_traits<Bi>::value_type num)
{
   if (first == last)
   {
      if (num)  Private::throwPartitionNotPossible();
      return;
   }

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
 *  next_combination()
 *  initial_combination()
 *
 *  Generates the next combination of  k  from  n  elements.
 *
 *  If no new combination can be created, false is returned.
 *
 *  next_combination() searches for the first element that can be moved to
 *  the next position (i.e, we search for the minimal  i  such that there is an
 *  element on position  i  and position  i + 1  is empty).
 *  This element is moved to the next position and all previous elements are
 *  returned to the leftmost positions.
 *
 *  There are two versions:
 *
 *    a) operating on a bool-array
 *    b) operating on the bits of an integer
 */

void initial_combination (bool* first, bool* last, int num);

template<typename T>
void initial_combination (T* pattern, int num)
{
   *pattern = lsBitsMask<T>(num);
}

bool next_combination (bool* first, bool* last);

template<typename T>
bool next_combination (T* pattern, int n)
{
   // make sure there are at least two stacks and some elements

   if (n <= 1 || *pattern == 0)  return false;
   const T last = T(1) << (n - 1);

   // Skip leading empty stacks

   T i = *pattern & ~(*pattern - 1);

   // track end-position for redistribution

   T end = 1;

   // Now we have found an element on position i.  We need to go on and find an
   // empty position

   for(;;)
   {
      // Check for (0,...,0, 1,...,1), i.e., the end

      if (i == last)  return false;
      i <<= 1;
      if ((*pattern & i) == 0)  break;
      end <<= 1;
   }

   // Now element  i-1  can be moved to position  i  and everything before that
   // can be moved back to the leftmost positions.

   *pattern = (*pattern & ~(i - 1)) | i | (end - 1);

   return true;
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

