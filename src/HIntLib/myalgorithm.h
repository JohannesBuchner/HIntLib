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

#ifndef MYALGORITHM_H
#define MYALGORITHM_H 1

#include <algorithm>


namespace HIntLib
{

/**
 *  Function deleting an object
 *
 *  Used by purge()
 */

template<class T>
// inline
T* delete_ptr (T* p) { delete p; return 0; }


/**
 *  purge()
 *
 *  Calls delete() on each member of the container
 */

template<class In>
inline
void purge (In first, In last)
{
   while (first != last)  delete *first++;

   // for_each (first, last, delete_ptr);     // does not compile :-((
}


/**
 *  next_partition ()
 *
 *  Generates the next partion.
 *
 *  We search for the first stack that is not empty.
 *  All elements but one return to first stack. The remaining elemnt moves on
 *  to the next stack.
 */

template<class Bi>
inline
bool next_partition (Bi first, Bi last)
{
   // make sure there are at least two stacks

   if (first + 1 >= last)  return false;
   
   // Check if the first stack is not empty
   // This special case is treated seperately

   if (*first > 0)
   {
      -- *first;      // move one element to the second stack
      ++ *(first+1);
      return true;
   }

   // search for the first non-empty stack

   for (Bi i = first + 1; i < last-1; ++i)
   {
      if (*i > 0)   // find postion that is not empty
      {
         *first = *i - 1;   // all but one elment goes back to stack # 0
         *i = 0;

         ++ *(i+1);  // the remaining element moves on to the next stack

         return true;
      };
   }

   return false;
}


/**
 *  Returns the difference between two unsigned integers of arbitrary size
 */

#define MYALGO_H_XXX(type) \
inline \
type diff (unsigned type a, unsigned type b) \
{ \
   return static_cast<type>(a) - static_cast<type>(b); \
}

MYALGO_H_XXX(char)
MYALGO_H_XXX(int)
MYALGO_H_XXX(long)
#ifdef HAVE_LONG_LONG_INT
MYALGO_H_XXX(long long)
#endif

#undef MYALGO_H_XXX

}  // namespace HIntLib

#endif

