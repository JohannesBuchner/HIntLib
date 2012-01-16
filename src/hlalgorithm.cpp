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

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/hlalgorithm.h>
#include <HIntLib/exception.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

namespace L = HIntLib;
namespace P = HIntLib::Private;


/**
 *  initial_combination()
 *
 *  Initialize the selection to  (X,...,X,-,...,-).
 */

void
L::initial_combination (bool* first, bool* last, int num)
{
   std::fill (first, first + num, true);
   std::fill (first + num, last, false);
}


/**
 *  next_combination()
 *
 *  Searches for the first element that can be moved to the next position
 *  (i.e, we search for the minimal  i  such that there is an element on
 *  position  i  and position  i + 1  is empty).
 *  This element is moved to the next position and all previous elements are
 *  returned to the leftmost positions.
 */

bool
L::next_combination (bool* first, bool* last)
{
   // make sure there are at least two stacks

   if (first + 1 >= last)  return false;

   // Skip leading empty stacks

   bool* i = first;
   while (! *i)
   {
      // If there are no elements, we are done

      if (++i == last)  return false;
   }

   // track end-position for redistribution

   bool* end = first;

   // Now we have found an element on position i.  We need to go on and find an
   // empty position

   for(;;)
   {
      // Check for (0,...,0, 1,...,1), i.e., the end

      if (++i == last)  return false;
      if (! *i)  break;
      ++end;
   }

   // Now element  i-1  can be moved to position  i  and everything before that
   // can be moved back to the leftmost positions.

   *i = true;

   std::fill (first, end, true);
   std::fill (end, i, false);

   return true;
}


/**
 *  throwPartitionNotPossible()
 */

void
P::throwPartitionNotPossible()
{
   throw OtherException ("Partitions of a positive number on zero slots "
                         "do not exist!");
}


