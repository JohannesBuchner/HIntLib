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

#include <iostream>

#include <HIntLib/hlalgorithm.h>
#include <HIntLib/hlmath.h>
#include <HIntLib/array.h>
#include <HIntLib/output.h>

#include "test.h"

using std::cout;
using std::cerr;
using std::endl;

const char options[] = "";
const char option_msg[] = "";
const char testProgramParameters[] = "slots items [OPTION]...";
const char testProgramUsage[] =
   "Lists all combinations of _items_ from _slots_ items.\n\n";
const char testProgramName[] = "test_combination";
const int  testProgramCopyright = 2006;


bool opt(int, const char*)
{
   return false;
}

void test(int argc, char** argv)
{
   if (argc != 2)  usage("Invalid number of arguments!");

   const int slots = atoi (argv [0]);
   const int items = atoi (argv [1]);
   
   if (slots < 0)  usage("Number of slots must be non-negative!");

   if (slots > 32)  usage("Sorry, but more than 32 slots are not allowed!");

   if (items < 0 || items > slots)
   {
      usage("Number of items must be between 0 and number of slots!");
   }

   NORMAL printHeader (cout);

   // Storage for the combinations

   L::Array<int> partition (slots);
   L::Array<bool> combination (slots);
   L::u32 combination32;

   // Initialize combinations

   L::initial_partition (&partition[0], &partition[slots], 1, items);
   L::initial_combination (&combination[0], &combination[slots], items);
   L::initial_combination (&combination32, items);

   int combinationCount = 0;
   L::u32 lastCode = 0;

   // Create all combinations

   for(;;)
   {
      // Print combination

      DEB1
      {
         cout << "   [";
         for (int i = 0; i < slots; ++i)
         {
            // WHITE CIRCLE, BLACK CIRCLE
            cout << Wgl4Ascii (
                  (combination[i] ? "\xe2\x97\x8f" : "\xe2\x97\x8b"),
                  (combination[i] ? "X" : "."));
         }
         cout << "]   " << combination32 << '\n';
      }

      // Make sure that there are exactly _slots_ entries set

      int count = 0;
      for (int i = 0; i < slots; ++i)  if (combination[i])  ++count;
      if (count != items)  error ("Wrong number of itmes in combination");

      // Convert partition result to number

      L::u32 partitionCode = 0;
      for (int i = slots - 1; i >= 0; --i)
      {
         const int x = partition[i];
         if (x != 0 && x != 1)  error ("Value in partition must be 0 or 1");
         partitionCode = (partitionCode << 1) | L::u32(x);
      }

      // Convert combination result to number

      L::u32 combinationCode = 0;
      for (int i = slots - 1; i >= 0; --i)
      {
         const bool x = combination[i];
         combinationCode = (combinationCode << 1) | L::u32(x);
      }

      // Compare all three codes

      if (partitionCode != combination32 || combinationCode != combination32)
      {
         cout << "Code from combination: " << combinationCode << "\n"
                 "Code from partition:   " << partitionCode   << "\n"
                 "Combination code:      " << combination32   << '\n';
         error ("Codes do not match");
      }

      // Make sure, the combinations are enumerated in increasing order

      if (lastCode != 0 && combination32 <= lastCode)
      {
         error ("Combination code did not increase");
      }
      lastCode = combination32;

      // Count this combination

      ++combinationCount;

      // Get next combination

      const bool nextPartition
         = L::next_partition   (&partition[0], &partition[slots], 1);
      const bool nextCombination
         = L::next_combination (&combination[0], &combination[slots]);
      const bool nextCombination32
         = L::next_combination (&combination32, slots);

      // Make sure, all next_X()-functions agree

      const int result = int (nextPartition) + int (nextCombination)
                                             + int (nextCombination32);

      if (result != 0 && result != 3)
      {
         error ("next_X() functions do not agree");
      }

      if (! result)  break;
   }

   /**
    * The number of partitions is given by
    *
    *    (slots)
    *    (items)
    *
    * Take _item_ items and _slot_ positions.
    */

   DEB1 cout << '\n';
   NORMAL cout << "Number of combinations = " << combinationCount << endl;
   if (combinationCount != L::choose (slots, items))
   {
      error ("Number of generated combinations is wrong");
   }
}

