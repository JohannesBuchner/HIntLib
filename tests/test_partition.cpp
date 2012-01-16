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


#include <iostream>
#include <iomanip>

#include <HIntLib/hlalgorithm.h>
#include <HIntLib/hlmath.h>
#include <HIntLib/array.h>

#include "test.h"

using std::cout;
using std::cerr;
using std::endl;

const char* options = "r:";

int RESTRICTION = -1;

bool opt(int c, const char* s)
{
   switch (c)
   {
   case 'r':  RESTRICTION = atoi (s); return true;
   }

   return false;
}

void usage()
{
   cerr <<
      "Usage: test_partition slots items [OPTION]...\n\n"
      "Lists all partitions of _items_ items partitioned into _slots_ slots.\n\n"
      << option_msg <<
      "  -r n   Restrict to partitions with no more than n items per slot\n"
      "\n";

   exit (1);
}

void test(int argc, char** argv)
{
   if (argc != 2)  usage();

   const int slots = atoi (argv [0]);
   const int items = atoi (argv [1]);
   
   if (slots < 1 || items < 0)  usage();

   if (RESTRICTION < 0)  RESTRICTION = items;

   if (items > RESTRICTION * slots)
   {
      cerr << "Restriction too low!" << endl;
      usage();
   }

   L::Array<int> partition (slots);
   L::initial_partition (&partition[0], &partition[slots], items);

   int partitionCount = 0;
   int restrictedPartitionCount = 0;

   // Create normal partitions

   do
   {
      int count = 0;
      bool restricted = true;
      for (int i = 0; i < slots; ++i)
      {
         int x = partition[i];
         DEB1 cout << x << ' ';
         count += x;
         if (x < 0)  error ("Value < 0!");
         if (x > RESTRICTION)  restricted = false;
      }
      DEB1 cout << endl;
      partitionCount ++;
      if (restricted)  restrictedPartitionCount++;

      if (count != items)  error ("Total number of items wrong!");
   }
   while (L::next_partition (&partition[0], &partition[slots]));

   /**
    * The number of partitions is given by
    *
    *    (slots + items - 1)
    *    (      items      )
    *
    * Take _item_ items and _slot_-1 separators. In how many different ways can
    * these  slots+items-1  things be arranged?
    */

   NORMAL cout << "# partitions = " << partitionCount << endl;
   if (partitionCount != L::choose (slots + items - 1, items))
   {
      error ("Number of partitions is wrong!");
   }
   DEB1 cout << "# restricted partitions = "
             << restrictedPartitionCount << endl;

   // Create restricted partitions

   L::initial_partition (&partition[0], &partition[slots], RESTRICTION, items);

   partitionCount = 0;

   do
   {
      int count = 0;
      for (int i = 0; i < slots; ++i)
      {
         int x = partition[i];
         DEB1 cout << x << ' ';
         count += x;
         if (x < 0)  error ("Value < 0!");
         if (x > RESTRICTION)  error ("Restriction not obeyed!");
      }
      DEB1 cout << endl;
      partitionCount ++;

      if (count != items)  error ("Total number of items wrong!");
   }
   while (L::next_partition (&partition[0], &partition[slots], RESTRICTION));

   // Compare number with previous count
 
   NORMAL cout << "# restricted partitions = " << partitionCount << endl;
   if (partitionCount != restrictedPartitionCount)
   {
      error ("Number of partitions is wrong!");
   }
}

