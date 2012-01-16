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
 *  create_prime
 *
 *  Creates the source file prime.cpp
 */

#include <iostream>
#include <iomanip>
#include <string.h>
#include <stdlib.h>

#include "defaults.h"
#include "eratosthenes.h"

using std::cout;
using std::setw;

int main (int argc, const char** argv)
{
   // Get mode from the command line

   if (argc != 2)  exit(1);

   enum {H, CPP} mode;

   if (strcmp(argv[1], "h") == 0)
      mode = H;
   else if (strcmp(argv[1], "cpp") == 0)
      mode = CPP;
   else
      exit (1);

   // Define and initialize the field

   const unsigned int SIZE = HINTLIB_PRIME_TABLE_SIZE;

   bool field [2 * SIZE];

   HIntLib::eratosthenes (field, 2 * SIZE, unsigned (0));

   // Print header

   cout <<
      "/*******************************************************/\n"
      "/***   This file is program-generated!               ***/\n"
      "/***                                                 ***/\n"
      "/***   Do not change!!!                              ***/\n"
      "/***                                                 ***/\n"
      "/***   Update " __FILE__ " to update this file.  ***/\n"
      "/*******************************************************/\n\n";

   if (mode == CPP)
   {
      cout <<
         "#ifdef __GNUG__\n"
         "#pragma implementation\n"
         "#endif\n"
         "\n"
         "#include <HIntLib/prime.h>\n"
         "\n"
         "const unsigned short HIntLib::Prime::nextPrimeArray [MAX_N+1] =\n"
         "{\n";

      // Create isPrime array 

      for (unsigned i = 0; i <= SIZE; i++)
      {
         unsigned j = i;

         while (! field [j]) ++j;

         cout << setw (6) << j << ",   // " << i << '\n';
      }

      cout <<
         "};\n\n"

      // Create nthPrime array

         "const unsigned short HIntLib::Prime::nthPrimeArray [NUM_PRIMES] =\n"
         "{\n";
   }

   unsigned count = 0;

   for (unsigned i = 0; i < SIZE; i++)
   {
      if (field [i])
      {
         if (mode == CPP)  cout << setw(6) << i << ",  // " << count << '\n';
         count++;
      }
   }

   if (mode == CPP)  cout << "};\n\n";

   if (mode == H)
   { 
      cout 
      << "#ifndef PRIME_H\n"
         "#error \"prime_generated.h must not be included directly!\"\n"
         "#endif\n"
         "\n"
         "static const unsigned MAX_N = " << SIZE << ";\n"
         "static const unsigned NUM_PRIMES = " << count << ";\n\n";
   }

   return 0;
}

