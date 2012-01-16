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

#include <HIntLib/polynomial2.h>
#include <HIntLib/eratosthenes.h>
#include <HIntLib/array.h>

#include "test.h"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;

using HIntLib::u32;


u32 SIZE = 1000;

template<class P>
void doIt (const u32 S, const P& dummy)
{
   L::Array<bool> f (S);

   L::eratosthenes (f, S, dummy);

   for (u32 i = 0; i < S; i++)
   {
      P x (i);

      DEB1  cout << setw(5) << i << setw(25) << P (i);

      if (f[i])
      {
         DEB1  cout << " is irreducible\n";

         for (u32 j = 2; j < S; ++j)
         {
            if (P(j) == x)  continue;

            if (x % P(j) == P(0))  error ("not irreducible");
         }
      }
      else
      {
         DEB1  cout << '\n';

         if (x == P(1))  continue;

         bool found = false;

         for (u32 j = 2; j < S; ++j)
         {
            if (P(j) == x)  continue;

            if (x % P(j) == P(0))  found = true;
         }

         if (! found)  error ("irreducible");
      }
   }
}

const char* options = "n:";

bool opt (int c, const char* s)
{
   switch (c)
   {
   case 'n':  SIZE = atoi(s) + 1; return true;
   }

   return false;
}

void usage()
{
   cerr <<
      "Usage: test_eratosthenes [OPTION]...\n\n"
      "Tests Eratosthenes' sieve.\n"
      "The output is compared to a manual divisability check.\n"
      "The test is performed for integers (u32) and polynomials.\n\n"
      << option_msg <<
      "  -n n   Size of the sieve (default = 1000).\n"
      "\n";

   exit (1);
}

void test (int argc, char **)
{
   if (argc)  usage();

   NORMAL  cout << "Testing Eratosthenes for 'u32'..." << endl;
   doIt (SIZE, u32 (0));
   NORMAL  cout << "Testing Eratosthenes for 'Polynomial2'..." << endl;
   doIt (SIZE, L::Polynomial2 (0));
}

