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

#include <HIntLib/prime.h>
#include <HIntLib/gcd.h>

#include "test.h"

using std::cout;
using std::cerr;
using std::setw;
using std::endl;
using std::flush;

using HIntLib::u32;

u32 SIZE = 1000;

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
      "Usage: test_prime [OPTION]...\n\n"
      "Tests if the functions in Prime work:"
              " test(), next(), nth(), eulerPhi()\n\n"
      << option_msg <<
      "  -n n   Largest number to test for (default = 1000).\n"
      "\n";

   exit (1);
}

void test (int argc, char**)
{
   if (argc)  usage();
   
   NORMAL  cout << "Searching for prime numbers..." << endl;

   u32 n = 0;
   u32 last = 0;

   for (u32 i = 0; i < SIZE; ++i)
   {
      DEB1  cout << setw(5) << i;
      DEB2  cout << flush;

      bool p = L::Prime::test(i);

      NORMAL
      {
         DEB1  cout << setw(10) << (p ? " prime" : "not prime");
         else if (p)  cout << setw(5) << i << endl;
         DEB2  cout << flush;
      }

      if (p)  // prime
      {
         DEB2  cout << "." << flush;
         for (u32 j = 2; j < i; ++j)
         {
            if (i % j == 0)  error ("Should not be prime");
         }

         DEB2  cout << "." << flush;
         try
         {
            if (i != L::Prime::nth(n++))  error ("nth() failed");
         }
         catch (L::PrimeNumberNth &) {}

         DEB2  cout << "." << flush;
         for (u32 j = last+1; j <= i; ++j)
         {
            if (L::Prime::next(j) != i)  error ("next() failed");
         }

         DEB2  cout << "." << flush;
         if (i != L::Prime::eulerPhi(i) + 1)  error ("eulerPhi() failed");

         last = i;
         DEB2  cout << "." << flush;
      }
      else   // not prime
      {
         if (i > 1)
         {
            u32 numDivisor = 0;
            u32 numRelPrime = 0;

            for (u32 j = 2; j < i; ++j)
            {
               if (i % j == 0)  ++numDivisor;
               if (L::gcd (i, j) == 1)  ++numRelPrime;
            }

            if (numDivisor == 0)  error ("No divisors found");

            if (L::Prime::eulerPhi(i) != numRelPrime+1)
               error ("eulerPhi() failed");
         }
      }

      DEB1  cout << endl;
   }
}

