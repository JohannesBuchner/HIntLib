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

#include <HIntLib/polynomial2.h>
#include <HIntLib/mymath.h>

#include "test.h"

using std::cerr;
using std::cout;
using std::endl;

typedef L::Polynomial2<L::u32> P;

unsigned SIZE = 1000;

void powers()
{
   DEB1   cout << "\n\n";
   NORMAL cout << "Checking powers...\n";
   DEB1   cout << "\n";

   const int S = 100;

   for (int i = 0; i < S; i++)
   {
      P x (i);
      DEB1 cout << "  Checking: " << x << endl;

      for (int j = 0; j < 5; ++j)
      {
         P xx = L::powInt (x, j);
         DEB2 cout << "    (" << x << ") ^ " << j << " = " << xx << endl;

         P y (1);
         for (int k = 0; k < j; ++k)  y *= x;

         if (xx != y)  error("x^k wrong");
      }
   }
}


void primitive()
{
   DEB1   cout << "\n\n";
   NORMAL cout << "Searching for primitive polynomials...\n";
   DEB1   cout << "\n";

   P lastPrim (0);
   int found = 0;
   int expected = 0;

   const unsigned S = SIZE;

   for (unsigned i = 0; i < S; ++i)
   {
      P p (i);

      if (p.isPrimitive ())
      {
         const int deg = p.degree();

         if (deg > lastPrim.degree())
         {
            if (found != expected)  error ("Wrong number found");

            lastPrim = p;
            found = 0;
            expected = HIntLib::Prime::eulerPhi ((1u << deg) - 1) / deg;

            DEB1  cout << "  Degree " << deg
                       << "    expected number = " << expected << endl;
         }

         ++ found;
         DEB2  cout << "    " << p << endl;
      }
   }
}


const char* options = "n:";

bool opt (int c, const char* s)
{
   switch (c)
   {
   case 'n':  SIZE = atoi (s); return true;
   }

   return false;
}

void usage()
{
   cerr <<
      "Usage: test_polynomial [OPTION]...\n\n"
      "Tests if the type Polynomial2<u32> behaves as expected.\n\n"
      << option_msg <<
      "  -n n   Uses tables of size _n_ for the tests (default = 1000).\n"
      "\n";

   exit (1);
}

void test (int argc, char**)
{
   if (argc)  usage();

   powers();
   primitive();
}

