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
#include <iomanip>

#include <HIntLib/tparameter.h>
#include <HIntLib/generatormatrixgen.h>
#include <HIntLib/shiftnet.h>

#include "test.h"

using namespace HIntLib;

using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::setw;

// Option processing

int MAX_M = 0;
int MIN_M = 1;

const char options[] = "m:i:";
const char option_msg[] =
   "  base   Base of the shift net\n"
   "  -m     Maximal m\n"
   "  -i     Minimal m\n";
const char testProgramParameters[] = "[OPTIONS] base";
const char testProgramUsage[] =
    "Test whether the shifnets are generated correctly.\n\n";
const char testProgramName[] = "test_shiftnet";
const int  testProgramCopyright = 2003;

bool opt(int c, const char* s)
{
   switch (c)
   {
   case 'm':  MAX_M = parseInt (s); return true;
   case 'i':  MIN_M = parseInt (s); return true;
   }

   return false;
}


/**
 *  Main program
 */

void test (int argc, char** argv)
{
   // Determine parameters

   if (argc != 1)  usage("Invalid number of arguments!");

   int base = atoi (argv[0]);

   if (! MAX_M)  MAX_M = maxMForShiftNet (base);

   // Check parameters
   
   if (MIN_M < 1)  usage("Minimal m must be positive!");

   if (MAX_M > int (maxMForShiftNet (base)))
   {
      usage("Maximal m is too large!");
   }

   if (MIN_M > MAX_M)
   {
      usage("Minimal m must not be larger than maximal m!");
   }

   // Do Tests

   NORMAL printHeader(cout);
   if (verbose >= 0)
   {
      cout << "Testing shift nets in base b = " << base
           << (verbose >= 2 ? ":" : "...") << endl;
   }

   DEB1
   {
      cout << setw (3) << "m" << setw(11) << "t" << "   Test results\n";
   }

   for (int m = MIN_M; m <= MAX_M; ++m)
   {
      DEB1  cout << setw (3) << m << flush;

      GeneratorMatrixGen<unsigned char> gm (base, m, m, m);
      initShiftNet (gm);

      const unsigned expected = optimalShiftNetT (base, m);

      DEB1  cout << setw(11) << expected << flush;

      if (confirmT(gm, expected))
      {
         DEB1 cout << setw(7) << "ok" << flush;
      }
      else
      {
         DEB1 cout << setw(7) << "failed" << flush;
         error ("Actual t is too high");
      }

      if (expected > 0)
      {
         if (confirmT(gm, expected - 1))
         {
            DEB1 cout << setw(7) << "failed" << flush;
            error ("Actual t is too low");
         }
         else
         {
            DEB1 cout << setw(7) << "ok" << flush;
         }
      }

      DEB1  cout << '\n';
   }
}


