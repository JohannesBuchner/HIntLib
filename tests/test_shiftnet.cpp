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
#include <stdlib.h>

#include <HIntLib/tparameter.h>
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

const char* options = "m:";

bool opt(int c, const char* s)
{
   switch (c)
   {
   case 'm':  MAX_M =atoi (s); return true;
   }

   return false;
}

void usage()
{
   cerr <<
      "Usage: test_shiftnet base [OPTION] ...\n\n"
      "Test if the shifnets are generated correctly.\n\n"
      "  base   Base of the shift net\n"
      << option_msg <<
      "  -m     Maximal m\n"
      "\n";

   exit (1);
}


/**
 *  Main program
 */

void test (int argc, char** argv)
{
   if (argc != 1)  usage();

   int base = atoi (argv[0]);

   if (! MAX_M)  MAX_M = maxMForShiftNet (base);

   NORMAL cout << "Testing shift nets in base b = " << base << endl;
   DEB1
   {
      cout << setw (3) << "m"
           << setw (6) << "exp t"
           << setw (6) << "act t" << endl;
   }

   for (int m = 1; m <= MAX_M; ++m)
   {
      DEB1  cout << setw (3) << m << flush;

      GeneratorMatrixGen<unsigned char> gm (base, m, m, m);
      initShiftNet (gm);

      unsigned expected = optimalShiftNetT (base, m);

      DEB1  cout << setw(6) << expected << flush;

      unsigned actual = tParameter (gm);

      DEB1  cout << setw(6) << actual << endl;

      if (actual != expected)
      {
         error ("Actual t not equal to expected t!");
      }
   }
}


