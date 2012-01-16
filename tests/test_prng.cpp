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
 *  test_prng
 *
 *  Exercises the PRNGs
 *
 *  No statistical tests are performed.
 *
 *  The sole purpose of this test program is to ensure that the template code
 *  complies and links.
 */

#include <iostream>
#include <iomanip>

#include <HIntLib/mersennetwister.h>
#include <HIntLib/builtinprng.h>
#include <HIntLib/lcg_pow2.h>
#include <HIntLib/lcg_prime.h>

#include "test.h"

using std::cout;
using std::endl;
using namespace HIntLib;


/**
 *  Command line arguments
 */

const char* options = "";

bool opt(int, const char*) { return false; }

void usage()
{
   std::cerr <<
      "Usage: test_prng [OPTION]...\n\n"
      << option_msg <<
      "\n";

   exit (1);
}

template<class T>
void doTest (const T*)
{
   T g;
   typedef typename T::ReturnType R;

   R max = g.getMax();
   real range = g.getRange();
   real res   = g.getResolution();

   DEB1 cout << "Max=" << max << "  Range=" << range << "  Resolution=" << res
             << endl;

   if (range < max || range > real(max) * 1.000001 + 1.)
   {
      error ("getRange() broken!");
   }

   if (L::abs (1.0 / res - range) > .0001)
   {
      error ("getResolution() broken!");
   }

   real xx = 0;
   unsigned yy = 0;
   R zz = 0;

   for (unsigned i = 0; i < 100000; ++i)
   {
      // (0,1)

      real x = g.getReal();
      if (x <= 0.0 || x >= 1.0)  error ("getReal() failed!");
      xx += x;

      // {0,...,n}

      int y = g (37);
      if (y < 0 || y >= 37)  error ("operator()(int) failed!");
      yy += y;

      // {0,...,MAX}

      R z = g();
      if (z > max)  error ("operator()() failed!");
      zz += z;
   }
}


/**
 *  Main program
 */

void test(int argc, char**)
{
   if (argc)  usage();

   // BuiltInPRNG

   NORMAL cout << "BuiltInPRNG" << endl;
   doTest (static_cast<BuiltInPRNG*>(0));

   // Mersenne Twister

   NORMAL cout << "MersenneTwister" << endl;
   doTest (static_cast<MersenneTwister*>(0));

   // LCG_Pow2<>

   NORMAL cout << "LCG_Pow2_69069_0" << endl;
   doTest (static_cast<LCG_Pow2_69069_0*>(0));

   NORMAL cout << "LCG_Pow2_69069" << endl;
   doTest (static_cast<LCG_Pow2_69069*>(0));

   NORMAL cout << "LCG_Pow2_LavauxJanssens32_0" << endl;
   doTest (static_cast<LCG_Pow2_LavauxJanssens32_0*>(0));

   NORMAL cout << "LCG_Pow2_LavauxJanssens32" << endl;
   doTest (static_cast<LCG_Pow2_LavauxJanssens32*>(0));

   NORMAL cout << "LCG_Pow2_LavauxJanssens48_0" << endl;
   doTest (static_cast<LCG_Pow2_LavauxJanssens48_0*>(0));

   NORMAL cout << "LCG_Pow2_LavauxJanssens48" << endl;
   doTest (static_cast<LCG_Pow2_LavauxJanssens48*>(0));

   NORMAL cout << "LCG_Pow2_Haynes_0" << endl;
   doTest (static_cast<LCG_Pow2_Haynes_0*>(0));

   NORMAL cout << "LCG_Pow2_Haynes" << endl;
   doTest (static_cast<LCG_Pow2_Haynes*>(0));

   // LCG_Prime

   NORMAL cout << "LCG_Prime_IMSL" << endl;
   doTest (static_cast<LCG_Prime_IMSL*>(0));

   NORMAL cout << "LCG_Prime_Fishman" << endl;
   doTest (static_cast<LCG_Prime_Fishman*>(0));

   NORMAL cout << "LCG_Prime_69621" << endl;
   doTest (static_cast<LCG_Prime_69621*>(0));

   NORMAL cout << "LCG_Prime_Lecuyer" << endl;
   doTest (static_cast<LCG_Prime_Lecuyer*>(0));

   // LCG_Combined
   
   NORMAL cout << "LCG_Combined_Lecuyer" << endl;
   doTest (static_cast<LCG_Combined_Lecuyer*>(0));
}

