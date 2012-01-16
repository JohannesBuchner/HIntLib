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
#include <string>

#include "test.h"

#ifdef HINTLIB_HAVE_SSTREAM
  #include <sstream>
#else
  #include <HIntLib/fallback_sstream.h>
#endif

using std::string;
using std::cout;

/*
 *  Global variables
 */

// Userdefined constants

unsigned SIZE = 100;
bool FLUSH = false;
unsigned W = 15;

string selectedCategory;
string selectedName;
string selectedType;

// counters

unsigned nilpotentsCounter;
unsigned unitsCounter;
unsigned primitivesCounter;
unsigned lastNorm;
int lastDegree;


/**
 *  Comamnd line arguments
 */

const char* options = "n:fc:m:t:";

bool opt (int c, const char* s)
{
   switch (c)
   {
   case 'n':  SIZE = atoi (s); return true;
   case 'f':  FLUSH = true; return true;
   case 'c':  selectedCategory = s; return true;
   case 'm':  selectedName = s; return true;
   case 't':  selectedType = s; return true;
   }

   return false;
}

void usage()
{
   std::cerr <<
      "Usage: test_arithmetic [OPTION]...\n\n"
      "Tests all kinds of rings and fields.\n\n"
      << option_msg <<
      "  -n n   Uses tables of size _n_ for the tests (default = 1000).\n"
      "  -f     Flush the output stream at regular times.\n"
      "  -c str Perform only tests for given algebra category.\n"
      "  -m str Perform only tests for given structure name.\n"
      "  -t str Perform only tests for given structure type.\n"
      "\n";

   exit (1);
}


/**
 *  performTest()
 */

bool performTest (const string& cat, const string& name, const string& type)
{
   if (((   selectedCategory == ""
         || cat .find (selectedCategory) != string::npos) &&
        (   selectedName == ""
         || name.find (selectedName)     != string::npos) &&
        (   selectedType == ""
         || type.find (selectedType)     != string::npos)))
   {
      NORMAL
      {
         cout << "Testing " << cat
              << " \"" << name << "\" (" << type << ")..." << std::endl;
      }
      else if (verbose == 0)
      {
         cout << '.' << std::flush;
      }

      return true;
   }

   return false;
}


/**
 * printNumberOrInf ()
 */

void printNumberOrInf (unsigned n)
{
   if (n)  cout << n;
   else    cout << "inf";
}


/**
 *  checkCounter()
 */

void checkCounter (unsigned expected, unsigned actual, bool all, const char* s)
{
   if (   (  all && actual != expected)
       || (! all && expected && actual > expected))
   {
      std::ostringstream ss;
      ss << s << " wrong (expected: " << expected << ", got: " << actual << ")";
      error (ss.str().c_str());
   }
}


/**
 *  checkIndex()
 */

void checkIndex (unsigned size, unsigned index, const char* s)
{
   if (size && index >= size)
   {
      std::ostringstream ss;
      ss << "index() of " << s << "()-result invalid";
      error (ss.str().c_str());
   }
}


/**
 *  checkInfiniteInFinite()
 */

void checkInfiniteInFinite (unsigned size, unsigned n, const char* s)
{
   if (size)
   {
      if (n == 0)
      {
         std::ostringstream ss;
         ss << "Infinite number of " << s << " in set of size " << size;
         error (ss.str().c_str());
      }
      else if (n > size)
      {
         std::ostringstream ss;
         ss << n << " " << s << " in set of size " << size;
         error (ss.str().c_str());
      }
   }
}



/**
 *  testXXXX()
 *
 *  Test routines defined in seperate compilation units
 */

namespace HIntLib
{
   void testComplex();
   void testReal();
   void testInteger();
   void testGF2();
   void testPoly2();
   void testGF2Vec();
   void testLookupField (unsigned);
   void testModularArithmetic (unsigned);
   void testModularArithmeticShort (unsigned);

} // namespace HIntLib


/**
 *  test()
 *
 *  Call the test routines defined above
 */

void test (int argc, char**)
{
   if (argc)  usage();

   HIntLib::testComplex();
   HIntLib::testReal();
   HIntLib::testInteger();
   HIntLib::testGF2();
   HIntLib::testPoly2();
   HIntLib::testGF2Vec();

   const unsigned bases [] =
      { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16, 17, 18,
       32,  64, 128, 256, // all prime powers <= 256
       27,  81, 243,
       25, 125,
       49,
       121,
       169,
       43, 251,          // some more primes
       91, 213, 255 };   // some more non-primes

   for (unsigned i = 0; i < sizeof(bases) / sizeof(unsigned); ++i)
   {
      const unsigned size = bases [i];

      HIntLib::testLookupField (size);
      HIntLib::testModularArithmetic (size);
   }

   HIntLib::testModularArithmeticShort (1999);

   if (verbose == 0)  cout << std::endl;
}

