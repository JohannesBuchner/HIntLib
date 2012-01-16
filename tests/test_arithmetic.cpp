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

#include <HIntLib/output.h>

#include "test_arithmetic.h"

/*
 *  Global variables
 */

// Userdefined constants

unsigned SIZE = 100;
bool FLUSH = false;
unsigned W = 15;

STRING selectedCategory;
STRING selectedName;
std::string selectedType;

// counters

unsigned nilpotentsCounter;
unsigned unitsCounter;
unsigned primitivesCounter;
unsigned lastNorm;
int lastDegree;

/**
 *  Comamnd line arguments
 */

const char option_msg[] =
      "  -n n   Uses tables of size _n_ for the tests (default = 1000).\n"
      "  -f     Flush the output stream at regular times.\n"
      "  -c str Perform only tests for given algebra category.\n"
      "  -m str Perform only tests for given structure name.\n"
      "  -t str Perform only tests for given structure type.\n";
const char options[] = "n:fc:m:t:";
const char testProgramParameters[] = "[OPTION]...";
const char testProgramUsage[] = "Tests all kinds of rings and fields.\n\n";
const char testProgramName[] = "test_arithmetic";
const int  testProgramCopyright = 2003;


bool opt (int c, const char* s)
{
   switch (c)
   {
   case 'n':  SIZE = parseInt(s); return true;
   case 'f':  FLUSH = true; return true;
   case 'c':
      {
#ifdef USE_WCHAR
         STRINGSTREAM ss;
         ss << s;
         selectedCategory = ss.str();
#else
         selectedCategory = s;
#endif
         return true;
      }
   case 'm':
      {
#ifdef USE_WCHAR
         STRINGSTREAM ss;
         ss << s;
         selectedName = ss.str();
#else
         selectedName = s;
#endif
         return true;
      }
   case 't':  selectedType = s; return true;
   }

   return false;
}


/**
 *  performTest()
 */

bool performTest (
      const STRING& cat, const STRING& name, const std::string& type)
{
   if (((   selectedCategory.empty()
         || cat .find (selectedCategory) != STRING::npos) &&
        (   selectedName.empty()
         || name.find (selectedName)     != STRING::npos) &&
        (   selectedType.empty()
         || type.find (selectedType)     != std::string::npos)))
   {
      NORMAL
      {
         COUT << "Testing " << cat << ' ';
         doubleQuote (COUT, name.c_str());
         COUT << " (" << type.c_str() << ")..." << std::endl;
      }
      else if (verbose == 0)
      {
         COUT << '.' << std::flush;
      }

      return true;
   }

   return false;
}


/**
 * printInfinity ()
 */

void printInfinity ()
{
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 3
   COUT << L"\x221e";
# else
   COUT << L"inf";
# endif
#else
   COUT << Wgl4Ascii ("\xe2\x88\x9e", "inf");
#endif
}


/**
 * printNumberOrInf ()
 */

void printNumberOrInf (unsigned n)
{
   if (n)  COUT << n;
   else printInfinity();
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
   if (argc)  usage("Too many arguments!");

#ifdef USE_WCHAR
#ifdef HINTLIB_STREAMS_SUPPORT_LOCALE
   COUT.imbue (std::cout.getloc());
#endif
#endif

   NORMAL printHeader (COUT);

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

   if (verbose == 0)  COUT << std::endl;
}

