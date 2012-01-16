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

using std::string;

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
   return ((   selectedCategory == ""
            || cat .find (selectedCategory) != string::npos) &&
           (   selectedName == ""
            || name.find (selectedName)     != string::npos) &&
           (   selectedType == ""
            || type.find (selectedType)     != string::npos));
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
}

