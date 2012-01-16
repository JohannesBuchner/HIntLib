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
#include <stdlib.h>

#include <HIntLib/defaults.h>

#ifdef HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/hintlib_limits.h>
#endif

#include <HIntLib/bitop.h>
#include <HIntLib/old_bitop.h>

#include "test.h"

using std::cout;
using std::cerr;
using std::endl;

namespace O = L::NormalBitOp;

using L::ms1;
using L::ls0;
using L::ls1;
using L::bit;
using L::leastSignificantBits;
using std::numeric_limits;


/**
 *  test_ms1()
 */

template<class T>
inline
void test_ms1 (T x)
{
   int i = ms1 (x);

   // compare with old routine

   if (i != O::ms1 (x))
   {
      cout << "ms1() != O::ms1() for " << x << endl;
      error();
   }

   // range

   if (i < -1 || i > int (sizeof(T) * 8) - 1)
   {
      cout << "ms1(): " << i << " out of range for " << x << endl;
      error();
   }

   // Is it a one?

   if (x == 0)
   {
      if (i != -1)
      {
         cout << "ms1(0) != -1" << endl;
         error();
      }
   }
   else
   {
      if (! bit (x, i))
      {
         cout << "ms1()-bit is not 1 for " << x << endl;
         error();
      }
   }

   if (i + 1 < int (sizeof(T) * 8))
   {
      T mask = T(-1) << (i+1);   // all bits above the msb

      if (x & mask)
      {
         cout << "ms1(): bits above msb are set for " << x << endl;
         error();
      }
   }
}


/**
 *  test_ls0()
 */

template<class T>
inline
void test_ls0 (T x)
{
   unsigned i = ls0 (x);

   if (i != O::ls0 (x))
   {
      cout << "ls0() != O::ls0() for " << x << endl;
      error();
   }

   // range
 
   if (i > sizeof(T) * 8)
   {
      cout << "ls0(): " << i << " out of range for " << x << endl;
      error();
   }

   // Is it 0?

   if (i < sizeof(T)*8 && bit (x, i))
   {
      cout << "ls0()-bit is not 0 for " << x << endl;
      error();
   }

   // All lower bits 1?

   if (i > 0)
   {
      if (leastSignificantBits (x + 1, i))
      {
         cout << "ls0: There are 0s below ls0() for " << x << endl;
         error();
      }
   }
}


/**
 *  test_ls1()
 */

template<class T>
inline
void test_ls1 (T x)
{
   int i = ls1 (x);

   if (i != O::ls1 (x))
   {
      cout << "ls1() != O::ls1() for " << x << endl;
      error();
   }
 
   // range
 
   if (i < -1 || i > int (sizeof(T) * 8) - 1)
   {
      cout << "ls1(): " << i << " out of range for " << x << endl;
      error();
   }
 
   // Is it 1?
 
   if (x == 0)
   {
      if (i != -1)
      {
         cout << "ls1(0) != -1" << endl;
         error();
      }
   }
   else
   {
      if (! bit (x, i))
      {
         cout << "ls1()-bit is not 1 for " << x << endl;
         error();
      }
   }
 
   // All lower bits 0?
 
   if (i > 0)
   {
      if (leastSignificantBits (x, i))
      {
         cout << "ls1: There are 1s below ls1() for " << x << endl;
         error();
      }
   }
}


/**
 *  Loop through all number of type T
 *
 *  Call all tests for each number
 */

template<class T>
inline
void test (const T &dummy)
{
   if (numeric_limits<T>::digits < 20)
   {
      T i = T();

      do
      {
         test_ms1 (i);
         test_ls0 (i);
         test_ls1 (i);

      } while (++i);
   }
   else
   {
      for (T i = 0; i < 0x133333ull; ++i)
      {
         test_ms1 (i);
         test_ls0 (i);
         test_ls1 (i);
      }
   
      T i = numeric_limits<T>::max() - 0x133333;
   
      do
      {
         test_ms1 (i);
         test_ls0 (i);
         test_ls1 (i);
    
      } while (++i);
   }
}

const char* options = "";

bool opt(int, const char*) { return false; }

void usage()
{
   cerr <<
      "Usage: test_bitop [OPTION]...\n\n"
      << option_msg <<
      "\n";

   exit (1);
}

void test(int argc, char**)
{
   if (argc)  usage();
   
   NORMAL  cout << "Testing type unsigned char" << endl;

   test (static_cast<unsigned char> (0));

   NORMAL  cout << "Testing type unsigned short" << endl;

   test (static_cast<unsigned short> (0));

   NORMAL  cout << "Testing type unsigned int" << endl;

   test (static_cast<unsigned> (0));

   NORMAL  cout << "Testing type unsigned long" << endl;

   test (static_cast<unsigned long> (0));

#ifdef HAVE_UNSIGNED_LONG_LONG_INT
   NORMAL  cout << "Testing type unsigned long long" << endl;

   test (static_cast<unsigned long long> (0));
#endif
}

