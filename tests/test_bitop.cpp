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

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_LIMITS
#  include <limits>
#else
#  include <HIntLib/fallback_limits.h>
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
using L::lsBits;
using L::lsBitsMask;
using L::popCount;
using L::parity;
using std::numeric_limits;


/**
 *  test_ms1()
 */

template<class T>
inline
void test_ms1 (T x)
{
   const int i = ms1 (x);

   // compare with old routine

   if (i != O::ms1 (x))
   {
      cout << "ms1() != O::ms1() for " << x << endl;
      error();
   }

   // range

   if (i < -1 || i > std::numeric_limits<T>::digits - 1)
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

   if (i + 1 < std::numeric_limits<T>::digits)
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
   const int i = ls0 (x);

   if (i != int(O::ls0 (x)))
   {
      cout << "ls0() != O::ls0() for " << x << endl;
      error();
   }

   // range
 
   if (i > int(std::numeric_limits<T>::digits))
   {
      cout << "ls0(): " << i << " out of range for " << x << endl;
      error();
   }

   // Is it 0?

   if (i < std::numeric_limits<T>::digits && bit (x, i))
   {
      cout << "ls0()-bit is not 0 for " << x << endl;
      error();
   }

   // All lower bits 1?

   if (lsBits (x, i) != lsBitsMask<T> (i))
   {
      cout << "ls0: There are 0s below ls0() for " << x << endl;
      error();
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
 
   if (i < -1 || i > std::numeric_limits<T>::digits - 1)
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
 
   if (lsBits (x, i))
   {
      cout << "ls1: There are 1s below ls1() for " << x << endl;
      error();
   }
}


/**
 *  test_popcount()
 */

template<class T>
inline
void test_popcount (T x)
{
   int i = popCount (x);

   // range
 
   if (i < 0 || i > std::numeric_limits<T>::digits)
   {
      cout << "popCount(" << x << ") = " << i << " is out of range!\n";
      error();
   }
 
   // Count digits

   int counter = 0;
   {
      T xx = x;
      while (xx)
      {
         if (xx & 1)  ++counter;
         xx >>= 1;
      }
   }

   if (counter != i)
   {
      cout << "popCount(" << x << ") = " << i << " should be "
           << counter << "!\n";
      error();
   }
}


/**
 *  test_parity()
 */

template<class T>
inline
void test_parity (T x)
{
   int i = parity (x);

   // range
 
   if (i < 0 || i > 1)
   {
      cout << "parity(" << x << ") = " << i << " should be in {0,1}!\n";
      error();
   }
 
   // Count digits

   int par = 0;
   {
      T xx = x;
      while (xx)
      {
         if (xx & 1)  par ^= 1;
         xx >>= 1;
      }
   }

   if (par != i)
   {
      cout << "parity(" << x << ") = " << i << " should be " << par << "!\n";
      error();
   }

   // Compare with popcount()

   if (i != (popCount(x) & 1))
   {
      cout << "parity(" << x << ") = " << i
           << " does not match popCount() & 1!\n";
      error();
   }
}


/**
 *  Loop through all number of type T
 *
 *  Call all tests for each number
 */

template<class T>
inline
void test (const T &)
{
   if (numeric_limits<T>::digits < 22)
   {
      T i = T();

      do
      {
         DEB1  cout << i;
         test_ms1 (i);
         test_ls0 (i);
         test_ls1 (i);
         test_popcount (i);
         test_parity (i);
         DEB1 cout << '\n';

      } while (++i);
   }
   else
   {
      for (L::u32 i = 0; i < 0x133333ull; ++i)
      {
         T ii (i);
         DEB1  cout << i;
         test_ms1 (ii);
         test_ls0 (ii);
         test_ls1 (ii);
         test_popcount (ii);
         test_parity (ii);
         DEB1 cout << '\n';
   
         ii = numeric_limits<T>::max() - i;
         DEB1  cout << ii;
         test_ms1 (ii);
         test_ls0 (ii);
         test_ls1 (ii);
         test_popcount (ii);
         test_parity (ii);
         DEB1 cout << '\n';
      }
   }
}

const char options[] = "";
const char option_msg[] = "";
const char testProgramParameters[] = "[OPTION]...";
const char testProgramUsage[] = "";
const char testProgramName[] = "test_bitop";
const int  testProgramCopyright = 2003;

bool opt(int, const char*) { return false; }

void test(int argc, char**)
{
   if (argc)  usage("Too many arguments!");

   NORMAL printHeader (cout);
   
   NORMAL  cout << "Testing type unsigned char" << endl;

   test (static_cast<unsigned char> (0));

   NORMAL  cout << "Testing type unsigned short" << endl;

   test (static_cast<unsigned short> (0));

   NORMAL  cout << "Testing type unsigned int" << endl;

   test (static_cast<unsigned> (0));

   NORMAL  cout << "Testing type unsigned long" << endl;

   test (static_cast<unsigned long> (0));

#ifdef HINTLIB_HAVE_UNSIGNED_LONG_LONG_INT
   NORMAL  cout << "Testing type unsigned long long" << endl;

   test (static_cast<unsigned long long> (0));
#endif
}

