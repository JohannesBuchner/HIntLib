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

#include <HIntLib/mymath.h>

#include "test.h"

using std::cout;
using std::endl;

const char* options = "";

bool opt(int, const char*) { return false; }

void usage()
{
   std::cerr <<
      "Usage: test_abs [OPTION]...\n\n"
      << option_msg <<
      "\n";

   exit (1);
}

using namespace HIntLib;

// Constants

const       float mFloat  = -17.001f;
const       float pFloat  =  17.001f;
const      double mDouble = -17.001;
const      double pDouble =  17.001;
const long double mLong
   = - static_cast<long double> (17001) / static_cast<long double> (1000);
const long double pLong
   =   static_cast<long double> (17001) / static_cast<long double> (1000);
#if HINTLIB_REAL > 3
const real mReal = -static_cast<real>(17001) / static_cast<real>(1000);
const real pReal =  static_cast<real>(17001) / static_cast<real>(1000);
#endif

// determining the type of an expression

template<typename T> int type (const T&)  { return -1; }

template<> int type (const float&)  { return 1; }
template<> int type (const double&)  { return 2; }
template<> int type (const long double&)  { return 3; }


/**
 *  Main program
 */

void test(int argc, char**)
{
   if (argc)  usage();

   cout << std::setprecision (25);

   // float

   NORMAL  cout << "Testing float..." << endl;

   int t = type (abs (mFloat));

   DEB1 cout << "x = " << mFloat << "   abs(x) = " << abs (mFloat)
             << "   expected = " << pFloat
             << "   type is " << t << endl;

   if (abs (mFloat) != pFloat)  error ("accuracy failed!");
   if (t != 1)  error ("type failed!");
   
   // double

   NORMAL  cout << "Testing double..." << endl;

   t = type (abs (mDouble));
   
   DEB1 cout << "x = " << mDouble << "   abs(x) = " << abs (mDouble)
             << "   expected = " << pDouble
             << "   type is " << t << endl;

   if (abs (mDouble) != pDouble)  error ("accuracy failed!");
   if (t != 2)  error ("type failed!");
   
   // long double

   NORMAL  cout << "Testing long double..." << endl;

   t = type (abs (mLong));
   
   DEB1 cout << "x = " << mLong << "   abs(x) = " << abs (mLong)
             << "   expected = " << pLong
             << "   type is " << t << endl;

   if (abs (mLong) != pLong)
   {
      if (double (abs (mLong)) != double (pLong))  error ("accuracy failed!");
      else
      {
         cout << "*** Warning: Accuracy of long double only like double!"
              << endl;
      }
   }
   if (t != 3)  error ("type failed!");
   
   // real

#if HINTLIB_REAL > 3
   NORMAL  cout << "Testing real..." << endl;
   DEB1 cout << "x = " << mReal << "   fabsl(x) = " << abs (mReal)
             << "   expected = " << pReal << endl;

   if (abs (mReal) != pReal)  error ("accuracy failed!");
#endif
}

