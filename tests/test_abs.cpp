/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 03,04,05  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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

#include <HIntLib/hlmath.h>

#include "test.h"

using std::cout;
using std::endl;

const char options[] = "";
const char option_msg[] = "";
const char testProgramParameters[] = "[OPTION]...";
const char testProgramUsage[] = "";
const char testProgramName[] = "test_abs";
const int  testProgramCopyright = 2003;

bool opt(int, const char*) { return false; }

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
   if (argc)  usage("Too many arguments!");

   NORMAL printHeader (cout);

   cout << std::setprecision (25);

   // float

   NORMAL  cout << "Testing float..." << endl;

   int t = type (L::abs (mFloat));

   DEB1 cout << "x = " << mFloat << "   abs(x) = " << L::abs (mFloat)
             << "   expected = " << pFloat
             << "   type is " << t << endl;

   if (L::abs (mFloat) != pFloat)  error ("accuracy failed");
   if (t != 1)  error ("type failed");
   
   // double

   NORMAL  cout << "Testing double..." << endl;

   t = type (L::abs (mDouble));
   
   DEB1 cout << "x = " << mDouble << "   abs(x) = " << L::abs (mDouble)
             << "   expected = " << pDouble
             << "   type is " << t << endl;

   if (L::abs (mDouble) != pDouble)  error ("accuracy failed");
   if (t != 2)  error ("type failed");
   
   // long double

   NORMAL  cout << "Testing long double..." << endl;

   t = type (L::abs (mLong));
   
   DEB1 cout << "x = " << mLong << "   abs(x) = " << L::abs (mLong)
             << "   expected = " << pLong
             << "   type is " << t << endl;

   if (L::abs (mLong) != pLong) error ("accuracy failed");
   if (t != 3)  error ("type failed");
}

