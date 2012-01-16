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

#ifdef __GNUG__
#pragma implementation
#endif

#define HINTLIB_LIBRARY_OBJECT

/**
 *  polynomialbase.cpp
 *
 *  Defines all non-template members of PolynomialRing<> base classes.
 */

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_OSTREAM
  #include <ostream>
#else
  #include <iostream>
#endif

#include <HIntLib/polynomialbase.h>

namespace P = HIntLib::Private;


/**
 *  funnySum ()
 *
 *  Return the first  n  terms of the sum
 *
 *      1 + base + base^2 + base^3 + ...
 */

unsigned
P::funnySum (int n, unsigned base)
{
   if (n <= 0)  return 0;

   unsigned res = 0;
   unsigned k = 1;

   for (int i = 0; i < n; ++i)
   {
      res += k;
      k *= base;
   }
   return res;
}


/**
 *  printVariable ()
 */

void
P::PRB::printVariable (std::ostream& o) const
{
   o << var;
}


/**
 *  printVariableWithBrackets ()
 *
 *  This method ist NOT width()-save!
 */

void
P::PRB::printVariableWithBrackets (std::ostream& o) const
{
   o << '[' << var << ']';
}


/**
 *  printVariablePow ()
 */

void
P::PRB::printVariablePow (std::ostream& o, unsigned i) const
{
   o << var;

   switch (i)
   {
   case 1:  break;
   case 2:  o << '\262'; break;
   case 3:  o << '\263'; break;
   default: o << '^' << i;
   }
}


/**
 *  squareBeatsLinear
 *
 *  Performing a full search over all degree-2 polynomials is an advantage if
 *  and only if
 *
 *     deg > squareBeatsLinear[p]
 */

const int P::PRB::squareBeatsLinear [] =
{
       -1,
       -1,
   1,  // p = 2
   1,  // p = 3
   1,  // p = 4
   1,  // p = 5
       -1,
   16, // p = 7
   17, // p = 8
   20, // p = 9
       -1,
   23, // p = 11
       -1,
   28, // p = 13
       -1,
       -1,
   35  // p = 16
};

