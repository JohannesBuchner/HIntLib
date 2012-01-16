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

/**
 *  KahanAdd
 *
 *  Compensated (Kahan) Summation
 *
 *  This class acts like a simple real when it is used for summation (wiht +=).
 *  However, its much more stable numerically.
 *
 *  See Krommer, Ueberhuber. Computational Integration. p 333 vor details
 */

#ifndef HINTLIB_KAHAN_ADD_H
#define HINTLIB_KAHAN_ADD_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/defaults.h>

namespace HIntLib
{

class KahanAdd
{
public:

   // Constructors

   KahanAdd (real init = 0.0) : sum (init), error (0.0) {}
   KahanAdd (const KahanAdd &k) : sum (k.sum), error (k.error) {}

   // Assignment

   KahanAdd& operator= (real x) { sum = x; error = 0.0; return *this; }
   KahanAdd& operator= (const KahanAdd &k)
      { sum = k.sum; error = k.error; return *this; }

   // Summation

   KahanAdd& operator+= (real x);

   // Conversion to real

   operator real () const  { return sum; }

private:
   volatile real sum;
   real error;
};


/******* Implementation *******/

inline
KahanAdd& KahanAdd::operator+= (real x)
{
   real adjustedX = x - error;

   // newSum MUST be truncated to a double. If it remains in the full internal
   // precision of the i387, the whole algorithm breaks.
   // 
   // volatile is used to force the compile to store newSum in memory and read
   // it back in the next statement.
   // This works with gcc 2.95.2 even if -O3 is turned on

   volatile real newSum = sum + adjustedX;
   error = (newSum - sum) - adjustedX;
   sum = newSum;

   return *this;
}

}   // namespace HIntLib

#endif
