/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schürer <rudolf.schuerer@sbg.ac.at>
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

#ifndef HINTLIB_COUNTER_H
#define HINTLIB_COUNTER_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/array.h>

namespace HIntLib
{

class Counter
{
public:
   Counter (unsigned _digits, unsigned _base)
      : digits (_digits), baseMinusOne (_base - 1),
        index (digits, 0) {}

   unsigned operator[] (unsigned i) const  { return index[i]; }
   bool next();

private:
   const unsigned digits;
   const unsigned baseMinusOne;
   Array<unsigned> index;
};

class CounterMixedBase
{
public:
   template<typename T>
   CounterMixedBase (T begin, T end)
      : digits (end - begin), baseMinusOne (digits),
        index (digits, 0)
   {
      unsigned* p = baseMinusOne.begin();
      while (begin != end)  *p++ = *begin++ - 1;
   }

   unsigned operator[] (unsigned i) const  { return index[i]; }
   bool next();

private:
   const unsigned digits;
   Array<unsigned> baseMinusOne;
   Array<unsigned> index;
};

}  // namespace HIntLib

#endif

