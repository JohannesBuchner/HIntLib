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

#ifndef HINTLIB_COUNTER_H
#define HINTLIB_COUNTER_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <iosfwd>

#include <HIntLib/array.h>

namespace HIntLib
{

/**
 *  CounterBase
 *
 *  The common base class for all counters
 */

class CounterBase
{
protected:
   CounterBase (unsigned _digits) : digits (_digits), index (digits, 0) {}

public:
   unsigned operator[] (unsigned i) const  { return index[i]; }
   const unsigned* array() const  { return index.begin(); }
   unsigned getNumDigits() const  { return digits; }

protected:
   void reset();
   
   const unsigned digits;
   Array<unsigned> index;
};

std::ostream& operator<< (std::ostream &, const CounterBase &);
#ifdef HINTLIB_BUILD_WCHAR
std::wostream& operator<< (std::wostream &, const CounterBase &);
#endif


/**
 *  Counter
 */

class Counter : public CounterBase
{
public:
   Counter (unsigned _digits, unsigned _base)
      : CounterBase (_digits), baseMinusOne (_base - 1) {}

   bool next();
   using CounterBase::reset;

private:
   const unsigned baseMinusOne;
};


/**
 *  CounterMixedBase
 */

class CounterMixedBase : public CounterBase
{
public:
   template<typename T>
   CounterMixedBase (T begin, T end)
      : CounterBase (end - begin), baseMinusOne (digits)
   {
      unsigned* p = baseMinusOne.begin();
      while (begin != end)  *p++ = *begin++ - 1;
   }

   bool next();
   using CounterBase::reset;

private:
   Array<unsigned> baseMinusOne;
};


/**
 *  GrayCounter
 */

class GrayCounter : public CounterBase
{
public:
   GrayCounter (unsigned _digits, unsigned _base)
      : CounterBase (_digits),
        baseMinusOne(_base - 1), direction(_digits, 1) {}

   bool next();
   unsigned nextDigit();
   unsigned nextDigit(unsigned*);
   void reset();

private:
   const unsigned baseMinusOne;
   Array<int> direction;
};

}  // namespace HIntLib

#endif

