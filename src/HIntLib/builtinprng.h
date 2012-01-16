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

/**
 *  BuiltIn PRNG
 *
 *  Makes the built-in PRNG available
 */

#ifndef HINTLIB_BUILTINPRNG_H
#define HINTLIB_BUILTINPRNG_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_CSTDLIB
  #include <cstdlib>
  #define HINTLIB_SLN std::
#else
  #include <stdlib.h>
  #define HINTLIB_SLN
#endif

#ifdef HINTLIB_HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/fallback_limits.h>
#endif


namespace HIntLib
{

class BuiltInPRNG
{
private:

   #if HINTLIB_STATIC_WORKS == 1
      HINTLIB_IMPORT static const real RANGE;
      HINTLIB_IMPORT static const real RESOLUTION;
   #else
      static real RANGE;
      static real RESOLUTION;
   #endif

   HINTLIB_IMPORT static bool inUse;

   // Don't copy, dont assign!!!

   BuiltInPRNG (const BuiltInPRNG &);
   BuiltInPRNG& operator= (const BuiltInPRNG&);

public:

   // Types

   typedef int ReturnType;

   // Information about the generator

   static ReturnType         getMax() { return RAND_MAX; }
   static const real&      getRange() { return RANGE; }
   static const real& getResolution() { return RESOLUTION; }

   // Return a random number

   u32 operator() ()  { return HINTLIB_SLN rand(); }  // {0,...,getMax()}
   int operator() (int max);              // {0,...,max-1}
   real getReal();

   // Initialize Generator

   void init (unsigned seed)  { HINTLIB_SLN srand (seed); }

   // Constructors

   BuiltInPRNG (unsigned start = 0);
   ~BuiltInPRNG ()  { inUse = false; }
};


/**
 *  getReal()
 *
 *  Returns a random real from [0,1]
 */

inline
real BuiltInPRNG::getReal()
{
   return (real(operator()()) + real(0.5)) * RESOLUTION;
}


inline
int BuiltInPRNG::operator() (int max)
{
   if (std::numeric_limits<int>::digits + std::numeric_limits<unsigned>::digits
    <= std::numeric_limits<u64>::digits
    && RAND_MAX + 1 != 0)
   {
      return int ((u64(operator()()) * max) / (RAND_MAX + 1));
   }
   else
   {
      return int (getReal() * max);
   }
}

} // namespace HIntLib

#undef HINTLIB_SLN

#endif

