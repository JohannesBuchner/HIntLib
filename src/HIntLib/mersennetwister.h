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
 *  MersenneTwister
 *
 *  The Mersenne Twister Random Number Generator
 *
 *  This implementation is based on
 *     [1] Makoto Matusmoto and Takuji Nishimura. Mersenne Twister: A 623-
 *         dimensionally equidistributed uniform pseudorandom number
 *         generator. ACT Transactions on Modeling and Computer Simulation.
 *         8(1):3-30, January 1998.
 *     [2] Source codes due to Takuji Nishimura, Shawn Cokus,
 *         and Richard J. Wagner.
 */

#ifndef HINTLIB_MERSENNETWISTER_H
#define HINTLIB_MERSENNETWISTER_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_CSTDDEF
  #include <cstddef>
  #define HINTLIB_SDN ::std::
#else
  #include <stddef.h>
  #define HINTLIB_SDN ::
#endif


namespace HIntLib
{

class MersenneTwister
{
private:

   static const unsigned N = 624;

   u32 *next;
   u32 state [N];

   void reload();

   #if HINTLIB_STATIC_WORKS == 1
      static const real RANGE;
      static const real RESOLUTION;
   #else
      static real RANGE;
      static real RESOLUTION;
   #endif

public:

   // Types

   typedef u32 ReturnType;

   // Information about the generator

   static ReturnType         getMax() { return ~u32(); }
   static const real&      getRange() { return RANGE; }
   static const real& getResolution() { return RESOLUTION; }

   static HINTLIB_SDN size_t getStateSize()
      { return N * sizeof(u32) + sizeof(HINTLIB_SDN ptrdiff_t); }

   // Return a random number

   u32 operator() ();         // {0,...,getMax()}
   int operator() (int max);  // {0,...,max-1}
   real getReal();            // (0,1)

   // Save and restore state

   void saveState (void *p) const;
   void restoreState (const void *p);

   // Initialize Generator

   void init (unsigned = 4357u);

   // These two initializers are here for comapatibility only and should not
   // be used

   void initOld (unsigned);
   void initCokus (unsigned start = 4357u)  { initOld (start / 2); }

   // Constructors

   MersenneTwister (unsigned start = 4357u);
   MersenneTwister (const MersenneTwister &);
   MersenneTwister& operator= (const MersenneTwister &);
};


inline
u32 MersenneTwister::operator() ()
{
   if (next == state + N)
   {
      reload ();

      next = state; // Leaving this outside reload() allows better optimization
   }

   u32 y  = *next++;

   y ^= (y >> 11);
   y ^= (y <<  7) & 0x9D2C5680u;
   y ^= (y << 15) & 0xEFC60000u;
   return (y ^ (y >> 18));
}

inline
int MersenneTwister::operator() (int max)
{
   return int((u64(operator()()) * max) >> 32);
}


/**
 *  getReal()
 *
 *  Returns a random real from [0,1]
 */

inline
real MersenneTwister::getReal()
{
   return (real(operator()()) + real(0.5)) * RESOLUTION;
}

} // namespace HIntLib

#undef HINTLIB_SDN

#endif

