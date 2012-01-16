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

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/mersennetwister.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <algorithm>

#ifdef HINTLIB_HAVE_CSTDDEF
  #define HINTLIB_SDN std::
#else
  #define HINTLIB_SDN
#endif

#ifdef HINTLIB_HAVE_CSTRING
  #include <cstring>
  #define HINTLIB_SSN std::
#else
  #include <string.h>
  #define HINTLIB_SSN
#endif

#include <HIntLib/lcg_pow2.h>

namespace L = HIntLib;
using L::real;
using L::u32;

#if HINTLIB_STATIC_WORKS == 1
   const real L::MersenneTwister::RANGE = real(~u32()) + real(1);
   const real L::MersenneTwister::RESOLUTION = real(1) / RANGE;
#else
   real L::MersenneTwister::RANGE;
   real L::MersenneTwister::RESOLUTION;
#endif


/**
 * initOld()
 */

void L::MersenneTwister::initOld (unsigned start)
{
   LCG_Pow2_69069_0 lcg (start);  // Create an LCG for seeding state

   for (u32* s = state; s != state + N; ++s)  *s = lcg();

   reload();

   next = state;
}


/**
 *  init()
 *
 *  Standard initialization routine. It is also used in
 *   - the original Makoto Matsumoto and Takuji Nishimura code
 *   - Wagern's implementation
 */

void L::MersenneTwister::init (unsigned start)
{
   LCG_Pow2_69069 lcg (start);

   for (u32* s = state; s != state + N; ++s)
   {
      *s  = lcg() & 0xffff0000ul;
      *s |= lcg() >> 16;
   }

   reload();

   next = state;
}

namespace {

   const size_t M = 397;          // a period parameter

   // some funny bit-shufflers used by reload()

   inline u32 hiBit  (u32 x)         { return x & 0x80000000u; }
   inline u32 loBit  (u32 x)         { return x & 0x00000001u; }
   inline u32 loBits (u32 x)         { return x & 0x7FFFFFFFu; }
   inline u32 mixBits (u32 x, u32 y) { return hiBit (x) | loBits (y); }

   inline u32 twist (u32 m, u32 s0, u32 s1)
      { return m ^ (mixBits(s0,s1) >> 1) ^ (loBit(s1) ? 0x9908B0DFul : 0ul); }
}


/**
 *  reload()
 */

void L::MersenneTwister::reload()
{
         u32* p0 = state;
   const u32* p2 = state + 2;
   const u32* pM = state + M;

   u32 s0 = state[0];
   u32 s1 = state[1];

   for (unsigned j = (N-M)/2+1; --j;)
   {
      *p0++ = twist(*pM++, s0, s1); s0=s1; s1=*p2++;
      *p0++ = twist(*pM++, s0, s1); s0=s1; s1=*p2++;
   }

   *p0++ = twist(*pM, s0, s1); s0=s1; s1=*p2++;    // same as above

   pM = state;

   for (unsigned j = (M-1)/2+1; --j;)
   {
      *p0++ = twist (*pM++, s0, s1); s0=s1; s1=*p2++;
      *p0++ = twist (*pM++, s0, s1); s0=s1; s1=*p2++;
   }

   *p0 = twist (*pM, s0, state[0]);
}


/**
 *  saveState()
 */

void L::MersenneTwister::saveState (void *p) const
{
   HINTLIB_SSN memcpy (p, state, N * sizeof(u32));

   HINTLIB_SDN ptrdiff_t left = (state + N) - next;

   HINTLIB_SSN memcpy (static_cast<char*>(p) + sizeof(state),
                       &left, sizeof (HINTLIB_SDN ptrdiff_t));
}


/**
 * restoreState()
 */

void L::MersenneTwister::restoreState (const void *p)
{
   HINTLIB_SSN memcpy (state, p, N * sizeof(u32));

   HINTLIB_SDN ptrdiff_t left;

   HINTLIB_SSN memcpy (&left, static_cast<const char*>(p) + sizeof(state),
                       sizeof(HINTLIB_SDN ptrdiff_t));

   next = (state + N) - left;
}

/**
 *  Constructor
 */

L::MersenneTwister::MersenneTwister (unsigned start)
   : next(state)
{
   #if HINTLIB_STATIC_WORKS == 0
      RANGE = real(~u32()) + real(1);
      RESOLUTION = real(1) / RANGE;
   #endif
   init(start);
}

/**
 *  Copy-Constructor
 */

L::MersenneTwister::MersenneTwister (const MersenneTwister &mt)
  : next ((state + N) - ((mt.state + N) - mt.next))
{
   #if HINTLIB_STATIC_WORKS == 0
      RANGE = real(~u32()) + real(1);
      RESOLUTION = real(1) / RANGE;
   #endif

   std::copy (mt.state, mt.state + N, state);
}


/**
 *  Assignment
 */

L::MersenneTwister& L::MersenneTwister::operator= (const MersenneTwister &mt)
{
   if (this != &mt)
   {
      std::copy (mt.state, mt.state + N, state);

      next = (state + N) - ((mt.state + N) - mt.next);
   }

   return *this;
}

#undef HINTLIB_SDN
#undef HINTLIB_SSN


