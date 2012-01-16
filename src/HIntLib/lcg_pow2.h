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
 *  LCG_Pow2
 *
 *  Linear Congruential Generator using the formula
 *
 *      x_{n+1} = x_n * a + c    (mod m=2^e)
 *
 *  Two distinct cases have to be considered:
 *
 *  **** c relative prime to m, i.e. c is odd
 *
 *  The maximal periode is is m, which is achieved iff   a mod 4 = 1
 *
 *  **** c = 0
 *
 *  The maximal achievable periode of this generator is
 *
 *     1        for e = 1
 *     2        for e = 2
 *     2^(e-2)  for e > 2
 *
 *  i.e., for common values of e only one quarter of of the available numbers
 *  are used.
 *  i.e., the two lower order bits of x are constant
 *
 *  This maximal periode is achieved if
 *     i)  x_0 relative prim to m, i.e. x_0 is odd
 *    ii)  a is a primitive element modulo m
 *
 *  Seeding the generator with an odd start value is guaranteed by the program.
 *  Supplying a primitive element as multiplier, however, is up to the user.
 *
 *  a is a primitive element modulo 2^e, if
 *     i)  e = 1  and a is odd
 *    ii)  e = 2  and  a mod 4 = 3
 *   iii)  e = 3  and  a mod 8 = 3, 5, or 7
 *    iv)  e > 3  and  a mod 8 = 3 or 5
 *
 *  Therefore, most common multipliers have  a mod 8 = 5, which makes them
 *  usable for c=0 as well as c odd.
 *
 *  a should be between 0.01 m and 0.99 m.
 *
 *  In addition to that, a should pass the spectral test.
 */

#ifndef LCG_POW2_H
#define LCG_POW2_H 1
 
#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/mymath.h>
#include <HIntLib/exception.h>


namespace HIntLib
{

/*
 *  T  used to hold x_n.  Should be either unsigned or u64 with at least e bits
 *  a  the multiplier
 *  e  defines the modulus m = 2^e
 */

template<class T, T a, unsigned e, T c = T(1)>
class LCG_Pow2
{
private:
   T state;

   // mask contains exactly  e  1-digits

   static const T mask;

   // This function is defend nowhere.
   // It is called for invalid template arguments, causing a compile-time error

   void illegalArgument();

public:
 
   // Types
 
   typedef T ReturnType;

   // Constants

   static const T A = a;
   static const T C = c;

   static const T MAX;
   static const real RANGE;
   static const real RESOLUTION;

   // Constants
 
   // Information about the generator
 
   static T getMax()  { return ~T() >> (std::numeric_limits<T>::digits - e); }
   static real getRange()      { return RANGE; }
   static real getResolution() { return RESOLUTION; }
 
   static size_t getStateSize() { return sizeof(T); }
 
   // Return a random number
 
   T operator() ();
   int operator() (int);
   real getReal();
 
   // Save and restore state
 
   void saveState (void *p) const    { *static_cast<T*>(p) = state; }
   void restoreState (const void *p) { state = *static_cast<const T*>(p); }
 
   // Initialize Generator

   void init (unsigned start)  { state = odd(c)  ?   start  :  2 * start + 1; }

   // Constructors
 
   LCG_Pow2 (unsigned start = a / 3, bool force = false);
};

template<class T, T a, unsigned e, T c>
const T LCG_Pow2<T,a,e,c>::mask = ~T() >> (std::numeric_limits<T>::digits - e);

template<class T, T a, unsigned e, T c>
const T LCG_Pow2<T,a,e,c>::MAX  = ~T() >> (std::numeric_limits<T>::digits - e);

template<class T, T a, unsigned e, T c>
const real LCG_Pow2<T,a,e,c>::RANGE = real(mask) + real(1);

template<class T, T a, unsigned e, T c>
const real LCG_Pow2<T,a,e,c>::RESOLUTION = real(1) / RANGE;


template<class T, T a, unsigned e, T c>
inline
LCG_Pow2<T,a,e,c>::LCG_Pow2 (unsigned start, bool force)
{
   // Do a number of sanity checks.  All this is optimized away at compile time

   if (std::numeric_limits<T>::is_signed)  illegalArgument();
   if (e > unsigned (std::numeric_limits<T>::digits))  illegalArgument();
   if (a >= ~T() >> (std::numeric_limits<T>::digits - e))  illegalArgument();

   if (! force)
   {
      if (c == 0)
      {
         switch (e)
         {
         case 0:  illegalArgument(); break;
         case 1:  if (even(a)) illegalArgument(); break;
         case 2:  if (a % 4 != 3) illegalArgument(); break;
         case 3:  if (a%8 != 3 && a%8 != 5 && a%8 != 7) illegalArgument(); break;
         default: if (a % 8 != 3 && a % 8 != 5) illegalArgument();
         }
      }
      else if (odd (c))
      {
         if (a % 4 != 1) illegalArgument();
      }
      else illegalArgument();
   }

   // a should not be too large or too small.

   // if ((a < MAX/100) || (MAX-a < MAX/100))  illegalArgument();

   init(start);
}


/**
 *  operator()
 */

template<class T, T a, unsigned e, T c>
inline
T LCG_Pow2<T,a,e,c>::operator() ()
{
   // Heuristic to figure out if (hopefully compiletime) recalculation of
   // mask is faster than memory lookup

   if (std::numeric_limits<T>::digits <= std::numeric_limits<unsigned>::digits)
   {
      T result = state & (~T() >> (std::numeric_limits<T>::digits - e));
      state = state * a + c;
      return result;
   }
   else
   {
      T result = state & mask;
      state = state * a + c;
      return result;
   }
}

/**
 *  getReal()
 *
 *  Returns a random real from [0,1]
 */
 
template<class T, T a, unsigned e, T c>
inline
real LCG_Pow2<T,a,e,c>::getReal()
{
   if (odd (c))
      return (real(operator()()) + real(0.5)) * RESOLUTION;
   else
   {
      return real(operator()()) * RESOLUTION;
   }
}
 
/**
 *  operator() (int)
 *
 *  Returns a number from {0,1,...,upperBound-1}
 *
 *  If the product of unsigned (the function parameter) and x fits into a u64,
 *  we use integer arithmetic.  If not, we fall back to floating point.
 *
 *  At least GCC is able to determine the correct algorithm at compile time.
 */

template<class T, T a, unsigned e, T c>
inline
int LCG_Pow2<T,a,e,c>::operator() (int upperBound)
{
   if (int(e) + std::numeric_limits<unsigned>::digits
             <= std::numeric_limits<u64>::digits)
   {
      return int ((u64(operator()()) * upperBound) >> e);
   }
   else
   {
      return int (getReal() * upperBound);
   }
}


/**
 *  Some common generators of this type
 *
 *  See D.E.Knuth, The Art of Computer Programming, Vol II, third ed, 3.3.4
 */

// The one propsed in the ANSI C reference
// DO NOT USE !!!

typedef LCG_Pow2<u32,1103515245,32,12345> LCG_Pow2_Ansi_C;

// A reminder of the good old days, when 35 was a comman word-size.
// Due to O. Taussky.  Line 11 in Knuth's Spectral Test list
// DO NOT USE !!!

typedef LCG_Pow2<u64,1220703125,35,0> LCG_Pow2_Taussky_0;
typedef LCG_Pow2<u64,1220703125,35,1> LCG_Pow2_Taussky;

// The infamous RANDU.  One of the worst LCGs that have ever been used
// Line 12 in Knuth's Spectral Test list
// DO NOT USE !!!

typedef LCG_Pow2<u32,65539u,31,0> LCG_Pow2_RANDU_0;
// typedef LCG_Pow2<u32,65539u,31,1> LCG_Pow2_RANDU;  // a is invalid for c=1

// Two suggested multipliers for 2^32
// Lines 13 and 14 in Knuth's Spectral Test list

typedef LCG_Pow2<u32,1812433253,32,0> LCG_Pow2_BoroshNiederreiter_0;
typedef LCG_Pow2<u32,1812433253,32,1> LCG_Pow2_BoroshNiederreiter;
typedef LCG_Pow2<u32,1566083941,32,0> LCG_Pow2_Waterman_0;
typedef LCG_Pow2<u32,1566083941,32,1> LCG_Pow2_Waterman;

// This one is very common and actually not too bad.
// G. Marsaglia calls it "a candidate for the best of all multipliers"!
// (Line 15)
// It's used to seed the state-array of MersenneTwister
// DO NOT USE !!!

typedef LCG_Pow2<u32,69069,32,0> LCG_Pow2_69069_0;
typedef LCG_Pow2<u32,69069,32,1> LCG_Pow2_69069;

// Two multiplies due to M. Lavaux and F. Janssens found in a computer search
// for spectrally good multipliers having a very high mu_2 (Lines 16 and 23)

typedef LCG_Pow2<u32,1664525,32,0> LCG_Pow2_LavauxJanssens32_0;
typedef LCG_Pow2<u32,1664525,32,1> LCG_Pow2_LavauxJanssens32;
typedef LCG_Pow2<u64,31167285ull,48,0>  LCG_Pow2_LavauxJanssens48_0;
typedef LCG_Pow2<u64,31167285ull,48,1>  LCG_Pow2_LavauxJanssens48;

// This one is used in the CRAY X-MP library (Line 22)

typedef LCG_Pow2<u64,44485709377909ull,48,0> LCG_Pow2_Cray_0;
typedef LCG_Pow2<u64,44485709377909ull,48,1> LCG_Pow2_Cray;

// An excellent one due to C. E. Haynes for a 64-bit generator (Line 26)

typedef LCG_Pow2<u64,6364136223846793005ull,64,0> LCG_Pow2_Haynes_0;
typedef LCG_Pow2<u64,6364136223846793005ull,64,1> LCG_Pow2_Haynes;

}  // namespace HIntLib

#endif

