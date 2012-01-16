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
 *  LCG_Prime
 *
 *  Linear Congruential Generator using the formula
 *
 *      x_{n+1} = x_n * a    (mod p)   with p prime
 *
 *  The implementation provided here is pretty fast.  However, only certain
 *  values can be used as multiplier:
 *
 *  Once  q = m / a
 *   are  r = m % a have been calculated, r must be less than q.
 *
 *  The maximum period length (p-1) is obtained, iff
 *    i)  x_0  relative primt to p, i.e. x_0 is not 0
 *   ii)  a is a primitive element modulo m, i.e.
 *           -) a != 0
 *           -) a^((p-1)/q) != 1 (mod p)  for any prime divisor q of p-1
 */

#ifndef HINTLIB_LCG_PRIME_H
#define HINTLIB_LCG_PRIME_H 1
 
#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/fallback_limits.h>
#endif


namespace HIntLib
{

/*
 *  T  used to hold x_n.  Should be either unsigned or u64 with at least e bits
 *  a  the multiplier
 *  m  modulus (prime)
 */

template<class T, T a, T m>
class LCG_Prime
{
private:
   T state;

   static const T q = m / a;
   static const T r = m % a;

   // This function is defend nowhere.
   // It is called for invalid template arguments, causing a compile-time error

   void invalidArgument();

public:
 
   // Types
 
   typedef T ReturnType;

   // Constants

   static const T    A = a;
   static const T    M = m;
   static const T    C = 0;
 
   static const T    MAX = m - 1;
   static const real RANGE;
   static const real RESOLUTION;

   // Information about the generator
 
   static T getMax()           { return MAX; }
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

   void init (unsigned start)  { state = start + 1; }

   // Constructors
 
   LCG_Prime (unsigned start = m - a);
};

template<class T, T a, T m>
const real LCG_Prime<T,a,m>::RANGE = real(m);

template<class T, T a, T m>
const real LCG_Prime<T,a,m>::RESOLUTION = real(1) / RANGE;


/**
 *  LCG_Prime()
 */

template<class T, T a, T m>
inline
LCG_Prime<T,a,m>::LCG_Prime (unsigned start)
{
   if (! std::numeric_limits<T>::is_signed)  invalidArgument();
   if (a == 0 || a >= m || r >= q)  invalidArgument();

   // There is no easy way for checking if m is prime, nor for making sure that
   // a is a primitive element (mod p)
 
   init(start);
}


/**
 * operator()
 */

template<class T, T a, T m>
inline
T LCG_Prime<T,a,m>::operator() ()
{
   state = a * (state % q) - r * (state / q);
   if (state < 0)  state += m;
   return state;
}


/**
 *  getReal()
 */

template<class T, T a, T m>
inline
real LCG_Prime<T,a,m>::getReal()
{
   return real(operator()()) * RESOLUTION;
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

template<class T, T a, T m>
inline
int LCG_Prime<T,a,m>::operator() (int upperBound)
{
   if (   std::numeric_limits<T>::digits + std::numeric_limits<unsigned>::digits
       <= std::numeric_limits<u64>::digits)
   {
      return int ((u64(operator()()) * upperBound) / m);
   }
   else
   {
      return int (getReal() * upperBound);
   }
}


/**
 *  LCG Combined
 *
 *  Combines two LCG_Prime generators to one large generator
 */

template<class T, T a1, T m1, T a2, T m2>
class LCG_Combined
{
private:

   typedef LCG_Prime<T,a1,m1> LCG1;
   typedef LCG_Prime<T,a2,m2> LCG2;
   LCG1 lcg1;
   LCG2 lcg2;

   // This function is defend nowhere.
   // It is called for invalid template arguments, causing a compile-time error

   void invalidArgument();

public:
 
   // Types
 
   typedef T ReturnType;

   // Constants

   static const T    M = m1 > m2 ? m1 : m2;
   static const T    MAX = (m1 > m2 ? m1 : m2) - 1;
   static const real RANGE;
   static const real RESOLUTION;

   // Information about the generator
 
   static T getMax()           { return MAX; }
   static real getRange()      { return RANGE; }
   static real getResolution() { return RESOLUTION; }
 
   static size_t getStateSize()
      { return LCG1::getStateSize() + LCG2::getStateSize(); }
 
   // Return a random number
 
   T operator() ();
   int operator() (int);
   real getReal();
 
   // Save and restore state
 
   void saveState (void *p) const
      { lcg1.saveState(p); lcg2.saveState (p+lcg1.getStateSize()); }
   void restoreState (const void *p)
      { lcg1.restoreState(p); lcg2.restoreState (p+lcg1.getStateSize()); }
 
   // Initialize Generator

   void init (unsigned start)
      { lcg1.init (start); lcg2.init (LCG2::M / (start+1));
        operator()(); operator()(); operator()(); }

   // Constructors
 
   LCG_Combined (unsigned start = 0);
};


template<class T, T a1, T m1, T a2, T m2>
const real LCG_Combined<T,a1,m1,a2,m2>::RANGE = real(M);

template<class T, T a1, T m1, T a2, T m2>
const real LCG_Combined<T,a1,m1,a2,m2>::RESOLUTION = real(1) / RANGE;


/**
 *  LCG_Combined()
 */

template<class T, T a1, T m1, T a2, T m2>
inline
LCG_Combined<T,a1,m1,a2,m2>::LCG_Combined (unsigned start)
{
   if (m1 == m2)  invalidArgument();
   init (start);
}


/**
 * operator()
 */

template<class T, T a1, T m1, T a2, T m2>
inline
T LCG_Combined<T,a1,m1,a2,m2>::operator() ()
{
   // T x = lcg1.operator()() - lcg2.operator()();
   T x = lcg2.operator()();
   return (x > 0) ? x : x + MAX;
}


/**
 *  getReal()
 */

template<class T, T a1, T m1, T a2, T m2>
inline
real LCG_Combined<T,a1,m1,a2,m2>::getReal()
{
   return real(operator()()) * RESOLUTION;
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

template<class T, T a1, T m1, T a2, T m2>
inline
int LCG_Combined<T,a1,m1,a2,m2>::operator() (int upperBound)
{
   if (   std::numeric_limits<T>::digits + std::numeric_limits<unsigned>::digits
       <= std::numeric_limits<u64>::digits)
   {
      return int ((u64(operator()()) * upperBound) / M);
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

// The "Minimal Standard Generator" (Park and Miller),
//     proposed by Lewis, Goodman, and Miller 1969
// This one is ok, however a is almost too small
// Used in the IMSL subroutine library since 1971
// Line 19 in TACP
// DO NOT USE !!!

typedef LCG_Prime<long, 7*7*7*7*7, (1ul << 31) - 1> LCG_Prime_ISML;

// The best multiplier for m=(2^31)-1 allowing this implementation technique
// Line 20 in TACP

typedef LCG_Prime<long, 48271, (1ul << 31) - 1> LCG_Prime_Fishman;

// And here is another good one (discussed in Park and Miller)

typedef LCG_Prime<long, 69621, (1ul <<31) - 1> LCG_Prime_69621;

// Another spectral-test success story due to P. L'Ecuyer
// Line 21 in TACP

typedef LCG_Prime<long, 40692, (1ul << 31) - 249> LCG_Prime_Lecuyer;

// The combined generator proposed by L'Ecuyer, using two LCG_Primes from above

typedef LCG_Combined<long, 48271, (1ul << 31) - 1, 40692, (1ul << 31) - 249>
   LCG_Combined_Lecuyer;

}  // namespace HIntLib

#endif

