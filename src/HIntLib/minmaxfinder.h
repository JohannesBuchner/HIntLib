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
 *  MinFinder
 *  MaxFinder
 *  MinMaxFinder
 *
 *  Determines the minimum/maximum of a sequence of values
 */

#ifndef HINTLIB_MIN_MAX_FINDER_H
#define HINTLIB_MIN_MAX_FINDER_H 1
 
#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#ifdef HINTLIB_HAVE_LIMITS
#  include <limits>
#else
#  include <HIntLib/fallback_limits.h>
#endif

namespace HIntLib
{

/**
 *  getUpperBound()
 *
 *  Returns a number that can be used as a starting value for searching for a
 *  minimum.
 *  It is the largest number of this type, or infinity, it it can be expressed.
 */
 
template <class T>
inline
T getUpperBound (const T&)
{
   return std::numeric_limits<T>::has_infinity
      ? std::numeric_limits<T>::infinity()
      : std::numeric_limits<T>::max();
}
 
 
/**
 *  getLowerBound()
 *
 *  Returns a number that can be used as a starting value for searching for a
 *  maximum.
 *  It is the smallest number of this type, of -infinity, if it can be
 *  expressed.
 */
 
template <class T>
inline
T getLowerBound (const T&)
{
   if (! std::numeric_limits<T>::is_signed)
      return T(0);
   else if (std::numeric_limits<T>::has_infinity)
      return -std::numeric_limits<T>::infinity();
   else if (std::numeric_limits<T>::is_integer)
      return std::numeric_limits<T>::min();
   else
      return -std::numeric_limits<T>::max();
}


/**
 *  MinFinder
 *
 *  Determines the minimum of a sequence
 */

template<class T>
class MinFinder
{
public:
   MinFinder() : minimum(getUpperBound(T(0))) {}
   MinFinder(T x)  : minimum(x) {}

   bool operator<< (T x)
   {
      if (x < minimum)  { minimum = x; return true; }
      return false;
   }

   T getMinimum() const  { return minimum; }

   void reset()     { minimum = getUpperBound(T(0)); }
   void reset(T x)  { minimum = x; }

private:
   T minimum;
};


/**
 *  MaxFinder
 *
 *  Determines the maximum of a sequence
 */

template<class T>
class MaxFinder
{
public:
   MaxFinder() : maximum(getLowerBound(T(0))) {}
   MaxFinder(T x)  : maximum(x) {}

   bool operator<< (T x)
   {
      if (x > maximum)  { maximum = x; return true; }
      return false;
   }

   T getMaximum() const  { return maximum; }

   void reset()     { maximum = getLowerBound(T(0)); }
   void reset(T x)  { maximum = x; }


private:
   T maximum;
};


/**
 *  MinMaxFinder
 *
 *  Determines the minimum and maximum of a sequence
 */

template<class T>
class MinMaxFinder : public MinFinder<T>,
                     public MaxFinder<T>
{
public:
   MinMaxFinder() {}
   MinMaxFinder(T x)  : MinFinder<T>(x), MaxFinder<T>(x) {}

   void operator<< (T x)
      { MinFinder<T>::operator<<(x); MaxFinder<T>::operator<<(x); }

   T getRange() const  { return this->getMaximum() - this->getMinimum(); }

   void reset()     { MinFinder<T>::reset();  MaxFinder<T>::reset(); }
   void reset(T x)  { MinFinder<T>::reset(x); MaxFinder<T>::reset(x); }
};

}  // namespace HIntLib

#endif

