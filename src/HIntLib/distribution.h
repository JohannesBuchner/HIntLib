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
 *  distribution.h
 *
 *  Provides a number of probability distributions to be used with PRNGs
 */

#ifndef HINTLIB_DISTRIBUTION_H
#define HINTLIB_DISTRIBUTION_H 1
 
#ifdef __GNUG__
#pragma interface
#endif


#include <HIntLib/hypercube.h>
#include <HIntLib/shiftscale.h>

namespace HIntLib
{

/**
 *  Uniform distributions in [a,b]
 */

template<class PRNG>
inline
real uniform(PRNG &prng)       //  [0,1]
{
   return prng.getReal();
}

template<class PRNG>
inline
real uniform(PRNG &prng, real max)   //   [0,max]
{
   return uniform(prng) * max;
}

template<class PRNG>
inline
real uniform(PRNG &prng, real min, real max)  //   [min,max]
{
   return uniform(prng) * (max - min) + min;
}


/**
 *  Uniform distribution in {a, a+1,...,b-1, b}
 */

template<class PRNG>
inline
int equidist(PRNG &prng, int max)     //  {0,1,...,max-2,max-1}
{
   return prng(max);
}

template<class PRNG>
inline
int equidist(PRNG &prng, int min, int max)  //  {min,min+1,...,max-2,max-1}
{
   return prng(max-min) + min;
}

/**
 *  Uniform distributions in  C^s
 */

template<class PRNG>
inline
void uniform (PRNG &prng, const Hypercube &h, real *p)
{
   const unsigned dim = h.getDimension();

   for (unsigned i = 0; i != dim; ++i)
   {
      p[i] = uniform (prng, h.getLowerBound(i), h.getUpperBound(i));
   }
}


/**
 *  UniformCube
 *
 *  Uniform distribution in C^s
 *
 *  This class should be used if a number of samples have to be drawn from
 *  a given Hypercube.  A linear transformation is precalculated to allow
 *  fast calculation of the sampling point based on the unprocessed integer
 *  output of the random number generator.
 */

template<class PRNG>
class UniformCube
{
public:
   UniformCube (const PRNG &g, const Hypercube &h)
      : dim (h.getDimension()), ss (h, g.getRange())  {}

   void operator() (PRNG &g, real *p)
   {
      for (unsigned i = 0; i != dim; ++i)  p[i] = ss[i] (g());
   }

private:
   const unsigned dim;
   const ShiftScale ss;
};

} // namespace HIntLib

#endif

