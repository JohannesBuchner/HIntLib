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
 *  Hypercube
 *
 *  A n-dimensional hypercube.
 *
 *  (C) 2001 by Rudolf Schürer
 */

#if defined __GNUG__ && ! defined PARALLEL
#pragma implementation
#endif

#include <HIntLib/hypercube.h>

#include <HIntLib/mymath.h>

namespace L = HIntLib;

using L::approx;


/*****************************************************************************
 *  Constructors
 *****************************************************************************/

L::Hypercube::Hypercube (unsigned _dim)
   : dim (_dim), data (2 * _dim, 0.5), volume (1.0) {}

L::Hypercube::Hypercube (unsigned _dim, const real a [], const real b [])
   : dim (_dim), data (2 * _dim)
{
   for (unsigned i = 0; i < dim; ++i)
   {
      center (i) =     (b [i] + a [i]) * 0.5;
      width (i)  = abs (b [i] - a [i]) * 0.5;
   }

   calcVolume ();
}

L::Hypercube::Hypercube (unsigned _dim, real a, real b)
   : dim (_dim), data (2 * _dim)
{
   const real c =     (b + a) * 0.5;
   const real w = abs (b - a) * 0.5;

   for (unsigned i = 0; i < dim; ++i)
   {
      center (i) = c;
      width (i)  = w;
   }

   calcVolume ();
}


/*****************************************************************************
 *  Member Functions
 *****************************************************************************/

/**
 *  Assignement operator
 *
 */

L::Hypercube& L::Hypercube::operator= (const Hypercube &h)
{
   if (&h != this)
   {
      checkDimensionEqual (*this, h);

      volume = h.volume;

      data.assign (h.data, 2*dim);
   }

   return *this;
}


/**
 *  Comparison
 */

bool L::Hypercube::operator== (const Hypercube &h) const
{
   for (unsigned i = 0; i < dim; ++i)
   {
      if (! approx (getCenter(i), h.getCenter(i))
       || ! approx (getWidth (i), h.getWidth (i)))  return false;
   }

   return true;
}


/*****************************************************************************
 *  Functions working with Hypercubes
 *****************************************************************************/

/**
 *  Build union of two hypercubes (if the result is a Hypercube)
 *
 *  The union is successful, if both cubes have a common surface
 *
 *  If the union is successful, the first cube contains the union,
 *  and true is returned.
 *
 *  If the union can not be performed, false is returned.
 */

const L::real FACTOR = 10.0;

bool L::unite (Hypercube &h, const Hypercube &hh)
{
  if (h.getDimension() != hh.getDimension())  return false;

  int found = -1;
  bool which = false;

  for (unsigned i = 0; i < h.getDimension(); ++i)
  {
    // Check for identity

    if (   approx (h.getCenter(i), hh.getCenter(i), FACTOR)
        && approx (h.getWidth (i), hh.getWidth (i), FACTOR))
    {
       continue;
    }

    // Dont find a fitting face twice

    if (found >= 0)  return false;

    // check if h is below hh

    if (approx (h.getUpperBound(i), hh.getLowerBound(i), FACTOR))
    {
       found = i;
       which = true;

       continue;
    }

    // Check if h is above hh

    if (approx (h.getLowerBound(i), hh.getUpperBound(i), FACTOR))
    {
       found = i;
       which = false;

       continue;
    }

    // Does not fit

    return false;
  }

  if (found >= 0)
  {
     if (which)
        h.set (found,  h.getLowerBound(found), hh.getUpperBound(found));
     else
        h.set (found, hh.getLowerBound(found),  h.getUpperBound(found));

     return true;
  }
  else
  {
     return false;
  }
}


/**
 *  isUnitCube()
 *
 *  Tests, if a given cube is the standard unit cube [0,1]x...x[0,1].
 */

bool L::isUnitCube (const Hypercube &h)
{
   for (unsigned i = 0; i < h.getDimension(); ++i)
   {
      if (! approx (h.getCenter(i), .5) || ! approx (h.getWidth(i), .5))
         return false;
   }

   return true;
}


/**
 *  isPointInside()
 *
 *  Returns true if a point is strictly inside the cube
 *
 *  If border cased are important, consider using whereIsPoint()
 */

bool L::isPointInside (const Hypercube &h, const real x [])
{
   for (unsigned i = 0; i < h.getDimension(); ++i)
   {
      if (x [i] < h.getLowerBound (i) || x [i] > h.getUpperBound (i))
         return false;
   }

   return true;
}


/**
 *  whereIsPoint
 *
 *  Returns
 *    BORDER,  if the point is close to the border of the cube
 *    INSIDE,  if it is inside this border region
 *    OUTSIDE, if it is outside
 */

L::Hypercube::Location L::whereIsPoint (const Hypercube &h, const real x[])
{

   // XXX eps should be calculated for each dimension

   const real eps
      = std::numeric_limits<real>::epsilon() * h.getDiameter(0) * 3.0;

   for (unsigned i = 0; i < h.getDimension(); ++i)
   {
      if (h.getLowerBound(i) - x[i] > eps || x[i] - h.getUpperBound(i) > eps)
         return Hypercube::OUTSIDE;
   }

   for (unsigned i = 0; i < h.getDimension(); ++i)
   {
      if (x[i] - h.getLowerBound(i) <= eps || h.getUpperBound(i) - x[i] <= eps)
         return Hypercube::BORDER;
   }

   return Hypercube::INSIDE;
}



