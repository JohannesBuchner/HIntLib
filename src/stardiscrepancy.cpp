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

#define HINTLIB_LIBRARY_OBJECT

#include <algorithm>

#include <HIntLib/discrepancy.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/hlmath.h>
#include <HIntLib/exception.h>
#include <HIntLib/array.h>

namespace L = HIntLib;

using std::max;
using L::real;
using L::Index;


namespace
{

/*
 *  starDiscrepancyDim1()
 *
 *  Calculates the star-discrepancy of n points in [0,1]^1.
 *
 *  We use the formula in Niederreiter 1992, Theorem 2.6.
 *
 *  Runtime is  O(n log n)  for sorting plus
 *              O(n)  for finding the maximum.
 */

real starDiscrepancyDim1 (const real* p, Index n)
{
   L::Array<real> points (p, n);
   std::sort (points.begin(), points.begin() + n);

   if (points[0] < 0 || points[n - 1] > 1)  throw L::FIXME (__FILE__, __LINE__);

   real m = 0;
   const real twoN = 2 * n;

   for (Index i = 0; i < n; ++i)
   {
      m = max (m, points[i] - real (2 * i + 1) / twoN);
   }

   return 1 / twoN + m;
}


/**
 *  Point2
 *
 *  A point in R^2
 */

struct Point2
{
   Point2 (real _x, real _y) : x(_x), y(_y) {}
   Point2 () : x(0), y(0) {}
   real x;
   real y;
};


/**
 *  used for sorting by 1) x-coordinate and 2) y-coordinate
 */

bool operator< (const Point2& p1, const Point2& p2)
{
   if (p1.x < p2.x)  return true;
   if (p1.x > p2.x)  return false;
   return p1.y < p2.y;
}


real starDiscrepancyDim2 (const real* p, Index n)
{
   L::Array<Point2> points (n + 1);

   for (Index i = 0; i < n; ++i)  points[i] = Point2 (p[2 * i], p[2 * i + 1]);
   points [n] = Point2 (1, 1);   // Add (1,1) as last point
   std::sort (points.begin(), points.begin() + n);

   real disc = 0;

   // Create all relevant cubes. Avoid duplicates

   real lastX = -1;
   for (Index i = 0; i <= n; ++i)
   {
      const real x = points[i].x;

      if (x <= lastX)  continue;
      lastX = x;
      real lastY = -1.;

      for (Index j = 0; j <= n; ++j)
      {
         const real y = points[j].y;

         if (y <= lastY)  continue;
         lastY = y;

         // Volume of the cube

         const real volume = x * y;

         // Count points inside the cube

         Index countOpen = 0;
         Index countClosed = 0;

         // for  k > i

         for (Index k = i; k < n; ++k)
         {
            if (points[k].x > x || points[k].y > y)  break;
            ++countClosed;
         }
         
         // for  k <= i  with  points[k].x == x

         int k = i-1;
         for (; k >= 0; --k)
         {
            if (points[k].x < x)  break;
            ++countClosed;
         }

         // for  k < i  with  points[k].x < x

         for (; k >= 0; --k)
         {
            if (points[k].y <= y)
            {
               ++countClosed;
               if (points[k].y < y)  ++countOpen;
            }
         }
         
         // Update discrepancy if a new maximal value was found

         disc = max (disc,
                  max (HINTLIB_MN abs(volume - (real (countClosed) / n)),
                       HINTLIB_MN abs(volume - (real (countOpen)   / n))));
      }
   }

   return disc;
}

}  // anonymous namespace


real
L::starDiscrepancy (const real* points, unsigned dim, Index n)
{
   checkDimensionNotZero (dim);
   checkDimensionLeq<2> (dim);

   return (dim == 1)
      ? starDiscrepancyDim1 (points, n)
      : starDiscrepancyDim2 (points, n);
}

