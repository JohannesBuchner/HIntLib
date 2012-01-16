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
 *  CubePartitioner
 *
 *  Splits a Hypercube into a given number of (almost) equally sized
 *  sub-Hypercubes and performs an action on each of them
 */

#ifdef __GNUG__
#pragma implementation
#endif

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/cubepartitioner.h>

#include <HIntLib/hypercube.h>
#include <HIntLib/hlmath.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;

namespace {

inline
int getNumSplits (
   unsigned k,             // dimension
   L::Index n,             // remaining points
   L::Index avgNumSplits,  // approximated number of splits
   const L::Array<L::Index>& splitThresholds)
{
   // By default, choose avgNumSplits
   // Usue avgNumSplits+1 only if using avgNumSplits is not enough to
   //    to split remaining dimensions avgNumSplit+1 times.
   // Using the lower split values in the first dimensions speeds up the
   //    algorithm.

   return (splitThresholds [k] >= n) ? avgNumSplits : avgNumSplits + 1;
}

} // namespace


namespace L = HIntLib;

void L::CubePartitioner::operator() (
   const Hypercube &h, Index numSections, DistanceType distanceType)
{
   const int dim = h.getDimension ();

   if (numSections < 1 || dim < 1)  return;

   // Hypercube hh represents the current section.

   Hypercube hh (h);

   // The average number of splits along each axis
   // We perform exactly avgNumSplits or avgNumSplits+1 splits

   Index avgNumSplits = Index (HINTLIB_MN floor (
         HINTLIB_MN pow (double (numSections), 1.0 / dim)));

   // do not underestimate this number due to rounding errors in pow()

   if (   HINTLIB_MN pow(double(avgNumSplits+1), dim)
       <= std::numeric_limits<Index>::max())
   {
      if (powInt (avgNumSplits+1, dim) <= numSections)  avgNumSplits++;
   }

   // Create a table containing (avgNumSplits+1)^(dim-i-1) * avgNumSplits
   // Used to determine if a coordinate has to be split avgNumSplit or
   //    avgNumSplit+1 times. (see getNumSplits())

   Array<Index> splitThresholds (dim);

   {
      Index x = avgNumSplits;

      for (int i = 0; i != dim; ++i)
      {
         splitThresholds [dim-i-1] = x;

         // make sure we do not overflow

         if (   double(x) * double(avgNumSplits+1) + .5
             >= double(std::numeric_limits<Index>::max()))
         {
            for (int j = i+1; j != dim; ++j)
               splitThresholds [dim-j-1] = std::numeric_limits<Index>::max();

            break;
         }

         // multiply x

         x *= (avgNumSplits + 1);
      }
   }

   // XXX Right now we do NOT use Index consistently.  Therefor we have to
   // back off if avgNumSplits is too large

   if (avgNumSplits >= Index(std::numeric_limits<int>::max()))
   {
      throw InternalError (__FILE__, __LINE__);
   }

   // # of direct splits in a certain dimension

   Array<int> numSplits (dim, 0);
   numSplits [0] = getNumSplits (0, numSections,
                                 avgNumSplits, splitThresholds);

   // Initialize coordinate 0 of hh, to make hh.move() down below work

   hh.set (0, h.getLowerBound (0) - h.getDiameter (0) / numSplits [0],
              h.getLowerBound (0));

   // Active split for each dimension

   Array<int> curSplit (dim, -1);

   // Total number of siblings
   //   curNumLeaves [0] = n
   //   curNumSubsectsion [dim] = 1

   Array<int> curNumLeaves (dim+1);
   curNumLeaves [0] = numSections;
   curNumLeaves [dim] = 1;

   // Number of subcubes already processed on this level
   // Used to calculate split-plane positions for EQUIVOLUME

   Array<int> curNumProcessedLeaves (dim);
   curNumProcessedLeaves [0] = 0;


   // Main loop

   for (;;)
   {
      // k stores the level in the tree we are working right now
      // k=0 means root
      // k=dim-1 means leave nodes

      // Search for the level where we have to switch branch

      int k = dim - 1;

      while (k >= 0 && ++curSplit [k] == numSplits [k])  k--;

      // Terminate when the whole tree is processed

      if (k < 0)  break;

      // Do initialization at all sub-levels k+1...

      for (int i = k + 1; i < dim; ++i)
      {
         // Calculate total number of subsections at this branch

         curNumLeaves [i] = curNumLeaves [i-1] / numSplits [i-1];

         // Adjust, if parent numSections is not divisible by numSplits

         if (curNumLeaves [i-1] % numSplits [i-1] > curSplit [i-1])
            curNumLeaves [i] ++;

         // Determine number of splits

         numSplits [i] = getNumSplits (i, curNumLeaves [i],
                                 avgNumSplits, splitThresholds);

         // Start with first branch

         curSplit [i] = curNumProcessedLeaves [i] = 0;
      }

      // Set coordinates in hh correctly
      // coordinate k is incremented
      // coordinates k... are reset to the lower border

      if (distanceType == EQUIDISTANT)
      {
         hh.move (k, h.getDiameter (k) / numSplits [k]);

         for (int i = k + 1; i < dim; ++i)
         {
             hh.set (i,
                h.getLowerBound (i),
                h.getLowerBound (i) + h.getDiameter (i) / numSplits [i]);
         }
      }
      else  // distanceType == EQUIVOLUME
      {
         hh.set (k,
            h.getLowerBound (k) + h.getDiameter (k)
              *  curNumProcessedLeaves[k]
                         / curNumLeaves[k],
           h.getLowerBound (k) + h.getDiameter (k)
              * (curNumProcessedLeaves[k] + curNumLeaves[k+1])
                         / curNumLeaves[k]);

         curNumProcessedLeaves [k] += curNumLeaves [k+1];

         for (int i = k + 1; i < dim; ++i)
         {
            hh.set (i,
               h.getLowerBound (i),
               h.getLowerBound (i)
             + h.getDiameter (i) * curNumLeaves[i+1] / curNumLeaves[i]);

            curNumProcessedLeaves [i] += curNumLeaves [i+1];
         }
      }

      // Performe job on this sub-cube

      action (hh);
   }
}


