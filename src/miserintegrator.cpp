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
 *  MISER Integrator
 *
 *  This implementation is based on:
 *    [1] W.H. Press, G.R. Farrar. 1990, Computers in Physics, vol. 4,
 *        pp. 190 - 195.
 *    [2] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P.
 *        Flannery. Numerical Recipes in C - The Art of Scientific Computing.
 *        second edition. Cambrdige University Press. Chapter 7.8.
 */

#ifdef __GNUG__
#pragma implementation
#endif

#include <algorithm>

#include <HIntLib/miserintegrator.h>

#include <HIntLib/mymath.h>
#include <HIntLib/statistic.h>
#include <HIntLib/hypercube.h>
#include <HIntLib/function.h>
#include <HIntLib/pointset.h>
#include <HIntLib/distribution.h>
#include <HIntLib/exception.h>

#include <iostream>
#include <iomanip>

using std::min;
using std::max;


namespace L = HIntLib;

namespace {

   using namespace L;

/**
 * MiserImp
 *
 *  One instance of this class is created for every call to
 *  MiserIntegrator::integrate().
 *
 *  Its main purpose is to provide a scope for the "global" variables during
 *  the recursion.
 */

class MiserImp
{
public:

   MiserImp (MiserIntegrator &, PRNG*, PointSet*, unsigned dim, Function &f);

   void recur (Hypercube &h, Index numPoints, EstErr &ee);

private:

   MiserIntegrator &miser;
   PRNG* mc;
   PointSet* ps;

   const unsigned dim;

   Function &f;

   Array<real> point;

   Array<MinMaxFinder<real> > lowerMinMax, upperMinMax;
};


/*
 *  Initialization
 */

MiserImp::MiserImp (MiserIntegrator &i, PRNG* _mc, PointSet* _ps,
                    unsigned dim, Function &f)
   : miser(i),
     mc(_mc), ps(_ps),
     dim(dim), f(f),
     point(dim),
     lowerMinMax (dim),
     upperMinMax (dim)
{}


/**
 *  The main method of MISER Integration
 */

void MiserImp::recur (
   Hypercube &h,
   Index numPoints,
   EstErr &ee)
{
   // If there are too few points left, do plain Monte Carlo

   if (numPoints < miser.MIN_POINTS)
   {
      StatisticVar<real> stat;

      ps->setCube (&h);
      ps->integrate (point, f, numPoints, stat);

      ee.set (stat.getMean(), stat.getVariance() / stat.getCount());

      return;
   }

   // Reset MinMaxFinders

   for (unsigned d = 0; d != dim; ++d)
   {
      lowerMinMax[d].reset();
      upperMinMax[d].reset();
   }

   // Determine number of presampling points

   Index numPointsPre
      = min(miser.MAX_PRESAMPLING_POINTS,
            max(miser.MIN_PRESAMPLING_POINTS,
                Index(numPoints * miser.PRESAMPLING_RATE)));

   // Make sure, an even number of points remains

   if ((numPoints - numPointsPre) % 2)  ++numPointsPre;

   // Do presampling

   UniformCube<PRNG> cubeSampler (*mc, h);

   for (Index i = 0; i != numPointsPre; ++i)
   {
      cubeSampler (*mc, point);

      const real value = f(point);

      for (unsigned d = 0; d != dim; ++d)
      {
         if (point[d] <= h.getCenter(d)) lowerMinMax[d] << value;
         else                            upperMinMax[d] << value;
      }
   }

   // Choose dimension for bisection

   int splitDim = -1;
   MinFinder<real> m;
   real upperVariance = 1.0;
   real lowerVariance = 1.0;

   for (unsigned d = 0; d != dim; ++d)
   {
      if (lowerMinMax[d].getRange() > 0.0 && upperMinMax[d].getRange() > 0.0)
      {
         real lower = pow (lowerMinMax[d].getRange(), 2.0 / 3.0);
         real upper = pow (upperMinMax[d].getRange(), 2.0 / 3.0);
         real sum = lower + upper;

         if (m << sum)
         {
            splitDim = d;
            lowerVariance = lower;
            upperVariance = upper;
         }
      }
   }

   // Did we find a valid dimension?  If not, take a random one

   if (splitDim < 0)  splitDim = mc->equidist (dim);

   // Calculate point budgets

   numPoints -= numPointsPre + 2 * miser.FREE_POINTS;

   Index numPointsLower = miser.FREE_POINTS
      + Index(numPoints * lowerVariance / (lowerVariance+upperVariance));

   if (odd (numPointsLower))  ++numPointsLower;

   const Index numPointsUpper =
      2 * miser.FREE_POINTS + numPoints - numPointsLower;

   // do recursion

   Hypercube hUpper (h, splitDim);

   recur (h, numPointsLower, ee);

   EstErr eeUpper;

   recur (hUpper, numPointsUpper, eeUpper);

   ee.set ((ee.getEstimate() + eeUpper.getEstimate()) / 2.0,
           (ee.getError() + eeUpper.getError()) / 4.0);
}

}  // end namespace


/**
 *  MiserIntegrator()
 *
 *  Initialize all parameters to defaults
 */

L::MiserIntegrator::MiserIntegrator (PRNG* _mc, PointSet* _ps)
   : mc (_mc), ps (_ps),
     MIN_POINTS(200),
     FREE_POINTS(25),
     MIN_PRESAMPLING_POINTS(40),
  // MAX_PRESAMPLING_POINTS(numeric_limits<Index>::max()),
     MAX_PRESAMPLING_POINTS(10000),    // original
     PRESAMPLING_RATE(0.1)
{}


/**
 *  MiserIntegrator::integrate()
 *
 *  Creates an object of type MiserImp and uses it for the actual calculations
 */

L::MiserIntegrator::Status L::MiserIntegrator::integrate (
   Function &f,
   const Hypercube &h,
   Index maxEval,
   real reqAbsError, real reqRelError,
   EstErr &ee)
{
   checkDimension(h, f);

   // We need an upper bound to the number of sampling points

   if (maxEval == 0)
   {
      #ifdef CRAY
         ee.set (0.0, 0.0);
         return ERROR;
      #else
         throw MaxEvaluationsRequired();
      #endif
   }

   MiserImp imp (*this, mc, ps, h.getDimension(), f);

   Hypercube hh (h);

   imp.recur (hh, maxEval, ee);

   ee.set (ee.getEstimate() * h.getVolume(),
           sqrt(ee.getError()) * h.getVolume());

   Status status = checkRequestedError (ee, reqAbsError, reqRelError);

   return (status == ERROR) ? MAX_EVAL_REACHED : status;
}

