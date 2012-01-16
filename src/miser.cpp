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

#include <HIntLib/miser.h>

#include <HIntLib/mymath.h>
#include <HIntLib/statistic.h>
#include <HIntLib/hypercube.h>
#include <HIntLib/integrand.h>
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

class MiserImp : public Job
{
public:

   MiserImp (Miser &, PointSet*, PointSet*, unsigned dim, Integrand &f);

   void recur (Hypercube &h, Index numPoints, EstErr &ee);
   void operator() (const real*);

private:

   Miser &miser;
   PointSet* presamplePointSet;
   PointSet* samplePointSet;

   const unsigned dim;

   Integrand &f;
   Point point;
   Array<StatisticMinMax<real> > leftStatistic, rightStatistic;

   const real* center;
   unsigned nextDefaultSplitDim;
};


/*
 *  Initialization
 */

MiserImp::MiserImp (Miser &i, PointSet* presample, PointSet* sample,
                    unsigned dim, Integrand &f)
   : miser(i),
     presamplePointSet (presample), samplePointSet (sample),
     dim(dim), f(f),
     point(dim),
     leftStatistic (dim),
     rightStatistic (dim),
     nextDefaultSplitDim (0)
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

      samplePointSet->setCube (&h);
      samplePointSet->integrate (point, f, numPoints, stat);

      ee.set (stat.getMean(), stat.getVariance() / stat.getCount());

      return;
   }

   // Determine number of presampling points

   Index numPointsPre
      = min(miser.MAX_PRESAMPLING_POINTS,
            max(miser.MIN_PRESAMPLING_POINTS,
                Index(numPoints * miser.PRESAMPLING_RATE)));

   // Make sure, an even number of points remains

   if (odd (numPoints - numPointsPre))  ++numPointsPre;

   // Reset MinMaxFinders

   for (unsigned d = 0; d != dim; ++d)
   {
      leftStatistic [d].reset();
      rightStatistic[d].reset();
   }

   // Do presampling

   center = h.getCenter();
   presamplePointSet->setCube (&h);
   presamplePointSet->doJob (point, *this, numPointsPre);

   // Choose dimension for bisection

   int splitDim = ++nextDefaultSplitDim % dim;
   MinFinder<real> m;
   real  leftVariance = 1.0;
   real rightVariance = 1.0;

   for (unsigned d = 0; d != dim; ++d)
   {
      if (   leftStatistic [d].getCount() >= 10
          && rightStatistic[d].getCount() >= 10
          && (   leftStatistic [d].getRange() > 0
              || rightStatistic [d].getRange() > 0))
      {
         real left  = pow ( leftStatistic [d].getRange(), real (2.0 / 3.0));
         real right = pow (rightStatistic [d].getRange(), real (2.0 / 3.0));
         real sum = left + right;

         if (m << sum)
         {
            splitDim = d;
             leftVariance = left;
            rightVariance = right;
         }
      }
   }

   // Don't allow too large a variance difference

   if ( leftVariance > 50 * rightVariance) rightVariance =  leftVariance / 50;
   if (rightVariance > 50 *  leftVariance)  leftVariance = rightVariance / 50;

   if (leftVariance <= 0.0 || rightVariance <= 0.0)
   {
       leftVariance = 1.0;
      rightVariance = 1.0;
   }
   
   // Calculate point budgets

   numPoints -= numPointsPre + 2 * miser.FREE_POINTS;

   Index numPointsLeft = miser.FREE_POINTS
      + Index(numPoints * leftVariance / (leftVariance + rightVariance));

   if (odd (numPointsLeft))  ++numPointsLeft;

   const Index numPointsRight =
      2 * miser.FREE_POINTS + numPoints - numPointsLeft;

   // do recursion

   Hypercube hRight (h, splitDim);

   EstErr eeLeft;
   recur (h, numPointsLeft, eeLeft);

   EstErr eeRight;
   recur (hRight, numPointsRight, eeRight);

   ee.set ((eeLeft.getEstimate() + eeRight.getEstimate()) / 2.0,
           (eeLeft.getError()    + eeRight.getError()) / 4.0);
}

void MiserImp::operator() (const real* p)
{
   const real value = f (p);

   for (unsigned d = 0; d != dim; ++d)
   {
      if (point[d] <= center[d])  leftStatistic [d] << value;
      else                       rightStatistic [d] << value;
   }
}

}  // end namespace


/**
 *  MiserIntegrator()
 *
 *  Initialize all parameters to defaults
 */

L::Miser::Miser (PointSet* presample, PointSet* sample)
   : presamplePointSet (presample), samplePointSet (sample)
{
   defaults();
}

L::Miser::Miser (PointSet* ps)
   : presamplePointSet (ps), samplePointSet (ps)
{
   defaults();
}

void L::Miser::defaults()
{
   MIN_POINTS = 200;
   FREE_POINTS = 25;
   MIN_PRESAMPLING_POINTS = MIN_POINTS / 4;
   MAX_PRESAMPLING_POINTS = std::numeric_limits<Index>::max();
   // MAX_PRESAMPLING_POINTS = 10000;
   PRESAMPLING_RATE = 0.10;
}


/**
 *  MiserIntegrator::integrate()
 *
 *  Creates an object of type MiserImp and uses it for the actual calculations
 */

L::Miser::Status L::Miser::integrate (
   Integrand &f,
   const Hypercube &h,
   Index maxEval,
   real reqAbsError, real reqRelError,
   EstErr &ee)
{
   checkDimension(h, f);

   // We need an upper bound to the number of sampling points

   if (maxEval == 0)
   {
      #ifdef HINTLIB_NO_EXCEPTIONS
         ee.set (0.0, 0.0);
         return ERROR;
      #else
         throw MaxEvaluationsRequired();
      #endif
   }

   MiserImp imp (*this, presamplePointSet, samplePointSet, h.getDimension(), f);

   Hypercube hh (h);

   imp.recur (hh, maxEval, ee);

   ee.set (ee.getEstimate() * h.getVolume(),
           sqrt(ee.getError()) * h.getVolume());

   Status status = checkRequestedError (ee, reqAbsError, reqRelError);

   return (status == ERROR) ? MAX_EVAL_REACHED : status;
}

