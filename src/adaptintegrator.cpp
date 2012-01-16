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
 *  AdaptIntegrator
 */

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <memory>

#ifdef HINTLIB_HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/fallback_limits.h>
#endif

#include <HIntLib/adaptintegrator.h>

#include <HIntLib/cubaturerule.h>
#include <HIntLib/regioncollection.h>
#include <HIntLib/exception.h>
#include <HIntLib/cubepartitioner.h>

#define DEB 0

#if DEB >= 1
#include <iostream>
using std::cout;
#endif

namespace L = HIntLib;

/**
 *  AdaptIntegrator
 *
 *  Strategy:
 *
 *     1) Use a CubePartition to divide the integration domain a predefined
 *        number of subcubes (default 1)
 *     2) Use a RegionCollection to subdivide the region with the largest error
 *        until the precision goal is met or we run out of integrand
 *        evaluations.
 */

namespace {

   using namespace HIntLib;

   /**
    * This CubePartitioner stores each subcube in a Region Collection
    */

   class MyCubePartitioner : public CubePartitioner
   {
   public:
      MyCubePartitioner (RegionCollection& _rc, Integrand& _f, EmbeddedRule& _r)
         : rc (_rc), f(_f), r(_r) {}

      virtual void action (const Hypercube& h);

   private:
      RegionCollection& rc;
      Integrand &f;
      EmbeddedRule &r;
   };

   void
   MyCubePartitioner::action (const Hypercube& h)
   {
      rc.push (new Region (h, f, r));
   }

} // namespace


/**
 * setXXX()
 */

void
L::AdaptIntegrator::setMinPercEval (double perc)
{
   if (perc < 0.0 || perc > 1.0)  throw 1;

   minPercEval = perc;
}

void
L::AdaptIntegrator::setPercEvalInitialRegions (double perc)
{
   if (perc < 0.0 || perc > 1.0)  throw 1;

   percEvalInitialRegions = perc;
}


/**
 * integrate()
 */

L::Integrator::Status
L::AdaptIntegrator::integrate (
   Integrand& f, const Hypercube& h, Index maxEval,
   real reqAbsError, real reqRelError, EstErr& ee)
{
#if DEB >= 1
      cout << "in AdaptIntegrator::integrate ("
           << maxEval << ',' << reqAbsError << ',' << reqRelError << ")\n";
#endif

   checkDimension (h, f);

   // Create Cubature Rule and Region Collection

   std::auto_ptr<EmbeddedRule> rule (getRule (h.getDimension()));
   const Index rulePoints = rule->getNumPoints ();

#if DEB >= 1
   cout << "   Points in rule: " << rulePoints << '\n';
#endif

   // Before we start the calculation, we have to determine three parameters:

   Index initialRegions;  // The number of initial regions (split non-adaptive)
   Index minIter;         // The minimum number of adaptive subdivisions
   Index maxIter;         // The maximum number of adaptive subdivisions

   // How to determine these three parameters depends whether  _maxEval_ is
   // given

   if (maxEval)   // with  maxEval
   {
      // Determine minimum number of evaluations.

      Index minEval = minNumEval;

      if (minPercEval >= 1.)
      {
         minEval = std::max (minEval, maxEval);
      }
      else if (minPercEval > 0.)
      {
         minEval = std::max (minEval, Index (maxEval * minPercEval));
      }

      // Determine number of initial regions

      initialRegions = numInitialRegions;

      {
         Index n = std::max (Index(1), numEvalInitialRegions);

         if (percEvalInitialRegions >= 1.)
         {
            n = std::max (n, maxEval);
         }
         else if (percEvalInitialRegions > 0.)
         {
            n = std::max (n, Index (maxEval * percEvalInitialRegions));
         }

         initialRegions = std::max (initialRegions, (n - 1) / rulePoints + 1);
      }

      // make sure we are allowed to spend that many points on initial regions

      Index initialPoints = initialRegions * rulePoints;

      if (initialPoints > maxEval)
      {
         initialRegions = maxEval / rulePoints;  // round down
         initialPoints  = initialRegions * rulePoints;
      }

      // adjust maxEval and minEval

      maxEval -= initialPoints;

      if (minEval > initialPoints)  minEval -= initialPoints;
      else minEval = 0;

      // Determine minimum and maximum number of iterations

      maxIter = maxEval / (2 * rulePoints);  // round down

      // Determine minimum number of iterations

      minIter = minEval  ?  (minEval - 1) / (2 * rulePoints) + 1 // round up
                         :  0;

      if (minIter > maxIter)  minIter = maxIter;
   }
   else  //  no  maxEval
   {
      // Determine the number of inital regions

      initialRegions = std::max (Index(1), numInitialRegions);

      if (numEvalInitialRegions > 0)
      {
         initialRegions = std::max (
               initialRegions,
               (numEvalInitialRegions - 1) / rulePoints + 1);  // round up
      }

      // Determine the minimum number of iterations

      minIter = minNumEval;

      if (minIter <= rulePoints * initialRegions)  minIter = 0;
      else minIter -= rulePoints * initialRegions;

      if (minIter)
      {
         minIter = (minIter - 1) / (2 * rulePoints) + 1;  // round up
      }
      
      // Set maximum number of iterations to max()

      maxIter = std::numeric_limits<Index>::max();
   }

   // Do initial split

#if DEB >= 1
      cout << "   " << initialRegions << " initial region(s)\n"
           << "   minIter/maxIter = " << minIter << '/' << maxIter << '\n';
#endif

   // We need at least one initial region

   if (initialRegions == 0)
   {
      #ifdef HINTLIB_NO_EXCEPTIONS
         ee.set (0.0, 0.0);
         return ERROR;
      #else
         throw NoEvaluationsPossible(maxEval);
      #endif
   }

   // Do initial split

   RegionCollection rc;

   if (initialRegions == 1)
   {
      rc.push (new Region (h, f, *rule));
   }
   else
   {
      MyCubePartitioner cp (rc, f, *rule);
      cp (h, initialRegions, cp.EQUIDISTANT);
   }

   // Do refinements

   Status status = ERROR;

   for (Index i = 0; i != maxIter; ++i)
   {
      status = checkRequestedError (rc.getEstErr (), reqAbsError, reqRelError);

      if (status != ERROR && i >= minIter)  break;

#if DEB >= 2
      :cout << "      In iteration " << (i + 1) << '\n';
#endif
      
      rc.refine (f, *rule);
   }

   rc.killQueueAndUpdateResult ();

   ee = rc.getEstErr ();

   return (status == ERROR) ? MAX_EVAL_REACHED : status;
}

