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
 *  AdaptIntegrator
 *
 */

#ifdef __GNUG__
#pragma implementation
#endif

#include <memory>

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/fallback_limits.h>
#endif

#include <HIntLib/adaptintegrator.h>

#include <HIntLib/embeddedrule.h>
#include <HIntLib/regioncollection.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;

L::Integrator::Status L::AdaptIntegrator::integrate (
   Function &f, const Hypercube &h, Index max,
   real reqAbsError, real reqRelError, EstErr &ee) 
{
   checkDimension (h, f);

   std::auto_ptr<EmbeddedRule> rule (getRule (h.getDimension()));

   // Determine number of iterations

   Index maxEvaluations;

   if (max)
   {
      maxEvaluations = max / rule->getNumPoints ();

      if (maxEvaluations < 1)
      {
         #ifdef HINTLIB_NO_EXCEPTIONS
            ee.set (0.0, 0.0);
            return ERROR;
         #else
            throw NoEvaluationsPossible(max); 
         #endif
      }

      maxEvaluations = (maxEvaluations - 1) / 2;
   }
   else
   {
      maxEvaluations = std::numeric_limits<Index>::max();
   } 

   // Create a Region Collection

   RegionCollection rc;

   // Push whole region

   rc.push (new Region (h, f, *rule));

   // Do refinements

   Status status = ERROR;

   for (Index i = 0; i != maxEvaluations; ++i)
   {
      status = checkRequestedError (rc.getEstErr (), reqAbsError, reqRelError);

      if (status != ERROR)  break;

      rc.refine (f, *rule);
   }

   rc.killQueueAndUpdateResult ();

   ee = rc.getEstErr ();

   return (status == ERROR) ? MAX_EVAL_REACHED : status;
}

