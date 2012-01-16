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

#ifndef HINTLIB_ADAPTINTEGRATOR_H
#define HINTLIB_ADAPTINTEGRATOR_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/rulebasedintegrator.h>

namespace HIntLib
{
   class EmbeddedRuleFactory;

   class AdaptIntegrator: public EmbeddedRuleBasedIntegrator
   {
   public:
      AdaptIntegrator (const EmbeddedRuleFactory *fac)
         : EmbeddedRuleBasedIntegrator(fac),
           minNumEval(0),
           minPercEval(0),
           numInitialRegions(0),
           numEvalInitialRegions(0),
           percEvalInitialRegions(0)
         {}

      virtual
      Status integrate (
         Integrand &f, const Hypercube &h, Index maxEval,
         real reqAbsError, real reqRelError, EstErr &ee);

      // Set a minimum number of integrand evaluations

      void setMinNumEval (Index n)  { minNumEval = n; }
      void setMinPercEval (double);

      // Set the number of (points spent in the) initial regions which are
      // created by subdividing the integration domain in a regular way.

      void setNumInitialRegions (Index n)    { numInitialRegions = n; }
      void setNumEvalInitialRegions (Index n) { numEvalInitialRegions = n; }
      void setPercEvalInitialRegions (double);
      
   private:
      Index minNumEval;
      double minPercEval;

      Index numInitialRegions;
      Index numEvalInitialRegions;
      double percEvalInitialRegions;
   };
}  // namespace HIntLib

#endif

