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
 *  CompRuleIntegrator
 *  CompRuleIntegratorEst
 *
 *  This integrator splits the cube into a number of subcubes and applies a
 *  basic rule to each of them.
 *
 *  The *Est version used EmbeddedRules and produeces an error estimate.
 */

#ifdef __GNUG__
#pragma implementation
#endif

#include <memory>

#include <HIntLib/compruleintegrator.h>

#include <HIntLib/cubaturerule.h>
#include <HIntLib/cubepartitioner.h>
#include <HIntLib/kahanadd.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;
using L::Integrator;


namespace {

   using namespace HIntLib;

   class MyCubePartitioner : public CubePartitioner
   {
   public:
      MyCubePartitioner (Integrand &f, CubatureRule &rule)
         : f (f), rule (rule), estimate (0.0) {}

      virtual void action (const Hypercube &h);

   private:
      Integrand &f;
      CubatureRule &rule;

   public:
      KahanAdd estimate;
   };

   void MyCubePartitioner::action (const Hypercube &h)
   {
      estimate += rule.eval (f, h);
   }

} // namespace


Integrator::Status L::CompRuleIntegrator::integrate (
   Integrand &f, const Hypercube &h, Index maxEval, real, real, EstErr &ee)
{
   checkDimension(h, f);

   if (! maxEval)
   {
      #ifdef HINTLIB_NO_EXCEPTIONS
         ee.set (0.0, 0.0);
         return ERROR;
      #else
         throw MaxEvaluationsRequired();
      #endif
   }

   // Calculate the number of sections 

   std::auto_ptr<CubatureRule> rule (getRule(h.getDimension()));

   Index numSections = maxEval / rule->getNumPoints();

   // We need at least one section to apply the basic rule once

   if (numSections == 0)
   {
      #ifdef HINTLIB_NO_EXCEPTIONS
         ee.set (0.0, 0.0);
         return ERROR;
      #else
         throw NoEvaluationsPossible (maxEval);
      #endif
   }

   // Use CubePartitioner to perform action() an every sub-cube

   MyCubePartitioner cp (f, *rule);

   cp (h, numSections, cp.EQUIDISTANT);

   // Copy to final result

   ee.setNoErr (cp.estimate);

   return MAX_EVAL_REACHED;
}


/****************** With estimate ******************/

namespace {

   using namespace HIntLib;

   class MyErrCubePartitioner : public CubePartitioner
   {
   public:
      MyErrCubePartitioner (Integrand &f, EmbeddedRule &rule)
         : f (f), rule (rule), estimate (0.0), error (0.0) {}
 
      virtual void action (const Hypercube &h);
 
   private:
      Integrand &f;
      EmbeddedRule &rule;
 
   public:
      KahanAdd estimate;
      KahanAdd error;
   };
 
   void MyErrCubePartitioner::action (const Hypercube &h)
   {
      EstErr ee (0.0, 0.0);

      rule.evalError (f, h, ee);       // Evaluate the EmbeddedRule
 
      estimate += ee.getEstimate ();   // ...and update totals
      error    += ee.getError ();
   }
} // namespace

Integrator::Status L::CompRuleIntegratorErr::integrate (
   Integrand &f, const Hypercube &h, Index maxEval,
   real reqAbsError, real reqRelError, EstErr &ee)
{
   checkDimension(h, f);

   if (!maxEval)
   {
      #ifdef HINTLIB_NO_EXCEPTIONS
         ee.set (0.0, 0.0);
         return ERROR;
      #else
         throw MaxEvaluationsRequired();
      #endif
   }

   std::auto_ptr<EmbeddedRule> rule (getRule(h.getDimension()));

   // Calculate the number of sections
 
   Index numSections = maxEval / rule->getNumPoints ();
 
   // We need at least one section to apply the basic rule once
 
   if (numSections == 0)
   {
      #ifdef HINTLIB_NO_EXCEPTIONS
         ee.set (0.0, 0.0);
         return ERROR;
      #else
         throw NoEvaluationsPossible (maxEval);
      #endif
   }

   // Use CubePartitioner to perform action an every sub-cube
 
   MyErrCubePartitioner cp (f, *rule);
 
   cp (h, numSections, cp. EQUIDISTANT);
 
   // Copy to final result

   ee.set (cp.estimate, cp.error);
 
   // Return status depending on estimated error
 
   Status status = checkRequestedError (ee, reqAbsError, reqRelError);

   return (status == ERROR) ? MAX_EVAL_REACHED : status;
}

