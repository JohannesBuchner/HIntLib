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
 *  RuleIntegrator
 *
 *  This integrator applies a basic rule once to estimate value and result
 *  of an integral.
 */

#define HINTLIB_LIBRARY_OBJECT

#include <memory>

#include <HIntLib/ruleintegrator.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/cubaturerule.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;

L::Integrator::Status L::RuleIntegrator::integrate (
   Integrand &f, const Hypercube &h, Index maxEval,
   real reqAbsError, real reqRelError, EstErr &ee)
{
   checkDimension(h, f);
   checkTerminationCriteria (maxEval, reqAbsError, reqRelError, false);

   std::auto_ptr<CubatureRule> rule (getRule(h.getDimension()));

   if (maxEval < rule->getNumPoints ())
   {
      #ifdef HINTLIB_NO_EXCEPTIONS
         ee.set (0.0, 0.0);
         return ERROR;
      #else
         throw NoEvaluationsPossible (maxEval);
      #endif
   }

   ee.setNoErr (rule->eval(f, h));

   return MAX_EVAL_REACHED;
}

L::Integrator::Status L::RuleIntegratorErr::integrate (
   Integrand &f, const Hypercube &h, Index maxEval,
   real reqAbsError, real reqRelError, EstErr &ee)
{
   checkDimension(h, f);
   checkTerminationCriteria (maxEval, reqAbsError, reqRelError, false);

   std::auto_ptr<EmbeddedRule> rule (getRule(h.getDimension()));

   if (maxEval < rule->getNumPoints ())
   {
      #ifdef HINTLIB_NO_EXCEPTIONS
         ee.set (0.0, 0.0);
         return ERROR;
      #else
         throw NoEvaluationsPossible (maxEval);
      #endif
   }

   rule->evalError (f, h, ee);

   Status status = checkRequestedError (ee, reqAbsError, reqRelError);

   return (status == ERROR) ? MAX_EVAL_REACHED : status;
}

