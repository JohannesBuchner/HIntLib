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

#ifndef HINTLIB_RULEINTEGRATOR_H
#define HINTLIB_RULEINTEGRATOR_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/rulebasedintegrator.h>

namespace HIntLib
{
   class RuleIntegrator : public CubatureRuleBasedIntegrator
   {
   public:
      RuleIntegrator (const CubatureRuleFactory *fac)
         : CubatureRuleBasedIntegrator(fac) {}

      virtual
      Status integrate (
         Integrand &, const Hypercube &, Index maxEval,
         real reqAbsError, real reqRelError, EstErr &);
   };

   class RuleIntegratorErr : public EmbeddedRuleBasedIntegrator
   {
   public:
      RuleIntegratorErr (const EmbeddedRuleFactory *fac)
         : EmbeddedRuleBasedIntegrator(fac) {}

      virtual
      Status integrate (
         Integrand &, const Hypercube &, Index maxEval,
         real reqAbsError, real reqRelError, EstErr &);
   };
}  // namespace HIntLib
 
#endif

