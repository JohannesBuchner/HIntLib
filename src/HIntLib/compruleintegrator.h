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

#ifndef COMPRULEINTEGRATOR_H
#define COMPRULEINTEGRATOR_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/rulebasedintegrator.h>

namespace HIntLib
{
   class CompRuleIntegrator : public CubatureRuleBasedIntegrator
   {
   public:
      CompRuleIntegrator (const CubatureRuleFactory *fac)
         : CubatureRuleBasedIntegrator(fac) {}

      virtual
      Status integrate (
         Function &, const Hypercube &, Index maxEval, real, real, EstErr &);
   };

   class CompRuleIntegratorErr : public EmbeddedRuleBasedIntegrator
   {
   public:
      CompRuleIntegratorErr (const EmbeddedRuleFactory *fac)
         : EmbeddedRuleBasedIntegrator(fac) {}

      Status integrate (
         Function &, const Hypercube &, Index maxEval,
         real reqAbsError, real reqRelError, EstErr &ee);
   };
}  // namespace HIntLib
 
#endif
