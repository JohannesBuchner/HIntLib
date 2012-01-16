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
 *  RuleBasedIntegrator
 */

#ifndef RULEBASEDINTEGRATOR_H
#define RULEBASEDINTEGRATOR_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/integrator.h>

namespace HIntLib
{
   class CubatureRule;
   class EmbeddedRule;
   class CubatureRuleFactory;
   class EmbeddedRuleFactory;

   class CubatureRuleBasedIntegrator: public Integrator
   {
   protected:
      CubatureRuleBasedIntegrator (const CubatureRuleFactory *);

      CubatureRule* getRule (unsigned) const;

   public:
      virtual ~CubatureRuleBasedIntegrator();

   private:
      CubatureRuleFactory *factory;
   };

   class EmbeddedRuleBasedIntegrator: public Integrator
   {
   protected:
      EmbeddedRuleBasedIntegrator (const EmbeddedRuleFactory *);

      EmbeddedRule* getRule (unsigned) const;

   public:
      virtual ~EmbeddedRuleBasedIntegrator();

   private:
      EmbeddedRuleFactory *factory;
   };
}  // namespace HIntLib

#endif

