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

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/rulebasedintegrator.h>
#include <HIntLib/cubaturerule.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#pragma implementation "cubaturerule.h"
#endif


namespace L = HIntLib;

/**
 *  CubatureRuleBasedIntegrator
 */

L::CubatureRuleBasedIntegrator::CubatureRuleBasedIntegrator
   (const CubatureRuleFactory *fac)
: factory(fac->clone()) {}

L::CubatureRuleBasedIntegrator::~CubatureRuleBasedIntegrator()
{
   delete factory;
}

L::CubatureRule* L::CubatureRuleBasedIntegrator::getRule (unsigned dim) const
{
   return factory->create(dim);
}


/**
 *  EmbeddedRuleBasedIntegrator
 */

L::EmbeddedRuleBasedIntegrator::EmbeddedRuleBasedIntegrator
   (const EmbeddedRuleFactory *fac)
: factory(fac->clone()) {}

L::EmbeddedRuleBasedIntegrator::~EmbeddedRuleBasedIntegrator()
{
   delete factory;
}

L::EmbeddedRule* L::EmbeddedRuleBasedIntegrator::getRule (unsigned dim) const
{
   return factory->create(dim);
}


