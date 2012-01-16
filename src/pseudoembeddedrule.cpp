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


#ifdef __GNUG__
#pragma implementation
#endif

#include <algorithm>

#include <HIntLib/pseudoembeddedrule.h>

#include <HIntLib/mymath.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;
using L::Index;
using L::real;


/**
 *  getNumPoints()
 */

inline
Index L::PseudoEmbeddedRule::getNumPoints() const
{
   return r1->getNumPoints() + r2->getNumPoints() + fd.getNumPoints();
}

Index L::PseudoDoubleEmbeddedRule::getNumPoints() const
{
   return PseudoEmbeddedRule::getNumPoints() + r3->getNumPoints();
}


/**
 *  isAllPointsInside()
 */

inline
bool L::PseudoEmbeddedRule::isAllPointsInside() const
{
   return r1->isAllPointsInside() && r2->isAllPointsInside();
}

bool L::PseudoDoubleEmbeddedRule::isAllPointsInside() const
{
   return PseudoEmbeddedRule::isAllPointsInside() && r3->isAllPointsInside();
}

/**
 *  getSumAbsWeight()
 */

real L::PseudoEmbeddedRule::getSumAbsWeight() const
{
   return r1->getSumAbsWeight();
}


/**
 *  evalError()
 */

unsigned L::PseudoEmbeddedRule::evalError (Integrand &f,
   const Hypercube &h, EstErr &ee)
{
   real result1 = r1->eval(f, h);

   real result2 = r2->eval(f, h);

   ee.set (result1, abs (result1 - result2));

   return fd (f, h);
}


unsigned L::PseudoDoubleEmbeddedRule::evalError (Integrand &f,
   const Hypercube &h, EstErr &ee)
{
   unsigned split = PseudoEmbeddedRule::evalError (f, h, ee);

   real result3 = r3->eval (f, h);

   // Set estimate from r1 and estimated error in result
 
   ee.set (ee.getEstimate(),
           std::max (ee.getError(), abs (ee.getEstimate() - result3)));

   return split;
} 


/**
 *  PseudoEmbeddedRuleFactory
 */

L::PseudoEmbeddedRuleFactory* L::PseudoEmbeddedRuleFactory::clone() const
{
   return new PseudoEmbeddedRuleFactory (factory1->clone(), factory2->clone());
}
 
L::PseudoEmbeddedRule* L::PseudoEmbeddedRuleFactory::create (unsigned dim)
{
   return new PseudoEmbeddedRule (dim, factory1.get(), factory2.get());
}


/**
 *  PseudoDoubleEmbeddedRuleFactory
 */

L::PseudoDoubleEmbeddedRuleFactory*
L::PseudoDoubleEmbeddedRuleFactory::clone() const
{
   return new PseudoDoubleEmbeddedRuleFactory (
         factory1->clone(), factory2->clone(), factory3->clone());
}
 
L::PseudoDoubleEmbeddedRule*
L::PseudoDoubleEmbeddedRuleFactory::create (unsigned dim)
{
   return new PseudoDoubleEmbeddedRule (
         dim, factory1.get(), factory2.get(), factory3.get());
}

