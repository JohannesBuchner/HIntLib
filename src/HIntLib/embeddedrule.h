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

#ifndef EMBEDRULE_H
#define EMBEDRULE_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/cubaturerule.h>
#include <HIntLib/esterr.h>


namespace HIntLib
{

class Function;

/**
 *  Embedded Rule
 *
 *  Abstract base class for Integration rules for hyperrectangles
 */

class EmbeddedRule : public CubatureRule
{
public:
   // Evaluates the rule on a given hyper-rectangle for a given function.
   // Estimates the integral and an error estimate

   virtual unsigned evalError (Function &, const Hypercube &, EstErr &ee) = 0;

   // Discard error from evalError() to get plain eval() done 

   virtual real eval (Function &, const Hypercube &);
};

inline
real EmbeddedRule::eval (Function &f, const Hypercube &h)
{
   EstErr ee;
   evalError (f, h, ee);
   return ee.getEstimate ();
}


/**
 *  Embedded Rule Factory
 *
 *  Specialization of CubatureRuleFactory  producing EmbeddedRules
 */

class EmbeddedRuleFactory : public CubatureRuleFactory
{
public:
   EmbeddedRuleFactory() {}
   virtual ~EmbeddedRuleFactory() {}
   virtual EmbeddedRule* create (unsigned) = 0;
   virtual EmbeddedRuleFactory* clone() const = 0;
private:
   EmbeddedRuleFactory (const EmbeddedRuleFactory&);
   EmbeddedRuleFactory& operator= (const EmbeddedRuleFactory&);
};

}  // namespace HIntLib

#endif

