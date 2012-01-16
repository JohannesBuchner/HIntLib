/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schürer <rudolf.schuerer@sbg.ac.at>
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

#ifndef HINTLIB_CUBATURERULE_H
#define HINTLIB_CUBATURERULE_H 1

#ifdef __GNUG__
#pragma interface
// Implementation in rulebasedintegrator.cpp
#endif

#include <HIntLib/defaults.h>
#include <HIntLib/hypercube.h>
#include <HIntLib/esterr.h>

namespace HIntLib
{
   class Integrand;

/**
 *  Cubature Rule
 *
 *  Abstract base class for Integration rules for hyperrectangles
 */

class CubatureRule
{
public:

   CubatureRule() {}
   virtual ~CubatureRule() {}

   // Evaluates the rule on a given hyper-rectangle for a given function.
   // Estimates the integral and an error estimate

   virtual real eval (Integrand &, const Hypercube &) = 0;

   // Here are some simple methods to query (static) information about a
   // certain rule

   virtual unsigned getDimension()      const = 0;
   virtual Index    getNumPoints()      const = 0;
   virtual unsigned getDegree()         const = 0;
   virtual bool     isAllPointsInside() const = 0;
   virtual real     getSumAbsWeight()   const = 0;

private:

   // Do not copy. Do not assign

   CubatureRule (const CubatureRule&);
   CubatureRule& operator= (const CubatureRule&);
};


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

   virtual unsigned evalError (Integrand &, const Hypercube &, EstErr &ee) = 0;

   // Discard error from evalError() to get plain eval() done 

   virtual real eval (Integrand &, const Hypercube &);
};

inline
real EmbeddedRule::eval (Integrand &f, const Hypercube &h)
{
   EstErr ee;
   evalError (f, h, ee);
   return ee.getEstimate ();
}


/**
 *  Cubature Rule Factory
 *
 *  Interface to factory opbjects creating CubatureRules and EmbeddedRules
 *  in arbitrary dimenstion.
 *
 *  Two methods are provided:
 *
 *  create(unsigned)
 *     Returns a pointer to a new CubatureRule (allocated on free store) for
 *     the given dimension.
 *
 *  clone()
 *     Returns a copy of the factory (also allocated on free store).
 *
 *  Factories cannot be copied or assigned.
 */

class CubatureRuleFactory
{
public:
   CubatureRuleFactory() {}
   virtual ~CubatureRuleFactory() {}
   virtual CubatureRule* create (unsigned) = 0;
   virtual CubatureRuleFactory* clone() const = 0;
private:
   CubatureRuleFactory (const CubatureRuleFactory&);
   CubatureRuleFactory& operator= (const CubatureRuleFactory&);
};


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

