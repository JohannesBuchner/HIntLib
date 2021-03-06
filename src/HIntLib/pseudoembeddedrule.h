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


#ifndef HINTLIB_PSEUDO_EMBED_RULE_H
#define HINTLIB_PSEUDO_EMBED_RULE_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <memory>

#include <HIntLib/cubaturerule.h>
#include <HIntLib/fourthdiff.h>


namespace HIntLib
{

/**
 *  PseudoEmbeddedRule
 */

class PseudoEmbeddedRule : public EmbeddedRule
{
public:
   PseudoEmbeddedRule (
      int d,
      CubatureRuleFactory *fac1,
      CubatureRuleFactory *fac2)
   : r1 (fac1->create(d)), r2 (fac2->create(d)), fd(d) {}

   virtual int evalError (Integrand &, const Hypercube &, EstErr &);

   virtual int   getDimension()      const  { return fd.getDimension(); }
   virtual Index getNumPoints()      const;
   virtual int   getDegree()         const  { return r1->getDegree(); }
   virtual bool  isAllPointsInside() const;
   virtual real  getSumAbsWeight()   const;

protected:
   const std::auto_ptr<CubatureRule> r1;
   const std::auto_ptr<CubatureRule> r2;
   FourthDiff fd;
};


/**
 *  PseudoDoubleEmbeddedRule
 */

class PseudoDoubleEmbeddedRule : public PseudoEmbeddedRule
{
public:
   PseudoDoubleEmbeddedRule (
      int dim,
      CubatureRuleFactory *fac1,
      CubatureRuleFactory *fac2,
      CubatureRuleFactory *fac3)
   : PseudoEmbeddedRule(dim, fac1, fac2), r3 (fac3->create(dim)) {}

   virtual int evalError (Integrand &f, const Hypercube &, EstErr &);

   virtual Index getNumPoints() const;
   virtual bool isAllPointsInside() const;

private:
   const std::auto_ptr<CubatureRule> r3;
};


/**
 *  PseudoEmbeddedRuleFactory
 *
 *  Allows the creation of PseudoEmbeddedRules with arbitrary dimensions
 */

class PseudoEmbeddedRuleFactory : public EmbeddedRuleFactory
{
public:
   PseudoEmbeddedRuleFactory (CubatureRuleFactory *fac1,
                              CubatureRuleFactory *fac2)
   : factory1 (fac1), factory2 (fac2) {}

   virtual PseudoEmbeddedRuleFactory* clone() const;
   virtual PseudoEmbeddedRule* create (int dim);

private:
   const std::auto_ptr<CubatureRuleFactory> factory1;
   const std::auto_ptr<CubatureRuleFactory> factory2;
};


/**
 *  PseudoDoubleEmbeddedRuleFactory
 *
 *  Allows the creation of PseudoDoubleEmbeddedRules with arbitrary dimensions
 */

class PseudoDoubleEmbeddedRuleFactory : public EmbeddedRuleFactory
{
public:
   PseudoDoubleEmbeddedRuleFactory (
      CubatureRuleFactory *fac1,
      CubatureRuleFactory *fac2,
      CubatureRuleFactory *fac3)
   : factory1(fac1), factory2(fac2), factory3(fac3) {}

   virtual PseudoDoubleEmbeddedRuleFactory* clone() const;
   virtual PseudoDoubleEmbeddedRule* create (int dim);

private:
   const std::auto_ptr<CubatureRuleFactory> factory1;
   const std::auto_ptr<CubatureRuleFactory> factory2;
   const std::auto_ptr<CubatureRuleFactory> factory3;
};

}  // namespace HIntLib

#endif

