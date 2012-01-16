/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2007  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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
 *  rule5mustardlynessblatt.h
 *
 *  Cubature rule of degree 5 with  2^dim + 2*dim + 1  points.
 *  All points are inside the hypercube.
 *
 *  This rule was published in
       D. Mustard, J. N. Lyness, and J. M. Blatt.  Numerical quadrature in n
          dimensions.  Comp. J. 6, 75-87, 1963.
 *
 *  It is also presented in
 *     A. H. Stoud.  Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 5-5
 */

#ifndef HINTLIB_RULE_5_MUSTARD_LYNESS_BLATT_H
#define HINTLIB_RULE_5_MUSTARD_LYNESS_BLATT_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/cubaturerule.h>
#include <HIntLib/orbitrule.h>


namespace HIntLib
{
   class CubatureRuleFactory;

   class Rule5MustardLynessBlatt : public CubatureRule, private OrbitRule
   {
   public:
      Rule5MustardLynessBlatt (int dim);

      virtual real eval (Integrand &, const Hypercube &);

      virtual int   getDimension()      const  { return dim; }
      virtual Index getNumPoints()      const;
      virtual int   getDegree()         const  { return 5; }
      virtual bool  isAllPointsInside() const  { return true; }
      virtual real  getSumAbsWeight()   const;

      static CubatureRuleFactory* getFactory();

   private:
      // Dimension dependent constants

      const real weight0;
      const real weight1;
   };

}  // namespace HIntLib

#endif

