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
 *  rule5stroud2.h
 *
 *  Cubature rule of degree 5 with  2^dim + 2*dim  points.
 *  All points are inside the hypercube for  2 <= dim <= 5.
 *
 *  This rule was published in
 *     A. H. Stroud., Some fifth degree integration formulas for symmetric
 *        regions. Math. Comp. 20, 90-97, 1966.
 *
 *  It is also presented in
 *     A. H. Stoud.  Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 5-4
 */

#ifndef HINTLIB_RULE_5_STROUD_2_H
#define HINTLIB_RULE_5_STROUD_2_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/cubaturerule.h>
#include <HIntLib/orbitrule.h>


namespace HIntLib
{
   class CubatureRuleFactory;

   class Rule5Stroud2 : public CubatureRule, private OrbitRule
   {
   public:
      Rule5Stroud2 (int dim);

      virtual real eval (Integrand &, const Hypercube &);

      virtual int   getDimension()      const  { return dim; }
      virtual Index getNumPoints()      const;
      virtual int   getDegree()         const  { return 5; }
      virtual bool  isAllPointsInside() const;
      virtual real  getSumAbsWeight()   const  { return 1.; }

      static CubatureRuleFactory* getFactory();

   private:
      // Dimension dependent constants

      const real alpha;
      const real beta;
      const real weightAlpha;
      const real weightBeta;

      // Array for temporary data

      Point scaledBeta;
   };

}  // namespace HIntLib

#endif

