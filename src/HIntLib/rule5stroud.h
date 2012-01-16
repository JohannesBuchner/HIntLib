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
 *  rule5stroud.h
 *
 *  Cubature rule of degree 5 with 3*dim^2 + 3*dim + 1 points.
 *  All points are inside the hypercube.
 *
 *  This rule was published in
 *     A. H. Stroud., Extensions of Symmetric Integration Formulas.
 *        Math Comput, V 22 (1968), 271-274
 *
 *  It is also presented in
 *     A. H. Stoud.  Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 5-3
 */

#ifndef HINTLIB_RULE5STROUD_H
#define HINTLIB_RULE5STROUD_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/cubaturerule.h>
#include <HIntLib/orbitrule.h>


namespace HIntLib
{
   class CubatureRuleFactory;

   class Rule5Stroud : public CubatureRule, private OrbitRule
   {
   public:
      Rule5Stroud (unsigned dim);

      virtual real eval (Integrand &, const Hypercube &);

      virtual unsigned getDimension()      const  { return dim; }
      virtual Index    getNumPoints()      const;
      virtual unsigned getDegree()         const  { return 5; }
      virtual bool     isAllPointsInside() const  { return true; }
      virtual real     getSumAbsWeight()   const;

      static CubatureRuleFactory* getFactory();

   private:
      // Arrays for temporary data

      Point aR;
      Point aMinusR;
      Point aS;
      Point aMinusS;
      Point aT;
      Point aMinusT;

      // Dimension dependent constants

      const real b0;
      const real b2;
      const real b4;
   };

}  // namespace HIntLib

#endif

