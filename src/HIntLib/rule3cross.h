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
 *  Rule 3 Cross
 *
 *     (r,0,...,0)FS      1 / 2*dim
 *
 *     with  r = sqrt(n/3)
 *
 *  Cubature rule of degree 3 with 2n points.
 *  There are points outside the cube for n > 3.
 *
 *  It is also presented in
 *     A. H. Stoud. Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 3-2
 */

#ifndef RULE_3_CROSS_H
#define RULE_3_CROSS_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/cubaturerule.h>
#include <HIntLib/orbitrule.h>


namespace HIntLib
{
   class CubatureRuleFactory;

   class Rule3Cross : public CubatureRule, private OrbitRule
   {
   public:
      Rule3Cross (unsigned dim);

      virtual real eval (Function &, const Hypercube &);

      virtual unsigned getDimension()   const { return dim; }
      virtual Index getNumPoints ()     const { return numR0_0fs(); }
      virtual unsigned getDegree ()     const { return 3; }
      virtual bool isAllPointsInside () const { return dim <= 3; }
      virtual real getSumAbsWeight ()   const { return 1.0; }

      static CubatureRuleFactory* getFactory();
 
   private:
      const real r;
      const real oneOver2Dim;
   };

}  // namespace HIntLib

#endif

