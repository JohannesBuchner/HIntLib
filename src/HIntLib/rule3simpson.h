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
 *  Rule 3 Simpson
 *
 *  Cubature rule of degree 3 with 3^n points.
 *  All points (except the center) are on the border of the hypercube.
 *
 *  It is also presented in
 *     A. H. Stoud. Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 3-6
 */

#ifndef HINTLIB_RULE3SIMPSON_H
#define HINTLIB_RULE3SIMPSON_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/cubaturerule.h>
#include <HIntLib/orbitrule.h>


namespace HIntLib
{
   class CubatureRuleFactory;

   class Rule3Simpson : public CubatureRule, private OrbitRule
   {
   public:
      Rule3Simpson (unsigned dim);

      virtual real eval (Integrand &, const Hypercube &);

      virtual unsigned getDimension()   const { return dim; }
      virtual Index getNumPoints ()     const { return num3powS(); }
      virtual unsigned getDegree ()     const { return 3; }
      virtual bool isAllPointsInside () const { return true; }
      virtual real getSumAbsWeight ()   const { return 1.0; }

      static CubatureRuleFactory* getFactory();
 
   private:
      const real oneDivTwoPowDim;
   };

}  // namespace HIntLib

#endif

