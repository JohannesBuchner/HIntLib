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
 *  rule3octahedron.h
 *
 *  Cubature rule of degree 3 with 2 * dim points.
 *  All points are inside the hypercube. All weights are positive.
 *
 *  This rule was published in
 *     A. H. Stroud., Remarks on the disposition of points in numerical
 *        integration. Math Tables Aids Comput, V 11 (1957), 257-261
 *
 *  It is also presented in
 *     A. H. Stoud.  Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 3-1
 */

#ifndef RULE3OCTAHEDRON_H
#define RULE3OCTAHEDRON_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/cubaturerule.h>
#include <HIntLib/array.h>


namespace HIntLib
{
   class CubatureRuleFactory;

   class Rule3Octahedron : public CubatureRule
   {
   public:
      Rule3Octahedron (unsigned dim);

      virtual real eval (Function &, const Hypercube &);

      virtual unsigned getDimension()  const  { return dim; }
      virtual Index getNumPoints()     const  { return 2 * dim; }
      virtual unsigned getDegree()     const  { return 3; }
      virtual bool isAllPointsInside() const  { return true; }
      virtual real getSumAbsWeight()   const  { return 1.0; }

      static CubatureRuleFactory* getFactory();

   private:

      const unsigned dim;
      const unsigned dim2;

      // Dimension dependent constants

      Array<real> r; 
      Array<real> p;
   };

}  // namespace HIntLib

#endif

