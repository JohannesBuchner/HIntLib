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
 *  rule7phillips.h
 *
 *  Cubature rule of degree 7
 *  All points are inside the hypercube.
 *
 *  This rule was published in
 *     G.M.Phillips:
 *     Numerical Integration over an N-Dimensional Rectangular Region
 *     The Computer Journal 10 (1967), 297-299.
 *
 *  It is also presented in
 *     A. H. Stoud. Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 7-1
 */

#ifndef HINTLIB_RULE7PHILLIPS_H
#define HINTLIB_RULE7PHILLIPS_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/cubaturerule.h>
#include <HIntLib/orbitrule.h>


namespace HIntLib
{
   class CubatureRuleFactory;

   class Rule7Phillips : public CubatureRule, private OrbitRule
   {
   public:
      Rule7Phillips (unsigned dim);

      virtual real eval (Integrand &, const Hypercube &);

      virtual unsigned getDimension()       const  { return dim; }
      virtual Index    getNumPoints ()      const;
      virtual unsigned getDegree ()         const  { return 7; }
      virtual bool     isAllPointsInside () const  { return true; }
      virtual real     getSumAbsWeight ()   const;

      static CubatureRuleFactory* getFactory();

   private:

      // Arrays for temporary data

      Point aCD;

      // Dimension dependent constants

      const real eN;
      const real temp;
      const real rB2;
      const real weightB1, weightB2, weightC, weightA;
   };

}  // namespace HIntLib

#endif

