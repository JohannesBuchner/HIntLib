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

/**
 *  Rule 2 Thacher
 *
 *  Cubature rule of degree 2 with  2*dim + 1  points.
 *  All points are inside the hypercube.
 *
 *  It is also presented in
 *     A. H. Stoud. Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 2-2
 */

#ifndef HINTLIB_RULE_2_THACHER_H
#define HINTLIB_RULE_2_THACHER_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/cubaturerule.h>
#include <HIntLib/array.h>


namespace HIntLib
{
   class CubatureRuleFactory;

   class Rule2Thacher : public CubatureRule
   {
   public:
      Rule2Thacher (int dim);

      virtual real eval (Integrand &, const Hypercube &);

      virtual int   getDimension()      const  { return dim; }
      virtual Index getNumPoints()      const  { return 2 * dim + 1; }
      virtual int   getDegree()         const  { return 2; }
      virtual bool  isAllPointsInside() const  { return true; }
      virtual real  getSumAbsWeight()   const;

      static CubatureRuleFactory* getFactory();

   private:
      const int dim;
      Point p;
   };

}  // namespace HIntLib

#endif

