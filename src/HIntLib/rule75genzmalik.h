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
 *  Rule 7-5 Genz
 *
 *  An embedded cubatare rule of degree 7 (embedded rule degree 5) due to
 *  A.Genz and Malik
 */

#ifndef HINTLIB_RULE_7_5_GENZ_MALIK_H
#define HINTLIB_RULE_7_5_GENZ_MALIK_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/cubaturerule.h>
#include <HIntLib/orbitrule.h>


namespace HIntLib
{
   class EmbeddedRuleFactory;

   class Rule75GenzMalik : public EmbeddedRule, private OrbitRule
   {
   public:

      Rule75GenzMalik (int dim);

      virtual int evalError (Integrand &f, const Hypercube &, EstErr &);

      virtual int   getDimension()      const  { return dim; }
      virtual Index getNumPoints()      const;
      virtual int   getDegree()         const  { return 7; }
      virtual bool  isAllPointsInside() const  { return true; }
      virtual real  getSumAbsWeight()   const;

      static EmbeddedRuleFactory* getFactory();

   private:

      // Arrays for temporary data

      Point widthLamda;

      // Dimension depending constants

      const real weight1,weight3, weight5;
      const real weightE1, weightE3;
   };

}  // namespace HIntLib

#endif

