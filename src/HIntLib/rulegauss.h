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
 *  Rule Gauss
 */

#ifndef HINTLIB_RULE_GAUSS_H
#define HINTLIB_RULE_GAUSS_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/productrule.h>

namespace HIntLib
{
   class CubatureRuleFactory;

   class RuleGauss : public ProductRule
   {
   public:
      RuleGauss (int dim, int numPoints);

      // The parameters are all trivially calculated

      virtual int getDegree()   const { return numPoints * 2 - 1; }
      virtual real getSumAbsWeight() const { return 1.; }

      // getFactory() needs to know the number of abscissas

      static CubatureRuleFactory* getFactory (int numPoints);

   private:
      Point offsets;
      Point weights;
   };

}  // namespace HIntLib

#endif

