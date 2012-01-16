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
 *  Product Rule
 *
 *  Common base class for all CubatureRules which are tensor product rules of
 *  a one-dimensional rule with weights and abscissas stored in an array.
 */

#ifndef HINTLIB_PRODUCT_RULE_H
#define HINTLIB_PRODUCT_RULE_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/cubaturerule.h>
#include <HIntLib/array.h>
#include <HIntLib/counter.h>

namespace HIntLib
{
   class ProductRule : public CubatureRule
   {
   public:
      virtual real eval (Integrand &, const Hypercube &);

      // The parameters are all trivially calculated

      virtual int   getDimension() const { return dim; }
      virtual Index getNumPoints() const;
      virtual bool  isAllPointsInside() const;

   protected:
      ProductRule (int _dim, int _numPoints);

      void init (const real* _offsets, const real* _weights)
      {
         offsets = _offsets;
         weights = _weights;
      }

      const int dim;
      const int numPoints;

   private:
      Point p;
      GrayCounter counter;
      const real* offsets;
      const real* weights;
   };

}  // namespace HIntLib

#endif

