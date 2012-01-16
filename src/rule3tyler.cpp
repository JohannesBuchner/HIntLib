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
 *  Rule 3 Tyler
 *
 *     (0,...,0)FS        (3-n)/3
 *     (1,0,...,0)FS      1/6
 *
 *  Cubature rule of degree  3  with  2n + 1  points.
 *  All points are inside the cube.
 *
 *  It is also presented in
 *     A. H. Stoud. Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 3-3
 */

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/rule3tyler.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/exception.h>
#include <HIntLib/hypercube.h>


namespace L = HIntLib;


/**
 *  Constructor
 *
 *  The constructor is used primarily to initialize all the dimension dependent
 *  constatns and to allocate (dimension dependent) memory
 */

L::Rule3Tyler::Rule3Tyler (int dim)
   : OrbitRule (dim), b0 (real (3 - dim) / real(3))
{
   checkDimensionNotZero (dim);
}

namespace
{
   const L::real b1 = L::real(1) / L::real(6);
}


/**
 *  getSumAbsWeight()
 */

L::real
L::Rule3Tyler::getSumAbsWeight() const
{
   // Multiply each weight with the number of sampling points it is used for.
   // Don't forget to take the absolute value for weights that meight be
   // negative

   return num0_0() * abs(b0) + numR0_0fs() * b1;
}


/**
 *  eval()
 *
 *  Do the actual function evaluation
 */

L::real
L::Rule3Tyler::eval (Integrand &f, const Hypercube &h)
{
   const real* center = h.getCenter();

   setCenter (center);

   return h.getVolume() *
   (   b0 * eval0_0    (f)
     + b1 * evalR0_0fs (f, center, h.getWidth())
   );
}

L::real
L::Rule3TylerDim3::eval (Integrand &f, const Hypercube &h)
{
   const real* center = h.getCenter();
   const real* width  = h.getWidth();

   real result = 0;
   real point[3];

   point[0] = center[0];
   point[1] = center[1];
   point[2] = center[2] + width[2];

   result += f(point);

   point[2] = center[2] - width[2];

   result += f(point);

   point[2] = center[2];
   point[1] = center[1] + width[1];

   result += f(point);

   point[1] = center[1] - width[1];

   result += f(point);

   point[1] = center[1];
   point[0] = center[0] + width[0];

   result += f(point);

   point[0] = center[0] - width[0];

   result += f(point);

   return result * h.getVolume() / 6;
}


/**
 *  Declare a CubatureRuleFactory which returns the proper class depending on
 *  the dimension.
 */

namespace
{
   class Rule3TylerCubatureRuleFactory : public L::CubatureRuleFactory
   {
   public:
      virtual L::CubatureRule* create (int dim);
      virtual Rule3TylerCubatureRuleFactory* clone () const;
   };

   L::CubatureRule*
   Rule3TylerCubatureRuleFactory::create (int dim)
   {
      return (dim == 3)
         ? static_cast<L::CubatureRule*> (new L::Rule3TylerDim3())
         : static_cast<L::CubatureRule*> (new L::Rule3Tyler(dim));
   }

   Rule3TylerCubatureRuleFactory*
   Rule3TylerCubatureRuleFactory::clone () const
   {
      return new Rule3TylerCubatureRuleFactory();
   }
}


/**
 *  getFactory()
 */

L::CubatureRuleFactory* L::Rule3Tyler::getFactory()
{
   return new Rule3TylerCubatureRuleFactory();
}

