/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002  Rudolf Sch�rer <rudolf.schuerer@sbg.ac.at>
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <HIntLib/rule3tyler.h>

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/mymath.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;
using L::real;


/**
 *  Constructor
 *
 *  The constructor is used primarily to initialize all the dimension dependent
 *  constatns and to allocate (dimension dependent) memory
 */

L::Rule3Tyler::Rule3Tyler (unsigned dim)
: OrbitRule (dim), b0 ((3.0 - dim) / 3.0), b1 (1.0 / 6.0)
{
   checkDimensionNotZero (dim);
}


/**
 *  getSumAbsWeight()
 */

real L::Rule3Tyler::getSumAbsWeight() const
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

real L::Rule3Tyler::eval (Function &f, const Hypercube &h)
{
   const real* center = h.getCenter();

   Scaler scaler (h.getWidth(), 1.0);
   setCenter (center);

   return h.getVolume() *
   (   b0 * eval0_0    (f)
     + b1 * evalR0_0fs (f, center, scaler)
   );
}

L::CubatureRuleFactory* L::Rule3Tyler::getFactory()
{
   return new DefaultCubatureRuleFactory<L::Rule3Tyler> ();
}
