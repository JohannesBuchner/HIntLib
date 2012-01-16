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

#ifdef __GNUG__
#pragma implementation
#endif

#include <HIntLib/rule3cross.h>

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/mymath.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;


/**
 *  The constructor is used primarily to initialize all the dimension dependent
 *  constatns and to allocate (dimension dependent) memory
 */

L::Rule3Cross::Rule3Cross (unsigned dim)
   : OrbitRule (dim),
     r (sqrt (dim / real (3.0))),
     oneOver2Dim (real (1.0) / (2 * dim))
{
   checkDimensionNotZero (dim);
}


/**
 *  Do the actual function evaluation
 */

L::real L::Rule3Cross::eval (Integrand &f, const Hypercube &h)
{
   const real* center = h.getCenter();

   Scaler scaler (h.getWidth(), r);
   setCenter (center);
   return h.getVolume() * evalR0_0fs (f, center, scaler) * oneOver2Dim;
}

L::CubatureRuleFactory* L::Rule3Cross::getFactory()
{
   return new DefaultCubatureRuleFactory<L::Rule3Cross> ();
}

