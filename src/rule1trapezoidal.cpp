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
 *  Rule 1 Trapezoidal
 *
 *  Cubature rule of degree 1 with 2^n points.
 *  All points are inside the hypercube.
 *
 *  A product rule based on the one-dimensional 2-point trapezoidal rule
 *
 *  It is presented in
 *     A. H. Stoud. Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 1-2
 */

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/rule1trapezoidal.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/exception.h>
#include <HIntLib/hypercube.h>


namespace L = HIntLib;


/**
 *  The constructor is used primarily to initialize all the dimension dependent
 *  constatns and to allocate (dimension dependent) memory
 */

L::Rule1Trapezoidal::Rule1Trapezoidal (int dim)
   : OrbitRule (dim), a (dim), oneDivTwoPowDim (real(1) / (Index(1) << dim))
{
   checkDimensionNotZero (dim);
   checkDimensionLeq<std::numeric_limits<Index>::digits - 1> (dim);
}


/**
 *  Do the actual function evaluation
 */

L::real L::Rule1Trapezoidal::eval (Integrand &f, const Hypercube &h)
{
   return h.getVolume() * evalR_Rfs (f, h.getCenter(), h.getWidth())
                        * oneDivTwoPowDim;
}

L::CubatureRuleFactory* L::Rule1Trapezoidal::getFactory()
{
   return new DefaultCubatureRuleFactory<L::Rule1Trapezoidal> ();
}

