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
 *  Rule 2 Ionescu
 *
 *  Cubature rule of degree 2 with  6  points in dimension 2.
 *  All points are inside the hypercube.
 *
 *  It is also presented in
 *     A. H. Stoud. Approximate Calculation of Multiple Integrals (1971)
 *  as formula C_2: 2-1
 */

#ifdef __GNUG__
#pragma implementation
#endif                                                                          

#include <HIntLib/rule2ionescu.h>

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/exception.h>
#include <HIntLib/integrand.h>


namespace L = HIntLib;

using L::real;


/**
 *  Constructor
 */

L::Rule2Ionescu::Rule2Ionescu (unsigned dim)
{
   checkDimensionNotZero (dim);

   if (dim != 2)  throw InvalidDimension (dim);
}


/**
 *  Do the actual function evaluation
 */

namespace
{
   const real w1_12 = real (1.0) / real (12.0);
   const real w7_12 = real (7.0) / real (12.0);
   const real wm1_6 = real (-1.0) / real (6.0);
}

real L::Rule2Ionescu::eval (Integrand &f, const Hypercube &h)
{
   real sum = 0.0;
   real p [2];

   p[0] = h.getUpperBound (0);
   p[1] = h.getLowerBound (1);
   sum += w1_12 * f(p);
   
   p[0] = h.getLowerBound (0);
   sum += 0.25 * f(p);
   
   p[1] = h.getUpperBound (1);
   sum += w1_12 * f(p);
   
   p[0] = h.getUpperBound (0);
   sum = sum + w7_12 * f(p)

   // evaluate first derivative
   // Scale by  getDiameter() / 2.0  to normalize  h  to [-1, 1]
             + wm1_6 * (  f.derivative (p, 0) * h.getDiameter (0)
                        + f.derivative (p, 1) * h.getDiameter (1));

   return h.getVolume() * sum;
}


L::CubatureRuleFactory* L::Rule2Ionescu::getFactory()
{
   return new DefaultCubatureRuleFactory<L::Rule2Ionescu> ();
}

