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
 *  ruleMidpoint.cpp
 *
 *  Cubature rule of degree 1 with 1 point (the center of the hypercube).
 *
 *  It is presented in
 *     A. H. Stoud. Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 1-1
 */

#ifdef __GNUG__
#pragma implementation
#endif

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/rule1midpoint.h>

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/integrand.h>

namespace L = HIntLib;

/**
 *  Calculate the volume of the hypercube and multiply with integrand value
 *  in the center.
 */

L::real L::Rule1Midpoint::eval (Integrand &f, const Hypercube &h)
{
   return h.getVolume() * f(h.getCenter());
}

/**
 *  getFactory()
 */

L::CubatureRuleFactory* L::Rule1Midpoint::getFactory()
{
   return new DefaultCubatureRuleFactory<L::Rule1Midpoint> ();
}

