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
 *  TestIntegrand
 *
 *  (C) 2001 by Rudolf Schuerer
 */

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/testintegrand.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#pragma implementation "integrand.h"
#endif

#include <HIntLib/hypercube.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;

using L::TestIntegrand;
using L::Hypercube;
using L::real;

/**
 *  Integrand
 */

real L::Integrand::derivative (const real [], int)
{
   throw DerivativeNotSupported ();
}

/**
 *  DomainCheckerIntegrand
 */

L::DomainCheckerIntegrand::DomainCheckerIntegrand (
   TestIntegrand *pf, const Hypercube *ph)
: TestIntegrand (pf->getDimension()), f(*pf), h(*ph), allInside(true)
{
   checkDimensionEqual (h.getDimension(), f.getDimension());
}

real L::DomainCheckerIntegrand::operator() (const real p [])
{
   allInside = allInside && isPointInside (h, p);

   return f(p);
}

real L::DomainCheckerIntegrand::derivative (const real p [], int a)
{
   allInside = allInside && isPointInside (h, p);

   return f.derivative (p, a);
}

