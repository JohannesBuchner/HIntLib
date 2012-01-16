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
 *  TestFunction
 *
 *  (C) 2001 by Rudolf Schürer
 */

#ifdef __GNUG__
#pragma implementation
#pragma implementation "function.h"
#endif
 
#include <HIntLib/testfunction.h>

#include <HIntLib/hypercube.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;

using L::TestFunction;
using L::Hypercube;
using L::real;

/**
 *  DomainCheckerFunction
 */

L::DomainCheckerFunction::DomainCheckerFunction (
   TestFunction *pf, const Hypercube *ph)
: TestFunction (pf->getDimension()), f(*pf), h(*ph), allInside(true)
{
   checkDimensionEqual (h.getDimension(), f.getDimension());
}

real L::DomainCheckerFunction::operator() (const real p [])
{
   allInside &= isPointInside (h, p);

   return f(p);
}


