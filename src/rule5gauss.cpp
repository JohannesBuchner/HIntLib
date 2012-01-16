/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schürer <rudolf.schuerer@sbg.ac.at>
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
 *  Rule 5 Gauss
 *
 *  Cubature rule of degree 5 with 3^n points.
 *  All points are inside the hypercube.
 *
 *  It is also presented in
 *     A. H. Stoud. Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 5-9
 */

#ifdef __GNUG__
#pragma implementation
#endif

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/rule5gauss.h>

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/hlmath.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;

using L::real;

namespace
{
#if HINTLIB_STATIC_WORKS == 1
   const real r = HINTLIB_MN sqrt (real(3) / real(5));
#else
   real r;
#endif

   const real w0 = real(8) / real(9);
   const real w1 = real(5) / real(9);
}


/**
 *  The constructor is used primarily to initialize all the dimension dependent
 *  constatns and to allocate (dimension dependent) memory
 */

L::Rule5Gauss::Rule5Gauss (unsigned dim)
   : OrbitRule (dim), a (dim), oneDivTwoPowDim (real(1) / (Index(1) << dim))
{
   checkDimensionNotZero (dim);

   unsigned maxDim = digitsRepresentable(Index(3));

   if (dim > maxDim)  throw DimensionTooHigh (dim, maxDim);

   #if HINTLIB_STATIC_WORKS == 0
      r  = HINTLIB_MN sqrt (real(3) / real(5));
   #endif
}


/**
 *  Do the actual function evaluation
 */

real L::Rule5Gauss::eval (Integrand &f, const Hypercube &h)
{
   const real* width = h.getWidth();

   for (unsigned d = 0; d < dim; ++d)  a[d] = width[d] * r;

   return h.getVolume() * eval3powS (f, h.getCenter(), a, w0, w1)
        * oneDivTwoPowDim;
}

L::CubatureRuleFactory* L::Rule5Gauss::getFactory()
{
   return new DefaultCubatureRuleFactory<L::Rule5Gauss> ();
}

