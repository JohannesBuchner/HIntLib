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
 *  Rule 2 Thacher
 *
 *  Cubature rule of degree 2 with  2*dim+1  points.
 *  All points are inside the hypercube.
 *
 *  It is also presented in
 *     A. H. Stoud. Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 2-2
 */

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/rule2thacher.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/hlmath.h>
#include <HIntLib/exception.h>
#include <HIntLib/integrand.h>


namespace L = HIntLib;

using L::real;

namespace
{
#if HINTLIB_STATIC_WORKS == 1
   const real r  = HINTLIB_MN sqrt (real(3)) / real(6);
   const real r2 = HINTLIB_MN sqrt (real(3)) / real(3);
#else
   real r, r2;
#endif
}

/**
 *  Constructor
 */

L::Rule2Thacher::Rule2Thacher (unsigned _dim)
: dim(_dim), p(_dim)
{
   checkDimensionNotZero (dim);

#if HINTLIB_STATIC_WORKS != 1
   r  = HINTLIB_MN sqrt (real(3)) / real(6);
   r2 = HINTLIB_MN sqrt (real(3)) / real(3);
#endif
}

/**
 *  Do the actual function evaluation
 */

real L::Rule2Thacher::eval (Integrand &f, const Hypercube &h)
{
   // Sample all points

   const real* center = h.getCenter();
   const real* width  = h.getWidth();
   real sum = 0.0;

   for (unsigned i = 0; i < dim; ++i)  p [i] = center[i] + r2 * width[i];

   sum += f(p);

   for (unsigned j = 0; j < dim; ++j)
   {
      for (unsigned i = 0; i < j; ++i)      p [i] = center[i] + r * width[i];
      p[j] = center[j] + width [j];
      for (unsigned i = j+1; i < dim; ++i)  p [i] = center[i] + r * width[i];

      sum -= r * f(p);

      p[j] = center[j] - width [j];

      sum += r * f(p);
   }

   return h.getVolume() * sum;
}

real L::Rule2Thacher::getSumAbsWeight() const
{
   return real(1) + real (dim) * r2;
}


L::CubatureRuleFactory* L::Rule2Thacher::getFactory()
{
   return new DefaultCubatureRuleFactory<L::Rule2Thacher> ();
}

