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
 *  rule3octahedron.cpp
 *
 *  Cubature rule of degree 3 with 2 * dim points.
 *  All points are inside the hypercube. All weights are positive.
 *
 *  This rule was published in
 *     A. H. Stroud., Remarks on the disposition of points in numerical
 *        integration. Math Tables Aids Comput, V 11 (1957), 257-261
 *
 *  It is also presented in
 *     A. H. Stoud.  Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 3-1                                                          
 */

#ifdef __GNUG__
#pragma implementation
#endif                                                                          

#include <HIntLib/rule3octahedron.h>

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/hlmath.h>
#include <HIntLib/exception.h>
#include <HIntLib/integrand.h>


namespace L = HIntLib;

using L::real;

namespace
{
#if HINTLIB_STATIC_WORKS == 1
   const real sqrt2by3    = HINTLIB_MN sqrt (real(2) / real(3));
   const real oneDivSqrt3 = real(1) / HINTLIB_MN sqrt (real(3));
#else
   real sqrt2by3, oneDivSqrt3;
#endif
}


/**
 *  The constructor is used primarily to initialize all the dimension dependent
 *  constatns and to allocate (dimension dependent) memory
 */

L::Rule3Octahedron::Rule3Octahedron (unsigned d)
   : dim(d), dim2(d * 2), r(d * d * 2), p(d)
{
   checkDimensionNotZero (dim);

#if HINTLIB_STATIC_WORKS == 0
   sqrt2by3    = HINTLIB_MN sqrt (real(2) / real(3));
   oneDivSqrt3 = real(1) / HINTLIB_MN sqrt (real(3));
#endif

   // Initialze r

   for (unsigned i = 0; i < dim2; ++i)
   {
      for (unsigned k = 0; k < dim / 2; ++k)
      {
         real t = (2*k + 1) * (i+1) * M_PI / dim;

         r [i*dim + 2*k] =     sqrt2by3 * HINTLIB_MN cos (t);

         r [i*dim + 2*k + 1] = sqrt2by3 * HINTLIB_MN sin (t);
      }

      if (dim % 2)
      {
         r [i*dim + (dim-1)] = (i % 2) ? oneDivSqrt3 : - oneDivSqrt3;
      }
   }
}


/**
 *  Do the actual function evaluation
 */

real L::Rule3Octahedron::eval (Integrand &f, const Hypercube &h)
{
   // Sample all points

   real sum = 0.0;
   const real* center = h.getCenter();
   const real* width  = h.getWidth();

   for (unsigned i = 0; i < dim2; i++)
   {
      for (unsigned k = 0; k < dim; k++)
      {
         p [k] = center [k] + r [i * dim + k] * width [k];
      }

      sum += f(p);
   }

   return h.getVolume() * sum / real (dim2);
}

L::CubatureRuleFactory* L::Rule3Octahedron::getFactory()
{
   return new DefaultCubatureRuleFactory<L::Rule3Octahedron> ();
}


