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
 *  rule5stroud.cpp
 *
 *  Cubature rule of degree 5 with 3*dim^2 + 3*dim + 1 points.
 *  All points are inside the hypercube.
 *
 *  This rule was published in
 *     A. H. Stroud., Extensions of Symmetric Integration Formulas.
 *        Math Comput, V 22 (1968), 271-274
 *
 *  It is also presented in
 *     A. H. Stoud.  Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 5-3
 */                                                                             

#ifdef __GNUG__
#pragma implementation
#endif

#include <HIntLib/rule5stroud.h>

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/mymath.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;
using L::real;
using L::Index;

namespace
{
   const real b1 = 25.0 / 168.0;
   const real b3 =  5.0 /  48.0;

   #if HINTLIB_STATIC_WORKS == 1
   const real r = sqrt ( 7.0                / 15.0);
   const real s = sqrt ((7.0 + sqrt (24.0)) / 15.0);
   const real t = sqrt ((7.0 - sqrt (24.0)) / 15.0);
   #else
   real r,s,t;
   #endif
}

/**
 *  The constructor is used primarily to initialize all the dimension dependent
 *  constatns and to allocate (dimension dependent) memory
 */

L::Rule5Stroud::Rule5Stroud (unsigned dim)
: OrbitRule (dim),

  aR      (dim),
  aMinusR (dim),
  aS      (dim),
  aMinusS (dim),
  aT      (dim),
  aMinusT (dim),

  b0 ((5 * sqr(dim) - 15 * dim + 14) / 14.0),
  b2 ((50.0 - 25.0 * dim) / 168.0),
  b4 ((10.0 -  5.0 * dim) /  48.0)
{
   #if HINTLIB_STATIC_WORKS == 0
      r = sqrt ( 7.0                / 15.0);
      s = sqrt ((7.0 + sqrt (24.0)) / 15.0);
      t = sqrt ((7.0 - sqrt (24.0)) / 15.0);
   #endif

   checkDimensionNotZero (dim);
   checkDimensionGeq<2> (dim);
}



Index L::Rule5Stroud::getNumPoints () const
{
   return num0_0 () + 2 * numRR0_0s () + 3 * numR0_0fs () + 2 * numRS0_0s ();
}

real L::Rule5Stroud::getSumAbsWeight() const
{
   // Multiply each weight with the number of sampling points it is used for.
   // Don't forget to take the absolute value for weights that meight be
   // negative
 
   return     num0_0    () * abs (b0)
      + 2.0 * numRR0_0s () * b1
      +       numR0_0fs () * (abs (b2) + 2.0 * abs (b4))
      + 2.0 * numRS0_0s () * b3;
}
 

/**
 *  Do the actual function evaluation
 */

L::real L::Rule5Stroud::eval (Function &f, const Hypercube &h)
{
   // Calculate center and width of the rectange. Adjust volume accordingly

   const real* width = h.getWidth();
   for (unsigned int i = 0; i < dim; i++)
   {
      real w = width [i];

      aR      [i] =   w * r;
      aMinusR [i] = - w * r;
      aS      [i] =   w * s;
      aMinusS [i] = - w * s;
      aT      [i] =   w * t;
      aMinusT [i] = - w * t;
   }

   // Evaluate integrand

   const real* c = h.getCenter();
   setCenter (c);

   return h.getVolume() *
   (
        b0 * eval0_0 (f)
      + b1 * (  evalRR0_0s (f, c, aR)
              + evalRR0_0s (f, c, aMinusR)
             )
      + b2 * evalR0_0fs (f, c, aR)
      + b3 * (  evalRS0_0s (f, c, aS, aMinusT)
              + evalRS0_0s (f, c, aMinusS, aT)
             )
      + b4 * (  evalR0_0fs (f, c, aS)
              + evalR0_0fs (f, c, aT)
             )
   );
}

/**
 *  getFactory()
 */

L::CubatureRuleFactory* L::Rule5Stroud::getFactory()
{
   return new DefaultCubatureRuleFactory<L::Rule5Stroud> ();
}



