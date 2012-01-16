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

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/rule5stroud.h>

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/hlmath.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;
using L::real;
using L::Index;

namespace
{
   const real b1 = real(25) / real(168);
   const real b3 = real (5) / real (48);

#if HINTLIB_STATIC_WORKS == 1
   const real r = HINTLIB_MN sqrt ( real(7) / real(15));
   const real s = HINTLIB_MN sqrt ((real(7) + HINTLIB_MN sqrt (real(24)))
                                  / real(15));
   const real t = HINTLIB_MN sqrt ((real(7) - HINTLIB_MN sqrt (real(24)))
                                  / real(15));
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

     b0 (real(5 * sqr(dim) - 15 * dim + 14) / real(14)),
     b2 (real(50 - 25 * int (dim)) / real(168)),
     b4 (real(10 -  5 * int (dim)) / real (48))
{
   checkDimensionNotZero (dim);
   checkDimensionGeq<2> (dim);

#if HINTLIB_STATIC_WORKS == 0
   r = HINTLIB_MN sqrt (real(7) / real(15));
   s = HINTLIB_MN sqrt ((real(7) + HINTLIB_MN sqrt (real(24))) / real(15));
   t = HINTLIB_MN sqrt ((real(7) - HINTLIB_MN sqrt (real(24))) / real(15));
#endif
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

   return   num0_0    () * abs (b0)
      + 2 * numRR0_0s () * b1
      +     numR0_0fs () * (abs (b2) + 2 * abs (b4))
      + 2 * numRS0_0s () * b3;
}


/**
 *  Do the actual function evaluation
 */

L::real L::Rule5Stroud::eval (Integrand &f, const Hypercube &h)
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



