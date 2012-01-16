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
 *  rule5hammer.cpp
 *
 *  Cubature rule of degree 5 with 2*dim^2 + 1 points.
 *  All points are inside the hypercube.
 *
 *  This rule was published in
 *     P. C. Hammer and A. H. Stroud., Numerical Evaluation of Multiple
 *        Integrals. Math Tables Aids Comput, V 12 (1958), 272-280
 *
 *  It is also presented in
 *     A. H. Stoud.  Approximate Calculation of Multiple Integrals (1971)
 *  as formula Cn: 5-2
 */                                                                             

#ifdef __GNUG__
#pragma implementation
#endif

#include <HIntLib/rule5hammer.h>

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/hlmath.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;
using L::real;
using L::Index;


namespace
{
   const real b2 = real(25) / real(324);

   #if HINTLIB_STATIC_WORKS == 1
      const real r  = HINTLIB_MN sqrt (real(3) / real(5));
   #else
      real r;
   #endif
}

/**
 *  The constructor is used primarily to initialize all the dimension dependent
 *  constatns and to allocate (dimension dependent) memory
 */

L::Rule5Hammer::Rule5Hammer (unsigned d)
  : OrbitRule (d),
    lamda (d),
    b0 ((real(25) * sqr(d) - real(115) * d + real(162)) / real(162)),
    b1 ((real(70) - real(25) * d) / real(162))
{
   #if HINTLIB_STATIC_WORKS == 0
      r  = HINTLIB_MN sqrt (real(3) / real(5));
   #endif

   checkDimensionNotZero (d);
   checkDimensionGeq<2> (d);
}




Index L::Rule5Hammer::getNumPoints() const
{
   return num0_0 () + numR0_0fs () + numRR0_0fs ();
}                                             

real L::Rule5Hammer::getSumAbsWeight() const
{
   // Multiply each weight with the number of sampling points it is used for.
   // Don't forget to take the absolute value for weights that meight be
   // negative
 
   return num0_0    () * abs(b0)
        + numR0_0fs () * abs(b1)
        + numRR0_0fs() * b2;
}
 

/**
 *  Do the actual function evaluation
 */

L::real L::Rule5Hammer::eval (Integrand &f, const Hypercube &h)
{
   // Calculate lamdas

   const real* width  = h.getWidth();
   const real* center = h.getCenter();

   for (unsigned i = 0; i < dim; i++)  lamda[i] = width[i] * r;

   setCenter (center);

   return h.getVolume() *
   (
        b0 * eval0_0     (f)
      + b1 * evalR0_0fs  (f, center, lamda)
      + b2 * evalRR0_0fs (f, center, lamda)
   );
}

L::CubatureRuleFactory* L::Rule5Hammer::getFactory()
{
   return new DefaultCubatureRuleFactory<L::Rule5Hammer> ();
}


