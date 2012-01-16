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
 *  Rule 7-5 Genz
 *
 *  An embedded cubatare rule of degree 7 (embedded rule degree 5) due to
 *  A.Genz and Malik
 */


#ifdef __GNUG__
#pragma implementation
#endif

#include <HIntLib/rule75genzmalik.h>

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/hlmath.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;
using L::real;
using L::Index;

namespace
{
   /**
    *  Various dimension independent constants
    */

#if HINTLIB_STATIC_WORKS == 1
   const real lamda2 = HINTLIB_MN sqrt (real(9) / real(70));
   const real lamda4 = HINTLIB_MN sqrt (real(9) / real(10));
   const real lamda5 = HINTLIB_MN sqrt (real(9) / real(19));
#else
   real lamda2, lamda4, lamda5;
#endif

   const real weight2  = real(980) / real (6561);
   const real weight4  = real(200) / real(19683);

   const real weightE2 = real(245) / real (486);
   const real weightE4 = real (25) / real (729);
}


/**
 *  The constructor is used primarily to initialize all the dimension dependent
 *  constatns and to allocate (dimension dependent) memory
 */

L::Rule75GenzMalik::Rule75GenzMalik (unsigned dim)
   : OrbitRule (dim),
     widthLamda (dim),
     weight1 (real(12824 - 9120*int(dim) + 400*sqr(int(dim))) / real(19683)),
     weight3 (real( 1820 -  400*int(dim)) / real(19683)),
     weight5 (real(6859) / real(19683) / real(Index(1) << dim)),
     weightE1(real(  729 -  950*int(dim) +  50*sqr(int(dim))) / real(729)),
     weightE3(real(  265 -  100*int(dim)) / real(1458))
{
   checkDimensionNotZero (dim);
   checkDimensionGeq<2> (dim);
   checkDimensionLeq<std::numeric_limits<Index>::digits - 1> (dim);

#if HINTLIB_STATIC_WORKS == 0
   lamda2 = HINTLIB_MN sqrt (real(9) / real(70));
   lamda4 = HINTLIB_MN sqrt (real(9) / real(10));
   lamda5 = HINTLIB_MN sqrt (real(9) / real(19));
#endif
}


/**
 *  Return the sum of the absolute values of the weights
 */

real L::Rule75GenzMalik::getSumAbsWeight () const
{
   // Multiply each weight with the number of sampling points it is used for.
   // Don't forget to take the absolute value for weights that meight be
   // negative

   return num0_0     () * abs (weight1)
        + numR0_0fs  () * (weight2 + abs (weight3))
        + numRR0_0fs () * weight4
        + numR_Rfs   () * weight5;
}

Index L::Rule75GenzMalik::getNumPoints () const
{
   return num0_0 () + 2 * numR0_0fs () + numRR0_0fs () + numR_Rfs ();
}                       

/**
 *  Do the actual function evaluation
 */

unsigned L::Rule75GenzMalik::evalError
   (Integrand &f, const Hypercube &h, EstErr &ee)
{
   // Initialize

   const real* center = h.getCenter();
   const real* width  = h.getWidth();

   setCenter (center);

   Scaler scalerLamda2 (width, lamda2);

   for (unsigned i = 0; i != dim; ++i)  widthLamda [i] = width [i] * lamda4; 

   // Evaluate function in the center, in f (lamda2,0,...,0) and
   // f (lamda3=lamda4, 0,...,0)
   // Estimate dimension with largest error

   real sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;

   unsigned dimDiffMax =
      evalR0_0fs4d (f, center, sum1, scalerLamda2, sum2,
                    static_cast<real*> (widthLamda), sum3);

   // Calculate sum4 for f (lamda4, lamda4, 0, ...,0)

   real sum4 = evalRR0_0fs (f, center, widthLamda);

   // Calculate sum5 for f (lamda5, lamda5, ..., lamda5)

   for (unsigned i = 0; i != dim; ++i)  widthLamda [i] = width [i] * lamda5;

   real sum5 = evalR_Rfs (f, center, widthLamda);

   // Calculate fifth and seventh order result

   real result =
      h.getVolume() * (weight1 * sum1 + weight2 * sum2 + weight3 * sum3 +
                       weight4 * sum4 + weight5 * sum5);
   real res5th =
      h.getVolume() * (weightE1 * sum1 + weightE2 * sum2 +
                       weightE3 * sum3 + weightE4 * sum4);

   // calculate error

   ee.set (result, abs (res5th - result));

   return dimDiffMax;
}


/**
 *  getFactory()
 */

L::EmbeddedRuleFactory* L::Rule75GenzMalik::getFactory()
{
   return new DefaultEmbeddedRuleFactory<L::Rule75GenzMalik> ();
}

