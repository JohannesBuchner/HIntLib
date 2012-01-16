/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002  Rudolf Sch�rer <rudolf.schuerer@sbg.ac.at>
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
 *  rule9stenger.cpp
 *
 *  Cubature rule of degree 9
 *  All points are inside the hypercube.
 *
 *  This rule was published in
 *     F. Stenger: Numerical Integration in N Dimensions
 *     M.S. theses, Univ. of Alberta, 1963
 *
 *  Listet as 9-1 in Stroud71
 */

#ifdef __GNUG__
#pragma implementation
#endif

#include <HIntLib/rule9stenger.h>

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/mymath.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;
using L::real;
using L::Index;


/**
 *  The constructor is used primarily to initialize all the dimension dependent
 *  constatns and to allocate (dimension dependent) memory
 */

L::Rule9Stenger::Rule9Stenger (unsigned dim)
   : OrbitRule (dim), aU (dim), aV (dim)
{
   // Check for valid dimension
   //   numRRRR0_0fs() becomes invalid for dim>216

   checkDimensionNotZero (dim);
   checkDimensionGeq<4> (dim);
   checkDimensionLeq<216> (dim);

   real u2 = 5.0 / 9.0 + sqrt (40.0 / 567.0);
   real v2 = 5.0 / 9.0 - sqrt (40.0 / 567.0);
   rU = sqrt (u2);
   rV = sqrt (v2);

   real u2MinV2 = u2 - v2;
   real u2DivV2 = u2 / v2;
   real u4 = sqr(u2); real v4 = sqr(v2);
   real u8 = sqr(u4); real v8 = sqr(v4);

   real dimMin1 = dim - 1.0;
   real dimMin2 = dim - 2.0;
   real dimMin3 = dim - 3.0;
 
   weightF = 1.0 / (525.0 * u2 * v2 * sqr(u2MinV2));
 
   weightH = (24.0 - 5.0 * dim) / (3240.0 * v8);
 
   weightI = (3.0 - 5.0 * v2) / (2160.0 * dimMin3 * u4 * u2 * u2MinV2);
 
   weightJ = (1.0  - 1296.0 * u8 * weightI) / (1296.0 * v8);
 
   weightE =   (3.0 - 5.0 * u2) / (-180.0 * v4 * u2MinV2)
             - weightF * u2DivV2
             - 2.0 * dimMin2 * (weightH + dimMin3 * weightJ);
 
   weightD =   (3.0 - 5.0 * v2) / (180.0 * u4 * u2MinV2)
             - weightF / u2DivV2 - 2.0 * dimMin2 * dimMin3 * weightI;
 
   weightC =   (3.0 - 5.0 * u2) / (-30.0 * v2 * u2MinV2)
             - 2 * dimMin1 * (weightE + weightF + dimMin2 *
                                (weightH + 2.0 / 3.0 * dimMin3 * weightJ));
 
   weightB =   (3.0 - 5.0 * v2) / (30.0 * u2 * u2MinV2)
             - 2.0 * dimMin1 * (weightD + weightF
                                  + 2.0 / 3.0 * dimMin2 * dimMin3 * weightI);
 
   weightA =
        1.0
      - 2.0 * dim * (  weightB + weightC
                     + dimMin1 * (  weightD + weightE + 2.0 * weightF
                                  +   1.0 / 3.0 * dimMin2
                                    * (  2 * weightH
                                       + dimMin3 * (weightI + weightJ)
                                      )
                                 )
                    ); 
#if 0
   real volume = (1 << dim);

   cout << "u " << rU << endl;
   cout << "v " << rV << endl;
   cout << "A " << weightA * volume << endl;
   cout << "B " << weightB * volume << endl;
   cout << "C " << weightC * volume << endl;
   cout << "D " << weightD * volume << endl;
   cout << "E " << weightE * volume << endl;
   cout << "F " << weightF * volume << endl;
   cout << "H " << weightH * volume << endl;
   cout << "I " << weightI * volume << endl;
   cout << "J " << weightJ * volume << endl;
#endif
}


Index L::Rule9Stenger::getNumPoints () const
{
   return   2 * numRRRR0_0fs () + numRRR0_0fs ()
          + 2 * numRR0_0fs ()   + numRS0_0fs ()
          + 2 * numR0_0fs ()    + num0_0 ();
}


real L::Rule9Stenger::getSumAbsWeight() const
{
   // Multiply each weight with the number of sampling points it is used for.
   // Don't forget to take the absolute value for weights that meight be
   // negative
#if 1
   return   numRRRR0_0fs () * (abs (weightI) + abs (weightJ))
          + numRRR0_0fs ()  *  abs (weightH)
          + numRR0_0fs ()   * (abs (weightD) + abs (weightE))
          + numRS0_0fs ()   *  abs (weightF)
          + numR0_0fs ()    * (abs (weightB) + abs (weightC))
          + num0_0 ()       *  abs (weightA);
#endif
#if 0
   return   numRRRR0_0fs () * ((weightI) + (weightJ))
          + numRRR0_0fs ()  *  (weightH)
          + numRR0_0fs ()   * ((weightD) + (weightE))
          + numRS0_0fs ()   *  (weightF)
          + numR0_0fs ()    * ((weightB) + (weightC))
          + num0_0 ()       *  (weightA); 
#endif
}
 

/**
 *  Do the actual function evaluation
 */

real L::Rule9Stenger::eval (Function &f, const Hypercube &h)
{
   // Calculate offsets from the center for points U and V

   const real* width = h.getWidth();
   for (unsigned i = 0; i != dim; ++i)
   {
      aU [i] = width [i] * rU;
      aV [i] = width [i] * rV;
   }

   // Evaluate integrand

   const real* c = h.getCenter();
   setCenter (c);

   return h.getVolume() *
   (
        weightJ * evalRRRR0_0fs (f, c, aV)
      + weightI * evalRRRR0_0fs (f, c, aU)
      + weightH * evalRRR0_0fs  (f, c, aV)
      + weightF * evalRS0_0fs   (f, c, aU, aV)
      + weightE * evalRR0_0fs   (f, c, aV)
      + weightD * evalRR0_0fs   (f, c, aU)
      + weightC * evalR0_0fs    (f, c, aV)
      + weightB * evalR0_0fs    (f, c, aU)
      + weightA * eval0_0       (f)
   ); 
}


/**
 *  getFactory()
 */

L::CubatureRuleFactory* L::Rule9Stenger::getFactory()
{
   return new DefaultCubatureRuleFactory<L::Rule9Stenger> ();
}
