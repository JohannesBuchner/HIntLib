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

// algebraic structures to be tested

#include <HIntLib/polynomial.h>
#include <HIntLib/quotientfield.h>
#include <HIntLib/factorring.h>

#include <HIntLib/gf2.h>

// The actual test templates

#include "test_arithmetic_tests.h"

namespace HIntLib
{

void testGF2()
{
   // GF2, GF2[x], GF2(x), GF2(x)[y], GF(x,y)

   GF2 gf2;
   doTests (gf2, "GF2");

   typedef PolynomialRing<GF2> P;
   P gf2poly (gf2);
   doTests (gf2poly, "PolynomialRing<GF2>");

   typedef QuotientField<P> Q;
   Q gf2quot (gf2poly);
   doTests (gf2quot, "QuotientField<PolynomialRing<GF2> >");

   typedef PolynomialRing<Q> PQ;
   PQ gf2quotpoly (gf2quot, 'y');
   doTests (gf2quotpoly, "PolynomialRing<QuotientField<PolynomialRing<GF2> > >");

   const unsigned char coeff [] = {0, 1, 1};
   const P::type p1 (coeff, coeff + 3);
      
   PQ::type p2 = gf2quotpoly.x(3);
   p2[2] = gf2quot.one();
   p2[0] = gf2quot.makeElement(p1);

   typedef FactorField<PQ> FF;
   FF gf2ff (gf2quotpoly, p2);
   doTests (gf2ff, "FactorField<PolynomialRing<QuotientField<PolynomialRing<GF2> > > >");
}

} // namespace HIntLib

