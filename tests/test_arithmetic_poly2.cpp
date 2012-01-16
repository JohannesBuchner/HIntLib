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

#include <HIntLib/quotientfield.h>

#include <HIntLib/polynomial2.h>
#include <HIntLib/gf2vectorspace.h>

// The actual test templates

#include "test_arithmetic_tests.h"

namespace HIntLib
{

void testPoly2()
{
   // GF[2], GF2(x)  using PolynomialRing2<>

   Polynomial2Ring<u32> gf2poly;
   doTests (gf2poly, "Polynomial2Ring<u32>");
   doTests (QuotientField<Polynomial2Ring<u32> > (gf2poly),
              "QuotientField<Polynomial2Ring<u32> >");
}

void testGF2Vec()
{
   // GF2^n

   for (unsigned d = 1; d <= 10; ++d)
   {
      doTests (GF2VectorSpace<u32> (d), "GF2VectorSpace<u32>");
   }
}

} // namespace HIntLib

