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

#include <HIntLib/realfield.h>

// The actual test templates

#include "test_arithmetic_tests.h"

namespace HIntLib
{

void testComplex()
{
   // C, C[x], C(x)

   ComplexField<HIntLib::real> complexField;
   doTests (complexField, "ComplexField<real>");

   PolynomialRing<ComplexField<> > complexPolynomials (complexField);
   doTests (complexPolynomials, "PolynomialRing<ComplexField<real> >");

   doTests (QuotientField<PolynomialRing<ComplexField<HIntLib::real> > >
         (complexPolynomials),
           "QuotientField<PolynomialRing<ComplexField<real> > >");
}

void testReal()
{
   // R, R[x], R(x)

   RealField<HIntLib::real> realField;
   doTests (realField, "RealField<real>");

   PolynomialRing<RealField<> > realPolynomials (realField);
   doTests (realPolynomials, "PolynomialRing<RealField<real> >");

   doTests (QuotientField<PolynomialRing<RealField<HIntLib::real> > >
         (realPolynomials),
           "QuotientField<PolynomialRing<RealField<real> > >");
}

} // namespace HIntLib

