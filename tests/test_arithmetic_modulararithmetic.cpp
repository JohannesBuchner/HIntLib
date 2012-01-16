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

// algebraic structures to be tested

#include <HIntLib/polynomial.h>

#include <HIntLib/integerring.h>
#include <HIntLib/galoisfield.h>

// The actual test templates

#include "test_arithmetic_tests.h"

namespace HIntLib
{

void testModularArithmetic (unsigned size)
{
   typedef unsigned char u8;
   typedef IntegerRing<> Z;
   Z z;

   unsigned maxDim32 = logInt (std::numeric_limits<unsigned>::max(), size);

   // Only prime numbers give a field

   if (Prime::test (size))
   {
      typedef ModularArithmeticField<u8> F;
      F f (size);
      doTests (f, "ModularArithmeticField<>");

      typedef PolynomialRing<F> PF;
      PF pf (f);
      doTests (pf, "PolynomialRing<ModularArithmeticField<>");
      
      // Factor fields modulo irreducible polynomials

      for (unsigned j = 1; j <= std::min (8u, maxDim32); ++j)
      {
         doTests (GaloisField<u8> (size, j), "GaloisField<>");
      }

      // Factor rings modulo arbitrary polynomials

      for (unsigned deg = 1; deg <= std::min (5u, maxDim32); ++deg)
      {
         for (unsigned j = 0;
              j < std::min (4u, (size - 1) * powInt (size,deg)); ++j)
         {
            FactorRing<PF> fpf (pf, pf.element (powInt (size, deg) + j));
            doTests (fpf,
               "FactorRing<PolynomialRing<ModularArithmeticField<> > >");
         }
      }

      FactorField<Z> f2 (z, size);
      doTests (f2, "FactorField<IntegerRing<> >");
   }

   // A ring can be created based on prime and non-prime elements

   {
      typedef ModularArithmeticRing<u8> R;
      R r (size);
      doTests (r, "ModularArithmeticRing<>");

      doTests (PolynomialRing<R> (r),
               "PolynomialRing<ModularArithmeticRing<> >");

      FactorRing<Z> r2 (z, size);
      doTests (r2, "FactorRing<IntegerRing<> >");
   }
}


void testModularArithmeticShort (unsigned size)
{
   ModularArithmeticField<unsigned short> f (size);
   doTests (f, "ModularArithmeticField<unsigned short>");
}

} // namespace HIntLib

