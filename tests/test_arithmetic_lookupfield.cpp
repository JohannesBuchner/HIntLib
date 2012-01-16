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
#include <HIntLib/quotientfield.h>

#include <HIntLib/lookupfield.h>
#include <HIntLib/onedimvectorspace.h>

// The actual test templates

#include "test_arithmetic_tests.h"

namespace HIntLib
{

void testLookupField (unsigned size)
{
   typedef unsigned char u8;

   int maxDim   = logInt (256u, size);
   int maxDim32 = logInt (std::numeric_limits<unsigned>::max(), size);

   unsigned p;
   int e;

   if (Prime::isPrimePower (size, p, e))
   {
      // general case

      {
         typedef LookupField<u8> F;
         LookupGaloisField<u8> f (size);
         doTests (f, "LookupField<>");

         typedef PolynomialRing<F> PF;
         PF pf (f);
         doTests (pf, "PolynomialRing<LookupField<> >");

         doTests (PolynomialRing<PF> (pf, 'y'),
                  "PolynomialRing<PolynomialRing<LookupField<> > >");

         doTests (QuotientField<PF> (pf),
                 "QuotientField<PolynomialRing<LookupField<> > >");

         doTests (OneDimVectorSpace<F> (f),
                  "OneDimVectorSpace<LookupField<> >");

         for (int j = 1; j <= maxDim; ++j)
         {
            doTests (LookupVectorSpace<u8,u8>(f, j), "LookupVectorSpace<>");
         }
      }
     
      // optimized version for primes

      if (e == 1)
      {
         typedef LookupFieldPrime<u8> F;
         LookupGaloisFieldPrime<u8> f (size);
         doTests (f, "LookupFieldPrime<>");

         doTests (PolynomialRing<F> (f),
                  "PolynomialRing<LookupFieldPrime<> >");
         doTests (OneDimVectorSpace<F> (f),
                  "OneDimVectorSpace<LookupFieldPrime<> >");
      }
     
      // optimized version for powers of 2

      if (p == 2)
      {
         typedef LookupFieldPow2<u8> F;
         LookupGaloisFieldPow2<u8> f (size);
         doTests (f, "LookupFieldPow2<>");

         doTests (PolynomialRing<F> (f),
                  "PolynomialRing<LookupFieldPow2<> >");
         doTests (OneDimVectorSpace<F> (f),
                  "OneDimVectorSpace<LookupFieldPow2<> >");

         for (int j = 1; j <= maxDim; ++j)
         {
            doTests (LookupVectorSpacePow2<u8,u8>(f,j),
                     "LookupVectorSpacePow2<>");
         }

         for (int j = 1; j <= maxDim32; ++j)
         {
            doTests (VectorSpacePow2<u32>(f,j), "VectorSpacePow2<>");
         }
      }
   }
}

} // namespace HIntLib

