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

#include <HIntLib/integerring.h>

// The actual test templates

#include "test_arithmetic_tests.h"

namespace HIntLib
{

void testInteger()
{
   // Z, Z[x], Z[x,y], Q, Q[x], Q(x)

   typedef IntegerRing<> Z;
   Z z;
   doTests (z, "IntegerRing<int>");

   typedef PolynomialRing<Z> PZ;
   PZ pz (z);
   doTests (pz, "PolynomialRing<IntegerRing<int> >");

   doTests (PolynomialRing<PZ> (pz, 'y'),
            "PolynomialRing<PolynomialRing<IntegerRing<int> > >");

   typedef QuotientField<Z> Q;
   Q q (z);
   doTests (q, "QuotientField<IntegerRing<int> >");

   typedef PolynomialRing<Q> PQ;
   PQ pq (q);
   doTests (pq, "PolynomialRing<QuotientField<IntegerRing<int> > >");

   doTests (QuotientField<PQ> (pq),
      "QuotientField<PolynomialRing<QuotientField<IntegerRing<int> > > >");

   // Q(i)

   PQ pq_i (q, 'i');
   PQ::type p2 = pq_i.x(2);
   p2[0] = q.makeElement (1);
   doTests (FactorField<PQ> (pq_i, p2),
      "FactorField<PolynomialRing<QuotientField<IntegerRing<int> > > >");

   // Q(Sqrt(2))

   PQ pq_r (q, 'r');
   PQ::type p3 = pq_r.x(2);
   p3[0] = q.makeElement(-2);
   doTests (FactorField<PQ> (pq_r, p3),
      "FactorField<PolynomialRing<QuotientField<IntegerRing<int> > > >");
}

} // namespace HIntLib

