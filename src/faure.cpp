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

#define HINTLIB_LIBRARY_OBJECT

#include <algorithm>

#include <HIntLib/faure.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/generatormatrixgen.h>
#include <HIntLib/lookupfield.h>

namespace L = HIntLib;

/**
 *  Faure
 *
 *  A Generator Matrix  initialized according to Faure.
 *
 *  For details, see
 *
 *  [1] Henry Faure.  Discr'epance de suites associ'ees `a un syst`eme de
 *      num'eration (en dimension $s$).  Acta Arithmetica, 41: 337-351, 1982.
 *  [2] Bennett L. Fox.  Algorithm 647: Implementation and Relative Efficiency
 *      of Quasirandom Sequence Generators.  ACM TOMS, 12(4):362-376, 1986.
 *
 *  The algorithm used here is quite different from [2], because we calculate
 *  the full Generator Matrix apriori, instead of delaying some calculations
 *  until the sequence is generated.
 */


/**
 *  init Faure ()
 */

void L::initFaure (GeneratorMatrixGen<unsigned char> &gm)
{
   // Matrix d=0 ist the identity matrix

   if (gm.getDimension() < 1)  return;
   gm.makeIdentityMatrix (0);

   if (gm.getDimension() < 2)  return;

   LookupGaloisField<unsigned char> a (gm.getBase());

   // Create Pascal's Triangle in Matrix d=1

   for (int r = 0; r < gm.getM(); ++r)  gm.setd (1, r, 0, a.one());

   for (int b = 1; b < gm.getPrec(); ++b)
   {
      unsigned char x = 0;

      for (int r = b; r < gm.getM(); ++r)
      {
         a.addTo (x, gm (1, r-1, b-1));
         gm.setd (1,r,b, x);
      }
   }

   // initialize Faure-matrices for d =  2,...,dim-1

   for (int d = 2; d < gm.getDimension(); ++d)
   {
      for (int r = 0; r < gm.getM(); ++r)
      {
         for (int b = 0; b < std::min (r, gm.getPrec()); ++b)
         {
            gm.setd (d,r,b, a.mul (gm (1,r,b), a.power (d, r-b)));
         }

         if (r < gm.getPrec())  gm.setd (d, r, r, a.one());
      }
   }
}

