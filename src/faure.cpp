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

#ifdef __GNUG__
#pragma implementation
#endif

#include <HIntLib/faure.h>

#include <HIntLib/modulararithmetic.h>
#include <HIntLib/prime.h>
// #include <HIntLib/exception.h>

namespace L = HIntLib;

/**
 *  Faure
 *
 *  A  Heap Allocated Generator Matrix  initialized according to Faure.
 *
 *  For details, see
 *
 *  [1] Henry Faure.  Discrépance de suites associées à un système de
 *      numération (en dimension $s$).  Acta Arithmetica, 41: 337-351, 1982.
 *  [2] Bennett L. Fox.  Algorithm 647: Implementation and Relative Efficiency
 *      of Quasirandom Sequence Generators.  ACM TOMS, 12(4):362-376, 1986.
 *
 *  The algorithm used here is quite different from [2], becuase we calculate
 *  the full Generator Matrix apriori, instead of delaying some calculations
 *  until the sequence is actually generated.
 */

/**
 *  Consturctor
 */

L::Faure::Faure (unsigned _dim)
   : HeapAllocatedGeneratorMatrixGen<unsigned char> (Prime::next(_dim), _dim)
{
   init (*this);
}

L::Faure::Faure (unsigned _dim, unsigned _m, unsigned _prec)
   : HeapAllocatedGeneratorMatrixGen<unsigned char>
       (Prime::next (_dim), _dim, _m, _prec)
{
   init (*this);
}


/**
 *  init()
 */

void L::Faure::init (MutableGeneratorMatrixGen<unsigned char> &gm)
{
   // Matrix d=0 ist the identity matrix

   if (gm.getDimension() < 1)  return;
   gm.makeIdentityMatrix (0);

   ModularIntegerRing<unsigned char> a (gm.getBase());

   if (gm.getDimension() < 2)  return;

   // Create Pascal's Triangle in Matrix d=1

   gm.set (1, 0, 0, a.one ());

   for (unsigned r = 1; r < gm.getM(); ++r)
   {
      gm.set (1, r, 0, a.one());
      if (r < gm.getPrecision())  gm.set (1, r, r, a.one());

      for (unsigned b = 1; b < std::min (r, gm.getPrecision()); ++b)
      {
         gm.set (1, r, b, a.add (gm (1, r-1, b-1), gm (1, r-1, b)));
      }
   }
   
   // initialize Faure-matrices for d =  2,...,dim-1
 
   for (unsigned d = 2; d < gm.getDimension(); ++d)
   {
      for (unsigned r = 0; r < gm.getM(); ++r)
      {
         for (unsigned b = 0; b < std::min (r, gm.getPrecision()); ++b)
         {
            unsigned char x = gm (1, r, b);
            for (unsigned i = 0; i < (r-b); ++i)  a.mulBy (x, d);
            gm.set (d, r, b, x);
         }

         if (r < gm.getPrecision())  gm.set (d, r, r, a.one());
      }
   }
}

