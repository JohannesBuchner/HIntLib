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

#include <HIntLib/niederreitermatrixgen.h>

#include <HIntLib/modulararithmetic.h>
#include <HIntLib/lookupfield.h>
#include <HIntLib/prime.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;

/**
 *   Constructor
 *
 *   For each dimension, determine an irreducible, monic polynomial and call
 *   init() to initialize the initialize the corresponding matrix.
 */

template<class A>
L::NiederreiterMatrixGen<A>::NiederreiterMatrixGen
   (const A& a, unsigned _dim, unsigned _m, unsigned _prec)
   : HeapAllocatedGeneratorMatrixGen<typename A::type>
      (a.size(), _dim, _m, _prec)
{
   init (*this, a);
}


/**
 *  init()  -  whole matrix
 */

template<class A>
void L::NiederreiterMatrixGen<A>::init
   (MutableGeneratorMatrixGen<typename A::type> &gm, A a)
{
   // Find _dim_ monic, irreducible polynomials

   PolyRing poly (a);
   Poly p;
   unsigned i = 0;

   for (unsigned d = 0; d < gm.getDimension(); ++d)
   {
      do
      {
         p = poly.element (i++);
      } while (! (poly.isMonic (p) && poly.isPrime (p)));

      init (gm, a, d, p);
   }
}


/**
 *   Calculation of vectors v according to Bratley-Fox based on a given
 *   polynomial.
 */

template<class A>
void L::NiederreiterMatrixGen<A>::init
   (MutableGeneratorMatrixGen<typename A::type> &gm,
    A a, unsigned d, const Poly &irred)
{
   T v [100];
   PolyRing poly (a);

   // cout << "Using " << irred << " for d=" << d << endl;

   const int degree = irred.degree ();

   Poly newPoly = poly.one();

   int u = 0;

   for (unsigned j = 0; j < gm.getPrec(); j++)
   {
      // cout << "  j=" << j << endl;
      // Do we need a new v?

      if (u == 0)
      {
         Poly oldPoly = newPoly;
         unsigned oldDegree = oldPoly.degree ();

         // calculate polyK+1 from polyK

         poly.mulBy (newPoly, irred);
         int newDegree = newPoly.degree ();
         // cout << "    newPolynomial: " << newPoly << endl;

         // kj can be set to any value between 0 <= kj < newDegree

         const unsigned kj = oldDegree   // proposed by BFN
                             // newDegree - 1    // standard, bad???
                             // 0 
                             // (newDegree > 3) ? 3 : oldDegree 
         ;

         std::fill (&v[0], &v[kj], a.zero());  // Set leading v's to 0

         v[kj] = a.one();                     // the next one is 1

         if (kj < oldDegree)
         {
            T term = oldPoly [kj];

            for (unsigned r = kj + 1; r < oldDegree; ++r)
            {
               v [r] = a.one (); // 1 is arbitrary. Could be 0, too

               a.addTo (term, a.mul (oldPoly [r], v [r]));
            }

            // set v[] not equal to -term

            v [oldDegree] = a.sub (a.one(), term);

            for (int r = oldDegree + 1; r < newDegree; ++r) v [r] = a.one(); //or 0
         }
         else
         {
            for (int r = kj + 1; r < newDegree; ++r) v [r] = a.one(); // or 0..
         }

         // All other elements are calculated by a recursion parameterized
         // by polyK

         for (unsigned r = 0; r < 100u - newDegree; ++r)
         {
            T term = a.zero();

            for (int i = 0; i < newDegree; ++i)
            {
               a.addTo (term, a.mul (newPoly [i], v [r+i]));
            }
            v [newDegree + r] = a.neg (term);
         } 
      }

      // Set data in ci

      // cout << "v =";
      for (unsigned r = 0; r < gm.getM(); ++r)
      {
         // cout << " " << v[r+u];
         gm.setd (d,r,j, v[r+u]);
      }
      // cout << endl;

      if (++u == degree) u = 0;
   } 
}


/**
 *  Niederreiter Matirx PP
 */

/**
 *  init ()
 */

void L::NiederreiterMatrixPP::init (unsigned char prime, unsigned power)
{
   if (power == 0)  throw GaloisFieldExponent();
   else if (power == 1)
   {
      ModularArithField<unsigned char> a (prime);
      NiederreiterMatrixGen<ModularArithField<unsigned char> >::init (*this, a);
   }
   else  // general case (power > 1)
   {
      LookupGaloisField<unsigned char> field (prime, power);
      NiederreiterMatrixGen<LookupField<unsigned char> >
         ::init (*this, field);
   }
}

void L::NiederreiterMatrixPP::init (unsigned char size)
{
   if (Prime::test (size))  init (size, 1);
   else
   {
      LookupGaloisField<unsigned char> field (size);
      NiederreiterMatrixGen<LookupField<unsigned char> >
         ::init (*this, field);
   }
}

namespace HIntLib
{
   template class NiederreiterMatrixGen<ModularArithField<unsigned char> >;
   template class NiederreiterMatrixGen<LookupField<unsigned char> >;
}

