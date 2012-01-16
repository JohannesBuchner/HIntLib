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

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/niederreitermatrixgen.h>

#include <HIntLib/generatormatrixgen.h>
#include <HIntLib/polynomial.h>
#include <HIntLib/lookupfield.h>
#include <HIntLib/prime.h>
#include <HIntLib/array.h>

namespace L = HIntLib;


/**
 *  init Niederreiter ()
 *
 *  Initialize whole matrix, given a certain arithmetic
 */

template<class A>
void L::initNiederreiter (GeneratorMatrixGen<typename A::type> &gm, A a)
{
   typedef typename A::type T;
   typedef PolynomialRing<A> PolyRing;
   typedef typename PolyRing::type Poly;

   // Find _dim_ monic, irreducible polynomials

   PolyRing poly (a);
   typename PolyRing::PrimeGenerator ig (poly);

   for (unsigned d = 0; d < gm.getDimension(); ++d)
   {
      initNiederreiter (gm, a, d, ig.next());
   }
}


/**
 *  init Niederreiter ()
 *
 *  Initialized a single coordiante with a certain irreducible polynomial
 *  over a given arithmetic.
 */

template<class A>
void L::initNiederreiter
  (GeneratorMatrixGen<typename A::type> &gm, A a, unsigned d,
   const typename PolynomialRing<A>::type &irred)
{
   typedef typename A::type T;
   typedef PolynomialRing<A> PolyRing;
   typedef typename PolyRing::type Poly;

   const int degree = irred.degree ();

   const unsigned vSize
      = std::max (gm.getM() + degree - 1,   // these elements are copied to gm
                  (gm.getTotalPrec() + degree + 1));  // used in the loop

   Array<T> v (vSize);
   PolyRing poly (a);

   // cout << "Using " << irred << " for d=" << d << endl;

   Poly newPoly = poly.one();

   int u = 0;

   for (unsigned j = 0; j < gm.getTotalPrec(); j++)
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

         std::fill (&v[0], &v[kj], T());  // Set leading v's to 0

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

         for (unsigned r = 0; r < vSize - newDegree; ++r)
         {
            T term = T();

            for (int i = 0; i < newDegree; ++i)
            {
               a.addTo (term, a.mul (newPoly [i], v [r+i]));
            }
            v [newDegree + r] = a.neg (term);
         }
      }

      // Set data in ci

      for (unsigned r = 0; r < gm.getM(); ++r)  gm.setd (d,r,j, v[r+u]);

      if (++u == degree) u = 0;
   }
}


/**
 *  init Niederreiter ()
 *
 *  Figures out the proper arithmetic to initialize the Niederreiter Matrix
 */

void L::initNiederreiter (GeneratorMatrixGen<unsigned char> &gm)
{
   const unsigned base = gm.getBase();

   if (Prime::test (base))
   {
      LookupGaloisFieldPrime<unsigned char> field (base);
      initNiederreiter (
            gm, static_cast<LookupFieldPrime<unsigned char>&> (field));
   }
   else
   {
      LookupGaloisField<unsigned char> field (base);
      initNiederreiter (gm, static_cast<LookupField<unsigned char>&> (field));
   }
}


#if 0
/**
 *  This function fills the passed array of polynomials with the first MAX_DIM
 *  irreducible polynomials.
 *
 *  No primitive polynomials are required, so a simple sieve can be used.
 */

void createIrredPolys (Poly* irredPolys)
{
   // Use the following map for Eratosthenes' sieve

   unsigned mapSize = 100;
   bool* map = 0;

   unsigned int count;

   do
   {
      mapSize *= 2;

      // dump old table and create a new one

      delete[] map;
      map = new bool [mapSize];

      // find irreducible polynomials

      eratosthenes (map, mapSize, Poly (0));

      // Count irreducible polynomials in this map

      count = 0;

      for (unsigned i = 0; i < mapSize; i++) if (map [i]) count++;

   } while (count < NiederreiterMatrix::MAX_DIM);  // do we have enough?

   // Copy the first MAX_DIM irreducible polynomials to irredPolys

   bool* p = map;

   for (unsigned i = 0; i < NiederreiterMatrix::MAX_DIM; i++)
   {
      while (! *p) p++;

      irredPolys [i] = Poly (p++ - map);
   }

   // Free the memory used for the map

   delete[] map;
}
#endif

