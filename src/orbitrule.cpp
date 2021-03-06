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

/**
 *  Orbit Rule
 *
 *  Baseclass for all Cubature Rules using some kind of orbits as abscissa sets
 */

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/orbitrule.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/counter.h>

namespace L = HIntLib;


/**
 *  eval RRR 0_0 fs ()
 */

L::real
L::OrbitRule::evalRRR0_0fs (Integrand &f, const real*c, const real* r)
{
   real sum = 0;

   for (int i = 0; i != dim - 2; ++i)
   {
      p[i] = c[i] - r[i];

      for (int j = i + 1; j != dim - 1; ++j)
      {
         p[j] = c[j] - r[j];

         for (int k = j + 1; k != dim; ++k)
         {
            // Process 2^3 Points defined by changing the sign of r in
            // position i, j, and k.
            // This is done is a Gray Code ordering, so only one coordinate
            // has to be recalculated for every point

                                                p[k] = c[k]-r[k]; sum += f(p);
            p[i] = c[i]+r[i];                                     sum += f(p);
                              p[j] = c[j]+r[j];                   sum += f(p);
            p[i] = c[i]-r[i];                                     sum += f(p);
                                                p[k] = c[k]+r[k]; sum += f(p);
            p[i] = c[i]+r[i];                                     sum += f(p);
                              p[j] = c[j]-r[j];                   sum += f(p);
            p[i] = c[i]-r[i];                                     sum += f(p);

            p[k] = c[k];    // Done with k  ->  Restore p [k]
         }
         p[j] = c[j];     // Done with j  ->  Restore p [j]
      }
      p[i] = c[i];        // Done with i  ->  Restore p [i]
   }

   return sum;
}


/**
 *  eval RRRR 0_0 fs ()
 */

L::real
L::OrbitRule::evalRRRR0_0fs (Integrand &f, const real*c, const real* r)
{
   real sum = 0;

   for (int i = 0; i != dim - 3; ++i)
   {
      p[i] = c[i] - r[i];

      for (int j = i + 1; j != dim - 2; ++j)
      {
         p[j] = c[j] - r[j];

         for (int k = j + 1; k != dim - 1; ++k)
         {
            p[k] = c[k] - r[k];

            for (int l = k + 1; l != dim; ++l)
            {
               // Process 2^4=16 Points defined by changing the sign of r in
               // position i, j, k and l.
               // This is done in a Gray Code ordering, so only one coordinate
               // has to be recalculated for each point

               p[l] = c[l] - r[l]; sum += f(p);
               p[i] = c[i] + r[i]; sum += f(p);
               p[j] = c[j] + r[j]; sum += f(p);
               p[i] = c[i] - r[i]; sum += f(p);
               p[k] = c[k] + r[k]; sum += f(p);
               p[i] = c[i] + r[i]; sum += f(p);
               p[j] = c[j] - r[j]; sum += f(p);
               p[i] = c[i] - r[i]; sum += f(p);
               p[l] = c[l] + r[l]; sum += f(p);
               p[i] = c[i] + r[i]; sum += f(p);
               p[j] = c[j] + r[j]; sum += f(p);
               p[i] = c[i] - r[i]; sum += f(p);
               p[k] = c[k] - r[k]; sum += f(p);
               p[i] = c[i] + r[i]; sum += f(p);
               p[j] = c[j] - r[j]; sum += f(p);
               p[i] = c[i] - r[i]; sum += f(p);

               p[l] = c[l]; // Done with l  ->  Restore p [l]
            }
            p[k] = c[k];    // Done with k  ->  Restore p [k]
         }
         p[j] = c[j];       // Done with j  ->  Restore p [j]
      }
      p[i] = c[i];          // Done with i  ->  Restore p [i]
   }

   return sum;
}


/**
 *  eval R_R fs ()
 *
 *  Evaluate the integral on all 2^n points (+/-r,...+/-r)
 *
 *  A gray-code ordering is used to minimize the number of coordinate updates
 *  in p.
 */

L::real
L::OrbitRule::evalR_Rfs (Integrand &f, const real* c, const real* r)
{
   real sum = 0;

   // This bitmap is used to store the signes of the current point
   // 0 bit = +
   // 1 bit = -

   Index signs = 0;

   // We start with the point where r is ADDed in every coordinate
   //   (This implies signs=0)

   for (int i = 0; i != dim; ++i)   p [i] = c [i] + r [i];

   // Loop through the points in gray-code ordering

   for (Index i = 0; ; ++i)
   {
      // evaluate function at the current point

      sum += f (p);

      // determine, which sign has to be flipped

      int d = ls0(i);

      // Do we have this coordinate? If not, we are done.

      if (d >= dim)  break;

      // Swith the sign for this coordinate and remember it

      Index mask = Index(1) << d;

      signs ^= mask;

      // Use the new sign to determine if we do an addition or an subtraction

      p [d] = (signs & mask) ? c [d] - r [d] : c [d] + r [d];
   }

   return sum;
}


/**
 *  eval 3 pow S ()
 *
 *
 *  Evaluate  f  at the 3^dim points  (r_i_1,...,r_i_dim) with
 *
 *     r_0 = -r
 *     r_1 = 0
 *     r_2 = +r
 *
 *  and i_1,...,i_dim ranging independently over the integers 0, 1, and 2.
 *
 *  The implementation using Counter is not necessarily the most efficient;
 *  enumerating the abscissas in base-3 Gray-code order would be faster.
 */

L::real
L::OrbitRule::eval3powS (
      Integrand &f, const real* c, const real* r, real w0, real w1)
{
   real sum = 0;
   GrayCounter counter (dim, 3);
   real weights [std::numeric_limits<Index>::digits];
   weights[dim] = 1.0;
   setCenter (c);
   int pos = dim - 1;

   for(;;)
   {
      for (int i = pos; i >= 0; --i)
      {
         weights[i] = (counter[i] ? w1 : w0) * weights[i + 1];
      }

      sum += weights[0] * f(p);

      unsigned newDigit;
      pos = counter.nextDigit(&newDigit);

      if (pos == dim)  break;

      switch (newDigit)
      {
      case 0: p[pos] = c[pos];          break;
      case 1: p[pos] = c[pos] - r[pos]; break;
      case 2: p[pos] = c[pos] + r[pos]; break;
      }
   }

   return sum;
}

