/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration
 *
 *  Copyright (C) 2007  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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

#include <HIntLib/productrule.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/hypercube.h>
#include <HIntLib/exception.h>
#include <HIntLib/integrand.h>

namespace L = HIntLib;


/**
 *  Constructor
 */

L::ProductRule::ProductRule (int _dim, int _numPoints)
   : dim (_dim),
     numPoints (_numPoints),
     p (_dim),
     counter (_dim, _numPoints)
{
   checkDimensionNotZero (dim);

   if (numPoints > 1)
   {
      int maxDim = logInt(std::numeric_limits<Index>::max(), Index(numPoints));

      if (dim > maxDim)  throw DimensionTooHigh (dim, maxDim);
   }
}


/**
 *  getNumPoints()
 */

L::Index
L::ProductRule::getNumPoints() const
{
   return powInt (Index(numPoints), dim);
}


/**
 *  isAllPointsInside()
 */

bool
L::ProductRule::isAllPointsInside() const
{
   for (int i = 0; i < numPoints; ++i)
   {
      if (offsets[i] < -1 || offsets[i] > 1)  return false;
   }

   return true;
}


/**
 *  eval()
 */

L::real
L::ProductRule::eval (Integrand &f, const Hypercube &h)
{
   const real* center = h.getCenter();
   const real* width  = h.getWidth();
   
   // numPoints == 1 must be treated separately.
   // In this case  dim  may be arbitrarily large, therefore the array
   // "weightProducts" set up below would overflow

   if (numPoints == 1)
   {
      for (int d = 0; d != dim; ++d)
      {
         p[d] = center[d] + width[d] * offsets[0];
      }

      return f(p) * HINTLIB_MN pow (weights[0], int(dim)) * h.getVolume();
   }

   // If we get here, numPoints >=2   and we need some non-trivial loops

   real sum = 0;

   switch (dim)
   {
      case 1:
         for (int i = 0; i < numPoints; ++i)
         {
            p[0] = center[0] + offsets[i] * width[0];
            sum += f (p) * weights[i];
         }
         break;

      case 2:
      {
         for (int i = 0; i < numPoints; ++i)
         {
            p[0] = center[0] + offsets[i] * width[0];
            real sum1 = 0;

            for (int j = 0; j < numPoints; ++j)
            {
               p[1] = center[1] + offsets[j] * width[1];
               sum1 += f (p) * weights[j];
            }

            sum += weights[i] * sum1;
         }
         break;
      }

      case 3:
      {
         for (int i = 0; i < numPoints; ++i)
         {
            p[0] = center[0] + offsets[i] * width[0];
            real sum1 = 0;

            for (int j = 0; j < numPoints; ++j)
            {
               p[1] = center[1] + offsets[j] * width[1];
               real sum2 = 0;

               for (int k = 0; k < numPoints; ++k)
               {
                  p[2] = center[2] + offsets[k] * width[2];
                  sum2 += f (p) * weights[k];
               }
               
               sum1 += weights[j] * sum2;
            }

            sum += weights[i] * sum1;
         }
         break;
      }

      default:
      {
         counter.reset();
         real weightProducts [std::numeric_limits<Index>::digits];
         weightProducts[dim] = 1.0;
         int pos = dim - 1;

         for (int d = 0; d != dim; ++d)
         {
            p[d] = center[d] + width[d] * offsets[0];
         }

         for(;;)
         {
            // Recalculate the weigth

            {
               real w = weightProducts[pos + 1];

               for (int i = pos; i >= 0; --i)
               {
                  weightProducts[i] = w *= weights[counter[i]];
               }
            }

            // evaluate f

            sum += weightProducts[0] * f(p);

            // get next Gray-Code

            unsigned newDigit;
            pos = counter.nextDigit(&newDigit);

            if (pos == dim)  break;

            // change a single coordinate in order to obtain the next point

            p[pos] = center[pos] + width[pos] * offsets[newDigit];
         }
      }
   }

   return sum * h.getVolume();
}


