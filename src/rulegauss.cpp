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

#include <HIntLib/rulegauss.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/hypercube.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;

/**
 *  Constructor
 *
 *  Calculate the abscissas and the weights.
 *
 *  This routine is based on Press et al., Numerical Recipies in C, Sect. 4.5
 */

L::RuleGauss::RuleGauss (int _dim, int _numPoints)
   : ProductRule(_dim, _numPoints),
     offsets(numPoints), weights(numPoints)
{
   init (offsets, weights);

   for (int i = 1; 2 * i - 1 <= numPoints; ++i)
   {
      // Calculate an approximation

      real z = HINTLIB_MN cos (Constants<real>::pi()
                                * (i - .25) / (numPoints + .5));
      real pp,z1;

      // Refine using Newton iteration

      do
      {
         real p1 = 1;
         real p2 = 0;
         for (int j = 1; j <= numPoints; ++j)
         {
            const real p3 = p2;
            p2 = p1;
            p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
         }
         pp = numPoints * (z * p1 - p2) / (z * z - 1.0);
         z1 = z;
         z = z1 - p1 / pp;
      }
      while (std::abs(z - z1) > std::numeric_limits<real>::epsilon());

      // store offset

      offsets[numPoints - i] = - (offsets[i - 1] = z);

      // calcualte and store weigth

      weights[i - 1] = weights[numPoints - i] = 1.0 / ((1.0 - z*z) * pp * pp);
   }
}


/**
 *  Declare a CubatureRuleFactory which is able to store the requested number
 *  of points.
 */

namespace
{
   class RuleGaussCubatureRuleFactory : public L::CubatureRuleFactory
   {
   public:
      RuleGaussCubatureRuleFactory (int _numPoints) : numPoints(_numPoints) {}
      
      virtual L::RuleGauss* create (int dim);
      virtual RuleGaussCubatureRuleFactory* clone () const;
   private:
      int numPoints;
   };

   L::RuleGauss*
   RuleGaussCubatureRuleFactory::create (int dim)
   {
      return new L::RuleGauss (dim, numPoints);
   }

   RuleGaussCubatureRuleFactory*
   RuleGaussCubatureRuleFactory::clone () const
   {
      return new RuleGaussCubatureRuleFactory (numPoints);
   }
}


/**
 *  getFactory()
 */

L::CubatureRuleFactory*
L::RuleGauss::getFactory (int numPoints)
{
   return new RuleGaussCubatureRuleFactory (numPoints);
}


