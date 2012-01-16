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

#include <HIntLib/rule5stroud2.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/exception.h>
#include <HIntLib/hypercube.h>


namespace L = HIntLib;

/**
 *  The constructor is used primarily to initialize all the dimension dependent
 *  constatns and to allocate (dimension dependent) memory
 */

L::Rule5Stroud2::Rule5Stroud2 (int dim)
   : OrbitRule (dim),
     alpha (HINTLIB_MN sqrt((real(5) * dim + 4) / 30)),
     beta  (HINTLIB_MN sqrt((real(5) * dim + 4) / (real(15) * dim - 12))),
     weightAlpha (real(40) / sqr(real(5) * dim + 4)),
     weightBeta  (sqr((real(5) * dim - 4) / (real(5) * dim + 4))
                    / (Index(1) << dim)),
     scaledBeta (dim)
{
   checkDimensionNotZero (dim);
   checkDimensionLeq<std::numeric_limits<Index>::digits - 1> (dim);
}


/**
 *  getNumPoints()
 */

L::Index
L::Rule5Stroud2::getNumPoints () const
{
   return numR0_0fs() + numR_Rfs();
}


/**
 *  isAllPointsInside()
 */

bool
L::Rule5Stroud2::isAllPointsInside() const
{
   return 2 <= dim && dim <= 5;
}


/**
 *  eval()
 */

L::real
L::Rule5Stroud2::eval (Integrand &f, const Hypercube &h)
{
   const real* center = h.getCenter();
   const real* width  = h.getWidth();

   setCenter(center);

   for (int d = 0; d < dim; ++d)  scaledBeta[d] = width[d] * beta;

   return h.getVolume() *
      (  weightAlpha * evalR0_0fs (f, center, Scaler(width, alpha))
       + weightBeta  * evalR_Rfs  (f, center, scaledBeta));
}


/**
 *  getFactory()
 */

L::CubatureRuleFactory*
L::Rule5Stroud2::getFactory()
{
   return new DefaultCubatureRuleFactory<L::Rule5Stroud2> ();
}


