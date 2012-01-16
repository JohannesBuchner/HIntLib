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

#include <HIntLib/rule5mustardlynessblatt.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/defaultcubaturerulefactory.h>
#include <HIntLib/exception.h>
#include <HIntLib/hypercube.h>


namespace L = HIntLib;

namespace
{
#if HINTLIB_STATIC_WORKS == 1
   const L::real alpha = HINTLIB_MN sqrt (L::real(2) / L::real(5));
#else
   L::real alpha;
#endif

   const L::real weightAlpha = L::real(5) / L::real(18);
}


/**
 *  The constructor is used primarily to initialize all the dimension dependent
 *  constatns and to allocate (dimension dependent) memory
 */

L::Rule5MustardLynessBlatt::Rule5MustardLynessBlatt (int dim)
   : OrbitRule (dim),
     weight0 ((8 - 5 * real(dim)) / 9),
     weight1 (real(1) / 9 / (Index(1) << dim))
{
   checkDimensionNotZero (dim);
   checkDimensionLeq<std::numeric_limits<Index>::digits - 1> (dim);

#if HINTLIB_STATIC_WORKS == 0
   alpha = HINTLIB_MN sqrt (real(2) / real(5));
#endif
}


/**
 *  getNumPoints()
 */

L::Index
L::Rule5MustardLynessBlatt::getNumPoints () const
{
   return num0_0() + numR0_0fs() + numR_Rfs();
}


/**
 *  getSumAbsWeight()
 */

L::real
L::Rule5MustardLynessBlatt::getSumAbsWeight() const
{
   return weight0 * num0_0() + weightAlpha * numR0_0fs() + weight1 * numR_Rfs();
}


/**
 *  eval()
 */

L::real
L::Rule5MustardLynessBlatt::eval (Integrand &f, const Hypercube &h)
{
   const real* center = h.getCenter();
   const real* width  = h.getWidth();

   setCenter(center);

   return h.getVolume() *
      (  weight0     * eval0_0    (f)
       + weightAlpha * evalR0_0fs (f, center, Scaler(width, alpha))
       + weight1     * evalR_Rfs  (f, center, width));
}


/**
 *  getFactory()
 */

L::CubatureRuleFactory*
L::Rule5MustardLynessBlatt::getFactory()
{
   return new DefaultCubatureRuleFactory<L::Rule5MustardLynessBlatt> ();
}


