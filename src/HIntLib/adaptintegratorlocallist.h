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
 *  AdaptIntegratorLocalList
 *
 *  Parallel Adaptive Integration with local region collections
 */

#ifndef HINTLIB_ADAPTINTEGRATORLOCALLIST_H
#define HINTLIB_ADAPTINTEGRATORLOCALLIST_H 1

#include <HIntLib/defaults.h>

#ifndef HINTLIB_PARALLEL
#error AdaptIntegratorLocalList can only be used in PARALLEL mode
#endif

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/rulebasedintegrator.h>

namespace HIntLib
{
   class EmbeddedRuleFactory;
   class RegionCollection;

   class AdaptIntegratorLocalList : public EmbeddedRuleBasedIntegrator
   {
   public:
      AdaptIntegratorLocalList (
         const EmbeddedRuleFactory *,
         MPI_Comm = MPI_COMM_WORLD,
         unsigned frequency = 5000, unsigned maxDimension = 0);

      virtual
      Status integrate (
         Integrand &f, const Hypercube &h, Index maxEval,
         real reqAbsError, real reqRelError, EstErr &ee);

   private:
      MPI_Comm comm;
      enum TAGS {ERROR_CMP, REGIONS};

      unsigned frequency;

      unsigned dimension;
      int rank;
      int nodes;

      int size  [8];  // Size of mesh in a given direction
      int coord [8];  // position of itself in the mesh
      int succ  [8];
      int pred  [8];

      void doExchange2   (RegionCollection &rc, unsigned dir, unsigned dim); 
      void doExchange    (RegionCollection &rc, unsigned dir, unsigned dim);
      void doExchangeOld (RegionCollection &rc, unsigned dir, unsigned dim);
   };
}  // namespace HIntLib

#endif

