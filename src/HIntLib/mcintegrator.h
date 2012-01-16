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

/**
 *  MCIntegrator
 *  MCIntegratorNoLB    (when compiled with PARALLEL defined)
 *
 *  The basic Monte Carlo integration routine
 *
 *  In Parllel mode, static load balancing ist provided
 */

// Check different flags, depending on parallel mode

#if defined(HINTLIB_PARALLEL) && !defined(HINTLIB_MCINTEGRATOR_MPI_H) || !defined(HINTLIB_PARALLEL) && !defined(HINTLIB_MCINTEGRATOR_H)


#include <HIntLib/integrator.h>


// Define Name macro and set flag according to parallel mode

#ifdef HINTLIB_PARALLEL
   #define HINTLIB_MCINTEGRATOR_MPI_H 1
   #define HINTLIB_NAME(x) x##StaticLB
#else
   #define HINTLIB_MCINTEGRATOR_H 1
   #define HINTLIB_NAME(x) x
#endif

/**
 *  MCIntegratorBase(NoLB)
 *
 *  The plain-vanilla Monte Carlo integration routine.
 *
 *  Depending on PARALLEL, a sequential algorithm or a parallel one using
 *  static load balancing is generated.
 */

namespace HIntLib
{
   class MultiPointSet;
   class PointSet;

   class HINTLIB_NAME(MCIntegrator) : public Integrator
   {
   private:
#if HINTLIB_PARALLEL
      typedef MultiPointSet PS;
#else
      typedef PointSet PS;
#endif

   public:
      HINTLIB_NAME(MCIntegrator) (PS* _ps) : ps(_ps) {}

      virtual
      Status integrate (
         Integrand &, const Hypercube &, Index maxEval,
         real reqAbsError, real reqRelError, EstErr &ee);

   private:
      PS* ps;
   };
}  // namespace HIntLib

#undef HINTLIB_NAME

#endif

