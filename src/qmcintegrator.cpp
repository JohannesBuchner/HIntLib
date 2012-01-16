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
 *  QMCIntegratorBase
 *
 *  Protected base class for all QMCIntegrators.
 *
 *  It provides a template function for integration with any Prseudo Random
 *  Number Generator
 */

#if defined __GNUG__ && ! defined HINTLIB_PARALLEL
#pragma implementation
#endif

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/qmcintegrator.h>

#include <HIntLib/staticloadbalancer.h>
#include <HIntLib/pointset.h>
#include <HIntLib/exception.h>
#include <HIntLib/statistic.h>
#include <HIntLib/integrand.h>
#include <HIntLib/hypercube.h>


#ifdef HINTLIB_PARALLEL
   #define NAME(x) x##StaticLB
#else
   #define NAME(x) x
#endif

namespace L = HIntLib;

L::Integrator::Status L::NAME(QMCIntegrator)::integrate (
   Integrand &f,
   const Hypercube &h,
   Index n,
   real, real,
   EstErr &ee)
{
   checkDimension(h, f);

   if (n == 0)
   {
      #ifdef HINTLIB_NO_EXCEPTIONS
         ee.set (0.0, 0.0);
         return ERROR;
      #else
         throw MaxEvaluationsRequired();
      #endif
   }

   // Decrease the number of points, to get a "good" number

   ps->setCube (&h);
   ps->randomize (getSeed());
   n = ps->getOptimalNumber (n, h);
   Statistic<> stat;
   Point point (h.getDimension());

   #ifdef HINTLIB_PARALLEL
      StaticLoadBalancer slb (n);

      ps->integratePartition (
         point, f, slb.getRange(), slb.getBegin(), slb.getEnd(), stat);

      stat.reduce();

      if (slb.getRank() > 0) return WRONG_NODE;
   #else
      ps->integrate(point, f, n, stat);
   #endif

   ee.setNoErr (stat.getMean() * h.getVolume());

   return MAX_EVAL_REACHED;
}

#undef NAME

