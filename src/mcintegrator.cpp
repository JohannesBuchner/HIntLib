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
 *  MCIntegratorBase
 *
 *  Protected base class for all MCIntegrators.
 *
 *  It provides a template function for integration with any Prseudo Random
 *  Number Generator
 */

#if defined __GNUG__ && ! defined HINTLIB_PARALLEL
#pragma implementation
#endif

#include <HIntLib/mcintegrator.h>

#include <HIntLib/staticloadbalancer.h>
#include <HIntLib/exception.h>
#include <HIntLib/pointset.h>
#include <HIntLib/function.h>
#include <HIntLib/hypercube.h>

#ifdef HINTLIB_PARALLEL
   #define HINTLIB_NAME(x) x##StaticLB
#else
   #define HINTLIB_NAME(x) x
#endif


namespace L = HIntLib;

L::Integrator::Status L::HINTLIB_NAME(MCIntegrator)::integrate (
   Function &f, const Hypercube &h, Index n,
   real reqAbsError, real reqRelError, EstErr &ee)
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

   StaticLoadBalancer slb (n);

   StatisticVar<> stat;

   Array<real> point (h.getDimension());

   ps->setCube (&h);
   #ifdef HINTLIB_PARALLEL
   ps->select (slb.getRank(), slb.getSize());
   #endif
   ps->integrate (point, f, slb.getRange(), stat);
   
   // In parallel mode, collect results from all nodes

   #ifdef HINTLIB_PARALLEL
      stat.reduce();

      if (slb.getRank() > 0) return WRONG_NODE;
   #endif

   ee.set (stat.getMean() * h.getVolume(),
           stat.getStdDevSample() * h.getVolume() / sqrt(real(n)));

   Status status = checkRequestedError (ee, reqAbsError, reqRelError);

   return (status == ERROR) ? MAX_EVAL_REACHED : status;
} 

#undef HINTLIB_NAME

