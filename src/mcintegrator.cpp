/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration
 *
 *  Copyright (C) 2002,03,04,05,06  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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

#include <HIntLib/defaults.h>

#if defined HINTLIB_USE_INTERFACE_IMPLEMENTATION && ! defined HINTLIB_PARALLEL
#pragma implementation
#endif

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/mcintegrator.h>

#include <HIntLib/staticloadbalancer.h>
#include <HIntLib/pointset.h>
#include <HIntLib/integrand.h>
#include <HIntLib/hypercube.h>
#include <HIntLib/exception.h>

#ifdef HINTLIB_PARALLEL
#  define HINTLIB_NAME(x) x##StaticLB
#else
#  define HINTLIB_NAME(x) x
#endif


namespace L = HIntLib;


#ifndef HINTLIB_PARALLEL

/**
 *  setMinPercEval()
 */

L::MCIntegrator &
L::MCIntegrator::setMinPercEval (double perc)
{
   checkPercentageRange (perc);

   minPercEval = perc;

   return *this;
}


/**
 *  MCIntegratorJob
 *
 *  A ReportingJob that does integration while checking the estimated error.
 */

namespace HIntLib
{
   class MCIntegratorJob : public ReportingJob
   {
   public:
      MCIntegratorJob(Integrand&, real, real, Index);

      bool operator() (const real*);
      const StatisticVar<>* getStatistic() const  { return &stat; }
      Index getN() const  { return n; }

   private:
      StatisticVar<> stat;
      Integrand& f;
      Index n;
      Index nextCheck;
      const real reqAbsError;
      const real reqRelError;
   };
}

inline
L::MCIntegratorJob::MCIntegratorJob(
      Integrand& _f, real _reqAbsError, real _reqRelError, Index minN)
   : f(_f), n(0), nextCheck(minN),
     reqAbsError(_reqAbsError), reqRelError(_reqRelError)
{}

bool
L::MCIntegratorJob::operator() (const real* point)
{
   // Take sample and record it

   stat << f(point);

   // Is it time for checking the termination criteria?

   if (++n < nextCheck)  return false;

   // Estimate the current integration error

   const real error = stat.getStdDevSample() / HINTLIB_MN sqrt(real(n));
   
   // Check absolute error

   real factor;

   if (reqAbsError)
   {
      if (error < reqAbsError)  return true;
      factor = L::sqr(error / reqAbsError);
   }
   else
   {
      factor = 14.0;
   }

   // Check relative error

   if (reqRelError)
   {
      const real reqRelAbsError = abs(stat.getMean()) * reqRelError;
      if (error  < reqRelAbsError)  return true;
      factor = std::min(factor, L::sqr(error / reqRelAbsError));
   }

   // Correct factor if it is too high or too low

   // std::cout << "factor = " << factor << '\n';

   if (factor > 14.0)  factor = 8.0;
   else if (factor < 1.01)  factor = 1.01;
   else if (factor > 2.0)  factor = (factor - 2.0) / 2.0 + 2.0;

   nextCheck = Index (factor * n);

   // std::cout << n << " -> " << nextCheck << '\n';

   return false;
}
#endif


/**
 *  integrate()
 */

L::Integrator::Status
L::HINTLIB_NAME(MCIntegrator)::integrate (
   Integrand &f, const Hypercube &h, Index n, 
   real reqAbsError, real reqRelError, EstErr &ee)
{
   checkDimension(h, f);
   checkTerminationCriteria (n, reqAbsError, reqRelError,
#ifdef HINTLIB_PARALLEL
         true
#else
         false
#endif
      );

   Point point (h.getDimension());

   ps->randomize (getSeed());
   ps->setCube (&h);

#ifndef HINTLIB_PARLLEL
   if (reqAbsError == 0 && reqRelError == 0)
#endif
   {
      StaticLoadBalancer slb (n);

      StatisticVar<> stat;

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
              stat.getStdDevSample() * h.getVolume()
                    / HINTLIB_MN sqrt(real(n)));
   }
#ifndef HINTLIB_PARALLEL
   else
   {
      // Determine the minimum number of integrand evaluations

      Index minN;

      if (minNumEval == 0 && minPercEval == 0)
      {
         // The default is 3 times the dimension of the problem

         minN = 3 * h.getDimension();
      }
      else
      {
         minN = minNumEval;
         if (n != 0)
         {
            if (minPercEval >= 1)  minN = std::max(minN, n);
            else if(minPercEval > 0)
            {
               minN = std::max(minN, Index(n * minPercEval));
            }
         }
      }

      MCIntegratorJob job (f, reqAbsError / h.getVolume(), reqRelError, minN);

      ps->doJobRep (point, job,
                    (n == 0) ? std::numeric_limits<Index>::max() : n);

      ee.set (job.getStatistic()->getMean() * h.getVolume(),
              job.getStatistic()->getStdDevSample() * h.getVolume()
                    / HINTLIB_MN sqrt(real(job.getN())));
   }
#endif

   Status status = checkRequestedError (ee, reqAbsError, reqRelError);

   return (status == ERROR) ? MAX_EVAL_REACHED : status;
}

#undef HINTLIB_NAME

