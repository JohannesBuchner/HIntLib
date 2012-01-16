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
 *  Monte Carlo integration
 *
 *  Templates for executing the inner loop of Monte Carlo integration routines
 */

#ifndef HINTLIB_MC_INTEGRATION_H
#define HINTLIB_MC_INTEGRATION_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/hypercube.h>
#include <HIntLib/distribution.h>
#include <HIntLib/pointset.h>
#include <HIntLib/cubepartitioner.h>
#include <HIntLib/integrand.h>


namespace HIntLib
{

class Hypercube;

/**
 *  mcIntegration()
 */

template<class Gen, class Stat>
inline
void mcIntegration (
   real point [], Gen &gen, const Hypercube &h, Integrand &f,
   typename Stat::CounterType n, Stat &stat)
{
   UniformCube<Gen> uc (gen, h);

   while (stat.getCount() != n)
   {
      uc (gen, point);

      stat << f(point);
   }
}

/**
 *  mcDoJob
 */

template<class Gen>
inline
void mcDoJob (real point [], Gen &gen, const Hypercube &h, Job &job, Index n)
{
   UniformCube<Gen> uc (gen, h);

   for (Index i = 0; i < n; ++i)
   {
      uc (gen, point);
      job (point);
   }
}

/**
 *  mcDoJobRep
 */

template<class Gen>
inline
bool mcDoJobRep (
   real point [], Gen &gen, const Hypercube &h, ReportingJob &job, Index n)
{
   UniformCube<Gen> uc (gen, h);

   for (Index i = 0; i < n; ++i)
   {
      uc (gen, point);
      if (job (point))  return true;
   }

   return false;
}

namespace Private
{
   /**
    *  Integration CP
    */

   template<class Gen, class Stat>
   class IntegrationCP : public CubePartitioner
   {
   protected:
      IntegrationCP (Integrand &f, real* p, Gen& gen)
         : f(f), p(p), gen (gen) {}
 
   protected:
      Integrand &f;
      real *p;
      Gen &gen;
 
   public:
      Stat stat;
   };
   
   /**
    *  Stratified Integration CP
    */

   template<class Gen, class Stat>
   class StratifiedIntegrationCP : public IntegrationCP<Gen,Stat>
   {
   public:
      StratifiedIntegrationCP (Integrand &f, real* point, Gen& gen)
         : IntegrationCP<Gen,Stat> (f, point, gen) {}
 
      void action (const Hypercube &h)
      {
         uniform (this->gen, h, this->p);
         this->stat << this-> f(this->p);
      }
   };

   /**
    *  AntitheticIntegrationCP
    */

   template<class Gen, class Stat>
   class AntitheticIntegrationCP : public IntegrationCP<Gen,Stat>
   {
   public:
      AntitheticIntegrationCP (Integrand &f, real* point, Gen& gen)
         : IntegrationCP<Gen,Stat> (f, point, gen) {}
 
      void action (const Hypercube &h)
      {
         uniform (this->gen, h, this->p);
         this->stat << this-> f(this->p);
 
         for (unsigned i = 0; i < h.getDimension(); ++i)
         {
            this->p[i] = 2 * h.getCenter(i) - this->p[i];
         }
         this->stat << this->f (this->p);
      }
   };

}  // namespace Private


/**
 *  stratifedIntegration()
 *
 *  Basic routine for stratified integration.  Should only be called directly
 *  by other integration routines.
 */

template<class Gen, class Stat>
inline
void stratifiedIntegration (
   real* point, Gen &gen, const Hypercube &h, Integrand &f, Index n,
   Stat &stat)
{
   Private::StratifiedIntegrationCP<Gen,Stat> cp (f, point, gen);
   cp (h, n, cp.EQUIVOLUME);
   stat = cp.stat;
}


/**
 *  antitheticIntegration()
 */

template<class Gen, class Stat>
inline
void antitheticIntegration (
   real* point, Gen &gen, const Hypercube &h, Integrand &f, Index n,
   Stat &stat)
{
   // Use CubePartitioner to perform action() an every sub-cube
 
   Private::AntitheticIntegrationCP<Gen,Stat> cp (f, point, gen);
 
   cp (h, n / 2, cp.EQUIVOLUME);
 
   stat = cp.stat;
 
   // if n is odd, draw one additional sample
 
   if (n & 1)
   {
      uniform (gen, h, point);
      stat << f(point);
   }
}

namespace Private
{
   /**
    *  Do Job CP
    */

   template<class Gen>
   class DoJobCP : public CubePartitioner
   {
   protected:
      DoJobCP (Job &_job, real* _point, Gen& _gen)
         : job(_job), point(_point), gen(_gen) {}
 
   protected:
      Job &job;
      real *point;
      Gen &gen;
   };
   
   /**
    *  Stratified do Job CP
    */

   template<class Gen>
   class StratifiedDoJobCP : public DoJobCP<Gen>
   {
   public:
      StratifiedDoJobCP (Job& job, real* point, Gen& gen)
         : DoJobCP<Gen> (job, point, gen) {}
 
      void action (const Hypercube &h)
      {
         uniform (this->gen, h, this->point);
         this->job (this->point);
      }
   };

   /**
    *  AntitheticIntegrationCP
    */

   template<class Gen>
   class AntitheticDoJobCP : public DoJobCP<Gen>
   {
   public:
      AntitheticDoJobCP (Job& job, real* point, Gen& gen)
         : DoJobCP<Gen> (job, point, gen) {}
 
      void action (const Hypercube &h)
      {
         uniform (this->gen, h, this->point);
         this->job (this->point);
 
         for (unsigned i = 0; i < h.getDimension(); ++i)
         {
            this->point[i] = 2 * h.getCenter(i) - this->point[i];
         }
         this->job (this->point);
      }
   };

}  // namespace Private


/**
 *  stratifedDoJob()
 */

template<class Gen>
inline
void stratifiedDoJob (
   real* point, Gen &gen, const Hypercube &h, Job &job, Index n)
{
   Private::StratifiedDoJobCP<Gen> cp (job, point, gen);
   cp (h, n, cp.EQUIVOLUME);
}


/**
 *  antitheticDoJob()
 */

template<class Gen>
inline
void antitheticDoJob (
   real* point, Gen &gen, const Hypercube &h, Job &job, Index n)
{
   // Use CubePartitioner to perform action() an every sub-cube
 
   Private::AntitheticDoJobCP<Gen> cp (job, point, gen);
   cp (h, n / 2, cp.EQUIVOLUME);

   // if n is odd, draw one additional sample
 
   if (n & 1)
   {
      uniform (gen, h, point);
      job (point);
   }
}

}  // namespace HIntLib

#endif

