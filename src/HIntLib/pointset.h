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


#ifndef POINT_SET_H
#define POINT_SET_H 1
 
#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/defaults.h>
#include <HIntLib/statistic.h>


namespace HIntLib
{
   class Hypercube;
   class Function;

   /**
    *  Job
    *
    *  Something that can be performed for every point in a Point Set
    */

   class Job
   {
   public:
      virtual void operator() (const real*) = 0;
      virtual ~Job() {}
   };


   /**
    *  Reporting Job
    *
    *  Something that can be performed for points in a Point Set
    *
    *  It returns false if it should be applied to further points,
    *  true if the algorithm can terminate.
    */

   class ReportingJob
   {
   public:
      virtual bool operator() (const real*) = 0;
      virtual ~ReportingJob() {}
   };


   /**
    *  Point Set
    *
    *  Abstract concept of a set/sequence of points.
    *
    *  Allows Function evaluation and Job application to all points.
    */

   class PointSet
   {
   public:
      typedef Statistic<real,real,Index> Stat;
      typedef StatisticVar<real,real,Index> StatVar;

      virtual ~PointSet() {}

      virtual Index getOptimalNumber (Index, const Hypercube &) = 0;
      virtual void setCube (const Hypercube *) = 0;

      virtual void integrate (real *, Function &, Index, Stat&) = 0;
      virtual void integrate (real *, Function &, Index, StatVar&) = 0;

      virtual void doJob    (real *,          Job &, Index) = 0;
      virtual bool doJobRep (real *, ReportingJob &, Index) = 0;
   };


   /**
    *  Partitionable Point Set
    *
    *  A Point Set that can be partitioned, i.e. arbitrary blocks can be
    *  extracted.
    */

   class PartitionablePointSet : virtual public PointSet
   {
   public:
      void integrate (real *, Function &, Index, Stat&);
      void integrate (real *, Function &, Index, StatVar&);

      virtual void integratePartition (
         real *, Function &, Index, Index, Index, Stat&) = 0;
      virtual void integratePartition (
         real *, Function &, Index, Index, Index, StatVar&) = 0;

      void doJob (real *, Job &, Index);

      virtual void doJobPartition (real *, Job &, Index, Index, Index) = 0;
   };


   /**
    *  Multi Point Set
    *
    *  A set of Point Sets.  Different sets can be selected with select().
    */

   class MultiPointSet : virtual public PointSet
   {
   public:
      virtual void select (unsigned i, unsigned num) = 0;
   };


   /**
    *  PRNG
    *
    *  Abstract base class for Pseudo Random Number Generators
    */

   class PRNG
   {
   public:
      virtual ~PRNG () {}
      
      // Initialize Generator
 
      virtual void init (unsigned start) = 0;

      // Return a random number

      virtual real uniform() = 0;        // (0,1)
      virtual real uniform(real) = 0;
      virtual real uniform(real, real) = 0;
 
      virtual int  equidist (int ub) = 0;           // 0,1,2,...,ub-1
      virtual int  equidist (int min, int ub) = 0;  // min,...,ub-1

      virtual void uniform (Hypercube &h, real* p) = 0;

      // Save and restore state
 
      virtual size_t getStateSize() const = 0;
      virtual void saveState (void *) const = 0;
      virtual void restoreState (const void *) = 0;

      // low level

      virtual real getRange() const = 0;
      virtual real operator()() = 0;
   };

}  // namespace HIntLib

#endif

