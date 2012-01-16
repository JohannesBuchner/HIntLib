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


#ifndef HINTLIB_QMCROUTINESIMP_H
#define HINTLIB_QMCROUTINESIMP_H 1
 
#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/pointset.h>
#include <HIntLib/qmcintegration.h>

// implementing *.cpp file: halton.cpp

namespace HIntLib
{

template<class Gen>
class QMCPointSetBase : public PartitionablePointSet
{
protected:
   const Index skip;
   const Hypercube* h;

public:

   QMCPointSetBase (Index s) : skip (s), h (0) {}

   virtual void setCube (const Hypercube *_h)  { h = _h; }
   virtual void randomize (unsigned)  {}

   Index getOptimalNumber (Index max, const Hypercube &h)
   {
      return Gen::getOptimalNumber (max + skip, h) - skip;
   }

   void doJobPartition
     (real *point, Job &job, Index num, Index begin, Index end)
   {
      Gen gen (*h);
      qmcDoJob (point, gen, job, begin, end);
   }

   bool doJobRep (real *point, ReportingJob &job, Index num)
   {
      Gen gen (*h);
      return qmcDoJob (point, gen, job, 0, num);
   }
};

template<class Gen, class Sum>
class QMCPointSet : public QMCPointSetBase<Gen>
{
public:
   QMCPointSet (Index s = 0) : QMCPointSetBase<Gen> (s) {}

   void integratePartition
      (real* point, Function &f, Index num, Index begin, Index end,
       PointSet::Stat& stat)
   {
      Gen gen (*h);
      Statistic<real, Sum, Index> s;
      qmcIntegration (point, gen, f, begin + skip, end + skip, s);
      stat = s;
   }

   void integratePartition
      (real* point, Function &f, Index num, Index begin, Index end,
       PointSet::StatVar& stat)
   {
      Gen gen (*h);
      StatisticVar<real, Sum, Index> s;
      qmcIntegration (point, gen, f, begin + skip, end + skip, s);
      stat = s;
   }
};

// Specialization for sum=real

template<class Gen>
class QMCPointSet<Gen,real> : public QMCPointSetBase<Gen>
{
public:
   QMCPointSet (Index s = 0) : QMCPointSetBase<Gen> (s) {}

   void integratePartition
      (real *point, Function &f, Index num, Index begin, Index end,
       PointSet::Stat& stat)
   {
      Gen gen (*h);
      qmcIntegration (point, gen, f, begin + skip, end + skip, stat);
   }

   void integratePartition
      (real *point, Function &f, Index num, Index begin, Index end,
       PointSet::StatVar& stat)
   {
      Gen gen (*h);
      qmcIntegration (point, gen, f, begin + skip, end + skip, stat);
   }
};

}  // namespace HIntLib

#endif

