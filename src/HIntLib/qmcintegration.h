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
 *  Monte Carlo integration
 *
 *  Templates for executing the inner loop of Monte Carlo integration routines
 */

#ifndef QMCINTEGRATION_H
#define QMCINTEGRATION_H 1

#include <HIntLib/hypercube.h>
#include <HIntLib/statistic.h>
#include <HIntLib/function.h>
#include <HIntLib/pointset.h>


namespace HIntLib
{

   template<class Gen, class Stat>
   inline
   void qmcIntegration (
      real point [], Gen &gen, Function &f,
      Index begin, Index end, Stat &stat)
   {
      if (begin < end)
      {
         gen.first (point, begin);

         do
         {
            stat << f(point);
            gen.next(point);
         }
         while (gen.getIndex() != end);
      }
   }

   template<class Gen>
   inline
   void qmcDoJob (real point [], Gen &gen, Job &job, Index begin, Index end)
   {
      if (begin < end)
      {
         gen.first (point, begin);

         do
         {
            job (point);
            gen.next(point);
         }
         while (gen.getIndex() != end);
      }
   }

   template<class Gen>
   inline
   bool qmcDoJob (
     real point [], Gen &gen, ReportingJob &job, Index begin, Index end)
   {
      if (begin < end)
      {
         gen.first (point, begin);

         do
         {
            if (job (point))  return true;
            gen.next(point);
         }
         while (gen.getIndex() != end);
      }

      return false;
   }

}  // namespace HIntLib

#endif

