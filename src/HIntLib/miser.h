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
 *  MISER Integrator
 *
 *  This implementation is based on:
 *
 *    [1] W.H. Press, G.R. Farrar. Recursive Stratified Sampling for
 *        Multidimensional Monte Carlo Integration. Computers in Physics,
 *        vol. 4(2), pp. 190 - 195, 1990.
 *    [2] W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery.
 *        Numerical Recipes in C - The Art of Scientific Computing.
 *        second edition. Cambrdige University Press. Chapter 7.8.
 */

#ifndef HINTLIB_MISER_H
#define HINTLIB_MISER_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/integrator.h>

namespace HIntLib
{
   class PRNG;
   class PointSet;

   class Miser : public Integrator
   {
   public:
      Miser (PointSet* presample, PointSet* sample);
      Miser (PointSet*);
      void defaults();

      Status integrate (
         Function &f, const Hypercube &h, Index maxEval,
         real reqAbsError, real reqRelError,
         EstErr &ee);

   private:
      PointSet* presamplePointSet;
      PointSet* samplePointSet;

   public:
      Index MIN_POINTS;
      Index FREE_POINTS;
      Index MIN_PRESAMPLING_POINTS;
      Index MAX_PRESAMPLING_POINTS;
      real PRESAMPLING_RATE;
   };
}  // namespace HIntLib

#endif
