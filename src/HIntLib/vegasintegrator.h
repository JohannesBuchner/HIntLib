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
 *  VEGAS Integrator
 *
 *  An adaptive Monte Carlo integration routine
 *
 *  This implementation is based on
 *    [1] Lepage, L.P. 1978, Journal of Computational Physics, vol. 27,
 *        pp. 192 - 203.
 *    [2] Lepage, L.P. 1980, "VEGAS: An Adaptive Multidimensional Integration
 *        Program", Publication CLNS-80/447, Cornell Univeristy.
 *    [3] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P.
 *        Flannery. Numerical Recipes in C - The Art of Scientific Computing.
 *        second edition. Cambrdige University Press. Chapter 7.8.
 */

#ifndef VEGASINTEGRATOR_H
#define VEGASINTEGRATOR_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/integrator.h>

namespace HIntLib
{
   class PRNG;

   class VegasIntegrator : public Integrator
   {
   public:
      VegasIntegrator (PRNG* _mc) : mc(_mc) {}
      
      Status integrate (
         Function &f, const Hypercube &h, Index maxEval,
         real, real, EstErr &ee);
   private:
      PRNG* mc;
   };
}  // namespace HIntLib

#endif
