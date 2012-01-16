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
 *  AdaptIntegrator
 *
 */

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/integrator.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/exception.h>
#include <HIntLib/hypercube.h>
#include <HIntLib/integrand.h>


namespace L = HIntLib;

void
L::Integrator::checkDimension (const Hypercube &h, const Integrand &f)
{
   L::checkDimensionEqual (h.getDimension(), f.getDimension());
}

void
L::Integrator::checkTerminationCriteria (
      Index maxEval, real reqAbsError, real reqRelError, bool maxEvalRequired)
{
   if (reqAbsError < 0)  throw RequestedErrorNegative (reqAbsError);

   if (reqRelError < 0)  throw RequestedErrorNegative (reqRelError);

   if (maxEvalRequired && maxEval == 0)  throw MaxEvaluationsRequired();

   if (reqAbsError == 0 && reqRelError == 0 && maxEval == 0)
   {
      throw TerminationCriterionMissing();
   }
}


