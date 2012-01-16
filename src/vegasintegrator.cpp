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

#ifdef __GNUG__
#pragma implementation
#endif

#include <algorithm>

#include <HIntLib/vegasintegrator.h>

#include <HIntLib/mymath.h>
#include <HIntLib/kahanadd.h>
#include <HIntLib/hypercube.h>
#include <HIntLib/function.h>
#include <HIntLib/exception.h>
#include <HIntLib/pointset.h>

using std::min;

namespace L = HIntLib;


L::VegasIntegrator::Status L::VegasIntegrator::integrate (
   Function &f,
   const Hypercube &h,
   Index maxEval,
   real reqRelError, real reqAbsError,
   EstErr &ee)
{
   checkDimension(h, f);

   // We need an upper bound to the number of sampling points

   if (maxEval == 0)
   {
      #ifdef CRAY
         ee.set (0.0, 0.0);
         return ERROR;
      #else
         throw MaxEvaluationsRequired ();
      #endif
   }

   // Constants

   static const unsigned MAX_ITER     = 6;
   static const unsigned MAX_SECTIONS = 100;

   // Calculate an appropriate number of iterations, sections and point budgets

   // Overflow precautions if Index is > 32 bits! XXX 

   const unsigned numPointsPreTotal
       = unsigned (min (maxEval / 4, Index(MAX_ITER) * MAX_SECTIONS * 200));

   const unsigned numIter = min (MAX_ITER, numPointsPreTotal / 800);
   const unsigned numPointsPre = numIter ? numPointsPreTotal / numIter : 0;
   const unsigned numPointsFinal = maxEval - numPointsPre * numIter;

   const unsigned numSections = std::max(3u, min(100u,
                                    unsigned(sqrt(real(numPointsPre)))));

#if 0
cerr << "  Sections: " << numSections << "  Iterations: " << numIter << endl
     << "  Presampling Points: " << numPointsPre
     << "  Final Points: " << numPointsFinal << endl
     << "  Points Pre Total: " << numPointsPreTotal
     << endl;
#endif

   // Other initialization

   const real ALPHA = 1.5;

   const unsigned dim = h.getDimension();

   const real volume = h.getVolume();

   Array<real> _sectionUB     (numSections * dim);
   Array<real> sectionUBCopy  (numSections - 1);
   // XXX Array<real> _sectionInt    (numSections * dim);
   Array<real> _sectionSqrSum (numSections * dim);

   Array<real> r (numSections);
   Array<unsigned> currentSection (dim);

   Array<real> point (dim);    // The point f is evaluated at

   #define sectionUB(d,s)     (_sectionUB    [(d) * numSections + (s)])
   /// XXX #define sectionInt(d,s)    (_sectionInt   [(d) * numSections + (s)])
   #define sectionSqrSum(d,s) (_sectionSqrSum[(d) * numSections + (s)])

   // Initialize sectionLB.  In the beginning sections are equidistributed.

   for (unsigned d = 0; d < dim; ++d)
   {
      for (unsigned i = 0; i < numSections; ++i)
      {
         sectionUB(d,i) = real(i+1) / real (numSections);
      }
   }

   real integral = 0.0;
   real stddev   = 0.0;
   real chi2a    =-1.0;

   real sumIntegral      = 0.0;
   real sumChi           = 0.0;
   real sumVarianceRecip = 0.0;


   // **** Iteration Loop ****************

   unsigned numPoints = numPointsPre;

   for (unsigned iter = 0; ; ++iter)
   {
      // for the final run, use more points

      if (iter == numIter)  numPoints = numPointsFinal;

      // these values depend on the number of points

      const real pointVolume = volume / real(numPoints); 

      const real dv2g = sqr(numPoints*volume)
                   / sqr(real(numPoints)) / (numPoints-1.0);

      // Reset sectionInt and sectionSqrSum

      for (unsigned d = 0; d < dim; ++d)
      {
         for (unsigned j = 0; j < numSections; ++j)
         {
            // XXX sectionInt(d,j) = 
            sectionSqrSum(d,j) = 0.0;
         }
      }

      KahanAdd sum        = 0.0;
      KahanAdd sumSquares = 0.0;

      // ************* Sample Loop

      for (unsigned i = 0; i < numPoints; ++i)
      {

         // Set up point

         real weight = pointVolume;

         for (unsigned d = 0; d < dim; ++d)
         {
            const real random = mc->uniform (real(numSections));
            const unsigned j = currentSection [d] = unsigned (random);

            real diameter;
            real coordinate;

            // Set x to the correct value in [0,1]

            if (j == 0)
            {
               diameter = sectionUB(d,0);
               coordinate = diameter * random;
            }
            else
            {
               diameter = sectionUB(d,j) - sectionUB(d,j-1);
               coordinate = sectionUB(d,j-1) + diameter * (random - real(j));
            }

            weight *= diameter * real(numSections);

            // Scale x to the Hypercube h

            point[d] = h.getLowerBound(d) + coordinate * h.getDiameter(d);
         }

         // Evalute f.  Update integral

         const real value  = weight * f(point);
         const real square = sqr(value);

         sum        += value;
         sumSquares += square;


         // Update integral for each section

         for (unsigned d = 0; d < dim; ++d)
         {
            // XXX sectionInt    (d,currentSection[d]) += value;
            sectionSqrSum (d,currentSection[d]) += square;
         }
      }

      // Update global integral

      real variance = sqrt(sumSquares * numPoints);
           variance = (variance - sum) * (variance + sum);

           variance *= dv2g;

      const real varianceRecip = 1.0 / variance;

      sumIntegral      += varianceRecip * sum;
      sumChi           += varianceRecip * sqr(sum);
      sumVarianceRecip += varianceRecip;

      integral = sumIntegral / sumVarianceRecip;
      stddev   = sqrt (1.0 / sumVarianceRecip);

      // Do chi^2 test to check consistency

      if (iter >= 1)
      {
         chi2a = (sumChi - sumIntegral * integral) / real(iter);
         if (chi2a < 0.0) chi2a = 0.0;
      }

#if 0
      cout << "Iteration " << iter << ": " << sum << " +/- " << sqrt(variance) << endl
           << "      Total: " << integral << " +/- " << stddev
           << "  chi = " << chi2a << endl;
#endif

      // If this was the last iteration, we can stop here

      if (iter == numIter)  break;

      // For each dimension, update density function

      for (unsigned d = 0; d < dim; ++d)
      {
         real sum;

         // Smoothen the sectionSqrSum distribution
         //   (replace its value by the average of itself and its neighbors)
         // Also, sum up the results

         {
            real x2 = sectionSqrSum(d, 0);             // j = 0
            real x3 = sectionSqrSum(d, 1);

            sum = sectionSqrSum(d,0) = (x2 + x3) / 2.0;

            for (unsigned j = 1; j+1 < numSections; ++j)// j = 1..numSections-2
            {
               real x12 = x2 + x3;
               x2 = x3;
               x3 = sectionSqrSum(d, j+1);

               sum += sectionSqrSum(d,j) = (x12 + x3) / 3.0;
            }

            // j = numSections-1

            sum += sectionSqrSum(d,numSections-1) = (x2 + x3) / 2.0;
         }

         // Calculate correction factors for each section

         real rSum = 0.0;

         for (unsigned j = 0; j < numSections; ++j)
         {
            const real x = std::max(real(1e-30), sectionSqrSum(d,j));

            rSum += r[j] = pow ((1.0 - x/sum) / (log(sum) - log(x)), ALPHA); 
         } 

         // Apply these corrections to the section upper bounds

         {
            const real rAvg = rSum / real(numSections);
            real x = 0.0;

            for (unsigned j = 0, rIndex = 0; j+1 < numSections; ++j)
            {
               while (x < rAvg)  x += r[rIndex++];

               x -= rAvg;

               const real lb = (rIndex > 1) ? sectionUB(d,rIndex-2) : 0.0;
               const real ub =                sectionUB(d,rIndex-1);

               sectionUBCopy[j] = ub - (ub-lb) * x / r[rIndex-1];
            } 
         }

         for (unsigned j = 0; j+1 < numSections; ++j)
         {
            sectionUB(d,j) = sectionUBCopy[j];
         }

      }  // for(d=0..dim-1)

#if 0
      for (unsigned d = 0; d < dim; ++d)
      {
         for (unsigned j = 0; j < numSections; ++j)
         {
            cout << setw(10) << sectionSqrSum(d,j);
         }
         cout << endl;
      }
      for (unsigned d = 0; d < dim; ++d)
      {
         for (unsigned j = 0; j < numSections; ++j)
         {
            cout << setw(10) << sectionUB(d,j);
         }
         cout << endl;
      }
#endif
   }

   // return the result

   ee.set (integral, stddev);

   Status status = checkRequestedError (ee, reqAbsError, reqRelError);
   return (status == ERROR) ? MAX_EVAL_REACHED : status;
}

