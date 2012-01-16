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

#include <iomanip>

#include <HIntLib/vegas.h>

#include <HIntLib/mymath.h>
#include <HIntLib/kahanadd.h>
#include <HIntLib/hypercube.h>
#include <HIntLib/shiftscale.h>
#include <HIntLib/function.h>
#include <HIntLib/exception.h>
#include <HIntLib/pointset.h>

using std::min;

namespace L = HIntLib;
using namespace L;

namespace
{
   class VegasJob : public Job
   {
   public:
      VegasJob (const Hypercube &_h, Function &_f, unsigned _numSections);

      void operator() (const real* point);
      void reset (real volume);
      void updateDensityFunction (Index numPoints);
      
      real getSum () const  { return sum; }
      real getSumSquares () const  { return sumSquares; }

      real ALPHA;
      
   private:
      real& sectionUB (unsigned d, unsigned j)
         { return  _sectionUB [j * dim + d]; }
      real& sectionSqrSum (unsigned d, unsigned j)
         { return  _sectionSqrSum [j * dim + d]; }
      Index& sectionCount (unsigned d, unsigned j)
         { return  _sectionCount [j * dim + d]; }

      const unsigned numSections;
      const unsigned dim;
      ShiftScale ss;
      Function &f;
      KahanAdd sum, sumSquares;
      real avgWeight;

      Array<real> _sectionUB;
      Array<real> _sectionSqrSum;
      Array<Index>_sectionCount;
//XXX Array<real> _sectionInt;    (numSections * dim);
      Array<real> point;
      Array<real> sectionUBCopy;
      Array<real> r;
      Array<unsigned> currentSection;
   };


   /**
    *  Constructor
    */

   VegasJob::VegasJob (const Hypercube &h, Function &_f, unsigned _numSections)
      : ALPHA (1.5),
        numSections (_numSections),
        dim (h.getDimension()),
        ss (h),
        f (_f),
        // sum, sumSquares, and weight are initialized in reset()
        _sectionUB (numSections * dim),
        _sectionSqrSum (numSections * dim),
        _sectionCount (numSections * dim),
        point (dim),
        sectionUBCopy (numSections - 1),
        r (numSections),
        currentSection (dim)
   {
      for (unsigned d = 0; d < dim; ++d)
      {
         for (unsigned i = 0; i < numSections; ++i)
         {
            sectionUB(d,i) = real(i+1) / real (numSections);
         }
      }
   }


   /**
    *  reset()
    */

   void VegasJob::reset (real _avgWeight)
   {
      sum = 0.0;
      sumSquares = 0.0;
      avgWeight = _avgWeight;

      // Reset sectionInt and sectionSqrSum

      for (unsigned d = 0; d < dim; ++d)
      {
         for (unsigned j = 0; j < numSections; ++j)
         {
            // XXX sectionInt(d,j) = 
            sectionSqrSum(d,j) = 0.0;
            sectionCount (d,j) = 0;
         }
      }
   }

   
   /**
    *  operator()
    */

   void VegasJob::operator() (const real* uni)
   {
      real weight = avgWeight;

      for (unsigned d = 0; d < dim; ++d)
      {
         const unsigned j = currentSection [d] = unsigned (uni[d]);

         real diameter, coordinate;

         // transform x to coordinate

         if (j == 0)
         {
            diameter = sectionUB(d,0);
            coordinate = diameter * uni[d];
         }
         else
         {
            diameter = sectionUB(d,j) - sectionUB(d,j-1);
            coordinate = sectionUB(d,j-1) + diameter * (uni[d] - real(j));
         }

         weight *= diameter * real(numSections);

         point[d] = ss[d] (coordinate);    // scale to destination cube
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
         sectionSqrSum (d, currentSection [d]) += square;
         ++sectionCount(d, currentSection [d]);
      }
   }


   /**
    *   updateDensityFunction()
    */

   void VegasJob::updateDensityFunction(Index numPoints)
   {
      // For each dimension, update density function

      for (unsigned d = 0; d < dim; ++d)
      {
         // Scale  sectionSqrSum  according to the number of sample points
        
         for (unsigned j = 0; j < numSections; ++j)
         {
            sectionSqrSum (d,j) =
               sectionSqrSum (d,j) / real (sectionCount (d,j))
                                   * (numPoints / numSections);
         }

         // Smoothen the sectionSqrSum distribution
         //   (replace its value by the average of itself and its neighbors)
         // Also, sum up the results

         real sum;
         {
            real x2 = sectionSqrSum(d, 0);             // j = 0
            real x3 = sectionSqrSum(d, 1);

            sum = sectionSqrSum(d,0) = (x2 + x3) / 2.0;

            for (unsigned j = 1; j+1 < numSections; ++j)// j = 1..numSections-2
            {
               real x1 = x2;
               x2 = x3;
               x3 = sectionSqrSum(d, j+1);

               sum += sectionSqrSum(d,j) = (x1 + x2 + x3) / 3.0;
            }

            // j = numSections-1

            sum += sectionSqrSum(d,numSections-1) = (x2 + x3) / 2.0;
         }

         // Calculate correction factors for each section

         // const real minimum = sum / (numSections * 1000.0);
         const real minimum = 1e-30;
         real rSum = 0.0;

         for (unsigned j = 0; j < numSections; ++j)
         {
            const real x = std::max (minimum, sectionSqrSum (d,j));

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
      cout << "sectionSqrSum:" << endl;
      for (unsigned d = 0; d < dim; ++d)
      {
         cout << "Dim " << d << ": ";
         for (unsigned j = 0; j < numSections; ++j)
         {
            cout << setw(18) << sectionSqrSum(d,j);
         }
         cout << endl;
      }
      cout << "sectionUB:" << endl;
      for (unsigned d = 0; d < dim; ++d)
      {
         cout << "Dim " << d << ": ";
         for (unsigned j = 0; j < numSections; ++j)
         {
            cout << setw(18) << sectionUB(d,j);
         }
         cout << endl;
      }
#endif
   }
}


/**
 *  integrate()
 *
 *  Main integration routine of the VEGAS algorithm
 */

L::Vegas::Status L::Vegas::integrate (
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
      #ifdef HINTLIB_NO_EXCEPTIONS
         ee.set (0.0, 0.0);
         return ERROR;
      #else
         throw MaxEvaluationsRequired ();
      #endif
   }

   // Constants

   static const unsigned MAX_ITER     = 5;
   static const unsigned MIN_SECTIONS = 5;
   static const unsigned MAX_SECTIONS = 50;

   // Calculate an appropriate number of iterations, sections and point budgets

   // Overflow precautions if Index is > 32 bits! XXX 

   const unsigned numPointsPreTotal
       = unsigned (min (maxEval / 4, Index(MAX_ITER) * MAX_SECTIONS * 500));

   const unsigned numIter = min (MAX_ITER, numPointsPreTotal / 800);
   const unsigned numPointsPre = numIter ? numPointsPreTotal / numIter : 0;
   const unsigned numPointsFinal = maxEval - numPointsPre * numIter;

   const unsigned numSections = std::max(MIN_SECTIONS, min(MAX_SECTIONS,
                                    unsigned(sqrt(real(numPointsPre)))));

#if 0
cerr << "\n  Sections: " << numSections << "  Iterations: " << numIter << endl
     << "  Presampling Points: " << numPointsPreTotal
     << " (" << real (numPointsPreTotal) / real (maxEval) * 100.0 << " %)"
     << "  Final Points: " << numPointsFinal << endl
     << "  Points Pre Total: " << numPointsPreTotal
     << endl;
#endif

   // Other initialization

   real integral = 0.0;
   real stddev   = 0.0;
   real chi2a    =-1.0;

   real sumIntegral      = 0.0;
   real sumChi           = 0.0;
   real sumVarianceRecip = 0.0;

   Array<real> point (h.getDimension());
   Hypercube sampleSpace (h.getDimension(), 0.0, real(numSections));
   ps->setCube (&sampleSpace);
   VegasJob vegasJob (h, f, numSections);
   vegasJob.ALPHA = ALPHA;

   // **** Iteration Loop ****************

   for (unsigned iter = 0; ; ++iter)
   {
      // for the final run, use more points

      const unsigned numPoints =
         (iter == numIter) ? numPointsFinal : numPointsPre;

      // Sample

      vegasJob.reset (h.getVolume() / real(numPoints));
      ps->doJob (point, vegasJob, numPoints);

      // Update global integral

      real sum = vegasJob.getSum();

      const real dv2g = sqr(numPoints * h.getVolume())
                   / sqr(real(numPoints)) / (numPoints-1.0);

      real variance = sqrt(vegasJob.getSumSquares() * numPoints);
           variance = (variance - sum) * (variance + sum) * dv2g;

      if (combineResults && variance > 0.0)
      {
         const real varianceRecip = 1.0 / variance;

         sumIntegral      += varianceRecip * sum;
         sumChi           += varianceRecip * sqr(sum);
         sumVarianceRecip += varianceRecip;

         integral = sumIntegral / sumVarianceRecip;
         stddev   = sqrt (1.0 / sumVarianceRecip);
      }
      else
      {
         integral = sum;
         stddev   = sqrt (variance);
      }

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

      vegasJob.updateDensityFunction (numPoints);

   }  // all iterations perfromed

   // return the result

   ee.set (integral, stddev);

   Status status = checkRequestedError (ee, reqAbsError, reqRelError);
   return (status == ERROR) ? MAX_EVAL_REACHED : status;
}

