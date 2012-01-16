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
 *  integrator.h
 *
 *  Abstract base class for objects that can integrate a given Integrand over
 *  a certain Hypercube
 */

#ifndef HINTLIB_INTEGRATOR_H
#define HINTLIB_INTEGRATOR_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/esterr.h>


namespace HIntLib
{

class Integrand;
class Hypercube;


class Integrator
{
public:

   /*
    * Return value of integrate
    *
    * If more than one value is possilbe, the first one from the list is
    * returned
    */

   enum Status {
      ABS_ERROR_REACHED,  // Estimated error <= reqAbsError
      REL_ERROR_REACHED,  // Estimated error <= reqRelError
      MAX_EVAL_REACHED,   // #integ evaluations reached
      ERROR               // Non of these, but some result was returned
   #ifdef HINTLIB_PARALLEL
     ,WRONG_NODE          // In parallel mode, only node 0 returns a result
   #endif
   };

   // A virtual destructor is required by many child classes

   virtual ~Integrator () {}

   // calcualte integral

   real operator() (Integrand &f, const Hypercube &h,
                    Index maxEval,
                    real reqAbsError, real reqRelError)
   {
      EstErr ee;

      integrate (f, h, maxEval, reqAbsError, reqRelError, ee);

      return ee.getEstimate ();
   }

   // Calcualate integral and error

   virtual
   Status integrate (
      Integrand &, const Hypercube &, Index maxEval,
      real reqAbsError, real reqRelError, EstErr &ee) = 0;

   void randomize (unsigned _seed)  { seed = _seed; }

protected:

   // Normal constructor.  No copy constructor

   Integrator () : seed (0) {}

   // A number of convenience functions used by many child classes

   static
   Status checkRequestedError (const EstErr &, real reqAbsErr, real reqRelErr);

   static void checkDimension (const Hypercube &, const Integrand &);
   static void checkTerminationCriteria (Index, real, real, bool);

   unsigned getSeed ()  { return seed; }

private:

   Integrator (const Integrator&);             // Do not copy!
   Integrator& operator= (const Integrator&);  // Do not assign!

   unsigned seed;
};


inline
Integrator::Status Integrator::checkRequestedError (
   const EstErr &ee, real reqAbsError, real reqRelError)
{
   if (ee.getError ()    <= reqAbsError) return ABS_ERROR_REACHED;
   if (ee.getRelError () <= reqRelError) return REL_ERROR_REACHED;
 
   return ERROR;
}

}  // namespace HIntLib

#endif

