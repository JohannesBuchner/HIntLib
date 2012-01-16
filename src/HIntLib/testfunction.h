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
 *  TestFunction
 *
 *  (C) 2001 by Rudolf Schürer
 *
 *  Provides a number of specialized Functions that are used for the evaluation
 *  of integration routines.
 *
 *  TestFunction
 *
 *     Calcualates the exact value of the integral for a given hypercube.
 *     The allows to compare the approximation from the integration routine
 *     with the correct result.
 *
 *  The following classes are all subclasses of TestFunction:
 *
 *  SimpleFunction
 *
 *     Allows the cration of a TestFunction based to a C-style function
 *
 *  EvalCounterFunction
 *
 *     Counts the number of evaluations
 *
 *  DomainCheckerFunction
 *
 *     Checks if all calls to operator() use points inside a specified
 *     Hypercube
 *
 *  NodeTrackerFunction
 *
 *     Writes all points passed in operator() into a stream
 */ 

#ifndef TESTFUNCTION_H
#define TESTFUNCTION_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/function.h>

namespace HIntLib
{

class hypercube;

/**
 *  TestFunction
 *
 *  Abstract class that adds the capability of returning the exact result for
 *  the integral over a specified domain.
 */
 
class TestFunction : public Function
{
public:
   TestFunction (unsigned dim) : Function (dim) {}

   virtual real getExactResult (const Hypercube &) const = 0;
};


/**
 *  CStyleTestFunction
 *
 *  Creats a TestFunction based on a C-type function.
 *
 *  A pointer to the C-type function and the correct result has to be passed
 *  to the constructor.
 */

class CStyleTestFunction : public TestFunction
{
public:

   typedef real F (unsigned, const real*);

   CStyleTestFunction (unsigned dim, F *f, real exactResult)
      : TestFunction(dim), f(f), exactResult(exactResult) {}

   virtual real operator() (const real* p)  { return f (dim, p); }
   virtual real getExactResult (const Hypercube &) const { return exactResult; }

private:

   F * const f;
   const real exactResult;
};


/**
 *  EvalCounterFunction
 *
 *  Counts the number of calls to operator()
 */

class EvalCounterFunction : public TestFunction
{
public:

   EvalCounterFunction (TestFunction *f)
      : TestFunction (f->getDimension()), f (*f), counter (0) {}

   virtual real operator() (const real* p)  { ++counter; return f(p); }
   virtual real getExactResult (const Hypercube &h) const
      { return f.getExactResult (h); }

   Index getCounter() const  { return counter; }

private:

   TestFunction &f;
   Index counter;
};


/**
 *  DomainCheckerFunction
 *
 *  Checks if all calls to operator() are inside a given Hypercube
 */

class DomainCheckerFunction : public TestFunction
{
public:

   DomainCheckerFunction (TestFunction *, const Hypercube *);

   virtual real operator() (const real p []);
   virtual real getExactResult (const Hypercube &hh) const
      { return f.getExactResult (hh); }

   bool isAllPointsInside (void)  { return allInside; }

private:

   TestFunction &f;
   const Hypercube &h;
   bool allInside;
};

}  // namespace HIntLib

#endif

