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

#ifndef HINTLIB_TEST_INTEGRAND_H
#define HINTLIB_TEST_INTEGRAND_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/integrand.h>

namespace HIntLib
{

class Hypercube;

/**
 *  TestIntegrand
 *
 *  Abstract class that adds the capability of returning the exact result for
 *  the integral over a specified domain.
 */
 
class TestIntegrand : public Integrand
{
public:
   TestIntegrand (int dim) : Integrand (dim) {}

   virtual real getExactResult (const Hypercube &) const = 0;
};


/**
 *  CStyleTestIntegrand
 *
 *  Creats a TestIntegrand based on a C-type function.
 *
 *  A pointer to the C-type function and the correct result has to be passed
 *  to the constructor.
 */

class CStyleTestIntegrand : public TestIntegrand
{
public:

   typedef real F (int, const real*);

   CStyleTestIntegrand (int dim, F *f, real exactResult)
      : TestIntegrand (dim), f(f), exactResult(exactResult) {}

   virtual real operator() (const real* p)  { return f (dim, p); }
   virtual real getExactResult (const Hypercube &) const { return exactResult; }

private:

   F * const f;
   const real exactResult;
};


/**
 *  EvalCounterIntegrand
 *
 *  Counts the number of calls to operator()
 */

class EvalCounterIntegrand : public TestIntegrand
{
public:

   EvalCounterIntegrand (TestIntegrand *f)
      : TestIntegrand (f->getDimension()), f (*f), counter (0) {}

   virtual real operator() (const real* p)  { ++counter; return f(p); }
   virtual real derivative (const real* p, int a)
      { ++counter; return f.derivative (p, a); }
   virtual real getExactResult (const Hypercube &h) const
      { return f.getExactResult (h); }

   Index getCounter() const  { return counter; }

private:

   TestIntegrand &f;
   Index counter;
};


/**
 *  DomainCheckerIntegrand
 *
 *  Checks if all calls to operator() are inside a given Hypercube
 */

class DomainCheckerIntegrand : public TestIntegrand
{
public:

   DomainCheckerIntegrand (TestIntegrand *, const Hypercube *);

   virtual real operator() (const real p []);
   virtual real derivative (const real* p, int);
   virtual real getExactResult (const Hypercube &hh) const
      { return f.getExactResult (hh); }

   bool isAllPointsInside (void)  { return allInside; }

private:

   TestIntegrand &f;
   const Hypercube &h;
   bool allInside;
};

}  // namespace HIntLib

#endif

