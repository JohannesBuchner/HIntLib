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


#ifndef HINTLIB_MC_POINTSET_H
#define HINTLIB_MC_POINTSET_H 1
 
#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/pointset.h>
#include <HIntLib/distribution.h>
#include <HIntLib/mcintegration.h>

namespace HIntLib
{

/*****************************************************************************/
/******** PRNG  Implementation  **********************************************/
/*****************************************************************************/

template<class T>
class PRNGImp : public PRNG
{
protected:
   T mc;

   // Constructor

public:
   PRNGImp(unsigned start) : mc(start) {}
   PRNGImp() {}

   void init (unsigned start)       { mc.init (start); }

   // Forward all calls

   real uniform ()                  { return HIntLib::uniform (mc); }
   real uniform (real ub)           { return HIntLib::uniform (mc, ub); }
   real uniform (real lb, real ub)  { return HIntLib::uniform (mc, lb, ub); }

   int  equidist (int ub)           { return HIntLib::equidist (mc, ub); }
   int  equidist (int min, int ub)  { return HIntLib::equidist (mc, min, ub); }

   void uniform (Hypercube &h, real* p)  { HIntLib::uniform (mc, h, p); }

   size_t getStateSize () const      { return mc.getStateSize(); }
   void saveState (void *p) const    { mc.saveState (p); }
   void restoreState (const void *p) { mc.restoreState (p); }

   real getRange() const  { return mc.getRange(); }
   real operator()()      { return real (mc()); }
};


/*****************************************************************************/
/**********  General Base  ***************************************************/
/*****************************************************************************/

class MCPointSetBase_ : public MultiPointSet
{
public:
   bool doJobRep (real *, ReportingJob &, Index n);
   Index getOptimalNumber (Index n, const Hypercube &);
   void setCube (const Hypercube *);

protected:
   MCPointSetBase_ (Index _alignment = 1) : h (0), alignment (_alignment) {}

   const Hypercube* h;

private:
   const Index alignment;
};


/*****************************************************************************/
/****  Monte Carlo and related Point Sets    *********************************/
/*****************************************************************************/

template<class T>
class MCPointSetBase : public MCPointSetBase_,
                       public PRNGImp<T>
{
protected:
   MCPointSetBase(Index alignment, unsigned start)
      : MCPointSetBase_(alignment), PRNGImp<T> (start) {}
   MCPointSetBase(Index alignment) : MCPointSetBase_(alignment) {}

public:
   virtual void select (unsigned start, unsigned num) { this->mc.init (start); }
   virtual void randomize (unsigned seed) { this->mc.init (seed); }
};


/**
 *  Monte Carlo Point Set
 */

template<class T>
class MonteCarloPointSetBase : public MCPointSetBase<T>
{
public:
   MonteCarloPointSetBase(unsigned start) : MCPointSetBase<T> (1, start) {}
   MonteCarloPointSetBase()               : MCPointSetBase<T> (1) {}


   void doJob (real *point, Job &job, Index n)
   {
      mcDoJob (point, this->mc, *(this->h), job, n);
   }

   bool doJobRep (real *point, ReportingJob &job, Index n)
   {
      return mcDoJobRep (point, this->mc, *(this->h), job, n);
   }
};

template<class T, class Sum = real>
class MonteCarloPointSet : public MonteCarloPointSetBase<T>
{
public:
   MonteCarloPointSet (unsigned start) : MonteCarloPointSetBase<T>(start) {}
   MonteCarloPointSet () {}

   void integrate (real point [], Integrand &f, Index n, PointSet::Stat &stat)
   {
      Statistic<real, Sum> s;
      mcIntegration (point, this->mc, *(this->h), f, n, s);
      stat = s;
   }

   void integrate (real point [], Integrand &f, Index n, PointSet::StatVar &stat)
   {
      StatisticVar<real, Sum> s;
      mcIntegration (point, this->mc, *(this->h), f, n, s);
      stat = s;
   }
};


/**
 *  Stratified Point Set
 */

template<class T>
class StratifiedPointSetBase : public MCPointSetBase<T>
{
public:
   StratifiedPointSetBase(unsigned start) : MCPointSetBase<T> (1, start) {}
   StratifiedPointSetBase()               : MCPointSetBase<T> (1) {}


   void doJob (real *point, Job &job, Index n)
   {
      stratifiedDoJob (point, this->mc, *(this->h), job, n);
   }
};

template<class T, class Sum = real>
class StratifiedPointSet : public StratifiedPointSetBase<T>
{
public:
   StratifiedPointSet (unsigned start) : StratifiedPointSetBase<T>(start) {}
   StratifiedPointSet () {}

   void integrate (real point [], Integrand &f, Index n, PointSet::Stat &stat)
   {
      Statistic<real, Sum> s;
      stratifiedIntegration (point, this->mc, *(this->h), f, n, s);
      stat = s;
   }

   void integrate (real point [], Integrand &f, Index n, PointSet::StatVar &stat)
   {
      StatisticVar<real, Sum> s;
      stratifiedIntegration (point, this->mc, *(this->h), f, n, s);
      stat = s;
   }
};


/**
 *  Antithetic Point Set
 */

template<class T>
class AntitheticPointSetBase : public MCPointSetBase<T>
{
public:
   AntitheticPointSetBase(unsigned start) : MCPointSetBase<T> (2, start) {}
   AntitheticPointSetBase()               : MCPointSetBase<T> (2)        {}

   void doJob (real *point, Job &job, Index n)
   {
      antitheticDoJob (point, this->mc, *(this->h), job, n);
   }
};

template<class T, class Sum = real>
class AntitheticPointSet : public AntitheticPointSetBase<T>
{
public:
   AntitheticPointSet (unsigned start) : AntitheticPointSetBase<T>(start) {}
   AntitheticPointSet () {}

   void integrate (
      real point [], Integrand &f, Index n, PointSet::Stat &stat)
   {
      Statistic<real, Sum> s;
      antitheticIntegration (point, this->mc, *(this->h), f, n, s);
      stat = s;
   }

   void integrate (
      real point [], Integrand &f, Index n, PointSet::StatVar &stat)
   {
      StatisticVar<real, Sum> s;
      antitheticIntegration (point, this->mc, *(this->h), f, n, s);
      stat = s;
   }
};


/*****************************************************************************/
/****   Create    Point Sets     *********************************************/
/*****************************************************************************/


/**
 *  MonteCarlo PointSet Create
 */

template<class T>
class MonteCarloPointSetCreateBase : public MCPointSetBase_
{
protected:
   MonteCarloPointSetCreateBase (unsigned s) : start (s) {}

   unsigned start;

public:
   virtual void select (unsigned s, unsigned num) { start = s; }
   virtual void randomize (unsigned seed) { start = seed; }

   void doJob (real *point, Job &job, Index n)
   {
      T mc (start);
      mcDoJob (point, mc, *(this->h), job, n);
   }

   bool doJobRep (real *point, ReportingJob &job, Index n)
   {
      T mc (start);
      return mcDoJobRep (point, mc, *(this->h), job, n);
   }
};

template<class T, class Sum = real>
class MonteCarloPointSetCreate : public MonteCarloPointSetCreateBase<T>
{
public:

   MonteCarloPointSetCreate (unsigned s = 0)
      : MonteCarloPointSetCreateBase<T> (s) {}

   void integrate (
      real point [], Integrand &f, Index n, PointSet::StatVar &stat)
   {
      T m (this->start);
      StatisticVar<real, Sum> s;
      mcIntegration (point, m, *(this->h), f, n, s);
      stat = s;
   }

   void integrate (
      real point [], Integrand &f, Index n, PointSet::Stat &stat)
   {
      T m (this->start);
      Statistic<real, Sum> s;
      mcIntegration (point, m, *(this->h), f, n, s);
      stat = s;
   }
};

/**
 *  Stratified PointSet Create
 */

template<class T>
class StratifiedPointSetCreateBase : public MCPointSetBase_
{
protected:
   StratifiedPointSetCreateBase (unsigned s) : start (s) {}

   unsigned start;

public:
   virtual void select (unsigned s, unsigned num) { start = s; }
   virtual void randomize (unsigned seed) { start = seed; }

   void doJob (real *point, Job &job, Index n)
   {
      T mc (start);
      stratifiedDoJob (point, mc, *(this->h), job, n);
   }
};

template<class T, class Sum = real>
class StratifiedPointSetCreate : public StratifiedPointSetCreateBase<T>
{
public:

   StratifiedPointSetCreate (unsigned s = 0)
      : StratifiedPointSetCreateBase<T> (s) {}

   void integrate (
      real point [], Integrand &f, Index n, PointSet::StatVar &stat)
   {
      T m (this->start);
      StatisticVar<real, Sum> s;
      stratifiedIntegration (point, m, *(this->h), f, n, s);
      stat = s;
   }

   void integrate (
      real point [], Integrand &f, Index n, PointSet::Stat &stat)
   {
      T m (this->start);
      Statistic<real, Sum> s;
      stratifiedIntegration (point, m, *(this->h), f, n, s);
      stat = s;
   }
};


/**
 *  Antithetic PointSet Create
 */

template<class T>
class AntitheticPointSetCreateBase : public MCPointSetBase_
{
protected:
   AntitheticPointSetCreateBase (unsigned s) : MCPointSetBase_(2), start (s) {}

   unsigned start;

public:
   virtual void select (unsigned s, unsigned num) { start = s; }
   virtual void randomize (unsigned seed) { start = seed; }

   void doJob (real *point, Job &job, Index n)
   {
      T mc (start);
      antitheticDoJob (point, mc, *(this->h), job, n);
   }
};

template<class T, class Sum = real>
class AntitheticPointSetCreate : public AntitheticPointSetCreateBase<T>
{
public:

   AntitheticPointSetCreate (unsigned s = 0)
      : AntitheticPointSetCreateBase<T> (s) {}

   void integrate (
      real point [], Integrand &f, Index n, PointSet::StatVar &stat)
   {
      T m (this->start);
      StatisticVar<real, Sum> s;
      antitheticIntegration (point, m, *(this->h), f, n, s);
      stat = s;
   }

   void integrate (
      real point [], Integrand &f, Index n, PointSet::Stat &stat)
   {
      T m (this->start);
      Statistic<real, Sum> s;
      antitheticIntegration (point, m, *(this->h), f, n, s);
      stat = s;
   }
};
}  // namespace HIntLib

#endif

