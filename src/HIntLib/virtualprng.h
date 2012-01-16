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


#ifndef HINTLIB_PRNG_H
#define HINTLIB_PRNG_H 1
 
#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_CSTDDEF
  #include <cstddef>
  #define HINTLIB_SDN ::std::
#else
  #include <stddef.h>
  #define HINTLIB_SDN ::
#endif

namespace HIntLib
{

/**
 *  VirtualPRNG
 *
 *  Abstract base class for Pseudo Random Number Generators
 */

class VirtualPRNG
{
public:
   // Initialize Generator
 
   virtual void init (unsigned start) = 0;

   // Information about the generator
 
   virtual u32  getMax() const = 0;
   virtual real getResolution() const = 0;
 
   // Return a random number
 
   virtual u32 operator() () = 0;     // {0,...,getMax()}
   virtual int operator() (int) = 0;  // {0,...,max-1}
   virtual real getReal() = 0;        // (0,1)
 
   // Save and restore state
 
   virtual HINTLIB_SDN size_t getStateSize() const = 0;
   virtual void saveState (void *) const = 0;
   virtual void restoreState (const void *) = 0;
};


/**
 *  VirtualPRNGNew<T>
 *
 *  A T, implementing the PRNG interface
 */

template<class T>
class VirtualPRNGNew : public VirtualPRNG, private T
{
public:

   // Constructor

   VirtualPRNGNew(unsigned start) : T(start) {}
   VirtualPRNGNew() : T() {}

   // Forward all calls

   virtual u32 getMax () const               { return T::getMax(); }
   virtual real   getResolution () const     { return T::getResolution(); }

   virtual u32 operator() ()                 { return T::operator()(); }
   virtual int operator() (int max)          { return T::operator()(max); }
   virtual real getReal()                    { return T::getReal(); }

   virtual HINTLIB_SDN size_t getStateSize () const{ return T::getStateSize(); }
   virtual void saveState (void *p) const    { T::saveState (p); }
   virtual void restoreState (const void *p) { T::restoreState (p); }

   virtual void init (unsigned start)        { T::init (start); }
};


/**
 *  VirtualPRNGNew<T>
 *
 *  A PRNG, using a T* to implement its functionality
 */

template<class T>
class VirtualPRNGWrapper : public VirtualPRNG
{
public:

   // Constructor

   VirtualPRNGWrapper(T* pp) : p(pp) {}

   // Forward all calls

   virtual u32  getMax () const         { return p->getMax(); }
   virtual real getResolution () const  { return p->getResolution(); }

   virtual u32  operator() ()           { return p->operator()(); }
   virtual int  operator() (int max)    { return p->operator()(max); }
   virtual real getReal()               { return p->getReal(); }

   virtual HINTLIB_SDN size_t getStateSize () const{ return p->getStateSize(); }
   virtual void saveState (void *b) const    { p->saveState (b); }
   virtual void restoreState (const void *b) { p->restoreState (b); }

   virtual void init (unsigned start)        { p->init (start); }

private:

   T* const p;
};

}  // namespace HIntLib

#undef HINTLIB_SDN

#endif

