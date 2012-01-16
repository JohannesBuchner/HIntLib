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
 *  RegionCollection
 */

#ifndef HINTLIB_REGIONCOLLECTION_H
#define HINTLIB_REGIONCOLLECTION_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <queue>
#include <vector>

#include <HIntLib/region.h>


namespace HIntLib
{

class EmbeddedRule;
class Integrand;

class RegionCollection
   : private std::priority_queue<Region*, std::vector<Region*>,
                                 RegionErrorLess>,
     private EstErr
{
private:
   typedef std::priority_queue<Region*, std::vector <Region*>,
                               RegionErrorLess> Q;
   typedef Q::container_type::iterator I;

public:
   using Q::size_type;

   RegionCollection() : EstErr (0.0, 0.0) {}
   ~RegionCollection();

   using Q::size;
   using Q::empty;

   using EstErr::getEstimate;
   using EstErr::getError;
   using EstErr::getRelError;
   real getTopError () const  { return top()->getError(); }
   const EstErr & getEstErr () const { return *static_cast<const EstErr*>(this); }

   void push (Region*);
   Region*   pop ();
   using Q::top;

   void refine (Integrand &, EmbeddedRule &);

   void killQueueAndUpdateResult ();
};

inline
void RegionCollection::push (Region *r)
{
   operator+= (r->getEstErr());

   Q::push (r);
}

inline
Region* RegionCollection::pop ()
{
   Region* r = top();

   operator-= (r->getEstErr());

   Q::pop();

   return r;
}

inline
void RegionCollection::refine (Integrand &f, EmbeddedRule &rule)
{
   Region *r = pop();   // Get worst region
 
   push (new Region (*r, f, rule));  // Split Region
   push (r);
}


/**
 *  splitHypercube()
 *
 *  Splits a Hypercube into _total_ sub-cubes.
 *  Returns sub-cube # _index_
 */
 
void splitHypercube (Integrand &, Hypercube &, unsigned total, unsigned index);
void splitHypercube (Integrand &, const Hypercube &, Hypercube* cubes [],
                     int total);

/**
 *  storeSubcube()
 */

inline
void storeSubcube (
   RegionCollection &rc, Hypercube h, unsigned total, unsigned index,
   Integrand &function, EmbeddedRule &rule)
{
   splitHypercube (function, h, total, index);

   rc.push (new Region (h, function, rule));
}


}  // namespace HIntLib

#endif

