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

#ifndef REGIONCOLLECTION_H
#define REGIONCOLLECTION_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <stddef.h>
#include <queue>
#include <vector>

#include <HIntLib/region.h>


namespace HIntLib
{

class EmbeddedRule;
class Function;

class RegionCollection
   : private std::priority_queue<Region*, std::vector<Region*>,
                                 RegionErrorLess>
{
private:
   typedef std::priority_queue<Region*, std::vector <Region*>,
                               RegionErrorLess> Q;
   typedef Q::container_type::iterator I;

public:

   RegionCollection() : ee (0.0, 0.0) {}
   ~RegionCollection();

   size_t size () const  { return Q::size(); }
   bool empty () const   { return Q::empty(); }

   real getEstimate () const  { return ee.getEstimate(); }
   real getError () const     { return ee.getError(); }
   real getRelError () const  { return ee.getRelError(); }
   real getTopError () const  { return Q::top()->getError(); }
   const EstErr & getEstErr () const { return ee; }

   void push (Region*);
   Region*   pop ();
   Region*   top () const  { return Q::top(); }

   void refine (Function &, EmbeddedRule &);

   void killQueueAndUpdateResult ();

private:

   EstErr ee;
};

inline
void RegionCollection::push (Region *r)
{
   ee += r->getEstErr();

   Q::push (r);
}

inline
Region* RegionCollection::pop ()
{
   Region* r = Q::top();

   ee -= r->getEstErr();

   Q::pop();

   return r;
}

inline
void RegionCollection::refine (Function &f, EmbeddedRule &rule)
{
   Region *r = pop();   // Get worst region
 
   push (new Region (*r, f, rule));  // Split Region
   push (r);
}

}  // namespace HIntLib

#endif

