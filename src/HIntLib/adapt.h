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
 *  Adapt
 *
 */

#ifndef HINTLIB_ADAPT_H
#define HINTLIB_ADAPT_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/regioncollection.h>
#include <HIntLib/region.h>


namespace HIntLib
{

/**
 *  splitHypercube()
 *
 *  Splits a Hypercube into _total_ sub-cubes.
 *  Returns sub-cube # _index_
 */
 
void splitHypercube (Function &, Hypercube &, unsigned total, unsigned index);
void splitHypercube (Function &, const Hypercube &, Hypercube* cubes [],
                     int total);


/**
 *  storeSubcube()
 */

inline
void storeSubcube (
   RegionCollection &rc, Hypercube h, unsigned total, unsigned index,
   Function &function, EmbeddedRule &rule)
{
   splitHypercube (function, h, total, index);

   rc.push (new Region (h, function, rule));
}

}

#endif

