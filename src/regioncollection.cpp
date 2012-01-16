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


#ifdef __GNUG__
#pragma implementation
#endif

#include <HIntLib/regioncollection.h>

#include <HIntLib/myalgorithm.h>

namespace L = HIntLib;

void L::RegionCollection::killQueueAndUpdateResult (void)
{
   set (0.0, 0.0);
 
   for (I i = c.begin (); i != c.end (); ++i)
   {
      operator+= ((*i)->getEstErr ());
      delete *i;
   }
 
   c.clear ();
}

L::RegionCollection::~RegionCollection (void)
{
   purge (c.begin(), c.end());
}

