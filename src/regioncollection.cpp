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
#pragma implementation "fourthdiff.h"
#endif

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/regioncollection.h>

#include <HIntLib/fourthdiff.h>
#include <HIntLib/bitop.h>
#include <HIntLib/hlalgorithm.h>

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


/**
 *  splitHypercube()
 *
 *  Splits a Hypercube into _total_ sub-cubes.
 *  Returns sub-cube # _index_
 */

void L::splitHypercube (Integrand &f, Hypercube &h,
                     unsigned total, unsigned index)
{
   FourthDiff fd (h.getDimension ());

   for (unsigned int level = 0; ; ++level)
   {
      int mask1 = 1 << level;               //  ..00100..00
      int mask2 = (mask1 * 2) - 1;          //  ..00111..11

      if (((index ^ mask1) & mask2) >= total) break;

      unsigned int split = fd (f, h);

      if (index & mask1)
         h.cutRight (split);
      else
         h.cutLeft (split);
   }
}

void L::splitHypercube (Integrand &f, const Hypercube &h, Hypercube* cubes [],
                        int total)
{
   cubes [0] = new Hypercube (h);

   FourthDiff fd (h.getDimension ());

   int levels = ms1 (total - 1);

   for (int level = 0; level <= levels; ++level)
   {
      int mask = 1u << level;

      for (int i = 0; i < mask; i++)
      {
         if (level != levels || i + mask < total)
         {
            unsigned split = fd (f, *cubes [i]);

            cubes [i + mask] = new Hypercube (*cubes [i], split);
         }
      }
   }
}

