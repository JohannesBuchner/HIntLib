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
 *  Region
 */

#if defined __GNUG__ && ! defined HINTLIB_PARALLEL
#pragma implementation
#pragma implementation "esterr.h"
#endif

#include <HIntLib/region.h>

#include <HIntLib/output.h>

namespace L = HIntLib;
using std::ostream;

/**
 *  Print the data of a region
 */

L::Region::Region (Region &r, Integrand &f, EmbeddedRule &rule)
  : h (r.h, r.splitDim),
    numOfSplits (++r.numOfSplits)
{
     eval (f, rule);
   r.eval (f, rule);
}


/**
 *  operator<< ()
 *
 *  for Hypercube, Region, and EstErr
 */

ostream& L::operator<< (ostream &o, const EstErr &ee)
{
   Private::Printer ss (o);
   ss << ee.getEstimate () << "(+/-" << ee.getError () << ')';
   return o;
}

ostream & L::operator<< (ostream &o, const Region &r)
{
   Private::Printer ss (o);
   ss << r.getHypercube () << ' ' << r.getEstErr();
   return o;
}

ostream & L::operator<< (ostream &o, const Hypercube &h)
{
   Private::Printer ss (o);

   ss << '[' << h.getLowerBound (0) << ',' << h.getUpperBound (0) << ']';

   for (unsigned i = 1; i < h.getDimension (); ++i)
   {
      ss << "x[" << h.getLowerBound (i) << ',' << h.getUpperBound (i) << ']';
   }

   return o;
}

