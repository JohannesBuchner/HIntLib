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

#if defined __GNUG__ && ! defined HINTLIB_PARALLEL
#pragma implementation
#endif

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_OSTREAM
  #include <ostream>
#else
  #include <iostream>
#endif

#ifdef HINTLIB_HAVE_SSTREAM
  #include <sstream>
#else
  #include <HIntLib/fallback_sstream.h>
#endif

#include <HIntLib/region.h>
#include <HIntLib/integerring.h>
#include <HIntLib/modulararithmetic.h>

namespace L = HIntLib;

using std::ostream;
using std::ostringstream;


/**
 *  This file contains  operator<< (ostrem&, const X &)  implementations
 *  for various types.
 */


/**
 *  EstErr
 */

ostream& L::operator<< (ostream &o, const EstErr &ee)
{
   ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   ss << ee.getEstimate () << "(+/-" << ee.getError () << ')';
 
   return o << ss.str();
}


/**
 *  Hypercube
 */

ostream & L::operator<< (ostream &o, const Hypercube &h)
{
   ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   ss << '[' << h.getLowerBound (0) << ',' << h.getUpperBound (0) << ']';

   for (unsigned i = 1; i < h.getDimension (); ++i)
   {
      ss << "x[" << h.getLowerBound (i) << ',' << h.getUpperBound (i) << ']';
   }

   return o << ss.str();
}


/**
 *  Region
 */

ostream & L::operator<< (ostream &o, const Region &r)
{
   ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   ss << r.getHypercube () << ' ' << r.getEstErr();

   return o << ss.str();
}


/**
 *  ZRing
 *  RRing
 */

ostream & L::operator<< (ostream &o, const ZRing &)
{
   return o << 'Z';
}

ostream & L::operator<< (ostream &o, const RRing &)
{
   return o << 'R';
}


