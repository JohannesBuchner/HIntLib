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

/**
 *  Region
 */

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/region.h>

#if defined HINTLIB_USE_INTERFACE_IMPLEMENTATION && ! defined HINTLIB_PARALLEL
#pragma implementation
#pragma implementation "esterr.h"
#endif

#include <HIntLib/output.h>

namespace L = HIntLib;

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
 *  operator<<  for  EstErr
 */

std::ostream& L::operator<< (std::ostream &o, const EstErr &ee)
{
   Private::Printer ss (o);
   ss << ee.getEstimate () << '(';
#if defined(HINTLIB_UTF8_SELECT)
   // PLUS-MINUS SIGN
   HINTLIB_UTF8_SELECT(ss.utf8(), ss << "\xc2\xb1", ss << "\xb1")
#else
   ss << "+/- ";
#endif
   ss << ee.getError () << ')';
   return o;
}

#ifdef HINTLIB_BUILD_WCHAR
std::wostream& L::operator<< (std::wostream &o, const EstErr &ee)
{
   Private::WPrinter ss (o);
   ss << ee.getEstimate () << L'(';
#if HINTLIB_CHARACTER_SET >= 2
   ss << L'\x00b1';  // PLUS-MINUS SIGN
#else
   ss << L"+/- ";
#endif
   ss << ee.getError () << L')';
   return o;
}
#endif


/**
 *  operator<<  for  Region
 */

std::ostream & L::operator<< (std::ostream &o, const Region &r)
{
   Private::Printer ss (o);
   ss << r.getHypercube () << ' ' << r.getEstErr();
   return o;
}

#ifdef HINTLIB_BUILD_WCHAR
std::wostream & L::operator<< (std::wostream &o, const Region &r)
{
   Private::WPrinter ss (o);
   ss << r.getHypercube () << ' ' << r.getEstErr();
   return o;
}
#endif


/**
 *  operator<<  for  Hypercube
 */

std::ostream & L::operator<< (std::ostream &o, const Hypercube &h)
{
   Private::Printer ss (o);

   ss << '[' << h.getLowerBound (0) << ',' << h.getUpperBound (0) << ']';

   for (int i = 1; i < h.getDimension (); ++i)
   {
#if defined(HINTLIB_UTF8_SELECT)
      // MULTIPLICATION SIGN
      HINTLIB_UTF8_SELECT(ss.utf8(), ss << "\xc3\x97", ss << "\xd7")
#else
      ss << 'x';
#endif
      ss << '[' << h.getLowerBound (i) << ',' << h.getUpperBound (i) << ']';
   }

   return o;
}
#ifdef HINTLIB_BUILD_WCHAR
std::wostream & L::operator<< (std::wostream &o, const Hypercube &h)
{
   Private::WPrinter ss (o);

   ss << L'[' << h.getLowerBound (0) << L',' << h.getUpperBound (0) << L']';

   for (int i = 1; i < h.getDimension (); ++i)
   {
#if HINTLIB_CHARACTER_SET >= 2
      ss << L'\x00d7';  // MULTIPLICATION SIGN
#else
      ss << L'x';
#endif
      ss << L'[' << h.getLowerBound (i) << L','
         << h.getUpperBound (i) << L']';
   }

   return o;
}
#endif

