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

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/output.h>

namespace P = HIntLib::Private;


/**
 *  Printer
 */

P::Printer::Printer (std::ostream &_o)
   : o (_o)
{
   flags (o.flags());
   precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   imbue (o.getloc());
#endif
}

P::Printer::~Printer ()
{
   // o << str(); // string<> does not handle width properly on GNU C++ 2.95.x
   o << str().c_str();
}


