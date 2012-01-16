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

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/polynomial2base.h>
#include <HIntLib/gf2vectorspace.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#pragma implementation "gf2.h"
#pragma implementation "gf2vectorspace.h"
#endif

#include <HIntLib/output.h>

namespace L = HIntLib;
namespace P = HIntLib::Private;

/******************************  GF2  ****************************************/


std::ostream &
L::operator<< (std::ostream &o, const GF2 &)
{
   return o << "GF2";
}
#ifdef HINTLIB_BUILD_WCHAR
std::wostream &
L::operator<< (std::wostream &o, const GF2 &)
{
   return o << L"GF2";
}
#endif

void
L::GF2::print (std::ostream& o, type x)
{
   Private::Printer ss (o);
   ss << unsigned (x) << " (2)";
}
#ifdef HINTLIB_BUILD_WCHAR
void
L::GF2::print (std::wostream& o, type x)
{
   Private::WPrinter ss (o);
   ss << unsigned (x) << L" (2)";
}
#endif

void
L::GF2::printShort (std::ostream& o, type x)
{
   o << unsigned (x);
}
#ifdef HINTLIB_BUILD_WCHAR
void
L::GF2::printShort (std::wostream& o, type x)
{
   o << unsigned (x);
}
#endif

void
L::GF2::printSuffix (std::ostream& o)
{
   o << "(2)";
}
#ifdef HINTLIB_BUILD_WCHAR
void
L::GF2::printSuffix (std::wostream& o)
{
   o << L"(2)";
}
#endif


/**********************  GF 2 Vector Space Base  *****************************/


std::ostream &
P::operator<< (std::ostream &o, const GF2VectorSpaceBase &v)
{
   Private::Printer ss (o);
   ss << "(GF2)";
   ss.power (v.dimension());
   return o;
}
#ifdef HINTLIB_BUILD_WCHAR
std::wostream &
P::operator<< (std::wostream &o, const GF2VectorSpaceBase &v)
{
   Private::WPrinter ss (o);
   ss << L"(GF2)";
   ss.power (v.dimension());
   return o;
}
#endif


/**********************  Polynomial 2 Ring Base  *****************************/


std::ostream &
P::operator<< (std::ostream &o, const Polynomial2RingBase &r)
{
   Private::Printer ss (o);

   ss << "GF2[";
   r.printVariable (ss);
   ss << ']';

   return o;
}
#ifdef HINTLIB_BUILD_WCHAR
std::wostream &
P::operator<< (std::wostream &o, const Polynomial2RingBase &r)
{
   Private::WPrinter ss (o);

   ss << L"GF2[";
   r.printVariable (ss);
   ss << L']';

   return o;
}
#endif

void
P::Polynomial2RingBase::printSuffix (std::ostream& o)
{
   o << "(2)";
}
#ifdef HINTLIB_BUILD_WCHAR
void
P::Polynomial2RingBase::printSuffix (std::wostream& o)
{
   o << L"(2)";
}
#endif

void
P::Polynomial2RingBase::printVariable (std::ostream& o) const
{
   o << var;
}
#ifdef HINTLIB_BUILD_WCHAR
void
P::Polynomial2RingBase::printVariable (std::wostream& o) const
{
   o << wvar;
}
#endif


