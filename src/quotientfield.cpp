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

#include <HIntLib/quotientfield.tcc>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/integerring.h>
#include <HIntLib/realfield.h>
#include <HIntLib/lookupfield.h>
#include <HIntLib/polynomial2.h>
#include <HIntLib/polynomial.h>

void
HIntLib::Private::printDoubleStruckCapitalQ (std::ostream& o)
{
#if HINTLIB_CHARACTER_SET == 4 && defined (HINTLIB_UTF8_SELECT)
   HINTLIB_UTF8_SELECT(utf8Support(o),
      o << "\xe2\x84\x9a",  // DOUBLE-STRUCK CAPITAL Q
      o << "Q")
#else
   o << "Q";
#endif
}

#ifdef HINTLIB_BUILD_WCHAR
void
HIntLib::Private::printDoubleStruckCapitalQ (std::wostream& o)
{
#if HINTLIB_CHARACTER_SET == 4
   o << L"\x211a";  // DOUBLE-STRUCK CAPITAL Q
#else
   o << L"Q";
#endif
}
#endif

namespace HIntLib
{
   HINTLIB_INSTANTIATE_QUOTIENTFIELD_INT (IntegerRing<int>)
   HINTLIB_INSTANTIATE_QUOTIENTFIELD_POL (PolynomialRing<RealField<real> >)
   HINTLIB_INSTANTIATE_QUOTIENTFIELD_POL (PolynomialRing<ComplexField<real> >)
   HINTLIB_INSTANTIATE_QUOTIENTFIELD_POL
      (PolynomialRing<QuotientField<IntegerRing<int> > >)
   HINTLIB_INSTANTIATE_QUOTIENTFIELD_POL
      (PolynomialRing<LookupField<unsigned char> >)
   HINTLIB_INSTANTIATE_QUOTIENTFIELD_POL (PolynomialRing<GF2>)
   HINTLIB_INSTANTIATE_QUOTIENTFIELD_POL (Polynomial2Ring<u32>)
}

