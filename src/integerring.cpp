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

#include <HIntLib/integerring.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/output.h>

namespace L = HIntLib;
namespace P = HIntLib::Private;

/**
 *  opeator<<
 */

std::ostream& P::operator<< (std::ostream &o, const ZRing&)
{
#if HINTLIB_CHARACTER_SET == 4 && defined (HINTLIB_UTF8_SELECT)
   HINTLIB_UTF8_SELECT(Private::utf8Support(o),
      return o << "\xe2\x84\xa4",  // DOUBLE-STRUCK CAPITAL Z
      return o << "Z")
#else
   return o << "Z";
#endif
}

#ifdef HINTLIB_BUILD_WCHAR
std::wostream& P::operator<< (std::wostream &o, const ZRing&)
{
#if HINTLIB_CHARACTER_SET == 4
   return o << L"\x2124";  // DOUBLE-STRUCK CAPITAL Z
#else
   return o << L"Z";
#endif
}
#endif


/**
 *  element()
 */

template<typename T>
T L::IntegerRing<T>::element (unsigned i)
{
   return (i & 1)  ?  T(i/2 + 1)  :  -T(i/2);
}


/**
 *  index()
 */

template<typename T>
unsigned L::IntegerRing<T>::index (const T& x)
{
   return (x > 0)  ?  x*2 - 1  :  x * -2;
}


/**
 *  factor()
 */

template<typename T>
int L::IntegerRing<T>::factor (Factorization& f, type n)
{
   unit_type u = makeCanonical (n);
   PrimeDivisors pd (n);

   unsigned e;
   while (unsigned prime = pd.next(e))
   {
      f.push_back (std::make_pair (type(prime), e));
   }
   return u;
}

namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template X IntegerRing<X>::element(unsigned); \
   template unsigned IntegerRing<X>::index(const X&); \
   template int IntegerRing<X>::factor(Factorization&, type);

   HINTLIB_INSTANTIATE(int)
#undef HINTLIB_INSTANTIATE
}

