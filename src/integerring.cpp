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

#include <HIntLib/integerring.h>

#ifdef HINTLIB_HAVE_OSTREAM
  #include <ostream>
#else
  #include <iostream>
#endif

namespace L = HIntLib;
namespace P = HIntLib::Private;

/**
 *  opeator<<
 */

std::ostream& P::operator<< (std::ostream &o, const ZRing&)
{
   return o << "Z";
}


/**
 *  element()
 */

template<typename T>
T L::IntegerRing<T>::element (unsigned i)
{
   return odd(i)  ?  T(i/2 + 1)  :  -T(i/2);
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

