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

#ifndef BIT_VECTOR_RING_H
#define BIT_VECTOR_RING_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/defaults.h>

#ifdef HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/hintlib_limits.h>
#endif

namespace HIntLib
{

template <class T>
class BitVectorRing
{
public:

   typedef T type;

   BitVectorRing (unsigned bits)
      : mask (~T(-1) >> (std::numeric_limits<T>::digits - bits))
   {
      if (bits < 1 || bits >= unsigned (std::numeric_limits<T>::digits))
      {
         throw 1;
      }
   } 

   T zero() const  { return T(0); }
   T one()  const  { return mask; }

   T& inc (T& a) const  { return a ^= mask; }
   T& dec (T& a) const  { return a ^= mask; }

   T add (const T& a, const T& b) const  { return a | b; }
   T& addTo (T& a,    const T& b) const  { return a |= b; }

   T neg (const T& a) const  { return a; }
   T& negate (T& a) const    { return a; }

   T sub (const T& a, const T& b) const  { return a | b; }
   T& subFrom (T& a,  const T& b) const  { return a |= b; }

   T mul (const T& a, const T& b) const  { return a & b; }
   T& mulBy (T& a,    const T& b) const  { return a &= b; }

   T pow (const T& a, unsigned k) const  { return a; }

   T size() const  { return mask + 1; }

private:
   const T mask;
};


}  // namespace HIntLib

#endif

