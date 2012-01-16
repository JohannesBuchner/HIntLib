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


#include <algorithm>
#include <iostream>
#include <iomanip>

#include <HIntLib/precalculatedfield.h>

#include <HIntLib/exception.h>

namespace L = HIntLib;

using std::setw;

/**
 *  Printing
 */

std::ostream& L::operator<< (std::ostream &o, const PrecalculatedFieldBase &r)
{
   return o << "Field with " << r.size() << " elements";
}


/**
 *  setPointers()
 */

template<class T>
void L::PrecalculatedField<T>::setPointers ()
{
   addTable = shared + 1;
   mulTable = addTable + sqr(s);
   negTable = mulTable + sqr(s);
   recTable = negTable + s;
}


/**
 *  copy()
 */

template<class T>
void L::PrecalculatedField<T>::copy ()
{
   if (*shared)
   {
      --*shared;

      const T* old = addTable;
      shared   = new T [2 * s * (s+1) + 1];
      setPointers();

      *shared = 0;
      std::copy (old, old + 2 * s * (s+1), addTable);
   }
}


/**
 *  destroy()
 */

template<class T>
void L::PrecalculatedField<T>::destroy ()
{
   if (! *shared)
   {
      delete[] shared;
   }
   else --*shared;
}


/**
 *  Constructor
 */

template<class T>
L::PrecalculatedField<T>::PrecalculatedField (unsigned _s)
   : PrecalculatedFieldBase (_s)
{
   if (     std::numeric_limits<T>::is_signed
       || ! std::numeric_limits<T>::is_integer)
   {
      throw InvalidType ("PrecalculatedField");
   }

   shared   = new T [2 * s * (s+1) + 1];
   setPointers();

   *shared = 0;

   std::fill (addTable, addTable + 2 * s * (s+1), T(0));
} 


/**
 *  Copy Constructor
 */

template<class T>
L::PrecalculatedField<T>::PrecalculatedField (const PrecalculatedField<T> &r)
   : PrecalculatedFieldBase (r)
{
   shared = r.shared;
   setPointers();
   ++ *shared;
}

/**
 *  operator=
 */

template<class T>
L::PrecalculatedField<T> &
L::PrecalculatedField<T>::operator= (const PrecalculatedField<T> &r)
{
   if (shared != r.shared)
   {
      destroy();
      shared   = r.shared;
      setPointers();
      ++ *shared;
   }
   return *this;
}


/**
 *  operator==
 */

template<class T>
bool L::PrecalculatedField<T>::operator== (const PrecalculatedField<T> &r) const
{
   return    size() == r.size()
          && std::equal (addTable, addTable + 2 * s * (s+1), r.addTable);
}


/**
 *  times()
 */

template<class T>
T  L::PrecalculatedField<T>::times (T x, unsigned k) const
{
   T result = zero();
   k %= characteristic;

   for(;;)
   {
      if (k & 1)  addTo (result, x);
      if ((k >>= 1) == 0)  return result;
      addTo (x, x);
   }
}


/**
 *  dump()
 */

template<class T>
void L::PrecalculatedField<T>::dump (std::ostream &o) const
{
   o << setw (3) << "+";
   for (unsigned j = 0; j < s; ++j)  o << setw(3) << j;
   o << "   shared: " << unsigned (*shared) << '\n';

   for (unsigned i = 0; i < s; ++i)
   {
      o << setw(3) << i;

      for (unsigned j = 0; j < s; ++j)
      {
         o << setw(3) << unsigned (add (i,j));
      }
      o << '\n';
   }

   o << '\n' << setw (3) << "*";
   for (unsigned j = 0; j < s; ++j)  o << setw(3) << j;
   o << '\n' << setw (3) << "-1";
   for (unsigned j = 0; j < s; ++j)  o << setw(3) << unsigned (neg (j));
   o << '\n';

   for (unsigned i = 0; i < s; ++i)
   {
      o << setw(3) << i;

      for (unsigned j = 0; j < s; ++j)
      {
         o << setw(3) << unsigned (mul (i,j));
      }
      o << '\n';
   }
   o << '\n';
}
   

/**
 *  setAdd ()
 *  setMul ()
 *  setNeg ()
 *  setRecip ()
 */

template<class T>
void L::PrecalculatedField<T>::setAdd (T a, T b, T sum)
{
   if (a   >= s)  throw PrecalculatedFieldSet (a,   s);
   if (b   >= s)  throw PrecalculatedFieldSet (b,   s);
   if (sum >= s)  throw PrecalculatedFieldSet (sum, s);
   copy();
   addTable [a + s * b] = sum;
}

template<class T>
void L::PrecalculatedField<T>::setMul (T a, T b, T prod)
{
   if (a    >= s)  throw PrecalculatedFieldSet (a,    s);
   if (b    >= s)  throw PrecalculatedFieldSet (b,    s);
   if (prod >= s)  throw PrecalculatedFieldSet (prod, s);
   copy();
   mulTable [a + s * b] = prod;
}

template<class T>
void L::PrecalculatedField<T>::setNeg (T a, T minus)
{
   if (a     >= s)  throw PrecalculatedFieldSet (a,     s);
   if (minus >= s)  throw PrecalculatedFieldSet (minus, s);
   copy();
   negTable [a] = minus;
}

template<class T>
void L::PrecalculatedField<T>::setRecip (T a, T r)
{
   if (a >= s)  throw PrecalculatedFieldSet (a, s);
   if (r >= s)  throw PrecalculatedFieldSet (r, s);
   copy();
   recTable [a] = r;
}


namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template void PrecalculatedField<X>::copy (); \
   template void PrecalculatedField<X>::destroy (); \
   template PrecalculatedField<X>::PrecalculatedField (unsigned); \
   template PrecalculatedField<X>::PrecalculatedField \
        (const PrecalculatedField<X>&); \
   template PrecalculatedField<X> & \
      PrecalculatedField<X>::operator= (const PrecalculatedField<X> &); \
   template bool PrecalculatedField<X>::operator== \
        (const PrecalculatedField<X> &) const; \
   template void PrecalculatedField<X>::dump (std::ostream &) const; \
   template X    PrecalculatedField<X>::times (X, unsigned) const; \
   template void PrecalculatedField<X>::setAdd (X, X, X); \
   template void PrecalculatedField<X>::setMul (X, X, X); \
   template void PrecalculatedField<X>::setNeg (X, X); \
   template void PrecalculatedField<X>::setRecip (X, X);

   HINTLIB_INSTANTIATE (unsigned char);
   HINTLIB_INSTANTIATE (unsigned short);
#undef HINTLIB_INSTANTIATE
}

