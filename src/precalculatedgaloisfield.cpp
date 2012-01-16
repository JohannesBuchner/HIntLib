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

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_OSTREM
  #include <ostream>
#else
  #include <iostream>
#endif

#include <HIntLib/precalculatedfield.h>
#include <HIntLib/galoisfield.h>
#include <HIntLib/prime.h>

namespace L = HIntLib;


/**
 *  make Galois Field ()
 */

template<class T>
void L::makeGaloisField
   (L::PrecalculatedField<T> &r, unsigned base, unsigned exponent)
{
        if (exponent == 0)  throw GaloisFieldExponent ();
   else if (exponent == 1)
   {
      ModularIntegerField<T> field (base);
      L::copy (r, field);
   }
   else
   {
      GaloisField<T> gf (base, exponent);
      L::copy (r, gf);
   }
}


/**
 *  Constructor
 */

template<class T>
L::PrecalculatedGaloisField<T>::PrecalculatedGaloisField
   (unsigned prime, unsigned power)
   : PrecalculatedField<T> (powInt (prime, power))
{
   makeGaloisField (*this, prime, power);
   characteristic = prime;
}

template<class T>
L::PrecalculatedGaloisField<T>::PrecalculatedGaloisField (unsigned size)
   : PrecalculatedField<T> (size)
{
   unsigned power;
   Prime::factorPrimePower (size, characteristic, power);
   makeGaloisField (*this, characteristic, power);
}


/**
 *  copy ()
 */

template<class T, class A>
void L::copy (L::PrecalculatedField<T> & dest, const A src)
{
   typedef typename A::type AT;

   if (   dest.size() != src.size()
       || ! src.is0 (src.element(0))
       || ! src.is1 (src.element(1)))  throw PrecalculatedFieldCopy();

   for (unsigned i = 0; i < src.size(); ++i)
   {
      AT a = src.element (i);

      dest.setNeg (i, src.index (src.neg (a)));

      if (i > 0)  dest.setRecip (i, src.index (src.recip (a)));

      for (unsigned j = 0; j <= i; ++j)
      {
         AT b = src.element (j);

         unsigned sum = src.index (src.add (a, b));
         dest.setAdd (i, j, sum);
         dest.setAdd (j, i, sum);

         unsigned prod = src.index (src.mul (a, b));
         dest.setMul (i, j, prod);
         dest.setMul (j, i, prod);
      }
   }
}


namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template void copy (PrecalculatedField<X> &, const GaloisField<X>); \
   template void copy (PrecalculatedField<X> &, const ModularIntegerField<X>);\
   template void makeGaloisField (PrecalculatedField<X> &,unsigned,unsigned); \
   template class PrecalculatedGaloisField<X>;

   HINTLIB_INSTANTIATE (unsigned char)
   HINTLIB_INSTANTIATE (unsigned short)
#undef HINTLIB_INSTANTIATE
}

