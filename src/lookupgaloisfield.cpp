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

#include <HIntLib/lookupfield.h>
#include <HIntLib/galoisfield.h>
#include <HIntLib/prime.h>

namespace L = HIntLib;

/**
 *  make Galois Field ()
 */

template<class T>
void L::makeGaloisField
   (T &r, unsigned base, unsigned exponent)
{
        if (exponent == 0)  throw GaloisFieldExponent ();
   else if (exponent == 1)
   {
      FactorField<typename T::type> field (base);
      L::copy (r, field);
   }
   else
   {
      GaloisField<typename T::type> field (base, exponent);
      L::copy (r, field);
   }
}


/**
 *  Lookup Galois Field
 */

template<class T>
L::LookupGaloisField<T>::LookupGaloisField (unsigned prime, unsigned power)
   : LookupField<T> (powInt (prime, power))
{
   if (getRefCount())  makeGaloisField (*this, prime, power);
}

template<class T>
L::LookupGaloisField<T>::LookupGaloisField (unsigned size)
   : LookupField<T> (size)
{
   if (getRefCount())
   {
      unsigned prime, power;
      Prime::factorPrimePower (size, prime, power);
      makeGaloisField (*this, prime, power);
   }
}


/**
 *  Lookup Galois Field 2
 */

template<class T>
L::LookupGaloisFieldPow2<T>::LookupGaloisFieldPow2
   (unsigned prime, unsigned power)
   : LookupFieldPow2<T> (1 << power)
{
   if (prime != 2)  throw FIXME (__FILE__, __LINE__);
   if (getRefCount())  makeGaloisField (*this, prime, power);
}

template<class T>
L::LookupGaloisFieldPow2<T>::LookupGaloisFieldPow2 (unsigned size)
   : LookupFieldPow2<T> (size)
{
   unsigned prime, power;

   if (getRefCount())
   {
      Prime::factorPrimePower (size, prime, power);
      makeGaloisField (*this, prime, power);
   }
   else
   {
      prime = characteristic();
   }

   if (prime != 2)  throw FIXME (__FILE__, __LINE__);
}


/**
 *  Lookup Galois Field Prime
 */

template<class T>
L::LookupGaloisFieldPrime<T>::LookupGaloisFieldPrime
   (unsigned prime, unsigned power)
   : LookupFieldPrime<T> (prime)
{
   if (power != 1)  throw FIXME (__FILE__, __LINE__);
   if (getRefCount())  makeGaloisField (*this, prime, power);
}

template<class T>
L::LookupGaloisFieldPrime<T>::LookupGaloisFieldPrime (unsigned size)
   : LookupFieldPrime<T> (size)
{
   if (! Prime::test(size))  throw FIXME (__FILE__, __LINE__);
   if (getRefCount())  makeGaloisField (*this, size, 1);
}


/**
 *  copy ()
 */

namespace
{

template<class T, class A>
void copyMulTable (L::LookupFieldBase<T> & dest, const A src)
{
   typedef typename A::type AT;

   if (   dest.size() != src.size()
       || ! src.is0 (src.element(0))
       || ! src.is1 (src.element(1)))  throw L::LookupFieldCopy();

   dest.setCharacteristic (src.characteristic(), src.extensionDegree());

   for (unsigned i = 0; i < src.size(); ++i)
   {
      AT a = src.element (i);

      if (i > 0)  dest.setRecip (i, src.index (src.recip (a)));
      if (i > 0)  dest.setOrder (i, src.order (a));

      for (unsigned j = 0; j <= i; ++j)
      {
         AT b = src.element (j);

         unsigned prod = src.index (src.mul (a, b));
         dest.setMul (i, j, prod);
         dest.setMul (j, i, prod);
      }
   }
}

template<class T, class A>
void copyAddTable (L::LookupField<T> & dest, const A src)
{
   typedef typename A::type AT;

   for (unsigned i = 0; i < src.size(); ++i)
   {
      AT a = src.element (i);

      dest.setNeg (i, src.index (src.neg (a)));

      for (unsigned j = 0; j <= i; ++j)
      {
         AT b = src.element (j);

         unsigned sum = src.index (src.add (a, b));
         dest.setAdd (i, j, sum);
         dest.setAdd (j, i, sum);
      }
   }
}
}

template<class T, class A>
void L::copy (LookupField<T> & dest, const A src)
{
   copyMulTable (dest, src);
   copyAddTable (dest, src);
}

template<class T, class A>
void L::copy (LookupFieldMulOnly<T> & dest, const A src)
{
   copyMulTable (dest, src);
}

namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template void copy (LookupField<X> &, const GaloisField<X>); \
   template void copy (LookupField<X> &, const FactorField<X>);\
   template void copy (LookupFieldMulOnly<X> &, const GaloisField<X>); \
   template void copy (LookupFieldMulOnly<X> &, const FactorField<X>);\
   template void makeGaloisField(LookupGaloisField<X>&,unsigned,unsigned); \
   template void makeGaloisField(LookupGaloisFieldPow2<X>&,unsigned,unsigned);\
   template void makeGaloisField(LookupGaloisFieldPrime<X>&,unsigned,unsigned);\
   template class LookupGaloisField<X>; \
   template class LookupGaloisFieldPow2<X>; \
   template class LookupGaloisFieldPrime<X>;

   HINTLIB_INSTANTIATE (unsigned char)
#undef HINTLIB_INSTANTIATE
}

