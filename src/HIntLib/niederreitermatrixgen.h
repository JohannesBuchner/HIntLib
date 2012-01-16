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

#ifndef NIEDERREITER_MATRIX_GEN_H
#define NIEDERREITER_MATRIX_GEN_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/generatormatrix.h>
#include <HIntLib/polynomial.h>
#include <HIntLib/mymath.h>

namespace HIntLib
{

   template<class A>
   class NiederreiterMatrixGen
      : public HeapAllocatedGeneratorMatrixGen<typename A::type>
   {
   private:
      typedef typename A::type T;
      typedef PolynomialRingField<A> PolyRing;
      typedef typename PolyRing::type Poly;

   public:
      static void init (
         MutableGeneratorMatrixGen<typename A::type> &, A);
      static void init (
         MutableGeneratorMatrixGen<typename A::type> &,
         A, unsigned d, const Poly &);
      
      void init (const A &a , unsigned d, const Poly &p)
         {  init (*this, a, d, p); }

      NiederreiterMatrixGen
         (const A &a, unsigned _dim, unsigned _m, unsigned _prec);
   };

   class NiederreiterMatrixPP
      : public HeapAllocatedGeneratorMatrixGen<unsigned char>
   {
   public:
      NiederreiterMatrixPP (unsigned _dim, unsigned char size)
         : HeapAllocatedGeneratorMatrixGen<unsigned char> (size, _dim)
         {  init (size); }
      NiederreiterMatrixPP (unsigned _dim, unsigned char prime, unsigned power)
         : HeapAllocatedGeneratorMatrixGen<unsigned char>
                  (powInt (unsigned (prime), power), _dim)
         {  init (prime, power); }

      NiederreiterMatrixPP (
         unsigned _dim, unsigned _m, unsigned _prec,
         unsigned char size)
      : HeapAllocatedGeneratorMatrixGen<unsigned char> (size, _dim, _m, _prec)
         {  init (size); }
      NiederreiterMatrixPP (
         unsigned _dim, unsigned _m, unsigned _prec,
         unsigned char prime, unsigned power)
      : HeapAllocatedGeneratorMatrixGen<unsigned char>
               (powInt (unsigned (prime), power), _dim, _m, _prec)
         {  init (prime, power); }

   private:
      void init (unsigned char prime, unsigned power);
      void init (unsigned char size);
   };

}   // namespace HIntLib

#endif


