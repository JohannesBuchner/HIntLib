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

/**
 *  eratosthenes.h
 *
 *  Uses Eratosthenes' sieve to find primitive elements in a certain
 *  Euclidean ring.
 *
 *  The template parameter type T has to have the following properties.
 *
 *  - T is an Euclidean ring, i.e.
 *      - commutative
 *      - there is a One
 *      - a * b != 0 for all a!=0 and b!=0
 *      - Every element != Zero is a nonambiguous product of prime elements
 *
 *  - There must be a bijection f:{0,1,2,3,...,n}->T
 *    - A constructor T(unsigned int) must be defined to implement f
 *    - T(0) = Zero
 *    - T(1) = One
 *    - x = p_1 *...* p_n, p_i primitive  ==> f(p_i) < p(x) 
 *
 *  - There must be an operator > on T x T with the following property
 *       a * b > c   ==>   a1 * b1 != c	  for all  a1 with f(a1) >= f(a) 
 *                                        for all  b1 with f(b1) >= f(b)
 *       For all c there should be an x with   x*x > last and x*c > last
 *
 *    Example:
 *       T = unsinged int:  Normal > does the job
 *       T = Polynomial2:   > compares degree of Polynomials
 */

#ifndef HINTLIB_ERATOSTHENES_H
#define HINTLIB_ERATOSTHENES_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif


namespace HIntLib
{

template <class T>
void eratosthenes (bool* field, unsigned size, const T&)
{
   T last (size);

   // Initialize the field - Zero and One are never prime
 
   field [0] = field [1] = false;                                                
   // Assume that all other elements are prime
 
   for (unsigned i = 2; i < size; i++)  field [i] = true;
 
   // Sieve
 
   for (unsigned primeCode = 2; ; primeCode++)
   {
      // if it's not prime, skip it
 
      if (! field [primeCode]) continue;
 
      // Convert it
 
      T prime (primeCode);

      if (prime * prime > last) break;

 
      // Multiply it with every possible value larger than primeCode
 
      for (unsigned i = primeCode; ; i++)
      {
         T product = prime * T (i);

         if (product > last) break;

         unsigned a = static_cast<unsigned> (product);

         if (a < size) field [a] = false;
      }
   }
}

template <class A>
void eratosthenes (A a, bool* field, unsigned size)
{
   typename A::type last (size);

   field [0] = false;  // 0 is never prime 

   unsigned i = 1;

   while (a.isUnit (a.element(i)))  field[i++] = false;

   // Assume that all other elements are prime
 
   for (i = 2; i < size; i++)  field [i] = true;
 
   // Sieve
 
   for (i = 2; ; i++)
   {
      // if it's not prime, skip it
 
      if (! field [i]) continue;
 
      // Convert it
 
      typename A::type prime (i);

      if (a.norm(a.mul (prime,prime)) > a.norm(last)) break;
 
      // Multiply it with ervery possible value larger than primeCode
 
      for (unsigned j = i; ; j++)
      {
         typename A::type product = a.mul (prime * a.element(i));

         if (a.norm(product) > a.norm(last)) break;

         unsigned index = a.index(product);

         if (index < size) field [index] = false;
      }
   }
}

}  // namespace HIntLib

#endif
