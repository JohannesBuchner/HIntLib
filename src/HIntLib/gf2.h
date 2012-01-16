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

#ifndef HINTLIB_GF2_H
#define HINTLIB_GF2_H 1

#ifdef __GNUG__
#pragma interface
// Implmentaion in polynomial2.cpp
#endif

#include <iosfwd>

#include <HIntLib/algebra.h>


namespace HIntLib
{

/**
 *  GF2
 *
 *  The field with two elements - 0 and 1.
 */

class GF2 : public BitOpBasedAddition<unsigned char>
{
public:
   GF2 () {}
   GF2 (unsigned) {}

   typedef unsigned char type;
   typedef cyclic_tag algebra_category;
   typedef nopolynomial_tag polynomial_category;
   typedef char_two char_category;
   typedef nozerodivisor_tag zerodivisor_category;
   typedef finite_tag size_category;

   static unsigned size()  { return 2; }
   static unsigned modulus()  { return 2; }
   static unsigned extensionDegree()  { return 1; }

   static type one()  { return 1; }

   static bool is1 (type a)  { return a; }

   static type element  (unsigned i)  { return type(i); }
   static unsigned index  (type x)  { return x; }

   static type times (const type& a, unsigned k)  { return a & k; }

   static type mul   (const type& a, const type& b)  { return a & b; }
   static void mulBy (      type& a, const type& b)  { a &= b; }

   static type sqr (const type& a)  { return a; }
   static void square (type&)     {}
   static type power (const type& a, unsigned)  { return a; }

   static type recip   (const type&)  { return 1; }
   static void reciprocal (type&)  {}

   static type div   (const type& a, const type&)  { return a; }
   static void divBy (      type&, const type&)  {}

   static type    frobenius (const type& a)  { return a; }
   static type invFrobenius (const type& a)  { return a; }

   static unsigned additiveOrder (const type& a)  { return a + 1; }
   static unsigned order (const type&)  { return 1; }
   static bool isPrimitiveElement (const type& a)  { return a; }
   
   static void print (std::ostream&, type);
   static void printShort (std::ostream&, type);
   static void printShort (std::ostream& o, const type& a, PrintShortFlag)
      { printShort (o, a); }
   static void printSuffix (std::ostream&);
};

std::ostream& operator<< (std::ostream&, const GF2&);


} // namespace HIntLib

#endif

