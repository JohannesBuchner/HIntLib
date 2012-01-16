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

#ifndef HINTLIB_ALGEBRA_H
#define HINTLIB_ALGEBRA_H 1

#ifdef __GNUG__
#pragma interface
#endif

namespace HIntLib
{
   /**
    *  Types for polynomial_tag
    */

   struct nopolynomial_tag {};
   struct polynomial_tag : public nopolynomial_tag {};

   /**
    *  Types for algebra_tag
    */

   struct group_tag {};
   struct ring_tag : public group_tag      { typedef void* MAGIC;};
   struct domain_tag : public ring_tag     { typedef void** MAGIC;};

   struct ufd_tag : public domain_tag      {};
   struct euclidean_tag : public ufd_tag   {};
   struct integer_tag : public euclidean_tag {};

   struct field_tag : public domain_tag    { typedef void*** MAGIC;};
   struct rational_tag : public field_tag  { typedef unsigned*** MAGIC; };
   struct real_tag     : public field_tag  { typedef int *** MAGIC; };
   struct complex_tag  : public field_tag  { typedef char*** MAGIC; };

   struct gf_tag : public field_tag        { typedef void**** MAGIC;};
   struct cyclic_tag : public gf_tag       { typedef void***** MAGIC;};

   struct vectorspace_tag : public group_tag {};

   /**
    *  Addition based on bit operations
    */

   template<typename T>
   class BitOpBasedAddition
   {
   public:
      static bool is0 (const T& a)  { return !a; }

      static T  add (const T& a, const T& b)  { return a ^ b; }
      static T  sub (const T& a, const T& b)  { return a ^ b; }

      static T& addTo   (T& a, const T& b)  { return a ^= b; }
      static T& subFrom (T& a, const T& b)  { return a ^= b; }

      static T  neg    (const T& a)  { return a; }
      static T& negate (      T& a)  { return a; }

      static T times (const T& a, unsigned k)  { return odd (k) ? a : T(); }
      static unsigned additiveOrder (const T& a)  { return a ? 2 : 1; }
      static unsigned characteristic()  { return 2; }

   protected:
      BitOpBasedAddition () {}
   };

} // namespace HIntLib

#define HINTLIB_TRIVIAL_DOMAIN_MEMBERS \
   unsigned additiveOrder (const type& a) const \
      { return is0(a) ? 1 : characteristic(); }

#define HINTLIB_TRIVIAL_CYCLIC_MEMBERS \
   HINTLIB_TRIVIAL_DOMAIN_MEMBERS \
   static unsigned extensionDegree()  { return 1; }

#endif

