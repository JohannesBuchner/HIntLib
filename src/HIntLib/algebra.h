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
    *  Tags used in algebraic structures
    */

   // Tags for algebra_category

   struct group_tag {};
   struct ringfield_tag : public group_tag {};

   struct ringdomain_tag : public ringfield_tag {};
   struct ring_tag : public ringdomain_tag   { typedef ring_tag P; };
   struct domain_tag : public ringdomain_tag  { typedef domain_tag P; };
   struct ufd_tag : public domain_tag      { typedef ufd_tag P; };
   struct euclidean_tag : public ufd_tag {};
   struct integer_tag : public euclidean_tag {};

   struct field_tag       : public ringfield_tag
                        { typedef field_tag P; typedef P F; };
   struct gf_tag          : public field_tag { typedef gf_tag P; typedef P F; };
   struct cyclic_tag      : public gf_tag {};
   struct realcomplex_tag : public field_tag { typedef realcomplex_tag F; };
   struct real_tag        : public realcomplex_tag { typedef real_tag P; };
   struct complex_tag     : public realcomplex_tag { typedef complex_tag P; };
   struct numberfield_tag : public field_tag { typedef numberfield_tag F; };
   struct rational_tag    : public numberfield_tag { typedef rational_tag P; };
   struct funfield_tag    : public field_tag { typedef funfield_tag F; };
   struct ratfunfield_tag : public funfield_tag {};

   struct vectorspace_tag : public group_tag {};

   // Tags for polynomial_category

   struct nopolynomial_tag {};
   struct polynomial_tag {};

   // Tags for char_category

   struct char_any {};
   struct char_non    : public char_any {}; 
   struct char_exists : public char_any {};
   struct char_zero  : public char_exists {};
   struct char_prime : public char_exists {};
   struct char_two   : public char_prime {};

   // Tags for zerodivisor_cateogyr

   struct zerodivisor_tag {};
   struct nozerodivisor_tag {};

   // Tags for size_category

   struct finite_tag   { static const bool finite = true;  };
   struct infinite_tag { static const bool finite = false; };

   // Tags for primedetection_category

   struct primedetection_tag {};
   struct noprimedetection_tag {};

   // flags for printShort()

   enum PrintShortFlag
   {
      FIT_FOR_MUL = 1
   };
   
   /**
    *  Addition based on bit operations
    */

   template<typename T>
   class BitOpBasedAddition
   {
   public:
      typedef char_two char_category;

      static bool is0 (const T& a)  { return !a; }

      static T  add (const T& a, const T& b)  { return a ^ b; }
      static T  sub (const T& a, const T& b)  { return a ^ b; }

      static void addTo   (T& a, const T& b)  { a ^= b; }
      static void subFrom (T& a, const T& b)  { a ^= b; }

      static T  neg (const T& a)  { return a; }
      static void negate  (T&)  {}

      static T dbl (const T&)  { return 0; }
      static void times2 (T& a)  { a = 0; }
      static T times (const T& a, unsigned k)  { return odd (k) ? a : 0; }

      static unsigned additiveOrder (const T& a)  { return a ? 2 : 1; }
      static unsigned characteristic()  { return 2; }

   protected:
      BitOpBasedAddition () {}
   };

   /**
    *  destructiveAssign()
    *
    *  Assigns src to dst in the most efficient way, allowing even the
    *  destruction of src.
    *
    *  The default action is normal assignment. However, specializations may be
    *  provided using, e.g. swap() or similar routines.
    */

   template<typename T>
   inline
   void destructiveAssign (T& dst, T& src) { dst = src; }

} // namespace HIntLib

#define HINTLIB_TRIVIAL_DOMAIN_MEMBERS \
   unsigned additiveOrder (const type& a) const \
      { return is0(a) ? 1 : characteristic(); }

#define HINTLIB_TRIVIAL_CYCLIC_MEMBERS \
   HINTLIB_TRIVIAL_DOMAIN_MEMBERS \
   static unsigned extensionDegree()  { return 1; } \
   const type& frobenius (const type& a) const  { return a; } \
   const type& invFrobenius (const type& a) const  { return a; }

#endif

