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

#ifndef HINTLIB_ONE_DIM_VECTOR_SPACE_H
#define HINTLIB_ONE_DIM_VECTOR_SPACE_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/algebra.h>
#include <HIntLib/output.h>

namespace HIntLib
{
   template<typename A>
   class OneDimVectorSpace : private A
   {
   public:
      typedef vectorspace_tag algebra_category;
      typedef typename A::size_category size_category;

      typedef A scalar_algebra;
      typedef typename A::type type;
      typedef type scalar_type;
      typedef type& scalar_reference;

      OneDimVectorSpace (const A& a) : A (a) {}

      const A& getScalarAlgebra() const  { return *this; }

      unsigned dimension() const  { return 1; }

      using A::size;
      using A::element;
      using A::index;
      using A::printSuffix;
      using A::is0;
      using A::add;
      using A::addTo;
      using A::sub;
      using A::subFrom;
      using A::neg;
      using A::negate;
      using A::dbl;
      using A::times2;
      using A::times;
      using A::mul;
      using A::additiveOrder;

      void printShort (std::ostream&,   const type&) const;
      void printShort (std::ostream& o, const type& x, PrintShortFlag) const
         { printShort (o, x); }
      void print      (std::ostream&, const type&) const;
#ifdef HINTLIB_BUILD_WCHAR
      void printShort (std::wostream&,   const type&) const;
      void printShort (std::wostream& o, const type& x, PrintShortFlag) const
         { printShort (o, x); }
      void print      (std::wostream&, const type&) const;
#endif

      scalar_type coord (const type& x, unsigned) const  { return x; }
      scalar_reference coord  (type& x, unsigned) const  { return x; }

      template<typename I> void toCoord (const type& x, I p) const  {  *p = x; }
      template<typename I> void fromCoord (type& x, I p) const  { x = *p; }

      type& scale (type& x, const scalar_type& a) const { return mulBy (x, a); }
   };

   template<typename A>
   void OneDimVectorSpace<A>::print (std::ostream &o, const type& x) const
   {
      Private::Printer ss (o);
      ss << '(';
      A::print (ss, x);
      ss << ')';
   }
#ifdef HINTLIB_BUILD_WCHAR
   template<typename A>
   void OneDimVectorSpace<A>::print (std::wostream &o, const type& x) const
   {
      Private::WPrinter ss (o);
      ss << L'(';
      A::print (ss, x);
      ss << L')';
   }
#endif

   template<typename A>
   void OneDimVectorSpace<A>::printShort (std::ostream &o, const type& x) const
   {
      Private::Printer ss (o);
      ss << '(';
      A::printShort (ss, x);
      ss << ')';
   }
#ifdef HINTLIB_BUILD_WCHAR
   template<typename A>
   void OneDimVectorSpace<A>::printShort (std::wostream &o, const type& x) const
   {
      Private::WPrinter ss (o);
      ss << L'(';
      A::printShort (ss, x);
      ss << L')';
   }
#endif

   template<typename A>
   std::ostream& operator<< (std::ostream &o, const OneDimVectorSpace<A> &v)
   {
      Private::Printer ss (o);
      ss << v.getScalarAlgebra();
      ss.power(1);
      return o;
   } 
#ifdef HINTLIB_BUILD_WCHAR
   template<typename A>
   std::wostream& operator<< (std::wostream &o, const OneDimVectorSpace<A> &v)
   {
      Private::WPrinter ss (o);
      ss << v.getScalarAlgebra();
      ss.power(1);
      return o;
   } 
#endif

} // namespace HIntLib

#endif

