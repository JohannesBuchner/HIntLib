/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schürer <rudolf.schuerer@sbg.ac.at>
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
 *  Array
 *
 *  Wrapper Class for C-style arrays
 *
 *  This class introduces NO overhead compared to C style arrays
 *  However, it guarantees proper cleanup (memory deallocation), even in the
 *  case of exceptions.
 */

#ifndef HINTLIB_ARRAY_H
#define HINTLIB_ARRAY_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_CSTDDEF
  #include <cstddef>
  #define HINTLIB_SN std::
#else
  #include <stddef.h>
  #define HINTLIB_SN
#endif

#include <algorithm>


namespace HIntLib
{

template<class T, int offset = 0>
class Array
{
public:

   // Types

   typedef T value_type;
   typedef HINTLIB_SN size_t size_type;
   typedef HINTLIB_SN ptrdiff_t differencey_type;

   typedef T* iterator;
   typedef const T* const_iterator;

   typedef T* pointer;
   typedef const T* const_pointer;
   typedef T& reference;
   typedef const T& const_reference;

   // Iterators

   iterator begin ()             { return array; }
   const_iterator begin () const { return array; }

   // Element access

         reference operator[] (int n)       { return array [n-offset]; }
   const_reference operator[] (int n) const { return array [n-offset]; }

         reference front()       { return *array; }
   const_reference front() const { return *array; }

   operator pointer()             { return array; }
   operator const_pointer() const { return array; }
   operator const_pointer()       { return array; }

   // Constructor

   explicit Array (size_type n) { array = new T [n]; }

   Array (size_type n, const value_type& val);
   Array (const_pointer a, size_type n);

   void assign (const_pointer a, size_type n);

   // Destructor

   ~Array ()  { delete[] array; }

private:
   T* array;

   // Assignment and copying is impossible because we don't have any size
   // information

   Array (const Array&);
   Array& operator= (const Array&);
};

typedef Array<real> Point;

} // namespace HIntLib


template<class T, int offset>
inline
HIntLib::Array<T, offset>::Array (size_type n, const value_type& val)
{
   array = new T [n];
   std::fill (array, array+n, val);
}

template<class T, int offset>
inline
HIntLib::Array<T, offset>::Array (const_pointer a, size_type n)
{
   array = new T [n];
   std::copy (a, a+n, array);
}

template<class T, int offset>
inline
void HIntLib::Array<T, offset>::assign (const_pointer a, size_type n)
{
   const T* oldArray = array;

   array = new T [n];
   std::copy (a, a+n, array);

   delete[] oldArray;
}                                                                               

#undef HINTLIB_SN

#endif


