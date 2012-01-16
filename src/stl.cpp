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

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/defaults.h>

namespace L = HIntLib;

#include <vector>
#include <queue>
#include <algorithm>
#include <iostream>
#include <string>

#ifdef HINTLIB_HAVE_CSTDDEF
  #include <cstddef>
  #define HINTLIB_SN std::
#else
  #include <stddef.h>
  #define HINTLIB_SN
#endif

#include <HIntLib/region.h>

namespace std
{
// fill<>()

template void fill<>(bool*,          bool*,          const bool&);
template void fill<>(char*,          char*,          const char&);
template void fill<>(L::real*,       L::real*,       const L::real&);
template void fill<>(int*,           int*,           const int&);
template void fill<>(unsigned char*, unsigned char*, const unsigned char&);
template void fill<>(unsigned char*, unsigned char*, const int&);
template void fill<>(unsigned short*,unsigned short*,const unsigned short&);
template void fill<>(unsigned short*,unsigned short*,const int&);
template void fill<>(L::u32*,        L::u32*,        const L::u32&);
template void fill<>(L::u32*,        L::u32*,        const int&);
#ifdef HINTLIB_U32_NOT_EQUAL_U64
template void fill<>(L::u64*,        L::u64*,        const L::u64&);
template void fill<>(L::u64*,        L::u64*,        const int&);
#endif
template void fill<>(L::EstErr*,     L::EstErr*,     const L::EstErr&);
template void fill<>(L::Region**,    L::Region**,    L::Region *const&);

// vector<>

template class vector<unsigned char>;
template class vector<unsigned short>;
template class vector<int>;
template class vector<char>;
template class vector<L::Region*>;
template class priority_queue<L::Region*, vector<L::Region*>,
                                 L::RegionErrorLess>;

template L::Region** fill_n<>(L::Region**,HINTLIB_SN size_t,L::Region *const&);
template unsigned char* fill_n<>(unsigned char*,HINTLIB_SN size_t,const unsigned char&);
template unsigned short* fill_n<>(unsigned short*,HINTLIB_SN size_t,const unsigned short&);
template int* fill_n<>(int*,HINTLIB_SN size_t,const int&);

#ifdef HINTLIB_SGI
template void __push_heap<>(L::Region**,int,int,L::Region*,L::RegionErrorLess);

template void __adjust_heap<>(L::Region**,int,int,L::Region*,L::RegionErrorLess);
#endif

// Stream related instantiations

template istream& getline<>(istream&, string&, char);
template ostream& operator<<(ostream&,const string&);

#ifdef HINTLIB_SGI
template istream& istream::_M_get_num<>(unsigned int&);

template ostream& ostream::_M_put_num<>(long);
template ostream& ostream::_M_put_num<>(unsigned long);
#ifdef HINTLIB_HAVE_UNSIGNED_LONG_LONG_INT
template ostream& ostream::_M_put_num<>(unsigned long long);
#endif
template ostream& ostream::_M_put_num<>(const void*);
template ostream& ostream::_M_put_num<>(L::real);
#endif

#undef HINTLIB_SN

}  // namespace std

