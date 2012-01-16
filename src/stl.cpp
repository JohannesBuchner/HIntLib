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

#include <vector>
#include <queue>
#include <algorithm>

#include <HIntLib/defaults.h>

#ifdef HINTLIB_SGI
#include <iostream>
#endif

#include <HIntLib/region.h>

namespace L = HIntLib;

template class std::vector<L::Region*>;
template class std::priority_queue<L::Region*, std::vector<L::Region*>,
                                 L::RegionErrorLess>;
template L::Region** std::fill_n<>(L::Region**,unsigned,L::Region *const&);

template void std::fill<>(bool*,          bool*,          const bool&);
template void std::fill<>(double*,        double*,        const double&);
template void std::fill<>(int*,           int*,           const int&);
template void std::fill<>(unsigned char*, unsigned char*, const unsigned char&);
template void std::fill<>(unsigned char*, unsigned char*, const int&);
template void std::fill<>(unsigned short*,unsigned short*,const int&);
#ifdef HINTLIB_32BIT
template void std::fill<>(L::u32*,        L::u32*,        const L::u32&);
template void std::fill<>(L::u32*,        L::u32*,        const int&);
#endif
template void std::fill<>(L::u64*,        L::u64*,        const L::u64&);
template void std::fill<>(L::EstErr*,     L::EstErr*,     const L::EstErr&);
template void std::fill<>(L::Region**,    L::Region**,    L::Region *const&);


#ifdef HINTLIB_SGI
template class std::vector<unsigned char>;
template class std::vector<unsigned short>;
template class std::vector<int>;
template class std::vector<char>;

template void
   std::__push_heap<>(L::Region**,int,int,L::Region*,L::RegionErrorLess);

template void
   std::__adjust_heap<>(L::Region**,int,int,L::Region*,L::RegionErrorLess);

template std::istream& std::getline<>(std::istream&, std::string&, char);
template std::istream& std::istream::_M_get_num<>(unsigned int&);
template std::ostream& std::ostream::_M_put_num<>(L::u32);
template std::ostream& std::ostream::_M_put_num<>(long);
template std::ostream& std::ostream::_M_put_num<>(const void*);
template std::ostream& std::ostream::_M_put_num<>(double);

#endif


