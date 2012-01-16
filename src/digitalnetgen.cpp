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

#include <HIntLib/digitalnetgen.tcc>

namespace HIntLib
{

/**
 *  Digital Net
 */

DigitalNet::DigitalNet (unsigned base, unsigned _m)
   : m (std::min (
        _m,
        unsigned(logInt (std::numeric_limits<Index>::max(), Index(base))))),
     size (powInt (base, m))
{}


#define HINTLIB_INSTANTIATE(X) \
   template class DigitalNetGen<X, Index>; \
   template void DigitalNetGenNormal<X, Index>::updateX (); \
   template void DigitalNetGenGray<X, Index>::updateX (); \
   template void DigitalNetGenCyclicGray<X, Index>::updateX ();

HINTLIB_INSTANTIATE (OneDimVectorSpace<ModularArithField<unsigned char> >)
HINTLIB_INSTANTIATE (OneDimVectorSpace<LookupField<unsigned char> >)
HINTLIB_INSTANTIATE (OneDimVectorSpace<LookupFieldPow2<unsigned char> >)
HINTLIB_INSTANTIATE (OneDimVectorSpace<LookupFieldPrime<unsigned char> >)
HINTLIB_INSTANTIATE (OneDimVectorSpace<GF2>)
HINTLIB_INSTANTIATE (GF2VectorSpace<u32>)
typedef LookupVectorSpace    <unsigned char,unsigned char> X1;
typedef LookupVectorSpacePow2<unsigned char,unsigned char> X2;
HINTLIB_INSTANTIATE (X1)
HINTLIB_INSTANTIATE (X2)

#undef HINTLIB_INSTANTIATE
}

