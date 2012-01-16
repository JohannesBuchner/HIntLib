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

#ifndef HINTLIB_GCD_H
#define HINTLIB_GCD_H 1

#ifdef __GNUG__
#pragma interface
#endif

namespace HIntLib
{

/**
 *  powerMod()
 */
   
template<class A>
typename A::type
powerMod (const A&, const typename A::type&, unsigned,
                    const typename A::type&);

template<class A>
typename A::type
powerMod
   (const A&, const typename A::type&, unsigned, unsigned,
              const typename A::type&);


/**
 *  genGcd()
 *
 *  Calculates the greatest common divisor
 *
 *  A has to be an Euclidean Ring
 */

template<class A>
typename A::type
genGcd (const A&,
   const typename A::type&  u, const typename A::type&  v,
         typename A::type& mu,       typename A::type& mv);
 
template<class A>
typename A::type
genGcd (const A&,
   const typename A::type&  u, const typename A::type& v,
         typename A::type& mu);

template<class A>
typename A::type
genGcd (const A&, const typename A::type& u, const typename A::type& v);


/**
 *  genLcm ()
 */

template<class A>
inline
typename A::type
genLcm (const A&a, const typename A::type& u, const typename A::type& v)
{
   return a.mul (a.div (u, genGcd (a, u, v)), v);
}


/**
 *  genIsCoprime ()
 */

template<class A>
bool
genIsCoprime (const A&, const typename A::type&, const typename A::type&);

}  // namespace HIntLib

#endif
