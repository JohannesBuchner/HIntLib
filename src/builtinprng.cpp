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
 *  BuiltIn PRNG
 *
 *  Makes the built-in PRNG available
 */

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/builtinprng.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/exception.h>


namespace L = HIntLib;
using L::real;

#if HINTLIB_STATIC_WORKS == 1
   const real L::BuiltInPRNG::RANGE = real(RAND_MAX) + real(1);
   const real L::BuiltInPRNG::RESOLUTION = real(1) / RANGE;
#else
   real L::BuiltInPRNG::RANGE;
   real L::BuiltInPRNG::RESOLUTION;
#endif

bool L::BuiltInPRNG::inUse = false;

L::BuiltInPRNG::BuiltInPRNG (unsigned start)
{
#if HINTLIB_STATIC_WORKS == 0
   RANGE = real(RAND_MAX) + real(1);
   RESOLUTION = real(1) / RANGE;
#endif

   if (inUse)  throw BuiltInPRNGUsedTwice();

   inUse = true;

   init(start);
}

