/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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
 *  bitop.h
 *
 *  Some functions that perform bit operations on integer-like types
 */

#ifndef HINTLIB_BITOP_TCC
#define HINTLIB_BITOP_TCC 1

#include <algorithm>

#include <HIntLib/bitop.h>

namespace HIntLib
{

/**
 *  thread()
 *
 *  Given integers  ...a3a2a1a0 and ...b3b2b1b0, create ...b3a3b2a2b1a1b0a0
 */

template<class T>
T thread (T a, T b, unsigned na, unsigned nb)
{
   const T amask = (T(1) << na) - 1;
   const T bmask = (T(1) << nb) - 1;
   int shift = 0;

   T result = 0;

   while (a | b)
   {
      result |= (a & amask) << shift;
      shift += na;
      a >>= na;

      result |= (b & bmask) << shift;
      shift += nb;
      b >>= nb;
   }

   return result;
}

/**
 *  unthread()
 *
 *  Given an integer  ...a5a4a3a2a1a0, create ...a4a2a0 and ...a5a3a1.
 */

template<class T>
void unthread (T both, T& a, T& b, unsigned na, unsigned nb)
{
   a = b = 0;
   const T amask = (T(1) << na) - 1;
   const T bmask = (T(1) << nb) - 1;
   int ashift = 0;
   int bshift = 0;

   while (both)
   {
      a |= ((both & amask) << ashift);
      ashift += na;
      both >>= na;

      b |= ((both & bmask) << bshift);
      bshift += nb;
      both >>= nb;
   }
}

/**
 *  threadn()
 */

template<class T>
T threadn (T* indices, unsigned num)
{
   T result = 0;
   T mask = 1;

   for (;;)
   {
      bool work = false;

      for (unsigned i = 0; i < num; ++i)
      {
          if (indices[i] & 1)  result |= mask;
          mask <<= 1;
          if (indices[i] >>= 1)  work = true;
      }

      if (! work)  break;
   }

   return result;
}

/**
 *  unthreadn()
 */

template<class T>
void unthreadn (T all, T* indices, unsigned num)
{
   T mask = 1;

   std::fill (indices, indices + num, T(0));

   while (all)
   {
      for (unsigned i = 0; i < num; ++i)
      {
         if (all & 1)  indices[i] |= mask;
         all >>= 1;
      }

      mask <<= 1;
   }
}

/**
 * threadinf()
 */

template<class T>
T threadinf (T* indices, unsigned indexBound)
{
   // determine bound on depth
   
   T m = 0;
   for (unsigned i = 0; i < indexBound; ++i)
   {
      m = std::max (m, indices[i]);
   }
   const unsigned depthBound = unsigned (ms1 (m) + 1);
   
   T mask = 1;
   T result = 0;
   unsigned index = 0;
   unsigned depth = 0;
   
   while (index < indexBound || depth < depthBound)
   {
      if (index < indexBound && indices[index] & (T(1) << depth))
      {
         result |= mask;
      }

      if (! (mask <<= 1))  break;

      if (index == depth)
      {
         depth = 0;
         ++index;
      }
      else if (depth + 1 == index)
      {
         ++depth;
         index = 0;
      }
      else if (index < depth)
      {
         ++index;
      }
      else  // index > depth
      {
         ++depth;
      }
   }

   return result;
}


/**
 * unthreadinf()
 */

template<class T>
unsigned unthreadinf (T all, T* indices)
{
   unsigned index = 0;
   unsigned depth = 0;

   indices[0] = 0;
   
   while (all)
   {
      if (all & 1)  indices[index] |= (T(1) << depth);
      all >>= 1;

      if (index == depth)
      {
         depth = 0;
         indices[++index] = 0;
      }
      else if (depth + 1 == index)
      {
         ++depth;
         index = 0;
      }
      else if (index < depth)
      {
         ++index;
      }
      else  // index > depth
      {
         ++depth;
      }
   }

   unsigned num = std::max (index, depth) + 1;

   while (num && indices[num - 1] == 0)  --num;

   return num;
}


}  // namespace HIntLib

#endif

