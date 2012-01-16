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

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/shiftscale.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/exception.h>

namespace L = HIntLib;

/**
 *  Constructor
 */

L::ShiftScale::ShiftScale (int dim) : t (dim) {}

L::ShiftScale::ShiftScale (const Hypercube &h)
   : t (h.getDimension())
{
   set (h);
}

L::ShiftScale::ShiftScale (const Hypercube &h, real ub)
   : t (h.getDimension())
{
   set (h, ub);
}

L::ShiftScale::ShiftScale (const Hypercube &h, real lb, real ub)
   : t (h.getDimension())
{
   set (h, lb, ub);
}

L::ShiftScale::ShiftScale (const Hypercube &dest, const Hypercube &src)
   : t (dest.getDimension())
{
  set (dest, src);
}


/**
 *  set()
 */

void L::ShiftScale::set (const Hypercube &h)
{
   for (int i = 0; i != h.getDimension(); ++i)
   {
      t[i].set (h.getLowerBound(i), h.getDiameter(i));
   }
}

void L::ShiftScale::set (const Hypercube &h, real ub)
{
   for (int i = 0; i != h.getDimension(); ++i)
   {
      t[i].set (h.getLowerBound(i), h.getDiameter(i) / ub);
   }
}

void L::ShiftScale::set (const Hypercube &h, real lb, real ub)
{
   for (int i = 0; i != h.getDimension(); ++i)
   {
      t[i].set (h.getLowerBound(i) - lb * h.getDiameter(i) / (ub-lb),
                h.getDiameter(i) / (ub-lb));
   }
}

void L::ShiftScale::set (const Hypercube &dest, const Hypercube &src)
{
   checkDimensionEqual (dest.getDimension(), src.getDimension());

   for (int i = 0; i != dest.getDimension(); ++i)
   {
      t[i].set (
           dest.getLowerBound(i)
         - src.getLowerBound(i) * dest.getDiameter(i) / src.getDiameter(i),
         dest.getDiameter(i) / src.getDiameter(i));
   }
}

