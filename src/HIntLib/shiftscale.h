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

/**
 *  ShiftScale
 *
 *  Shifts and scales a point in space
 */

#ifndef SHIFTSCALE_H
#define SHIFTSCALE_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/hypercube.h>

namespace HIntLib
{

/**
 *  ShiftScale1
 *
 *  A 1-dimensional affine tranformation of the form
 *
 *     y(x) = sc * x + sh
 */

class ShiftScale1
{
public:
   ShiftScale1 (real sh, real sc) : shift(sh), scale(sc) {}

   ShiftScale1 () : shift(0), scale(0) {}

   void set (real sh, real sc)  { shift = sh; scale = sc; }

   real operator() (real x) const  { return shift + scale * x; }

private:
   real shift;
   real scale;
};


/**
 *  ShiftScale
 *
 *  Transforms each point of the cube
 *
 *         [0,1]^s   or   [0,ub]^s   or   [lb,ub]^s   of   h
 *
 *  into a given Hypercube.
 */

class ShiftScale
{
public:
   explicit ShiftScale (unsigned dim);
   
   ShiftScale (const Hypercube &);
   ShiftScale (const Hypercube &, real ub);
   ShiftScale (const Hypercube &, real lb, real ub);
   ShiftScale (const Hypercube &, const Hypercube &);

   void set (const Hypercube &);
   void set (const Hypercube &, real ub);
   void set (const Hypercube &, real lb, real ub);
   void set (const Hypercube &, const Hypercube &);

   const ShiftScale1& operator[] (unsigned i) const  { return t[i]; }
private:
   Array<ShiftScale1> t;
};


}  // namespace HIntLib

#endif

