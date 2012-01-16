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
 *  CubePartitioner
 *
 *  Splits a Hypercube into a given number of (almost) equally sized
 *  sub-Hypercubes and performs an action on each of them
 *
 *  The virtual method action() has to be overwritten to define the action that
 *  should be performed on each subcube.
 *
 *  DistanceType:
 *
 *    Determines the position of the planes used for splitting.
 *
 *    EQUIDISTANT:
 *      The planes are spaced equidistantly along the coordinate axis.
 *
 *    EQUIVOLUME:
 *      The planes are spaced such that all resulting subcubes have identical
 *      volume.
 */

#ifndef HINTLIB_CUBEPARTITIONER_H
#define HINTLIB_CUBEPARTITIONER_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/defaults.h>

namespace HIntLib
{

class Hypercube;

class CubePartitioner
{
public:

   enum DistanceType {EQUIDISTANT, EQUIVOLUME};

   void operator() (
      const Hypercube &, Index numSubcubes, DistanceType = EQUIVOLUME);

   virtual void action (const Hypercube &) = 0;

   virtual ~CubePartitioner () {};
};

}  // namespace HIntLib

#endif

