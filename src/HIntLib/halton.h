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
 *   Halton
 *
 *   QRNGenerator producing the Halton sequence based on the first dim
 *   prime numbers.
 */

#ifndef HINTLIB_HALTON_H
#define HINTLIB_HALTON_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/qrnsequencebase.h>
#include <HIntLib/shiftscale.h>


namespace HIntLib
{

class Halton : public QRNSequenceBase
{
public:
   Halton (const Hypercube &h);

   static Index getOptimalNumber (Index max, const Hypercube &)  { return max; }

   void first          (real* p, Index _n = 0)  { n = _n-1; next(p); }
   void firstDontScale (real* p, Index _n = 0)  { n = _n-1; nextDontScale(p); }
   void next          (real*);
   void nextDontScale (real*);
   Index getOptimalNumber (Index max) const  { return max; }

private:
   ShiftScale ss;
};

}  // namespace HIntLib

#endif

