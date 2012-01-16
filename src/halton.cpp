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
#pragma implementation "qmcroutinesimp.h"
#endif

#include <HIntLib/halton.h>

#include <HIntLib/qmcroutinesimp.h>

#include <HIntLib/hlmath.h>
#include <HIntLib/prime.h>

namespace L = HIntLib;

using L::real;
using L::Index;

/**
 *  Halton::Halton()
 */

L::Halton::Halton (const Hypercube &h)
   : QRNSequenceBase (h), ss(h)
{}

/**
 *  Create next point
 */

void L::Halton::next (real points[])
{
   ++n;

   for (unsigned i = getDimension(); i; )
   {
      --i;

      points [i] = ss[i] (radicalInverseFunction (n, Prime::nth(i)));
   }
}

void L::Halton::nextDontScale (real points[])
{
   ++n;

   for (unsigned i = getDimension(); i; )
   {
      --i;

      points [i] = radicalInverseFunction (n, Prime::nth(i));
   }
}


