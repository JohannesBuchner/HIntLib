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
 *  FourthDiff
 */

#ifndef HINTLIB_FOURTHDIFF_H
#define HINTLIB_FOURTHDIFF_H 1

#ifdef __GNUG__
#pragma interface
// Implementation in regioncollection.cpp
#endif

#include <algorithm>

#include <HIntLib/array.h>
#include <HIntLib/mymath.h>
#include <HIntLib/integrand.h>
#include <HIntLib/hypercube.h>
#include <HIntLib/minmaxfinder.h>


namespace HIntLib
{

/**
 *  FourthDiff
 *
 *  Determines the direction of the largest fourth divided difference
 *
 *  The 4th divided difference is an approximation of the 4th partial
 *  derivative in the given domain.
 */

class FourthDiff
{
public:
   FourthDiff (unsigned _dim) : dim (_dim), p (_dim) {};

   unsigned getDimension() const  { return dim; }
   unsigned getNumPoints() const  { return 4 * dim + 1; }

   unsigned operator() (Integrand &, const Hypercube &);
   unsigned operator() (Integrand &, const real* c, const real* w);

private:
   unsigned dim;
   Point p;
};

inline
unsigned FourthDiff::operator() (Integrand &f, const real* c, const real* w)
{
   std::copy (&c[0], &c[dim],  &p[0]);

   MaxFinder<real> mf;
   unsigned dimDiffMax = 0;
 
   const real center = f(p);

   for (unsigned d = 0; d < dim; ++d)
   {
      p[d] = c[d] -       w[d]; real f1  = f(p);
      p[d] = c[d] +       w[d];      f1 += f(p);
      p[d] = c[d] - 0.5 * w[d]; real f2  = f(p);
      p[d] = c[d] + 0.5 * w[d];      f2 += f(p);
      p[d] = c[d];
 
      if (mf << abs (f1 - 4.0*f2 + 6.0*center))  dimDiffMax = d;
   }
 
   return dimDiffMax;                                                           
}

inline
unsigned FourthDiff::operator() (Integrand &f, const Hypercube &h)
{
   return operator()(f, h.getCenter(), h.getWidth());
}

}  // namespace HIntLib

#endif

