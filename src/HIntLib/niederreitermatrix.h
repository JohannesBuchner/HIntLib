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
 *  The matrices defining the Niederreiter Digital Net
 */

#ifndef HINTLIB_NIEDERREITER_MATRIX_H
#define HINTLIB_NIEDERREITER_MATRIX_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/fallback_limits.h>
#endif

#include <HIntLib/generatormatrix2.h>

namespace HIntLib
{

class NiederreiterMatrix : public GeneratorMatrix2<u64>
{
public:
   NiederreiterMatrix() : GeneratorMatrix2<u64> (MAX_DIM)
   {
      setMatrix (v_mem[0]);
   }

   static const unsigned MAX_DIM = 200;

private:

   // Vectors used to calculate the sequence

   static const u64 v_mem [DEFAULT_M_BASE2][MAX_DIM];
   // static const unsigned t_s [MAX_DIM];
};

}  // namespace HIntLib

#endif

