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
 *  The matrices defining the Sobol Digtal Net
 */

#ifndef SOBOL_MATRIX_H
#define SOBOL_MATRIX_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/defaults.h>

#ifdef HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/hintlib_limits.h>
#endif

#include <HIntLib/generatormatrix.h>


namespace HIntLib
{

class SobolMatrix : public GeneratorMatrix2<u64>
{
public:
   SobolMatrix()
      : GeneratorMatrix2<u64> (MAX_DIM, MAX_LOG_N, PRECISION, v_mem[0])  {}

   unsigned getT (unsigned d) const  { return t_s [d]; }

   // We don't have data for creating the v array for dimension larger than 40

   static const unsigned MAX_DIM = 40;
   static const unsigned PRECISION = std::numeric_limits<real>::digits - 1;
   static const unsigned MAX_LOG_N =
      std::numeric_limits<Index>::digits > 48 ? 48 :
      std::numeric_limits<Index>::digits - 1;

private:

   // Vectors used to calculate the sequence

   static const u64 v_mem [MAX_LOG_N][MAX_DIM];
   static const unsigned t_s [MAX_DIM];
};

}  // namespace HIntLib

#endif

