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

#ifndef HINTLIB_FAURE_H
#define HINTLIB_FAURE_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/generatormatrix.h>

namespace HIntLib
{
   class Faure
      : public HeapAllocatedGeneratorMatrixGen<unsigned char>
   {
   public:
      static void init (MutableGeneratorMatrixGen<unsigned char> &);
      
      Faure (unsigned _dim);
      Faure (unsigned _dim, unsigned _m, unsigned _prec);

      Faure (unsigned _dim, unsigned _base)
      : HeapAllocatedGeneratorMatrixGen<unsigned char> (_base, 1, _dim)
      { init (*this); }

      Faure (unsigned _dim, unsigned _m, unsigned _prec, unsigned _base)
         : HeapAllocatedGeneratorMatrixGen<unsigned char>
            (_base, 1, _dim, _m, _prec)
      { init (*this); }
   };

}   // namespace HIntLib

#endif


