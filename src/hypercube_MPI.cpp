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
#pragma implementation "hypercube.h"
#endif
 
#include <HIntLib/mympi.h>
 
#include "hypercube.cpp"

namespace L = HIntLib;

MPI_Datatype L::Hypercube::getMPIDatatype () const
{
   MPI_Datatype type   = MPIType<real>::type;
   int          length = 2 * dim;
   MPI_Aint     disp;

   const real* p = data;

   MPI_Address (const_cast<real*> (p), &disp);

   MPI_Datatype newType;

   MPI_Type_struct (1, &length, &disp, &type, &newType);

   return newType;
}

L::Hypercube::Hypercube (unsigned _dim, RecvBuffer &b)
   : dim (_dim), data (2 * _dim)
{
   b >> *this;
}

