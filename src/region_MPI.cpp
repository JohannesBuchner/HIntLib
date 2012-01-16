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

#if defined __GNUG__
#pragma implementation "region.h"
#pragma implementation "esterr.h"
#endif

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/hlmpi.h>

#include "region.cpp"

namespace L = HIntLib;

MPI_Datatype L::EstErr::getMPIDatatype () const
{
   MPI_Datatype type   = MPIType<real>::type;
   int          length = 2;
   MPI_Aint     disp;

   EstErr *p = const_cast<EstErr*> (this);

   MPI_Address (p, &disp);

   MPI_Datatype newType;

   MPI_Type_struct (1, &length, &disp, &type, &newType);

   return newType;
}

#ifdef __GNUG__
inline
#endif
void L::Region::initAfterReceive ()
{
   ee.initAfterReceive ();
   h.initAfterReceive ();
}

MPI_Datatype L::Region::getMPIDatatype () const
{
   MPI_Datatype types  [4];
   int          length [4] = { 1, 1, 1, 1 };
   MPI_Aint     disp   [4];

   Region *p = const_cast<Region*> (this);

   types [0] = MPI_UNSIGNED; MPI_Address (&p->splitDim,    disp + 0);
   types [1] = MPI_UNSIGNED; MPI_Address (&p->numOfSplits, disp + 1);

   types [2] = h.getMPIDatatype ();
   MPI_Address (MPI_BOTTOM, disp + 2);

   types [3] = ee.getMPIDatatype ();
   MPI_Address (MPI_BOTTOM, disp + 3);

   MPI_Datatype newType;

   MPI_Type_struct (4, length, disp, types, &newType);

   MPI_Type_free (&types [2]);
   MPI_Type_free (&types [3]);

   return newType;
}

L::Region::Region (unsigned dim, RecvBuffer &b)
   : h (dim)
{
   MPI_Datatype t = getMPIDatatype ();

   MPI_Type_commit (&t);

   b.unpack (MPI_BOTTOM, 1, t);

   initAfterReceive ();

   MPI_Type_free (&t);
}

int L::Region::recv (int source, int tag, MPI_Comm comm, MPI_Status *status)
{
   MPI_Datatype t = getMPIDatatype ();

   MPI_Type_commit (&t);

   int res = MPI_Recv (MPI_BOTTOM, 1, t, source, tag, comm, status);

   initAfterReceive ();

   MPI_Type_free (&t);

   return res;
}

L::SendBuffer& L::operator<< (SendBuffer &b, const Region &r)
{
   return b << r.splitDim << r.numOfSplits << r.h << r.ee;
}

L::RecvBuffer& L::operator>> (RecvBuffer &b, Region &r)
{
   return b >> r.splitDim >> r.numOfSplits >> r.h >> r.ee;
}

int L::Region::send (int dest, int tag, MPI_Comm comm) const
{
   MPI_Datatype t = getMPIDatatype ();

   MPI_Type_commit (&t);

   int res = MPI_Send (MPI_BOTTOM, 1, t, dest, tag, comm);

   MPI_Type_free (&t);

   return res;
}

int L::Region::isend (int dest, int tag, MPI_Comm comm) const
{
   MPI_Request request;
   MPI_Datatype t = getMPIDatatype ();

   MPI_Type_commit (&t);

   int res = MPI_Isend (MPI_BOTTOM, 1, t, dest, tag, comm, &request);

   MPI_Type_free (&t);
   MPI_Request_free (&request);

   return res;
}


