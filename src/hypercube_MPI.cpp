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

#include <HIntLib/hlmpi.h>
#include <HIntLib/hypercube.h>

namespace L = HIntLib;

MPI_Datatype
L::Hypercube::getMPIDatatype () const
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

 
int
L::Hypercube::send (int dest, int tag, MPI_Comm comm) const
{
   const real* p = data;

   return MPI_Send (const_cast<real*> (p),
                    2 * dim, MPIType<real>::type, dest, tag, comm);
}

void
L::Hypercube::isend (int dest, int tag, MPI_Comm comm) const
{
   const real* p = data;

   MPI_Request request;
 
   MPI_Isend (const_cast<real*> (p), 2 * dim, MPIType<real>::type,
              dest, tag, comm, &request);

   MPI_Request_free (&request);
}
 
int
L::Hypercube::recv (int source, int tag, MPI_Comm comm, MPI_Status *status)
{
    const real* p = data;

    int res = MPI_Recv (const_cast<real*> (p),
                    2 * dim, MPIType<real>::type, source, tag, comm, status);
 
    calcVolume();

    return res;
}
 
L::Hypercube::Hypercube (unsigned dim,
                      int source, int tag, MPI_Comm comm, MPI_Status *status)
   : dim (dim), data (2 * dim)
{
   recv (source, tag, comm, status);
}

void
L::Hypercube::initAfterReceive ()
{
   calcVolume();
}

L::SendBuffer &
L::operator<< (SendBuffer &b, const Hypercube &h)
{
   const real* p = h.data;
   b.pack (p, 2 * h.dim, MPIType<real>::type);
   return b;
}

L::RecvBuffer &
L::operator>> (RecvBuffer &b, Hypercube &h)
{
   real* p = h.data;
   b.unpack (p, 2 * h.dim, MPIType<real>::type);
   h.initAfterReceive();
   return b;
}

