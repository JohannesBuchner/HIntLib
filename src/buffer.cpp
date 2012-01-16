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
 *  SendBuffer
 *  RecvBuffer
 */

#ifdef __GNUG__
#pragma implementation
#endif

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/hlmpi.h>

#include <HIntLib/buffer.h>

namespace L = HIntLib;


/**
 *  Send Buffer
 */

bool L::SendBuffer::pack (const void* buf, int n, MPI_Datatype t)
{
   if (! full)
   {
      int maxSize;

      MPI_Pack_size (n, t, comm, &maxSize);

      full = pos + maxSize > SIZE;

      if (! full)
      {
         MPI_Pack (const_cast<void*> (buf), n, t, data, SIZE, &pos, comm);
      }
   }

   return !full;
}


/**
 *  Recv Buffer
 */

void L::RecvBuffer::recv (int source, int tag)
{
   MPI_Recv (data, SIZE, MPI_PACKED, source, tag, comm, &status);
   MPI_Get_count (&status, MPI_PACKED, &bytesReceived);
}


/**
 *  bufferSendRecv()
 */

void L::bufferSendRecv (SendBuffer &out, int dest,   int destTag,
                        RecvBuffer &in,  int source, int sourceTag)
{
   MPI_Sendrecv (
      out.data, out.pos, MPI_PACKED, dest,   destTag,
      in.data,  in.SIZE, MPI_PACKED, source, sourceTag,
      in.comm, &in.status);
   MPI_Get_count (&in.status, MPI_PACKED, &in.bytesReceived);
}

