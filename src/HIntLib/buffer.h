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
 *  Buffer
 *
 *  Send- and Receive buffer for MPI
 */

#ifndef HINTLIB_BUFFER_H
#define HINTLIB_BUFFER_H 1

#ifndef HINTLIB_PARALLEL
#error "buffer.h can only be used in PARALLEL mode"
#endif

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/defaults.h>


namespace HIntLib
{

/**
 *  PackBuffer
 *
 *  Common base class for ReceiveBuffer and SendBuffer
 */

class PackBuffer
{
public:
   static const int SIZE = 128000;
protected:
   PackBuffer (MPI_Comm comm) : comm (comm), pos (0) {}

   MPI_Comm comm;
   int pos;
 
   char data [SIZE];
};


class RecvBuffer;


/**
 *  SendBuffer
 *
 *  A buffer for sending arbitrary data
 */

class SendBuffer : public PackBuffer
{
public:
   SendBuffer (MPI_Comm comm) : PackBuffer (comm), full (false) {}

   bool pack (const void* buf, int n, MPI_Datatype t);
   template <class T>
      SendBuffer& operator<< (T x)
         { pack (&x, 1, MPIType<T>::type); return *this; }

   operator bool ()  { return ! full; }

   int send (int dest, int tag);

friend void bufferSendRecv (SendBuffer &out, int dest,   int destTag,
                            RecvBuffer &in,  int source, int sourceTag);
private:
   bool full;
};

inline
int SendBuffer::send (int dest, int tag)
{
   return MPI_Send (data, pos, MPI_PACKED, dest, tag, comm);
}


/**
 *  ReceiveBuffer
 */

class RecvBuffer : public PackBuffer
{
public:
   RecvBuffer (MPI_Comm comm) : PackBuffer (comm), bytesReceived (0) {}
   RecvBuffer (MPI_Comm comm, int source, int tag);

   void recv (int source, int tag);

   void unpack (void* buf, int n, MPI_Datatype t);
   template <class T>
      RecvBuffer& operator>> (T &x)
         { unpack (&x, 1, MPIType<T>::type); return *this; }

   int numBytesLeft () const  { return bytesReceived - pos; }
   bool empty () const        { return numBytesLeft() <= 0; }
   operator bool () const     { return ! empty(); }

   const MPI_Status* getStatus () { return &status; }

friend void bufferSendRecv (SendBuffer &out, int dest,   int destTag,
                            RecvBuffer &in,  int source, int sourceTag);

private:
   int bytesReceived;
   MPI_Status status;
};

inline
RecvBuffer::RecvBuffer (MPI_Comm comm, int source, int tag)
   : PackBuffer (comm)
{
   recv (source, tag);
}

inline
void RecvBuffer::unpack (void* buf, int n, MPI_Datatype t)
{
   MPI_Unpack (data, SIZE, &pos, buf, n, t, comm);
}


}  // namespace

#endif

