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


#ifndef REGION_H
#define REGION_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <functional>
#include <iosfwd>

#include <HIntLib/esterr.h>
#include <HIntLib/embeddedrule.h>
#include <HIntLib/hypercube.h>
#ifdef PARALLEL
   #include <HIntLib/buffer.h>
#endif


namespace HIntLib
{

class Region
{
private:
   void eval (Function &f, EmbeddedRule &rule)
      { splitDim = rule.evalError (f, h, ee); }

public:

#ifdef PARALLEL
   // Create an empty (invalid) Region that can be filled with MPI_Receive

   Region (unsigned dim) : h (dim) {}
   Region (unsigned dim, RecvBuffer &b);
   Region (unsigned dim,
           int source, int tag, MPI_Comm comm, MPI_Status *status);
#endif

   // Create new Region

   Region (const Hypercube &h, Function &f, EmbeddedRule &rule)
      : h (h), numOfSplits (0) { eval (f, rule); }

   // Split an old Region

   Region (Region &, Function &, EmbeddedRule &);

   // get functions

   real getEstimate () const              { return ee.getEstimate (); }
   real getError () const                 { return ee.getError (); }
   const EstErr& getEstErr () const       { return ee; }
   const Hypercube& getHypercube () const { return h; }

#ifdef PARALLEL

   MPI_Datatype getMPIDatatype () const;
   int send (int dest, int tag, MPI_Comm comm) const;
   int isend (int dest, int tag, MPI_Comm comm) const;
   void initAfterReceive ();
#endif

private:

#ifdef PARALLEL
   int recv (int source, int tag, MPI_Comm comm, MPI_Status *status);
#endif

   Hypercube h;
   EstErr ee;

   unsigned splitDim;
   unsigned numOfSplits;

#ifdef PARALLEL
   friend SendBuffer& operator<< (SendBuffer &, const Region &);
   friend RecvBuffer& operator>> (RecvBuffer &, Region &);
#endif
};

std::ostream & operator<< (std::ostream &, const Region &);

#ifdef PARALLEL
SendBuffer& operator<< (SendBuffer &, const Region &);
RecvBuffer& operator>> (RecvBuffer &, Region &);
#endif

struct RegionErrorLess : public std::binary_function<Region*,Region*,bool>
{
   bool operator() (const Region *r1, const Region *r2)
   {
      return r1->getError () < r2->getError (); 
   }
};


/********** Implementation ****************/


#ifdef PARALLEL

inline
Region::Region (unsigned dim,
                int source, int tag, MPI_Comm comm, MPI_Status *status)
   : h (dim)
{
   recv (source, tag, comm, status);
}

#endif  // defined PARALLEL

}  // namespace HIntLib

#endif

