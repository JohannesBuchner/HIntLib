/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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


#ifndef HINTLIB_REGION_H
#define HINTLIB_REGION_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <functional>
#include <iosfwd>

#include <HIntLib/cubaturerule.h>
#include <HIntLib/hypercube.h>
#ifdef HINTLIB_PARALLEL
#  include <HIntLib/buffer.h>
#endif


namespace HIntLib
{

class Region
{
private:
   void eval (Integrand &f, EmbeddedRule &rule)
      { splitDim = rule.evalError (f, h, ee); }

public:

#ifdef HINTLIB_PARALLEL
   // Create an empty (invalid) Region that can be filled with MPI_Receive

   Region (int dim) : h (dim) {}
   Region (int dim, RecvBuffer &b);
   Region (int dim, int source, int tag, MPI_Comm comm, MPI_Status *status);
#endif

   // Create new Region

   Region (const Hypercube &h, Integrand &f, EmbeddedRule &rule)
      : h (h), numOfSplits (0) { eval (f, rule); }

   // Split an old Region

   Region (Region &, Integrand &, EmbeddedRule &);

   // get functions

   real getEstimate () const              { return ee.getEstimate (); }
   real getError () const                 { return ee.getError (); }
   const EstErr& getEstErr () const       { return ee; }
   const Hypercube& getHypercube () const { return h; }

#ifdef HINTLIB_PARALLEL

   MPI_Datatype getMPIDatatype () const;
   int send (int dest, int tag, MPI_Comm comm) const;
   int isend (int dest, int tag, MPI_Comm comm) const;
   void initAfterReceive ();
#endif

private:

#ifdef HINTLIB_PARALLEL
   int recv (int source, int tag, MPI_Comm comm, MPI_Status *status);
#endif

   Hypercube h;
   EstErr ee;

   int splitDim;
   Index numOfSplits;

#ifdef HINTLIB_PARALLEL
   friend SendBuffer& operator<< (SendBuffer &, const Region &);
   friend RecvBuffer& operator>> (RecvBuffer &, Region &);
#endif
};

std::ostream & operator<< (std::ostream &, const Region &);
#ifdef HINTLIB_BUILD_WCHAR
std::wostream & operator<< (std::wostream &, const Region &);
#endif

#ifdef HINTLIB_PARALLEL
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


#ifdef HINTLIB_PARALLEL

inline
Region::Region (int dim, int source, int tag, MPI_Comm comm, MPI_Status *status)
   : h (dim)
{
   recv (source, tag, comm, status);
}

#endif  // defined HINTLIB_PARALLEL

}  // namespace HIntLib

#endif

