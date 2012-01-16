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
 *  Hypercube
 *
 *  A n-dimensional hypercube.
 */

#ifndef HINTLIB_HYPERCUBE_H
#define HINTLIB_HYPERCUBE_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <iosfwd>

#include <HIntLib/mymath.h>
#include <HIntLib/array.h>
#ifdef HINTLIB_PARALLEL
   #include <HIntLib/buffer.h>
#endif


namespace HIntLib
{

class Hypercube
{

// These inline functions are used by other (public) functions
// So they have to be defined here

private:
   real* width  ()            { return &data [dim]; }
   real* center ()            { return &data [0]; }
   real& width  (unsigned i)  { return width()  [i]; }
   real& center (unsigned i)  { return center() [i]; }
   void calcVolume ();

public:

   /*
    *  Normal Constructors
    */

   // [0,1]^dim

   explicit Hypercube (unsigned dim);

   // [a_1,b_1]x..x[a_dim,b_dim]

   Hypercube (unsigned dim, const real a [], const real b []);

   // [a,b]^dim

   Hypercube (unsigned dim, real a, real b);

   // Copy constructor

   Hypercube (const Hypercube &);


   /*
    *  Create a new Hypercube be splitting another one
    */

   Hypercube (Hypercube &, unsigned dim);

   // Assignement

   Hypercube& operator= (const Hypercube&);

   // Comparision

   bool operator== (const Hypercube &) const;
   bool operator!= (const Hypercube &h) const  { return ! (*this == h); }

#ifdef HINTLIB_PARALLEL
   Hypercube (unsigned dim,
              int source, int tag, MPI_Comm comm, MPI_Status *status);
   Hypercube (unsigned dim, RecvBuffer &b);
   int send (int dest, int tag, MPI_Comm comm) const;
   void isend (int dest, int tag, MPI_Comm comm) const;
   int recv (int source, int tag, MPI_Comm comm, MPI_Status *status); 

   MPI_Datatype getMPIDatatype (void) const;

   void initAfterReceive ();

   friend SendBuffer& operator<< (SendBuffer &, const Hypercube &);
   friend RecvBuffer& operator>> (RecvBuffer &, Hypercube &);
#endif

   /*
    *  Get geometry information about a cube
    */

   unsigned getDimension () const       { return dim; }
         real getVolume () const        { return volume; }
   const real* getCenter () const       { return &data [0]; }
   const real* getWidth  () const       { return &data [dim]; }
         real  getCenter (unsigned i) const { return data [i]; }
         real  getWidth  (unsigned i) const { return data [i + dim]; }
         real  getUpperBound (unsigned i) const
                     { return getCenter (i) + getWidth (i); }
         real  getLowerBound (unsigned i) const
                     { return getCenter (i) - getWidth (i); }
         real  getDiameter (unsigned i) const { return 2.0 * getWidth (i); }

   // Places a point can be relative to a Hypercube

   enum Location {INSIDE, OUTSIDE, BORDER};

   /*
    *  Methods for changing a cube
    */

   void set (unsigned dim, real a, real b);

   void move (unsigned dim, real distance);

   void cutLeft  (unsigned dim);
   void cutRight (unsigned dim);

private:
   const unsigned dim;
   Array<real> data;

   real volume;
};


// A number of tests

bool isUnitCube    (const Hypercube &);
bool isPointInside (const Hypercube &, const real[]);
Hypercube::Location
     whereIsPoint  (const Hypercube &, const real[]);

// Printing a hypercube

std::ostream& operator<< (std::ostream & o, const Hypercube &h);

// build union of two cubes

bool unite (Hypercube &h, const Hypercube &hh);

void throwDimensionMismatch(unsigned dim1, unsigned dim2);
inline
void checkDimensionEqual (const Hypercube &h1, const Hypercube &h2)
{
   if (h1.getDimension() != h2.getDimension())
   {
      throwDimensionMismatch (h1.getDimension(), h2.getDimension());
   }
}


/****** Implementation ***************/

inline
void Hypercube::calcVolume (void)
{
   real v = 1.0;

   // Workaround another g++ Optimizer bug. :-(

   #ifdef __GNUG__
   for (volatile unsigned i = 0; i < dim; ++i)
   #else
   for (unsigned i = 0; i < dim; ++i)
   #endif
   {
      v *= getDiameter (i);
   }

   volume = v;
}

inline
void Hypercube::cutLeft (unsigned split)
{
   width  (split) /= 2.0;
   center (split) -= width (split);

   volume /= 2.0;
}

inline
void Hypercube::cutRight (unsigned split)
{
   width  (split) /= 2.0;
   center (split) += width (split);

   volume /= 2.0;
}

inline Hypercube::Hypercube (const Hypercube &h)
   : dim (h.dim), data (h.data, 2 * dim), volume (h.volume)
{}

inline Hypercube::Hypercube (Hypercube &h, unsigned split)
   : dim (h.dim), data (h.data, 2 * dim), volume (h.volume)
{
   h.cutLeft  (split);
     cutRight (split);
}

inline
void Hypercube::set (unsigned i, real a, real b)
{
   center() [i] =     (a + b) * .5;
   width()  [i] = abs (b - a) * .5;
   calcVolume ();
}
                                
inline
void Hypercube::move (unsigned i, real distance)
{
   center() [i] += distance;
}

}  // namespace HIntLib

#endif

