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

/**
 *  Hypercube
 *
 *  A n-dimensional hypercube.
 */

#ifndef HINTLIB_HYPERCUBE_H
#define HINTLIB_HYPERCUBE_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <iosfwd>

#include <HIntLib/hlmath.h>
#include <HIntLib/array.h>
#ifdef HINTLIB_PARALLEL
#  include <HIntLib/buffer.h>
#endif


namespace HIntLib
{

class Hypercube
{

// These inline functions are used by other (public) functions
// So they have to be defined here

private:
   real* width  ()       { return &data [dim]; }
   real* center ()       { return &data [0]; }
   real& width  (int i)  { return width()  [i]; }
   real& center (int i)  { return center() [i]; }
   void calcVolume ();

public:

   /*
    *  Normal Constructors
    */

   // [0,1]^dim

   explicit Hypercube (int dim);

   // [a_1,b_1]x..x[a_dim,b_dim]

   Hypercube (int dim, const real a [], const real b []);

   // [a,b]^dim

   Hypercube (int dim, real a, real b);

   // Copy constructor

   Hypercube (const Hypercube &);


   /*
    *  Create a new Hypercube be splitting another one
    */

   Hypercube (Hypercube &, int dim);

   // Assignement

   Hypercube& operator= (const Hypercube&);

#ifdef HINTLIB_PARALLEL
   Hypercube (int dim, int source, int tag, MPI_Comm comm, MPI_Status *status);
   Hypercube (int dim, RecvBuffer &b);
   int send (int dest, int tag, MPI_Comm comm) const;
   void isend (int dest, int tag, MPI_Comm comm) const;
   int recv (int source, int tag, MPI_Comm comm, MPI_Status *status); 

   MPI_Datatype getMPIDatatype () const;

   void initAfterReceive ();

   friend SendBuffer& operator<< (SendBuffer &, const Hypercube &);
   friend RecvBuffer& operator>> (RecvBuffer &, Hypercube &);
#endif

   /*
    *  Get geometry information about a cube
    */

   int getDimension () const  { return dim; }
         real getVolume () const       { return volume; }
   const real* getCenter () const      { return &data [0]; }
   const real* getWidth  () const      { return &data [dim]; }
         real  getCenter (int i) const { return data [i]; }
         real  getWidth  (int i) const { return data [i + dim]; }
         real  getUpperBound (int i) const
                     { return getCenter (i) + getWidth (i); }
         real  getLowerBound (int i) const
                     { return getCenter (i) - getWidth (i); }
         real  getDiameter (int i) const { return 2.0 * getWidth (i); }

   // Places a point can be relative to a Hypercube

   enum Location {INSIDE, OUTSIDE, BORDER};

   /*
    *  Methods for changing a cube
    */

   void set (int dim, real a, real b);

   void move (int dim, real distance);

   void cutLeft  (int dim);
   void cutRight (int dim);

private:
   const int dim;
   Array<real> data;

   real volume;
};


// Comparing two cubes

bool operator== (const Hypercube &, const Hypercube &);

// Some other tests

bool isUnitCube    (const Hypercube &);
bool isPointInside (const Hypercube &, const real[]);
Hypercube::Location
     whereIsPoint  (const Hypercube &, const real[]);


// Printing a hypercube

std::ostream& operator<< (std::ostream &, const Hypercube &);
#ifdef HINTLIB_BUILD_WCHAR
std::wostream& operator<< (std::wostream &, const Hypercube &);
#endif

// build union of two cubes

bool unite (Hypercube &, const Hypercube &);

void throwDimensionMismatch(int dim1, int dim2);
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
void Hypercube::calcVolume ()
{
   real v = 1.0;
   for (int i = 0; i < dim; ++i)  v *= getDiameter (i);
   volume = v;
}

inline
void Hypercube::cutLeft (int split)
{
   width  (split) /= 2.0;
   center (split) -= width (split);

   volume /= 2.0;
}

inline
void Hypercube::cutRight (int split)
{
   width  (split) /= 2.0;
   center (split) += width (split);

   volume /= 2.0;
}

inline Hypercube::Hypercube (const Hypercube &h)
   : dim (h.dim), data (h.data, 2 * dim), volume (h.volume)
{}

inline Hypercube::Hypercube (Hypercube &h, int split)
   : dim (h.dim), data (h.data, 2 * dim), volume (h.volume)
{
   h.cutLeft  (split);
     cutRight (split);
}

inline
void Hypercube::set (int i, real a, real b)
{
   center() [i] =     (a + b) * .5;
   width()  [i] = abs (b - a) * .5;
   calcVolume ();
}
                                
inline
void Hypercube::move (int i, real distance)
{
   center() [i] += distance;
}

}  // namespace HIntLib

#endif

