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
 *  Quasi Random Number Sequence Base
 *
 *  (c) 2000 by Rudolf Schuerer
 *
 *  This is the abstract base class for all Quasi Random Number Generators.
 *
 *  The following interface is provided:
 *
 *  Constructor QRNGenerator (unsigned int dim)
 *     Creates a generator for a given dimension
 *
 *  unsigned int getDimension()
 *     Returns the dimension of the sequence
 *
 *  T getOtimalNumber (T max)
 *     Returns the largest number less or equal to max, for which the sequence
 *     has "good" properties.
 *
 *  void first (real[], Index)
 *     Generates the first point of the sequence, starting with at a given
 *     Index.
 *
 *  void next (real [])
 *     Copies the next quasi random vector into the specified array.
 */
  
#ifndef HINTLIB_QRNSEQUENCEBASE_H
#define HINTLIB_QRNSEQUENCEBASE_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <HIntLib/hypercube.h>

namespace HIntLib
{

class QRNSequenceBase
{
public:
   QRNSequenceBase (const Hypercube &);

   unsigned getDimension() const  { return h.getDimension(); }
   const Hypercube& getRegion() const  { return h; }
   Index getIndex() const  { return n; }   // last generated point

protected:
   const Hypercube h;
   Index n;

private:
   QRNSequenceBase (const QRNSequenceBase &);   // Do not copy
                   // Assignment impossible due to const member
};

/**
 *  Digital Net
 *
 *  Basic things shared by all digital nets
 *
 *  It has a certain size, given by  _base_ ^ _m_, which must fit into Index.
 */

class DigitalNet
{
protected:
   DigitalNet (unsigned base, unsigned m);  // implemented in digitalnetgen.cpp

   const unsigned m;
   const Index size;
public:
   Index getSize() const  { return size; }
   unsigned getM() const  { return m; }
   typedef enum {TRUNCATE, CENTER, FULL} Truncation;
};

}  // namespace HIntLib

#endif

