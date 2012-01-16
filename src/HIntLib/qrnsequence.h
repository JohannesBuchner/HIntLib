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
 *  Quasi Random Number Sequence
 *
 *  This is the abstract base class for Quasi Random Number sequences.
 */
  
#ifndef HINTLIB_QRNSEQUENCE_H
#define HINTLIB_QRNSEQUENCE_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

namespace HIntLib
{

class Hypercube;

class QRNSequence
{
public:

   virtual ~QRNSequence () {}

   virtual unsigned getDimension() const = 0;
   virtual const Hypercube& getRegion() const = 0;
   virtual Index getIndex() const = 0;

   virtual void first          (real*, Index = 0) = 0;
   virtual void firstDontScale (real*, Index = 0) = 0;
   virtual void next          (real*) = 0;
   virtual void nextDontScale (real*) = 0;
   virtual Index getOptimalNumber (Index max) const = 0;
};

template<class Seq>
class QRNSequenceT : public QRNSequence
{
private:
   Seq seq;

public:
   QRNSequenceT (const Hypercube &h) : seq (h) {}
   ~QRNSequenceT () {}

   unsigned getDimension() const      { return seq.getDimension(); }
   const Hypercube& getRegion() const { return seq.getRegion(); }
   Index getIndex() const             { return seq.getIndex(); }
   void first          (real *p, Index n)  { seq.first          (p, n); }
   void firstDontScale (real *p, Index n)  { seq.firstDontScale (p, n); }
   void next          (real *p)  { seq.next (p); }
   void nextDontScale (real *p)  { seq.nextDontScale (p); }
   Index getOptimalNumber (Index max) const
      { return seq.getOptimalNumber(max); }
};

template<class Seq>
class QRNSequenceP : public QRNSequence
{
private:
   Seq* seq;

public:
   QRNSequenceP (Seq* s) : seq (s) {}
   ~QRNSequenceP () { delete seq; }

   unsigned getDimension() const      { return seq->getDimension(); }
   const Hypercube& getRegion() const { return seq->getRegion(); }
   Index getIndex() const             { return seq->getIndex(); }
   void first          (real *p, Index n)  { seq->first          (p, n); }
   void firstDontScale (real *p, Index n)  { seq->firstDontScale (p, n); }
   void next          (real *p)  { seq->next (p); }
   void nextDontScale (real *p)  { seq->nextDontScale (p); }
   Index getOptimalNumber (Index max) const
      { return seq->getOptimalNumber(max); }
};

}  // namespace HIntLib

#endif

