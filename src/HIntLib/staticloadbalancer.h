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
 *  Static Load Balancer
 *
 *  Distributes a task statically among all available processing nodes
 */

#ifndef HINTLIB_STATICLOADBALANCER_H
#define HINTLIB_STATICLOADBALANCER_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/defaults.h>


namespace HIntLib
{

class StaticLoadBalancerBase
{
public:

   enum OP {SUM, MIN, MAX};

   Index getBegin(void) const  { return begin; }
   Index getEnd  (void) const  { return end; }
   Index getRange(void) const  { return end - begin; }
   int   getRank (void) const  { return rank; }
   int   getSize (void) const  { return size; }

protected:

   StaticLoadBalancerBase(void) {};
   StaticLoadBalancerBase(int r, int s, Index b, Index e)
      : rank(r), size(s), begin(b), end(e) {}

   int rank,size;
   Index begin, end;

};

class DummyStaticLoadBalancer : public StaticLoadBalancerBase
{
public:
   DummyStaticLoadBalancer(Index n)
      : StaticLoadBalancerBase (0, 1, 0, n) {}
};

#ifdef HINTLIB_PARALLEL

   class MPIStaticLoadBalancer : public StaticLoadBalancerBase
   {
   public:
      MPIStaticLoadBalancer(Index n, MPI_Comm comm = MPI_COMM_WORLD);
   };

   typedef MPIStaticLoadBalancer StaticLoadBalancer;
#else
   typedef DummyStaticLoadBalancer StaticLoadBalancer;
#endif

}  // namespace HIntLib

#endif
 
