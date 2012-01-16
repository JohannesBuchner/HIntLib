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
 *  AdaptIntegratorLocalList
 *
 *  Parallel Adaptive Integration with local region collections
 */

#include <memory>

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/hlmpi.h>

#include <HIntLib/adaptintegratorlocallist.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/hlmath.h>
#include <HIntLib/regioncollection.h>
#include <HIntLib/cubaturerule.h>
#include <HIntLib/buffer.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;


L::AdaptIntegratorLocalList::AdaptIntegratorLocalList (
   const EmbeddedRuleFactory *fac,
   MPI_Comm _comm,
   unsigned aFrequency, unsigned maxDimension)
: EmbeddedRuleBasedIntegrator (fac),
  comm (_comm),
  frequency (aFrequency)
{
   // Determine number of nodes and own rank

   MPI_Comm_rank (comm, &rank);
   MPI_Comm_size (comm, &nodes);

   // Determine topology

   if (maxDimension == 1)
   {
      size [0] = nodes;
      dimension = 1;
   }
   else
   {
      for (int i = 1; i * i <= nodes; ++i)
      {
         if (nodes % i == 0)
         {
            size [1] = i;
            size [0] = nodes / i;
         }

         dimension = (size [1] == 1) ? 1 : 2;
      }

      switch (maxDimension)
      {
      case 0:
      case 8:
         if (nodes == 256)
         {
            for (int i = 0; i < 8; ++i)  size [i] = 2;
            dimension = 8;
            break;
         }
      case 7:
         if (nodes == 128)
         {
            for (int i = 0; i < 7; ++i)  size [i] = 2;
            dimension = 7;
            break;
         }
      case 6:
         if (nodes == 64)
         {
            for (int i = 0; i < 6; ++i)  size [i] = 2;
            dimension = 6;
            break;
         }
      case 5:
         if (nodes == 32)
         {
            size [0] = size [1] = size [2] = size [3] = size [4] = 2;
            dimension = 5;
            break;
         }
      case 4:
         if (nodes == 16)
         {
            size [0] = size [1] = size [2] = size [3] = 2;
            dimension = 4;
            break;
         }
         if (nodes == 24)
         {
            size [0] = 3;
            size [1] = size [2] = size [3] = 2;
            dimension = 4;
            break;
         }
      case 3:
         if (nodes == 8)
         {
            size [0] = size [1] = size [2] = 2;
            dimension = 3;
            break;
         }
         if (nodes == 18)
         {
            size [0] = size [1] = 3;
            size [2] = 2;
            dimension = 3;
            break;
         }
         if (nodes == 27)
         {
            size [0] = size [1] = size [2] = 3;
            dimension = 3;
            break;
         }
         if (nodes == 32)
         {
            size [0] = size [1] = 4;
            size [2] = 2;
            dimension = 3;
            break;
         }
      default: ;
      }
   }

   // Determine own position, successor and predesessor

   int x = 1;

   for (unsigned i = 0; i != dimension; ++i)
   {
      coord [i] = (rank % (x * size [i])) / x;

      succ [i] = rank + (coord [i] == size [i] - 1 ? 1 - size [i] :  1) * x;
      pred [i] = rank + (coord [i] == 0            ? size [i] - 1 : -1) * x;

      x *= size [i];
   }

#if 0
   for (unsigned int i = 0; i < dimension; i++)
   {
      cout << rank << "  " << i << "  " << coord [i] << "  " << size [i] << "  " << succ [i] << "  " << pred [i] << endl;
   }
#endif
}


inline
void L::AdaptIntegratorLocalList::doExchange2 (
   RegionCollection &rc, unsigned dir, unsigned dim)
{
   MPI_Status status;

   // Exchange total error and error of worst region

   struct
   {
      real errorTotal;
      real errorWorstRegion;
   } sendBuffer, receiveBuffer;

   sendBuffer.errorTotal       = rc.getError ();
   sendBuffer.errorWorstRegion = rc.getTopError ();

   MPI_Sendrecv (
      &sendBuffer,    2, MPIType<real>::type, succ [dir], ERROR_CMP,
      &receiveBuffer, 2, MPIType<real>::type, pred [dir], ERROR_CMP,
      comm, &status);

   if (sendBuffer.errorTotal == receiveBuffer.errorTotal)  return;

   // Do we have to send or to receive

   if (sendBuffer.errorTotal > receiveBuffer.errorTotal)
   {
      if (sendBuffer.errorTotal - 2 * sendBuffer.errorWorstRegion
          <= receiveBuffer.errorTotal)  return;

      // Create output buffer

      SendBuffer buffer (comm);

      // Pack regions into Buffer

      while (   rc.size () >= 2
             && receiveBuffer.errorTotal + 2 * rc.getTopError ()
              < rc.getError ())
      {
         Region *r = rc.top ();

         if (! (buffer << *r)) break;

         receiveBuffer.errorTotal += r->getError ();

         rc.pop ();
         delete r;
      }

      // Send it

      buffer.send (pred [dir], REGIONS);
   }
   else
   {
      if (receiveBuffer.errorTotal - 2 * receiveBuffer.errorWorstRegion
          <= sendBuffer.errorTotal)  return;

      // Receive data regions and store them in Region Collection

      RecvBuffer buffer (comm, succ [dir], REGIONS);

      while (! buffer.empty ())  rc.push (new Region (dim, buffer));
   }
}


inline
void L::AdaptIntegratorLocalList::doExchange (
   RegionCollection &rc, unsigned dir, unsigned dim)
{
   // Exchange total error

   real otherError;

   {
      MPI_Status status;
      real error = rc.getError ();

      MPI_Sendrecv (
         &error,      1, MPIType<real>::type, succ [dir], ERROR_CMP,
         &otherError, 1, MPIType<real>::type, pred [dir], ERROR_CMP,
         comm, &status);
   }

   // Create output buffer and pack regions into it

   SendBuffer outBuffer (comm);

   while (   rc.size () >= 2
          && otherError + 2 * rc.getTopError () < rc.getError ())
   {
      Region *r = rc.top ();

      if (! (outBuffer << *r))  break;

      otherError += r->getError ();

      rc.pop ();
      delete r;
   }

   // Exchange Regions

   RecvBuffer inBuffer (comm);

   bufferSendRecv (outBuffer, pred [dir], REGIONS,
                   inBuffer,  succ [dir], REGIONS);

   // Unpack all received regions and push onto RegionCollection

   while (! inBuffer.empty ())  rc.push (new Region (dim, inBuffer));
}


void L::AdaptIntegratorLocalList::doExchangeOld (
   RegionCollection &rc, unsigned dir, unsigned /* dim */)
{
   MPI_Status status;

   Region *p = rc.top ();

   Region old (*p);

   MPI_Datatype oldType = old.getMPIDatatype ();
   MPI_Datatype newType = p-> getMPIDatatype ();

   MPI_Type_commit (&oldType);
   MPI_Type_commit (&newType);

   MPI_Sendrecv (
      MPI_BOTTOM, 1, oldType, succ [dir], REGIONS,
      MPI_BOTTOM, 1, newType, pred [dir], REGIONS,
      comm, &status);

   p->initAfterReceive ();

   MPI_Type_free (&oldType);
   MPI_Type_free (&newType);

   rc.push (p);
}


L::Integrator::Status L::AdaptIntegratorLocalList::integrate (
   Integrand &f, const Hypercube &h, Index maxEvaluations,
   real reqAbsError, real reqRelError, EstErr &finalEE)
{
   checkDimension(h, f);
   checkTerminationCriteria (maxEvaluations, reqAbsError, reqRelError, true);

   std::auto_ptr<EmbeddedRule> rule (getRule (h.getDimension()));

   const Index rulePoints = rule->getNumPoints ();

   // Subtract evaluations for initial cube splitting

   const Index initialPoints =
        (1 + 4 * h.getDimension()) * nodes * (ms1 (nodes - 1) + 1)
      + rulePoints * nodes;

   if (initialPoints > maxEvaluations)
   {
#     ifdef HINTLIB_NO_EXCEPTIONS
         finalEE.set (0.0, 0.0);
         return ERROR;
#     else
         throw NoEvaluationsPossible(maxEvaluations);
#     endif
   }

   // Number of possible iterations

   const Index maxIter = (maxEvaluations - initialPoints)
                       / (rulePoints * nodes * 2);

   // Determine interval

   unsigned interval = 1 + frequency / rulePoints;

   // Create Region Collection and push initial hypercube

   RegionCollection rc;

   storeSubcube (rc, h, nodes, rank, f, *rule);

   unsigned dir = 0;

   for (Index i = 0; i != maxIter; ++i)
   {
      if (i > 0 && (i % (2 * ms1 (i+1)) == 0 || i % interval == 0))
      {
         dir = (dir + 1) % dimension;

         if (size [dir] > 1)
         {
            if (size [dir] == 2)
               doExchange2 (rc, dir, h.getDimension ());
            else if (size [dir] != 1)
               doExchange (rc, dir, h.getDimension ());
         }
      }

      // Refine a region

      rc.refine (f, *rule);
   }

   rc.killQueueAndUpdateResult ();

   // Collect results from all nodes

   EstErr ee = rc.getEstErr ();

   MPI_Reduce (&ee, &finalEE, 2, MPIType<real>::type, MPI_SUM, 0, comm);

   // Only node 0 returns a result

   if (rank != 0)  return WRONG_NODE;

   // Determine if we have met some accuracy requirement and return status

   Status status =
         checkRequestedError (finalEE, reqAbsError, reqRelError);

   return (status == ERROR) ? MAX_EVAL_REACHED : status;
}

