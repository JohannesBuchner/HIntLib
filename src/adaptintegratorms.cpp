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
 *  AdaptIntegratorMS
 *
 *  Parallel Adaptive Integration using a Master-Slave approach
 */

#ifdef __GNUG__
#pragma implementation
#endif

#include <memory>

#include <HIntLib/mympi.h>

#include <HIntLib/adaptintegratorms.h>

#include <HIntLib/cubaturerule.h>
#include <HIntLib/regioncollection.h>
#include <HIntLib/buffer.h>
#include <HIntLib/array.h>
#include <HIntLib/exception_MPI.h>
#include <HIntLib/bitop.h>

// #include "perf.h"

namespace L = HIntLib;
using L::Integrator;

enum TAGS { RESULT, CUBE, DONE };

L::AdaptIntegratorMS::AdaptIntegratorMS
   (const EmbeddedRuleFactory *fac, MPI_Comm _comm)
: EmbeddedRuleBasedIntegrator (fac), frequency (1000), comm (_comm)
{
   MPI_Comm_rank (comm, &rank);
   MPI_Comm_size (comm, &nodes);
}

inline
unsigned L::AdaptIntegratorMS::numWorkers() const
{
   return nodes - 1; 
}

void L::AdaptIntegratorMS::sendTerminationSignal ()
{
   for (unsigned w = 1; w <= numWorkers(); ++w)
   {
      // PA::send (PA::SS, 0, i);

      int dummy;
      MPI_Request request;
 
      MPI_Isend (&dummy, 0, MPI_INT, w, DONE, comm, &request);
      MPI_Request_free (&request);

      // PA::send (PA::SE, 0, i);
   }
}

/**
 *  sendTopRegionToSlave()
 *  sendTopRegionsToAllSlaves()
 */

inline
void L::AdaptIntegratorMS::sendTopRegionToSlave (
      RegionCollection &rc, int w, EstErr* localEE,
      RegionCollection::size_type n)
{
   // PA::send (PA::SS, 0, w);
   SendBuffer buffer (comm);

   for (RegionCollection::size_type i = 0; i < n; ++i)
   {
      Region *r = rc.top();

      if (! (buffer << *r))  break;

      localEE [w] += r->getEstErr();

      rc.pop(); 
      delete r;
   }
 
   buffer.send (w, CUBE);
   // PA::send (PA::SE, 0, w);
}

void L::AdaptIntegratorMS::sendTopRegionToAllSlaves (
       RegionCollection &rc, EstErr* localEE, RegionCollection::size_type n)
{
   Array<Region*> regions (n * numWorkers());

   for (unsigned j = 0; j < n * numWorkers(); ++j)  regions[j] = rc.pop();

   for (unsigned w = 0; w < numWorkers(); ++w)
   {
      // PA::send (PA::SS, 0, w);
      SendBuffer buffer (comm);

      for (unsigned i = 0; i < n; ++i)
      {
         Region *r = regions [i * numWorkers() + w];

         if (buffer << *r)
         {
            localEE [w] += r->getEstErr();

            delete r;
         }
         else
         {
            rc.push (r);
         }
      }
      buffer.send (w+1, CUBE);
      // PA::send (PA::SE, 0, w);
   }
}


/**
 *  recvRegionsFromSlave()
 *  recvRegionsFromAllSlaves()
 */

inline
int L::AdaptIntegratorMS::recvRegionsFromSlave (
   RegionCollection &rc, int workerMask, unsigned dim, EstErr* localEE)
{
   // PA::send (PA::RS, workerMask, 0);
   RecvBuffer buffer (comm, workerMask, RESULT);
   // PA::send (PA::RE, workerMask, 0);
 
   // PA::calc (PA::CS, 0);

   const int w = buffer.getStatus ()->MPI_SOURCE;
 
   buffer >> localEE[w];
 
   while (buffer)
   {
      Region *r = new Region (dim, buffer);
 
      localEE[w] -= r->getEstErr();
 
      rc.push (r);
   }

   // PA::calc (PA::CE, 0);

   return w;
}

void L::AdaptIntegratorMS::recvRegionsFromAllSlaves (
   RegionCollection &rc, unsigned dim, EstErr* localEE)
{
   for (unsigned w = 1; w <= numWorkers(); ++w)
   {
      recvRegionsFromSlave (rc, MPI_ANY_SOURCE, dim, localEE);
   }
}

/**
 *  Algorithm for MASTER in ASYNCHRONOUS mode
 */

Integrator::Status L::AdaptIntegratorMSAsync::master (
   const Hypercube &h,
   real reqAbsError, real reqRelError, EstErr &ee, Index maxIter)
{
   const unsigned dim = h.getDimension ();

   unsigned numActiveWorkers = numWorkers ();

   Array<EstErr> localEE (nodes, EstErr (0.0, 0.0));

   RegionCollection rc;

   bool shutDown = false;
   unsigned initialResultsMissing = numWorkers();
   Array<bool> initialMissing (numWorkers (), true);

   maxIter -= numWorkers();   // # of times we can send new Regions

   // Main loop

   Status status = ERROR;

   while (numActiveWorkers)
   {
      // Get result from slave

      const int worker =
         recvRegionsFromSlave (rc, MPI_ANY_SOURCE, dim, localEE);

      if (initialMissing [worker])
      {
          initialMissing [worker] = false;
          --initialResultsMissing;
      }

      // Check if termination criterions

      if (! shutDown)
      {         
         // Check requested errors

         EstErr totalEE = rc.getEstErr ();
         for (unsigned w = 1; w <= numWorkers(); ++w)  totalEE += localEE[w];

         status = checkRequestedError (totalEE, reqAbsError, reqRelError);

         shutDown = (status != ERROR && initialResultsMissing == 0)
                 || maxIter == 0;
      }

      if (! shutDown)
      {
         --maxIter;
         sendTopRegionToSlave (rc, worker, localEE, rc.size() / numWorkers());
      }
      else
      {
         --numActiveWorkers;
      }
   }

   sendTerminationSignal ();

   // Clean up Region Collection and recalculate result
      
   rc.killQueueAndUpdateResult ();
 
   ee = rc.getEstErr ();
   for (unsigned w = 1; w <= numWorkers(); ++w)  ee += localEE[w];

   return (status == ERROR) ? MAX_EVAL_REACHED : status; 
}


/**
 *  Algorithm for MASTER in INTERLEAVED mode
 */

Integrator::Status L::AdaptIntegratorMSInter::master (
   const Hypercube &h,
   real reqAbsError, real reqRelError, EstErr &ee, Index maxIter)
{
   const unsigned dim = h.getDimension ();
   maxIter /= numWorkers();
 
   Array<EstErr> localEE (nodes, EstErr (0.0, 0.0));
 
   RegionCollection rc;
 
   // Main loop
 
   for (Index i = 1; i != maxIter; ++i)
   {
      // Collect from and deal out to each node

      for (unsigned w = 1; w <= numWorkers(); ++w)
      {
         recvRegionsFromSlave (rc, w, dim, localEE);
         if (i < maxIter)
         {
            sendTopRegionToSlave (rc, w, localEE, (rc.size()/numWorkers()) + 1);
         }
      }

      // Check if error criterion is met
      
      EstErr totalEE = rc.getEstErr ();
      for (unsigned w = 1; w <= numWorkers(); ++w)  totalEE += localEE[w];

      Status status = checkRequestedError (totalEE, reqAbsError, reqRelError);

      // Shut down calculation estErr is acceptable

      if (status != ERROR)
      {
         // Collect all remainig results

         if (i < maxIter)
         {
            recvRegionsFromAllSlaves (rc, dim, localEE);

            // Check EstErr again

            totalEE = rc.getEstErr();
            for (unsigned w = 1; w <= numWorkers(); ++w)  totalEE += localEE[w];
            status = checkRequestedError (totalEE, reqAbsError, reqRelError);
         }

         if (status != ERROR)  break;

         // if the second check fails, we have to send regions to all nodes
         // and continue

         sendTopRegionToAllSlaves (rc, localEE, rc.size() / numWorkers());
      }
   }  // next iteration
 
   sendTerminationSignal();

   rc.killQueueAndUpdateResult ();
 
   ee = rc.getEstErr();
   for (unsigned w = 1; w <= numWorkers(); ++w)  ee += localEE[w];

   Status status = checkRequestedError (ee, reqAbsError, reqRelError);
 
   return (status == ERROR) ? MAX_EVAL_REACHED : status;
}


/**
 *  Algorithm for MASTER in SYNCHRONOUS mode
 */

Integrator::Status L::AdaptIntegratorMSSync::master (
   const Hypercube &h,
   real reqAbsError, real reqRelError, EstErr &ee, Index maxIter)
{ 
   const unsigned dim = h.getDimension ();
   maxIter /= numWorkers();
 
   Array<EstErr> localEE (nodes, EstErr (0.0, 0.0));
 
   RegionCollection rc;
 
   for (Index i = 1; ; ++i)
   {
      recvRegionsFromAllSlaves (rc, dim, localEE);

      // Check for absolute and relative error

      EstErr totalEE = rc.getEstErr();
      for (unsigned w = 1; w <= numWorkers(); ++w)  totalEE += localEE[w];
// cout << "Iteration: " << i << "   Error: " << totalEE.getError() << endl;
 
      Status status = checkRequestedError (totalEE, reqAbsError, reqRelError);
      if (status != ERROR || i == maxIter)  break;
 
      // Send new regions to slaves for processing
 
      sendTopRegionToAllSlaves (rc, localEE, rc.size() / numWorkers());
   }
 
   sendTerminationSignal ();
 
   // Clean up Region Collection and recalculate result
 
   rc.killQueueAndUpdateResult();
 
   ee = rc.getEstErr ();
   for (unsigned w = 1; w <= numWorkers(); ++w)  ee += localEE[w];
 
   Status status = checkRequestedError (ee, reqAbsError, reqRelError);
   return (status == ERROR) ? MAX_EVAL_REACHED : status;
}


/**
 *  Algorithm for the SLAVE nodes (rank >= 1)
 */

inline
void L::AdaptIntegratorMS::slave (
   Integrand &f, const Hypercube &h, EmbeddedRule& rule, unsigned interval)
{
   // PA::calc (PA::CS, rank);
   // Create Region Collection and push initial hypercube

   RegionCollection rc;

   storeSubcube (rc, h, numWorkers (), rank - 1, f, rule);

   // PA::calc (PA::CE, rank);

   for (;;)
   {
      // Refine _interval_ times

      // PA::calc (PA::CS, rank);
      for (unsigned i = 0; i < interval; ++i)  rc.refine (f, rule);
      // PA::calc (PA::CE, rank);

      // Send results to manager
 
      // PA::send (PA::SS, rank, 0);
      SendBuffer buffer (comm);

      buffer << rc.getEstErr();

      for (unsigned i = 0;
           i < std::min(RegionCollection::size_type(numSentBack), rc.size());
           ++i)
      {
         Region *region = rc.top();
 
         if (! (buffer << *region))  break;

         rc.pop(); 
         delete region;
      }

      buffer.send (0, RESULT);
      // PA::send (PA::SE, rank, 0);
 
      // Wait for response from master
 
      // PA::send (PA::RS, 0, rank);
      MPI_Status status;
      MPI_Probe (0, MPI_ANY_TAG, comm, &status);
 
      // What kind of message did we receive??
 
      switch (status.MPI_TAG)
      {
         case CUBE:   // Receive new Regions from master for processing
         {
            RecvBuffer buffer (comm, 0, CUBE);
            while (buffer)
            {
               rc.push (new Region (h.getDimension(), buffer));
            }

            break;
         }
 
         case DONE:   // We are done. -> return
         {
            int dummy;
            MPI_Recv (&dummy, 0, MPI_INT, 0, DONE, comm, &status);
            return;
         }
 
         default:  throw InternalError(__FILE__, __LINE__);
      }
      // PA::send (PA::RE, 0, rank);
   }
}


Integrator::Status L::AdaptIntegratorMS::integrate (
   Integrand &f, const Hypercube &h, Index maxEvaluations,
   real reqAbsError, real reqRelError, EstErr &ee) 
{
    checkDimension (h, f);

   // PA::calc (PA::CS, rank);

   // We need at least one slave for this algorithm
 
   if (numWorkers() == 0)
   {
      #ifdef HINTLIB_NO_EXCEPTIONS
         ee.set (0.0, 0.0);
         return ERROR;
      #else
         throw TooFewNodes();
      #endif
   }

   // Figure out load balancing interval

   std::auto_ptr<EmbeddedRule> rule (getRule (h.getDimension()));

   interval = 1 + frequency / rule->getNumPoints();

   // Determine the number of possible iterations
 
   Index maxIter = std::numeric_limits<Index>::max();

   if (maxEvaluations)
   {
      const unsigned dim = h.getDimension ();
      const unsigned w = numWorkers ();

      // Determine the number of points we need to get started

      Index initialPoints
         = (1 + 4 * dim) * (w - 1 + w * ms1(w-1) + 1)  // Splitting
         + rule->getNumPoints () * w;                  // Initial estimate

      // Make sure we can do at least one iteration

      if (maxEvaluations <
          Index (2 * interval * rule->getNumPoints() * w) + initialPoints)
      {
         #ifdef HINTLIB_NO_EXCEPTIONS
            ee.set (0.0, 0.0);
            return ERROR;
         #else
            throw NoEvaluationsPossible(maxEvaluations);
         #endif
      }

      maxIter = maxEvaluations - initialPoints;
      maxIter /= (2 * interval * rule->getNumPoints());
   }

   numSentBack = numWorkers();

   // PA::calc (PA::CE, rank);

   if (rank == 0)
   {
      return master (h, reqAbsError, reqRelError, ee, maxIter);
   }
   else
   {
      slave (f, h, *rule, interval);
      return WRONG_NODE;
   }
}
             

