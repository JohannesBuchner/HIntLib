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
 *  Statistic
 */

#ifndef STATISTIC_H
#define STATISTIC_H 1
 
#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/defaults.h>
#include <HIntLib/minmaxfinder.h>
#include <HIntLib/mymath.h>


namespace HIntLib
{

/**
 *  Statistic
 *
 *  Basic statistic class.
 *
 */

template<class Data = real, class Sum = Data, class CT = Index>
class Statistic
{
public:

   // Types

   typedef CT   CounterType;
   typedef Data DataType;
   typedef Sum  SumType;

   // Constructor

   Statistic() : sum(0), n(0) {}
   Statistic(real sum, CT n) : sum(sum), n(n) {}

   // Gather data

   void operator<< (Data x)  { sum += x; ++n; }

   // Get results

   Data  getSum()   const  { return Data(sum); }
   CT    getCount() const  { return n; }
   real  getMean()  const  { return real(sum) / real(n); }

   // assignment

   template<class Sum1>
   const Statistic<Data,Sum,CT>& operator= (const Statistic<Data,Sum1,CT>& s)
   {
      sum = s.getSum(); n = s.getCount(); return *this;
   }

   void reset()  { sum = 0; n = 0; }

   #ifdef PARALLEL
      void reduce (MPI_Comm = MPI_COMM_WORLD);
   #endif

private:
   Sum sum;
   CT n;
};

#ifdef PARALLEL

template<class Data, class Sum, class CT>
inline
void Statistic<Data,Sum,CT>::reduce (MPI_Comm comm)
{
   Data locSum = sum;
   Data gloSum;
   CT   gloN;
   
   MPI_Reduce (&locSum, &gloSum, 1, MPIType<Data>::type, MPI_SUM, 0, comm);
   MPI_Reduce (&n,      &gloN,   1, MPIType<CT>::type,   MPI_SUM, 0, comm);

   int rank;
   MPI_Comm_rank (comm, &rank);

   if (rank == 0)
   {
      sum = gloSum;
      n   = gloN;
   }
}

#endif


/**
 *  StatisticMinMax
 */

template<class Data = real, class Sum = Data, class CT = Index>
class StatisticMinMax : public Statistic<Data,Sum,CT>,
                        public MinMaxFinder<Data>
{
public:

   void operator<< (Data x)
   {
      Statistic<Data,Sum,CT>::operator<< (x); 
      MinMaxFinder<Data>::operator<< (x);
   }

   void reset()
   {
      Statistic<Data,Sum,CT>::reset();
      MinMaxFinder<Data>::reset();
   }
};


/**
 *  StatisticVar
 */

template<class Data = real, class Sum = Data, class CT = Index>
class StatisticVar : public Statistic<Data,Sum,CT>
{
public:

   StatisticVar() : Statistic<Data,Sum,CT>(), squares(0) {}
   StatisticVar(real sum, real squares, CT n)
      : Statistic<Data,Sum,CT>(sum, n), squares(squares) {}

   // Gather data

   void operator<< (Data x)
   {
      Statistic<Data,Sum,CT>::operator<< (x);
      squares += sqr(x);
   }

   // Return results

   real getVariance() const 
   {
      return real(squares) / real(getCount()) - sqr(getMean());
   }

   real getVarianceSample() const
   {
      return   (real(squares) - real(getCount()) * sqr(getMean()))
             / real(getCount()-1);
   }

   real getStdDev()       const  { return sqrt(getVariance()); }
   real getStdDevSample() const  { return sqrt(getVarianceSample()); }

   Data getSumSquares()    const { return Data(squares); }
   real getGeometricMean() const { return sqrt(real(squares)); }

   // assignment

   template<class Sum1>
   const StatisticVar<Data,Sum,CT>& operator=
      (const StatisticVar<Data,Sum1,CT>& s)
   {
      Statistic<Data,Sum,CT>::operator=(s); squares = s.getSumSquares();
      return *this;
   }

   void reset()  { Statistic<Data,Sum,CT>::reset(); squares = 0; }

   #ifdef PARALLEL
      void reduce (MPI_Comm = MPI_COMM_WORLD);
   #endif

private:
   Sum squares;
};

#ifdef PARALLEL

template<class Data, class Sum, class CT>
inline
void StatisticVar<Data,Sum,CT>::reduce (MPI_Comm comm)
{
   Statistic<Data,Sum,CT>::reduce(comm);

   Data locSquares = squares;
   Data gloSquares;
   
   MPI_Reduce (&locSquares, &gloSquares, 1, MPIType<Data>::type, MPI_SUM, 0, comm);

   int rank;
   MPI_Comm_rank (comm, &rank);

   if (rank == 0)  squares = gloSquares;
}

#endif

} // namespace HIntLib

#endif

