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


#ifndef HINTLIB_HLMPI_H
#define HINTLIB_HLMPI_H 1

#include <mpi.h>

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/fallback_limits.h>
#endif

#define HINTLIB_PARALLEL 1

namespace HIntLib
{

class MPI
{
public:
   MPI (int *argc, char ***argv);
   ~MPI ();
   int getRank()  { return rank; }

private:
   int rank;
};

/**
 * class MPIType
 *
 * Make the appropriate MPIType availabel for each builit-in datatype.
 *
 * A template mechanism similar to <limits> is used.
 */

template<class T> class MPIType {};

template<> class MPIType<char>
   { public: static const MPI_Datatype type =
       std::numeric_limits<char>::is_signed ? MPI_CHAR : MPI_UNSIGNED_CHAR; };
template<> class MPIType<signed char>
   { public: static const MPI_Datatype type = MPI_CHAR; };
template<> class MPIType<short>
   { public: static const MPI_Datatype type = MPI_SHORT; };
template<> class MPIType<int>
   { public: static const MPI_Datatype type = MPI_INT; };
template<> class MPIType<long>
   { public: static const MPI_Datatype type = MPI_LONG; };
#ifdef HINTLIB_HAVE_LONG_LONG_INT
template<> class MPIType<long long>
   { public: static const MPI_Datatype type = MPI_LONG_LONG_INT; };
#endif
template<> class MPIType<unsigned char>
   { public: static const MPI_Datatype type = MPI_UNSIGNED_CHAR; };
template<> class MPIType<unsigned short>
   { public: static const MPI_Datatype type = MPI_UNSIGNED_SHORT; };
template<> class MPIType<unsigned>
   { public: static const MPI_Datatype type = MPI_UNSIGNED; };
template<> class MPIType<unsigned long>
   { public: static const MPI_Datatype type = MPI_UNSIGNED_LONG; };
#ifdef HINTLIB_HAVE_UNSIGNED_LONG_LONG_INT
template<> class MPIType<unsigned long long>
   { public: static const MPI_Datatype type = MPI_LONG_LONG_INT; }; // XXX
#endif
template<> class MPIType<float>
   { public: static const MPI_Datatype type = MPI_FLOAT; };
template<> class MPIType<double>
   { public: static const MPI_Datatype type = MPI_DOUBLE; };
#ifdef HINTLIB_HAVE_LONG_DOUBLE
template<> class MPIType<long double>
   { public: static const MPI_Datatype type = MPI_LONG_DOUBLE; };
#endif

}  // namespace HIntLib

#endif

