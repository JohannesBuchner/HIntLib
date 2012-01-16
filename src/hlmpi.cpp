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
 *  MPI
 */

#ifdef __GNUG__
#pragma implementation
#endif

#include <HIntLib/hlmpi.h>

#include <HIntLib/exception_MPI.h>

namespace L = HIntLib;

L::MPI::MPI(int *argc, char ***argv)
{
   int initialized;
   int status = MPI_Initialized (&initialized);
   if (status != MPI_SUCCESS || initialized)  throw MPIError(status);

   status = MPI_Init (argc, argv);
   if (status != MPI_SUCCESS)  throw MPIError(status);

   status = MPI_Comm_rank (MPI_COMM_WORLD, &rank);
   if (status != MPI_SUCCESS)
   {
      MPI_Finalize();
      throw MPIError(status);
   }
}

L::MPI::~MPI ()
{
   int status = MPI_Finalize();
   if (status != MPI_SUCCESS)  throw MPIError(status);
}

