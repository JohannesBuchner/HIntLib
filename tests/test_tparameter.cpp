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

#include <iostream>
#include <iomanip>
#include <memory>

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_CSTDLIB
  #include <cstdlib>
  #define HINTLIB_SLN std::
#else
  #include <stdlib.h>
  #define HINTLIB_SLN
#endif

#include <HIntLib/make.h>
#include <HIntLib/tparameter.h>
#include <HIntLib/array.h>
#include <HIntLib/hlalgorithm.h>

#include "test.h"

using namespace std;
using namespace HIntLib;

// Option processing

const char* options = "erbm:d:t";

bool ADD_EQUI = false;
bool RESTRICTED = false;
bool NO_BOUNDS = false;
bool PRINT_THICKNESS = false;
const unsigned MIN_M = 1;
const unsigned MIN_S = 1;
unsigned MAX_S = 55 + 1;
unsigned MAX_M = 25 + 1;

bool opt(int c, const char* s)
{
   switch (c)
   {
   case 'e':  ADD_EQUI = true; return true;
   case 'r':  RESTRICTED = true; return true;
   case 'b':  NO_BOUNDS = true; return true;
   case 'm':  MAX_M = HINTLIB_SLN atoi (s) + 1; return true;
   case 'd':  MAX_S = HINTLIB_SLN atoi (s) + 1; return true;
   case 't':  PRINT_THICKNESS = true; return true;
   }

   return false;
}

void usage()
{
   cerr <<
      "Usage: test_tparameter [OPTION] net...\n\n"
      "Tests if determining the t-parameter works correctly.\n\n"
      "  net    # of the GeneratorMatrix or filename of libseq-format matrix\n"
      << option_msg <<
      "  -e     Add equidistributed coordinate\n"
      "  -r     Test restricted t_parameter\n"
      "  -b     Do tests without optimized bounds\n"
      "  -m n   Maximum m (default = 25)\n"
      "  -d n   Maximum dimension (default = 55)\n"
      "  -t     Print thickness instead of t-parameter\n"
      "\n";

   exit (1);
}


// Other global names

typedef GeneratorMatrixGen<unsigned char> Matrix;


/**
 *  determineT ()
 */

void determineT (
      unsigned s, unsigned m, const Matrix& matrix, int* t_matrix,
      const Matrix* lowDimMatrix)
{
   int ub = m;  // trivial upper bound

   // t can only increase by one compared to m-1
   // With some thought, you see that it also works for the
   // equidistributed case.

   if (t_matrix [(m-1) * MAX_S + s] >= 0) ub = t_matrix [(m-1) * MAX_S + s] + 1;

   int lb = 0;  // trivial lower bound
   TOption opts = DEFAULT;

   // Do we know t for a lower-dimensional sub-matrix?

   if (lowDimMatrix && lowDimMatrix->getM() >= m
                    && lowDimMatrix->getTotalPrec() >=m)
   {
      GMCopy copy; copy.dim(s-1).m(m).totalPrec(m).equi(ADD_EQUI);
      Matrix g1 (*lowDimMatrix, copy);
      Matrix g2 (matrix,        copy);
      if (g1 == g2)
      {
         // t cannot decrease compared to s-1

         lb = std::max (lb, t_matrix [m * MAX_S + (s-1)]);
         opts = LOWER_DIM_OK;
      }
   }

   GMCopy copy;
   copy.dim(s).m(m).totalPrec(m).equi(ADD_EQUI);

   Matrix g (matrix, copy);

   int t = t_matrix [m * MAX_S + s] = tParameter (g, lb, ub, opts);

   if (NO_BOUNDS)
   {
      if (t != tParameter (g, 0, m))  error ("Bound optimization broken!");
   }

   // check if restricted results make sense

   if (RESTRICTED)
   {
      Array<int> t_restricted (MAX_M + 1);

      for (unsigned i = 0; i <= m; ++i)
      {
          t_restricted [i] = tParameterRestricted (g, 0, m, i, DEFAULT);
      }

      if (t_restricted[0] != 0)  error ("t_res(0) != 0");
      if (t_restricted[m] != t)  error ("t_res(m) != t");
      for (unsigned i = 0; i < m; ++i)
      {
         if (t_restricted[i] > t_restricted[i+1])
         {
            error ("t_res(i) > t_res(i+1)");
         }
      }

      NORMAL
      {
         cout << setw(3) << s << setw(3) << m << setw(4) << t;
         for (unsigned i = 0; i <= m; ++i)
         {
            cout << setw(3) << t_restricted[i];
         }
         cout << endl;
      }
   }
   else
   {
      NORMAL  cout << setw(3) << (PRINT_THICKNESS ? m - t : t) << flush;
   }
}


/**
 *  Main program
 */

void test (int argc, char** argv)
{
   if (argc != 1)  usage();

   char* endptr; 
   int matrix = HINTLIB_SLN strtol (argv[0], &endptr, 10);

   Array<Matrix*> matrices (MAX_S, 0);    // The fullsized matrix for each s
   Array<int> t_matrix (MAX_S * MAX_M, -1);

   if (endptr != argv[0] && *endptr == '\0')
   { 
      NORMAL cout << Make::getGeneratorMatrixGenName (matrix) << "\n\n";
      NORMAL if (! RESTRICTED)  cout << "   s=";

      for (unsigned s = MIN_S; s < MAX_S; ++s)
      {
         try
         {
            matrices [s] = Make::generatorMatrixGen (matrix, s - ADD_EQUI);
            NORMAL if (! RESTRICTED)  cout << setw(3) << s << flush;
         }
         catch (InvalidDimension &)
         {
            NORMAL if (! RESTRICTED)  cout << setw(3) << "";
         }
      }
      NORMAL
      {
         if (! RESTRICTED)  cout << "\n";
         else cout << setw(3) << "s" << setw(3) << "m" << setw(4) << "t"
                   << endl;
      }
   }
   else   // reading matrix from file
   { 
      std::auto_ptr<Matrix> gm (loadLibSeq (argv[0]));
      NORMAL cout << "Matrix from file '" << argv[0] << "'\n\n";
      NORMAL if (! RESTRICTED)  cout << "   s=";

      for (unsigned s = MIN_S; s < MAX_S; ++s)
      {
         unsigned ss = s - ADD_EQUI;

         if (ss <= gm->getDimension())
         {
            GMCopy copy; copy.dim(ss);
            matrices [s] = new Matrix (*gm, copy);
            NORMAL if (! RESTRICTED)  cout << setw(3) << s << flush;
         }
         else
         {
            NORMAL if (! RESTRICTED)  cout << setw(3) << "";
         }
      }
      NORMAL
      {
         if (! RESTRICTED)  cout << "\n";
         else cout << setw(3) << "s" << setw(3) << "m" << setw(4) << "t"
                   << endl;
      }
   }

   for (unsigned m = MIN_M; m < MAX_M; ++m)
   {
      NORMAL if (! RESTRICTED)  cout << "m=" << setw(2) << m << ' ';

      for (unsigned s = MIN_S; s < MAX_S; ++s)
      {
         Matrix* matrix = matrices [s];

         if (matrix && matrix->getM() >= m && matrix->getTotalPrec() >= m)
         {
            determineT (s, m, *matrix, t_matrix, matrices[s-1]);
         }
         else
         {
            NORMAL  if (! RESTRICTED)  cout << setw(3) << "";
         }
      }
      NORMAL  if (! RESTRICTED)  cout << endl;
   }

   purge (&matrices[0], &matrices[MAX_S]);
}

