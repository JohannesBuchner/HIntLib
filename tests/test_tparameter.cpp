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
#include <stdlib.h>

#include <HIntLib/make.h>
#include <HIntLib/generatormatrix.h>
#include <HIntLib/array.h>
#include <HIntLib/myalgorithm.h>

#include "test.h"

using namespace std;
using namespace HIntLib;

const char* options = "erbm:d:";

bool ADD_EQUI = false;
bool RESTRICTED = false;
bool NO_BOUNDS = false;
unsigned MAX_S = 55 + 1;
unsigned MAX_M = 25 + 1;

bool opt(int c, const char* s)
{
   switch (c)
   {
   case 'e':  ADD_EQUI = true; return true;
   case 'r':  RESTRICTED = true; return true;
   case 'b':  NO_BOUNDS = true; return true;
   case 'm':  MAX_M = atoi (s) + 1; return true;
   case 'd':  MAX_S = atoi (s) + 1; return true;
   }

   return false;
}

void usage()
{
   cerr <<
      "Usage: test_tparameter [OPTION] net...\n\n"
      "Tests if determining the t-parameter works correctly.\n\n"
      "  net    # of the GeneratorMatrix\n"
      << option_msg <<
      "  -e     Add equidistributed coordinate\n"
      "  -r     Test restricted t_parameter\n"
      "  -b     Do tests without optimized bounds\n"
      "  -m n   Maximum m (default = 25)\n"
      "  -d n   Maximum dimension (default = 55)\n"

      "\n";

   exit (1);
}

void test (int argc, char** argv)
{
   if (argc != 1)  usage();

   int matrix = atoi (argv[0]);

   const unsigned MIN_M = 1;
   const unsigned MIN_S = 1;

   typedef GeneratorMatrixGen<unsigned char> Matrix;
   typedef GeneratorMatrixGenCopy<unsigned char> MatrixCopy;

   Array<Matrix*> matrices (MAX_S, 0);    // The fullsized matrix for each s
   Array<int> t_matrix (MAX_S * MAX_M, -1);
   Array<int> t_restricted (MAX_M + 1);

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
      else cout << setw(3) << "s" << setw(3) << "m" << setw(4) << "t" << endl;
   }

   for (unsigned m = MIN_M; m < MAX_M; ++m)
   {
      NORMAL if (! RESTRICTED)  cout << "m=" << setw(2) << m << ' ';

      for (unsigned s = MIN_S; s < MAX_S; ++s)
      {
         if (matrices [s] && matrices[s]->getM() >= m)
         {
            int ub = m;  // trivial upper bound

            // t can only increase by one compared to m-1
            // With some thought, you see that it also works for the
            // equidistributed case.

            if (t_matrix [(m-1) * MAX_S + s] >= 0)
            {
               ub = t_matrix [(m-1) * MAX_S + s] + 1;
            }

            int lb = 0;  // trivial lower bound
            bool dimOpt = false;

            // Do we know t for a lower-dimensional sub-matrix?

            if (matrices [s-1])
            {
               GMCopy copy; copy.dim(s-1).m(m).totalPrec(m).equi(ADD_EQUI);
               MatrixCopy g1 (*matrices[s-1], copy);
               MatrixCopy g2 (*matrices[s],   copy);
               if (g1 == g2)
               {
                  // t cannot decrease compared to s-1

                  lb = std::max (lb, t_matrix [m * MAX_S + (s-1)]);
                  dimOpt = true;
               }
            }

            int t;
            int tt = 0;

            // Use the optimal implementaion depending on b and m

            GMCopy copy; copy.dim(s).m(m).totalPrec(m).equi(ADD_EQUI);

            if (matrices[s]->getBase() > 2)
            {
               GeneratorMatrixGenCopy<unsigned char> g (*matrices[s], copy);
               t = t_parameter (g, lb, ub, dimOpt);

               if (NO_BOUNDS)  tt = t_parameter (g, 0, m);

               if (RESTRICTED)
               {
                  for (unsigned i = 0; i <= m; ++i)
                  {
                     t_restricted [i] = t_parameter2 (g, 0, m, i);
                  }
               }
            }
            else if (m <= unsigned (std::numeric_limits<u32>::digits))
            {
               GeneratorMatrix2Copy<u32> g (*matrices[s], copy);
               t = t_parameter (g, lb, ub, dimOpt);

               if (NO_BOUNDS)  tt = t_parameter (g, 0, m);

               if (RESTRICTED)
               {
                  for (unsigned i = 0; i <= m; ++i)
                  {
                     t_restricted [i] = t_parameter2 (g, 0, m, i);
                  }
               }
            }
            else
            {
               GeneratorMatrix2Copy<u64> g (*matrices[s], copy);
               t = t_parameter (g, lb, ub, dimOpt);

               if (NO_BOUNDS)  tt = t_parameter (g, 0, m);

               if (RESTRICTED)
               {
                  for (unsigned i = 0; i <= m; ++i)
                  {
                     t_restricted [i] = t_parameter2 (g, 0, m, i);
                  }
               }
            }

            t_matrix [m * MAX_S + s] = t;

            if (NO_BOUNDS)
            {
               if (t != tt)  error ("Bound optimzation broken!");
            }

            // check if restricted results make sense

            if (RESTRICTED)
            {
               if (t_restricted[0] != 0)  error ("t_res(0) != 0");
               if (t_restricted[m] != t)  error ("t_res(m) != t");
               for (unsigned i = 0; i < m; ++i)
               {
                  if (t_restricted[i] > t_restricted[i+1])
                  {
                     error ("t_res(i) > t_res(i+1)");
                  }
               }
            }

            NORMAL
            {
               if (! RESTRICTED)
               {
                  cout << setw(3) << t << flush;
               }
               else
               {
                  cout << setw(3) << s << setw(3) << m << setw(4) << t;
                  for (unsigned i = 0; i <= m; ++i)
                  {
                     cout << setw(3) << t_restricted[i];
                  }
                  cout << endl;

               }
            }
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

