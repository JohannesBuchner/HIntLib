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

#include <HIntLib/mymath.h>
#include <HIntLib/array.h>
#include <HIntLib/linearalgebra.h>
#include <HIntLib/mersennetwister.h>

#include "test.h"

using namespace HIntLib;

using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::setw;

const char* options = "";

bool opt (int c, const char*)
{
   switch (c)
   {
   // case 'm':  MAX_M =atoi (s); return true;
   }

   return false;
}

void usage()
{
   cerr <<
      "Usage: test_linearalgebra q size\n\n"
      "Exercises functions like isLinearIndependent(), matrixRank(), and\n"
      "matrixInverse().\n\n"
      "  q      Use arithmetic over F_q\n"
      "  size   Size of the matrices\n"
      << option_msg <<
      "\n";

   exit (1);
}

void printSquareMatrix (const unsigned char* m, unsigned size)
{
   for (unsigned row = 0; row < size; ++row)
   {
      for (unsigned col = 0; col < size; ++col)
      {
         cout << setw (3) << int (m[row * size + col]);
      }
      cout << '\n';
   }
}

void test (int argc, char** argv)
{
   if (argc != 2)  usage();
      
   const int b    = atoi(argv[0]);
   const int size = atoi(argv[1]);

   if (b < 2 || size < 1)  usage();

   const unsigned size2 = sqr (size);
   MersenneTwister mt;
   std::auto_ptr<LinearAlgebra> la (makeLinearAlgebra(b));

   NORMAL cout << "Creating random matrix..." << endl;

   Array<unsigned char> m (size2);
   for (unsigned i = 0; i < size2; ++i)  m[i] = mt(b);

   DEB1
   {
      cout << "\n";
      printSquareMatrix (m, size);
      cout << endl;
   }

   NORMAL cout << "Determining rank..." << flush;

   Array<unsigned char> copy (m.begin(), size2);
   unsigned rank = la->matrixRank (copy.begin(), size, size);

   NORMAL cout << rank << endl;

   NORMAL cout << "Checking for singularity..." << flush;

   std::copy (m.begin(), m.begin() + size2, copy.begin());
   bool indep = la->isLinearlyIndependent (copy.begin(), size, size);

   NORMAL cout << (indep ? "independent" : "singular") << endl;

   if ((rank == unsigned (size)) != indep)
   {
      error ("matrixRank() and isLinearlyIndependent() do not match!");
   }

   NORMAL cout << "Determining inverse..." << flush;

   std::copy (m.begin(), m.begin() + size2, copy.begin());
   bool invResult = la->matrixInverse (copy.begin(), size);

   NORMAL cout << (invResult ? "ok" : "singular") << endl;

   if (invResult != indep)
   {
      error ("matrixInverse() singularity detection broken!");
   }

   if (invResult)
   {
      DEB1
      {
         cout << "\n";
         printSquareMatrix (copy, size);
         cout << endl;
      }

      NORMAL cout << "Calculating product..." << endl;

      Array<unsigned char> product (size2);

      la->matrixMul (product.begin(), m.begin(), copy.begin(), size);
   
      DEB1
      {
         cout << "\n";
         printSquareMatrix (product, size);
         cout << endl;
      }

      for (int row = 0; row < size; ++row)
      {
         for (int col = 0; col < size; ++col)
         {
            if ((row == col) != product[row * size + col])
            {
               error ("m times m^-1 not equal to identity matrix!");
            }
         }
      }
   }
}

