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

#include <HIntLib/hlmath.h>
#include <HIntLib/array.h>
#include <HIntLib/linearalgebra.h>
#include <HIntLib/linearalgebragen.h>
#include <HIntLib/linearalgebra2.h>
#include <HIntLib/mersennetwister.h>

#include "test.h"

using namespace HIntLib;

using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::setw;

const char* options = "x:n:";

int SEED = 0;
int NUM_TESTS = 1;

bool opt (int c, const char* s)
{
   switch (c)
   {
   case 'x':  SEED = atoi (s); return true;
   case 'n':  NUM_TESTS = atoi (s); return true;
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
      "  -x seed Seed for random number generator [0]\n"
      "  -n num Number of tests [1]\n"
      << option_msg <<
      "\n";

   exit (1);
}


/**
 *  printMatrix()
 *  printSquareMatrix()
 */

void printMatrix (const unsigned char* m, unsigned numRows, unsigned numCols)
{
   cout << '\n';

   for (unsigned row = 0; row < numRows; ++row)
   {
      for (unsigned col = 0; col < numCols; ++col)
      {
         cout << setw (3) << int (*m++);
      }
      cout << '\n';
   }
   cout << endl;
}

inline
void printSquareMatrix (const unsigned char* m, unsigned size)
{
   printMatrix (m, size, size);
}


/**
 *  doTest()
 */

void doTest (const int b, const unsigned char* m, const int size)
{
   const bool base2 = (b == 2) && (size <= 64);
   const unsigned size2 = sqr (size);
   std::auto_ptr<LinearAlgebra> la (LinearAlgebra::make (b));

   Array<unsigned char> temp (size2);
   Array<unsigned char> copy (size2);
   Array<u64> m2    (base2 ? size : 1);
   Array<u64> copy2 (base2 ? size : 1);
   Array<u64> temp2 (base2 ? size : 1);

   DEB1 printSquareMatrix (m, size);

   // packMatrix() and unpackMatrix()

   if (base2)
   {
      NORMAL cout << "packing and unpacking..." << endl;

      packMatrix (m, size, m2.begin(), m2.begin() + size);
      unpackMatrix (m2.begin(), m2.begin() + size, size, temp.begin());

      DEB1 printSquareMatrix (temp.begin(), size);

      if (! std::equal (temp.begin(), temp.begin() + size2, m))
      {
         error ("packMatrix() or unpackMatrix() failed!");
      }
   }

   // matrixMul(), vectorMatrixMul()

   NORMAL cout << "Multiplying matrix..." << flush;

   for (int i = 0; i < size; ++i)
   {
      la->vectorMatrixMul (m + i * size, m, size, size,
                           temp.begin() + i * size);
   }
   la->matrixMul (m, m, size, copy.begin());

   if (! std::equal (temp.begin(), temp.begin() + size2, copy.begin()))
   {
      printSquareMatrix (temp.begin(), size);
      printSquareMatrix (copy.begin(), size);
      
      error ("squareMatrixMul() and vectorMatrixMul() inconsistent");
   }

   if (base2)
   {
      for (int i = 0; i < size; ++i)
      {
         temp2[i] = vectorMatrixMul (m2[i], m2.begin(), m2.begin() + size);
      }

      unpackMatrix (temp2.begin(), temp2.begin() + size, size, temp.begin());

      if (! std::equal (temp.begin(), temp.begin() + size2, copy.begin()))
      {
         printSquareMatrix (temp.begin(), size);
         printSquareMatrix (copy.begin(), size);
         
         error ("vectorMatrixMul() for base 2 is different");
      }
   }
   
   // matrixMul(), matrixVectorMul()
   
   matrixTranspose (m, size, size, copy.begin());
   for (int i = 0; i < size; ++i)
   {
      la->matrixVectorMul (m, copy.begin() + i * size, size, size,
                           temp.begin() + i * size);
   }
   matrixTranspose (temp.begin(), size, size, copy.begin());
   la->matrixMul (m, m, size, temp.begin());

   if (! std::equal (temp.begin(), temp.begin() + size2, copy.begin()))
   {
      printSquareMatrix (temp.begin(), size);
      printSquareMatrix (copy.begin(), size);
      
      error ("squareMatrixMul() and matrixVectorMul() inconsistent");
   }

   if (base2)
   {
      matrixTranspose (m, size, size, copy.begin());
      packMatrix (copy.begin(), size, copy2.begin(), copy2.begin() + size);

      for (int i = 0; i < size; ++i)
      {
         temp2[i] = matrixVectorMul (m2.begin(), m2.begin() + size, copy2[i]);
      }

      unpackMatrix (temp2.begin(), temp2.begin() + size, size, copy.begin());
      matrixTranspose (copy.begin(), size, size, temp.begin());
      la->matrixMul (m, m, size, copy.begin());

      if (! std::equal (temp.begin(), temp.begin() + size2, copy.begin()))
      {
         printSquareMatrix (temp.begin(), size);
         printSquareMatrix (copy.begin(), size);
         
         error ("matrixVectorMul() for base 2 is different");
      }
   }

   // matrixRank()

   NORMAL cout << "Determining rank..." << flush;

   std::copy (m, m + size2, copy.begin());
   unsigned rank = la->matrixRank (copy.begin(), size, size);

   NORMAL cout << rank << endl;

   if (base2)
   {
      std::copy (m2.begin(), m2.begin() + size, copy2.begin());
      if (matrixRank (copy2.begin(), copy2.begin() + size) != rank)
      {
         error ("matrixRank() for base 2 is different!");
      }
   }

   // numLinearlyIndependentVectors()

   NORMAL cout << "Determining number of initial independent vectors..."
               << flush;

   std::copy (m, m + size2, copy.begin());
   unsigned numLi
      = la->numLinearlyIndependentVectors (copy, size, size);

   NORMAL cout << numLi << endl;

   if (numLi > rank)  error ("rank lower than number of independent vectors");

   // isLinearlyIndependent()

   NORMAL cout << "Checking for singularity..." << flush;

   std::copy (m, m + size2, copy.begin());
   bool indep = la->isLinearlyIndependent (copy, size, size);

   NORMAL cout << (indep ? "independent" : "singular") << endl;

   if ((rank == unsigned (size)) != indep)
   {
      error ("matrixRank() and isLinearlyIndependent() do not match");
   }

   if ((numLi == unsigned (size)) != indep)
   {
      error ("numLinearlyIndependentVectors() and isLinearlyIndependent() "
             "do not match!");
   }

   if (base2)
   {
      std::copy (m2.begin(), m2.begin() + size, copy2.begin());
      if (isLinearlyIndependent (copy2.begin(), copy2.begin() + size) != indep)
      {
         error ("isLinearlyIndependent() for base 2 is different!");
      }
   }

   // basisSupplement()

   NORMAL cout << "Determining basis supplement..." << flush;

   std::copy (m, m + size2, copy.begin());
   unsigned bs = la->basisSupplement (copy, size, size);

   NORMAL cout << " " << bs << " vectors valid" << endl;
   DEB1 printSquareMatrix (copy, size);

   if (bs != numLi)
   {
      error ("basisSupplement() and numLinearlyIndepenedentVectors() "
             "inconsistent");
   }

   if (! std::equal (copy.begin(), &copy[bs * size], m))
   {
      error ("Not a supplement to given basis");
   }

   if (! la->isLinearlyIndependent (copy, size, size))
   {
      error ("Result of basisSupplement() is not a basis");
   }

   // nullSpace()

   NORMAL cout << "Determining null space..." << flush;

   std::copy (m, m + size2, copy.begin());
   unsigned dimNS = la->nullSpace (copy, size, size, temp);

   NORMAL cout << " dim = " << dimNS << endl;

   DEB1  printMatrix (temp, dimNS, size);

   if (dimNS != size - rank)  error ("rank and dim(nullspace) do not match");

   for (unsigned i = 0; i < dimNS; ++i)
   {
      la->matrixVectorMul (m, &temp[i * size], size, size, copy);
      if (! la->isZeroMatrix (copy, size, 1))
      {
         error ("kernel base vector not in kernel");
      }
   }

   if (base2)
   {
      // normal

      std::copy (m2.begin(), m2.begin() + size, copy2.begin());

      unsigned dimNS2 =
         nullSpace (copy2.begin(), copy2.begin() + size, size, temp2.begin());

      if (dimNS != dimNS2)
      {
         cout << "Result: " << dimNS2 << '\n';
         error ("nullSpace() for base 2 yields different dimension!");
      }
      else
      {
         unpackMatrix (temp2.begin(), temp2.begin() + dimNS, size,
                       copy.begin());

         if (! std::equal (temp.begin(), temp.begin() + dimNS * size,
                           copy.begin()))
         {
            error ("nullSpace() for base 2 is different");
            printMatrix (temp.begin(), dimNS, size);
            printMatrix (copy.begin(), dimNS, size);
         }

         for (unsigned i = 0; i < dimNS; ++i)
         {
            if (matrixVectorMul (m2.begin(), m2.begin() + size, temp2 [i]) != 0)
            {
               error ("kernel base vector not in kernel");
            }
         }
      }

      // transposed

      matrixTranspose (m, size, size, copy.begin());
      packMatrix (copy.begin(), size, copy2.begin(), copy2.begin() + size);
      dimNS2 =
         nullSpaceT (copy2.begin(), copy2.begin() + size, temp2.begin());

      if (dimNS != dimNS2)
      {
         cout << "Result: " << dimNS2 << '\n';
         error ("nullSpaceT() for base 2 yields different dimension!");
      }
      else
      {
         unpackMatrix (temp2.begin(), temp2.begin() + dimNS, size,
                       copy.begin());

         if (! std::equal (temp.begin(), temp.begin() + dimNS * size,
                           copy.begin()))
         {
            error ("nullSpaceT() for base 2 is different");
            printMatrix (temp.begin(), dimNS, size);
            printMatrix (copy.begin(), dimNS, size);
         }

         for (unsigned i = 0; i < dimNS; ++i)
         {
            if (matrixVectorMul (m2.begin(), m2.begin() + size, temp2 [i]) != 0)
            {
               error ("kernel base vector not in kernel");
            }
         }
      }


   }

   // matrixInverse()

   NORMAL cout << "Determining inverse..." << flush;

   std::copy (m, m + size2, copy.begin());
   bool invResult = la->matrixInverse (copy.begin(), size);

   NORMAL cout << (invResult ? "ok" : "singular") << endl;

   if (invResult != indep)
   {
      error ("matrixInverse() singularity detection broken");
   }

   if (invResult)
   {
      DEB1 printSquareMatrix (copy, size);

      NORMAL cout << "Calculating product..." << endl;

      Array<unsigned char> product (size2);

      la->matrixMul (m, copy.begin(), size, product.begin());
   
      DEB1 printSquareMatrix (product, size);

      if (! la->isIdentityMatrix (product, size))
      {
         error ("m times m^-1 not equal to identity matrix");
      }
   }

   NORMAL cout << endl;
}

void test (int argc, char** argv)
{
   if (argc != 2)  usage();
      
   const int b    = atoi(argv[0]);
   const int size = atoi(argv[1]);
   const int size2 = size * size;

   if (b < 2 || size < 1)  usage();

   Array<unsigned char> m (size2, 0);
   doTest (b, m.begin(), size);

   for (int i = 0; i < size; ++i)  m[i * (size + 1)] = 1;
   doTest (b, m.begin(), size);

   MersenneTwister mt;
   mt.init (SEED);

   for (int i = 0; i < NUM_TESTS; ++i)
   {
      for (int i = 0; i < size2; ++i)  m[i] = mt(b);
      doTest (b, m.begin(), size);
   }
}


