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

#include <HIntLib/array.h>
#include <HIntLib/generatormatrix2.h>
#include <HIntLib/generatormatrixgen.h>
#include <HIntLib/generatormatrixvec.h>
#include <HIntLib/sobolmatrix.h>
#include <HIntLib/mersennetwister.h>

#include "test.h"

using namespace std;
using namespace HIntLib;

unsigned DIM = 40;

const char* options = "d:";

bool opt (int c, const char* s)
{
   switch (c)
   {
      case 'd':  DIM = atoi (s); return true;
   }

   return false;
}

void usage()
{
   cerr <<
      "Usage: test_genmatrix [OPTION]...\n\n"
      "Checks GeneratorMatrix2 and GeneratorMatrixGen.\n\n"
      << option_msg <<
      "  -d dim Dimension of the matrices (default = 40).\n"
      "\n";

   exit (1);
}


typedef MutableGeneratorMatrixGen<unsigned char> MGen;
typedef MutableGeneratorMatrixVec<unsigned char> MVec;
typedef MutableGeneratorMatrix2<u64> M2;

typedef HeapAllocatedGeneratorMatrixGen<unsigned char> HMGen;
typedef HeapAllocatedGeneratorMatrixVec<unsigned char> HMVec;
typedef HeapAllocatedGeneratorMatrix2<u64> HM2;

typedef GeneratorMatrixGenCopy<unsigned char> CMGen;
typedef GeneratorMatrixVecCopy<unsigned char> CMVec;
typedef GeneratorMatrix2Copy<u64> CM2;

MersenneTwister mt;

void initMatrix (GeneratorMatrix &m)
{
   mt.init (42);

   for (unsigned d = 0; d < m.getDimension(); ++d)
   {
      for (unsigned r = 0; r < m.getM(); ++ r)
      {
         for (unsigned b = 0; b < m.getTotalPrec(); ++b)
         {
            m.setDigit (d, r, b, mt (m.getBase()));
         }
      }
   }
}

bool verifyMatrix (GeneratorMatrix &m)
{
   mt.init (42);

   for (unsigned d = 0; d < m.getDimension(); ++d)
   {
      for (unsigned r = 0; r < m.getM(); ++ r)
      {
         for (unsigned b = 0; b < m.getTotalPrec(); ++b)
         {
            if (int (m.getDigit (d, r, b)) != mt (m.getBase()))
            {
               return false;
            }
         }
      }
   }

   return true;
}

void testOneMatrix (GeneratorMatrix &m)
{
   initMatrix (m);

   if (! verifyMatrix (m))  error ("getDigit()/setDigit() broken!");

   if (! (m == m))  error ("operator==() broken!");
}

template<class T>
void testTwoMatrices (GeneratorMatrix& m1, T& m2)
{
   initMatrix (m1);

   assign (m1, m2);

   if (! verifyMatrix (m2))
   {
      error ("assign() broken!");
      DEB1
      {
         m1.dump (cout);
         m2.dump (cout);
      }
   }

   if (! (m1 == m2))  error ("operator==() broken!");

   for (unsigned i = 0; i < 5; ++i)
   {
      unsigned d = mt (m2.getDimension());
      unsigned r = mt (m2.getM());
      unsigned b = mt (m2.getTotalPrec());
      unsigned char ori = m2.getDigit (d,r,b);
      m2.setDigit (d,r,b, (ori+1) % m2.getBase());

      if (m2 == m1)
      {
         error ("operator==() does not detect error");
      }
      m2.setDigit (d,r,b, ori);
   }
}

void doTests (unsigned base)
{
   NORMAL cout << "Allocating matrices" << endl;

   unsigned numMatrices
      = digitsRepresentable(static_cast<unsigned char>(base)) + 1
      + (base == 2);
   DEB1 cout << "  Number of matrices: " << numMatrices << endl;

   Array<GeneratorMatrix*> matrices (numMatrices);

   for (unsigned vec = 0;
        vec < digitsRepresentable(static_cast<unsigned char>(base)); ++vec)
   {
      matrices [vec] = new HMVec (base, vec + 1, DIM);
   }

   if (base == 2)
   {
      matrices [numMatrices - 2] = new HMGen (base, DIM);
      matrices [numMatrices - 1] = new HM2 (DIM);
   }
   else
   {
      matrices [numMatrices - 1] = new HMGen (base, DIM);
   }

   NORMAL cout << "Testing Gen, base=" << base << endl;
   
   HMGen mgen (base, DIM);
   testOneMatrix (mgen);

   for (unsigned i = 0; i < numMatrices; ++i)
   {
      DEB1 cout << "  ...in combination with matrix " << i << endl;
      testTwoMatrices (*matrices[i], mgen);
      testTwoMatrices (*matrices[i], static_cast<GeneratorMatrix&>(mgen));
   }

   for (unsigned vec = 1;
        vec <= digitsRepresentable(static_cast<unsigned char>(base)); ++vec)
   {
      NORMAL cout << "Testing Vec, vec=" << vec << ", base=" << base << endl;
   
      HMVec mvec (base, vec, DIM);
      testOneMatrix (mvec);

      for (unsigned i = 0; i < numMatrices; ++i)
      {
         DEB1 cout << "  ...in combination with matrix " << i << endl;
         testTwoMatrices (*matrices[i], mvec);
         testTwoMatrices (*matrices[i], static_cast<GeneratorMatrix&>(mvec));
      }
   }

   if (base == 2)
   {
      NORMAL cout << "Testing 2, base=" << base << endl;
   
      HM2 m2 (DIM);
      testOneMatrix (m2);

      for (unsigned i = 0; i < numMatrices; ++i)
      {
         DEB1 cout << "  ...in combination with matrix " << i << endl;

         testTwoMatrices (*matrices[i], m2);
         testTwoMatrices (*matrices[i], static_cast<GeneratorMatrix&>(m2));
      }
   }
}


void oldTests ()
{
   SobolMatrix sobol;
   
   typedef GeneratorMatrixGenCopy8 GMGen;
   typedef GeneratorMatrixVecCopy8 GMVec;
   typedef GeneratorMatrix2Copy<u32> GM2;

   NORMAL  cout << "Creating Gen copy..." << endl;
   GMGen sobolGen (sobol);
   DEB2  sobolGen.dump (cout);

   NORMAL  cout << "Creatin base-2 copy..." << endl;
   GeneratorMatrix2Copy<u64> sobol2 (sobolGen);
   DEB2  sobol2.dump (cout);

   NORMAL  cout << "Testing: Equal in base 2" << endl;

   if (! (sobol == sobol2))
   {
      error ("Testing: Equal in base 2      failed!");
   }

   GMGen sobolGenAgain (sobol2);

   NORMAL  cout << "Testing: Equal in gen-base" << endl;

   if (! (sobolGen == sobolGenAgain))
   {
      error ("Testing: Equal in gen-basee      failed!");
   }

   NORMAL  cout << "Creating 32-bit copy of sobol..." << endl;
   GM2   sobol32 (sobol);
   DEB2  sobol32.dump(cout);

   // m < prec

   GMCopy copy1;
   copy1.dim(10).m(13).totalPrec(20).vec(5);

   NORMAL  cout << "Creating truncated base-2   matrix, m < prec..." << endl;
   GM2   sobol2truncatedA (sobol32, copy1);
   DEB2  sobol2truncatedA.dump (cout);
   
   NORMAL  cout << "Creating truncated gen-base matrix, m < prec..." << endl;
   GMGen sobolGenTruncatedA (sobolGen, copy1);
   DEB2  sobolGenTruncatedA.dump (cout);

   DEB1 cout << "Comparing..." << endl;
   if (! (sobol2truncatedA == sobolGenTruncatedA))
   {
      error ("Truncated matrices, m < prec  are not equal!");
   }

   NORMAL  cout << "Creating truncated vec-base matrix, m < prec..." << endl;
   GMVec sobolVecTruncatedA (sobolGen, copy1);
   DEB2  sobolVecTruncatedA.dump (cout);

   DEB1 cout << "Comparing..." << endl;
   if (! (sobol2truncatedA == sobolVecTruncatedA))
   {
      error ("Truncated matrices, m < prec  are not equal!");
   }

   // m > prec

   GMCopy copy2;
   copy2.dim(10).m(20).totalPrec(13).vec(5);

   NORMAL  cout << "Creating truncated base-2   matrix, m > prec..." << endl;
   GM2 sobol2truncatedB (sobol32, copy2);
   DEB2  sobol2truncatedB.dump (cout);
   
   NORMAL  cout << "Creating truncated gen-base matrix, m > prec..." << endl;
   GMGen sobolGenTruncatedB (sobolGen, copy2);
   DEB2  sobolGenTruncatedB.dump (cout);

   DEB1  cout << "Comparing..." << endl;
   if (! (sobol2truncatedB == sobolGenTruncatedB))
   {
      error ("Truncated matrices, m > prec  are not equal!");
   }

   NORMAL  cout << "Creating truncated vec-base matrix, m > prec..." << endl;
   GMVec sobolVecTruncatedB (sobolGen, copy2);
   DEB2  sobolVecTruncatedB.dump (cout);

   DEB1  cout << "Comparing..." << endl;
   if (! (sobol2truncatedB == sobolVecTruncatedB))
   {
      error ("Truncated matrices, m > prec  are not equal!");
   }

   // m < prec, equi

   copy1.equi();
   
   NORMAL
      cout << "Creating truncated base-2   matrix, m < prec, equi..." << endl;
   GM2 sobol2truncatedEquiA (sobol32, copy1);
   DEB2  sobol2truncatedEquiA.dump (cout);

   NORMAL
      cout << "Creating truncated gen-base matrix, m < prec, equi..." << endl;
   GMGen sobolGenTruncatedEquiA (sobolGen, copy1);
   DEB2  sobolGenTruncatedEquiA.dump (cout);

   DEB1  cout << "Comparing..." << endl;
   if (! (sobol2truncatedEquiA == sobolGenTruncatedEquiA))
   {
      error ("Truncated matrices, m < prec, equi  are not equal!");
   }

   NORMAL
      cout << "Creating truncated vec-base matrix, m < prec, equi..." << endl;
   GMVec sobolVecTruncatedEquiA (sobolGen, copy1);
   DEB2  sobolVecTruncatedEquiA.dump (cout);

   DEB1  cout << "Comparing..." << endl;
   if (! (sobol2truncatedEquiA == sobolVecTruncatedEquiA))
   {
      error ("Truncated matrices, m < prec, equi  are not equal!");
   }

   // m > prec, equi

   copy2.equi();

   NORMAL
      cout << "Creating truncated base-2   matrix, m > prec, equi..." << endl;
   GM2   sobol2truncatedEquiB (sobol32, copy2);
   DEB2  sobol2truncatedEquiB.dump (cout);

   NORMAL
      cout << "Creating truncated gen-base matrix, m > prec, equi..." << endl;
   GMGen sobolGenTruncatedEquiB (sobolGen, copy2);
   DEB2  sobolGenTruncatedEquiB.dump (cout);

   DEB1 cout << "Comparing..." << endl;
   if (! (sobol2truncatedEquiB == sobolGenTruncatedEquiB))
   {
      error ("Truncated matrices, m > prec, equi  are not equal!");
   }

   NORMAL
      cout << "Creating truncated vec-base matrix, m > prec, equi..." << endl;
   GMVec sobolVecTruncatedEquiB (sobolGen, copy2);
   DEB2  sobolVecTruncatedEquiB.dump (cout);

   DEB1 cout << "Comparing..." << endl;
   if (! (sobol2truncatedEquiB == sobolVecTruncatedEquiB))
   {
      error ("Truncated matrices, m > prec, equi  are not equal!");
   }

   // prepareForGrayCode

   NORMAL cout << "Testing prepareForGrayCode() ..." << endl;

   GM2 grayCodeA (sobol2truncatedA);
   grayCodeA.prepareForGrayCode();

   if (grayCodeA == sobol2truncatedA)
   {
      error ("prepareForGrayCode() was null operation");
   }

   NORMAL cout << "Testing restoreFromGrayCode() ..." << endl;

   grayCodeA.restoreFromGrayCode();

   if (grayCodeA != sobol2truncatedA)
   {
      error ("restoreFromGrayCode() does not undo prepareForGrayCode()");
   }
}


void test (int argc, char**)
{
   if (argc)  usage();

   doTests (2);
   doTests (3);
   doTests (4);
   doTests (5);
   doTests (7);
   doTests (9);
   doTests (11);
   doTests (13);

   oldTests();
}


