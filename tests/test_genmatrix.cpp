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

#include <HIntLib/sobolmatrix.h>
#include <HIntLib/niederreitermatrix.h>

#include "test.h"

using namespace std;
using namespace HIntLib;


const char* options = "";

bool opt(int, const char*) { return false; }

void usage()
{
   cerr <<
      "Usage: test_genmatrix [OPTION]...\n\n"
      "Checks GeneratorMatrix2 and GeneratorMatrixGen.\n\n"
      << option_msg <<
      "\n";

   exit (1);
}

void test (int argc, char**)
{
   if (argc)  usage();

   NORMAL  cout << "Creating Sobol's matrix" << endl;
   
   SobolMatrix sobol;
   DEB1  sobol.dump (cout);

   NORMAL  cout << "Testing: sobol==sobol" << endl;

   if (! (sobol == sobol))
   {
      error ("Testing: sobol==sobol    failed!");
   }

   typedef GeneratorMatrixGenCopy8 GMgen;
   typedef GeneratorMatrix2Copy<u32> GM2;

   NORMAL  cout << "Creating Gen copy..." << endl;
   GeneratorMatrixGenCopy8 sobolGen (sobol);
   DEB1  sobolGen.dump (cout);

   NORMAL  cout << "Creatin base-2 copy..." << endl;
   GeneratorMatrix2Copy<u64> sobol2 (sobolGen);
   DEB1  sobol2.dump (cout);

   NORMAL  cout << "Testing: Equal in base 2" << endl;

   if (! (sobol == sobol2))
   {
      error ("Testing: Equal in base 2      failed!");
   }

   GMgen sobolGenAgain (sobol2);

   NORMAL  cout << "Testing: Equal in gen-base" << endl;

   if (! (sobolGen == sobolGenAgain))
   {
      error ("Testing: Equal in gen-basee      failed!");
   }

   NORMAL  cout << "Creating 32-bit copy of sobol..." << endl;
   GM2 sobol32 (sobol);

   // m < prec

   NORMAL  cout << "Creating truncated base-2   matrix, m < prec..." << endl;
   GM2 sobol2truncatedA (sobol32, 10, 13, 20);
   DEB1  sobol2truncatedA.dump (cout);
   
   NORMAL  cout << "Creating truncated gen-base matrix, m < prec..." << endl;
   GeneratorMatrixGenCopy8 sobolGenTruncatedA (sobolGen, 10, 13, 20);
   DEB1  sobolGenTruncatedA.dump (cout);

   NORMAL cout << "Comparing both matrices..." << endl;
   if (! (sobol2truncatedA == sobolGenTruncatedA))
   {
      error ("Truncated matrices, m < prec  are not equal!");
   }

   // m > prec

   NORMAL  cout << "Creating truncated base-2   matrix, m > prec..." << endl;
   GM2 sobol2truncatedB (sobol32, 10, 20, 13);
   DEB1  sobol2truncatedB.dump (cout);
   
   NORMAL  cout << "Creating truncated gen-base matrix, m > prec..." << endl;
   GMgen sobolGenTruncatedB (sobolGen, 10, 20, 13);
   DEB1  sobolGenTruncatedB.dump (cout);

   NORMAL  cout << "Comparing both matrices..." << endl;
   if (! (sobol2truncatedB == sobolGenTruncatedB))
   {
      error ("Truncated matrices, m > prec  are not equal!");
   }

   // m < prec, equi

   NORMAL
      cout << "Creating truncated base-2   matrix, m < prec,  equi..." << endl;
   GM2 sobol2truncatedEquiA (sobol32, 10, 13, 20, true);
   DEB1  sobol2truncatedEquiA.dump (cout);

   NORMAL
      cout << "Creating truncated gen-base matrix, m < prec, equi..." << endl;
   GMgen sobolGenTruncatedEquiA (sobolGen, 10, 13, 20, true);
   DEB1  sobolGenTruncatedEquiA.dump (cout);

   NORMAL  cout << "Comparing both matrices..." << endl;
   if (! (sobol2truncatedEquiA == sobolGenTruncatedEquiA))
   {
      error ("Truncated matrices, m < prec, equi  are not equal!");
   }

   // m > prec, equi

   NORMAL
      cout << "Creating truncated base-2   matrix, m > prec,  equi..." << endl;
   GM2 sobol2truncatedEquiB (sobol32, 10, 20, 13, true);
   DEB1  sobol2truncatedEquiB.dump (cout);

   NORMAL
      cout << "Creating truncated gen-base matrix, m > prec, equi..." << endl;
   GMgen sobolGenTruncatedEquiB (sobolGen, 10, 20, 13, true);
   DEB1  sobolGenTruncatedEquiB.dump (cout);

   NORMAL cout << "Comparing both matrices..." << endl;
   if (! (sobol2truncatedEquiB == sobolGenTruncatedEquiB))
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

