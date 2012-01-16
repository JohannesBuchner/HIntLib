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
 *  Creates a C file containing an array with the Niederreiter v vectors.
 *
 *  The algorithm used to create this array follows Bratley, Fox, and
 *  Niederreiter; Implementation and Tests of Low-Discrepancy Sequences in
 *  ACM TOMS Vol 2 No 3 (1992), pp 195 - 213
 */

#include <iostream>
#include <iomanip>
#include <algorithm>

// This is not really a library object. However, since library files are
// linked static, we do not require DLL-export names.
#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/niederreitermatrix.h>
#include <HIntLib/niederreitermatrixgen.h>
#include <HIntLib/generatormatrixgen.h>

using std::cout;
using std::setw;

using namespace HIntLib;


/**
 *  Main function
 *
 */

int main()
{
   // Create matrix using the general base algorithm

   GeneratorMatrixGen<unsigned char> matrix (
         2,
         NiederreiterMatrix::MAX_DIM,
         GeneratorMatrix::DEFAULT_M_BASE2,
         GeneratorMatrix2<u64>::CORR_DEFAULT_TOTALPREC_BASE2);
   initNiederreiter (matrix);

   // Convert to base 2

   GeneratorMatrix2<u64> matrixBase2 (matrix);

   // Generate output file

   cout << "/***********************************************************/\n"
           "/***   This file is program-generated!                   ***/\n"
           "/***                                                     ***/\n"
           "/***   Do not change!!!                                  ***/\n"
           "/***                                                     ***/\n"
           "/***   Update " __FILE__            " to update this file. */\n"
           "/***********************************************************/\n"
           "\n"
           "#ifdef __GNUG__\n"
           "#pragma implementation\n"
           "#endif\n"
           "\n"
           "#define HINTLIB_LIBRARY_OBJECT\n"
           "\n"
           "#include <HIntLib/niederreitermatrix.h>\n"
           "\n"
           "namespace L = HIntLib;\n"
           "\n"
           "const L::u64 "
              "L::NiederreiterMatrix::v_mem [DEFAULT_M_BASE2][MAX_DIM] =\n";

   matrixBase2.CArrayDump (cout);

#if 0
   cout << "\n"
           "const unsigned L::NiederreiterMatrix::t_s [MAX_DIM] =\n"
           "{\n";

   unsigned sum = 0;

   for (unsigned i = 0; i < NiederreiterMatrix::MAX_DIM; ++i)
   {
      // sum += irredPolys [i].degree() - 1;

      cout << setw(6) << sum << ",   // t_s [" << i + 1 << "]\n";
   }

   cout << "};"
#endif

   cout << "\n\n";

   return 0;
}


