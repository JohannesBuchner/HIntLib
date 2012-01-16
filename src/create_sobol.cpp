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
 *   Creates a C file (sobol.cpp) containing the v array for Sobol's sequence.
 *
 *   The setup routine for the v array is a modified version of a program from
 *   Wolfgang Ch. Schmid
 */

#include <iostream>
#include <iomanip>

#include <HIntLib/polynomial2.h>

#include <HIntLib/sobolmatrix.h>

using std::cout;
using std::setw;

using namespace HIntLib;

typedef Polynomial2<u32> P;


/*
 *  List of MAX_DIM primitive polynomials on the field {0,1}, ordered by degree.
 *
 *  Polynomials are not sorted to allow property A and A'.  However, all
 *  possible polynomials of lower degree are listed befor the first polynomial
 *  of the next higher degree.
 */

const P polynomials [SobolMatrix::MAX_DIM] =
{
   P (2),    // deg 1   x       this one is NOT primitive!!!
   P (3),    //         x+1
   P (7),    // deg 2   x²+x+1
   P (11),   // deg 3   x³+x+1
   P (13),   //         x³+x²+1
   P (19),   // deg 4   x^4+x+1
   P (25),   //         x^4+x³+1
   P (37),   // deg 5   x^5+x²+1
   P (59),
   P (47),
   P (61),
   P (55),
   P (41),
   P (67),   // deg 6
   P (97),
   P (91),
   P (109),
   P (103),
   P (115),
   P (131),  // deg 7
   P (193),
   P (137),
   P (145),
   P (143),
   P (241),
   P (157),
   P (185),
   P (167),
   P (229),
   P (171),
   P (213),
   P (191),
   P (253),
   P (203),
   P (211),
   P (239),
   P (247),
   P (285),  // deg 8
   P (369),
   P (299)
};


const unsigned MAX_DEGREE = 8;  // Maximum degree of polynomials


/**
 *  Initial values for the recursion
 *
 *  We don't know how to calculate these values, so it's impossible to use 
 *  this class for dimensions >40. :-( 
 *
 *  The values listed here are taken from Bratley and Fox, ACM TOMS 14.1 (1988)
 *  pages 88-100
 */
  
const u64 vinit [SobolMatrix::MAX_DIM][MAX_DEGREE] =
{  //  x   x²   x³  x^4  x^5  x^6  x^7  x^8
   {   1,   0,   0,   0,   0,   0,   0,   0},
   {   1,   0,   0,   0,   0,   0,   0,   0},
   {   1,   1,   0,   0,   0,   0,   0,   0},
   {   1,   3,   7,   0,   0,   0,   0,   0},
   {   1,   1,   5,   0,   0,   0,   0,   0},
   {   1,   3,   1,   1,   0,   0,   0,   0},
   {   1,   1,   3,   7,   0,   0,   0,   0},
   {   1,   3,   3,   9,   9,   0,   0,   0},
   {   1,   3,   7,  13,   3,   0,   0,   0},
   {   1,   1,   5,  11,  27,   0,   0,   0},
   {   1,   3,   5,   1,  15,   0,   0,   0},
   {   1,   1,   7,   3,  29,   0,   0,   0},
   {   1,   3,   7,   7,  21,   0,   0,   0},
   {   1,   1,   1,   9,  23,  37,   0,   0},
   {   1,   3,   3,   5,  19,  33,   0,   0},
   {   1,   1,   3,  13,  11,   7,   0,   0},
   {   1,   1,   7,  13,  25,   5,   0,   0},
   {   1,   3,   5,  11,   7,  11,   0,   0},
   {   1,   1,   1,   3,  13,  39,   0,   0},
   {   1,   3,   1,  15,  17,  63,  13,   0},
   {   1,   1,   5,   5,   1,  27,  33,   0},
   {   1,   3,   3,   3,  25,  17, 115,   0},
   {   1,   1,   3,  15,  29,  15,  41,   0},
   {   1,   3,   1,   7,   3,  23,  79,   0},
   {   1,   3,   7,   9,  31,  29,  17,   0},
   {   1,   1,   5,  13,  11,   3,  29,   0},
   {   1,   3,   1,   9,   5,  21, 119,   0},
   {   1,   1,   3,   1,  23,  13,  75,   0},
   {   1,   3,   3,  11,  27,  31,  73,   0},
   {   1,   1,   7,   7,  19,  25, 105,   0},
   {   1,   3,   5,   5,  21,   9,   7,   0},
   {   1,   1,   1,  15,   5,  49,  59,   0},
   {   1,   1,   1,   1,   1,  33,  65,   0},
   {   1,   3,   5,  15,  17,  19,  21,   0},
   {   1,   1,   7,  11,  13,  29,   3,   0},
   {   1,   3,   7,   5,   7,  11, 113,   0},
   {   1,   1,   5,   3,  15,  19,  61,   0},
   {   1,   3,   1,   1,   9,  27,  89,   7},
   {   1,   1,   3,   7,  31,  15,  45,  23},
   {   1,   3,   3,   9,   9,  25, 107,  39}
};


// Prototypes

class SobolM : public  HeapAllocatedGeneratorMatrix2<u64>
{
public:
   SobolM (unsigned _dim, unsigned _m, unsigned _vec);
};


/**
 *  Main function
 */

int main (void)
{
   try
   {
      SobolM m (SobolMatrix::MAX_DIM,
                SobolMatrix::MAX_LOG_N,
                SobolMatrix::TOTAL_PREC);

      cout << "/*******************************************************/\n"
              "/***   This file is program-generated!               ***/\n"
              "/***                                                 ***/\n"
              "/***   Do not change!!!                              ***/\n"
              "/***                                                 ***/\n"
              "/***   Update " __FILE__ " to update this file.  ***/\n"
              "/*******************************************************/\n"
              "\n"
              "#ifdef __GNUG__\n"
              "#pragma implementation\n"
              "#endif\n"
              "\n"
              "#include <HIntLib/sobolmatrix.h>\n"
              "\n"
              "namespace L = HIntLib;\n"
              "\n"
              "const L::u64 L::SobolMatrix::v_mem [MAX_LOG_N][MAX_DIM] =\n";

      m.CArrayDump (cout);

      cout << "\n"
              "const unsigned L::SobolMatrix::t_s [MAX_DIM] =\n"
              "{\n";

      unsigned sum = 0;

      for (unsigned i = 0; i < SobolMatrix::MAX_DIM; ++i)
      {
         sum += polynomials [i].degree() - 1;
   
         cout << setw(6) << sum << ",   // t_s [" << i + 1 << "]\n";
      }

      cout << "};\n\n";
 
      return 0;
   } 
   catch (std::exception &e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      throw;
   }
}



/**
 *   Calculation of vectors v according to Bratley-Fox
 */

SobolM::SobolM (unsigned _dim, unsigned _m, unsigned _totalPrec)
   : HeapAllocatedGeneratorMatrix2<u64> (_dim, _m, _totalPrec)
{
   // Initialize v for all possible dimensions

   for (unsigned d = 0; d < dim; ++d)
   {
      P p = polynomials [d];

      const int deg = p.degree();

      // Initializ the first elements of v from the table

      for (int i = 0; i < deg; ++i)
      {
         setv(d,i, vinit [d][i] << (totalPrec - 1 - i));
      }

      // Calculate the following elements according to the recurrency

      for (unsigned i = deg; i < m; ++i)
      {
         BaseType x = operator()(d, i-deg) >> deg;

         for (int j = 1; j <= deg; ++j)
         {
            if (p[deg - j])  x ^= operator()(d, i-j);
         }

         setv (d, i, x);
      }
   }  //  for dim
}

