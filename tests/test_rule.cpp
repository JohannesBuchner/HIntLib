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
 *  test_rule
 *
 *  Performes various tests with CubatureRule and EmbeddedRule
 *
 *      Ensures correctness for rules
 */

#include <iostream>
#include <iomanip>
#include <memory>
#include <stdlib.h>

#include <HIntLib/cubaturerule.h>
#include <HIntLib/make.h>
#include <HIntLib/testfunction.h>
#include <HIntLib/mymath.h>
#include <HIntLib/mersennetwister.h>
#include <HIntLib/distribution.h>

#include "test.h"

using namespace std;
using namespace HIntLib;

MersenneTwister mt;


/**
 *  User-tuneable constants
 */

const real ACCEPTABLE_ERROR = 2000.0;
      Index    MAX_NUM_POINTS = 3000;
      unsigned NUM_TESTS = 2;
const unsigned FIRST_TEST_DIM =  1;
      unsigned  LAST_TEST_DIM = 100;

      unsigned FIRST_RULE = 1;
      unsigned  LAST_RULE = 100;



/**
 *  TestPolynomial
 *
 *  A TestFunction, modelling a random multi-variate polynomial with a certain
 *  maximal degree.
 */

class TestPolynomial : public TestFunction
{
public:

   TestPolynomial (unsigned dim, unsigned degree);
   ~TestPolynomial ()  { delete [] exponents; }

   real getExactResult (const Hypercube &) const;
   real operator() (const real *);

private:

   // Static member with constructor used to initialize static data

   static
   struct Dummy
   {
      Dummy ();
   } dummy;

   // Number of monomials in each polynomial

   static const unsigned NUM_MON = 130;

   const unsigned degree;

   real coef [NUM_MON];
   unsigned char *exponents;   // [NUM_MON][dim]

   unsigned char& getExponent (unsigned i, unsigned d)
      { return exponents [d + i * dim]; }
   unsigned char getExponent (unsigned i, unsigned d) const
      { return exponents [d + i * dim]; }

   /*
    * Static monomial table
    */

   static const unsigned MAX_DEGREE = 10; // Highest degree in use

   static const unsigned NUM_MONOMIALS;   // Size of the following three tables

   static const unsigned char monomials [][MAX_DEGREE]; // monomials
   static unsigned monomialsDim [];                     // # differnt variables
   static unsigned monomialsDeg [];                     // degree

   // Formatted output

   friend ostream& operator<< (ostream&, const TestPolynomial&);
};


const unsigned TestPolynomial::MAX_DEGREE;


/**
 *  List with all polynomials up to a certain degree, without listing
 *  permutations of variables twice.
 *
 *  Each entry contains the exponents of x_1,...,x_MAX_DEGREE
 */

/**
 *  The following data ensures that no monomials have been forgotten in the
 *  table.
 *
 *  Number of entries for a given degree n equals the number of unordered
 *  partitions x_1+...+x_k=n of the number n, denoted as p(n).
 *
 *  This number can be calculated using the numbers
 *      p_k(n) := #{unordered partitions of n with <= k summands
 *
 *  Using p_k(n), p(n) can be calculated as p_n(n).
 *
 *  p_k(n) can be calculated by the recursion
 *
 *      p_k (n) = p_k(n-k) + p_(k-1)(n)
 *      p_1 (n) = 1
 *
 *  n  k 1  2  3  4  5  6  7  8  9 10
 *   0   1
 *   1   1
 *   2   1  2
 *   3   1  2  3
 *   4   1  3  4  5
 *   5   1  3  5  6  7
 *   6   1  4  7  9 10 11
 *   7   1  4  8 11 13 14 15
 *   8   1  5 10 15 18 20 21 22
 *   9   1  5 12 18 23 26 28 29 30
 *  10   1  6 14 23 30 35 38 40 41 42
 *  11   1  6 16 27 37 44 49 52 54 55 56
 *
 *  More details can be found in my script "Kombinatorik" from R.Wolf.
 */

const unsigned char TestPolynomial::monomials [][MAX_DEGREE] =
{
   { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},   // 1          deg 0   #=1
   { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},   // x          deg 1   #=1
   { 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},   // x²         deg 2   #=2
   { 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},   // xy
   { 3, 0, 0, 0, 0, 0, 0, 0, 0, 0},   // x³         deg 3   #=3
   { 2, 1, 0, 0, 0, 0, 0, 0, 0, 0},   // x²y
   { 1, 1, 1, 0, 0, 0, 0, 0, 0, 0},   // xyz
   { 4, 0, 0, 0, 0, 0, 0, 0, 0, 0},   // x^4        deg 4   #=5
   { 3, 1, 0, 0, 0, 0, 0, 0, 0, 0},   // x³y
   { 2, 2, 0, 0, 0, 0, 0, 0, 0, 0},   // x²y²
   { 2, 1, 1, 0, 0, 0, 0, 0, 0, 0},   // x²yz
   { 1, 1, 1, 1, 0, 0, 0, 0, 0, 0},   // xyzu
   { 5, 0, 0, 0, 0, 0, 0, 0, 0, 0},   // x^5        deg 5   #=7
   { 4, 1, 0, 0, 0, 0, 0, 0, 0, 0},   // x^4 y
   { 3, 2, 0, 0, 0, 0, 0, 0, 0, 0},   // x³y²
   { 3, 1, 1, 0, 0, 0, 0, 0, 0, 0},   // x³yz
   { 2, 2, 1, 0, 0, 0, 0, 0, 0, 0},   // x²y²z
   { 2, 1, 1, 1, 0, 0, 0, 0, 0, 0},
   { 1, 1, 1, 1, 1, 0, 0, 0, 0, 0},
   { 6, 0, 0, 0, 0, 0, 0, 0, 0, 0},   //            deg 6   #=11
   { 5, 1, 0, 0, 0, 0, 0, 0, 0, 0},
   { 4, 2, 0, 0, 0, 0, 0, 0, 0, 0},
   { 4, 1, 1, 0, 0, 0, 0, 0, 0, 0},
   { 3, 3, 0, 0, 0, 0, 0, 0, 0, 0},
   { 3, 2, 1, 0, 0, 0, 0, 0, 0, 0},
   { 3, 1, 1, 1, 0, 0, 0, 0, 0, 0},
   { 2, 2, 2, 0, 0, 0, 0, 0, 0, 0},
   { 2, 2, 1, 1, 0, 0, 0, 0, 0, 0},
   { 2, 1, 1, 1, 1, 0, 0, 0, 0, 0},
   { 1, 1, 1, 1, 1, 1, 0, 0, 0, 0},
   { 7, 0, 0, 0, 0, 0, 0, 0, 0, 0},  // deg 7   #=15
   { 6, 1, 0, 0, 0, 0, 0, 0, 0, 0},
   { 5, 2, 0, 0, 0, 0, 0, 0, 0, 0},
   { 5, 1, 1, 0, 0, 0, 0, 0, 0, 0},
   { 4, 3, 0, 0, 0, 0, 0, 0, 0, 0},
   { 4, 2, 1, 0, 0, 0, 0, 0, 0, 0},
   { 4, 1, 1, 1, 0, 0, 0, 0, 0, 0},
   { 3, 3, 1, 0, 0, 0, 0, 0, 0, 0},
   { 3, 2, 2, 0, 0, 0, 0, 0, 0, 0},
   { 3, 2, 1, 1, 0, 0, 0, 0, 0, 0},
   { 3, 1, 1, 1, 1, 0, 0, 0, 0, 0},
   { 2, 2, 2, 1, 0, 0, 0, 0, 0, 0},
   { 2, 2, 1, 1, 1, 0, 0, 0, 0, 0},
   { 2, 1, 1, 1, 1, 1, 0, 0, 0, 0},
   { 1, 1, 1, 1, 1, 1, 1, 0, 0, 0},
   { 8, 0, 0, 0, 0, 0, 0, 0, 0, 0},  // deg 8   #=22
   { 7, 1, 0, 0, 0, 0, 0, 0, 0, 0},
   { 6, 2, 0, 0, 0, 0, 0, 0, 0, 0},
   { 6, 1, 1, 0, 0, 0, 0, 0, 0, 0},
   { 5, 3, 0, 0, 0, 0, 0, 0, 0, 0},
   { 5, 2, 1, 0, 0, 0, 0, 0, 0, 0},
   { 5, 1, 1, 1, 0, 0, 0, 0, 0, 0},
   { 4, 4, 0, 0, 0, 0, 0, 0, 0, 0},
   { 4, 3, 1, 0, 0, 0, 0, 0, 0, 0},
   { 4, 2, 2, 0, 0, 0, 0, 0, 0, 0},
   { 4, 2, 1, 1, 0, 0, 0, 0, 0, 0},
   { 4, 1, 1, 1, 1, 0, 0, 0, 0, 0},
   { 3, 3, 2, 0, 0, 0, 0, 0, 0, 0},
   { 3, 3, 1, 1, 0, 0, 0, 0, 0, 0},
   { 3, 2, 2, 1, 0, 0, 0, 0, 0, 0},
   { 3, 2, 1, 1, 1, 0, 0, 0, 0, 0},
   { 3, 1, 1, 1, 1, 1, 0, 0, 0, 0},
   { 2, 2, 2, 2, 0, 0, 0, 0, 0, 0},
   { 2, 2, 2, 1, 1, 0, 0, 0, 0, 0},
   { 2, 2, 1, 1, 1, 1, 0, 0, 0, 0},
   { 2, 1, 1, 1, 1, 1, 1, 0, 0, 0},
   { 1, 1, 1, 1, 1, 1, 1, 1, 0, 0},
   { 9, 0, 0, 0, 0, 0, 0, 0, 0, 0},  // deg 9   #=30
   { 8, 1, 0, 0, 0, 0, 0, 0, 0, 0},
   { 7, 2, 0, 0, 0, 0, 0, 0, 0, 0},
   { 7, 1, 1, 0, 0, 0, 0, 0, 0, 0},
   { 6, 3, 0, 0, 0, 0, 0, 0, 0, 0},
   { 6, 2, 1, 0, 0, 0, 0, 0, 0, 0},
   { 6, 1, 1, 1, 0, 0, 0, 0, 0, 0},
   { 5, 4, 0, 0, 0, 0, 0, 0, 0, 0},
   { 5, 3, 1, 0, 0, 0, 0, 0, 0, 0},
   { 5, 2, 2, 0, 0, 0, 0, 0, 0, 0},
   { 5, 2, 1, 1, 0, 0, 0, 0, 0, 0},
   { 5, 1, 1, 1, 1, 0, 0, 0, 0, 0},
   { 4, 4, 1, 0, 0, 0, 0, 0, 0, 0},
   { 4, 3, 2, 0, 0, 0, 0, 0, 0, 0},
   { 4, 3, 1, 1, 0, 0, 0, 0, 0, 0},
   { 4, 2, 2, 1, 0, 0, 0, 0, 0, 0},
   { 4, 2, 1, 1, 1, 0, 0, 0, 0, 0},
   { 4, 1, 1, 1, 1, 1, 0, 0, 0, 0},
   { 3, 3, 3, 0, 0, 0, 0, 0, 0, 0},
   { 3, 3, 2, 1, 0, 0, 0, 0, 0, 0},
   { 3, 3, 1, 1, 1, 0, 0, 0, 0, 0},
   { 3, 2, 2, 2, 0, 0, 0, 0, 0, 0},
   { 3, 2, 2, 1, 1, 0, 0, 0, 0, 0},
   { 3, 2, 1, 1, 1, 1, 0, 0, 0, 0},
   { 3, 1, 1, 1, 1, 1, 1, 0, 0, 0},
   { 2, 2, 2, 2, 1, 0, 0, 0, 0, 0},
   { 2, 2, 2, 1, 1, 1, 0, 0, 0, 0},
   { 2, 2, 1, 1, 1, 1, 1, 0, 0, 0},
   { 2, 1, 1, 1, 1, 1, 1, 1, 0, 0},
   { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0},
   {10, 0, 0, 0, 0, 0, 0, 0, 0, 0},  // deg 10   #=42
   { 9, 1, 0, 0, 0, 0, 0, 0, 0, 0},
   { 8, 2, 0, 0, 0, 0, 0, 0, 0, 0},
   { 8, 1, 1, 0, 0, 0, 0, 0, 0, 0},
   { 7, 3, 0, 0, 0, 0, 0, 0, 0, 0},
   { 7, 2, 1, 0, 0, 0, 0, 0, 0, 0},
   { 7, 1, 1, 1, 0, 0, 0, 0, 0, 0},
   { 6, 4, 0, 0, 0, 0, 0, 0, 0, 0},
   { 6, 3, 1, 0, 0, 0, 0, 0, 0, 0},
   { 6, 2, 2, 0, 0, 0, 0, 0, 0, 0},
   { 6, 2, 1, 1, 0, 0, 0, 0, 0, 0},
   { 6, 1, 1, 1, 1, 0, 0, 0, 0, 0},
   { 5, 5, 0, 0, 0, 0, 0, 0, 0, 0},
   { 5, 4, 1, 0, 0, 0, 0, 0, 0, 0},
   { 5, 3, 2, 0, 0, 0, 0, 0, 0, 0},
   { 5, 3, 1, 1, 0, 0, 0, 0, 0, 0},
   { 5, 2, 2, 1, 0, 0, 0, 0, 0, 0},
   { 5, 2, 1, 1, 1, 0, 0, 0, 0, 0},
   { 5, 1, 1, 1, 1, 1, 0, 0, 0, 0},
   { 4, 4, 2, 0, 0, 0, 0, 0, 0, 0},
   { 4, 4, 1, 1, 0, 0, 0, 0, 0, 0},
   { 4, 3, 3, 0, 0, 0, 0, 0, 0, 0},
   { 4, 3, 2, 1, 0, 0, 0, 0, 0, 0},
   { 4, 3, 1, 1, 1, 0, 0, 0, 0, 0},
   { 4, 2, 2, 2, 0, 0, 0, 0, 0, 0},
   { 4, 2, 2, 1, 1, 0, 0, 0, 0, 0},
   { 4, 2, 1, 1, 1, 1, 0, 0, 0, 0},
   { 4, 1, 1, 1, 1, 1, 1, 0, 0, 0},
   { 3, 3, 3, 1, 0, 0, 0, 0, 0, 0},
   { 3, 3, 2, 2, 0, 0, 0, 0, 0, 0},
   { 3, 3, 2, 1, 1, 0, 0, 0, 0, 0},
   { 3, 3, 1, 1, 1, 1, 0, 0, 0, 0},
   { 3, 2, 2, 2, 1, 0, 0, 0, 0, 0},
   { 3, 2, 2, 1, 1, 1, 0, 0, 0, 0},
   { 3, 2, 1, 1, 1, 1, 1, 0, 0, 0},
   { 3, 1, 1, 1, 1, 1, 1, 1, 0, 0},
   { 2, 2, 2, 2, 2, 0, 0, 0, 0, 0},
   { 2, 2, 2, 2, 1, 1, 0, 0, 0, 0},
   { 2, 2, 2, 1, 1, 1, 1, 0, 0, 0},
   { 2, 2, 1, 1, 1, 1, 1, 1, 0, 0},
   { 2, 1, 1, 1, 1, 1, 1, 1, 1, 0},
   { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
};

const unsigned TestPolynomial::NUM_MONOMIALS
     = sizeof (monomials) / sizeof (monomials [0]);

unsigned TestPolynomial::monomialsDim [NUM_MONOMIALS];
unsigned TestPolynomial::monomialsDeg [NUM_MONOMIALS];
TestPolynomial::Dummy TestPolynomial::dummy;


/**
 *  Initialize static members monomialsDeg[] and monomialsDim[]
 */

TestPolynomial::Dummy::Dummy ()
{
   for (unsigned i = 0; i < NUM_MONOMIALS; ++i)
   {
      unsigned sum = 0;
      unsigned dim = 0;

      for (unsigned j = 0; j < MAX_DEGREE; ++j)
      {
         if (monomials [i][j])
         {
            sum += monomials [i][j];
            dim = j + 1;
         }
      }

      monomialsDim [i] = dim;
      monomialsDeg [i] = sum;
   }
}

TestPolynomial::TestPolynomial (unsigned dim, unsigned degree)
   : TestFunction (dim), degree (degree),
     exponents (new unsigned char [NUM_MON * dim])
{
   unsigned int m = 0;

   for (unsigned i = 0; i < NUM_MON; ++i)
   {
      coef [i] = uniform (mt, -10.0, 10.0);  // Choose koefficient from [-10,10]

      // Copy monomomial

      for (unsigned int j = 0; j < min (dim, MAX_DEGREE); ++j)
      {
         getExponent(i, j) = monomials [m][j];
      }

      // Fill remaining exponents with 0

      for (unsigned int j = MAX_DEGREE; j < dim; ++j)  getExponent(i, j) = 0;

      // Shuffle Exponents

      random_shuffle (&getExponent (i,0), &getExponent (i,dim), mt);

      // Find next monomial with acceptable degree and dimension

      do
      {
         m = (m + 1) % NUM_MONOMIALS;
      }
      while (monomialsDeg [m] > degree || monomialsDim [m] > dim);
   }
}

real TestPolynomial::getExactResult (const Hypercube &h) const
{
   real sum = 0.0;

   for (unsigned i = 0; i < NUM_MON; ++i)
   {
      real prod = coef [i];

      for (unsigned d = 0; d < dim; ++d)
      {
         unsigned e = getExponent (i, d) + 1;

         prod *= (  powInt(h.getUpperBound (d), e)
                  - powInt(h.getLowerBound (d), e)) / e;
      }

      sum += prod;
   }

   return sum;
}

/**
 *  Evaluate polynomail on a given point x
 */

real TestPolynomial::operator() (const real *x)
{
   real sum = 0.0;

   for (unsigned i = 0; i < NUM_MON; ++i)   // all monomials in polynomial
   {
      real prod = coef [i];

      for (unsigned d = 0; d < dim; ++d)
      {
         for (unsigned j = getExponent (i, d); j > 0; --j)  prod *= x [d];
      }

      sum += prod;
   }

   return sum;
}

ostream& operator<< (ostream& os, const TestPolynomial &p)
{
   for (unsigned i = 0; i < TestPolynomial::NUM_MON; ++i)
   {
      os << (p.coef [i] > 0.0 ? " + " : " - ") << abs (p.coef [i]);

      for (unsigned d = 0; d < p.dim; ++d)
      {
         if (p.getExponent (i, d) > 0)
         {
            os << " x";

            if (p.dim > 1)  os << '_' << d;

            switch (p.getExponent (i, d))
            {
            case 1: break;
            case 2: os << '²'; break;
            case 3: os << '³'; break;
            default: os << '^' << int (p.getExponent (i, d));
            }
         }
      }
   }

   return os;
}


bool checkDegree (CubatureRule& r, unsigned degree, unsigned numTests)
{
   const unsigned dim = r.getDimension ();

   const Hypercube h (dim, -0.2, 1.7);

   for (unsigned i = 0; i < numTests; ++i)
   {
      TestPolynomial p (r.getDimension (), degree);

      DEB2  cout << "Using polynomial " << p << endl;

      EvalCounterFunction f1 (&p);
      DomainCheckerFunction f (&f1, &h);

      real result = r.eval (f, h);
      real correct = f.getExactResult (h);

      if (! f.isAllPointsInside() && r.isAllPointsInside())
      {
         error ("\nAbscissas outside the cube encountered");
      }
      else if (f.isAllPointsInside() && ! r.isAllPointsInside())
      {
         error ("\nisAllPointsInside() should return true");
      }

      if (f1.getCounter () != r.getNumPoints ())
      {
         cout << "\ngetNumPoints() wrong ("
              << r.getNumPoints () << "!=" << f1.getCounter ()
              << ")!!!" << endl;
         error();
      }

// cout << setprecision (15);
// cout << setw (20) << result << setw (20) << correct << endl;

      if (  abs ((result - correct) / max(abs(correct), real(1.0)))
          >  numeric_limits<real>::epsilon() * ACCEPTABLE_ERROR
            * (r.getSumAbsWeight () + real (dim)))
      {
         DEB1
         {
            cout << "Error factor: "
                 << abs((result - correct) / max(abs(correct), real(1.0)))
                    / (numeric_limits<real>::epsilon()
                        * (r.getSumAbsWeight () + real (dim)))
                 << " exceeds "
                 << ACCEPTABLE_ERROR
                 << endl;
         }

         return false;
      }
   }

   return true;
}

void testRule (CubatureRule& r, unsigned dim)
{
   if (dim != r.getDimension())
   {
      cout << "\ngetDimension() is broken (" << dim << "!=" << r.getDimension()
           << ")!!!" << endl;
      error();
   }

   Index num_points = r.getNumPoints();

   if (num_points > MAX_NUM_POINTS)
   {
      NORMAL  cout << "Rule has " << num_points << " (> " << MAX_NUM_POINTS
                   << ") abscissas.  Test skipped." << endl;
      return;
   }
   else if (verbose == 0)  cout << ' ' << dim << flush;

   unsigned degree = r.getDegree ();

   bool allright = true;

   for (unsigned d = 0; d <= degree; ++d)
   {
      if (checkDegree (r, d, NUM_TESTS))
      {
         NORMAL  cout << "Degree " << d << " passed." << endl;
      }
      else
      {
         cout << "\nDegree " << d << " failed!!!" << endl;
         error();
         allright = false;
      }
   }

   if (allright)
   {
      if (checkDegree (r, degree + 1, 2 * NUM_TESTS + 5))
      {
         error ("\nDegree seems to be higher");
      }
      else
      {
         NORMAL  cout << "Degree seems to be correct." << endl;
      }
   }
}

const char* options = "n:p:d:r:";

bool opt(int c, const char* s)
{
   switch (c)
   {
      case 'n':  NUM_TESTS      = atoi (s); return true;
      case 'p':  MAX_NUM_POINTS = atoi (s); return true;
      case 'd':  LAST_TEST_DIM  = atoi (s); return true;
      case 'r':  LAST_RULE = FIRST_RULE = atoi (s); return true;
   }

   return false;
}

void usage()
{
   cerr <<
      "Usage: test_rule [OPTION]...\n\n"
      "Tests if cubature rule work correctly.\n\n"
      << option_msg <<
      "  -n n   Number of testruns to be performed (default = 2)\n"
      "  -p n   Tests requiring more that _p_ points are skipped"
                     " (default = 3000)\n"
      "  -d n   Maximal dimension (default = 100)\n"
      "  -r n   Tests only rule number _n_\n"
      "\n";

   exit (1);
}

void test(int argc, char**)
{
   if (argc)  usage();
   for (unsigned rule = FIRST_RULE; rule <= LAST_RULE; ++rule)
   {
      const char* name;
      try
      {
         name = Make::getCubatureRuleFactoryName (rule);
      }
      catch (Make::CubatureRuleDoesNotExist &)
      {
         continue;
      }
      
      cout << "Testing Rule #" << rule << " (" << name << ") ..." << endl;

      for (unsigned dim = FIRST_TEST_DIM; dim <= LAST_TEST_DIM; ++dim)
      {
         NORMAL  cout << "Dimension: " << dim << endl;

         try
         {
            auto_ptr<CubatureRuleFactory> f (Make::cubatureRuleFactory(rule));

            auto_ptr<CubatureRule> r (f->create(dim));
            testRule (*r, dim);
         }
         catch (InvalidDimension &e)
         {
            NORMAL  cout << "Can not create rule in dimension "
                         << e.getDimension() << "!" << endl;
         }
      }

      if (verbose == 0)  cout << endl;
   }
}


