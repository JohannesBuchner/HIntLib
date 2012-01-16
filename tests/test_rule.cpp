/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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
 *  Ensures correctness for rules
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

#include <HIntLib/cubaturerule.h>
#include <HIntLib/make.h>
#include <HIntLib/testintegrand.h>
#include <HIntLib/hlmath.h>
#include <HIntLib/mersennetwister.h>
#include <HIntLib/distribution.h>
#include <HIntLib/output.h>

#include "test.h"

using namespace std;
using namespace HIntLib;

MersenneTwister mt;

const char* epsilonString;
const char* deltaString;

/**
 *  User-tuneable constants
 */

      Index    MAX_NUM_POINTS = 3000; // do not check rules with more points
      unsigned NUM_TESTS = 3;         // number of test-polynomials / dimesnion
const unsigned FIRST_TEST_DIM =  1;
      unsigned  LAST_TEST_DIM = 100;
      unsigned FIRST_RULE = 1;
      unsigned  LAST_RULE = 100;
const unsigned NUM_MON = 130;        // Number of monomials in each polynomial
const real     MAX_COEFF = 10.0;     // Maximal coefficient in each polynomial


/**
 *  epsilonPower ()
 *
 *  Returns  (1+epsilon)^power , with epsilong the accuracy of real.
 */

real epsilonPowerMinusOne (int power)
{
   const real e = numeric_limits<real>::epsilon();

   switch (power)
   {
   case 0:  return 0;
   case 1:  return e;
   case 2:  return 2 * e;
   default:
      return
         power * e
       + (power * (power-1) / 2) * e * e
       + (power * (power-1) * (power-2) / 6) * e * e *e;
   }
}

inline
real epsilonPower (int power)
{
   return epsilonPowerMinusOne (power) + real(1);
}


/**
 *  TestPolynomial
 *
 *  A TestIntegrand, modelling a random multi-variate polynomial with a certain
 *  maximal degree.
 */

class TestPolynomial : public TestIntegrand
{
public:

   TestPolynomial (unsigned dim, unsigned degree);

   real operator() (const real *);
   real derivative (const real *, unsigned);

   real getExactResult (const Hypercube &) const;
   real getPartialResult (const Hypercube &, unsigned) const;
   real getErrorOfExactResult (const Hypercube&) const;
   real getErrorOfEvaluation (const Hypercube&) const;

private:

   // Static member with constructor used to initialize static data

   static
   struct Dummy
   {
      Dummy ();
   } dummy;

   const unsigned degree;

   real coef [NUM_MON];

   Array<unsigned char> exponents;  // [NUM_MON][dim]

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
 *  More details can be found in the lecture notes "Kombinatorik" from R.Wolf.
 */

const unsigned char TestPolynomial::monomials [][MAX_DEGREE] =
{
   { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},   // 1          deg 0   #=1
   { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},   // x          deg 1   #=1
   { 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},   // x^2        deg 2   #=2
   { 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},   // x y
   { 3, 0, 0, 0, 0, 0, 0, 0, 0, 0},   // x^3        deg 3   #=3
   { 2, 1, 0, 0, 0, 0, 0, 0, 0, 0},   // x^2 y
   { 1, 1, 1, 0, 0, 0, 0, 0, 0, 0},   // x y z
   { 4, 0, 0, 0, 0, 0, 0, 0, 0, 0},   // x^4        deg 4   #=5
   { 3, 1, 0, 0, 0, 0, 0, 0, 0, 0},   // x^3 y
   { 2, 2, 0, 0, 0, 0, 0, 0, 0, 0},   // x^2 y^2
   { 2, 1, 1, 0, 0, 0, 0, 0, 0, 0},   // x^2 y z
   { 1, 1, 1, 1, 0, 0, 0, 0, 0, 0},   // x y z u
   { 5, 0, 0, 0, 0, 0, 0, 0, 0, 0},   // x^5        deg 5   #=7
   { 4, 1, 0, 0, 0, 0, 0, 0, 0, 0},   // x^4 y
   { 3, 2, 0, 0, 0, 0, 0, 0, 0, 0},   // x^3 y^2
   { 3, 1, 1, 0, 0, 0, 0, 0, 0, 0},   // x^3 y z
   { 2, 2, 1, 0, 0, 0, 0, 0, 0, 0},   // x^2 y^2 z
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
         if (! monomials [i][j])  break;

         sum += monomials [i][j];
         dim = j + 1;
      }

      monomialsDim [i] = dim;
      monomialsDeg [i] = sum;
   }
}


/**
 *  Constructor
 */

TestPolynomial::TestPolynomial (unsigned dim, unsigned degree)
   : TestIntegrand (dim), degree (degree), exponents (NUM_MON * dim)
{
   unsigned int m = 0;

   for (unsigned i = 0; i < NUM_MON; ++i)
   {
      coef [i] = uniform (mt, -MAX_COEFF, MAX_COEFF);

      // Copy monomomial

      for (unsigned j = 0; j < min (dim, MAX_DEGREE); ++j)
      {
         getExponent(i, j) = monomials [m][j];
      }

      // Fill remaining exponents with 0

      for (unsigned j = MAX_DEGREE; j < dim; ++j)  getExponent(i, j) = 0;

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


/**
 *  getExactResult()
 */

real TestPolynomial::getExactResult (const Hypercube &h) const
{
   real sum = 0;

   for (unsigned i = 0; i < NUM_MON; ++i)
   {
      real prod = coef[i];

      for (unsigned d = 0; d < dim; ++d)
      {
         int e = getExponent (i, d) + 1;

         prod *= (e == 1)
            ? h.getDiameter (d)
            : (  HINTLIB_MN pow (h.getUpperBound (d), e)
               - HINTLIB_MN pow (h.getLowerBound (d), e)) / real (e);
      }

      sum += prod;
   }

   return sum;
}


/**
 *  getPartialResult()
 *
 *  Returns the integral, considering only monimials with degree _deg_.
 */

real TestPolynomial::getPartialResult (const Hypercube &h, unsigned deg) const
{
   real sum = 0;

   for (unsigned i = 0; i < NUM_MON; ++i)
   {
      real prod = coef [i];

      unsigned thisDegree = 0;
      for (unsigned d = 0; d < dim; ++d)  thisDegree += getExponent (i, d);
      if (thisDegree != deg)  continue;

      for (unsigned d = 0; d < dim; ++d)
      {
         int e = getExponent (i, d) + 1;

         prod *= (e == 1)
            ? h.getDiameter (d)
            : (  HINTLIB_MN pow (h.getUpperBound (d), e)
               - HINTLIB_MN pow (h.getLowerBound (d), e)) / real (e);
      }

      sum += prod;
   }

   return sum;
}


/**
 *  getErrorOfExactResult ()
 *
 *  Returns an upper bound on the rounding error in  getExactResult()
 */

real TestPolynomial::getErrorOfExactResult (const Hypercube &h) const
{
   const real onePlusEpsilon = 1.0 + numeric_limits<real>::epsilon();
   real sum = 0;

   for (unsigned i = 0; i < NUM_MON; ++i)
   {
      real realProd = 1;
      real prod     = onePlusEpsilon;
      real factor = coef [i];

      for (unsigned d = 0; d < dim; ++d)
      {
         int e = getExponent (i, d) + 1;

         if (e == 1)
         {
            factor *= h.getDiameter(d);
            prod   *= sqr (onePlusEpsilon);
         }
         else
         {
            real v1 = HINTLIB_MN pow (h.getUpperBound (d), e);
            real v2 = HINTLIB_MN pow (h.getLowerBound (d), e);

            realProd *= v1 - v2;

            real ePower = epsilonPower (e);

            if (v1 > v2)
            {
               v1 *= (v1 > 0.0) ?       ePower : 1.0 / ePower;
               v2 *= (v2 > 0.0) ? 1.0 / ePower :       ePower;
            }
            else
            {
               v1 *= (v1 > 0.0) ? 1.0 / ePower :       ePower;
               v2 *= (v2 > 0.0) ?       ePower : 1.0 / ePower;
            }

            prod *= v1 - v2;
            
            factor /= real (e);
         }
      }

      sum += L::abs (realProd - prod) * L::abs (factor);
   }

   return sum;
}


/**
 *  Evaluate polynomail on a given point x
 */

real TestPolynomial::operator() (const real *x)
{
   real sum = 0;

   for (unsigned i = 0; i < NUM_MON; ++i)   // all monomials in polynomial
   {
      real prod = coef [i];

      for (unsigned d = 0; d < dim; ++d)
      {
         const int exponent = getExponent (i,d);
         if (exponent)  prod *= HINTLIB_MN pow (x[d], exponent);
      }

      sum += prod;
   }

   return sum;
}


/**
 *  getErrorOfEvaluation ()
 *
 *  Returns an upper bound the round-off error encountered during a integrand
 *  evaluation
 */

real TestPolynomial::getErrorOfEvaluation (const Hypercube &h) const
{
   Array<real> maxCoordinate (h.getDimension());

   for (unsigned i = 0; i < h.getDimension(); ++i)
   {
      maxCoordinate [i]
         = max (L::abs (h.getLowerBound(i)), L::abs (h.getUpperBound(i)));
   }

   real sum = 0.0;

   for (unsigned i = 0; i < NUM_MON; ++i)   // all monomials in polynomial
   {
      real prod = L::abs (coef [i]);

      unsigned degreeOfMonomial = 0;

      for (unsigned d = 0; d < dim; ++d)
      {
         degreeOfMonomial += getExponent (i, d);
         prod *= HINTLIB_MN pow (maxCoordinate [d], getExponent (i, d));
      }

      sum += prod * (epsilonPowerMinusOne (degreeOfMonomial + 1));
   }

   return sum;
}

real TestPolynomial::derivative (const real *x, unsigned a)
{
   real sum = 0.0;

   for (unsigned i = 0; i < NUM_MON; ++i)   // all monomials in polynomial
   {
      if (getExponent (i, a) == 0)  continue;

      real prod = coef [i];

      for (unsigned d = 0; d < dim; ++d)
      {
         if (d == a)
         {
            prod = prod * getExponent (i, a)
                        * HINTLIB_MN pow (x[d], getExponent (i,a) - 1);
         }
         else
         {
            prod *= HINTLIB_MN pow (x[d], getExponent (i, d));
         }
      }

      sum += prod;
   }

   return sum;
}

#if 0
real TestPolynomial::derivative (const real *x, unsigned a, unsigned b)
{
   real sum = 0.0;

   if (a == b)
   {
      for (unsigned i = 0; i < NUM_MON; ++i)   // all monomials in polynomial
      {
         if (getExponent (i, a) <= 1)  continue;

         real prod = coef [i];

         for (unsigned d = 0; d < dim; ++d)
         {
            if (d == a)
            {
               unsigned exponent = getExponent (i, a);
               prod = exponent * (exponent - 1)
                    * prod * HINTLIB_MN pow (x[d], exponent - 2);
            }
            else
            {
               prod *= HINTLIB_MN pow (x[d], getExponent (i, d));
            }
         }

         sum += prod;
      }
   }
   else   // a != b
   {
      for (unsigned i = 0; i < NUM_MON; ++i)   // all monomials in polynomial
      {
         if (getExponent (i, a) == 0)  continue;
         if (getExponent (i, b) == 0)  continue;

         real prod = coef [i];

         for (unsigned d = 0; d < dim; ++d)
         {
            if (d == a || d == b)
            {
               prod = prod * getExponent (i,d)
                           * HINTLIB_MN pow (x[d], getExponent (i,d) - 1);
            }
            else
            {
               prod *= HINTLIB_MN pow (x[d], getExponent (i, d));
            }
         }

         sum += prod;
      }

   }

   return sum;
}
#endif

ostream& operator<< (ostream& os, const TestPolynomial &p)
{
   HIntLib::Private::Printer ss (os);
   ss << setprecision (3);

   for (unsigned i = 0; i < NUM_MON; ++i)
   {
      ss << ' ';
      if (p.coef [i] > 0.0)  ss << '+';
      else ss.minusSign();
      ss << ' ' << L::abs (p.coef [i]);
      bool spacePrinted = false;

      for (unsigned d = 0; d < p.dim; ++d)
      {
         if (p.getExponent (i, d) > 0)
         {
            if (! spacePrinted)
            {
               ss << ' ';
               spacePrinted = true;
            }

            ss << 'x';

            if (p.dim > 1)  ss.subscript (d);

            const int exponent = p.getExponent(i,d);

            if (exponent > 1)  ss.power (exponent);
         }
      }
   }

   return os;
}


/**
 *  check Degree ()
 *
 *  Deterimines, if a given CubatureRule _r_ has at least a certain degree
 *  _degree_.
 *
 *  _numTests_ random polynomials of degree _degree_ are generated and the
 *  approximation produced by _r_ is compared to the exact result.
 *
 *  The main difficulaty is to distinguish between round-off errors and real
 *  failures in the degree of the CubatureRule.
 */

bool checkDegree (
   CubatureRule& r, unsigned degree, unsigned numTests, bool disproof)
{
   const unsigned dim = r.getDimension ();

   const Hypercube h (dim, -0.2, 1.7);

   for (unsigned i = 0; i < numTests; ++i)
   {
      TestPolynomial p (r.getDimension (), degree);

      DEB3  cout << "     Using polynomial " << p << endl;

      EvalCounterIntegrand f1 (&p);
      DomainCheckerIntegrand f (&f1, &h);

      const real result = r.eval (f, h);
      const real correct = f.getExactResult (h);

      // check isAllPointsInside()

      if (! f.isAllPointsInside() && r.isAllPointsInside())
      {
         error ("\nAbscissas outside the cube encountered");
      }
      else if (f.isAllPointsInside() && ! r.isAllPointsInside())
      {
         error ("\nisAllPointsInside() should return true");
      }

      // check getNumPoints()
      
      if (f1.getCounter () != r.getNumPoints ())
      {
         cout << "\ngetNumPoints() wrong ("
              << r.getNumPoints () << " != " << f1.getCounter ()
              << ")!!!" << endl;
         error();
      }

      // check Result()

      const real error = L::abs (result - correct);

      const real expectedExactError = p.getErrorOfExactResult (h);
      const real expectedApproxError
         = p.getErrorOfEvaluation (h) * r.getSumAbsWeight() * h.getVolume()
         * epsilonPower (2 * h.getDimension()) * 2.0;

      const real threshold = expectedExactError + expectedApproxError;
      const real problem = L::abs (p.getPartialResult (h, degree));

      DEB2
      {
         cout.setf (ostream::left, ostream::adjustfield);
         cout << setprecision (10)
              << "     Q(f)=" << setw (13) << result
              << "I(f)="   << setw (13) << correct
              << setprecision (4)
              << deltaString << '=' << setw(10) << error
              << epsilonString << '=' << setw(10) << threshold
              << epsilonString << "_exact=" << setw(10) << expectedExactError
              << epsilonString << "_approx=" << setw(10) << expectedApproxError
              << "probl.=" << problem << endl;
      }
      
      if (error > 25.0 * threshold)
      {
         DEB1 cout << "     Error is " << error << ", exceeding " << threshold
                   << " by a factor of " << error / threshold << ".\n";
         return false;
      }

      if (disproof && threshold > problem / 1000.0)
      {
         DEB1 cout << "Cannot determine degree due to numerical instabilities!"
                   << endl;
         if (verbose == 0)  cout << '-';
         return false;
      }
   }

   return true;
}


/**
 *  test Rule()
 *
 *  Ensures that a given CubatureRule _r_ workes correctly, primarilly by
 *  calling checkDegree().
 */

void testRule (CubatureRule& r, unsigned dim)
{
   // Reset random number generator

   mt.init (12345);

   if (dim != r.getDimension())
   {
      cout << "\ngetDimension() is broken (" << dim << "!=" << r.getDimension()
           << ")!!!" << endl;
      error();
   }

   Index num_points = r.getNumPoints();

   if (num_points > MAX_NUM_POINTS)
   {
      NORMAL  cout << "   Rule has " << num_points << " (> " << MAX_NUM_POINTS
                   << ") abscissas.  Test skipped." << endl;
      return;
   }
   else if (verbose == 0)
   {
      cout << dim << ' ' << flush;
   }

   unsigned degree = r.getDegree ();

   bool allright = true;

   for (unsigned d = 0; d <= degree; ++d)
   {
      if (checkDegree (r, d, NUM_TESTS, false))
      {
         NORMAL  cout << "   Degree " << d << " passed." << endl;
      }
      else
      {
         cout << "\n   Degree " << d << " failed!!!" << endl;
         error();
         allright = false;
      }
   }

   if (allright)
   {
      if (checkDegree (r, degree + 1, 2 * NUM_TESTS + 5, true))
      {
         error ("\n   Degree seems to be higher");
      }
      else
      {
         NORMAL  cout << "   Degree seems to be correct." << endl;
      }
   }
}


/**
 *  Command line arguments
 */

/**
 *  Verbosity levels
 *
 *  -sss  -2    Print nothing at all
 *  -ss   -1    Print only the names of the evaluated rules
 *  -s     0    Names and dimensions
 *  normal 1     + Result for each degree
 *  -v     2     + Error
 *  -vv    3     + Error statistic for each polynomial
 *  -vvv   4     + each evaluated polynomial
 */

const char options[] = "n:p:d:r:";
const char option_msg[] =
   "  -n n   Number of testruns to be performed (default = 3)\n"
   "  -p n   Tests requiring more that _p_ points are skipped"
         " (default = 3000)\n"
   "  -d n   Maximal dimension (default = 100)\n"
   "  -r n   Tests only rule number _n_\n";
const char testProgramParameters[] = "[OPTION]...";
const char testProgramUsage[] =
   "Tests whether cubature rules work correctly.\n\n";
const char testProgramName[] = "test_rule";
const int  testProgramCopyright = 2003;

bool opt(int c, const char* s)
{
   switch (c)
   {
      case 'n':  NUM_TESTS      = HINTLIB_SLN atoi (s); return true;
      case 'p':  MAX_NUM_POINTS = HINTLIB_SLN atoi (s); return true;
      case 'd':  LAST_TEST_DIM  = HINTLIB_SLN atoi (s); return true;
      case 'r':  LAST_RULE = FIRST_RULE = HINTLIB_SLN atoi (s); return true;
   }

   return false;
}


/**
 *  Main program
 */

void test(int argc, char**)
{
   if (argc)  usage("Too many arguments!");

   if (verbose >= 0)  printHeader (cout);

   // GREEK CAPITAL LETTER DELTA
   deltaString = Wgl4Ascii ("\xce\x94", "error");

   // GREEK SMALL LETTER EPSILON
   epsilonString = Wgl4Ascii ("\xce\xb5", "eps");

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
      
      if (verbose >= -1)
      {
         cout << "Testing Rule #" << rule << " ";
         doubleQuote (cout, name);
         cout << (verbose == -1 ? "..." : ":") << endl;
      }

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
            NORMAL  cout << "   Cannot create rule in dimension "
                         << e.getDimension() << "!" << endl;
         }
      }

           if (verbose == 0)  cout << endl;
      else if (verbose >= 1)  cout << '\n';
   }
}


