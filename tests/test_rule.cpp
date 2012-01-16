/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002,03,04,05,06,07  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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
int NUM_TESTS      =   3;         // number of test-polynomials / dimesnion
int FIRST_TEST_DIM =   1;
int LAST_TEST_DIM  = 100;
int FIRST_RULE     =   1;
int LAST_RULE      = 200;

const int  NUM_MON = 130;        // Number of monomials in each polynomial
const real MAX_COEFF = 10.0;     // Maximal coefficient in each polynomial
bool STANDARD_UNIT_CUBE = false;


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

   TestPolynomial (int dim, int degree);

   real operator() (const real *);
   real derivative (const real *, int);

   real getExactResult (const Hypercube &) const;
   real getPartialResult (const Hypercube &, int) const;
   real getErrorOfExactResult (const Hypercube&) const;
   real getErrorOfEvaluation (const Hypercube&) const;

private:

   // Static member with constructor used to initialize static data

   static
   struct Dummy
   {
      Dummy ();
   } dummy;

   const int degree;

   real coef [NUM_MON];

   Array<unsigned char> exponents;  // [NUM_MON][dim]

   unsigned char& getExponent (int i, int d)
      { return exponents [d + i * dim]; }
   unsigned char getExponent (int i, int d) const
      { return exponents [d + i * dim]; }

   /*
    * Static monomial table
    */

   enum { MAX_DEGREE = 10 }; // Highest degree in use
   static const int NUM_MONOMIALS;   // Size of the following three tables

   static const unsigned char monomials [][MAX_DEGREE]; // monomials
   static int monomialsDim [];                     // # differnt variables
   static int monomialsDeg [];                     // degree

   // Formatted output

   friend ostream& operator<< (ostream&, const TestPolynomial&);
};


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
   { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},

   // the only rules requiring degrees >= 11 are only two-dimensional
   // Therefore only the monomials (d-i,i,0,...,0) are required

   {11, 0, 0, 0, 0, 0, 0, 0, 0, 0},  // deg 11  (incomplete, only dim <= 3)
   {10, 1, 0, 0, 0, 0, 0, 0, 0, 0},
   { 9, 2, 0, 0, 0, 0, 0, 0, 0, 0},
   { 9, 1, 1, 0, 0, 0, 0, 0, 0, 0},
   { 8, 3, 0, 0, 0, 0, 0, 0, 0, 0},
   { 8, 2, 1, 0, 0, 0, 0, 0, 0, 0},
   { 7, 4, 0, 0, 0, 0, 0, 0, 0, 0},
   { 7, 3, 1, 0, 0, 0, 0, 0, 0, 0},
   { 7, 2, 2, 0, 0, 0, 0, 0, 0, 0},
   { 6, 5, 0, 0, 0, 0, 0, 0, 0, 0},
   { 6, 4, 1, 0, 0, 0, 0, 0, 0, 0},
   { 6, 3, 2, 0, 0, 0, 0, 0, 0, 0},
   {12, 0, 0, 0, 0, 0, 0, 0, 0, 0},  // deg 12  (incomplete, only dim <= 3)
   {11, 1, 0, 0, 0, 0, 0, 0, 0, 0},
   {10, 2, 0, 0, 0, 0, 0, 0, 0, 0},
   {10, 1, 1, 0, 0, 0, 0, 0, 0, 0},
   { 9, 3, 0, 0, 0, 0, 0, 0, 0, 0},
   { 9, 2, 1, 0, 0, 0, 0, 0, 0, 0},
   { 8, 4, 0, 0, 0, 0, 0, 0, 0, 0},
   { 8, 3, 1, 0, 0, 0, 0, 0, 0, 0},
   { 8, 2, 2, 0, 0, 0, 0, 0, 0, 0},
   { 7, 5, 0, 0, 0, 0, 0, 0, 0, 0},
   { 7, 4, 1, 0, 0, 0, 0, 0, 0, 0},
   { 7, 3, 2, 0, 0, 0, 0, 0, 0, 0},
   { 6, 6, 0, 0, 0, 0, 0, 0, 0, 0},
   { 6, 5, 1, 0, 0, 0, 0, 0, 0, 0},
   { 6, 4, 2, 0, 0, 0, 0, 0, 0, 0},
   { 6, 3, 3, 0, 0, 0, 0, 0, 0, 0},
   {13, 0, 0, 0, 0, 0, 0, 0, 0, 0},  // deg 13  (incomplete, only dim <= 3)
   {12, 1, 0, 0, 0, 0, 0, 0, 0, 0},
   {11, 2, 0, 0, 0, 0, 0, 0, 0, 0},
   {11, 1, 1, 0, 0, 0, 0, 0, 0, 0},
   {10, 3, 0, 0, 0, 0, 0, 0, 0, 0},
   {10, 2, 1, 0, 0, 0, 0, 0, 0, 0},
   { 9, 4, 0, 0, 0, 0, 0, 0, 0, 0},
   { 9, 3, 1, 0, 0, 0, 0, 0, 0, 0},
   { 9, 2, 2, 0, 0, 0, 0, 0, 0, 0},
   { 8, 5, 0, 0, 0, 0, 0, 0, 0, 0},
   { 8, 4, 1, 0, 0, 0, 0, 0, 0, 0},
   { 8, 3, 2, 0, 0, 0, 0, 0, 0, 0},
   { 7, 6, 0, 0, 0, 0, 0, 0, 0, 0},
   { 7, 5, 1, 0, 0, 0, 0, 0, 0, 0},
   { 7, 4, 2, 0, 0, 0, 0, 0, 0, 0},
   { 7, 3, 3, 0, 0, 0, 0, 0, 0, 0},
   {14, 0, 0, 0, 0, 0, 0, 0, 0, 0},  // deg 14  (incomplete, only dim <= 3)
   {13, 1, 0, 0, 0, 0, 0, 0, 0, 0},
   {12, 2, 0, 0, 0, 0, 0, 0, 0, 0},
   {12, 1, 1, 0, 0, 0, 0, 0, 0, 0},
   {11, 3, 0, 0, 0, 0, 0, 0, 0, 0},
   {11, 2, 1, 0, 0, 0, 0, 0, 0, 0},
   {10, 4, 0, 0, 0, 0, 0, 0, 0, 0},
   {10, 3, 1, 0, 0, 0, 0, 0, 0, 0},
   {10, 2, 2, 0, 0, 0, 0, 0, 0, 0},
   { 9, 5, 0, 0, 0, 0, 0, 0, 0, 0},
   { 9, 4, 1, 0, 0, 0, 0, 0, 0, 0},
   { 9, 3, 2, 0, 0, 0, 0, 0, 0, 0},
   { 8, 6, 0, 0, 0, 0, 0, 0, 0, 0},
   { 8, 5, 1, 0, 0, 0, 0, 0, 0, 0},
   { 8, 4, 2, 0, 0, 0, 0, 0, 0, 0},
   { 8, 3, 3, 0, 0, 0, 0, 0, 0, 0},
   { 7, 7, 0, 0, 0, 0, 0, 0, 0, 0},
   { 7, 6, 1, 0, 0, 0, 0, 0, 0, 0},
   { 7, 5, 2, 0, 0, 0, 0, 0, 0, 0},
   { 7, 4, 3, 0, 0, 0, 0, 0, 0, 0},
   {15, 0, 0, 0, 0, 0, 0, 0, 0, 0},  // deg 15  (incomplete, only dim <= 3)
   {14, 1, 0, 0, 0, 0, 0, 0, 0, 0},
   {13, 2, 0, 0, 0, 0, 0, 0, 0, 0},
   {13, 1, 1, 0, 0, 0, 0, 0, 0, 0},
   {12, 3, 0, 0, 0, 0, 0, 0, 0, 0},
   {12, 2, 1, 0, 0, 0, 0, 0, 0, 0},
   {11, 4, 0, 0, 0, 0, 0, 0, 0, 0},
   {11, 3, 1, 0, 0, 0, 0, 0, 0, 0},
   {11, 2, 2, 0, 0, 0, 0, 0, 0, 0},
   {10, 5, 0, 0, 0, 0, 0, 0, 0, 0},
   {10, 4, 1, 0, 0, 0, 0, 0, 0, 0},
   {10, 3, 2, 0, 0, 0, 0, 0, 0, 0},
   { 9, 6, 0, 0, 0, 0, 0, 0, 0, 0},
   { 9, 5, 1, 0, 0, 0, 0, 0, 0, 0},
   { 9, 4, 2, 0, 0, 0, 0, 0, 0, 0},
   { 9, 3, 3, 0, 0, 0, 0, 0, 0, 0},
   { 8, 7, 0, 0, 0, 0, 0, 0, 0, 0},
   { 8, 6, 1, 0, 0, 0, 0, 0, 0, 0},
   { 8, 5, 2, 0, 0, 0, 0, 0, 0, 0},
   { 8, 4, 3, 0, 0, 0, 0, 0, 0, 0},
   {16, 0, 0, 0, 0, 0, 0, 0, 0, 0},  // deg 16  (incomplete, only dim <= 3)
   {15, 1, 0, 0, 0, 0, 0, 0, 0, 0},
   {14, 2, 0, 0, 0, 0, 0, 0, 0, 0},
   {14, 1, 1, 0, 0, 0, 0, 0, 0, 0},
   {13, 3, 0, 0, 0, 0, 0, 0, 0, 0},
   {13, 2, 1, 0, 0, 0, 0, 0, 0, 0},
   {12, 4, 0, 0, 0, 0, 0, 0, 0, 0},
   {12, 3, 1, 0, 0, 0, 0, 0, 0, 0},
   {12, 2, 2, 0, 0, 0, 0, 0, 0, 0},
   {11, 5, 0, 0, 0, 0, 0, 0, 0, 0},
   {11, 4, 1, 0, 0, 0, 0, 0, 0, 0},
   {11, 3, 2, 0, 0, 0, 0, 0, 0, 0},
   {10, 6, 0, 0, 0, 0, 0, 0, 0, 0},
   {10, 5, 1, 0, 0, 0, 0, 0, 0, 0},
   {10, 4, 2, 0, 0, 0, 0, 0, 0, 0},
   {10, 3, 3, 0, 0, 0, 0, 0, 0, 0},
   { 9, 7, 0, 0, 0, 0, 0, 0, 0, 0},
   { 9, 6, 1, 0, 0, 0, 0, 0, 0, 0},
   { 9, 5, 2, 0, 0, 0, 0, 0, 0, 0},
   { 9, 4, 3, 0, 0, 0, 0, 0, 0, 0},
   { 8, 8, 0, 0, 0, 0, 0, 0, 0, 0},
   { 8, 7, 1, 0, 0, 0, 0, 0, 0, 0},
   { 8, 6, 2, 0, 0, 0, 0, 0, 0, 0},
   { 8, 5, 3, 0, 0, 0, 0, 0, 0, 0},
   { 8, 4, 4, 0, 0, 0, 0, 0, 0, 0},
   {17, 0, 0, 0, 0, 0, 0, 0, 0, 0},  // deg 17  (incomplete, only dim <= 3)
   {16, 1, 0, 0, 0, 0, 0, 0, 0, 0},
   {15, 2, 0, 0, 0, 0, 0, 0, 0, 0},
   {15, 1, 1, 0, 0, 0, 0, 0, 0, 0},
   {14, 3, 0, 0, 0, 0, 0, 0, 0, 0},
   {14, 2, 1, 0, 0, 0, 0, 0, 0, 0},
   {13, 4, 0, 0, 0, 0, 0, 0, 0, 0},
   {13, 3, 1, 0, 0, 0, 0, 0, 0, 0},
   {13, 2, 2, 0, 0, 0, 0, 0, 0, 0},
   {12, 5, 0, 0, 0, 0, 0, 0, 0, 0},
   {12, 4, 1, 0, 0, 0, 0, 0, 0, 0},
   {12, 3, 2, 0, 0, 0, 0, 0, 0, 0},
   {11, 6, 0, 0, 0, 0, 0, 0, 0, 0},
   {11, 5, 1, 0, 0, 0, 0, 0, 0, 0},
   {11, 4, 2, 0, 0, 0, 0, 0, 0, 0},
   {11, 3, 3, 0, 0, 0, 0, 0, 0, 0},
   {10, 7, 0, 0, 0, 0, 0, 0, 0, 0},
   {10, 6, 1, 0, 0, 0, 0, 0, 0, 0},
   {10, 5, 2, 0, 0, 0, 0, 0, 0, 0},
   {10, 4, 3, 0, 0, 0, 0, 0, 0, 0},
   { 9, 8, 0, 0, 0, 0, 0, 0, 0, 0},
   { 9, 7, 1, 0, 0, 0, 0, 0, 0, 0},
   { 9, 6, 2, 0, 0, 0, 0, 0, 0, 0},
   { 9, 5, 3, 0, 0, 0, 0, 0, 0, 0},
   { 9, 4, 4, 0, 0, 0, 0, 0, 0, 0},
   {18, 0, 0, 0, 0, 0, 0, 0, 0, 0},  // deg 18  (incomplete, only dim <= 3)
   {17, 1, 0, 0, 0, 0, 0, 0, 0, 0},
   {16, 2, 0, 0, 0, 0, 0, 0, 0, 0},
   {16, 1, 1, 0, 0, 0, 0, 0, 0, 0},
   {15, 3, 0, 0, 0, 0, 0, 0, 0, 0},
   {15, 2, 1, 0, 0, 0, 0, 0, 0, 0},
   {14, 4, 0, 0, 0, 0, 0, 0, 0, 0},
   {14, 3, 1, 0, 0, 0, 0, 0, 0, 0},
   {14, 2, 2, 0, 0, 0, 0, 0, 0, 0},
   {13, 5, 0, 0, 0, 0, 0, 0, 0, 0},
   {13, 4, 1, 0, 0, 0, 0, 0, 0, 0},
   {13, 3, 2, 0, 0, 0, 0, 0, 0, 0},
   {12, 6, 0, 0, 0, 0, 0, 0, 0, 0},
   {12, 5, 1, 0, 0, 0, 0, 0, 0, 0},
   {12, 4, 2, 0, 0, 0, 0, 0, 0, 0},
   {12, 3, 3, 0, 0, 0, 0, 0, 0, 0},
   {11, 7, 0, 0, 0, 0, 0, 0, 0, 0},
   {11, 6, 1, 0, 0, 0, 0, 0, 0, 0},
   {11, 5, 2, 0, 0, 0, 0, 0, 0, 0},
   {11, 4, 3, 0, 0, 0, 0, 0, 0, 0},
   {10, 8, 0, 0, 0, 0, 0, 0, 0, 0},
   {10, 7, 1, 0, 0, 0, 0, 0, 0, 0},
   {10, 6, 2, 0, 0, 0, 0, 0, 0, 0},
   {10, 5, 3, 0, 0, 0, 0, 0, 0, 0},
   {10, 4, 4, 0, 0, 0, 0, 0, 0, 0},
   { 9, 9, 0, 0, 0, 0, 0, 0, 0, 0},
   { 9, 8, 1, 0, 0, 0, 0, 0, 0, 0},
   { 9, 7, 2, 0, 0, 0, 0, 0, 0, 0},
   { 9, 6, 3, 0, 0, 0, 0, 0, 0, 0},
   { 9, 5, 4, 0, 0, 0, 0, 0, 0, 0},
};

const int TestPolynomial::NUM_MONOMIALS
     = sizeof (monomials) / sizeof (monomials [0]);

int TestPolynomial::monomialsDim [NUM_MONOMIALS];
int TestPolynomial::monomialsDeg [NUM_MONOMIALS];
TestPolynomial::Dummy TestPolynomial::dummy;


/**
 *  Initialize static members monomialsDeg[] and monomialsDim[]
 */

TestPolynomial::Dummy::Dummy ()
{
   for (int i = 0; i < NUM_MONOMIALS; ++i)
   {
      int sum = 0;
      int dim = 0;

      for (int j = 0; j < MAX_DEGREE; ++j)
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

TestPolynomial::TestPolynomial (int dim, int degree)
   : TestIntegrand (dim), degree (degree), exponents (NUM_MON * dim)
{
   int m = NUM_MONOMIALS;

   for (int i = 0; i < NUM_MON; ++i)
   {
      // Find next monomial with acceptable degree and dimension

      do
      {
         if (m == 0)  m = NUM_MONOMIALS;
         --m;
      }
      while (monomialsDeg [m] > degree || monomialsDim [m] > dim);

      // Determine coefficient

      coef [i] = uniform (mt, -MAX_COEFF, MAX_COEFF);

      // Copy monomomial

      for (int j = 0; j < min (dim, int(MAX_DEGREE)); ++j)
      {
         getExponent(i, j) = monomials [m][j];
      }

      // Fill remaining exponents with 0

      for (int j = MAX_DEGREE; j < dim; ++j)  getExponent(i, j) = 0;

      // Shuffle Exponents

      random_shuffle (&getExponent (i,0), &getExponent (i,dim), mt);
   }
}


/**
 *  getExactResult()
 */

real TestPolynomial::getExactResult (const Hypercube &h) const
{
   real sum = 0;

   for (int i = 0; i < NUM_MON; ++i)
   {
      real prod = coef[i];

      for (int d = 0; d < dim; ++d)
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

real TestPolynomial::getPartialResult (const Hypercube &h, int deg) const
{
   real sum = 0;

   for (int i = 0; i < NUM_MON; ++i)
   {
      real prod = coef [i];

      int thisDegree = 0;
      for (int d = 0; d < dim; ++d)  thisDegree += getExponent (i, d);
      if (thisDegree != deg)  continue;

      for (int d = 0; d < dim; ++d)
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

   for (int i = 0; i < NUM_MON; ++i)
   {
      real realProd = 1;
      real prod     = onePlusEpsilon;
      real factor = coef [i];

      for (int d = 0; d < dim; ++d)
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
            
            factor /= real(e);
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

   for (int i = 0; i < NUM_MON; ++i)   // all monomials in polynomial
   {
      real prod = coef [i];

      for (int d = 0; d < dim; ++d)
      {
         const int exponent = getExponent (i,d);
         if (exponent)  prod *= HINTLIB_MN pow (x[d], exponent);
      }

      sum += prod;
   }

   DEB4
   {
      cout << "      f(";
      for (int d = 0; d < dim; ++d)
      {
         cout << x[d];
         if (d + 1 < dim)  cout << ',';
      }
      cout << ") = " << sum << '\n';
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

   for (int i = 0; i < h.getDimension(); ++i)
   {
      maxCoordinate [i]
         = max (L::abs (h.getLowerBound(i)), L::abs (h.getUpperBound(i)));
   }

   real sum = 0.0;

   for (int i = 0; i < NUM_MON; ++i)   // all monomials in polynomial
   {
      real prod = L::abs (coef [i]);

      int degreeOfMonomial = 0;

      for (int d = 0; d < dim; ++d)
      {
         degreeOfMonomial += getExponent (i, d);
         prod *= HINTLIB_MN pow (maxCoordinate [d], getExponent (i, d));
      }

      sum += prod * (epsilonPowerMinusOne (degreeOfMonomial + 1));
   }

   return sum;
}

real TestPolynomial::derivative (const real *x, int a)
{
   real sum = 0.0;

   for (int i = 0; i < NUM_MON; ++i)   // all monomials in polynomial
   {
      if (getExponent (i, a) == 0)  continue;

      real prod = coef [i];

      for (int d = 0; d < dim; ++d)
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
real TestPolynomial::derivative (const real *x, int a, int b)
{
   real sum = 0.0;

   if (a == b)
   {
      for (int i = 0; i < NUM_MON; ++i)   // all monomials in polynomial
      {
         if (getExponent (i, a) <= 1)  continue;

         real prod = coef [i];

         for (int d = 0; d < dim; ++d)
         {
            if (d == a)
            {
               int exponent = getExponent (i, a);
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
      for (int i = 0; i < NUM_MON; ++i)   // all monomials in polynomial
      {
         if (getExponent (i, a) == 0)  continue;
         if (getExponent (i, b) == 0)  continue;

         real prod = coef [i];

         for (int d = 0; d < dim; ++d)
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

   for (int i = 0; i < NUM_MON; ++i)
   {
      ss << ' ';
      if (p.coef [i] > 0.0)  ss << '+';
      else ss.minusSign();
      ss << ' ' << L::abs (p.coef [i]);
      bool spacePrinted = false;

      for (int d = 0; d < p.dim; ++d)
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

bool checkDegree (CubatureRule& r, int degree, int numTests, bool disproof)
{
   const int dim = r.getDimension ();

   const Hypercube h (dim, STANDARD_UNIT_CUBE ? -1 : -.2,
                           STANDARD_UNIT_CUBE ?  1 : 1.7);

   for (int i = 0; i < numTests; ++i)
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
         cout << setprecision(10)
              << "     Q(f)=" << setw(13) << result
              << "I(f)="   << setw(13) << correct
              << setprecision(4)
              << deltaString << '=' << setw(10) << error
              << epsilonString << '=' << setw(10) << threshold
              << epsilonString << "_exact=" << setw(10) << expectedExactError
              << epsilonString << "_approx=" << setw(10) << expectedApproxError
              << "probl.=" << problem << endl;
      }
      
      if (error > 25. * threshold)
      {
         DEB1 cout << "     Error is " << error << ", exceeding " << threshold
                   << " by a factor of " << error / threshold << ".\n";
         return false;
      }

      if (disproof && threshold > problem / 1400.)
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

void testRule (CubatureRule& r, int dim)
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

   const int degree = r.getDegree ();

#if 0
   DEB1
   {
      cout << "   Degree " << degree << ", " << num_points
           << " points, sum |w| = " << r.getSumAbsWeight() << ", "
           << (r.isAllPointsInside() ? "" : "not ") << "all points inside.\n";
   }
#endif

   bool allright = true;

   for (int d = 0; d <= degree; ++d)
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
 *  -vvvv  5     + each integrand evaluation
 */

const char options[] = "n:p:d:r:u";
const char option_msg[] =
   "  -n n   Number of testruns to be performed (default = 3)\n"
   "  -p n   Tests requiring more that _p_ points are skipped"
         " (default = 3000)\n"
   "  -d [nmin-]nmax   Range of dimensions (default nmin = 1, nmax = 100)\n"
   "  -r n   Tests only rule number _n_\n"
   "  -u     Use standard unit cube [-1,1]^s\n";
const char testProgramParameters[] = "[OPTION]...";
const char testProgramUsage[] =
   "Tests whether cubature rules work correctly.\n\n";
const char testProgramName[] = "test_rule";
const int  testProgramCopyright = 2003;

bool opt(int c, const char* s)
{
   switch (c)
   {
      case 'n':  NUM_TESTS      = parseInt (s); return true;
      case 'p':  MAX_NUM_POINTS = parseInt (s); return true;
      case 'd':  parseRange (s, 1, FIRST_TEST_DIM, LAST_TEST_DIM); return true;
      case 'r':  LAST_RULE = FIRST_RULE = parseInt (s); return true;
      case 'u':  STANDARD_UNIT_CUBE = true; return true;
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

   for (int rule = FIRST_RULE; rule <= LAST_RULE; ++rule)
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

      for (int dim = FIRST_TEST_DIM; dim <= LAST_TEST_DIM; ++dim)
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


