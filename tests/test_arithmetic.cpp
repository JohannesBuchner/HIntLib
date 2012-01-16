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
#include <algorithm>

#include <HIntLib/polynomial2.h>
#include <HIntLib/integerring.h>
#include <HIntLib/galoisfield.h>

#include "test.h"

using std::cerr;
using std::cout;
using std::endl;
using std::setw;

using namespace HIntLib;

unsigned SIZE = 100;
unsigned W = 15;


/**
 *   test Ring ()
 */

template<class R>
void testRing (const R &r)
{
   typedef typename R::type T;

   const T zero = r.zero();
   if (r.size() && r.index(zero) >= r.size())  error ("index(zero) invalid");

   T one = r.one();
   if (r.size() && r.index(one) >= r.size())  error ("index(one) invalid");

   unsigned S = SIZE;
   if (r.size() > 0)  S = std::min (S, r.size());

   for (unsigned i = 0; i < S; i++)
   {
      T a = r.element (i);

      DEB1 cout << "  Checking: " << setw(W) << a;

      // element() and index()

      DEB2 cout << ", index = " << i;
      if (r.index(a) != i)  error ("index(a) != i");

      // is0() and is1()

      bool is0 = r.is0(a);
      bool is1 = r.is1(a);
      DEB2
      {
         if (is0)  cout << ", is zero";
         if (is1)  cout << ", is one";
      }
      if (is0 != (a == zero))  error ("is0() fails");
      if (is1 != (a == one ))  error ("is1() fails");
      
      // zero

      if (r.add (a, zero) != a)  error ("a+0 != a");

      // one

      if (r.mul (a, one) != a)  error ("a*1 != a");

      // additive inverse

      T neg = r.neg(a);
      DEB2 cout << ", -a = " << neg;

      if (r.size() && r.index(neg) >= r.size())  error ("-a invalid");
      if (! r.is0 (r.add (a, neg)))  error ("a + (-a) != 0");

      {
         // negate()

         T x = a;
         r.negate (x);
         if (x != neg)  error ("negate(a) != -a");
      }

      DEB1 cout << endl;
   }

   S = unsigned(sqrt(double(SIZE)));
   if (r.size() > 0)  S = std::min (S, r.size());

   for (unsigned i = 0; i < S; i++)
   {
      const T a = r.element (i);

      DEB3 cout << "  Checking: " << setw(W) << a;

      {  // check times()

         DEB3 cout << ", ";

         T x = zero;
         for (unsigned l = 0; l < SIZE / S; ++l)
         {
            DEB3 cout << r.times (a, l) << '/';
            if (x != r.times (a, l))  error ("a+...+a != times (a,k)");
            r.addTo (x, a);
         }
      }

      {  // check power()

         DEB3 cout << ", ";
         T x = a;
         for (unsigned l = 1; l < SIZE / S; ++l)
         {
            DEB3 cout << r.power (a, l) << '/';
            if (x != r.power (a, l))  error ("a*...*a != power (a,k)");
            try
            {
               r.mulBy (x, a);
            }
            catch (Overflow &)
            {
               break;
            }
         }
      }

      DEB3 cout << endl;

      for (unsigned j = 0; j < S; j++)
      {
         T b = r.element (j);

         DEB1 cout << "  Checking: " << setw(W) << a << setw(W) << b;

         // Make sure, == make sense

         if ((a == b) != (i == j))  error ("a==b incorrect");
         if ((a != b) != (i != j))  error ("a!=b incorrect");
         
         // Make sure, operators dont leave domain

         T add = r.add (a,b);
         T sub = r.sub (a,b);
         T mul = r.mul (a,b);

         if (r.size())
         {
            if (r.index (add) >= r.size())  error ("a+b invalid");
            if (r.index (sub) >= r.size())  error ("a-b invalid");
            if (r.index (mul) >= r.size())  error ("a*b invalid");
         }
         
         // Commutativity

         DEB2 cout << ", a + b = " << add << ", a * b = " << mul;
         if (add != r.add (b,a))  error ("a+b != b+a");
         if (mul != r.mul (b,a))  error ("a*b != b*a");

         // other functions

         T negB = r.neg (b);
         if (r.add (a, negB) != sub)  error ("sub() fails");

         T x;

         x = a;
         r.addTo (x, b);
         if (x != add) error ("addTo() fails");

         x = a;
         r.subFrom (x, b);
         if (x != sub) error ("subFrom() fails");

         x = a;
         r.mulBy (x, b);
         if (x != mul) error ("mulBy() fails");

         DEB1 cout << endl;
      }
   }

   S = unsigned(pow(double(SIZE), 1.0 / 3.0));
   if (r.size() > 0)  S = std::min (S, r.size());

   for (unsigned i = 0; i < S; i++)
   {
      for (unsigned j = 0; j < S; j++)
      {
         for (unsigned k = 0; k < S; ++k)
         {
            T a = r.element (i);
            T b = r.element (j);
            T c = r.element (k);

            DEB1 cout << "  Checking: " << a << ", " << b << ", " << c << endl;

            // Associativity

            if (r.add(r.add(a,b), c) != r.add(a, r.add(b,c)))
                  error("(a+b) + c != a + (b+c)");
            if (r.mul(r.mul(a,b), c) != r.mul(a, r.mul(b,c)))
                  error("(a*b) * c != a * (b*c)");

            // Distributivity

            if (r.mul (a, r.add(b,c)) != r.add (r.mul(a,b), r.mul(a,c)))
                  error("a*(b+c) != (a*b)+(a*c)");
            if (r.mul (r.add(a,b), c) != r.add (r.mul(a,c), r.mul(b,c)))
               error("(a+b)*c != (a*c)+(b*c)");
         }
      }
   }

}


/**
 *  test UFD ()
 */

template<class R>
void testUFD (const R &r)
{
   typedef typename R::type T;

   testRing (r);

   unsigned S = SIZE;
   if (r.size() > 0)  S = std::min (S, r.size());

   for (unsigned i = 0; i < S; i++)
   {
      T a = r.element (i);

      DEB1 cout << "  Checking: " << setw(W) << a;

      bool is0     = r.is0 (a);
      bool isUnit  = r.isUnit (a);
      bool isPrime = r.isPrime (a);

      DEB2
      {
         if (is0)      cout << ", is zero";
         if (isUnit)   cout << ", is unit";
         if (isPrime)  cout << ", is prime";
      }

      if (isUnit)
      {
         T recip = r.unitRecip (a);
         DEB2  cout << ", a^-1 = " << recip;
         if (! r.is1 (r.mul (a, recip)))
         {
            error ("Reciprocal of unit is wrong!");
         }
      }

      if (is0 + isUnit + isPrime > 1)
      {
         error ("Classifiaction ambiguous!");
      }
      DEB1 cout << endl;
   }

   S = unsigned(sqrt(double(SIZE)));
   if (r.size() > 0)  S = std::min (S, r.size());

   for (unsigned i = 0; i < S; i++)
   {
      for (unsigned j = 0; j < S; j++)
      {
         T a = r.element (i);
         T b = r.element (j);
         DEB1 cout << "  Checking: " << setw(W) << a << setw(W) << b;

         if (! r.is0(a) && ! r.is0(b))
         {
            T prod = r.mul (a,b);
            DEB2 cout << ", a*b = " << prod;

            bool aIsUnit = r.isUnit (a);
            bool bIsUnit = r.isUnit (b);

            if (aIsUnit && bIsUnit)
            {
               if (! r.isUnit (prod))
                  error ("Product of units must be a unit!");
            }
            else
            {
               bool aIsPrime = r.isPrime (a);
               bool bIsPrime = r.isPrime (b);

               if ((aIsUnit && bIsPrime) || (aIsPrime && bIsUnit))
               {
                  if (! r.isPrime (prod))
                     error ("Product of unit and prime must be prime!");
               }
               else
               {
                  if (r.isUnit (prod) || r.isPrime (prod) || r.is0 (prod))
                  {
                     error ("product of two non-units must be composit!");
                  }
               }
            }
         }

         if (! r.is0 (b))
         {
            T q, re;
            r.div (a, b, q, re);
            DEB2 cout << ", quotient = " << q << ", remainder = " << re;
            if (r.add (r.mul (b, q), re) != a)  error ("div() fails");
         }

         DEB1 cout << endl;
      }
   }
}


/**
 *  test Field ()
 */

template<class F>
void testField (const F &f)
{
   typedef typename F::type T;

   testUFD (f);

   unsigned S = unsigned(sqrt(double(SIZE)));
   if (f.size() > 0)  S = std::min (S, f.size());

   for (unsigned i = 0; i < S; i++)
   {
      T a = f.element (i);

      DEB1 cout << "  Checking: " << setw(W) << a;

      if (! f.is0 (a))
      {
         // multiplicative inverse

         if (! f.isUnit (a))  error ("Each elemnt != 0 must be a unit!");
         T r = f.recip(a);
         DEB2 cout << ", a^-1 = " << r;
         if (f.size() && f.index(r) >= f.size())  error ("1/b invalid");
         if (! f.is1 (f.mul (a, r)))  error ("a * (a^-1) != 1");
      }

      DEB1 cout << endl;
   }

   S = unsigned(sqrt(double(SIZE)));
   if (f.size() > 0)  S = std::min (S, f.size());

   for (unsigned i = 0; i < S; i++)
   {
      for (unsigned j = 0; j < S; j++)
      {
         T a = f.element (i);
         T b = f.element (j);

         DEB1 cout << "  Checking: " << a << ", " << b << endl;

         if (! f.is0(b))
         {
            // Make sure, operators dont leave domain

            if (f.size())
            {
               if (f.index (f.div (a,b)) >= f.size())  error ("a/b invalid");
            }

            // other functions

            T r = f.recip (b);
            if (f.mul (a, r) != f.div (a, b))  error ("div() fails");

            T x;

            x = a;
            f.divBy (x, b);
            if (x != f.div (a,b)) error ("divBy() fails");
         }
      }
   }
}


const char* options = "n:";

bool opt (int c, const char* s)
{
   switch (c)
   {
   case 'n':  SIZE = atoi (s); return true;
   }

   return false;
}

void usage()
{
   cerr <<
      "Usage: test_arithmetic [OPTION]...\n\n"
      "Tests all kinds of rings and fields.\n\n"
      << option_msg <<
      "  -n n   Uses tables of size _n_ for the tests (default = 1000).\n"
      "\n";

   exit (1);
}

void test (int argc, char**)
{
   if (argc)  usage();

   // Integers

   IntegerRing<int> intRing;
   NORMAL cout << "Testing Euclidian Ring " << intRing << "..." << endl;
   testUFD (intRing);

   // Polynomials over the integers

   PolynomialRing<IntegerRing<int> > polyInt (intRing);
   NORMAL cout << "Testing           Ring " << polyInt << "..." << endl;
   testRing (polyInt);

   // Polynomials over GF(2)

   Polynomial2Ring<u32> poly2;
   NORMAL cout << "Testing Euclidian Ring " << poly2 << "..." << endl;
   testUFD (poly2);

   // Modular arithmetic

   for (unsigned i = 2; i <= 30; ++i)
   {
      if (Prime::test (i))
      {
         ModularIntegerField<unsigned short> field (i);
         NORMAL cout << "Testing          Field " << field << "..." << endl;
         testField (field);

         PolynomialRingField<ModularIntegerField<unsigned short> > poly (field);
         NORMAL cout << "Testing Euclidian Ring " << poly << "..." << endl;
         testUFD (poly);
         
         unsigned short coef [] = {1, 0, 0};
         PolynomialRingField<ModularIntegerField<unsigned short> >::type
            p (sizeof(coef)/sizeof(unsigned short), coef);
         ModularArithmeticRing<PolynomialRingField<
            ModularIntegerField<unsigned short> > > mod (poly, p);
         NORMAL cout << "Testing           Ring " << mod << "..." << endl;
         testRing (mod);

         for (unsigned j = 1; j <= 4; ++j)
         {
            GaloisField<unsigned short> gf (i, j);
            NORMAL cout << "Testing          Field " << gf << "..." << endl;
            testField (gf);
         }
      }
      else
      {
         ModularIntegerRing<unsigned short> modInt (i);
         NORMAL cout << "Testing           Ring " << modInt << "..." << endl;
         testRing (modInt);

         PolynomialRing<ModularIntegerRing<unsigned short> > polyModInt(modInt);
         NORMAL cout << "Testing           Ring " << polyModInt << "..." <<endl;
         testRing (polyModInt);
      }
   }

#if 1
   ModularIntegerRing<unsigned short> modInt (213);
   NORMAL cout << "Testing           Ring " << modInt << "..." << endl;
   testRing (modInt);

   ModularIntegerField<unsigned short> gf1999 (1999);
   NORMAL cout << "Testing          Field " << gf1999 << "..." << endl;
   testField (gf1999);
#endif

}

