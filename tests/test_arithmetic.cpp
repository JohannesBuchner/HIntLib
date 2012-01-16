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
#include <HIntLib/lookupfield.h>
#include <HIntLib/onedimvectorspace.h>
#include <HIntLib/array.h>

#include "test.h"

using std::cerr;
using std::cout;
using std::endl;
using std::setw;

using namespace HIntLib;

unsigned SIZE = 100;
unsigned W = 15;


/*************  Test procedure for rings  ************************************/


/**
 *   Ring
 */

template<class R>
void doTest (const R &r, ring_tag)
{
   // doTest (r, XXXX);

   const typename R::type one = r.one();
   if (r.size() && r.index(one) >= r.size())  error ("index(one) invalid");
   DEB2 { cout << "  One:  "; r.printShort (cout, one); cout << endl; }
}

template<class R>
void doTest (const R &r, const typename R::type& a, ring_tag)
{
   // doTest (r, a, XXXX);

   // is1(), one()

   const typename R::type one = r.one();

   bool is1 = r.is1(a);
   DEB2 if (is1)  cout << ", is one";
   if (is1 != (a == one ))  error ("is1() fails");
      
   // neutral

   if (r.mul (a, one) != a)  error ("a*1 != a");

   // check power()

   if (sqr(r.index(a)) < SIZE)
   {
      DEB3 cout << ", power = ";

      unsigned SS = unsigned (sqrt(double (SIZE)));

      typename R::type x = a;
      for (unsigned l = 1; l < SS; ++l)
      {
         DEB3 { r.printShort (cout, r.power (a, l)); cout << '/'; }
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

}

template<class R>
void doTest
   (const R& r, const typename R::type& a, const typename R::type& b, ring_tag)
{
   // doTest (r, a, b, XXXX);

   // mul()

   const typename R::type mul = r.mul (a,b);

   DEB2
   {
      cout << ", a*b = ";
      r.printShort (cout, mul);
   }

   if (r.size())  if (r.index (mul) >= r.size())  error ("a*b invalid");

   // Commutativity

   if (mul != r.mul (b,a))  error ("a*b != b*a");

   typename R::type x = a;
   r.mulBy (x, b);
   if (x != mul) error ("mulBy() fails");
}

template<class R>
void doTest
   (const R& r, const typename R::type& a,
                const typename R::type& b,
                const typename R::type& c, ring_tag)
{
   // doTest (r, a, b, c, XXXX);

   // Associativity of multiplication

   if (r.mul(r.mul(a,b), c) != r.mul(a, r.mul(b,c)))
         error("(a*b) * c != a * (b*c)");

   // Distributivity

   if (r.mul (a, r.add(b,c)) != r.add (r.mul(a,b), r.mul(a,c)))
         error("a*(b+c) != (a*b)+(a*c)");
   if (r.mul (r.add(a,b), c) != r.add (r.mul(a,c), r.mul(b,c)))
      error("(a+b)*c != (a*c)+(b*c)");
}


/**
 *  Domain
 */

template<class R>
void doTest
  (const R &r, const typename R::type& a, const typename R::type& b, domain_tag)
{
   doTest (r, a, b, ring_tag());

   if (r.is0(a) || r.is0(b))  return;
   if (r.is0(r.mul (a, b)))  error ("a * b = 0!");
}


/**
 *  UFD
 */

template<class R>
void doTest (const R &r, const typename R::type& a, ufd_tag)
{
   typedef typename R::type T;

   doTest (r, a, domain_tag());

   bool is0     = r.is0 (a);
   bool isUnit  = r.isUnit (a);
   bool isPrime = r.isPrime (a);
   bool isIrred = r.isIrreducible (a);
   bool isComp  = r.isComposit (a);

   DEB2
   {
      if (isUnit)   cout << ", unit";
      if (isPrime)  cout << ", prime";
      if (isComp)   cout << ", composit";
   }

   if (isIrred != isPrime)  error ("isIrreducible() != isPrime()");

   if (is0 + isUnit + isPrime + isComp != 1)
   {
      error ("Classifiaction ambiguous!");
   }

   if (isUnit)
   {
      T recip = r.unitRecip (a);
      if (r.size() && r.index(recip) >= r.size())  error ("b^-1 invalid");

      DEB2
      {
         cout << ", a^-1 = ";
         r.printShort (cout, recip);
      }

      if (! r.is1 (r.mul (a, recip)))
      {
         error ("Reciprocal of unit is wrong!");
      }
   }
}

template<class R>
void doTest
  (const R &r, const typename R::type& a, const typename R::type& b, ufd_tag)
{
   doTest (r, a, b, domain_tag());

   if (r.is0(a) || r.is0(b))  return;

   typename R::type prod = r.mul (a,b);

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
         if (! r.isComposit (prod))
         {
            error ("product of two non-units must be composit!");
         }
      }
   }
}


/**
 *  Euclidean Ring
 */

static unsigned lastNorm = 0;

template<class R>
void doTest (const R &r, const typename R::type& a, euclidean_tag)
{
   typedef typename R::type T;

   doTest (r, a, ufd_tag());

   unsigned norm = r.norm (a);

   DEB2
   {
      cout << ", norm = " << norm;
   }

   if (r.is0(a) != (norm == 0))  error ("norm(a)==0 iff a==0 failed!");
   if (r.isUnit(a) != (norm == 1))  error ("norm(a)==1 iff isUnit(a) failed!");
   if (norm < lastNorm)  error ("norm() decreased!");
   lastNorm = norm;

   if (r.is0 (a))  return;

   unsigned size = r.size();
   unsigned num = r.numOfRemainders (a);

   DEB2
   {
      cout << ", numRem: ";
      if (num) cout << num; else cout << "infinite";
   }

   if (! size)  return;

   if (num == 0 || num >= size)  error ("numOfRemainders() >= size()!");
}

template<class R>
void doTest
   (const R &r, const typename R::type& a, const typename R::type& b,
    euclidean_tag)
{
   doTest (r, a, b, ufd_tag());

   lastNorm = 0;

   // norm

   unsigned normA  = r.norm(a);
   unsigned normB  = r.norm(b);
   unsigned normAB = r.norm(r.mul(a,b));

   if (normA * normB < normAB)
   {
      error ("norm(a)norm(b) < norm(ab)!");
   }

   if (r.is0 (b))  return;

   if (normA > normAB)
   {
      error ("norm(ab) < norm(a)!");
   }

   // div()

   typename R::type q, re;
   r.div (a, b, q, re);

   // Make sure, operators don't leave domain

   if (r.size())
   {
      if (r.index (q)  >= r.size())  error ("quot() invalid!");
      if (r.index (re) >= r.size())  error ("rem() invalid!");
   }

   DEB2
   {
      cout << ", quotient = ";
      r.printShort (cout, q);
      cout << ", remainder = ";
      r.printShort (cout, re);
   }

   if (r.add (r.mul (b, q), re) != a)  error ("div() fails!");

   // rem() and quot()

   typename R::type x = r.rem(a,b);
   if (x != re)  error ("rem() inconsistent with div()!");

   x = r.quot(a,b);
   if (x != q)  error ("quot() inconsistent with div()!");

   if (r.norm (re) >= r.norm (b))  error ("norm(rem()) >= norm(b)!");
}


/**
 *  Field
 */

template<class R>
void doTest (const R &f, field_tag)
{
   doTest (f, euclidean_tag());

   unsigned c = f.characteristic();
   unsigned s = f.size();
   DEB2 cout << "  Char: " << c << endl;

   if (s == 0 && c != 0) error ("characteristic != 0 in infinite field!");

   if (s != 0)
   {
      unsigned power;
      unsigned prime;
      if (! Prime::isPrimePower (s, prime, power))
      {
         error ("Size of finite field is not a prime power!");
      }
      else
      {
         if (c != prime)  error ("characterisitc wrong!");
      }
   }
}

template<class R>
void doTest (const R &f, const typename R::type& a, field_tag)
{
   typedef typename R::type T;

   doTest (f, a, euclidean_tag());

   if (f.is0 (a))  return;

   // multiplicative inverse

   if (! f.isUnit (a))  error ("Each elemnt != 0 must be a unit!");
   typename R::type r = f.recip(a);

   if (r != f.unitRecip(a))  error ("recip() != unitRecip()");
   if (! f.is1 (f.mul (a, r)))  error ("a * (a^-1) != 1");

   // times characteristic

   if (! f.is0 (f.times (a, f.characteristic())))
   {
      error ("times(a,characteristic) != 0!");
   }
}

template<class R>
void doTest
   (const R &f, const typename R::type& a, const typename R::type& b,
    field_tag)
{
   doTest (f, a, b, euclidean_tag());

   if (f.is0(b))  return;

   typename R::type q = f.div(a, b);
   typename R::type r = f.recip (b);
   if (f.mul (a, r) != q)  error ("div() fails");

   typename R::type x;

   x = a;
   f.divBy (x, b);
   if (x != q) error ("divBy() fails");

   if (!f.is0 (f.rem (a,b)))  error ("rem(a,b) is not 0!");
   if (q != f.quot (a,b))  error ("div(a,b) != quot(a,b)");
}


/**
 *  Cyclic
 */

template<class R>
void doTest (const R &f, cyclic_tag)
{
   doTest (f, field_tag());

   unsigned c = f.characteristic();
   unsigned s = f.size();

   if (s == 0) error ("Infinite field can not have cyclic additive group!");

   if (s != c || ! Prime::test (s))  error ("Additive group not cyclic!");
}

/*************  Test procedure for vector spaces  ****************************/


template<class V>
void doTest (const V& v, vectorspace_tag)
{
   typedef typename V::scalar_algebra A;
   typedef typename V::type T;
   typedef typename V::scalar_type ST;
   const A alg = v.getScalarAlgebra();
   const unsigned dim = v.dimension();

   DEB2
   {
      cout << "  Dim:  " << dim << endl
           << "  Alg:  " << alg << endl;
   }

   if (alg.size() == 0 && v.size() != 0)
   {
      error("Finite vector space over infinte field!");
   }
   if (alg.size() != 0 && v.size() == 0)
   {
      error("Infinite vector space over finte field!");
   }
   if (alg.size() != 0 && v.size() != 0)
   {
      if (v.size() != std::numeric_limits<unsigned>::max())
      {
         if (v.size() != powInt (alg.size(), dim))
            error ("dim (F^n) != (dim F)^n");
      }
   }

   Array<ST> coords (dim);
   T x = v.zero();
   v.toCoord (x, &coords[0]);

   for (unsigned i = 0; i < dim; ++i)
   {
      if (! alg.is0 (coords[i]))  error ("zero() has non-zero coordinates!");
   }
}

template<class V>
void doTest (const V& v, const typename V::type& a, vectorspace_tag)
{
   typedef typename V::scalar_algebra A;
   typedef typename V::type T;
   typedef typename V::scalar_type ST;
   const A alg = v.getScalarAlgebra();
   const unsigned dim = v.dimension();

   // toCoords()
   
   Array<ST> coords (dim);
   v.toCoord (a, &coords[0]);

   T aCopy = a;

   DEB2
   {
      cout << ", coord = ";
      for (unsigned i = 0; i < dim; ++i)
      {
         if (i > 0)  cout << ',';
         alg.printShort (cout, coords[i]);
      }
   }

   for (unsigned i = 0; i < dim; ++i)
   {
      if (alg.size() && alg.index(coords[i]) >= alg.size())
         error ("Coordinate invalid!");
      if (coords[i] != v.coord(a,i))  error("coord() and toCoords() mismatch!");
      typename V::scalar_reference ref = v.coord(aCopy, i);
      if (coords[i] != ref)  error("coord() and toCoords() mismatch!");
   }

   // fromCoords()

   T b;
   v.fromCoord (b, &coords[0]);
   if (a != b)  error ("fromCoord() broken!");

   // coords()

   b = v.zero();
   for (unsigned i = 0; i < dim; ++i)  v.coord(b, i) = coords[i];
   if (a != b)  error ("Assignment to coord() broken!");

   b = v.zero();
   for (int i = dim - 1; i >= 0; --i)  v.coord(b, i) = coords[i];
   if (a != b)  error ("Assignment to coord() broken!");

   // Multiplication with 1

   b = v.mul (a, alg.one());

   if (b != a)  error ("1v != v!");

   // Multiplication with 0

   b = v.mul (a, alg.zero());

   if (! v.is0 (b))  error ("0v != 0!");

   //  neg()

   b = v.neg (a);
   for (unsigned i = 0; i < dim; ++i)
   {
      if (v.coord(b,i) != alg.neg (v.coord(a,i)))
      {
         error ("neg() for vector does not match neg() for coordinates!");
      }
   }

   // mul()

   if (sqr (v.index(a)) < SIZE)
   {
      unsigned SS = unsigned (sqrt(double(SIZE)));
      if (alg.size())  SS = std::min (SS, alg.size());

      for (unsigned l = 0; l < SS; ++l)
      {
         ST scal = alg.element (l);

         b = v.mul (a, scal);

         for (unsigned i = 0; i < dim; ++i)
         {
            if (v.coord(b,i) != alg.mul (v.coord(a,i), scal))
            {
               error ("mul() for vector does not match mul() for coordinates!");
            }
         }
      }
   }

   // Associativity and Distributivity (scalar, scalar, vector)

   if (cube (v.index(a)) < SIZE)
   {
      unsigned SS = unsigned (pow(double(SIZE), 1.0/3.0));
      if (alg.size())  SS = std::min (SS, alg.size());

      for (unsigned l1 = 0; l1 < SS; ++l1)
      {
         ST scal1 = alg.element (l1);

         for (unsigned l2 = 0; l2 < SS; ++l2)
         {
            ST scal2 = alg.element (l2);

            if (   v.mul (a, alg.mul (scal1, scal2))
                != v.mul (v.mul (a, scal1), scal2))
            {
               error ("(ab)v != a(bv)!");
            }

            if (   v.mul (a, alg.add (scal1, scal2))
                != v.add (v.mul(a, scal1), v.mul(a, scal2)))
            {
               error ("(a+b)v != av + bv!");
            }
         }
      }
   }
}

template<class V>
void doTest
   (const V& v, const typename V::type& a, const typename V::type& b,
    vectorspace_tag)
{
   typedef typename V::scalar_algebra A;
   typedef typename V::type T;
   typedef typename V::scalar_type ST;
   const A alg = v.getScalarAlgebra();
   const unsigned dim = v.dimension();

   // add()

   T sum = v.add (a, b);

   for (unsigned i = 0; i < dim; ++i)
   {
      if (v.coord(sum,i) != alg.add (v.coord(a,i), v.coord(b,i)))
      {
         error ("add() for vector does not match add() for coordinates!");
      }
   }

   // Distributivity (scalar, vector, vector)

   if (cube (v.index(a)) < SIZE && cube (v.index(b)) < SIZE)
   {
      unsigned SS = unsigned (pow(double(SIZE), 1.0/3.0));
      if (alg.size())  SS = std::min (SS, alg.size());

      for (unsigned l = 0; l < SS; ++l)
      {
         ST scal = alg.element (l);

         if (   v.mul (v.add(a,b), scal)
             != v.add (v.mul (a, scal), v.mul (b, scal)))
         {
            error ("a (v1+v2) != (a v1) + (b + v2)!");
         }
      }
   }
}

template<class V>
void doTest
   (const V& v, const typename V::type& a,
                const typename V::type& b,
                const typename V::type& c, vectorspace_tag)
{}

/*************  Test procedure for groups  ***********************************/

template<typename A>
void doTests(A r, group_tag)
{
   typedef typename A::algebra_category cat;
   typedef typename A::type T;

   const T zero = r.zero();
   if (r.size() && r.index(zero) >= r.size())  error ("index(zero) invalid");
   DEB2 { cout << "  Zero: "; r.printShort (cout, zero); cout << endl; }

   doTest (r, cat());

   unsigned S = SIZE;
   if (r.size() > 0)  S = std::min (S, r.size());

   // One variable

   for (unsigned i = 0; i < S; i++)
   {
      T a = r.element (i);

      DEB1 { cout << "  Checking: " << setw(W); r.printShort (cout, a); }

      // element() and index()

      DEB2 cout << ", index = " << i;
      if (r.index(a) != i)  error ("index(a) != i");

      // is0()

      bool is0 = r.is0(a);
      DEB2  if (is0)  cout << ", is zero";
      if (is0 != (a == zero))  error ("is0() fails");
      
      // zero

      if (r.add (a, zero) != a)  error ("a+0 != a");

      // additive inverse

      T neg = r.neg(a);
      DEB2
      {
         cout << ", -a = ";
         r.printShort (cout, neg);
      }

      if (r.size() && r.index(neg) >= r.size())  error ("-a invalid");
      if (! r.is0 (r.add (a, neg)))  error ("a + (-a) != 0");

      {
         // negate()

         T x = a;
         r.negate (x);
         if (x != neg)  error ("negate(a) != -a");
      }

      // check times()

      if (sqr (i) < SIZE)
      {
         DEB3 cout << ", times = ";

         unsigned SS = unsigned (sqrt(double(SIZE)));

         T x = zero;
         for (unsigned l = 0; l < SS; ++l)
         {
            T sum = r.times(a,l);
            DEB3 { r.printShort (cout, sum); cout << '/'; }
            if (x != sum)  error ("a+...+a != times (a,k)");
            r.addTo (x, a);
         }
      }

      doTest (r, a, cat());

      DEB1 cout << endl;
   }

   // two variables

   S = unsigned(sqrt(double(SIZE)));
   if (r.size() > 0)  S = std::min (S, r.size());

   for (unsigned i = 0; i < S; i++)
   {
      const T a = r.element (i);

      for (unsigned j = 0; j < S; j++)
      {
         const T b = r.element (j);

         DEB1
         {
            cout << "  Checking: "
                 << setw(W); r.printShort (cout, a);
            cout << setw(W); r.printShort (cout, b);
         }

         // Make sure, == makes sense

         if ((a == b) != (i == j))  error ("a==b incorrect");
         if ((a != b) != (i != j))  error ("a!=b incorrect");
         
         // Make sure, operators dont leave domain

         T add = r.add (a,b);
         T sub = r.sub (a,b);

         DEB2
         {
            cout << ", a+b = ";
            r.printShort (cout, add);
            cout << ", a-b = ";
            r.printShort (cout, sub);
         }

         if (r.size())
         {
            if (r.index (add) >= r.size())  error ("a+b invalid");
            if (r.index (sub) >= r.size())  error ("a-b invalid");
         }
   
         // Commutativity

         if (add != r.add (b,a))  error ("a+b != b+a");

         // derived operations

         T negB = r.neg (b);
         if (r.add (a, negB) != sub)  error ("sub() fails");

         T x;

         x = a;
         r.addTo (x, b);
         if (x != add) error ("addTo() fails");

         x = a;
         r.subFrom (x, b);
         if (x != sub) error ("subFrom() fails");

         doTest (r, a, b, cat());

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

            DEB1
            {
               cout << "  Checking: "
                    << setw(W); r.printShort (cout, a);
               cout << setw(W); r.printShort (cout, b);
               cout << setw(W); r.printShort (cout, c);
            }

            // Associativity

            if (r.add(r.add(a,b), c) != r.add(a, r.add(b,c)))
                  error("(a+b) + c != a + (b+c)");

            doTest (r, a, b, c, cat());

            DEB1 cout << endl;
         }
      }
   }
}


/*************  All algebraic structures  ************************************/

const char* typeName (ring_tag)        { return "ring"; }
const char* typeName (domain_tag)      { return "integral domain"; }
const char* typeName (ufd_tag)         { return "UFD"; }
const char* typeName (euclidean_tag)   { return "Euclidean ring"; }
const char* typeName (field_tag)       { return "field"; }
const char* typeName (cyclic_tag)      { return "cyclic field"; }
const char* typeName (vectorspace_tag) { return "vector space"; }

const char* polyName (nopolynomial_tag) { return ""; }
const char* polyName (  polynomial_tag) { return ", polynomials"; }

template<typename A>
void doTests(A r, const char* type)
{
   typedef typename A::algebra_category cat;
   typedef typename A::polynomial_category poly;
   typedef typename A::type T;

   NORMAL cout << "Testing " << typeName(cat()) << polyName(poly())
               << " \"" << r << "\" (" << type << ")..." << endl;

   DEB2
   {
      cout << "  Size: ";
      if (r.size() == 0) cout << "infinite"; else cout << r.size();
      cout << endl;
   }
   
   doTests (r, cat());
}
   

/************  Command line arguments  ***************************************/

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


/************  Main Program  *************************************************/

/**
 *  test()
 *
 *  Main program
 *
 *  Executes all test cases
 */

void test (int argc, char**)
{
   if (argc)  usage();

   typedef unsigned char u8;

   {  // Integers and integer polynomials

      IntegerRing<int> intRing;
      doTests (intRing, "IntegerRing<int>");
      doTests (PolynomialRing<IntegerRing<int> > (intRing),
               "PolynomialRing<IntegerRing<int> >");
   }
   {  // GF2 and GF2 polynomials

      GF2 gf2;
      doTests (gf2, "GF2");
      doTests (PolynomialRing<GF2> (gf2), "PolynomialRing<GF2>");
   }
   doTests (Polynomial2Ring<u32>(), "Polynomial2Ring<u32>");

   for (unsigned d = 1; d <= 10; ++d)
   {
      doTests (GF2VectorSpace<u32> (d), "GF2VectorSpace<u32>");
   }

   // Algebraic structures with  n  elements

   for (unsigned i = 2; i <= 17; ++i)
   {
      // Finite field of the given size

      unsigned p;
      unsigned e;
      if (Prime::isPrimePower (i, p, e))
      {
         unsigned maxDim   = digitsRepresentable (static_cast<u8>(i));
         unsigned maxDim32 = digitsRepresentable (static_cast<u32>(i));

         {
            LookupGaloisField<u8> gf (i);
            doTests (gf, "LookupGaloisField<>");
            doTests (PolynomialRing<LookupField<u8> > (gf),
                     "PolynomialRing<LookupField<>>");
            doTests (OneDimVectorSpace<LookupField<u8> > (gf),
                     "OneDimVectorSpace<LookupField<>>");

            for (unsigned j = 1; j <= maxDim; ++j)
            {
               doTests (LookupVectorSpace<u8,u8>(gf, j), "LookupVectorSpace<>");
            }
         }
        
         // optimized version for primes

         if (Prime::test (i))
         {
            LookupGaloisFieldPrime<u8> gf (i);
            doTests (gf, "LookupGaloisFieldPrime<>");
            doTests (PolynomialRing<LookupFieldPrime<u8> > (gf),
                     "PolynomialRing<LookupFieldPrime<>>");
            doTests (OneDimVectorSpace<LookupFieldPrime<u8> > (gf),
                     "OneDimVectorSpace<LookupFieldPrime<>>");
         }
        
         // optimized version for powers of 2

         if (p == 2)
         {
            LookupGaloisFieldPow2<u8> gf (i);
            doTests (gf, "LookupGaloisFieldPow2<>");
            doTests (PolynomialRing<LookupFieldPow2<u8> > (gf),
                     "PolynomialRing<LookupFieldPow2<>>");
            doTests (OneDimVectorSpace<LookupFieldPow2<u8> > (gf),
                     "OneDimVectorSpace<LookupFieldPow2<>>");

            for (unsigned j = 1; j <= maxDim; ++j)
            {
               doTests (LookupVectorSpacePow2<u8,u8>(gf,j),
                        "LookupVectorSpacePow2<>");
            }

            for (unsigned j = 1; j <= maxDim32; ++j)
            {
               doTests (VectorSpacePow2<u32>(gf,j), "VectorSpacePow2<>");
            }
         }
      }

      if (Prime::test (i))
      {
         ModularArithField<unsigned short> field (i);
         doTests (field, "ModularArithField<>");

         PolynomialRing<ModularArithField<unsigned short> > poly (field);
         doTests (poly, "PolynomialRing<ModularArithField<>");
         
         unsigned short coef [] = {0, 0, 1};
         PolynomialRing<ModularArithField<unsigned short> >::type
            p (coef, coef + sizeof(coef)/sizeof(unsigned short));
         ModularArith<PolynomialRing<
            ModularArithField<unsigned short> > > mod (poly, p);
         doTests (mod, "ModularArithmetic<PolynomialRing<ModularArithField<>>>");

         for (unsigned j = 1; j <= 4; ++j)
         {
            doTests (GaloisField<u8> (i, j), "GaloisField<>");
         }
      }
      else
      {
         ModularArith<unsigned short> modInt (i);
         doTests (modInt, "ModularArith<>");
         doTests (PolynomialRing<ModularArith<unsigned short> > (modInt),
                  "PolynomialRing<ModularArith<>>");
      }
   }

   ModularArith<unsigned short> modInt (213);
   doTests (modInt, "ModularArith<>");

   ModularArithField<unsigned short> gf1999 (1999);
   doTests (gf1999, "ModularArithField<>");
}

