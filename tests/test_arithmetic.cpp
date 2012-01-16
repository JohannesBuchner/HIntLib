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
#include <utility>

#include <HIntLib/polynomial2.h>
#include <HIntLib/integerring.h>
#include <HIntLib/galoisfield.h>
#include <HIntLib/lookupfield.h>
#include <HIntLib/quotientfield.h>
#include <HIntLib/onedimvectorspace.h>
#include <HIntLib/array.h>

#include "test.h"

using std::cerr;
using std::cout;
using std::endl;
using std::setw;

using namespace HIntLib;

#ifdef HINTLIB_NAMESPACE_REL_OPS
using namespace std::rel_ops;
#else
using namespace std;
#endif

unsigned SIZE = 100;
bool FLUSH = false;
unsigned W = 15;


/**
 *  fl();
 */

inline void fl ()
{
   DEB2 if (FLUSH)  cout << std::flush;
}

/**
 *  ensure Equal Type()
 */

template<typename X>
inline
void ensureEqualType (const X*, const X*) {}


/*****************************************************************************/
/*************  Groups                                           *************/
/*****************************************************************************/


template<class G>
void doTest (const G& g, group_tag)
{
   // doTest (g, XXXX());

   typedef typename G::type T;

   const T zero = T();
   DEB2
   {
      cout << "  Zero: ";
      g.printShort (cout, zero);
   }

   unsigned ind = g.index(zero);

   DEB3  cout << " (ind = " << ind << ')';

   if (g.size() && ind >= g.size())  error ("index(zero) invalid");

   if (ind != 0 || ! g.is0 (g.element(0)))
   {
      error ("index(zero) != 0!");
   }
   DEB2  cout << endl;
}

template<class G>
void doTest (const G& g, const typename G::type& a, group_tag)
{
   // doTest (g, a, XXXX());

   typedef typename G::type T;

   // Copy and assignment

   T x (a);
   if (x != a)  error ("Copy constructor or operator==() broken!");
   x = a;
   if (x != a)  error ("Assignment, or operator==() broken!");
   
   const T zero = T();
   const unsigned s = g.size();

   // is0()

   bool is0 = g.is0(a);
   DEB2  if (is0)  cout << ", zero";
   if (is0 != (a == zero))  error ("is0() fails");
   
   // zero

   if (g.add (a, zero) != a)  error ("a+0 != a");

   // additive inverse

   T neg = g.neg(a);
   DEB2
   {
      cout << ", -a = ";
      g.printShort (cout, neg);
   }

   if (s && g.index(neg) >= s)  error ("-a invalid");
   if (! g.is0 (g.add (a, neg)))  error ("a + (-a) != 0");

   {
      // negate()

      T x = a;
      g.negate (x);
      if (x != neg)  error ("negate(a) != -a");
   }
   
   fl();

   // check times()

   if (sqr (g.index(a)) < SIZE)
   {
      DEB3 cout << ", times = ";

      unsigned SS = unsigned (sqrt(double(SIZE)));

      T x = zero;
      for (unsigned l = 0; l < SS; ++l)
      {
         T sum = g.times(a,l);
         DEB3 { g.printShort (cout, sum); cout << '/'; }
         if (x != sum)
         {
            error ("a+...+a != times (a,k)");
            cout << "a+...+a = ";
            g.printShort (cout, x);
            cout << ",  times (a," << l << ") = ";
            g.printShort (cout, sum);
            cout << ",  diff = ";
            g.printShort (cout, g.sub (x, sum));
            cout << endl;
         }
         g.addTo (x, a);
      }
   }

   fl();

   // additive order of the element

   if (sqr(s) < SIZE)
   {
      unsigned o = g.additiveOrder (a);
      DEB2
      {
         cout << ", add.order: ";
         if (o) cout << o; else cout << "inf";
      }

      if (o == 0)
      {
         if (s != 0)  error ("infinite additive order in finite group!");

         return;
      }

      if (s && o > s)  error ("additiveOrder() larger than size()");
      if (s && s % o != 0)  error ("additiveOrder() does not divide size()!");

      T x = zero;
      for (unsigned i = 1; i < o; ++i)
      {
         g.addTo (x, a);
         if (g.is0(x))  error ("additveOrder() is too high!");
      }
      g.addTo (x,a);
      if (! g.is0(x))  error ("additiveOrder() is too low!");
   }

   fl();
}

template<class G>
void doTest
   (const G& g, const typename G::type& a, const typename G::type& b, group_tag)
{
   // doTest (g, a, b, XXXX());

   typedef typename G::type T;

   // Make sure, operators dont leave domain

   T add = g.add (a,b);
   T sub = g.sub (a,b);

   DEB2
   {
      cout << ", a+b = ";
      g.printShort (cout, add);
      cout << ", a-b = ";
      g.printShort (cout, sub);
   }

   if (g.size())
   {
      if (g.index (add) >= g.size())  error ("a+b invalid");
      if (g.index (sub) >= g.size())  error ("a-b invalid");
   }
   
   // Commutativity

   if (add != g.add (b,a))  error ("a+b != b+a");

   // derived operations

   T negB = g.neg (b);
   if (g.add (a, negB) != sub)  error ("sub() fails");

   T x;

   x = a;
   g.addTo (x, b);
   if (x != add) error ("addTo() fails");

   x = a;
   g.subFrom (x, b);
   if (x != sub) error ("subFrom() fails");

   fl();
}

template<class G>
void doTest
   (const G& g, const typename G::type& a,
                const typename G::type& b,
                const typename G::type& c, group_tag)
{
   // doTest (g, a, b, c, XXXX());

   // Associativity of addition

   if (g.add(g.add(a,b), c) != g.add(a, g.add(b,c)))
   {
      error("(a+b) + c != a + (b+c)");
   }

   fl();
}


/*****************************************************************************/
/*************  Rings                                            *************/
/*****************************************************************************/


template<class R>
void doTest (const R &r, ring_tag)
{
   doTest (r, group_tag());

   const typename R::type one = r.one();
   if (r.size() && r.index(one) >= r.size())  error ("index(one) invalid");
   DEB2 { cout << "  One:  "; r.printShort (cout, one); cout << endl; }

   if (r.index(one) != 1 || ! r.is1(r.element(1)))  error ("index(one) != 1!");
}

template<class R>
void doTest (const R &r, const typename R::type& a, ring_tag)
{
   doTest (r, a, group_tag());

   typedef typename R::type T;

   // is1(), one()

   const T one = r.one();

   bool is1 = r.is1(a);
   DEB2 if (is1)  cout << ", one";
   if (is1 != (a == one ))  error ("is1() fails");
      
   // neutral

   if (r.mul (a, one) != a)  error ("a*1 != a");

   // check power()

   if (sqr(r.index(a)) < SIZE)
   {
      DEB3 cout << ", power = ";

      unsigned SS = unsigned (sqrt(double (SIZE)));

      if (SS > 10) SS = 7;

      T x = a;
      for (unsigned l = 1; l < SS; ++l)
      {
         try
         {
            T p = r.power (a, l);
         
            DEB3
            {
               r.printShort (cout, p);
               DEB4  { cout << '='; r.printShort (cout, x); }
               cout << '/';
            }
            if (x != p)
            {
               error ("a*...*a != power (a,k)");
               cout << "a*...*a = ";
               r.print (cout, x);
               cout << ",  power (a," << l << ") = ";
               r.print (cout, p);
               cout << endl;
            }
            r.mulBy (x, a);
         }
         catch (Overflow &)
         {
            break;
         }
      }
   }

   fl();
}

template<class R>
void doTest
   (const R& r, const typename R::type& a, const typename R::type& b, ring_tag)
{
   doTest (r, a, b, group_tag());

   // mul()

   typedef typename R::type T;

   const T mul = r.mul (a,b);

   DEB2
   {
      cout << ", a*b = ";
      r.printShort (cout, mul);
   }

   if (r.size())  if (r.index (mul) >= r.size())  error ("a*b invalid");

   // Commutativity

   if (mul != r.mul (b,a))  error ("a*b != b*a");

   T x = a;
   r.mulBy (x, b);
   if (x != mul) error ("mulBy() fails");

   fl();
}

template<class R>
void doTest
   (const R& r, const typename R::type& a,
                const typename R::type& b,
                const typename R::type& c, ring_tag)
{
   doTest (r, a, b, c, group_tag());

   // Associativity of multiplication

   if (r.mul(r.mul(a,b), c) != r.mul(a, r.mul(b,c)))
         error("(a*b) * c != a * (b*c)");

   // Distributivity

   if (r.mul (a, r.add(b,c)) != r.add (r.mul(a,b), r.mul(a,c)))
         error("a*(b+c) != (a*b)+(a*c)");
   if (r.mul (r.add(a,b), c) != r.add (r.mul(a,c), r.mul(b,c)))
      error("(a+b)*c != (a*c)+(b*c)");

   fl();
}


/*****************************************************************************/
/*************  Domains                                          *************/
/*****************************************************************************/


template<class R>
void doTest (const R &r, domain_tag)
{
   doTest (r, ring_tag());

   unsigned c = r.characteristic();
   DEB2
   {
      cout << "  Char: ";
      if (c) cout << c; else cout << "inf";
      cout << endl;
   }

   if (c != 0)
   {
      if (! Prime::test (c))
      {
         error ("Finite characteristic not a prime number!");
         return;
      }
   }
}

template<class R>
void doTest (const R &r, const typename R::type& a, domain_tag)
{
   doTest (r, a, ring_tag());

   typedef typename R::type T;

   unsigned o = r.additiveOrder (a);
  
   if (r.is0 (a))  return;

   // additive order must match characteristic

   if (o != r.characteristic())
   {
      error ("additiveOrder() != characteristic()!");
   }

   // multiplicative order of the element

   if (! r.is0 (a))
   {
      unsigned s = r.size();

      if (sqr(s) < SIZE)
      {
         unsigned ord = r.order (a);

         DEB2
         {
            cout << ", mult.order: ";
            if (ord) cout << ord; else cout << "inf";
         }

         if (ord == 0)
         {
            if (s != 0)
            {
               error ("infinite multiplicative order in finite ring!");
            }
         }
         else  // ord finite
         {
            if (s && ord > s-1)  error ("order() larger than size()-1");
            if (s && (s-1) % ord != 0)
            {
               error ("order() does not divide size()-1!");
            }

            T x = r.one();
            for (unsigned i = 1; i < ord; ++i)
            {
               r.mulBy (x, a);
               if (r.is1(x))  error ("order() is too high!");
            }
            r.mulBy (x,a);
            if (! r.is1(x))  error ("order() is too low!");
         }
      }
   }

   fl();
}

template<class R>
void doTest
  (const R &r, const typename R::type& a, const typename R::type& b, domain_tag)
{
   doTest (r, a, b, ring_tag());

   if (r.is0(a) || r.is0(b))  return;
   if (r.is0(r.mul (a, b)))  error ("a * b = 0!");

   fl();
}


/*****************************************************************************/
/*************  UFDs                                             *************/
/*****************************************************************************/


template<class R>
void doTest (const R &r, ufd_tag)
{
   typedef typename R::type T;
   typedef typename R::unit_type UT;

   doTest (r, domain_tag());

   unsigned numUnits = r.numUnits();
   unsigned size = r.size();

   DEB2
   {
      cout << "  # units: ";
      if (numUnits)  cout << numUnits; else cout << "inf";
      cout << endl;
   }

   if (size)
   {
      if (numUnits == 0 || numUnits >= size)
      {
         error ("More units than elements!");
      }
   }

   // make sure, unit_type works

   {
      UT x = r.unitElement(0);
      UT def = x;
      if (def != x)  error ("unit_type broken");
      def = x;
      if (def != x)  error ("unit_type broken");
   }

   DEB3 cout << "Units: ";
   bool first = true;

   unsigned ub = numUnits ? std::min (SIZE, numUnits) : SIZE;
   for (unsigned i = 0; i < ub; ++i)
   {
      UT unit = r.unitElement (i);
      T  unitElem = r.fromUnit (unit);

      DEB3
      {
         if (! first)  cout << ", ";
         r.printShort (cout, unitElem);
         first = false;
      }

      if (! r.isUnit (unitElem))
      {
         error ("result of unitElement() not isUnit()!");
      }
      
      if (r.unitIndex (unit) != i)
      {
         error ("unitElement() and unitIndex() do not match!");
      }

      if (size)
      {
         if (r.element(i + 1) != unitElem)
         {
            error ("element() and unitElement() do not match!");
         }
      }
   }

   DEB3 cout << endl;

   // Prime generator

   typename R::PrimeGenerator pg (r);

   DEB3 cout << "Primes: ";
   first = true;

   for (unsigned i = 0; i < SIZE; ++i)
   {
      T p = pg.next();

      DEB3
      {
         if (! first)  cout << ", ";
         r.printShort (cout, p);
         first = false;
      }

      if (! r.isCanonical (p))
      {
         return error ("PrimeGenerator::next() does not return canonical!");
      }

      if (! r.isPrime (p))
      {
         return error ("PrimeGenerator::next() does not return prime!");
      }
   }

   DEB3 cout << endl;
}

template<class R>
void doTest (const R &r, const typename R::type& a, ufd_tag)
{
   typedef typename R::type T;
   typedef typename R::unit_type UT;

   doTest (r, a, domain_tag());

   bool is0         = r.is0 (a);
   bool isUnit      = r.isUnit (a);
   bool isPrime     = r.isPrime (a);
   bool isComp      = r.isComposit (a);
   bool isCanonical = r.isCanonical (a);

   DEB2
   {
      if (isUnit)      cout << ", unit";
      if (isPrime)     cout << ", prime";
      if (isComp)      cout << ", composit";
      if (isCanonical) cout << ", canonical";
   }

   if (is0 + isUnit + isPrime + isComp != 1)
   {
      error ("Classifiaction ambiguous!");
   }

   if (is0 && ! isCanonical)  error ("0 should be canonical!");

   if (isUnit)
   {
      UT unit = r.toUnit (a);
      UT recip = r.unitRecip (unit);
      T recipElem = r.fromUnit (recip);

      DEB2
      {
         cout << ", a^-1 = ";
         r.printShort (cout, recipElem);
      }
      if (r.size() && r.index(recipElem) >= r.size())  error ("a^-1 invalid");

      if (! r.is1 (r.mul (a, recipElem)))
      {
         error ("Reciprocal of unit is wrong!");
      }
      if (! r.is1 (r.mulUnit (a, recip)))
      {
         error ("Reciprocal of unit is wrong!");
      }
      if (! r.is1 (r.fromUnit (r.mulUnit (unit, recip))))
      {
         error ("Reciprocal of unit is wrong!");
      }

      unsigned numUnits = r.numUnits();
      unsigned index = r.unitIndex (unit);
      DEB2 cout << ", unit index = " << index;
      if (numUnits && index >= numUnits)  error ("unitIndex() invalid!");
   }

   // make canonical

   {
      T copy = a;
      UT unit = r.makeCanonical (copy);
      T unitElem = r.fromUnit (unit);

      if (r.mulUnit (copy, unit) != a || r.mul (copy, unitElem) != a)
      {
         error ("makeCanonical() does not factor!");
      }

      if (! r.isUnit(unitElem))
      {
         error ("makeCanonical() does not return a unit!");
      }
      if (! r.isCanonical (copy))
      {
         error ("makeCanonical() does not create a canonical form!");
      }

      // Canonical form of units

      if (isUnit)
      {
         if (a != unitElem || ! r.is1 (copy))
         {
            error ("makeCanonical() of unit is broken!");
         }
      }

      // Canonical form of (non-)canonical forms

      if (isCanonical)
      {
         if (! r.is1(unitElem) || copy != a)
         {
            error ("makeCanonical() of canonical form invalid!");
         }
      }
      else // ! isCanonical
      {
         if (r.is1(unitElem) || copy == a)
         {
            error ("makeCanonical() for non-canonical form invalid!");
         }
      }
   }

   fl();
}

template<class R>
void doTest
  (const R &r, const typename R::type& a, const typename R::type& b, ufd_tag)
{
   doTest (r, a, b, domain_tag());

   if (r.is0(a) || r.is0(b))  return;

   typename R::type prod = r.mul (a,b);

   bool aIsUnit      = r.isUnit (a);
   bool bIsUnit      = r.isUnit (b);
   bool aIsCanonical = r.isCanonical (a);
   bool bIsCanonical = r.isCanonical (b);

   if (aIsCanonical && bIsCanonical)
   {
      if (! r.isCanonical (prod))
      {
         error ("Product of canonicals must be canonical!");
      }
   }

   if (aIsUnit && bIsUnit)
   {
      if (! r.isUnit (prod))  error ("Product of units must be a unit!");
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

   fl();
}


/*****************************************************************************/
/*************  Euclidean Ring                                   *************/
/*****************************************************************************/


unsigned lastNorm = 0;

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
      if (num) cout << num; else cout << "inf";
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

   if (! r.is0 (b))
   {
      if (normA > normAB)  error ("norm(ab) < norm(a)!");

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

      if (! (r.add (r.mul (b, q), re) == a))  error ("div() fails!");

      // rem() and quot()

      typename R::type x = r.rem(a,b);
      if (! (x == re))  error ("rem() inconsistent with div()!");

      x = r.quot(a,b);
      if (! (x == q))  error ("quot() inconsistent with div()!");

      if (r.norm (re) >= r.norm (b))  error ("norm(rem()) >= norm(b)!");
   }

   fl();
}


/*****************************************************************************/
/*************  Field                                            *************/
/*****************************************************************************/


template<class R>
void doTest (const R &f, field_tag)
{
   typedef typename R::type T;

   doTest (f, domain_tag());

   unsigned s = f.size();

   if (s == 0)  return;   // remaining stuff only for finite fields

   unsigned power;
   unsigned prime;
   if (! Prime::isPrimePower (s, prime, power))
   {
      error ("Size of finite field is not a prime power!");
   }
}

template<class R>
void doTest (const R &f, const typename R::type& a, field_tag)
{
   typedef typename R::type T;

   doTest (f, a, domain_tag());

   if (! f.is0 (a))
   {
      // multiplicative inverse

      T r = f.recip(a);
      DEB2
      {
         cout << ", a^-1 = ";
         f.printShort (cout, r);
      }

      if (! f.is1 (f.mul (a, r)))  error ("a * (a^-1) != 1");
   }

   fl();
}

template<class R>
void doTest
   (const R &f, const typename R::type& a, const typename R::type& b,
    field_tag)
{
   doTest (f, a, b, domain_tag());

   if (! f.is0(b))
   {
      typename R::type q = f.div(a, b);
      typename R::type r = f.recip (b);
      DEB2
      {
         cout << ", a/b = ";
         f.printShort (cout, q);
      }

      if (f.mul (a, r) != q)  error ("div() fails");

      typename R::type x;

      x = a;
      f.divBy (x, b);
      if (x != q) error ("divBy() fails");
   }

   fl();
}


/*****************************************************************************/
/*************  Finite Field                                     *************/
/*****************************************************************************/


unsigned numPrimitives = 0;

template<class R>
void doTest (const R &f, gf_tag)
{
   numPrimitives = 0;

   typedef typename R::type T;

   doTest (f, field_tag());

   unsigned s = f.size();
   unsigned c = f.characteristic();
   unsigned deg = f.extensionDegree();

   DEB2  cout << "  Degree: " << deg << endl;

   if (s == 0)
   {
      error ("size infinite in finite field!");
      return;
   }

   unsigned power;
   unsigned prime;
   if (! Prime::isPrimePower (s, prime, power))
   {
      error ("Size of finite field is not a prime power!");
   }

   if (c != prime)    error ("characteristic() and field size incompatible!");
   if (deg != power)  error ("extensionDegree() and field size incompatible!");
}

template<class R>
void doTest (const R &f, const typename R::type& a, gf_tag)
{
   typedef typename R::type T;

   doTest (f, a, field_tag());

   // isPrimitiveElement()

   if (! f.is0 (a))
   {
      unsigned s = f.size();

      if (sqr(s) < SIZE)
      {
         unsigned o = f.order (a);
         bool isPrimitiveElement = f.isPrimitiveElement(a);

         if (isPrimitiveElement)
         {
            DEB2  cout << ", primitive root";
            ++numPrimitives;
         }

         if (f.index (a) == s - 1)
         {
            if (numPrimitives != Prime::eulerPhi (s - 1))
            {
               error ("number of primitive elements is wrong!");
            }
         }
   
         if (isPrimitiveElement != (o == s - 1))
         {
            error ("isPrimitiveRoot() broken!");
         }
      }
   }

   fl();
}


/*****************************************************************************/
/*************  Cyclic Field                                     *************/
/*****************************************************************************/


template<class R>
void doTest (const R &f, cyclic_tag)
{
   doTest (f, gf_tag());

   unsigned c = f.characteristic();
   unsigned s = f.size();

   if (s != c || ! Prime::test (s))  error ("Additive group not cyclic!");
}


/*****************************************************************************/
/*************  Vectorspace                                      *************/
/*****************************************************************************/

template<class V>
void doTest (const V& v, vectorspace_tag)
{
   doTest (v, group_tag());
   
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

   ensureEqualType (static_cast<ST*>(0), static_cast<typename A::type*>(0));

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
      if (logInt (std::numeric_limits<unsigned>::max(), unsigned (alg.size()))
          <= int (v.dimension()))
      {
         if (v.size() != powInt (alg.size(), dim))
            error ("dim (F^n) != (dim F)^n");
      }
   }

   Array<ST> coords (dim);
   T x = T();
   v.toCoord (x, &coords[0]);

   for (unsigned i = 0; i < dim; ++i)
   {
      if (! alg.is0 (coords[i]))  error ("zero has non-zero coordinates!");
   }
}

template<class V>
void doTest (const V& v, const typename V::type& a, vectorspace_tag)
{
   doTest (v, a, group_tag());
   
   typedef typename V::scalar_algebra A;
   typedef typename V::type T;
   typedef typename V::scalar_type ST;
   const A& alg = v.getScalarAlgebra();
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

   b = T();
   for (unsigned i = 0; i < dim; ++i)  v.coord(b, i) = coords[i];
   if (a != b)  error ("Assignment to coord() broken!");

   b = T();
   for (int i = dim - 1; i >= 0; --i)  v.coord(b, i) = coords[i];
   if (a != b)  error ("Assignment to coord() broken!");

   // Multiplication with 1

   b = v.mul (a, alg.one());

   if (b != a)  error ("1v != v!");

   // Multiplication with 0

   b = v.mul (a, ST());

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

   fl();
}

template<class V>
void doTest
   (const V& v, const typename V::type& a, const typename V::type& b,
    vectorspace_tag)
{
   doTest (v, a, b, group_tag());
   
   typedef typename V::scalar_algebra A;
   typedef typename V::type T;
   typedef typename V::scalar_type ST;
   const A& alg = v.getScalarAlgebra();
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

   fl();
}


/*****************************************************************************/
/*************  Polynomials                                      *************/
/*****************************************************************************/

template<class R>
void doTest (const R&, nopolynomial_tag) {}

template<class R>
void doTest (const R& r, polynomial_tag)
{
   typedef typename R::coeff_algebra A;
   typedef typename R::type T;
   typedef typename R::coeff_type CT;

   const A alg = r.getCoeffAlgebra();
   DEB2  cout << "  Alg:  " << alg << endl;

   ensureEqualType (static_cast<CT*>(0), static_cast<typename A::type*>(0));
   ensureEqualType (static_cast<CT*>(0),
                    static_cast<typename T::coeff_type*>(0));

   // x()

   T x = r.x();
   DEB2
   {
      cout << "  x: ";
      r.print (cout, x);
      cout << endl;
   }

   if (x.degree() != 1 || ! alg.is1(x.lc()) || ! alg.is0 (x[0]))
   {
      error ("x() is invalid");
   }

   for (int i = 0; i < 5; ++i)
   {
      T xx = r.x (i);

      if (xx.degree() != i || ! alg.is1(xx.lc()))
      {
         error ("x(k) is invalid!");
      }
      for (int j = 0; j < i; ++j)
      {
         if (! alg.is0 (xx[j]))  error ("x(k) is invalid!");
      }
      
      T xxdiv1 = xx;
      xxdiv1.divByX (i);
      T xxdiv2 = xx;
      xxdiv2.divByX (i+1);

      if (! r.is1 (xxdiv1) || ! r.is0 (xxdiv2))  error ("divByX() broken!");

      T one = r.one();
      if (xx != one.mulByX (i))  error ("mulByX() broken!");
   }
}

template<class R>
void doTestPoly (const R&, const typename R::type&, ring_tag) {}

template<class R>
void doTestPoly (const R& r, const typename R::type& a, gf_tag)
{
   doTestPoly (r, a, field_tag());

   // isPrimitive

   bool isPrimitive = r.isPrimitive (a);

   DEB2 if (isPrimitive)  cout << ", primitive";
}

template<class R>
void doTest (const R&, const typename R::type&, nopolynomial_tag) {}

template<class R>
void doTest (const R& r, const typename R::type& a, polynomial_tag)
{
   typedef typename R::coeff_algebra A;
   typedef typename R::type T;
   typedef typename R::coeff_type CT;

   const A& alg = r.getCoeffAlgebra();

   // degree 

   int deg = a.degree();
   DEB2  cout << ", degree = " << deg;

   if (deg == -1)
   {
      if (! r.is0(a))  error ("Polynomial with deg=-1 is not 0!");
   }
   else
   {
      bool monic = r.isCanonical (a);
      DEB2  if (monic)  cout << ", monic";

      CT l = a.lc();
      DEB2
      {
         cout << ", lc = ";
         alg.printShort (cout, l);
      }

      if (monic != alg.is1(l))  error ("monic() is wrong!");
      if (l != a[a.degree()])  error ("lc() != leading coefficient!");
      if (alg.is0(l))  error ("lc() is zero!");
   }

   // Access

   Array<CT> coeff (deg + 1);
   a.toCoeff (&coeff[0]);

   DEB2 cout << ", coeff = /";
   for (int i = 0; i <= deg; ++i)
   {
      DEB2
      {
         alg.printShort (cout, a[i]);
         cout << '/';
      }
      
      if (coeff[i] != a[i])  error ("toCoeff() or operator[]() broken!");
   }

   // Constructor with iterators

   T b (&coeff[0], &coeff[deg+1]);

   if (! (a == b))  error ("Constructor with iterators broken!");

   // mulByX(), divByX()

   T div (a);
   T mul (a);
   
   for (int i = 0; i < 5; ++i)
   {
      T div1 (a);
      T mul1 (a);

      div1.divByX (i);
      mul1.mulByX (i);

      if (div != div1)  error ("divByX() broken!");
      if (mul != mul1)  error ("mulByX() broken!");

      div.divByX();
      mul.mulByX();
   }

   // derivative()

   T deriv = r.derivative (a);

   DEB2
   {
      cout << ", derivative = ";
      r.printShort (cout, deriv);
   }
   
   T deriv1 = r.derivative (r.neg (a));
   T deriv2 = r.neg (r.derivative(a));

   if (! (deriv1 == deriv2))
   {
      error ("derivative() inconsistent with negation!");
   }

   doTestPoly (r, a, typename A::algebra_category());

   // evaluate()

   if (sqr (r.index(a)) > SIZE)  return;

   unsigned ub = unsigned (sqrt(double(SIZE)));
   if (ub > alg.size())  ub = alg.size();
   
   for (unsigned i = 0; i < ub; ++i)
   {
      CT x = alg.element (i);
      CT eval = r.evaluate (a, x);

      CT res = CT();
      for (int i = a.degree(); i >= 0; --i)
      {
         res = alg.add (alg.mul (res, x), a[i]);
      }
      
      if (eval != res)  error ("evaluate() broken!");
   }

   fl();
}

template<class R>
void doTest (const R&, const typename R::type&,
                       const typename R::type&, nopolynomial_tag) {}

template<class R>
void doTest (const R& r, const typename R::type& a,
                         const typename R::type& b, polynomial_tag)
{
   typedef typename R::coeff_algebra A;
   typedef typename R::type T;
   typedef typename R::coeff_type CT;

   // derivative()

   T deriv1 = r.derivative (r.add (a,b));
   T deriv2 = r.add (r.derivative(a), r.derivative(b));

   if (! (deriv1 == deriv2))
   {
      error ("derivative() inconsistent with addition!");
   }

   deriv1 = r.derivative (r.mul (a,b));
   deriv2 = r.add (r.mul (a, r.derivative(b)), r.mul (b, r.derivative(a)));

   if (! (deriv1 == deriv2))
   {
      error ("derivative() inconsistent with multiplication!");
   }

   fl();
}


/*************  All algebraic structures  ************************************/

const char* typeName (group_tag)       { return "group"; }
const char* typeName (ring_tag)        { return "ring"; }
const char* typeName (domain_tag)      { return "integral domain"; }
const char* typeName (ufd_tag)         { return "UFD"; }
const char* typeName (euclidean_tag)   { return "Euclidean ring"; }
const char* typeName (integer_tag)     { return "Euclidean ring of integers"; }
const char* typeName (field_tag)       { return "field"; }
const char* typeName (rational_tag)    { return "field of rational numbers"; }
const char* typeName (real_tag)        { return "field of real numbers"; }
const char* typeName (complex_tag)     { return "field of complex numbers"; }
const char* typeName (gf_tag)          { return "finite field"; }
const char* typeName (cyclic_tag)      { return "cyclic field"; }
const char* typeName (vectorspace_tag) { return "vector space"; }

const char* polyName (nopolynomial_tag) { return ""; }
const char* polyName (  polynomial_tag) { return ", polynomials"; }


/*************  Main test procedure  *****************************************/

template<class A>
void doTests(A r, const char* type)
{
   typedef typename A::algebra_category cat;
   typedef typename A::polynomial_category pol;
   typedef typename A::type T;

   NORMAL cout << "Testing " << typeName(cat()) << polyName(pol())
               << " \"" << r << "\" (" << type << ")..." << endl;

   DEB1
   {
      cout << "  Size: ";
      if (r.size() == 0) cout << "inf"; else cout << r.size();
      cout << endl;
   }
   unsigned S = SIZE;
   if (r.size() > 0)  S = std::min (S, r.size());

   // No variable

   doTest (r, cat());
   doTest (r, pol());

   // One variable

   for (unsigned i = 0; i < S; i++)
   {
      const T a = r.element (i);

      DEB1 { cout << setw(W); r.printShort (cout, a); }

      // element() and index()

      unsigned ind = r.index(a);
      DEB2 cout << ", index = " << ind;
      fl();
      if (ind != i)  error ("index(a) != i");

      doTest (r, a, cat());
      doTest (r, a, pol());

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
            cout << setw(W); r.printShort (cout, a);
            cout << setw(W); r.printShort (cout, b);
         }

         // Make sure, == makes sense

         if ((a == b) != (i == j))  error ("a==b incorrect");
         if ((a != b) != (i != j))  error ("a!=b incorrect");
         
         doTest (r, a, b, cat());
         doTest (r, a, b, pol());

         DEB1 cout << endl;
      }
   }

   // three variables

   S = unsigned(pow(double(SIZE), 1.0 / 3.0));
   if (r.size() > 0)  S = std::min (S, r.size());

   for (unsigned i = 0; i < S; i++)
   {
      const T a = r.element (i);

      for (unsigned j = 0; j < S; j++)
      {
         const T b = r.element (j);

         for (unsigned k = 0; k < S; ++k)
         {
            const T c = r.element (k);

            DEB1
            {
               cout << setw(W); r.printShort (cout, a);
               cout << setw(W); r.printShort (cout, b);
               cout << setw(W); r.printShort (cout, c);
            }

            doTest (r, a, b, c, cat());

            DEB1 cout << endl;
         }
      }
   }
}


/************  Command line arguments  ***************************************/

const char* options = "n:f";

bool opt (int c, const char* s)
{
   switch (c)
   {
   case 'n':  SIZE = atoi (s); return true;
   case 'f':  FLUSH = true; return true;
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
      "  -f     Flush the output stream at regular times.\n"
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

   {  // Reals and real polynomials

      RealField<real> realField;
      doTests (realField, "RealField<real>");

      PolynomialRing<RealField<> > realPolynomials (realField);
      doTests (realPolynomials, "PolynomialRing<RealField<real> >");

      doTests (QuotientField<PolynomialRing<RealField<real> > >
            (realPolynomials),
              "QuotientField<PolynomialRing<RealField<real> > >");
   }

   {  // Integers and integer polynomials

      IntegerRing<> intRing;
      doTests (intRing, "IntegerRing<int>");
      doTests (PolynomialRing<IntegerRing<int> > (intRing, 'y'),
              "PolynomialRing<IntegerRing<int> >");

      QuotientField<IntegerRing<int> > rationalField (intRing);
      doTests (rationalField, "QuotientField<IntegerRing<int> >");

      doTests (PolynomialRing<QuotientField<IntegerRing<int> > >(rationalField),
              "PolynomialRing<QuotientField<IntegerRing<int> > >");
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

   const unsigned bases [] =
      { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16, 17, 18,
       32,  64, 128,
       27,  81, 243,
       25, 125,
       49,
       121,
       169,
       43, 91};

   for (unsigned ii = 0; ii < sizeof(bases)/sizeof(unsigned); ++ii)
   {
      const unsigned i = bases [ii];

      // Finite field of the given size

      unsigned p;
      unsigned e;
      if (Prime::isPrimePower (i, p, e))
      {
         unsigned maxDim   = digitsRepresentable (static_cast<u8>(i));
         unsigned maxDim32 = logInt (std::numeric_limits<unsigned>::max(), i);

         {
            LookupGaloisField<u8> gf (i);
            doTests (gf, "LookupGaloisField<>");
            PolynomialRing<LookupField<u8> > polys (gf);
            doTests (polys, "PolynomialRing<LookupField<> >");
            doTests (QuotientField<PolynomialRing<LookupField<u8> > > (polys),
                    "QuotientField<PolynomialRing<LookupField<> > >");
            doTests (OneDimVectorSpace<LookupField<u8> > (gf),
                     "OneDimVectorSpace<LookupField<> >");

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
                     "PolynomialRing<LookupFieldPrime<> >");
            doTests (OneDimVectorSpace<LookupFieldPrime<u8> > (gf),
                     "OneDimVectorSpace<LookupFieldPrime<> >");
         }
        
         // optimized version for powers of 2

         if (p == 2)
         {
            LookupGaloisFieldPow2<u8> gf (i);
            doTests (gf, "LookupGaloisFieldPow2<>");
            doTests (PolynomialRing<LookupFieldPow2<u8> > (gf),
                     "PolynomialRing<LookupFieldPow2<> >");
            doTests (OneDimVectorSpace<LookupFieldPow2<u8> > (gf),
                     "OneDimVectorSpace<LookupFieldPow2<> >");

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
         FactorField<unsigned short> field (i);
         doTests (field, "FactorField<>");

         PolynomialRing<FactorField<unsigned short> > poly (field);
         doTests (poly, "PolynomialRing<FactorField<>");
         
         unsigned short coef [] = {0, 0, 1};
         PolynomialRing<FactorField<unsigned short> >::type
            p (coef, coef + sizeof(coef)/sizeof(unsigned short));
         FactorRing<PolynomialRing<
             FactorField<unsigned short> > > mod (poly, p);
         doTests (mod,
            "FactorRing<PolynomialRing<FactorField<> > >");

         for (unsigned j = 1; j <= 4; ++j)
         {
            doTests (GaloisField<u8> (i, j), "GaloisField<>");
         }
      }
      else
      {
         FactorRing<unsigned short> modInt (i);
         doTests (modInt, "FactorRing<>");
         doTests (PolynomialRing<FactorRing<unsigned short> > (modInt),
                  "PolynomialRing<FactorRing<> >");
      }
   }

   FactorRing<unsigned short> modInt (213);
   doTests (modInt, "FactorRing<>");

   FactorField<unsigned short> gf1999 (1999);
   doTests (gf1999, "FactorField<>");
}

