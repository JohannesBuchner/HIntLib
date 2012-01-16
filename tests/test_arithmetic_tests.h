/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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

#include <iomanip>
#include <algorithm>
#include <utility>

#include <HIntLib/hlmath.h>

#ifdef HINTLIB_HAVE_SSTREAM
#  include <sstream>
#else
#  include <HIntLib/fallback_sstream.h>
#endif

#include <HIntLib/exception.h>
#include <HIntLib/prime.h>
#include <HIntLib/array.h>
#include <HIntLib/gcd.h>
#include <HIntLib/output.h>

#include "test_arithmetic.h"

using std::endl;
using std::setw;

using HIntLib::group_tag;
using HIntLib::ringfield_tag;
using HIntLib::ringdomain_tag;
using HIntLib::ring_tag;
using HIntLib::domain_tag;
using HIntLib::ufd_tag;
using HIntLib::euclidean_tag;
using HIntLib::integer_tag;
using HIntLib::field_tag;
using HIntLib::numberfield_tag;
using HIntLib::rational_tag;
using HIntLib::real_tag;
using HIntLib::complex_tag;
using HIntLib::gf_tag;
using HIntLib::cyclic_tag;
using HIntLib::funfield_tag;
using HIntLib::ratfunfield_tag;
using HIntLib::vectorspace_tag;

using HIntLib::finite_tag;
using HIntLib::infinite_tag;

using HIntLib::polynomial_tag;
using HIntLib::nopolynomial_tag;

using HIntLib::noprimedetection_tag;
using HIntLib::primedetection_tag;
using HIntLib::factor_tag;

using HIntLib::Overflow;
using HIntLib::Prime;
using HIntLib::Array;
using HIntLib::FIT_FOR_MUL;
using HIntLib::sqr;
using HIntLib::logInt;
using HIntLib::powInt;
using HIntLib::cube;


/**
 *  fl();
 */

inline void fl ()
{
   DEB2 if (FLUSH)  COUT << std::flush;
}

/**
 *  ensure Equal Type()
 */

template<typename X>
inline
void ensureEqualType (const X*, const X*) {}

/**
 *  Cache
 */

struct CacheBase
{
   // no elements

   unsigned size;
   unsigned numUnits;
   unsigned numNilpotents;
   unsigned characteristic;
   int extensionDegree;
   int dimension;
   unsigned baseSize;

   // one element

   unsigned index;         // reinitialized for two elements
   bool is0;               // reinitialized for two elements
   unsigned additiveOrder;
   unsigned order;
   int degree;
   bool isPrime;           // reinitialized for two elements
   bool isPrimitive;
   bool is1;
   bool isUnit;            // reinitialized for two elements
   bool isNilpotent;
   bool isCanonical;
   bool isPrimitiveElement;
   unsigned norm;          // reinitialized for two elements

   // two elements

   unsigned indexB;
   bool is0B;
   bool isUnitB;
};

template<class A>
struct Cache : public CacheBase
{
   typedef typename A::type T;

   Cache (const A& _a) : arith (_a) {}

   const A& operator()()  { return arith; }
   void print (const T& x) const  { arith.printShort (COUT, x); }

   // no elements

   T zero;
   T one;
   T x;

   // one element

   T neg;
   T dbl;
   T derivative;  // reinitialized for two elements
   T sqr;
   T recip;

   // two elements

   T add;
   T sub;
   T mul;

private:
   const A arith;
};


/*****************************************************************************/
/*************  size_category                                    *************/
/*****************************************************************************/


// no elements

inline
void doTest (CacheBase& c, finite_tag)
{
   if (c.size == 0)  error ("size() == 0 in finite_tag structure");
}

inline
void doTest (CacheBase& c, infinite_tag)
{
   if (c.size != 0)  error ("size() != 0 in infinite_tag structure");
}


/*****************************************************************************/
/*************  Groups                                           *************/
/*****************************************************************************/


template<class G>
void doTest (Cache<G>& c, group_tag)
{
   // doTest (c, XXXX());

   typedef typename G::type T;

   c.zero = T();

   DEB1
   {
      COUT << "  Zero: ";
      c().print (COUT, c.zero);
      COUT << endl;
   }

   if (c().index(c.zero) != 0)  error ("index(zero) invalid");

   if (! c().is0 (c().element(0)))
   {
      error ("element(0) != 0");
   }
}

template<class G>
void doTest (Cache<G>& c, const typename G::type& a, group_tag)
{
   // doTest (c, a, XXXX());

   typedef typename G::type T;

   // Copy and assignment

   {
      T zero;
      T copy (a);
      if (! (copy == a))  error ("Copy constructor or operator==() broken");
      zero = a;
      if (! (zero == a))  error ("Assignment, or operator==() broken");
   }
   
   // is0()

   c.is0 = c().is0(a);
   DEB2  if (c.is0)  COUT << ", zero";
   if (c.is0 != (a == c.zero))  error ("is0() or operator==() fails");
   
   // zero

   if (! (c().add (a, c.zero) == a))  error ("a+0 != a");

   // additive inverse

   c.neg = c().neg(a);
   DEB2
   {
      COUT <<
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 3
         L", \x2212" L"a = ";  // MINUS SIGN
# else
         L", -a = ";
# endif
#else
         Wgl4Ascii (", \xe2\x88\x92" "a = ", ", -a = ");
#endif
      c.print (c.neg);
   }

   checkIndex (c.size, c().index(c.neg), "neg");
   if (! c().is0 (c().add (a, c.neg)))  error ("a + (-a) != 0");

   {
      // negate()

      T x (a);
      c().negate (x);

      DEB4
      {
         COUT << ", negate(a) = ";
         c.print (x);
      }

      if (! (x == c.neg))  error ("negate(a) != -a");
   }
   
   fl();

   // double

   {
      c.dbl = c().dbl(a);
      DEB2
      {
         COUT << ", 2a = ";
         c.print (c.dbl);
      }

      if (! (c.dbl == c().add (a,a)))  error ("dbl(a) invalid");

      T d = a;
      c().times2 (d);
      if (! (d == c.dbl))
      {
         error ("times2() invalid");
         c.print (d);
      }

      fl();
   }

   // check times()

   if (sqr (c.index) < SIZE)
   {
      DEB3
      {
         COUT <<
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 4
            // LEFT-POINTING ANGLE BRACKET
            // RIGHT-POINTING ANGLE BRACKET
            L", \x2329na\x232a = \x2329";
# else
            L", <na> = <";
# endif
#else
            UnicodeAscii (", \xe2\x8c\xa9na\xe2\x8c\xaa = \xe2\x8c\xa9",
                          ", <na> = <");
#endif
      }

      unsigned SS = unsigned (HINTLIB_MN sqrt(double(SIZE)));

      T x = T();
      for (unsigned l = 0; l < SS; ++l)
      {
         T sum = c().times(a,l);
         DEB3 { c.print (sum); COUT << ','; }
         if (! (x == sum))
         {
            error ("a+...+a != times (a,k)");
            COUT << "a+...+a = ";
            c.print (x);
            COUT << ",  times (a," << l << ") = ";
            c.print (sum);
            COUT << ",  diff = ";
            c.print (c().sub (x, sum));
            COUT << endl;
         }
         c().addTo (x, a);
      }

      DEB3
      {
         COUT <<
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 4
            // RIGHT-POINTING ANGLE BRACKET
            L"...\x232a";
# else
            L"...>";
# endif
#else
            UnicodeAscii ("...\xe2\x8c\xaa", "...>");
#endif
      }
   }

   fl();

   // additive order of the element

   c.additiveOrder = c().additiveOrder (a);
   DEB2
   {
      COUT << ", add.order: ";
      printNumberOrInf (c.additiveOrder);
   }

   if (c.additiveOrder == 0)
   {
      if (c.size != 0)  error ("infinite additive order in finite group");
   }
   else
   {
      if (c.size && c.additiveOrder > c.size)
      {
         error ("additiveOrder() larger than size()");
      }
      if (c.size && c.size % c.additiveOrder != 0)
      {
         error ("additiveOrder() does not divide size()");
      }

      if (sqr(c.size) < SIZE)
      {
         T x = T();
         for (unsigned i = 1; i < c.additiveOrder; ++i)
         {
            c().addTo (x, a);
            if (c().is0(x))  error ("additveOrder() is too high");
         }
         c().addTo (x,a);
         if (! c().is0(x))  error ("additiveOrder() is too low");
      }
      else
      {
         if (! c().is0(c().times(a, c.additiveOrder)))
         {
            error ("additiveOrder() wrong");
         }
      }
   }

   fl();
}

inline
void counterTest (CacheBase&, bool, group_tag)  {}

template<class G>
void doTest (Cache<G>& c, const typename G::type& a,
                          const typename G::type& b, group_tag)
{
   // doTest (c, a, b, XXXX());

   typedef typename G::type T;

   // reinitialize one-element cache values

   c.is0B   = c().is0 (b);
   c.indexB = c().index (b);

   if (c.indexB == 0)
   {
      c.is0    = c().is0 (a);
      c.index  = c().index (a);
   }

   // Make sure that operators don't leave domain

   c.sub = c().sub (a,b);
   c.add = c().add (a,b);

   DEB2
   {
      COUT << ", a+b = ";
      c.print (c.add);
      COUT <<
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 3
         L", a\x2212" L"b = ";  // MINUS SIGN
# else
         L", a-b = ";
# endif
#else
         Wgl4Ascii (", a\xe2\x88\x92" "b = ", ", a-b = ");
#endif
      c.print (c.sub);
   }

   checkIndex (c.size, c().index (c.add), "add");
   checkIndex (c.size, c().index (c.sub), "add");
   
   // Commutativity

   if (! (c.add == c().add (b,a)))  error ("a+b != b+a");

   // derived operations

   if (! (c().add (a, c().neg(b)) == c.sub))  error ("sub() fails");

   T x;

   x = a;
   c().addTo (x, b);
   if (! (x == c.add)) error ("addTo() fails");

   x = a;
   c().subFrom (x, b);
   if (! (x == c.sub)) error ("subFrom() fails");

   fl();
}

template<class G>
void doTest
   (Cache<G>& cc, const typename G::type& a,
                  const typename G::type& b,
                  const typename G::type& c, group_tag)
{
   // doTest (cc, a, b, c, XXXX());

   // Associativity of addition

   if (! (cc().add(cc().add(a, b), c) == cc().add(a, cc().add(b,c))))
   {
      error("(a+b) + c != a + (b+c)");
   }

   fl();
}


/*****************************************************************************/
/*************  char_category                                    *************/
/*****************************************************************************/


inline void requireCharacteristic (HIntLib::char_exists) {}

// no element

inline
void doTest (CacheBase&, HIntLib::char_none)
{
   DEB1  COUT << "  Char: none" << endl;
}

template<class R>
inline
void doTest (Cache<R>& c, HIntLib::char_zero)
{
   DEB1
   {
      COUT << "  Char: ";
      printInfinity();
      COUT << endl;
   }

   c.characteristic = c().characteristic();

   if (c.characteristic != 0)
   {
      error ("Finite characteristic() in char_zero structure");
   }
}

template<class R>
inline
void doTest (Cache<R>& c, HIntLib::char_two)
{
   DEB1  COUT << "  Char: 2" << endl;

   c.characteristic = c().characteristic();

   if (c.characteristic != 2)
   {
      error ("Wrong characteristic() in char_two structure");
   }
}

template<class R>
inline
void doTest (Cache<R>& c, HIntLib::char_prime)
{
   c.characteristic = c().characteristic();

   DEB1  COUT << "  Char: prime = " << c.characteristic << endl;

   if (c.characteristic == 0)
   {
      error ("Infinite characteristic() in char_prime structure");
   }

   if (! Prime::test (c.characteristic))
   {
      error ("Finite characteristic not a prime number");
   }
}

// one element

template<class R>
inline
void doTest (Cache<R>&, const typename R::type&, HIntLib::char_none)
{}

template<class R>
inline
void doTest (Cache<R>& c, const typename R::type&, HIntLib::char_exists)
{
   // additive order must match characteristic

   if (! c.is0)
   {
      if (c.additiveOrder != c.characteristic)
      {
         error ("additiveOrder() != characteristic()");
      }
   }
}

// two elements

template<class R>
inline
void doTest (Cache<R>&, const typename R::type&,
                        const typename R::type&, HIntLib::char_any)
{}


/*****************************************************************************/
/*************  zerodivisor_category                             *************/
/*****************************************************************************/

inline void requireNozerodivisor (HIntLib::nozerodivisor_tag) {}

// no element

inline
void doTest (CacheBase&, HIntLib::zerodivisor_tag)
{
   DEB1  COUT << "  Zero divisors: yes" << endl;
}

inline
void doTest (CacheBase&, HIntLib::nozerodivisor_tag)
{
   DEB1  COUT << "  Zero divisors: no" << endl;
}

// one element

template<class R>
inline
void doTest (Cache<R>&, const typename R::type&, HIntLib::zerodivisor_tag) {}

template<class R>
inline
void doTest (Cache<R>& c, const typename R::type& a, HIntLib::nozerodivisor_tag)
{
   // multiplicative order of the element

   typedef typename R::type T;

   if (! c.is0)
   {
      c.order = c().order (a);

      DEB2
      {
         COUT << ", mult.order: ";
         printNumberOrInf (c.order);
      }

      if (c.order == 0)
      {
         if (c.size != 0)
         {
            error ("infinite multiplicative order in finite ring");
         }
      }
      else  // ord finite
      {
         if (c.size && c.order > c.size-1)  error ("order() larger than size()-1");
         if (c.size && (c.size-1) % c.order != 0)
         {
            error ("order() does not divide size()-1");
         }

         if (sqr(c.size) < SIZE)
         {
            T x = c.one;
            for (unsigned i = 1; i < c.order; ++i)
            {
               c().mulBy (x, a);
               if (c().is1(x))  error ("order() is too high");
            }
            c().mulBy (x,a);
            if (! c().is1(x))  error ("order() is too low");
         }
         else
         {
            if (! c().is1 (c().power (a, c.order)))  error ("order() wrong");
         }
      }
   }

   fl();
}

// two elements

template<class R>
inline
void doTest (
   Cache<R>&,
   const typename R::type&, const typename R::type&, HIntLib::zerodivisor_tag)
{}

template<class R>
inline
void doTest (
   Cache<R>& c, const typename R::type&,
                const typename R::type&, HIntLib::nozerodivisor_tag)
{
   // no zero divisors

   if (! c.is0 && ! c.is0B)
   {
      if (c().is0 (c.mul))  error ("a * b = 0");
   }
}


/*****************************************************************************/
/*************  Polynomials                                      *************/
/*****************************************************************************/

// no element

template<class R>
void doTest (Cache<R>&, nopolynomial_tag)
{
   DEB1 COUT << "  Polynomials: no" << endl;
}

template<class R>
void doTest (Cache<R>& c, polynomial_tag)
{
   lastDegree = -2;
   primitivesCounter = 0;

   typedef typename R::coeff_algebra A;
   typedef typename R::type T;
   typedef typename R::coeff_type CT;

   const A& alg = c().getCoeffAlgebra();
   DEB1 COUT << "  Polynomials: yes, coeff algebra: " << alg << endl;

   c.baseSize = alg.size();

   ensureEqualType (static_cast<CT*>(0), static_cast<typename A::type*>(0));
   ensureEqualType (static_cast<CT*>(0),
                    static_cast<typename T::coeff_type*>(0));

   // x()

   c.x = c().x();
   DEB1
   {
      COUT << "  var: ";
      c().print (COUT, c.x);
      COUT << endl;
   }

   if (c.x.degree() != 1 || ! alg.is1(c.x.lc()) || ! alg.is0 (c.x[0]))
   {
      error ("x() is invalid");
   }

   for (int i = 0; i < 5; ++i)
   {
      T xx = c().x (i);

      if (xx.degree() != i || ! alg.is1(xx.lc()))
      {
         error ("x(k) is invalid");
      }
      for (int j = 0; j < i; ++j)
      {
         if (! alg.is0 (xx[j]))  error ("x(k) is invalid");
      }
      
      T xxdiv1 = xx;
      xxdiv1.divByX (i);
      T xxdiv2 = xx;
      xxdiv2.divByX (i+1);

      if (! c().is1 (xxdiv1) || ! c().is0 (xxdiv2))  error ("divByX() broken");

      T one = c.one;
      if (! (xx == one.mulByX (i)))  error ("mulByX() broken");
   }
}

// one element

template<class R>
inline
void doTestPoly (Cache<R>&, const typename R::type&, ringfield_tag) {}

template<class R>
void doTestPoly (Cache<R>& c, const typename R::type& a, field_tag)
{
   doTestPoly (c, a, ringfield_tag());

   typedef typename R::type T;
   typedef typename R::unit_type UT;
   typedef typename R::Factorization F;

   // squarefreeFactor()

   if (! c.is0)
   {
      F f;
      UT unit = c().squarefreeFactor (f, a);

      DEB2
      {
         COUT << ", s.f.f: ";
         c.print (c().fromUnit (unit));
         for (typename F::const_iterator i = f.begin(); i != f.end(); ++i)
         {
            COUT <<
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 4
               L"\x22c5";  // DOT OPERATOR
# elif HINTLIB_CHARACTER_SET >= 2
               L"\x00b7";  // MIDDLE DOT
# else
               L" * ";
# endif
#else
# if HINTLIB_CHARACTER_SET >= 4 && defined(HINTLIB_UTF8_SELECT)
               (utf8 ? "\xe2\x8b\x85" : "\xb7");  // DOT OPERATOR / MIDDLE DOT
# elif HINTLIB_CHARACTER_SET >= 2 && defined(HINTLIB_UTF8_SELECT)
               (utf8 ? "\xc2\xb7" : "\xb7");      // MIDDLE DOT
# else
               " * ";
# endif
#endif
            c().printShort (COUT, i->first, FIT_FOR_MUL);
            if (i->second > 1)
            {
               HIntLib::Private::PrinterSelector<CHAR>::printer ss (COUT);
               ss.power (i->second);
            }
         } 
      }

      T prod = c().fromUnit (unit);

      for (typename F::const_iterator i = f.begin(); i != f.end(); ++i)
      {
         if (! c().isSquarefree (i->first)) error ("factor not square free");
         if (! c().isCanonical(i->first))   error ("factor not canonical");
         if (c().is0(i->first) || c().isUnit (i->first))
         {
            error ("factor 0 or unit");
         }

         for (typename F::const_iterator j = f.begin(); j != f.end(); ++j)
         {
            if (i == j)  continue;
            if (! genIsCoprime (c(), i->first, j->first))
            {
               error ("factors are not coprime");
            }
         }

         c().mulBy (prod, c().power (i->first, i->second));
      } 

      if (! (prod == a))
      {
         COUT << "\nDifference: ";
         c().print (COUT, c().sub (prod, a));
         
         error ("Product of factors wrong");
      }
   }

   fl();
}

template<class R>
void doTestPoly (Cache<R>& c, const typename R::type& a, gf_tag)
{
   doTestPoly (c, a, field_tag());

   if (c.degree != lastDegree)
   {
      if (c.degree > lastDegree)
      {
         const unsigned order = powInt (c.baseSize, lastDegree) - 1;

         const unsigned expected =
            lastDegree <= 0 ? 0 : Prime::eulerPhi (order) / lastDegree;

         if (primitivesCounter != expected)
         {
            error ("Number of primitive polynomials is wrong");
            COUT << "\nPrimitives found: " << primitivesCounter
                 << "  Expected: " << expected << endl;
         }
      }

      primitivesCounter = 0;
      lastDegree = c.degree;
   }

   fl();

   // isPrimitive

   c.isPrimitive = c().isPrimitive (a);

   if (c.isPrimitive)
   {
      DEB2 COUT << ", primitive";

      if (! c().isPrime(a))  error ("primitive without prime");

      ++primitivesCounter;
   }

   fl();
}

template<class R>
inline
void doTest (Cache<R>&, const typename R::type&, nopolynomial_tag) {}

template<class R>
void doTest (Cache<R>& c, const typename R::type& a, polynomial_tag)
{
   typedef typename R::coeff_algebra A;
   typedef typename R::type T;
   typedef typename R::coeff_type CT;

   const A& alg = c().getCoeffAlgebra();

   // degree 

   c.degree = a.degree();
   DEB2  COUT << ", degree = " << c.degree;
   fl();

   // isConstant

   bool constant = a.isConstant();

   DEB2 if (constant)  COUT << ", constant";

   if (constant != (c.degree <= 0))  error ("isConstant() is wrong");

   if (c.degree == -1)
   {
      if (! c.is0)  error ("Polynomial with deg=-1 is not 0");
   }
   else
   {
      // isMonic()

      bool monic = c().isMonic (a);
      DEB2  if (monic)  COUT << ", monic";

      // lc()

      CT lc = a.lc();
      DEB2
      {
         COUT << ", lc = ";
         alg.printShort (COUT, lc);
      }

      if (monic != alg.is1(lc))  error ("monic() is wrong");
      if (! (lc == a[a.degree()]))  error ("lc() != leading coefficient");
      if (alg.is0(lc))  error ("lc() is zero");

      // ct()

      CT ct = a.ct();
      DEB2
      {
         COUT << ", ct = ";
         alg.printShort (COUT, ct);
      }

      if (! (ct == a[0]))  error ("ct() != constant term");

      fl();
   }

   // Access

   Array<CT> coeff (c.degree + 1);
   a.toCoeff (&coeff[0]);

   DEB2 COUT << ", coeff = [";
   for (int i = 0; i <= c.degree; ++i)
   {
      DEB2
      {
         alg.printShort (COUT, a[i]);
         if (i < c.degree)  COUT << ';';
      }
      
      if (! (coeff[i] == a[i]))  error ("toCoeff() or operator[]() broken");
   }
   DEB2 COUT << ']';

   fl();

   // Constructor with iterators

   T b (&coeff[0], &coeff[c.degree + 1]);

   if (! (a == b))  error ("Constructor with iterators broken");

   // mulByX(), divByX()

   T div (a);
   T mul (a);
   
   for (int i = 0; i < 5; ++i)
   {
      T div1 (a);
      T mul1 (a);

      div1.divByX (i);
      mul1.mulByX (i);

      if (! (div == div1))  error ("divByX() broken");
      if (! (mul == mul1))  error ("mulByX() broken");

      div.divByX();
      mul.mulByX();
   }

   // derivative()

   c.derivative = c().derivative (a);

   DEB2
   {
      COUT << ", deriv. = ";
      c.print (c.derivative);
   }
   
   T deriv1 = c().derivative (c.neg);
   T deriv2 = c().neg (c.derivative);

   if (! (deriv1 == deriv2))
   {
      error ("derivative() inconsistent with negation");
   }

   fl();

   // evaluate()

   if (sqr (c.index) <= SIZE)
   {
      DEB3  COUT << ", evaluation: {";

      unsigned ub = unsigned (HINTLIB_MN sqrt(double(SIZE)));
      if (c.baseSize && ub > c.baseSize)  ub = c.baseSize;

      for (unsigned i = 0; i < ub; ++i)
      {
         CT x = alg.element (i);
         CT eval = c().evaluate (a, x);

         DEB3
         {
            COUT << "a(";
            alg.printShort (COUT, x);
            COUT << ")=";
            alg.printShort (COUT, eval);

            if (i + 1 < ub)  COUT << ", ";
         }

         CT res = CT();
         for (int j = c.degree; j >= 0; --j)
         {
            res = alg.add (alg.mul (res, x), a[j]);
         }
      
         if (! (eval == res))  error ("evaluate() broken");
      }

      DEB3 
      {
         if (c.baseSize == 0 || c.baseSize > ub)  COUT << ",...";
         COUT << '}';
      }
   }

   doTestPoly (c, a, typename A::algebra_category());
}

// two elements

template<class R>
inline
void doTest (Cache<R>&, const typename R::type&,
                        const typename R::type&, nopolynomial_tag) {}

template<class R>
void doTest (Cache<R>& c, const typename R::type& a,
                          const typename R::type& b, polynomial_tag)
{
   typedef typename R::coeff_algebra A;
   typedef typename R::type T;
   typedef typename R::coeff_type CT;

   // reinitialize one-element cache values

   if (c.indexB == 0)
   {
      c.derivative = c().derivative (a);
   }

   // derivative()

   const T derivB = c().derivative (b);

   if (! (c().derivative (c.add) == c().add (c.derivative, derivB)))
   {
      error ("derivative() inconsistent with addition");
   }

   if (! (   c().derivative (c.mul)
          == c().add (c().mul (a, derivB), c().mul (b, c.derivative))))
   {
      error ("derivative() inconsistent with multiplication");
   }
}


/*****************************************************************************/
/*************  Ring / Field                                     *************/
/*****************************************************************************/

template<class R>
void doTest (Cache<R>& c, ringfield_tag)
{
   doTest (c, group_tag());

   // Check one()

   c.one = c().one();

   DEB1
   {
      COUT << "  One:  ";
      c().print (COUT, c.one);
      COUT << endl;
   }

   if (c().index(c.one) != 1)  error ("index(one) != 1");

   if (! c().is1(c().element(1)))  error ("element(1) invalid");

   // call characteristic and zerodivisor tests

   doTest (c, typename R::polynomial_category());
   doTest (c, typename R::char_category());
   doTest (c, typename R::zerodivisor_category());
}

template<class R>
void doTest (Cache<R>& c, const typename R::type& a, ringfield_tag)
{
   doTest (c, a, group_tag());

   typedef typename R::type T;

   // is1(), one()

   c.is1 = c().is1(a);
   DEB2 if (c.is1)  COUT << ", one";
   if (c.is1 != (a == c.one ))  error ("is1() fails");
      
   // neutral

   if (! (c().mul (a, c.one) == a))  error ("a*1 != a");

   // square

   {
      try
      {
         c.sqr = c().sqr(a);
         DEB2
         {
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 2
            COUT << L", a\x00b2 = ";  // SUPERSCRIPT TWO
# else
            COUT << L", a^2 = ";
# endif
#else
            COUT << Latin1Ascii(", a\xc2\xb2 = ", ", a\xb2 = ", ", a^2 = ");
#endif
            c.print (c.sqr);
         }

         if (! (c.sqr == c().mul (a,a)))  error ("sqr(a) invalid");

         T s = a;
         c().square (s);
         if (! (s == c.sqr))  error ("square() invalid");

         fl();
      }
      catch (Overflow &) {}
   }

   // check power()

   if (sqr(c.index) < SIZE)
   {
      DEB3
      {
         COUT <<
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 4
            // LEFT-POINTING ANGLE BRACKET
            // SUPERSCRIPT LATIN SMALL LETTER N
            // RIGHT-POINTING ANGLE BRACKET
            L", \x2329" L"a\x207f\x232a = \x2329";
# elif HINTLIB_CHARACTER_SET >= 3
            // SUPERSCRIPT LATIN SMALL LETTER N
            L", <a\x207f> = <";
# else
            L", <na> = <";
# endif
#else
            UnicodeAscii (", \xe2\x8c\xa9", ", <")
         << Wgl4Ascii ("a\xe2\x81\xbf", "a^n")
         << UnicodeAscii ("\xe2\x8c\xaa = \xe2\x8c\xa9", "> = <");
#endif
      }

      unsigned SS = unsigned (HINTLIB_MN sqrt(double (SIZE)));

      if (SS > 10) SS = 7;

      T x = a;
      for (unsigned l = 1; l < SS; ++l)
      {
         try
         {
            T p = c().power (a, l);
         
            DEB3
            {
               c.print (p);
               DEB4  { COUT << '='; c.print (x); }
               COUT << ',';
            }
            if (! (x == p))
            {
               error ("a*...*a != power (a,k)");
               COUT << "a*...*a = ";
               c().print (COUT, x);
               COUT << ",  power (a," << l << ") = ";
               c().print (COUT, p);
               COUT << endl;
            }
            c().mulBy (x, a);
         }
         catch (Overflow &)
         {
            break;
         }
      }

      DEB3
      {
         COUT <<
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 4
            // RIGHT-POINTING ANGLE BRACKET
            L"...\x232a";
# else
            L"...>";
# endif
#else
            UnicodeAscii ("...\xe2\x8c\xaa", "...>");
#endif
      }
   }

   fl();

   doTest (c, a, typename R::polynomial_category());
   doTest (c, a, typename R::char_category());
   doTest (c, a, typename R::zerodivisor_category());
}

template<class R>
void doTest
   (Cache<R>& c, const typename R::type& a,
                 const typename R::type& b, ringfield_tag)
{
   doTest (c, a, b, group_tag());

   // mul()

   typedef typename R::type T;

   c.mul = c().mul (a,b);

   DEB2
   {
      COUT <<
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 4
         L", a\x22C5" L"b = ";  // DOT OPERATOR
# elif HINTLIB_CHARACTER_SET >= 2
         L", a\x00b7" L"b = ";  // MIDDLE DOT
# else
         L", a*b = ";
# endif
#else
# if HINTLIB_CHARACTER_SET >= 4  && defined(HINTLIB_UTF8_SELECT)
         (utf8 ? ", a\xe2\x8B\x85" "b = " : ", a\xb7" "b = "); // both
# elif HINTLIB_CHARACTER_SET >= 2  && defined(HINTLIB_UTF8_SELECT)
         (utf8 ? ", a\xc2\xb7" "b = " : ", a\xb7" "b = ");  // MIDDLE DOT
# else
         ", a*b = ";
# endif
#endif
      c.print (c.mul);
   }

   checkIndex (c.size, c().index (c.mul), "mul");

   // Commutativity of multiplication

   if (! (c.mul == c().mul (b,a)))  error ("a*b != b*a");

   T x = a;
   c().mulBy (x, b);
   if (! (x == c.mul))  error ("mulBy() fails");

   fl();

   doTest (c, a, b, typename R::polynomial_category());
   doTest (c, a, b, typename R::zerodivisor_category());
   doTest (c, a, b, typename R::char_category());
}

template<class R>
void doTest
   (Cache<R>& cc, const typename R::type& a,
                  const typename R::type& b,
                  const typename R::type& c, ringfield_tag)
{
   doTest (cc, a, b, c, group_tag());

   // Associativity of multiplication

   if (! (cc().mul(cc().mul(a, b), c) == cc().mul(a, cc().mul(b, c))))
   {
      error("(a*b) * c != a * (b*c)");
   }

   // Distributivity

   if (! (   cc().mul (a, cc().add(b, c))
          == cc().add (cc().mul(a, b), cc().mul(a,c))))
   {
      error("a*(b+c) != (a*b)+(a*c)");
   }

   if (! (   cc().mul (cc().add(a, b), c)
          == cc().add (cc().mul(a,c), cc().mul(b,c))))
   {
      error("(a+b)*c != (a*c)+(b*c)");
   }

   fl();
}


/*****************************************************************************/
/*************  Rings / Domain                                   *************/
/*****************************************************************************/


template<class R>
void doTest (Cache<R>& c, ringdomain_tag)
{
   doTest (c, ringfield_tag());

   typedef typename R::type T;
   typedef typename R::unit_type UT;

   // make sure, unit_type works

   {
      UT x = c().toUnit(c.one);
      UT def = x;
      if (! (def == x))  error ("unit_type broken");
      def = x;
      if (! (def == x))  error ("unit_type broken");
   }

   // isUnit(one), isUnit(zero)

   if (! c().isUnit (c.one))  error ("1 must be a unit");
   if (  c().isUnit (c.zero)) error ("0 must not be a unit");

   // initialize counter for units

   unitsCounter = 0;

   // numUnits()

   c.numUnits = c().numUnits();

   DEB1
   {
      COUT << "  # units: ";
      printNumberOrInf (c.numUnits);
      COUT << endl;
   }

   checkInfiniteInFinite (c.size, c.numUnits, "units");
}

template<class R>
void doTest (Cache<R>& c, const typename R::type& a, ringdomain_tag)
{
   doTest (c, a, ringfield_tag());

   typedef typename R::type T;
   typedef typename R::unit_type UT;

   // check properties of units

   c.isUnit = c().isUnit (a);

   DEB2 if (c.isUnit) COUT << ", unit";
   fl();

   if (c.is0 + c.isUnit > 1)  error ("Classifiaction ambiguous");

   if (c.isUnit)
   {
      ++unitsCounter;

      UT unit = c().toUnit (a);

      if (! (c().fromUnit(unit) == a))
      {
         error ("toUnit() and fromUnit() do not match");
      }

      UT recip = c().unitRecip (unit);
      T recipElem = c().fromUnit (recip);

      DEB2
      {
         COUT <<
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 3
            L", 1\x2215" L"a = ";  // DIVISION SLASH
# else
            L", 1/a = ";
# endif
#else
            Wgl4Ascii(", 1\xe2\x88\x95" "a = ", ", 1/a = ");
#endif
         c.print (recipElem);
      }

      checkIndex (c.size, c().index(recipElem), "unitRecip");

      if (! c().is1 (c().mul (a, recipElem)))
      {
         error ("Reciprocal of unit is wrong");
      }
      if (! c().is1 (c().mulUnit (a, recip)))
      {
         error ("Reciprocal of unit is wrong");
      }
      if (! c().is1 (c().fromUnit (c().mulUnit (unit, recip))))
      {
         error ("Reciprocal of unit is wrong");
      }
   }

   fl();
}

inline
void counterTest (CacheBase& c, bool all, ringdomain_tag)
{
   counterTest (c, all, ringfield_tag());

   // numUnits()

   checkCounter (c.numUnits, unitsCounter, all, "numUnits()");
}

template<class R>
void doTest (Cache<R>& c, const typename R::type& a,
                          const typename R::type& b, ringdomain_tag)
{
   doTest (c, a, b, ringfield_tag());

   typedef typename R::type T;
   typedef typename R::unit_type UT;

   // reinitialize one-element cache values

   if (c.indexB == 0)
   {
      c.isUnit = c().isUnit (a);
   }

   // product of units/non-units

   c.isUnitB = c().isUnit (b);
   const bool pIsUnit = c().isUnit (c.mul);

   if (c.isUnit && c.isUnitB)
   {
      if (! pIsUnit)  error ("Product of units must be a unit");
   }
   else
   {
      if (pIsUnit)  error ("Product with non-unit cannot be a unit");
   }

   // mulUnit(), mulByUnit()

   if (c.isUnit)
   {
      UT aUnit = c().toUnit (a);

      if (! (c().mulUnit (b, aUnit) == c.mul))
      {
         error("mulUnit(type,unit_type) broken");
      }

      T x (b);
      c().mulByUnit (x, aUnit);
      if (! (x == c.mul))  error ("mulByUnit(type&, unit_type) broken");

      if (c.isUnitB)
      {
         UT bUnit = c().toUnit (b);
         if (! (c().fromUnit(c().mulUnit (aUnit, bUnit)) == c.mul))
         {
            error ("mulUnit(unit_type, unit_type) broken");
         }

         UT xUnit (aUnit);
         c().mulByUnit (xUnit, bUnit);

         if (! (c().fromUnit(xUnit) == c.mul))
         {
            error ("mulByUnit(unit_type&, unit_type) broken");
         }
      }
   }
}


/*****************************************************************************/
/*************  Rings                                            *************/
/*****************************************************************************/

template<class R>
void doTest (Cache<R>& c, ring_tag)
{
   doTest (c, ringdomain_tag());

   typedef typename R::type T;

   // numNilpotent()

   c.numNilpotents = c().numNilpotents();
   nilpotentsCounter = 0;

   DEB1
   {
      COUT << "  # nilpotents: ";
      printNumberOrInf (c.numNilpotents);
      COUT << endl;
   }

   checkInfiniteInFinite (c.size, c.numNilpotents, "nilpotents");

   if (c.size && c.size % c.numNilpotents != 0)
   {
      error ("numNilpotents() does not divide size()");
   }
}

template<class R>
void doTest (Cache<R>& c, const typename R::type& a, ring_tag)
{
   doTest (c, a, ringdomain_tag());

   typedef typename R::type T;

   // isNilpotent()

   c.isNilpotent = c().isNilpotent (a);

   DEB2 if (c.isNilpotent) COUT << ", nilpotent";

   if (c.isUnit && c.isNilpotent)  error ("unit cannot be nilpotent");
   if (c.is0 && ! c.isNilpotent)  error ("0 must be nilpotent");

   if (c.isNilpotent)
   {
      ++nilpotentsCounter;

      T p (a);
      unsigned i = 0;

      while (! c().is0(p))
      {
         c().square (p);

         if (++i > 20)
         {
            error ("Element does not seem to be nilpotent");
            break;
         }
      }
   }
   else
   {
      T p (a);

      for (unsigned i = 0; i < 5; ++i)
      {
         c().square (p);

         if (c().is0 (p))
         {
            error ("Element is nilpotent");
            break;
         }
      }
   } 
}

inline
void counterTest (CacheBase& c, bool all, ring_tag)
{
   counterTest (c, all, ringdomain_tag());

   // numNilpotents()

   checkCounter (c.numNilpotents, nilpotentsCounter, all, "numNilpotents()");
}


/*****************************************************************************/
/*************  Prime Detection                                  *************/
/*****************************************************************************/

// no element

template<class R>
void doTest (Cache<R>&, noprimedetection_tag)
{
   DEB1  COUT << "  Prime detection: no" << endl;
}

template<class R>
void doTest (Cache<R>& c, primedetection_tag)
{
   DEB1  COUT << "  Prime detection: yes" << endl;

   // Prime generator

   typedef typename R::type T;
   typename R::PrimeGenerator pg (c());

   DEB3 COUT << "Primes: ";
   bool first = true;

   for (unsigned i = 0; i < SIZE; ++i)
   {
      T p = pg.next();

      DEB3
      {
         if (! first)  COUT << ", ";
         c.print (p);
         fl();
         first = false;
      }

      if (! c().isCanonical (p))
      {
         return error ("PrimeGenerator::next() does not return canonical");
      }

      if (! c().isPrime (p))
      {
         return error ("PrimeGenerator::next() does not return prime");
      }
   }

   DEB3 COUT << ",..." << endl;
}

template<class R>
void doTest (Cache<R>& c, factor_tag)
{
   doTest (c, primedetection_tag());
   DEB1  COUT << "  Factorization possible" << endl;
}

// one element

template<class R>
void doTest (Cache<R>&, const typename R::type&, noprimedetection_tag)
{}

template<class R>
void doTest (Cache<R>& c, const typename R::type& a, primedetection_tag)
{
   c.isPrime = c().isPrime (a);
   bool isComp  = c().isComposit (a);

   DEB2
   {
      if (c.isPrime) COUT << ", prime";
      if (  isComp)  COUT << ", composit";
   }

   if (c.is0 + c.isUnit + c.isPrime + isComp != 1)
   {
      error ("Classifiaction ambiguous");
   }

   fl();
}

template<class R>
void doTest (Cache<R>& c, const typename R::type& a, factor_tag)
{
   typedef typename R::type T;
   typedef typename R::unit_type UT;
   typedef typename R::Factorization F;

   doTest (c, a, primedetection_tag());

   // factor()

   if (! c.is0)
   {
      F f;
      UT unit = c().factor (f, a);

      DEB2
      {
         COUT << ", fact.: ";
         c.print (c().fromUnit (unit));
         for (typename F::const_iterator i = f.begin(); i != f.end(); ++i)
         {
            COUT << 
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 4
               L"\x22c5";  // DOT OPERATOR
# elif HINTLIB_CHARACTER_SET >= 2
               L"\x00b7";  // MIDDLE DOT
# else
               L" * ";
# endif
#else
# if HINTLIB_CHARACTER_SET >= 4 && defined(HINTLIB_UTF8_SELECT)
               (utf8 ? "\xe2\x8b\x85" : "\xb7");  // DOT OPERATOR / MIDDLE DOT
# elif HINTLIB_CHARACTER_SET >= 2 && defined(HINTLIB_UTF8_SELECT)
               (utf8 ? "\xc2\xb7" : "\xb7");      // MIDDLE DOT
# else
               " * ";
# endif
#endif
            c().printShort (COUT, i->first, FIT_FOR_MUL);
            
            if (i->second > 1)
            {
               HIntLib::Private::PrinterSelector<CHAR>::printer ss (COUT);
               ss.power (i->second);
            }
         } 
      }

      T prod = c().fromUnit (unit);

      for (typename F::const_iterator i = f.begin(); i != f.end(); ++i)
      {
         if (! c().isPrime (i->first)) error ("factor not prime");
         if (! c().isCanonical(i->first))   error ("factor not canonical");

         for (typename F::const_iterator j = f.begin(); j != f.end(); ++j)
         {
            if (i == j)  continue;
            if (! genIsCoprime (c(), i->first, j->first))
            {
               error ("factors are not coprime");
            }
         }

         c().mulBy (prod, c().power (i->first, i->second));
      } 

      if (! (prod == a))
      {
         COUT << "\nDifference: ";
         c().print (COUT, c().sub (prod, a));
         error ("Product of factors wrong");
      }
   }

   fl();
}

// two elements

template<class R>
void doTest
  (Cache<R>&, const typename R::type&, const typename R::type&,
   noprimedetection_tag)
{}

template<class R>
void doTest
  (Cache<R>& c, const typename R::type& a, const typename R::type& b,
   primedetection_tag)
{
   // reinitialize one-element cache values

   if (c.indexB == 0)
   {
      c.isPrime = c().isPrime (a);
   }

   if (c.is0 || c.is0B)  return;

   if (! c.isUnit || ! c.isUnitB)  // two units has already been checked
   {
      bool bIsPrime = c().isPrime (b);

      if ((c.isUnit && bIsPrime) || (c.isPrime && c.isUnitB))
      {
         if (! c().isPrime (c.mul))
            error ("Product of unit and prime must be prime");
      }
      else
      {
         if (! c().isComposit (c.mul))
         {
            error ("product of two non-units must be composit");
         }
      }
   }
}


/*****************************************************************************/
/*************  Domains                                          *************/
/*****************************************************************************/


template<class R>
void doTest (Cache<R>& c, domain_tag)
{
   doTest (c, ringdomain_tag());
   doTest (c, typename R::primedetection_category());

   requireCharacteristic (typename R::char_category());
   requireNozerodivisor  (typename R::zerodivisor_category());
}

template<class R>
void doTest (Cache<R>& c, const typename R::type& a, domain_tag)
{
   doTest (c, a, ringdomain_tag());
   doTest (c, a, typename R::primedetection_category());

   if (! c().isAssociate (a, a))
   {
      error ("isAssociate() must be reflexive");
   }

   if (! c.is0)
   {
      if (! c().isDivisor (a, a))
      {
         error ("isDivisor() must be reflexive");
      }
   }
}

template<class R>
void doTest (Cache<R>& c, const typename R::type& a,
                          const typename R::type& b, domain_tag)
{
   doTest (c, a, b, ringdomain_tag());
   doTest (c, a, b, typename R::primedetection_category());

   // isAssociate

   bool isAssociate = c().isAssociate (a, b);

   DEB2 if (isAssociate)
   {
      COUT <<
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 4
         L", a\x223C" "b";  // TILDE OPERATOR
# else
         L", a~b";
# endif
#else
         UnicodeAscii (", a\xe2\x88\xbc" "b", ", a~b");
#endif
   }

   if (isAssociate != c().isAssociate (b,a))
   {
      error ("isAssociate() must be symmetric");
   }

   // isDivisor

   typedef typename R::type T;

   if (! c.is0B)
   {
      bool isDivisor = c().isDivisor (a,b);
      bool isUnitB   = c().isUnit (b);

      DEB2 if (isDivisor)
      {
         COUT <<
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 4
            L", a\x2223" L"b";  // DIVIDES
# else
            L", a|b";
# endif
#else
            UnicodeAscii (", a\xe2\x88\xa3" "b", ", a|b");
#endif
      }

      if (! isAssociate && isDivisor && ! c.is0 && c().isDivisor (b,a))
      {
         error ("isDivisor() must be anti-symmetric");
      }

      T div = T();

      if (c().isDivisor (a,b,div) != isDivisor)
      {
         error ("isDivisor()-2 and isDivisor()-3 inconsistent");
      }

      if (isDivisor)
      {
         DEB2
         {
            COUT <<
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 3
               L", a\x2215" L"b = ";  // DIVISION SLASH
# else
               L", a/b = ";
# endif
#else
               Wgl4Ascii (", a\xe2\x88\x95" "b = ", ", a/b = ");
#endif
            c.print (div);
         }

         if (! (c().mul (b,div) == a))  error ("exact division broken");

         if (! (c().div (a,b) == div))  error ("div() broken");

         T div2 (a);
         c().divBy (div2,b);
         if (! (div2 == div))  error ("divBy() broken");

         if (isAssociate && ! c().isUnit (div))
         {
            error ("divisor for associates must be a unit");
         }
      }

      if (isUnitB && ! isDivisor)  error ("unit must be a divisor");
      if (c.isUnit && isUnitB && ! c().isUnit (div))
      {
         error ("unit/unit must be unit");
      }
      if (c.isUnit && ! isUnitB && isDivisor)
      {
         error ("unit has non-unit divisor");
      }

      fl();

      // isAssociate

      if (isAssociate)
      {
         if (! isDivisor)
         {
            error ("associates must divide each other");
         }
      }
   }
   else  // b = 0
   {
      if (c.is0 != isAssociate)  error ("isAssociate() for zero wrong");
   }
}

template<class R>
void doTest (Cache<R>& cc, const typename R::type& a,
                           const typename R::type& b,
                           const typename R::type& c, domain_tag)
{
   doTest (cc, a, b, c, ringdomain_tag());

   // Transitivity of isAssociate()

   if (   cc().isAssociate(a,b) && cc().isAssociate(b,c)
       && ! cc().isAssociate(a,c))
   {
      error ("isAssociate() must be transitive");
   }

   // Transitivity of isDivisor()

   if (! cc().is0(b) && ! cc().is0(c))
   {
      if (cc().isDivisor(a,b) && cc().isDivisor(b,c) && ! cc().isDivisor(a,c))
      {
         error ("isDivisor() must be transitive");
      }
   }
}


/*****************************************************************************/
/*************  UFDs                                             *************/
/*****************************************************************************/

template<class R>
void doTest (Cache<R>& c, ufd_tag)
{
   typedef typename R::type T;
   typedef typename R::unit_type UT;

   doTest (c, domain_tag());

   // enumerate units

   DEB3 COUT << "Units: ";
   bool first = true;

   unsigned ub = c.numUnits ? std::min (SIZE, c.numUnits) : SIZE;
   for (unsigned i = 0; i < ub; ++i)
   {
      UT unit = c().unitElement (i);
      T  unitElem = c().fromUnit (unit);

      DEB3
      {
         if (! first)  COUT << ", ";
         c.print (unitElem);
         first = false;
      }

      if (! c().isUnit (unitElem))
      {
         error ("result of unitElement() not isUnit()");
      }
      
      if (c().unitIndex (unit) != i)
      {
         error ("unitElement() and unitIndex() do not match");
      }

      if (c.size)
      {
         if (! (c().element(i + 1) == unitElem))
         {
            error ("element() and unitElement() do not match");
         }
      }
   }

   DEB3
   {
      if (c.numUnits == 0 || ub < c.numUnits)
      {
         COUT << ",...";
      }
      COUT << endl;
   }
}

template<class R>
void doTest (Cache<R>& c, const typename R::type& a, ufd_tag)
{
   typedef typename R::type T;
   typedef typename R::unit_type UT;

   doTest (c, a, domain_tag());

   c.isCanonical = c().isCanonical (a);

   DEB2  if (c.isCanonical) COUT << ", canonical";

   if (c.is0 && ! c.isCanonical)  error ("0 should be canonical");

   // make canonical

   {
      T copy = a;
      UT unit = c().makeCanonical (copy);
      T unitElem = c().fromUnit (unit);

      if (! (c().mulUnit (copy, unit) == a && c().mul (copy, unitElem) == a))
      {
         error ("makeCanonical() does not factor");
      }

      if (! c().isUnit(unitElem))
      {
         error ("makeCanonical() does not return a unit");
      }
      if (! c().isCanonical (copy))
      {
         error ("makeCanonical() does not create a canonical form");
      }

      // Canonical form of units

      if (c.isUnit)
      {
         if (! (a == unitElem && c().is1 (copy)))
         {
            error ("makeCanonical() of unit is broken");
         }
      }

      // Canonical form of (non-)canonical forms

      if (c.isCanonical)
      {
         if (! (c().is1(unitElem) && copy == a))
         {
            error ("makeCanonical() of canonical form invalid");
         }
      }
      else // ! isCanonical
      {
         if (c().is1(unitElem) || copy == a)
         {
            error ("makeCanonical() for non-canonical form invalid");
         }
      }
   }

   fl();

   // advanced unit properties

   if (c.isUnit)
   {
      UT unit = c().toUnit (a);

      unsigned index = c().unitIndex (unit);
      DEB2 COUT << ", unit index = " << index;
      checkIndex (c.numUnits, index, "unitIndex");
   }
}

template<class R>
void doTest
  (Cache<R>& c, const typename R::type& a, const typename R::type& b, ufd_tag)
{
   doTest (c, a, b, domain_tag());

   typedef typename R::type T;

   // reinitialize one-element cache values

   if (c.indexB == 0)
   {
      c.isCanonical = c().isCanonical (a);
   }

   // isCononical()
   
   const bool isCanonicalB = c().isCanonical (b);

   if (c.isCanonical && isCanonicalB)
   {
      if (! c().isCanonical (c.mul))
      {
         error ("Product of canonicals must be canonical");
      }
   }
}


/*****************************************************************************/
/*************  Euclidean Ring                                   *************/
/*****************************************************************************/


template<class R>
void doTest (Cache<R>& c, const typename R::type& a, euclidean_tag)
{
   typedef typename R::type T;

   doTest (c, a, ufd_tag());

   c.norm = c().norm (a);

   DEB2
   {
      COUT <<
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 4
         L", \x2016" L"a\x2016 = "  // DOUBLE VERTICAL LINE
# else
         L", |a| = "
# endif
#else
         UnicodeAscii (", \xe2\x80\x96" "a\xe2\x80\x96 = ", ", |a| = ")
#endif
         << c.norm;
   }

   if (c.is0 != (c.norm == 0))  error ("norm(a)==0 iff a==0 failed");
   if (c.isUnit != (c.norm == 1))  error ("norm(a)==1 iff isUnit(a) failed");
   if (c.numUnits && ((c.norm) < lastNorm)) error ("norm() decreased");
   lastNorm = c.norm;

   if (! c.is0)
   {
      unsigned num = c().numOfRemainders (a);

      DEB2
      {
         COUT << ", numRem: ";
         printNumberOrInf (num);
      }

      checkInfiniteInFinite (c.size, num, "remainders");
   }
}

template<class R>
void doTest (Cache<R>& c, const typename R::type& a,
                          const typename R::type& b, euclidean_tag)
{
   typedef typename R::type T;

   doTest (c, a, b, ufd_tag());

   lastNorm = 0;

   // reinitialize one-element cache values

   if (c.indexB == 0)
   {
      c.norm = c().norm (a);
   }

   // norm

   unsigned normB  = c().norm(b);
   unsigned normAB = c().norm(c.mul);

   if (c.norm * normB < normAB)  error ("norm(a)norm(b) < norm(ab)");

   // div(), quot(), and rem()

   if (! c().is0 (b))
   {
      if (c.norm > normAB)  error ("norm(ab) < norm(a)");

      // div()

      T q, re;
      c().div (a, b, q, re);

      // Make sure, operators don't leave domain

      checkIndex (c.size, c().index (q),  "quot");
      checkIndex (c.size, c().index (re), "rem");

      DEB2
      {
         COUT << ", quot. = ";
         c.print (q);
         COUT << ", rem. = ";
         c.print (re);
      }

      if (! (c().add (c().mul (b, q), re) == a))  error ("div() fails");

      // rem() and quot()

      if (! (c().rem(a,b) == re))  error ("rem() inconsistent with div()");

      T x = a;
      c().reduce (x, b);
      if (! (x == re))  error ("reduce() inconsistent with div()");

      if (! (c().quot(a,b) == q))  error ("quot() inconsistent with div()");

      x = a;
      c().quotient (x, b);
      if (! (x == q))  error ("quotient() inconsistent with div()");

      if (c().norm (re) >= c().norm (b))  error ("norm(rem()) >= norm(b)");

      fl();
   }

   // genGcd()

   {
      T ma, mb;
      T g = genGcd (c(), a, b, ma, mb);

      DEB2
      {
         COUT << ", gcd(a,b) = ";
         c.print (g);
         COUT << ", ma = ";
         c.print (ma);
         COUT << ", mb = ";
         c.print (mb);
      }

      if (c().is0 (g))
      {
         if (! c.is0 || ! c.is0B)  error ("gcd(a,b) is zero");
      }
      else
      {
         if (! c().isDivisor(a,g) || ! c().isDivisor(b,g))
         {
            error ("gcd(a,b) does not divide a or b");
         }
      }

      if (! (c().add (c().mul (a, ma), c().mul (b, mb)) == g))
      {
         error ("multiplicators wrong");
      }

      if (! (genGcd (c(), a, b) == g))  error ("genGcd()-3 broken");

      T ma2;
      if (! (genGcd (c(), a, b, ma2) == g && ma2 == ma))
      {
         error ("genGcd()-4 broken");
      }

      fl();
   }
}


/*****************************************************************************/
/*************  Field                                            *************/
/*****************************************************************************/


template<class R>
void doTest (Cache<R>& c, field_tag)
{
   requireCharacteristic (typename R::char_category());
   requireNozerodivisor  (typename R::zerodivisor_category());

   typedef typename R::type T;

   doTest (c, ringfield_tag());

   if (c.size)
   {
      unsigned prime;

      if (! Prime::isPrimePower (c.size, prime, c.extensionDegree))
      {
         error ("Size of finite field must be a prime power");
      }
      else
      {
         if (c.characteristic != prime)
         {
            error ("characteristic() and field size incompatible");
         }
      }
   }
}

template<class R>
void doTest (Cache<R>& c, const typename R::type& a, field_tag)
{
   typedef typename R::type T;

   doTest (c, a, ringfield_tag());

   if (! c.is0)
   {
      // multiplicative inverse

      c.recip = c().recip(a);
      DEB2
      {
         COUT <<
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 3
            L", 1\x2215" "a = ";  // DIVISION SLASH
# else
            L", 1/a = ";
# endif
#else
            Wgl4Ascii (", 1\xe2\x88\x95" "a = ", ", 1/a = ");
#endif
         c.print (c.recip);
      }

      if (! c().is1 (c().mul (a, c.recip)))  error ("a * (a^-1) != 1");

      {
         T rr = a;
         c().reciprocal (rr);

         if (! (rr == c.recip))  error ("reciprocal() broken");
      }

      fl();
   }
}

template<class R>
void doTest
   (Cache<R>& c, const typename R::type& a, const typename R::type& b,
    field_tag)
{
   typedef typename R::type T;

   doTest (c, a, b, ringfield_tag());

   if (! c.is0B)
   {
      T q = c().div(a, b);
      T r = c().recip (b);
      DEB2
      {
         COUT <<
#ifdef USE_WCHAR
# if HINTLIB_CHARACTER_SET >= 3
            L", a\x2215" "b = ";  // DIVISION SLASH
# else
            L", a/b = ";
# endif
#else
            Wgl4Ascii (", a\xe2\x88\x95" "b = ", ", a/b = ");
#endif
         c.print (q);
      }

      if (! (c().mul (a, r) == q))  error ("div() fails");

      T x = a;
      c().divBy (x, b);
      if (! (x == q)) error ("divBy() fails");

      fl();
   }
}


/*****************************************************************************/
/*************  Finite Field                                     *************/
/*****************************************************************************/

template<class R>
void doTest (Cache<R>& c, gf_tag)
{
   primitivesCounter = 0;

   typedef typename R::type T;

   doTest (c, field_tag());

   int extensionDegree = c().extensionDegree();

   DEB1  COUT << "  Degree: " << extensionDegree << endl;

   if (c.size == 0)
   {
      error ("size infinite in finite field");
   }
   else
   {
      if (c.extensionDegree != extensionDegree)
      {
         error ("extensionDegree() and field size incompatible");
      }
   }
}

template<class R>
void doTest (Cache<R>& c, const typename R::type& a, gf_tag)
{
   typedef typename R::type T;

   doTest (c, a, field_tag());

   // isPrimitiveElement()

   c.isPrimitiveElement = c().isPrimitiveElement(a);

   if (c.isPrimitiveElement)
   {
      DEB2  COUT << ", primitive element";
      ++primitivesCounter;
   }

   if (! c.is0)
   {
      if (c.isPrimitiveElement != (c.order == c.size - 1))
      {
         error ("isPrimitiveRoot() incompatible with order()");
      }
   }
   else  // a = 0
   {
      if (c.isPrimitiveElement)  error ("0 cannot be primitive element");
   }

   fl();
}

inline
void counterTest (CacheBase& c, bool all, gf_tag)
{
   counterTest (c, all, field_tag());

   // number of primitive elements

   checkCounter (Prime::eulerPhi (c.size - 1), primitivesCounter, all,
                 "number of primitive elements");
}


/*****************************************************************************/
/*************  Cyclic Field                                     *************/
/*****************************************************************************/


template<class R>
void doTest (Cache<R>& c, cyclic_tag)
{
   doTest (c, gf_tag());

   if (c.extensionDegree != 1)  error ("extensionDegree() != 1");
}


/*****************************************************************************/
/*************  Vectorspace                                      *************/
/*****************************************************************************/

template<class V>
void doTest (Cache<V>& c, vectorspace_tag)
{
   doTest (c, group_tag());
   
   typedef typename V::scalar_algebra A;
   typedef typename V::type T;
   typedef typename V::scalar_type ST;
   const A alg = c().getScalarAlgebra();

   c.dimension = c().dimension();
   c.baseSize  = alg.size();

   DEB1
   {
      COUT << "  Dim:  " << c.dimension << "\n"
              "  Alg:  " << alg << endl;
   }

   ensureEqualType (static_cast<ST*>(0), static_cast<typename A::type*>(0));

   if (c.baseSize == 0 && c.size != 0)
   {
      error("Finite vector space over infinte field");
   }
   if (c.baseSize != 0 && c.size == 0)
   {
      error("Infinite vector space over finte field");
   }
   if (c.baseSize != 0 && c.size != 0)
   {
      if (logInt (std::numeric_limits<unsigned>::max(), unsigned (c.baseSize))
          <= int (c.dimension))
      {
         if (c.size != powInt (c.baseSize, c.dimension))
            error ("dim (F^n) != (dim F)^n");
      }
   }

   Array<ST> coords (c.dimension);
   T x = T();
   c().toCoord (x, &coords[0]);

   for (int i = 0; i < c.dimension; ++i)
   {
      if (! alg.is0 (coords[i]))  error ("zero has non-zero coordinates");
   }
}

template<class V>
void doTest (Cache<V>& c, const typename V::type& a, vectorspace_tag)
{
   doTest (c, a, group_tag());
   
   typedef typename V::scalar_algebra A;
   typedef typename V::type T;
   typedef typename V::scalar_type ST;
   const A& alg = c().getScalarAlgebra();

   // toCoords()
   
   Array<ST> coords (c.dimension);
   c().toCoord (a, &coords[0]);

   T aCopy = a;

   DEB2
   {
      COUT << ", coord = [";
      for (int i = 0; i < c.dimension; ++i)
      {
         if (i > 0)  COUT << ';';
         alg.printShort (COUT, coords[i]);
      }
      COUT << ']';
   }

   for (int i = 0; i < c.dimension; ++i)
   {
      checkIndex (c.baseSize, alg.index(coords[i]), "coords");
      if (! (coords[i] == c().coord(a,i)))
      {
         error("coord() and toCoords() mismatch");
      }
      typename V::scalar_reference ref = c().coord(aCopy, i);
      if (! (coords[i] == ref))  error("coord() and toCoords() mismatch");
   }

   // fromCoords()

   T b;
   c().fromCoord (b, &coords[0]);
   if (! (a == b))  error ("fromCoord() broken");

   // coords()

   b = T();
   for (int i = 0; i < c.dimension; ++i)  c().coord(b, i) = coords[i];
   if (! (a == b))  error ("Assignment to coord() broken");

   b = T();
   for (int i = c.dimension - 1; i >= 0; --i)  c().coord(b, i) = coords[i];
   if (! (a == b))  error ("Assignment to coord() broken");

   // Multiplication with 1

   b = c().mul (a, alg.one());

   if (! (b == a))  error ("1v != v");

   // Multiplication with 0

   b = c().mul (a, ST());

   if (! c().is0 (b))  error ("0v != 0");

   //  neg()

   b = c.neg;
   for (int i = 0; i < c.dimension; ++i)
   {
      if (! (c().coord(b,i) == alg.neg (c().coord(a,i))))
      {
         error ("neg() for vector does not match neg() for coordinates");
      }
   }

   // mul()

   if (sqr (c.index) < SIZE)
   {
      unsigned SS = unsigned (HINTLIB_MN sqrt(double(SIZE)));
      if (c.baseSize)  SS = std::min (SS, c.baseSize);

      for (unsigned l = 0; l < SS; ++l)
      {
         ST scal = alg.element (l);

         b = c().mul (a, scal);

         for (int i = 0; i < c.dimension; ++i)
         {
            if (! (c().coord(b,i) == alg.mul (c().coord(a,i), scal)))
            {
               error ("mul() for vector does not match mul() for coordinates");
            }
         }
      }
   }

   // Associativity and Distributivity (scalar, scalar, vector)

   if (cube (c.index) < SIZE)
   {
      unsigned SS = unsigned (HINTLIB_MN pow(double(SIZE), 1.0/3.0));
      if (c.baseSize)  SS = std::min (SS, c.baseSize);

      for (unsigned l1 = 0; l1 < SS; ++l1)
      {
         ST scal1 = alg.element (l1);

         for (unsigned l2 = 0; l2 < SS; ++l2)
         {
            ST scal2 = alg.element (l2);

            if (! (   c().mul (a, alg.mul (scal1, scal2))
                   == c().mul (c().mul (a, scal1), scal2)))
            {
               error ("(ab)v != a(bv)");
            }

            if (! (   c().mul (a, alg.add (scal1, scal2))
                   == c().add (c().mul(a, scal1), c().mul(a, scal2))))
            {
               error ("(a+b)v != av + bv");
            }
         }
      }
   }

   fl();
}

template<class V>
void doTest
   (Cache<V>& c, const typename V::type& a, const typename V::type& b,
    vectorspace_tag)
{
   doTest (c, a, b, group_tag());
   
   typedef typename V::scalar_algebra A;
   typedef typename V::type T;
   typedef typename V::scalar_type ST;
   const A& alg = c().getScalarAlgebra();

   // add()

   for (int i = 0; i < c.dimension; ++i)
   {
      if (! (c().coord(c.add, i) == alg.add (c().coord(a,i), c().coord(b,i))))
      {
         error ("add() for vector does not match add() for coordinates");
      }
   }

   // Distributivity (scalar, vector, vector)

   if (cube (c.index) < SIZE && cube (c.indexB) < SIZE)
   {
      unsigned SS = unsigned (HINTLIB_MN pow(double(SIZE), 1.0/3.0));
      if (c.baseSize)  SS = std::min (SS, c.baseSize);

      for (unsigned l = 0; l < SS; ++l)
      {
         ST scal = alg.element (l);

         if (! (   c().mul (c.add, scal)
                == c().add (c().mul (a, scal), c().mul (b, scal))))
         {
            error ("a (v1+v2) != (a v1) + (b + v2)");
         }
      }
   }

   fl();
}


/*************  Names of algebraic structures  *******************************/

#ifdef USE_WCHAR
inline const CHAR* typeName (group_tag)       { return L"group"; }
inline const CHAR* typeName (ring_tag)        { return L"ring"; }
inline const CHAR* typeName (domain_tag)      { return L"integral domain"; }
inline const CHAR* typeName (ufd_tag)         { return L"UFD"; }
inline const CHAR* typeName (euclidean_tag)   { return L"Euclidean domain"; }
inline const CHAR* typeName (integer_tag)     { return L"Euclidean domain of integers"; }
inline const CHAR* typeName (field_tag)       { return L"field"; }
inline const CHAR* typeName (numberfield_tag) { return L"number field"; }
inline const CHAR* typeName (rational_tag)    { return L"field of rational numbers"; }
inline const CHAR* typeName (real_tag)        { return L"field of real numbers"; }
inline const CHAR* typeName (complex_tag)     { return L"field of complex numbers"; }
inline const CHAR* typeName (gf_tag)          { return L"finite field"; }
inline const CHAR* typeName (cyclic_tag)      { return L"cyclic field"; }
inline const CHAR* typeName (funfield_tag)    { return L"function field"; }
inline const CHAR* typeName (ratfunfield_tag) { return L"rational function field"; }
inline const CHAR* typeName (vectorspace_tag) { return L"vector space"; }
#else
inline const CHAR* typeName (group_tag)       { return "group"; }
inline const CHAR* typeName (ring_tag)        { return "ring"; }
inline const CHAR* typeName (domain_tag)      { return "integral domain"; }
inline const CHAR* typeName (ufd_tag)         { return "UFD"; }
inline const CHAR* typeName (euclidean_tag)   { return "Euclidean domain"; }
inline const CHAR* typeName (integer_tag)     { return "Euclidean domain of integers"; }
inline const CHAR* typeName (field_tag)       { return "field"; }
inline const CHAR* typeName (numberfield_tag) { return "number field"; }
inline const CHAR* typeName (rational_tag)    { return "field of rational numbers"; }
inline const CHAR* typeName (real_tag)        { return "field of real numbers"; }
inline const CHAR* typeName (complex_tag)     { return "field of complex numbers"; }
inline const CHAR* typeName (gf_tag)          { return "finite field"; }
inline const CHAR* typeName (cyclic_tag)      { return "cyclic field"; }
inline const CHAR* typeName (funfield_tag)    { return "function field"; }
inline const CHAR* typeName (ratfunfield_tag) { return "rational function field"; }
inline const CHAR* typeName (vectorspace_tag) { return "vector space"; }
#endif


/*************  Main test procedure  *****************************************/

template<class A>
void doTests(A r, const char* type)
{
   typedef typename A::algebra_category cat;
   typedef typename A::type T;

   STRINGSTREAM ss;
#ifdef HINTLIB_STREAMS_SUPPORT_LOCALE
   ss.imbue (COUT.getloc());
#endif
   ss << r;

   if (! performTest (typeName(cat()), ss.str(), type))  return;

   // Initialize cache

   Cache<A> c (r);

   c.size = c().size();
   unsigned S = SIZE;
   if (c.size > 0)  S = std::min (S, c.size);

   // No variable

   DEB1
   {
      COUT << "  Size: ";
      printNumberOrInf (c.size);
      COUT << endl;
   }

   doTest (c, typename A::size_category());
   doTest (c, cat());

   // One variable

   for (unsigned i = 0; i < S; i++)
   {
      T a;
      try
      {
         a = r.element (i);
      }
      catch (Overflow &)
      {
         break;
      }

      DEB1 { COUT << setw(W); c.print (a); }

      // element() and index()

      c.index = r.index(a);
      DEB2 COUT << ", index = " << c.index;
      fl();
      if (c.index != i)  error ("index(a) != i");

      doTest (c, a, cat());

      DEB1 COUT << endl;
   }

   counterTest (c, S == c.size, cat());

   // two variables

   S = unsigned(HINTLIB_MN sqrt(double(SIZE)));
   if (c.size > 0)  S = std::min (S, c.size);

   for (unsigned i = 0; i < S; i++)
   {
      const T a = r.element (i);

      for (unsigned j = 0; j < S; j++)
      {
         const T b = r.element (j);

         DEB1
         {
            COUT << setw(W); c.print (a);
            COUT << setw(W); c.print (b);
            fl();
         }

         // Make sure, == makes sense

         if ((a == b) != (i == j))  error ("a==b incorrect");

         doTest (c, a, b, cat());

         DEB1 COUT << endl;
      }
   }

   // three variables

   S = unsigned(HINTLIB_MN pow(double(SIZE), 1.0 / 3.0));
   if (c.size > 0)  S = std::min (S, c.size);

   for (unsigned i = 0; i < S; i++)
   {
      const T a = r.element (i);

      for (unsigned j = 0; j < S; j++)
      {
         const T b = r.element (j);

         for (unsigned k = 0; k < S; ++k)
         {
            const T d = r.element (k);

            DEB1
            {
               COUT << setw(W); c.print (a);
               COUT << setw(W); c.print (b);
               COUT << setw(W); c.print (d);
            }

            doTest (c, a, b, d, cat());

            DEB1 COUT << endl;
         }
      }
   }
}

