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

#ifndef HINTLIB_POLYNOMIAL_RATIONAL_TCC
#define HINTLIB_POLYNOMIAL_RATIONAL_TCC 1

#ifndef HINTLIB_POLYNOMIAL_FIELD_TCC
#include <HIntLib/polynomial_field.tcc>
#endif

#include <HIntLib/counter.h>


/********************  Polynomial Ring <rational_tag>  ***********************/


/**
 *  isPrime()
 *
 *  We use Kronecker's method for determining if there is a non-trivial
 *  factorization.
 *
 *  See Lidl/Niederreiter, Finite Fields, Exercise 1.30, and
 *      Knuth, TACP, vol 2, 4.6.2, p. 449-450 for details.
 */

template<class A>
bool
HIntLib::Private::PRBA_Rational<A>::isPrime (const type& p) const
{
   const int degree = p.degree();
   const A& aa (this->a);

   if (degree <= 0)  return false;
   if (degree == 1)  return true;  //  bx+a  is always prime

   // ...+ cx^2 + bx + 0  is never prime

   if (aa.is0 (p.ct()))  return false;

   // Make sure p is square-free

   if (! isSquarefree (p))  return false;
   
   // we have to check for divisors up to degree  s

   int s = degree / 2;

   // Integer algebra and type

   typedef typename A::base_algebra int_algebra;
   int_algebra intAlg = aa.getBaseAlgebra();
   typedef typename A::base_type int_type;

   // Integer polynomial algebra and type 

   typedef PolynomialRing<int_algebra> intPoly_algebra;
   intPoly_algebra intPolyAlg (intAlg);
   typedef typename intPoly_algebra::type intPoly_type;

   // Make a copy of p with integer coefficients
   
   int_type gcdNum = p.ct().numerator();
   int_type lcmDen = p.ct().denominator();

   for (int i = 1; i <= degree; ++i)
   {
      if (! aa.is0 (p[i]))
      {
         gcdNum = genGcd (intAlg, gcdNum, p[i].numerator());
         lcmDen = genLcm (intAlg, lcmDen, p[i].denominator());
      }
   }
   
   intPoly_type intp; intp.reserve (degree + 1);

   for (int i = degree; i >= 0; --i)
   {
      if (aa.is0 (p[i]))
      {
         intp.mulByX();
      }
      else
      {
         intp.mulAndAdd (
            intAlg.mul (intAlg.quot (p[i].numerator(), gcdNum),
                        intAlg.quot (lcmDen, p[i].denominator())));
      }
   }

#if 0
cout << endl
     << "degree = " << degree << endl
     << "s = " << s << endl
     << "p = ";
printShort (cout, p);
cout << endl << "pi= ";
intPolyAlg.printShort (cout, intp);
cout << endl << "gcdNum = ";
intAlg.printShort (cout, gcdNum);
cout << endl << "lcmDen = ";
intAlg.printShort (cout, lcmDen);
cout << endl;
#endif

   // Determine divisors of f(a)
   
   Array<std::vector<int_type> > divisors (s + 1);
   Array<int_type> nodes (s + 1);

   for (int num = 0, index = 0; num <= s; ++index)
   {
      const int_type ii = intAlg.element (index);
      const int_type value = intPolyAlg.evaluate (intp, ii);

      if (intAlg.is0 (value))  continue;

      nodes [num] = ii;

#if 0
cout << "node " << num << ": f(";
intAlg.printShort (cout, ii);
cout << ") = ";
intAlg.printShort (cout, value);
#endif

      divisors[num].push_back (intAlg.one());
      if (num) divisors[num].push_back (intAlg.neg (intAlg.one()));

      unsigned normValue = intAlg.norm (value);

      if (normValue > 1)
      {
         divisors[num].push_back (value);
         if (num) divisors[num].push_back (intAlg.neg (value));
   
         for (unsigned j = 2; ; ++j)
         {
            int_type candidate (j);

            if (intAlg.norm (intAlg.mul (candidate, candidate)) > normValue)
            {
               break;
            }

            if (intAlg.isDivisor (value, candidate))
            {
               divisors[num].push_back (candidate);
               if (num) divisors[num].push_back (intAlg.neg (candidate));

               int_type candidate2 = intAlg.div (value, candidate);

               if (candidate2 != candidate)
               {
                  divisors[num].push_back (candidate2);
                  if (num) divisors[num].push_back (intAlg.neg (candidate2));
               }
            }
         }
      }

#if 0
cout << "   divisors: ";
for (unsigned j = 0; j < divisors[num].size(); ++j)
{
   intAlg.printShort (cout, divisors[num][j]);
   cout << ", ";
}
cout << endl;
#endif

      ++num;
   }

   // create Basis for Lagrange polynomials

   Array<type> lagrangeBasis (s+1);
   
   for (int i = 0; i <= s; ++i)
   {
      // construct denominator

      int_type den = intAlg.one();

      for (int j = 0; j <= s; ++j)
      {
         if (i == j)  continue;

         intAlg.mulBy (den, intAlg.sub (nodes[i], nodes[j]));
      }

      type pp; pp.reserve (s+1);
      pp.mulAndAdd (aa.makeElement (intAlg.one(), den));

      // construct numerator

      for (int j = 0; j <= s; ++j)
      {
         if (i == j)  continue;

         type copy (pp);
         pp.mulByX();
         mulBy (copy, aa.makeElement (intAlg.neg (nodes[j])));
         addTo (pp, copy);
      }

      lagrangeBasis [i] = pp;
   }

   // Enumerate all s+1 - tupels

   Array<unsigned> numDivisors (s + 1);
   for (int i = 0; i <= s; ++i)  numDivisors[i] = divisors[i].size();
   CounterMixedBase counter (numDivisors.begin(), numDivisors.begin() + s + 1);

   do
   {
      // create Lagrange interpolation

      type candidate;

      for (int i = 0; i <= s; ++i)
      {
         addTo (candidate, mul (lagrangeBasis[i],
                                aa.makeElement (divisors[i][counter[i]])));
      }

      // is it a (non-trivial) divisor?

#if 0
cout << "candidate = ";
printShort (cout, candidate);
if (candidate.degree() >= 1)
{
   cout << "   rem = ";
   printShort (cout, rem (p, candidate));
}
cout << endl;
#endif

      if (! candidate.isConstant() && isDivisor (p, candidate))  return false;
   }
   while (counter.next());
   
   return true;
} 

/**
 *  next()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::NextRational,A,T> x)
{
   const Private::PRBA_Rational<A>* a = x.a;
   unsigned* np = x.np;

   for (;;)
   {
      P p = a->elementMonic ((*np)++);
      if (p.degree() >= 6)  throw Overflow();

      if (a->isPrime (p))
      {
         swap (p);
         return;
      }
   }
}


/*********************  Instantiations  **************************************/


// Instantiate PolynomialRingBase<rational_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_RATIONAL(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_FIELD(X) \
   template Polynomial<X::type>::Polynomial \
      (Private::PG<Private::NextRational,X >);\
   namespace Private { \
   template bool PRBA_Rational<X >::isPrime (const type&) const; \
   }

#endif

