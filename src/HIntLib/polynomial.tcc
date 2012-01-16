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

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_OSTREM
  #include <ostream>
#else
  #include <iostream>
#endif

#ifdef HINTLIB_HAVE_SSTREAM
  #include <sstream>
#else
  #include <HIntLib/fallback_sstream.h>
#endif

#include <HIntLib/polynomial.h>

#include <HIntLib/gcd.tcc>
#include <HIntLib/counter.h>
#include <HIntLib/linearalgebra.h>

namespace L = HIntLib;


/********************** Polynomial <> ****************************************/

/**
 *  Assignment
 */

template<typename T>
L::Polynomial<T>&
L::Polynomial<T>::operator= (const L::Polynomial<T>& p)
{
   if (this != &p) c = p.c;
   return *this;
}


/**
 *  divByX()
 *  mulByX()
 */

template<typename A>
L::Polynomial<A>&
L::Polynomial<A>::divByX (unsigned k)
{
   for (unsigned i = 0; i < k; ++i)  divByX();
   return *this;
}

template<typename A>
L::Polynomial<A>&
L::Polynomial<A>::mulByX (unsigned k)
{
   for (unsigned i = 0; i < k; ++i)  mulByX();
   return *this;
}


/**
 *  operator==()  (specialization for T = real)
 */

namespace HIntLib
{
template<typename T>
bool operator== (
      const Polynomial<Real<T> >& p1, const Polynomial<Real<T> >& p2)
{
   // determine magnitude of coefficients

   int deg1 = p1.degree();
   int deg2 = p2.degree();
   int deg = std::max (deg1, deg2);

   T mag = 0.0;

   for (int i = 0; i <= deg; ++i)
   {
      if (i <= deg1 && mag < abs (T(p1[i])))  mag = abs (T(p1[i]));
      if (i <= deg2 && mag < abs (T(p2[i])))  mag = abs (T(p2[i]));
   }

   mag *= (deg + 1) * std::numeric_limits<T>::epsilon() * 100.0;

   for (int i = 0; i <= deg; ++i)
   {
      T x1 = (i <= deg1) ? T(p1[i]) : T(0);
      T x2 = (i <= deg2) ? T(p2[i]) : T(0);
      if (abs (x1 - x2) > mag)  return false;   
   }

   return true;
}
}


/*******************  Polynomail Ring Base Base  *****************************/

/**
 *  x ()
 */

template<typename A>
typename L::PolynomialRingBB<A>::type
L::PolynomialRingBB<A>::x (unsigned k) const
{
   type p (k + 1);
   p.mulAndAdd (a.one());
   p.mulByX (k);
   return p;
}


/**
 *  is1 ()
 */

template<typename A>
bool
L::PolynomialRingBB<A>::is1 (const type &p) const
{
   return p.degree() == 0 && a.is1(p.lc());
}


/**
 *  isCanonical ()
 */

template<class A>
bool
L::PolynomialRingBB<A>::isCanonical (const type &p) const
{
   return p.degree() < 0 || a.is1(p.lc());
}


/**
 *  element()
 */

template<typename A>
typename L::PolynomialRingBB<A>::type
L::PolynomialRingBB<A>::element(unsigned i) const
{
   if (! i)  return type();

   unsigned s = a.size()  ?  a.size()  : 5;

   unsigned ii = i;
   unsigned num = 0;
   while (ii)
   {
      ++num;
      ii /= s;
   }
   
   type p = x (num - 1);

   unsigned j = 0;
   while (i)
   {
      p[j++] = a.element (i % s);
      i /= s;
   }
   return p;
}


/**
 *  index()
 */

template<typename A>
unsigned
L::PolynomialRingBB<A>::index (const type &p) const
{
   unsigned s = a.size()  ?  a.size()  : 5;

   unsigned n = 0;
   for (int i = p.degree(); i >= 0; --i)  n = s * n + a.index (p[i]);
   return n;
}


/**
 *  element Monic ()
 */

template<typename A>
typename L::PolynomialRingBB<A>::type
L::PolynomialRingBB<A>::elementMonic(unsigned ind) const
{
   if (! ind)  return type();
   if (ind == 1)  return one();

   const unsigned s = a.size()  ?  a.size()  : 5;

   ind -= 2;
   unsigned max = s;
   unsigned numCoeffs = 1;

   while (ind >= max)
   {
      ind -= max;
      max *= s;
      ++numCoeffs;
   }

   type p = x (numCoeffs);

   for (unsigned i = 0; i < numCoeffs; ++i)
   {
      p[i] = a.element (ind % s);
      ind /= s;
   }

   return p;
}


/**
 *  negate()
 */

template<typename A>
typename L::PolynomialRingBB<A>::type &
L::PolynomialRingBB<A>::negate (type& p) const
{ 
   for (int i = p.degree(); i >= 0; --i)  a.negate (p[i]);
   return p;
}


/**
 *  neg()
 */

template<typename A>
typename L::PolynomialRingBB<A>::type
L::PolynomialRingBB<A>::neg (const type& p) const
{
   type n (p.degree() + 1);
   for (int i = p.degree(); i >= 0; --i)  n.mulAndAdd (a.neg (p[i]));
   return n;
}


/**
 *  add()
 */

template<typename A>
typename L::PolynomialRingBB<A>::type
L::PolynomialRingBB<A>::add (const type& p1, const type& p2) const
{
   const int deg1 = p1.degree();
   const int deg2 = p2.degree();
   const int maxDeg = std::max (deg1, deg2); 

   type p (maxDeg + 1);
   int start;
   
   if (deg1 != deg2)
   {
      start = std::min (deg1, deg2);
      const type& greater = deg1 > deg2 ? p1 : p2; 

      for (int i = maxDeg; i > start; --i)  p.mulAndAdd (greater[i]);
   }
   else
   {
      start = deg1;

      while (start >= 0 && a.is0 (a.add (p1[start], p2[start])))  --start;
   }

   for (int i = start; i >= 0; --i)  p.mulAndAdd (a.add (p1[i], p2[i]));

   return p;
}


/**
 *  times()
 */

template<typename A>
typename L::PolynomialRingBB<A>::type
L::PolynomialRingBB<A>::times (const type& p, unsigned k) const
{
   type n (p.degree() + 1);
   bool nonZero = false;

   for (int i = p.degree(); i >= 0; --i)
   {
      coeff_type c = a.times (p[i], k);

      if (nonZero || ! a.is0 (c))
      {
         nonZero = true;
         n.mulAndAdd (c);
      }
   }
   return n;
}


/**
 *  derivative()
 */

template<typename A>
typename L::PolynomialRingBB<A>::type
L::PolynomialRingBB<A>::derivative (const type& p) const
{
   if (p.degree() < 0)  return type();
   
   type n (p.degree());
   bool nonZero = false;

   for (int i = p.degree(); i > 0; --i)
   {
      coeff_type c = a.times (p[i], i);

      if (nonZero || ! a.is0 (c))
      {
         nonZero = true;
         n.mulAndAdd (c);
      }
   }
   return n;
}


/**
 *  evaluate()
 */

template<typename A>
typename A::type
L::PolynomialRingBB<A>::evaluate (const type& p, const coeff_type& x) const
{
   if (is0 (p))  return coeff_type();
   if (a.is0 (x))  return p[0];

   coeff_type res = coeff_type();

   if (a.is1 (x))
   {
      for (int i = p.degree(); i >= 0; --i)  a.addTo (res, p[i]);
   }
   else
   {
      for (int i = p.degree(); i >= 0; --i)  res = a.add (a.mul (res, x), p[i]);
   }

   return res;
}


/**
 *  operator<<
 */

template<typename A>
std::ostream&
L::operator<< (std::ostream &o, const PolynomialRingBB<A> &a)
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   ss << a.getCoeffAlgebra() << '[' << a.getVar() << ']';

   return o << ss.str().c_str();
}

template<typename A>
void
L::PolynomialRingBB<A>::printShort (std::ostream& o, const type& p) const
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   bool output = false;

   for (int i = p.degree(); i >= 0; --i)
   {
      typename A::type coef = p[i];
      if (a.is0(coef))  continue;

      bool coefPrinted = false;

      if (a.is1(coef))
      {
         if (output)  ss << '+';
      }
      else if (a.is1(a.neg(coef)))
      {
         ss << '-';
      }
      else
      {
         if (output)  ss << '+';
         a.printShort (ss, coef);
         coefPrinted = true;
      }

      switch (i)
      {
         case 0:  if (! coefPrinted)  a.printShort (ss, a.one());
                  break;
         case 1:  ss << var; break;
         case 2:  ss << var << '\262'; break;
         case 3:  ss << var << '\263'; break;
         default: ss << var << '^' << i;
      }

      output = true;
   }

   if (! output)  a.printShort (ss, coeff_type());

   o << ss.str().c_str();
}

template<typename A>
void
L::PolynomialRingBB<A>::print (std::ostream& o, const type& p) const
{
   std::ostringstream ss;
   ss.flags (o.flags());
   ss.precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCAL
   ss.imbue (o.getloc());
#endif

   printShort (ss, p);
   ss << ' ';
   printSuffix (ss);

   o << ss.str().c_str();
}


/********************  Polynomial Ring <ring_tag>  ***************************/


/**
 *  additive Order()
 */

template<class A>
unsigned
L::PolynomialRingBase<A,L::ring_tag>::additiveOrder (const type& u) const
{
   unsigned n = 1;
   for (int i = 0; i <= u.degree(); ++i)  n = lcm (n, a.additiveOrder (u[i]));
   return n;
}


/**
 *  mul()
 */

template<typename A>
typename L::PolynomialRingBase<A,L::ring_tag>::type
L::PolynomialRingBase<A,L::ring_tag>::mul (
      const type& _p1, const type& _p2) const
{
   const type& p1 = _p1.degree() > _p2.degree()  ?  _p1 : _p2;
   const type& p2 = _p1.degree() > _p2.degree()  ?  _p2 : _p1;
   
   const int deg1 = p1.degree();
   const int deg2 = p2.degree();

   if (deg1 < 0 || deg2 < 0)  return type(0);

   type p (deg1 + deg2 + 1);
   bool nonZero = false;

   for (int k = deg1 + deg2; k >= 0; --k)
   {
      coeff_type c = coeff_type();

      for (int i = std::max (k - deg2, 0); i <= std::min (k, deg1); ++i)
      {
         a.addTo (c, a.mul (p1[i], p2[k - i]));
      }

      if (nonZero || ! a.is0 (c))
      {
         nonZero = true;
         p.mulAndAdd (c);
      }
   }

   return p;
}


/********************  Polynomial Ring <domain_tag> and subclasses  **********/


/**
 *  mul()
 *
 *  For polynomials over domains, the leading coefficients can not cancel out
 *  each other.
 */

template<typename A>
typename L::PolynomialRingBase<A,L::domain_tag>::type
L::PolynomialRingBase<A,L::domain_tag>::mul (
      const type& _p1, const type& _p2) const
{
   const type& p1 = _p1.degree() > _p2.degree()  ?  _p1 : _p2;
   const type& p2 = _p1.degree() > _p2.degree()  ?  _p2 : _p1;
   
   const int deg1 = p1.degree();
   const int deg2 = p2.degree();

   if (deg1 < 0 || deg2 < 0)  return type();

   type p (deg1 + deg2 + 1);

   for (int k = deg1 + deg2; k >= 0; --k)
   {
      coeff_type c = coeff_type();

      for (int i = std::max (k - deg2, 0); i <= std::min (k, deg1); ++i)
      {
         a.addTo (c, a.mul (p1[i], p2[k - i]));
      }

      p.mulAndAdd (c);
   }

   return p;
}


/**
 *  order()
 */

template<typename A>
unsigned
L::PolynomialRingBase<A,L::domain_tag>::order (const type& p) const
{
   switch (p.degree())
   {
   case -1: throwDivisionByZero();
   case  0: return a.order (p.lc());  // no x, no polynomial
   }

   // if there is an x, it will never go away if A is a domain
   
   return 0;
}


/********************  Polynomial Ring <field_tag> and subclasses  ***********/


/**
 *  div ()
 *
 *  Division of polynomials
 *
 *  See Knuth, TACP, vol 2, 4.6.1, Algo D
 */

template<class A>
void L::PolynomialRingBase<A,L::field_tag>::div
   (const type& u, const type& v, type& _q, type& _r) const
{
   const int vdeg = v.degree();
   if (vdeg == -1)  throwDivisionByZero();
   const int udeg = u.degree();

   if (udeg < vdeg)
   {
      _q = type();
      _r = u;
      return;
   }

   coeff_type qq = a.recip (v.lc());
   type uu (u);
   type q (udeg - vdeg + 1);

   for (int k = udeg-vdeg; k >= 0; --k)
   {
      coeff_type lc = a.mul (uu[vdeg+k], qq);
      q.mulAndAdd (lc);

      for (int j = vdeg + k - 1; j >= k; --j)
      {
         a.subFrom (uu[j], a.mul (lc, v[j-k]));
      }
   }

   _q = q;
   
   bool nonZero = false;
   type r (vdeg);
   for (int i = vdeg-1; i >= 0; --i)
   {
      if (nonZero || ! a.is0(uu[i]))
      {
         r.mulAndAdd (uu[i]);
         nonZero = true;
      }
   }
   _r = r;
}


/**
 *  quot()
 */

template<class A>
typename L::PolynomialRingBase<A,L::field_tag>::type
L::PolynomialRingBase<A,L::field_tag>::quot (const type& u, const type& v) const
{
   const int vdeg = v.degree();
   if (vdeg == -1)  throwDivisionByZero();
   const int udeg = u.degree();

   if (udeg < vdeg)  return type(0);

   coeff_type qq = a.recip (v.lc());
   type uu (u);
   type q (udeg - vdeg + 1);

   for (int k = udeg-vdeg; k >= 0; --k)
   {
      coeff_type lc = a.mul (uu[vdeg+k], qq);
      q.mulAndAdd (lc);

      for (int j = vdeg + k - 1; j >= k; --j)
      {
         a.subFrom (uu[j], a.mul (lc, v[j-k]));
      }
   }

   return q;
}


/**
 *  rem()
 */

template<class A>
typename L::PolynomialRingBase<A,L::field_tag>::type
L::PolynomialRingBase<A,L::field_tag>::rem (const type& u, const type& v) const
{
   const int vdeg = v.degree();
   if (vdeg == -1)  throwDivisionByZero();
   const int udeg = u.degree();

   if (udeg < vdeg)  return u;

   coeff_type qq = a.recip (v.lc());
   type uu (u);

   for (int k = udeg-vdeg; k >= 0; --k)
   {
      coeff_type lc = a.mul (uu[vdeg+k], qq);

      for (int j = vdeg + k - 1; j >= k; --j)
      {
         a.subFrom (uu[j], a.mul (lc, v[j-k]));
      }
   }
   
   bool nonZero = false;
   type r (vdeg);
   for (int i = vdeg-1; i >= 0; --i)
   {
      if (nonZero || ! a.is0(uu[i]))
      {
         r.mulAndAdd (uu[i]);
         nonZero = true;
      }
   }

   return r;
}


/**
 *  num Of Remainders()
 */

template<class A>
unsigned
L::PolynomialRingBase<A,L::field_tag>::numOfRemainders (const type& p) const
{
   return a.size()  ?  powInt (a.size(), p.degree()) : 0;
}


/**
 *  make Canonical ()
 */

template<typename A>
typename L::PolynomialRingBase<A,L::field_tag>::unit_type
L::PolynomialRingBase<A,L::field_tag>::makeCanonical (type &p) const
{
   if (is0 (p)) return a.one();

   unit_type l = p.lc();
   unit_type il = a.recip (l);

   p.lc() = a.one();
   for (int i = 0; i < p.degree(); ++i)  a.mulBy (p[i], il);

   return l;
}


/**
 *  mulUnit ()
 */

template<typename A>
typename L::PolynomialRingBase<A,L::field_tag>::type
L::PolynomialRingBase<A,L::field_tag>::mulUnit
   (const type &p, const unit_type& u) const
{
   type pp (p.degree() + 1);
   for (int i = p.degree(); i >= 0; --i)  pp.mulAndAdd (a.mul (p[i], u));
   return pp;
}


/**
 *  mulByUnit()
 */

template<typename A>
typename L::PolynomialRingBase<A,L::field_tag>::type &
L::PolynomialRingBase<A,L::field_tag>::mulByUnit
   (type &p, const unit_type& u) const
{
   for (int i = 0; i <= p.degree(); ++i)  a.mulBy (p[i], u);
   return p;
}

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
L::PolynomialRingBase<A,L::rational_tag>::isPrime (const type& p) const
{
   const int degree = p.degree();

   if (degree <= 0)  return false;
   if (degree == 1)  return true;  //  bx+a  is always prime

   if (a.is0 (p[0]))  return false; // ...+ cx^2 + bx + 0  is never prime

   // Make sure p is square-free

   if (! isUnit (genGcd (*static_cast<const PolynomialRing<A>*>(this),
                         p, derivative(p))))  return false;
   
   // we have to check for divisors up to degree  s

   int s = degree / 2;

   // Integer algebra and type

   typedef typename A::base_algebra int_algebra;
   int_algebra intAlg = a.getBaseAlgebra();
   typedef typename A::base_type int_type;

   // Integer polynomial algebra and type 

   typedef PolynomialRing<int_algebra> intPoly_algebra;
   intPoly_algebra intPolyAlg (intAlg);
   typedef typename intPoly_algebra::type intPoly_type;

   // Make a copy of p with integer coefficients
   
   int_type gcdNum = p[0].numerator();
   int_type lcmDen = p[0].denominator();

   for (int i = 1; i <= degree; ++i)
   {
      if (! a.is0 (p[i]))
      {
         gcdNum = genGcd (intAlg, gcdNum, p[i].numerator());
         lcmDen = genLcm (intAlg, lcmDen, p[i].denominator());
      }
   }
   
   intPoly_type intp (degree + 1);

   for (int i = degree; i >= 0; --i)
   {
      if (a.is0 (p[i]))
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
      int_type ii = intAlg.element (index);
      int_type value = intPolyAlg.evaluate (intp, ii);

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
            int_type candidate = j;

            if (intAlg.norm (intAlg.mul (candidate, candidate)) > normValue)
            {
               break;
            }

            if (intAlg.is0 (intAlg.rem (value, candidate)))
            {
               divisors[num].push_back (candidate);
               if (num) divisors[num].push_back (intAlg.neg (candidate));

               if (intAlg.quot (value, candidate) != candidate)
               {
                  int_type candidate2 = intAlg.quot (value, candidate);
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
      // denominator

      int_type den = intAlg.one();

      for (int j = 0; j <= s; ++j)
      {
         if (i == j)  continue;

         intAlg.mulBy (den, intAlg.sub (nodes[i], nodes[j]));
      }

      type pp (s + 1);
      pp.mulAndAdd (a.makeElement (intAlg.one(), den));

      for (int j = 0; j <= s; ++j)
      {
         if (i == j)  continue;

         mulBy (pp,
           type(2).mulAndAdd (a.one()).
                   mulAndAdd (a.makeElement (intAlg.neg (nodes[j]))));
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

      type candidate (s + 1);

      for (int i = 0; i <= s; ++i)
      {
         addTo (candidate, mulUnit (lagrangeBasis[i],
                                    a.makeElement (divisors[i][counter[i]])));
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

      if (candidate.degree() >= 1 && is0 (rem (p, candidate)))  return false;
   }
   while (counter.next());
   
   return true;
} 

/**
 *  next()
 */

template<class A>
typename L::PolynomialRingBase<A,L::rational_tag>::type
L::PolynomialRingBase<A,L::rational_tag>::PrimeGenerator::next()
{
   type p;

   do
   {
      p = alg.elementMonic (n++);
   }
   while (! alg.isPrime (p));

   return p;
}


/********************  Polynomial Ring <real_tag>  ***************************/


/**
 *  isPrime()
 */

template<class A>
bool
L::PolynomialRingBase<A,L::real_tag>::isPrime (const type& p) const
{
   const int degree = p.degree();

   if (degree == 1)  return true;  //  bx+a  is always prime

   // 0 is zero
   // constants are units
   // polynomials of degree >= 3 as well as  ax^2+bx  can always be factored
   //    over the reals

   if (degree <= 0 || degree >= 3 || a.is0(p[0]))  return false;

   // Calculate discriminant to determine if there is a (real) solution

   const coeff_type disc = sqr (p[1]) - coeff_type(4) * p[2] * p[0];
   
   return disc < 0.0 && ! a.is0 (disc);
} 


/**
 *  next()
 */

template<class A>
typename L::PolynomialRingBase<A,L::real_tag>::type
L::PolynomialRingBase<A,L::real_tag>::PrimeGenerator::next()
{
   unsigned nn = n++;

   type p (3);

   if (nn & 1)
   {
      // create polynomial  x^2 + bx + c
      // with arbitrary  b, and  c  such that  b^2 < 4c

      unsigned ib, ic;
      unthread (nn / 2, ib, ic);

      real b = alg.element (ib);  // arbitrary real number
      real c = sqr(b) / 4.0 + alg.element (ic * 2 + 1);  // (x > 0) + b^2 / 4

      p.mulAndAdd (1.0).mulAndAdd (b).mulAndAdd (c);
   }
   else
   {
      // create polynomial  x + c

      p.mulAndAdd (1.0).mulAndAdd (alg.element(nn / 2));
   }

   return p;
}


/********************  Polynomial Ring <complex_tag>  ************************/


/**
 *  next()
 */

template<class A>
typename L::PolynomialRingBase<A,L::complex_tag>::type
L::PolynomialRingBase<A,L::complex_tag>::PrimeGenerator::next()
{
   // create polynomial of the form   x + c  with  c  any complex number

   return type(2).mulAndAdd (alg.one()).mulAndAdd (alg.element(n++));
}


/********************  Polynomial Ring <gf_tag>  *****************************/


/**
 *  isPrimitive()
 *
 *  Determine if  p  is primitive, i.e.  x^m = 1 mod p  implies m >= 2^deg - 1
 *
 *  See Lidl/Niederreiter, Finite Fields, Theorem 3.18, and
 *      Knuth, TACP, vol 2, 3.2.2, p.30.
 */

template<class A>
bool L::PolynomialRingBase<A,L::gf_tag>::isPrimitive (const type& p) const
{
   const int deg = p.degree();

   if (deg <= 0 || ! isCanonical (p))  return false;

   // calculate  a0 = (-1)^deg * p[0]

   coeff_type a0 = p[0];
   if (odd (deg))  a.negate (a0);

   // i)  a0  must be a primitive element

   if (! a.isPrimitiveElement (a0))  return false;
   
   // ii) x^r mod p  must be equal to a0

   const unsigned r = (powInt (a.size(), deg) - 1) / (a.size() - 1);
   const type polyx = x();
   const type xPowR
      = powerMod (*static_cast<const PolynomialRing<A>*>(this), polyx, r, p);

   if (! isUnit (xPowR) || xPowR[0] != a0)  return false;

   // iii) For all divisors rr of r, degree of x^rr must be positive

   PrimeDivisors pd (r);

   while (unsigned prime = pd.next())
   {
      if (powerMod (*static_cast<const PolynomialRing<A>*>(this),
                    polyx, r / prime, p).degree() <= 0)
      {
         return false;
      }
   }

   return true;
}


/**
 *  is Prime ()
 *
 *  After checking for a number of simple cases, we use Berlekamp's algorithm
 *  to determine the number of non-trivial factors.
 *
 *  See Knuth, TACP, vol 2, 4.6.2, and
 *      Lidl/Niederreiter, Finite Fields, 4.1.
 */

template<class A>
bool L::PolynomialRingBase<A,L::gf_tag>::isPrime (const type& p) const
{
   const int degree = p.degree();

   if (degree == 1)  return true;       //  bx+a  is always prime
   if (degree <= 0 || a.is0(p[0]))  return false;  // zero, units or ...+0

   // There is a large chance for a linear divisor.
   // So try linear factors first is A is not too large.

   if (a.size() < 50)
   {
      type q (2);
      q.mulAndAdd (a.one()).mulAndAdd (coeff_type());

      for (unsigned i = 1; i < a.size(); ++i)   // x+1, x+2,...
      {
         q[0] = a.element (i);
         if (is0 (rem (p, q)))  return false;
      }

      if (degree <= 3)  return true;
   }

   // if the degree and field size is small, try all possible divisors

   if (degree <=  5 && a.size() <= 11    // max. 11^2 = 121 iterations
    || degree <=  7 && a.size() <= 5     // max.  5^3 = 125 iterations
    || degree <=  9 && a.size() <= 4     // max.  4^4 = 256 iterations
    || degree <= 11 && a.size() <= 5     // max.  3^5 = 243 iterations
    || degree <= 19 && a.size() <= 2)    // max.  2^8 = 256 iterations
   {
      // check all possible divisors with degree >= 2.
      // if  q  is a divisor, so is  k q.  Thus we only check divisors of the
      //    form ...+1.

      // try x^2+1, x^2+x+1, x^2+2x+1,...

      for (unsigned i = sqr(a.size()) + 1; ; i += a.size())
      {
         type q = element (i);
         if (2 * q.degree() > degree)  return true;

         if (is0 (rem (p, q)))  return false;
      }
   }

   // Make sure p is square-free

   if (! isUnit (genGcd (*static_cast<const PolynomialRing<A>*>(this),
                         p, derivative(p))))  return false;

   // Berlekamp's Algorithm

   type factor =
      (int (a.size()) < degree) ? x(a.size()) :
      powerMod (*static_cast<const PolynomialRing<A>*>(this), x(), a.size(), p);

   // The frist row contains only zeros, so we do not include it in the matrix

   Array<coeff_type> matrix ((degree - 1) * degree);
   type q = one();

   // Calculate the matrix  B - I

   for (int row = 0; row < degree - 1; ++row)
   {
      q = rem (mul (q, factor), p);

      for (int col = 0; col < degree; ++col)
      {
         matrix [row * degree + col] = 
               (col <= q.degree()) ? q[col] : coeff_type();
      }

      a.subFrom (matrix [row * degree + row + 1], a.one());
   }

   // degree - rank  gives the number of factors

   int rank = matrixRank (a, matrix.begin(), degree, degree - 1);
   return rank == degree - 1;
}


/**
 *  next()
 */

template<typename A>
typename L::PolynomialRingBase<A,L::gf_tag>::type
L::PolynomialRingBase<A,L::gf_tag>::PrimeGenerator::next()
{
   type p;

   do
   {
      p = alg.elementMonic (n++);
   }
   while (! alg.isPrime (p));
      
   return p;
}


/*********************  Instantiations  **************************************/

/**
 *  Instantiate Polynomial<>
 */

#define HINTLIB_INSTANTIATE_POLYNOMIAL(X) \
   template Polynomial<X >& Polynomial<X >::operator= (const Polynomial<X >&); \
   template Polynomial<X >& Polynomial<X >::divByX (unsigned); \
   template Polynomial<X >& Polynomial<X >::mulByX (unsigned);


/**
 *  Instantiate PolynomialRingBB<>
 */

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_BB(X) \
   template PolynomialRingBB<X >::type \
            PolynomialRingBB<X >::x(unsigned) const; \
   template bool PolynomialRingBB<X >::is1 (const type&) const; \
   template bool PolynomialRingBB<X >::isCanonical (const type&) const; \
   template PolynomialRingBB<X >::type \
            PolynomialRingBB<X >::element(unsigned) const; \
   template PolynomialRingBB<X >::type \
            PolynomialRingBB<X >::elementMonic(unsigned) const; \
   template unsigned PolynomialRingBB<X >::index (const type&) const; \
   template PolynomialRingBB<X >::type& \
            PolynomialRingBB<X >::negate (type&) const; \
   template PolynomialRingBB<X >::type \
            PolynomialRingBB<X >::neg (const type&) const; \
   template PolynomialRingBB<X >::type \
            PolynomialRingBB<X >::add (const type&, const type&) const; \
   template PolynomialRingBB<X >::type \
            PolynomialRingBB<X >::times (const type&, unsigned) const; \
   template PolynomialRingBB<X >::type \
            PolynomialRingBB<X >::derivative (const type&) const; \
   template X::type \
            PolynomialRingBB<X >::evaluate \
                (const type&, const coeff_type&) const; \
   template std::ostream& \
      operator<< (std::ostream &, const PolynomialRingBB<X > &); \
   template void \
            PolynomialRingBB<X >::print (std::ostream &, const type&) const; \
   template void \
            PolynomialRingBB<X >::printShort(std::ostream &, const type&) const;

// Instantiate PolynomialRingBase<ring_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_RING(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_BB(X) \
   template unsigned PolynomialRingBase<X,ring_tag>::additiveOrder (\
         const type& u) const; \
   template PolynomialRingBase<X,ring_tag>::type \
            PolynomialRingBase<X,ring_tag>::mul(const type&, const type&) const;

// Instantiate PolynomialRingBase<domain_tag> and subclasses

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_DOMAIN(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_BB(X) \
   template PolynomialRingBase<X,domain_tag>::type \
            PolynomialRingBase<X,domain_tag>::mul( \
                  const type&, const type&) const; \
   template unsigned \
            PolynomialRingBase<X,domain_tag>::order (const type&) const;


// Instantiate PolynomialRingBase<field_tag> and subclasses

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_FIELD(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_DOMAIN(X) \
   template void PolynomialRingBase<X,field_tag>::div \
        (const type&, const type&, type&, type&) const; \
   template PolynomialRingBase<X,field_tag>::type \
            PolynomialRingBase<X,field_tag>::quot \
               (const type&, const type&) const; \
   template PolynomialRingBase<X,field_tag>::type \
            PolynomialRingBase<X,field_tag>::rem \
               (const type&, const type&) const; \
   template unsigned PolynomialRingBase<X,field_tag>::numOfRemainders \
        (const type&) const; \
   template PolynomialRingBase<X,field_tag>::unit_type \
            PolynomialRingBase<X,field_tag>::makeCanonical (type&) const; \
   template PolynomialRingBase<X,field_tag>::type \
            PolynomialRingBase<X,field_tag>::mulUnit \
               (const type&, const unit_type&) const; \
   template PolynomialRingBase<X,field_tag>::type & \
            PolynomialRingBase<X,field_tag>::mulByUnit \
               (type&, const unit_type&) const; \
   HINTLIB_INSTANTIATE_GENGCD(PolynomialRing<X >)

// Instantiate PolynomialRingBase<rational_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_RATIONAL(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_FIELD(X) \
   template bool PolynomialRingBase<X,rational_tag>::isPrime \
               (const type&) const; \
   template PolynomialRingBase<X,rational_tag>::type \
            PolynomialRingBase<X,rational_tag>::PrimeGenerator::next();

// Instantiate PolynomialRingBase<real_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_REAL(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_FIELD(X) \
   template bool PolynomialRingBase<X,real_tag>::isPrime (const type&) const; \
   template PolynomialRingBase<X,real_tag>::type \
            PolynomialRingBase<X,real_tag>::PrimeGenerator::next();

// Instantiate PolynomialRingBase<complex_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_COMPLEX(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_FIELD(X) \
   template PolynomialRingBase<X,real_tag>::type \
            PolynomialRingBase<X,real_tag>::PrimeGenerator::next();

// Instantiate PolynomialRingBase<gf_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_GF(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_FIELD(X) \
   template bool PolynomialRingBase<X,gf_tag>::isPrimitive (const type&) const;\
   template bool PolynomialRingBase<X,gf_tag>::isPrime (const type&) const; \
   template PolynomialRingBase<X,gf_tag>::type \
            PolynomialRingBase<X,gf_tag>::PrimeGenerator::next();

