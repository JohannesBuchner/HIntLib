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

#ifdef __GNUG__
#pragma implementation
#endif

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

#include <HIntLib/integerring.h>
#include <HIntLib/modulararithmetic.h>
#include <HIntLib/precalculatedfield.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;


/********************** Polynomial <> ****************************************/

/**
 *  Constructor
 */

template<class T>
L::Polynomial<T>::Polynomial (unsigned size, const T* coef)
   : c (coef, coef+size)
{}


/**
 *  Assignment
 */

template<class T>
L::Polynomial<T>& L::Polynomial<T>::operator= (const L::Polynomial<T>& p)
{
   if (this != &p) c = p.c;
   return *this;
}


/**
 *  divByXPow()
 *  mulByXPow()
 */

template<class A>
L::Polynomial<A>& L::Polynomial<A>::divByXPow (unsigned k)
{
   for (unsigned i = 0; i < k; ++i)  divByX();
   return *this;
}

template<class A>
L::Polynomial<A>& L::Polynomial<A>::mulByXPow (unsigned k)
{
   for (unsigned i = 0; i < k; ++i)  mulByX();
   return *this;
}


/**
 *  Prints a polynomial to an output stream.
 */

template<class A>
std::ostream& L::operator<<
   (std::ostream& o, const L::Polynomial<A>& p)
{
   std::ostringstream ss;

   bool output = false;

   for (int i = p.degree(); i >= 0; --i)
   {
      if (output) ss << '+';
   
      ss << p[i];

      switch (i)
      {
      case 0:  break;
      case 1:  ss << 'x'; break;
      case 2:  ss << "x\262"; break;
      case 3:  ss << "x\263"; break;
      default: ss << "x^" << i;
      }

      output = true;
   }

   if (! output)  ss << '0';

   o << ss.str().c_str();

   return o;
}


/*******************  Polynomail Ring <>  ************************************/

/**
 *  one()
 */

template<class A>
typename L::PolynomialRing<A>::type L::PolynomialRing<A>::one() const
{
   P p (1);
   p.mulAndAdd (a.one());
   return p;
}


/**
 *  is1 ()
 */

template<class A>
bool L::PolynomialRing<A>::is1 (const P &p) const
{
   return p.degree() == 0 && a.is1(p.lc());
}


/**
 *  isMonic ()
 */

template<class A>
bool L::PolynomialRing<A>::isMonic (const P &p) const
{
   return p.degree() >= 0 && a.is1(p.lc());
}


/**
 *  element()
 */

template<class A>
typename L::PolynomialRing<A>::type
L::PolynomialRing<A>::element(unsigned i) const
{
   unsigned s = a.size()  ?  a.size()  : 5;

   unsigned ii = i;
   unsigned num = 0;
   while (ii)
   {
      ++num;
      ii /= s;
   }
   
   P p (num);
   p.mulByXPow (num);

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

template<class A>
unsigned L::PolynomialRing<A>::index (const P &p) const
{
   unsigned s = a.size()  ?  a.size()  : 5;

   unsigned n = 0;
   for (int i = p.degree(); i >= 0; --i)  n = s * n + a.index (p[i]);
   return n;
}


/**
 *  negate()
 */

template<class A>
typename L::PolynomialRing<A>::type & L::PolynomialRing<A>::negate (P& p) const
{ 
   for (int i = p.degree(); i >= 0; --i)  a.negate (p[i]);
   return p;
}


/**
 *  neg()
 */

template<class A>
typename L::PolynomialRing<A>::type L::PolynomialRing<A>::neg (const P& p) const
{
   P n (p.degree() + 1);
   for (int i = p.degree(); i >= 0; --i)  n.mulAndAdd (a.neg (p[i]));
   return n;
}


/**
 *  add()
 */

template<class A>
typename L::PolynomialRing<A>::type
L::PolynomialRing<A>::add (const P& p1, const P& p2) const
{
   const int deg1 = p1.degree();
   const int deg2 = p2.degree();
   const int maxDeg = std::max (deg1, deg2); 

   P p (maxDeg + 1);
   int start;
   
   if (deg1 != deg2)
   {
      p.mulByXPow (maxDeg + 1);
      start = std::min (deg1, deg2);
      const P& greater = deg1 > deg2 ? p1 : p2; 

      for (int i = maxDeg; i > start; --i)  p[i] = greater[i];
   }
   else
   {
      start = deg1;

      while (start >= 0 && a.is0 (a.add (p1[start], p2[start])))  --start;

      p.mulByXPow (start + 1);
   }

   for (int i = start; i >= 0; --i)  p[i] = a.add (p1[i], p2[i]);

   return p;
}


/**
 *  mul()
 */

template<class A>
typename L::PolynomialRing<A>::type
L::PolynomialRing<A>::mul (const P& _p1, const P& _p2) const
{
   const P& p1 = _p1.degree() > _p2.degree()  ?  _p1 : _p2;
   const P& p2 = _p1.degree() > _p2.degree()  ?  _p2 : _p1;
   
   const int deg1 = p1.degree();
   const int deg2 = p2.degree();

   if (deg1 < 0 || deg2 < 0)  return P(0);

   P p (deg1 + deg2 + 1);
   bool nonZero = false;

   for (int k = deg1 + deg2; k >= 0; --k)
   {
      T c = a.zero();

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


/**
 *  times()
 */

template<class A>
typename L::PolynomialRing<A>::type
L::PolynomialRing<A>::times (const P& p, unsigned k) const
{
   P n (p.degree() + 1);
   bool nonZero = false;

   for (int i = p.degree(); i >= 0; --i)
   {
      T c = a.times (p[i], k);

      if (nonZero || ! a.is0 (c))
      {
         nonZero = true;
         n.mulAndAdd (c);
      }
   }
   return n;
}


/**
 *  operator<<
 */

template<class A>
std::ostream&
L::operator<< (std::ostream &o, const PolynomialRing<A> &a)
{
   return o << a.arithmetic() << "[x]";
}


/********************  Polynomial Ring Field <> ******************************/


/**
 *  div
 *
 *  Division of polynomials
 *
 *  See D.E.Knuth, Art o Computer Programming, 4.6.1, Algo D
 */

template<class A>
void L::PolynomialRingField<A>::div
   (const P& u, const P& v, P& _q, P& _r) const
{
   const int vdeg = v.degree();
   if (vdeg == -1)  throw DivisionByZero();
   const int udeg = u.degree();

   if (udeg < vdeg)
   {
      _q = zero();
      _r = u;
      return;
   }

   T qq = a.recip (v.lc());
   P uu (u);
   P q (udeg - vdeg + 1);

   for (int k = udeg-vdeg; k >= 0; --k)
   {
      T lc = a.mul (uu[vdeg+k], qq);
      q.mulAndAdd (lc);

      for (int j = vdeg + k - 1; j >= k; --j)
      {
         a.subFrom (uu[j], a.mul (lc, v[j-k]));
      }
   }

   _q = q;
   
   bool nonZero = false;
   P r (vdeg);
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
 * is Prime()
 */

template<class A>
bool L::PolynomialRingField<A>::isPrime (const P& p) const
{
   if (p.degree() <= 0)  return false;
   if (p.degree() == 1)  return true;

   if (! a.size())  throw InternalError (__FILE__, __LINE__);

   for (unsigned i = a.size(); ; ++i)
   {
      P q = element (i);

      if (2 * q.degree() > p.degree())  return true;

      P x, r;
      div (p, q, x, r);
      if (is0(r))  return false;
   }
}


/**
 *  num Of Remainders()
 */

template<class A>
unsigned L::PolynomialRingField<A>::numOfRemainders (const P& p) const
{
   return a.size()  ?  powInt (a.size(), p.degree()) : 0;
}

namespace HIntLib
{
   // Instantiate Polynomial<>

#define HINTLIB_INSTANTIATE(X) \
   template Polynomial<X>::Polynomial (unsigned, const X*); \
   template Polynomial<X>& Polynomial<X>::operator= (const Polynomial<X>&); \
   template Polynomial<X>& Polynomial<X>::divByXPow (unsigned); \
   template Polynomial<X>& Polynomial<X>::mulByXPow (unsigned); \
   template std::ostream& operator<< (std::ostream &, const Polynomial<X> &);

   // template class Polynomial<X>;

   HINTLIB_INSTANTIATE (int)
   HINTLIB_INSTANTIATE (unsigned char)
   HINTLIB_INSTANTIATE (unsigned short)
#undef HINTLIB_INSTANTIATE

   // Instantiate PolynomialRing<>

#define HINTLIB_INSTANTIATE(X) \
   template PolynomialRing<X>::type PolynomialRing<X>::one() const; \
   template bool PolynomialRing<X>::is1 (const P&) const; \
   template bool PolynomialRing<X>::isMonic (const P&) const; \
   template PolynomialRing<X>::type \
            PolynomialRing<X>::element(unsigned) const; \
   template unsigned PolynomialRing<X>::index (const P&) const; \
   template PolynomialRing<X>::type& PolynomialRing<X>::negate (P&) const; \
   template PolynomialRing<X>::type PolynomialRing<X>::neg (const P&) const; \
   template PolynomialRing<X>::type \
            PolynomialRing<X>::add (const P&, const P&) const; \
   template PolynomialRing<X>::type \
            PolynomialRing<X>::mul (const P&, const P&) const; \
   template PolynomialRing<X>::type \
            PolynomialRing<X>::times (const P&, unsigned) const; \
   template std::ostream& \
      operator<< (std::ostream &, const PolynomialRing<X> &);

   HINTLIB_INSTANTIATE (IntegerRing<int>)
   HINTLIB_INSTANTIATE (ModularIntegerRing<unsigned char>)
   HINTLIB_INSTANTIATE (ModularIntegerRing<unsigned short>)
   HINTLIB_INSTANTIATE (ModularIntegerField<unsigned char>)
   HINTLIB_INSTANTIATE (ModularIntegerField<unsigned short>)
   HINTLIB_INSTANTIATE (PrecalculatedField<unsigned char>)
   HINTLIB_INSTANTIATE (PrecalculatedField<unsigned short>)
#undef HINTLIB_INSTANTIATE

   // Instantiate PolynomialRingField<>

#define HINTLIB_INSTANTIATE(X) \
   template void PolynomialRingField<X>::div \
      (const P&, const P&, P&, P&) const; \
   template bool PolynomialRingField<X>::isPrime (const P&) const; \
   template unsigned PolynomialRingField<X>::numOfRemainders (const P&) const;

   HINTLIB_INSTANTIATE (ModularIntegerField<unsigned char>)
   HINTLIB_INSTANTIATE (ModularIntegerField<unsigned short>)
   HINTLIB_INSTANTIATE (PrecalculatedField<unsigned char>)
   HINTLIB_INSTANTIATE (PrecalculatedField<unsigned short>)
#undef HINTLIB_INSTANTIATE
}

