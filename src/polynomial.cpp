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
#include <HIntLib/lookupfield.h>
#include <HIntLib/exception.h>
#include <HIntLib/polynomial2.h>

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
 *  divByXPow()
 *  mulByXPow()
 */

template<typename A>
L::Polynomial<A>&
L::Polynomial<A>::divByXPow (unsigned k)
{
   for (unsigned i = 0; i < k; ++i)  divByX();
   return *this;
}

template<typename A>
L::Polynomial<A>&
L::Polynomial<A>::mulByXPow (unsigned k)
{
   for (unsigned i = 0; i < k; ++i)  mulByX();
   return *this;
}


/**
 *  Prints a polynomial to an output stream.
 */

template<typename A>
std::ostream&
L::operator<< (std::ostream& o, const L::Polynomial<A>& p)
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

   return o << ss.str().c_str();
}


/*******************  Polynomail Ring <>  ************************************/

/**
 *  one()
 */

template<typename A>
typename L::PolynomialRingBB<A>::type
L::PolynomialRingBB<A>::one() const
{
   type p (1);
   p.mulAndAdd (a.one());
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
 *  isMonic ()
 */

template<class A>
bool
L::PolynomialRingBB<A>::isMonic (const type &p) const
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
   unsigned s = a.size()  ?  a.size()  : 5;

   unsigned ii = i;
   unsigned num = 0;
   while (ii)
   {
      ++num;
      ii /= s;
   }
   
   type p (num);
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
      p.mulByXPow (maxDeg + 1);
      start = std::min (deg1, deg2);
      const type& greater = deg1 > deg2 ? p1 : p2; 

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

template<typename A>
typename L::PolynomialRingBB<A>::type
L::PolynomialRingBB<A>::mul (const type& _p1, const type& _p2) const
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
      typename A::type c = a.zero();

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

template<typename A>
typename L::PolynomialRingBB<A>::type
L::PolynomialRingBB<A>::times (const type& p, unsigned k) const
{
   type n (p.degree() + 1);
   bool nonZero = false;

   for (int i = p.degree(); i >= 0; --i)
   {
      typename A::type c = a.times (p[i], k);

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
   if (p.degree() < 0)  return a.zero();

   coeff_type xx = x;
   coeff_type res = p[0];

   for (int i = 1; i <= p.degree(); ++i)
   {
      a.addTo (res, a.mul(x, p[i]));
      a.mulBy (xx, x);
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

   ss << a.getCoeffAlgebra() << "[x]";

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
         case 1:  ss << 'x'; break;
         case 2:  ss << "x\262"; break;
         case 3:  ss << "x\263"; break;
         default: ss << "x^" << i;
      }

      output = true;
   }

   if (! output)  a.printShort (ss, a.zero());

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


/********************  Polynomial Ring Field <> ******************************/


/**
 *  div
 *
 *  Division of polynomials
 *
 *  See D.E.Knuth, Art o Computer Programming, 4.6.1, Algo D
 */

template<class A>
void L::PolynomialRingBase<A,L::field_tag>::div
   (const type& u, const type& v, type& _q, type& _r) const
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

   typename A::type qq = a.recip (v.lc());
   type uu (u);
   type q (udeg - vdeg + 1);

   for (int k = udeg-vdeg; k >= 0; --k)
   {
      typename A::type lc = a.mul (uu[vdeg+k], qq);
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
   if (vdeg == -1)  throw DivisionByZero();
   const int udeg = u.degree();

   if (udeg < vdeg)  return type(0);

   typename A::type qq = a.recip (v.lc());
   type uu (u);
   type q (udeg - vdeg + 1);

   for (int k = udeg-vdeg; k >= 0; --k)
   {
      typename A::type lc = a.mul (uu[vdeg+k], qq);
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
   if (vdeg == -1)  throw DivisionByZero();
   const int udeg = u.degree();

   if (udeg < vdeg)  return u;

   typename A::type qq = a.recip (v.lc());
   type uu (u);

   for (int k = udeg-vdeg; k >= 0; --k)
   {
      typename A::type lc = a.mul (uu[vdeg+k], qq);

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
 * is Prime()
 */

template<class A>
bool L::PolynomialRingBase<A,L::field_tag>::isPrime (const type& p) const
{
   if (p.degree() == 1)  return true;       //  bx+a  is always prime
   if (p.degree() <= 0 || a.is0(p[0]))  return false;  // 0, a  or ...+0

   if (! a.size())  throw InternalError (__FILE__, __LINE__);

   // check all possible divisors starting with  x+1.
   // No need to check x, because we got rid of this case already.
   // if  q  is a divisor, so is  k q.  Thus we only check divisors with ...+1.

   for (unsigned i = a.size() + 1; ; i += a.size())   // x+1, 2x+1,...
   {
      type q = element (i);
      if (2 * q.degree() > p.degree())  return true;

      if (is0 (rem (p, q)))  return false;
   }
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

namespace HIntLib
{
   // Instantiate Polynomial<>

#define HINTLIB_INSTANTIATE(X) \
   template Polynomial<X>& Polynomial<X>::operator= (const Polynomial<X>&); \
   template Polynomial<X>& Polynomial<X>::divByXPow (unsigned); \
   template Polynomial<X>& Polynomial<X>::mulByXPow (unsigned); \
   template std::ostream& operator<< (std::ostream &, const Polynomial<X> &);

   HINTLIB_INSTANTIATE (int)
   HINTLIB_INSTANTIATE (unsigned char)
   HINTLIB_INSTANTIATE (unsigned short)
#undef HINTLIB_INSTANTIATE

   // Instantiate PolynomialRingBB<>

#define HINTLIB_INSTANTIATE(X) \
   template PolynomialRingBB<X >::type \
            PolynomialRingBB<X >::one() const; \
   template bool PolynomialRingBB<X >::is1 (const type&) const; \
   template bool PolynomialRingBB<X >::isMonic (const type&) const; \
   template PolynomialRingBB<X >::type \
            PolynomialRingBB<X >::element(unsigned) const; \
   template unsigned PolynomialRingBB<X >::index (const type&) const; \
   template PolynomialRingBB<X >::type& \
            PolynomialRingBB<X >::negate (type&) const; \
   template PolynomialRingBB<X >::type \
            PolynomialRingBB<X >::neg (const type&) const; \
   template PolynomialRingBB<X >::type \
            PolynomialRingBB<X >::add (const type&, const type&) const; \
   template PolynomialRingBB<X >::type \
            PolynomialRingBB<X >::mul (const type&, const type&) const; \
   template PolynomialRingBB<X >::type \
            PolynomialRingBB<X >::times (const type&, unsigned) const; \
   template X::type \
            PolynomialRingBB<X >::evaluate \
                (const type&, const coeff_type&) const; \
   template std::ostream& \
      operator<< (std::ostream &, const PolynomialRingBB<X > &); \
   template void \
            PolynomialRingBB<X >::print (std::ostream &, const type&) const; \
   template void \
            PolynomialRingBB<X >::printShort(std::ostream &, const type&) const;

   HINTLIB_INSTANTIATE (IntegerRing<int>)
   HINTLIB_INSTANTIATE (GF2)
   HINTLIB_INSTANTIATE (ModularArith<unsigned char>)
   HINTLIB_INSTANTIATE (ModularArith<unsigned short>)
   HINTLIB_INSTANTIATE (ModularArithField<unsigned char>)
   HINTLIB_INSTANTIATE (ModularArithField<unsigned short>)
   HINTLIB_INSTANTIATE (LookupField<unsigned char>)
   HINTLIB_INSTANTIATE (LookupFieldPow2<unsigned char>)
   HINTLIB_INSTANTIATE (LookupFieldPrime<unsigned char>)
#undef HINTLIB_INSTANTIATE

   // Instantiate PolynomialRingBase<field_tag>

#define HINTLIB_INSTANTIATE(X) \
   template void PolynomialRingBase<X,field_tag>::div \
      (const type&, const type&, type&, type&) const; \
   template PolynomialRingBase<X,field_tag>::type \
            PolynomialRingBase<X,field_tag>::quot \
               (const type&, const type&) const; \
   template PolynomialRingBase<X,field_tag>::type \
            PolynomialRingBase<X,field_tag>::rem \
               (const type&, const type&) const; \
   template bool PolynomialRingBase<X,field_tag>::isPrime (const type&) const; \
   template unsigned PolynomialRingBase<X,field_tag>::numOfRemainders (const type&) const;

   HINTLIB_INSTANTIATE (GF2)
   HINTLIB_INSTANTIATE (ModularArithField<unsigned char>)
   HINTLIB_INSTANTIATE (ModularArithField<unsigned short>)
   HINTLIB_INSTANTIATE (LookupField<unsigned char>)
   HINTLIB_INSTANTIATE (LookupFieldPow2<unsigned char>)
   HINTLIB_INSTANTIATE (LookupFieldPrime<unsigned char>)
#undef HINTLIB_INSTANTIATE
}

