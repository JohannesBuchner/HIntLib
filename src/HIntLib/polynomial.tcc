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

#ifndef HINTLIB_POLYNOMIAL_TCC
#define HINTLIB_POLYNOMIAL_TCC 1

#include <iostream>

#include <HIntLib/polynomial.h>

#include <HIntLib/output.h>
#include <HIntLib/exception.h>

/********************** Polynomial <> ****************************************/

namespace HIntLib
{

/**
 *  Copy constructor
 */

template<typename T>
Polynomial<T>::Polynomial (const P& p)
   : c (p.c)
{}


/**
 *  Destructor
 *
 *  Not inline because it includes the destruction of the vector<T>
 */

template<typename T>
Polynomial<T>::~Polynomial() {}


/**
 *  divByX()
 *  mulByX()
 */

template<typename A>
Polynomial<A>&
Polynomial<A>::divByX (unsigned k)
{
   for (unsigned i = 0; i < k; ++i)  divByX();
   return *this;
}

template<typename A>
Polynomial<A>&
Polynomial<A>::mulByX (unsigned k)
{
   for (unsigned i = 0; i < k; ++i)  mulByX();
   return *this;
}


/**
 *  operator==()
 */

// for equal()-bug

#ifdef HINTLIB_EQUAL_BUG
template<typename T>
bool
operator== (const Polynomial<T>& p1, const Polynomial<T>& p2)
{
   if (p1.numCoeff() != p2.numCoeff())  return false;

   for (typename Polynomial<T>::CDownI i1 = p1.fromLc(), i2 = p2.fromLc();
        i1 != p1.toA0(); ++i1, ++i2)
   {
      if (! (*i1 == *i2)) return false;
   }

   return true;
}
#endif

} // namespace HIntLib


/*******************  Polynomial Ring Base -- A  *****************************/

/**
 *  x ()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::Xn,A,T> x)
{
   const unsigned size = x.n;
   c.reserve (size);
   mulAndAdd (x.a->one());
   c.resize (size);
}

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::X,A,T> x)
{
   c.reserve (2);
   mulAndAdd (x.a->one());
   mulAndAdd (T());
}


/**
 *  index()
 */

template<typename A>
unsigned
HIntLib::Private::PRBA<A>::indexImp (const type &p, finite_tag) const
{
   const unsigned s = a.size();
   unsigned n = 0;
   for (int i = p.degree(); i >= 0; --i)  n = s * n + a.index (p[i]);
   return n;
}


template<typename A>
unsigned
HIntLib::Private::PRBA<A>::indexImp (const type &p, infinite_tag) const
{
   unsigned indices [10];   // sqrt(digits) is sufficient
   int num = p.numCoeff();
   if (num > 10)  num = 10;

   for (int i = 0; i < num; ++i)  indices [i] = a.index (p[i]);

   return threadinf (indices, num);
}


/**
 *  element()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::
Polynomial (Private::PG<Private::Element<finite_tag>,A,T> x)
{
   unsigned i = x.n;
   if (! i)  return;

   const A* a = x.a;
   const unsigned s = a->size();

   c.resize (logInt (i, s) + 1);

   for (unsigned j = 0; i; )
   {
      (*this)[j++] = a->element (i % s);
      i /= s;
   }
}

template<typename T> template<typename A>
HIntLib::Polynomial<T>::
Polynomial (Private::PG<Private::Element<infinite_tag>,A,T> x)
{
   unsigned i = x.n;
   const A* a = x.a;

   unsigned indices [10];  // sqrt(digits) is sufficient
   int num = unthreadinf (i, indices);
   reserve (num);
   for (int i = num - 1; i >= 0; --i)  mulAndAdd (a->element(indices[i]));
}


/**
 *  elementMonic ()
 *
 *  Enumerate monic polynomials with non-zero constant term
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::
Polynomial (Private::PG<Private::ElementMonic<finite_tag>,A,T> x)
{
   unsigned ind = x.n;
   const A* a = x.a;

   // elementMonic(0) = 1

   if (!ind)
   {
      mulAndAdd (a->one());
      return;
   }

   const unsigned s = a->size();

   // determine the number of coefficients

   --ind;
   unsigned max = s;
   unsigned numCoeffs = 1;

   while (ind >= max)
   {
      ind -= max;
      max *= s;
      ++numCoeffs;
   }

   // set leading coefficient

   c.reserve (numCoeffs + 1);
   mulAndAdd (a->one());

   // set remaining coefficients

   for (unsigned i = 0; i < numCoeffs; ++i)
   {
      max /= s;
      mulAndAdd (a->element (ind / max));
      ind %= max;
   }
}

template<typename T> template<typename A>
HIntLib::Polynomial<T>::
Polynomial (Private::PG<Private::ElementMonic<infinite_tag>,A,T> x)
{
   unsigned ind = x.n;
   const A* a = x.a;

   // elementMonic(0) = 1

   if (!ind)
   {
      mulAndAdd (a->one());
      return;
   }

   unsigned indices [10];  // sqrt(digits) is sufficient
   int num = unthreadinf (ind, indices);

   reserve (num + 1);
   mulAndAdd (a->one());
   mulAndAdd (a->element(indices[num - 1] - 1));  // num > 0, since i != 0
   for (int i = num - 2; i >= 0; --i)  mulAndAdd (a->element(indices[i]));
}


/**
 *  is1 ()
 */

template<typename A>
bool
HIntLib::Private::PRBA<A>::is1 (const type &p) const
{
   return p.numCoeff() == 1 && a.is1(p.lc());
}


/**
 *  isMonic ()
 */

template<class A>
bool
HIntLib::Private::PRBA<A>::isMonic (const type &p) const
{
   return p.degree() < 0 || a.is1(p.lc());
}


/**
 *  negate()
 */

template<typename A>
void
HIntLib::Private::PRBA<A>::negateImp (type& p, char_any) const
{ 
   const typename type::DownI end = p.toA0();
   for (typename type::DownI i = p.fromLc(); i != end; ++i)  a.negate(*i);
}


/**
 *  neg()
 */

template<typename T> template<typename A, typename X>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::Neg<X>,A,T> x)
{
   const A* a = x.a;
   const P* p = x.p;

   reserve (p->numCoeff());

   const CDownI end = p->toA0();
   for (CDownI i = p->fromLc(); i != end; ++i)  mulAndAdd (a->neg (*i));
}

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::Neg<char_two>,A,T> x)
   : c (x.p->c) {}


/**
 *  dbl()
 */

// If we don't have a characteristic

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::Dbl<char_none>,A,T> x)
{
   const A* a = x.a;
   const P* p = x.p;

   const CDownI end = p->toA0();
         CDownI i = p->fromLc();

   while (i != end)
   {
       coeff_type c = a->dbl (*i++);

       if (! a->is0 (c))
       {
          reserve (end - i + 1);
          mulAndAdd (c);
          break;
       }
   }

   while (i != end)  mulAndAdd (a->dbl (*i++));
}

// char = 0

template<typename T> template<typename A>
HIntLib::Polynomial<T>::
Polynomial (Private::PG<Private::Dbl<char_zero>,A,T> x)
{
   const A* a = x.a;
   const P* p = x.p;

   const CDownI end = x.p->toA0();
         CDownI i   = x.p->fromLc();

   reserve (end - i);
   while (i != end)  mulAndAdd (a->dbl (*i++));
}

// char finite

template<typename T> template<typename A>
HIntLib::Polynomial<T>::
Polynomial (Private::PG<Private::Dbl<char_prime>,A,T> x)
{
   const A* a = x.a;
   if (a->characteristic() == 2)  return;
   const P* p = x.p;

   const CDownI end = x.p->toA0();
         CDownI i   = x.p->fromLc();

   reserve (end - i);
   while (i != end)  mulAndAdd (a->dbl (*i++));
}

// char = 2 is handled inline


/**
 *  times2()
 */

// ring with zero divisors

template<typename A>
void
HIntLib::Private::PRBA<A>::times2Imp (type& p, char_none) const
{
   const typename type::DownI end = p.toA0();
         typename type::DownI i = p.fromLc();
   unsigned numToErase = 0;

   // start doubeling entries. Check for non-zero result

   while (i != end)
   {
       a.times2 (*i);

       if (! a.is0 (*i))
       {
          numToErase = i - p.fromLc();
          break;
       }
       ++i;
   }

   if (i != end)  // we found a non-zero element
   {
      ++i;
      while (i != end)  a.times2 (*i++);  // double the rest
      p.erase (numToErase);    // remove 0 coefficients
   }
   else
   {
      p.makeZero();
   }
}

// char = prime

template<typename A>
void
HIntLib::Private::PRBA<A>::times2Imp (type& p, char_prime) const
{
   if (a.characteristic() == 2)
   {
      p.makeZero();
   }
   else
   {
      const typename type::DownI end = p.toA0();
            typename type::DownI i   = p.fromLc();
      while (i != end)  a.times2 (*i++);
   }
}

// char = 0

template<typename A>
void
HIntLib::Private::PRBA<A>::times2Imp (type& p, char_zero) const
{
   const typename type::DownI end = p.toA0();
         typename type::DownI i   = p.fromLc();
   while (i != end)  a.times2 (*i++);
}

// char = 2 is inline


/**
 *  times()
 */

// characteristic does not exit

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (
      Private::PG<Private::Times<char_none>,A,T> x)
{
   unsigned n = x.n;
   if (!n)  return;

   const A* a = x.a;
   const P* p = x.p;

   const CDownI end = p->toA0();
         CDownI i = p->fromLc();

   while (i != end)
   {
       coeff_type c = a->times (*i++, n);

       if (! a->is0 (c))
       {
          reserve (end - i + 1);
          mulAndAdd (c);
          break;
       }
   }

   while (i != end)  mulAndAdd (a->times (*i++, n));
}

// char = 0, no cancelation can occure

template<typename T> template<typename A>
HIntLib::Polynomial<T>::
Polynomial (Private::PG<Private::Times<char_zero>,A,T> x)
{
   unsigned n = x.n;
   if (!n)  return;

   const A* a = x.a;

   const CDownI end = x.p->toA0();
         CDownI i   = x.p->fromLc();

   reserve (end - i);
   while (i != end)  mulAndAdd (a->times (*i++, n));
}

// char prime

template<typename T> template<typename A>
HIntLib::Polynomial<T>::
Polynomial (Private::PG<Private::Times<char_prime>,A,T> x)
{
   const A* a = x.a;
   unsigned ca = a->characteristic();
   unsigned n = x.n;

   if (ca && n >= ca)  n %= ca;  // if we can reduce n, do so

   if (!n)  return;   // 0 x = 0

   if (n == 1)   // 1 x = x
   {
      c = x.p->c;
      return;
   }

   const CDownI end = x.p->toA0();
         CDownI i   = x.p->fromLc();

   while (i != end)
   {
       coeff_type coeff = a->times (*i++, n);

       if (! a->is0 (coeff))
       {
          reserve (end - i + 1);
          mulAndAdd (coeff);
          break;
       }
   }

   while (i != end)  mulAndAdd (a->times (*i++, n));
}

// char = 2

template<typename T> template<typename A>
HIntLib::Polynomial<T>::
Polynomial (Private::PG<Private::Times<char_two>,A,T> x)
{
   if (x.n & 1)  c = x.p->c;
}



/**
 *  add()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::Add,A,T> x)
{
   const A* a = x.a;
   const P* p1;
   const P* p2;
   
   if (x.p1->numCoeff() > x.p2->numCoeff()) // the larger poly is p1
   {
      p1 = x.p1;
      p2 = x.p2;
   }
   else
   {
      p1 = x.p2;
      p2 = x.p1;
   }

   const int num1 = p1->numCoeff();
   const int num2 = p2->numCoeff();

   const CDownI end1 = p1->toA0();
         CDownI i1   = p1->fromLc();
         CDownI i2   = p2->fromLc();
   
   if (num1 != num2)
   {
      reserve (num1);
      const CDownI oldI1 = i1;
      std::copy (oldI1, i1 = end1 - num2, back_inserter(c));
   }
   else
   {
      while (i1 != end1)
      {
         coeff_type coeff = a->add (*i1++, *i2++);

         if (! a->is0 (coeff))
         {
            reserve ((end1 - i1) + 1);
            mulAndAdd (coeff);
            break;
         }
      }
   }

   while (i1 != end1)  mulAndAdd (a->add (*i1++, *i2++));
}


/**
 *  sub()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::Sub,A,T> x)
{
   const A* a = x.a;
   const P* p1 = x.p1;
   const P* p2 = x.p2;

   const int num1 = p1->numCoeff();
   const int num2 = p2->numCoeff();

   const CDownI end1 = p1->toA0();
         CDownI i1   = p1->fromLc();
         CDownI i2   = p2->fromLc();

   if (num1 > num2)   // p1 longer: copy extra coefficients
   {
      reserve (num1);
      const CDownI oldI1 = i1;
      std::copy (oldI1, i1 = end1 - num2, back_inserter(c));
   }
   else if (num1 < num2)   // p2 longer: copy negatives of extra coefficients
   {
      reserve (num2);

      const CDownI end = p2->toA0() - num1;
      while (i2 != end)  mulAndAdd (a->neg (*i2++));
   }
   else   // deg1 == deg2:  find first non-zero result
   {
      while (i1 != end1)
      {
         coeff_type coeff = a->sub (*i1++, *i2++);

         if (! a->is0 (coeff))
         {
            reserve ((end1 - i1) + 1);
            mulAndAdd (coeff);
            break;
         }
      }
   }

   while (i1 != end1)  mulAndAdd (a->sub (*i1++, *i2++));
}


/**
 *  addTo()
 */

template<typename A>
void
HIntLib::Private::PRBA<A>::addTo (type& p1, const type& p2) const
{ 
   const int num1 = p1.numCoeff();
   const int num2 = p2.numCoeff();

   if (num1 > num2)   // target is larger
   {
      const typename type:: DownI end1 = p1.toA0();
            typename type:: DownI i1   = p1.fromLc() + (num1 - num2);
            typename type::CDownI i2   = p2.fromLc();

      while (i1 != end1)  a.addTo (*i1++, *i2++);
   }
   else if (num1 == num2)  // equal size
   {
      const typename type:: DownI end1 = p1.toA0();
            typename type:: DownI i1   = p1.fromLc();
            typename type::CDownI i2   = p2.fromLc();

      unsigned numToErase = 0;

      while (i1 != end1)
      {
          a.addTo (*i1, *i2++);

          if (! a.is0 (*i1))
          {
             numToErase = i1 - p1.fromLc();
             break;
          }
          ++i1;
      }

      if (i1 != end1)  // we found a non-zero element
      {
         ++i1;
         while (i1 != end1)  a.addTo (*i1++, *i2++);  // add the rest
         p1.erase (numToErase);    // remove 0 coefficients
      }
      else
      {
         p1.makeZero();
      }
   }
   else  // num2 > num1
   {
      type pp (Private::PG<Private::Add,A> (a, p1, p2));
      p1.swap (pp);
   }
}


/**
 *  subFrom()
 */

template<typename A>
void
HIntLib::Private::PRBA<A>::subFrom (type& p1, const type& p2) const
{ 
   const int num1 = p1.numCoeff();
   const int num2 = p2.numCoeff();

   if (num1 > num2)   // target is larger
   {
      const typename type:: DownI end1 = p1.toA0();
            typename type:: DownI i1   = p1.fromLc() + (num1 - num2);
            typename type::CDownI i2   = p2.fromLc();

      while (i1 != end1)  a.subFrom (*i1++, *i2++);
   }
   else if (num1 == num2)  // equal size
   {
      const typename type:: DownI end1 = p1.toA0();
            typename type:: DownI i1   = p1.fromLc();
            typename type::CDownI i2   = p2.fromLc();

      unsigned numToErase = 0;

      while (i1 != end1)
      {
          a.subFrom (*i1, *i2++);

          if (! a.is0 (*i1))
          {
             numToErase = i1 - p1.fromLc();
             break;
          }
          ++i1;
      }

      if (i1 != end1)  // we found a non-zero element
      {
         ++i1;
         while (i1 != end1)  a.subFrom (*i1++, *i2++);  // sub the rest
         p1.erase (numToErase);    // remove 0 coefficients
      }
      else
      {
         p1.makeZero();
      }
   }
   else  // num2 > num1
   {
      type pp (Private::PG<Private::Sub,A> (a, p1, p2));
      p1.swap (pp);
   }
}


/**
 *  mul()
 */

// For polynomials over domains, the leading coefficients cannot cancel out

template<typename T> template<typename A>
HIntLib::Polynomial<T>::
Polynomial (Private::PG<Private::Mul<nozerodivisor_tag>,A,T> x)
{
   const A* a = x.a;
   const P* p1 = x.p1->degree() > x.p2->degree()  ?  x.p1 : x.p2;
   const P* p2 = x.p1->degree() > x.p2->degree()  ?  x.p2 : x.p1;
   
   const int deg1 = p1->degree();
   const int deg2 = p2->degree();

   if (deg1 < 0 || deg2 < 0)  return;

   reserve (deg1 + deg2 + 1);

   for (int k = deg1 + deg2; k >= 0; --k)
   {
      coeff_type c = coeff_type();

      for (int i = std::max (k - deg2, 0); i <= std::min (k, deg1); ++i)
      {
         a->addTo (c, a->mul ((*p1)[i], (*p2)[k - i]));
      }

      mulAndAdd (c);
   }
}

// if there are zero divisors, watch out for cancelation!

template<typename T> template<typename A>
HIntLib::Polynomial<T>::
Polynomial (Private::PG<Private::Mul<zerodivisor_tag>,A,T> x)
{
   const A* a = x.a;
   const P* p1 = x.p1->degree() > x.p2->degree()  ?  x.p1 : x.p2;
   const P* p2 = x.p1->degree() > x.p2->degree()  ?  x.p2 : x.p1;
   
   const int deg1 = p1->degree();
   const int deg2 = p2->degree();

   if (deg1 < 0 || deg2 < 0)  return;

   reserve (deg1 + deg2 + 1);
   bool nonZero = false;

   for (int k = deg1 + deg2; k >= 0; --k)
   {
      coeff_type c = coeff_type();

      for (int i = std::max (k - deg2, 0); i <= std::min (k, deg1); ++i)
      {
         a->addTo (c, a->mul ((*p1)[i], (*p2)[k - i]));
      }

      if (nonZero || ! a->is0 (c))
      {
         nonZero = true;
         mulAndAdd (c);
      }
   }
}


/**
 *  mul()  --  with coeff_type
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::MulCoeff,A,T> x)
{
   const A* a = x.a;
   const P* p = x.p;
   const T* u = x.u;
   
   reserve (p->numCoeff());

   const CDownI end = p->toA0();
   for (CDownI i = p->fromLc(); i != end; ++i)  mulAndAdd (a->mul (*i, *u));
}


/**
 *  mulBy()  --  with coeff type
 */

template<typename A>
void
HIntLib::Private::PRBA<A>::mulBy (type &p, const coeff_type& u) const
{
   const typename type::DownI end = p.toA0();
   for (typename type::DownI i = p.fromLc(); i != end; ++i)  a.mulBy (*i, u);
}


/**
 *  divByLinearFactor()
 *
 *  Given a polynomial  p  and a ring element  u, set  p := p / (X - u).
 */

template<typename A>
void
HIntLib::Private::PRBA<A>::divByLinearFactor (
      type &p, const coeff_type& u) const
{
   if (p.degree() < 1)  throw DivisionByZero();
   p.divByX();

   const typename type::DownI end = p.toA0();
   for (typename type::DownI i = p.fromLc() + 1; i != end; ++i)
   {
      a.addTo (*i, a.mul (u, *(i-1)));
   }
}


/**
 *  sqr()
 *
 *  For polynomials over domains, the leading coefficients cannot cancel out
 *  each other.
 */

// If there are zero divisors, watch out for cancelation

template<typename T> template<typename A, typename Y>
HIntLib::Polynomial<T>::
Polynomial (Private::PG<Private::Sqr<zerodivisor_tag,Y>,A,T> x)
{
   const A* a = x.a;
   const P* p = x.p;
   const int deg = p->degree();

   if (deg < 0)  return;

   reserve (2 * deg + 1);
   bool nonZero = false;

   for (int k = 2 * deg; k >= 0; --k)
   {
      coeff_type c = coeff_type();

      const int ub = (k + 1) / 2;
      for (int i = std::max (k - deg, 0); i < ub; ++i)
      {
         a->addTo (c, a->mul ((*p)[i], (*p)[k - i]));
      }

      a->times2 (c);

      if ((k & 1) == 0)  a->addTo (c, a->sqr ((*p)[k / 2]));

      if (nonZero || ! a->is0 (c))
      {
         nonZero = true;
         mulAndAdd (c);
      }
   }
}

// general case in an integral domain

template<typename T> template<typename A, typename Y>
HIntLib::Polynomial<T>::
Polynomial (Private::PG<Private::Sqr<nozerodivisor_tag,Y>,A,T> x)
{
   const A* a = x.a;
   const P* p = x.p;
   const int deg = p->degree();

   if (deg < 0)  return;

   reserve (2 * deg + 1);

   for (int k = 2 * deg; k >= 0; --k)
   {
      coeff_type c = coeff_type();

      const int ub = (k + 1) / 2;
      for (int i = std::max (k - deg, 0); i < ub; ++i)
      {
         a->addTo (c, a->mul ((*p)[i], (*p)[k - i]));
      }

      a->times2 (c);

      if ((k & 1) == 0)  a->addTo (c, a->sqr ((*p)[k / 2]));

      mulAndAdd (c);
   }
}

// char = 2  can be done in linear time

template<typename T> template<typename A>
HIntLib::Polynomial<T>::
Polynomial (Private::PG<Private::Sqr<nozerodivisor_tag,char_two>,A,T> x)
{
   const A* a = x.a;
   const P* p = x.p;
   const int deg = p->degree();

   if (deg < 0)  return;

   reserve (2 * deg + 1);

   mulAndAdd (a->sqr ((*p)[deg]));

   for (int k = deg - 1; k >= 0; --k)
   {
      mulByX();
      mulAndAdd (a->sqr ((*p)[k]));
   }
}



/**
 *  power()
 *
 *  This is Algorithm 1.2.3 (Left-Right Binary) from Herni Cohen, CANT.
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::Power,A,T> data)
{
   const A* a = data.a;
   unsigned exponent = data.n;

   if (! exponent)
   {
      mulAndAdd (a->one());
      return;
   }

   unsigned mask = 1u << ms1 (exponent);

   c = data.p->c;
   
   typedef typename A::zerodivisor_category ZC;
   typedef typename A::char_category CC;
   const Private::PG<Private::Mul<ZC>,A,T> mulData (*a, *this, *(data.p));
   const Private::PG<Private::Sqr<ZC,CC>,A,T> sqrData (*a, *this);
   
   while (mask >>= 1)
   {
      {
         P t (sqrData);
         swap (t);
      }

      if (exponent & mask)
      {
         P t (mulData);
         swap (t);
      }
   }
}


/**
 *  derivative()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::Derivative,A,T> x)
{
   const A* a = x.a;
   const P* p = x.p;

   if (p->numCoeff() <= 1)  return;
   
   reserve (p->degree());
   bool nonZero = false;

   for (int i = p->degree(); i > 0; --i)
   {
      coeff_type c = a->times ((*p)[i], i);

      if (nonZero || ! a->is0 (c))
      {
         nonZero = true;
         mulAndAdd (c);
      }
   }
}


/**
 *  evaluate()
 */

template<typename A>
typename A::type
HIntLib::Private::PRBA<A>::evaluate (const type& p, const coeff_type& x) const
{
   if (is0 (p))  return coeff_type();
   if (a.is0 (x))  return p.ct();

   coeff_type res = coeff_type();
   const typename type::CDownI end = p.toA0();

   if (a.is1 (x))
   {
      for (typename type::CDownI i = p.fromLc(); i != end; ++i)  
      {
         a.addTo (res, *i);
      }
   }
   else
   {
      for (typename type::CDownI i = p.fromLc(); i != end; ++i)  
      {
         res = a.add (a.mul (res, x), *i);
      }
   }
   return res;
}


/**
 *  operator<<
 */

template<typename A>
std::ostream&
HIntLib::Private::operator<< (std::ostream &o, const PRBA<A> &a)
{
   Printer ss (o);
   ss << a.getCoeffAlgebra();
   a.printVariableWithBrackets (ss);
   return o;
}

#ifdef HINTLIB_BUILD_WCHAR
template<typename A>
std::wostream&
HIntLib::Private::operator<< (std::wostream &o, const PRBA<A> &a)
{
   WPrinter ss (o);
   ss << a.getCoeffAlgebra();
   a.printVariableWithBrackets (ss);
   return o;
}
#endif


/**
 *  printShort()
 */

template<typename A>
void
HIntLib::Private::PRBA<A>::
printShort (std::ostream& o, const type& p, PrintShortFlag f) const
{
   // The zero-polynomial is a special case

   if (p.is0())
   {
      a.printShort (o, coeff_type(), f);  // print zero-element of coeff_algebra
      return;
   }

   if (a.is0 (p.lc()))  throw InternalError (__FILE__, __LINE__);

   // count non-zero terms

   int nonZeroTerms = 0;

   for (int i = 0; i <= p.degree(); ++i)
   {
      if (! a.is0 (p[i]))  ++nonZeroTerms;
      if (nonZeroTerms >= 2)  break;
   }
   
   Printer ss (o);

#ifdef HINTLIB_ENCODING_LOCALE
   const bool utf8 = ss.utf8();
#endif
   bool needPlus = o.flags() & o.showpos;
   PrintShortFlag flag = PrintShortFlag(f | FIT_FOR_MUL);

   if ((f & FIT_FOR_MUL) && nonZeroTerms >= 2)
   {
      if (needPlus)
      {
         ss << '+';
         needPlus = false;
      }
      ss << '(';
   }

   for (int i = p.degree(); i >= 0; --i)
   {
      // print coefficient

      const coeff_type coef = p[i];
      if (a.is0(coef))  continue;

      bool coefPrinted = false;

      if (a.is1(coef))
      {
         if (needPlus)  ss << '+';
      }
      else if (a.is1(a.neg(coef)))
      {
         ss.minusSign();
      }
      else
      {
         if (i == 0 && ! (f & FIT_FOR_MUL && nonZeroTerms >= 2)) flag = f;
         if (needPlus) ss.setf (ss.showpos);
         a.printShort (ss, coef, flag);
         coefPrinted = true;
      }

      // print x^i

      if (i == 0)
      {
         ss.unsetf (ss.showpos);
         if (! coefPrinted)  a.printShort (ss, a.one());
      }
      else
      {
         printVariable (ss);
         if (i >= 2)  ss.power (i);
      }

      needPlus = true;
   }

   if ((f & FIT_FOR_MUL) && nonZeroTerms >= 2)  ss << ')';
}

#ifdef HINTLIB_BUILD_WCHAR
template<typename A>
void
HIntLib::Private::PRBA<A>::
printShort (std::wostream& o, const type& p, PrintShortFlag f) const
{
   // The zero-polynomial is a special case

   if (p.is0())
   {
      a.printShort (o, coeff_type(), f);  // print zero-element of coeff_algebra
      return;
   }

   if (a.is0 (p.lc()))  throw InternalError (__FILE__, __LINE__);

   // count non-zero terms

   int nonZeroTerms = 0;

   for (int i = 0; i <= p.degree(); ++i)
   {
      if (! a.is0 (p[i]))  ++nonZeroTerms;
      if (nonZeroTerms >= 2)  break;
   }
   
   WPrinter ss (o);

   bool needPlus = o.flags() & o.showpos;
   PrintShortFlag flag = PrintShortFlag(f | FIT_FOR_MUL);

   if ((f & FIT_FOR_MUL) && nonZeroTerms >= 2)
   {
      if (needPlus)
      {
         ss << L'+';
         needPlus = false;
      }
      ss << L'(';
   }

   for (int i = p.degree(); i >= 0; --i)
   {
      // print coefficient

      const coeff_type coef = p[i];
      if (a.is0(coef))  continue;

      bool coefPrinted = false;

      if (a.is1(coef))
      {
         if (needPlus)  ss << L'+';
      }
      else if (a.is1(a.neg(coef)))
      {
         ss.minusSign();
      }
      else
      {
         if (i == 0 && ! (f & FIT_FOR_MUL && nonZeroTerms >= 2)) flag = f;
         if (needPlus) ss.setf (ss.showpos);
         a.printShort (ss, coef, flag);
         coefPrinted = true;
      }

      // print x^i

      if (i == 0)
      {
         ss.unsetf (ss.showpos);
         if (! coefPrinted)  a.printShort (ss, a.one());
      }
      else
      {
         printVariable (ss);
         if (i >= 2)  ss.power (i);
      }

      needPlus = true;
   }

   if ((f & FIT_FOR_MUL) && nonZeroTerms >= 2)  ss << L')';
}
#endif


/**
 *  print()
 */

template<typename A>
void
HIntLib::Private::PRBA<A>::print (std::ostream& o, const type& p) const
{
   Printer ss (o);

   printShort (ss, p);
   ss << ' ';
   printSuffix (ss);
}

#ifdef HINTLIB_BUILD_WCHAR
template<typename A>
void
HIntLib::Private::PRBA<A>::print (std::wostream& o, const type& p) const
{
   WPrinter ss (o);

   printShort (ss, p);
   ss << L' ';
   printSuffix (ss);
}
#endif


/*********************  Instantiations  **************************************/

/**
 *  Instantiate Polynomial<>
 */

#define HINTLIB_INSTANTIATE_POLYNOMIAL_NO_EQUAL(X) \
   template Polynomial<X >::Polynomial(const Polynomial<X >&); \
   template Polynomial<X >::~Polynomial(); \
   template Polynomial<X >& Polynomial<X >::divByX (unsigned); \
   template Polynomial<X >& Polynomial<X >::mulByX (unsigned);

#ifdef HINTLIB_EQUAL_BUG
#define HINTLIB_INSTANTIATE_POLYNOMIAL(X) \
   HINTLIB_INSTANTIATE_POLYNOMIAL_NO_EQUAL(X) \
   template bool operator== (const Polynomial<X >&, const Polynomial<X >&);
#else
#define HINTLIB_INSTANTIATE_POLYNOMIAL(X) \
   HINTLIB_INSTANTIATE_POLYNOMIAL_NO_EQUAL(X)
#endif


/**
 *  Instantiate Polynomial Ring Base
 */

#ifdef HINTLIB_BUILD_WCHAR
#define HINTLIB_INSTANTIATE_POLYNOMIALRING_BB_W(Y) \
   template std::wostream& operator<< (std::wostream &, const PRBA<Y > &); \
   template void PRBA<Y >::print (std::wostream &, const type&) const; \
   template void PRBA<Y >::printShort \
               (std::wostream &, const type&, PrintShortFlag) const;
#else
#define HINTLIB_INSTANTIATE_POLYNOMIALRING_BB_W(Y)
#endif

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_BB(Y) \
   template Polynomial<Y::type>::Polynomial(Private::PG<Private::X,Y >); \
   template Polynomial<Y::type>::Polynomial(Private::PG<Private::Xn,Y >); \
   template Polynomial<Y::type>::Polynomial \
      (Private::PG<Private::Element<Y::size_category>,Y >); \
   template Polynomial<Y::type>::Polynomial \
      (Private::PG<Private::ElementMonic<Y::size_category>,Y >); \
   template Polynomial<Y::type>::Polynomial \
      (Private::PG<Private::Neg<Y::char_category>,Y >); \
   template Polynomial<Y::type>::Polynomial \
      (Private::PG<Private::Dbl<Y::char_category>,Y >); \
   template Polynomial<Y::type>::Polynomial \
      (Private::PG<Private::Times<Y::char_category>,Y >); \
   template Polynomial<Y::type>::Polynomial(Private::PG<Private::Add,Y >); \
   template Polynomial<Y::type>::Polynomial(Private::PG<Private::Sub,Y >); \
   template Polynomial<Y::type>::Polynomial \
      (Private::PG<Private::Sqr<Y::zerodivisor_category, \
                                Y::char_category>,Y >); \
   template Polynomial<Y::type>::Polynomial \
      (Private::PG<Private::Mul<Y::zerodivisor_category>,Y >); \
   template Polynomial<Y::type>::Polynomial(Private::PG<Private::MulCoeff,Y >);\
   template Polynomial<Y::type>::Polynomial(Private::PG<Private::Power,Y >); \
   template Polynomial<Y::type>::Polynomial \
      (Private::PG<Private::Derivative,Y >); \
   namespace Private { \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_BB_W(Y) \
   template unsigned PRBA<Y >::indexImp (const type&, Y::size_category) const; \
   template bool PRBA<Y >::is1 (const type&) const; \
   template bool PRBA<Y >::isMonic (const type&) const; \
   template void PRBA<Y >::negateImp (type&, char_any) const; \
   template void PRBA<Y >::times2Imp (type&, Y::char_category) const; \
   template void PRBA<Y >::addTo (type&,const type&) const; \
   template void PRBA<Y >::subFrom (type&,const type&) const; \
   template void PRBA<Y >::mulBy (type&, const coeff_type&) const; \
   template void PRBA<Y >::divByLinearFactor (type &, const coeff_type&) const;\
   template Y::type PRBA<Y >::evaluate (const type&, const coeff_type&) const; \
   template std::ostream& operator<< (std::ostream &, const PRBA<Y > &); \
   template void PRBA<Y >::print (std::ostream &, const type&) const; \
   template void PRBA<Y >::printShort \
               (std::ostream &, const type&, PrintShortFlag) const; \
   }

#endif

