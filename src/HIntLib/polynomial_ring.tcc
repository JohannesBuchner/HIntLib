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

#ifndef HINTLIB_POLYNOMIAL_RING_TCC
#define HINTLIB_POLYNOMIAL_RING_TCC 1

#ifndef HINTLIB_POLYNOMIAL_TCC
#include <HIntLib/polynomial.tcc>
#endif


/********************  Polynomial Ring <ring_tag>  ***************************/


/**
 *  additive Order()
 */

template<class A>
unsigned
HIntLib::Private::PRBA_Ring<A>::additiveOrder (const type& u) const
{
   unsigned n = 1;

   const typename type::CDownI end = u.toA0();
         typename type::CDownI i   = u.fromLc();

   while (i != end)  n = lcm (n, a.additiveOrder (*i++));
   return n;
}


/**
 *  isUnit()
 */

template<class A>
bool
HIntLib::Private::PRBA_Ring<A>::isUnit (const type& p) const
{
   if (p.numCoeff() == 0 || ! a.isUnit (p.ct()))  return false;

   const typename type::CDownI end = p.toA0() - 1;
         typename type::CDownI i   = p.fromLc();

   while (i != end)  if (! a.isNilpotent (*i++))  return false;
   return true;
}


/**
 *  isNilpotent()
 *
 *  A polynomial is nilpotent if all its coefficients are nilpotent.
 */

template<class A>
bool
HIntLib::Private::PRBA_Ring<A>::isNilpotent (const type& p) const
{
   const typename type::CDownI end = p.toA0();
         typename type::CDownI i   = p.fromLc();

   while (i != end)  if (! a.isNilpotent (*i++))  return false;
   return true;
}


/**
 *  unitRecip()
 *
 *  A unit has the form  u = 1 - b, where b is nilpotent.
 *
 *  The inverse of u is  u^-1 = (1 + b + b^2 + ...)
 *
 *  See, e.g.: Wikipedia, keyword: Nilpotent
 */

template<typename T> template<typename CA>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::UnitRecipR,CA,T> x)
{
   typedef Private::PRBA_Ring<CA> A;
   const A* a = x.a;
   const CA& ca = a->getCoeffAlgebra();
   const P* u = x.p;

   P b (*u);

   typename CA::unit_type runit = ca.unitRecip (ca.toUnit (b.ct()));

   // Construct b

   a->mulBy (b, runit);
   a->negate (b);
   b.ct() = coeff_type();

   // Initialize result to  1 + b

   c = b.c;
   ca.addTo (ct(), ca.one());

   // add powers of  b  until these powers become 0

   P power (b);

   for (;;)
   {
      a->mulBy (power, b);
      if (a->is0 (power))  break;
      a->addTo (*this, power);
   }

   // multiply with unit

   a->mulBy (*this, runit);
}


/********************  Polynomial Ring <domain_tag> and subclasses  **********/


/**
 *  mulUnit ()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::MulUnitD,A,T> x)
{
   const A* a = x.a;
   const P* p = x.p;
   const typename A::unit_type* u = x.u;
   
   reserve (p->numCoeff());

   const P::CDownI end = p->toA0();
   for (P::CDownI i = p->fromLc(); i != end; ++i)
   {
      mulAndAdd (a->mulUnit (*i, *u));
   }
}


/**
 *  mulByUnit()
 */

template<typename A>
void
HIntLib::Private::PRBA_Domain<A>:: mulByUnit (type &p, const unit_type& u) const
{
   const typename type::DownI end = p.toA0();
   for (typename type::DownI i = p.fromLc(); i != end; ++i)
   {
      a.mulByUnit (*i, u);
   }
}


/**
 *  isDivisor()
 */

template<class A>
bool
HIntLib::Private::PRBA_Domain<A>::isDivisor (const type& u, const type& v) const
{
   if (v.is0())  throwDivisionByZero();
   if (&u == &v || u.is0())  return true;
   if (u.numCoeff() < v.numCoeff())  return false;

   // we could get rid of one iteration here, but its probably not worth it

   type uu (u);

   typedef typename A::type coeff_type;
   typedef typename type:: DownI  I;
   typedef typename type::CDownI CI;

          I ubegin = uu.fromLc();
   const  I uend   = uu.toA0();
   const CI vend   = v.toA0();

   for (;;)
   {
      coeff_type factor;
      if (! a.isDivisor (*ubegin, v.lc(), factor))  return false;

      I ui = ++ubegin;

      for (CI vi = v.fromLc() + 1; vi != vend; )
      {
         a.subFrom (*ui++, a.mul (factor, *vi++));
      }

      if (ui == uend)  break;
   }

   typename type::DownI rend = uu.toA0();
      
   for ( ; ubegin != rend; ++ubegin)
   {
      if (! a.is0(*ubegin))  return false;
   }

   return true;
}

template<class A>
bool
HIntLib::Private::PRBA_Domain<A>::
isDivisor (const type& u, const type& v, type& result) const
{
   if (v.is0())  throwDivisionByZero();
   if (&u == &v || u.is0())
   {
      result.makeZero();
      return true;
   }
   if (u.numCoeff() < v.numCoeff())  return false;

#if 0
   if (&result == &u)
   {
      typename type::DownI qend = polyReduceS (a, result, v);
      typename type::DownI rend = result.toA0();
      
      for (typename type::DownI i = qend; i != rend; ++i)
      {
         if (! a.is0(*i))  return false;
      }

      result.getC().erase (qend, rend);
      return true;
   } 
#endif
   else
   {
      type uu (u);

#if 0
      if (&result == &v)  // watch out for aliasing
      {
         type temp;
         temp.reserve (u.numCoeff() - v.numCoeff() + 1);
         bool ok = findLc (a, uu, polyReduce (a, uu, v, temp)) == uu.toA0();
         if (ok)  destructiveAssign (result, temp);
         return ok;
      }
      else
#endif
      {
         result.makeZero();
         result.reserve (u.numCoeff() - v.numCoeff() + 1);

         // return findLc (a, uu, polyReduce (a, uu, v, result)) == uu.toA0();
   typedef typename A::type coeff_type;
   typedef typename type:: DownI  I;
   typedef typename type::CDownI CI;

          I ubegin = uu.fromLc();
   const  I uend   = uu.toA0();
   const CI vend   = v.toA0();

   for (;;)
   {
      coeff_type factor;
      if (! a.isDivisor (*ubegin, v.lc(), factor))  return false;
      result.mulAndAdd (factor);

      I ui = ++ubegin;

      for (CI vi = v.fromLc() + 1; vi != vend; )
      {
         a.subFrom (*ui++, a.mul (factor, *vi++));
      }

      if (ui == uend)  break;
   }

   typename type::DownI rend = uu.toA0();
      
   for ( ; ubegin != rend; ++ubegin)
   {
      if (! a.is0(*ubegin))  return false;
   }

   return true;
      }
   }
}


/**
 *  div()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::DivDomain,A,T> x)
{
   const P* u = x.p1;
   const P* v = x.p2;
   
   if (v->is0())  throwDivisionByZero();
   if (u->is0())  return;

   int num = int (u->numCoeff()) - int (v->numCoeff()) + 1;
   if (num <= 0)  throwDivisionByZero();
   reserve (num);

   P uu (*u);

          DownI ubegin = uu.fromLc();
   const  DownI uend   = uu.toA0();
   const CDownI vend   = v->toA0();

   const A* a = x.a;

   for (;;)
   {
      coeff_type factor = a->div (*ubegin, v->lc());
      mulAndAdd (factor);

      DownI ui = ++ubegin;

      for (CDownI vi = v->fromLc() + 1; vi != vend; )
      {
         a->subFrom (*ui++, a->mul (factor, *vi++));
      }

      if (ui == uend)  break;
   }
}


/**
 *  divBy()
 */

template<typename A>
void
HIntLib::Private::PRBA_Domain<A>::divBy (type& u, const type& v) const
{
   if (v.is0())  throwDivisionByZero();
   if (u.is0())  return;
   if (u.numCoeff() < v.numCoeff())  throw DivisionByZero();

   typedef typename A::type coeff_type;
   typedef typename type:: DownI  I;
   typedef typename type::CDownI CI;

          I ubegin = u.fromLc();
   const  I uend   = u.toA0();
   const CI vend   = v.toA0();

   for (;;)
   {
      a.divBy (*ubegin, v.lc());

      I ui = ubegin + 1;

      for (CI vi = v.fromLc() + 1; vi != vend; )
      {
         a.subFrom (*ui++, a.mul (*ubegin, *vi++));
      }

      ++ubegin;

      if (ui == uend)  break;
   }

   u.getC().erase (ubegin, u.toA0());
}


/**
 *  order()
 *
 *  This is a copy of the identical routine for fields.
 *  Copying it literaly avoids multiple inhertance.
 */

template<typename A>
unsigned
HIntLib::Private::PRBA_Domain<A>::order (const type& p) const
{
   const unsigned num = p.numCoeff();

   if (! num)  throwDivisionByZero();  // p == 0

   // if there is an x, it will never go away because A is a domain
   if (num > 1)  return 0;

   return a.order (p.lc());   // no x, no polynomial
}


/********************  Polynomial Ring <ufd_tag>  ****************************/


/**
 *  isCanonical ()
 */

template<class A>
bool
HIntLib::Private::PRBA_UFD<A>::isCanonical (const type &p) const
{
   return p.degree() < 0 || a.isCanonical(p.lc());
}


/**
 *  make Canonical ()
 */

template<typename A>
typename HIntLib::Private::PRBA_UFD<A>::unit_type
HIntLib::Private::PRBA_UFD<A>::makeCanonical (type& p) const
{
   if (is0 (p)) return a.toUnit (a.one());

   unit_type l = a.makeCanonical (p.lc());
   unit_type il = a.unitRecip (l);

   const typename type::DownI end = p.toA0();
   for (typename type::DownI i = p.fromLc() + 1; i != end; ++i)
   {
      a.mulByUnit (*i, il);
   }

   return l;
}


/*********************  Instantiations  **************************************/


// Instantiate PolynomialRingBase<ring_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_RING(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_BB(X) \
   template Polynomial<X::type>::Polynomial \
      (Private::PG<Private::UnitRecipR,X >);\
   namespace Private { \
   template bool PRBA_Ring<X >::isUnit (const type& u) const; \
   template bool PRBA_Ring<X >::isNilpotent (const type& u) const; \
   template unsigned PRBA_Ring<X >::additiveOrder (const type& u) const; \
   }

// Instantiate PolynomialRingBase<domain_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_DOMAIN(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_BB(X) \
   template Polynomial<X::type>::Polynomial(Private::PG<Private::MulUnitD,X >);\
   template Polynomial<X::type>::Polynomial \
      (Private::PG<Private::DivDomain,X >);\
   namespace Private { \
   template void PRBA_Domain<X >::mulByUnit (type&, const unit_type&) const; \
   template unsigned PRBA_Domain<X >::order (const type&) const; \
   template void PRBA_Domain<X >::divBy (type&, const type&) const; \
   template bool PRBA_Domain<X >::isDivisor \
      (const type&, const type&) const; \
   template bool PRBA_Domain<X >::isDivisor \
      (const type&, const type&, type&) const; \
   }

// Instantiate PolynomialRingBase<ufd_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_UFD(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_DOMAIN(X) \
   namespace Private { \
   template bool PRBA_UFD<X >::isCanonical (const type&) const; \
   template PRBA_UFD<X >::unit_type \
            PRBA_UFD<X >::makeCanonical (type&) const; \
   }

#endif

