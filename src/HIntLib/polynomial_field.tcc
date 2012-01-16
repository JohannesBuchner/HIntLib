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

#ifndef HINTLIB_POLYNOMIAL_FIELD_TCC
#define HINTLIB_POLYNOMIAL_FIELD_TCC 1

#ifndef HINTLIB_POLYNOMIAL_TCC
#include <HIntLib/polynomial.tcc>
#endif

#include <HIntLib/gcd.tcc>


/********************  Polynomial Ring <field_tag> and subclasses  ***********/


/**
 *  order()
 *
 *  This is a copy of the identical routine for domains.
 *  Copying it literaly avoids multiple inhertance.
 */

template<typename A>
unsigned
HIntLib::Private::PRBA_Field<A>::order (const type& p) const
{
   const unsigned num = p.numCoeff();

   if (! num)  throwDivisionByZero();  // p == 0

   // if there is an x, it will never go away because A is a domain
   if (num > 1)  return 0;

   return this->a.order (p.lc());   // no x, no polynomial
}


/**
 *  polyReduce()
 *
 *  Reduces a polynomial  u  modulo  v.
 *
 *  The result is stored in the  deg v - 1  least significant coefficients of u.
 *
 *  Preconditions:  deg v != 0  and  deg u >= deg v
 *
 *  See Knuth, TACP, vol 2, 4.6.1, Algo D
 */

namespace HIntLib
{
namespace Private
{

template<typename A>
typename Polynomial<typename A::type>::DownI
polyReduce (const A& a, Polynomial<typename A::type>& u,
                  const Polynomial<typename A::type>& v)
{
   typedef typename A::type coeff_type;
   typedef Polynomial<coeff_type> P;
   typedef typename P:: DownI  I;
   typedef typename P::CDownI CI;

   const coeff_type q = a.recip (v.lc());

          I ubegin = u.fromLc();
   const  I uend   = u.toA0();
   const CI vend   = v.toA0();

   for (;;)
   {
      const coeff_type factor = a.mul (*ubegin, q);
      I ui = ++ubegin;

      for (CI vi = v.fromLc() + 1; vi != vend; )
      {
         a.subFrom (*ui++, a.mul (factor, *vi++));
      }

      if (ui == uend)  break;
   }

   return ubegin;
}

template<typename A>
typename Polynomial<typename A::type>::DownI
polyReduce (const A& a, Polynomial<typename A::type>& u,
                  const Polynomial<typename A::type>& v,
                        Polynomial<typename A::type>& result)

{
   typedef typename A::type coeff_type;
   typedef Polynomial<coeff_type> P;
   typedef typename P:: DownI  I;
   typedef typename P::CDownI CI;

   const coeff_type q = a.recip (v.lc());

          I ubegin = u.fromLc();
   const  I uend   = u.toA0();
   const CI vend   = v.toA0();

   for (;;)
   {
      const coeff_type factor = a.mul (*ubegin, q);
      result.mulAndAdd (factor);
      I ui = ++ubegin;

      for (CI vi = v.fromLc() + 1; vi != vend; )
      {
         a.subFrom (*ui++, a.mul (factor, *vi++));
      }

      if (ui == uend)  break;
   }

   return ubegin;
}

template<typename A>
typename Polynomial<typename A::type>::DownI
polyReduceS (const A& a, Polynomial<typename A::type>& u,
                   const Polynomial<typename A::type>& v)
{
   typedef typename A::type coeff_type;
   typedef Polynomial<coeff_type> P;
   typedef typename P:: DownI  I;
   typedef typename P::CDownI CI;

   const coeff_type q = a.recip (v.lc());

          I ubegin = u.fromLc();
   const  I uend   = u.toA0();
   const CI vend   = v.toA0();

   for (;;)
   {
      a.mulBy (*ubegin, q);
      I ui = ubegin + 1;

      for (CI vi = v.fromLc() + 1; vi != vend; )
      {
         a.subFrom (*ui++, a.mul (*ubegin, *vi++));
      }

      ++ubegin;

      if (ui == uend)  break;
   }

   return ubegin;
}

template<typename A>
typename Polynomial<typename A::type>::DownI
findLc (const A& a, Polynomial<typename A::type>& p,
        typename Polynomial<typename A::type>::DownI i)
{
   typename Polynomial<typename A::type>::DownI end = p.toA0();
   while (i != end && a.is0 (*i))  ++i;
   return i;
}

} // namespace Private
} // namespace HIntLib


/**
 *  div ()
 *
 *  Division of polynomials
 */

template<class A>
void
HIntLib::Private::PRBA_Field<A>::
div (const type& u, const type& v, type& q, type& r) const
{
   if (v.is0())  throwDivisionByZero();

   int num = int (u.numCoeff()) - int (v.numCoeff()) + 1;

   if (&u == &r)
   {
      q.makeZero();

      if (num <= 0)  return;

      q.reserve (num);
      r.getC().erase (r.fromLc(),
                      findLc (this->a, r, polyReduce (this->a, r, v, q)));
   }
   else if (&u == &q)
   {
      if (num <= 0)
      {
          destructiveAssign (r, q);
          q.makeZero();
      }

      const typename type::DownI qend   = polyReduceS (this->a, q, v);
      const typename type::DownI rbegin = findLc(this->a, q, qend);
      const typename type::DownI rend   = q.toA0();

      r.makeZero();
      r.reserve (rend - rbegin);
      std::copy (rbegin, rend, back_inserter(r.getC()));

      q.getC().erase (q.fromLc(), qend);
   }
   else
   {
      q.makeZero();

      if (num <= 0)
      {
         r = u;
         return;
      }

      type uu (u);

      q.reserve (num);
      const typename type::CDownI start
            = findLc (this->a, uu, polyReduce (this->a, uu, v, q));
      const typename type::CDownI end   = uu.toA0();

      r.makeZero();
      r.reserve (end - start);
      std::copy (start, end, back_inserter(r.getC()));
   }
}


/**
 *  quot()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::Quot,A,T> x)
{
   const P* u = x.p1;
   const P* v = x.p2;
   
   if (v->is0())  throwDivisionByZero();

   int num = int (u->numCoeff()) - int (v->numCoeff()) + 1;
   if (num <= 0)  return;
   reserve (num);

   P uu (*u);
   Private::polyReduce (*x.a, uu, *v, *this);
}


#if 0
/**
 *  div()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::Div,A,T> x)
{
   const P* u = x.p1;
   const P* v = x.p2;
   
   if (v->is0())  throwDivisionByZero();
   if (u->is0())  return;

   int num = int (u->numCoeff()) - int (v->numCoeff()) + 1;
   if (num <= 0)  throwDivisionByZero();
   reserve (num);

   P uu (*u);
   Private::polyReduce (*x.a, uu, *v, *this);
}
#endif


/**
 *  quotient()
 */

template<typename A>
void
HIntLib::Private::PRBA_Field<A>::quotient (type& u, const type& v) const
{
   if (v.is0())  throwDivisionByZero();
   if (u.numCoeff() < v.numCoeff())
   {
      u.makeZero();
      return;
   }

   u.getC().erase (polyReduceS (this->a, u, v), u.toA0());
}


#if 0
/**
 *  divBy()
 */

template<typename A>
void
HIntLib::Private::PRBA_Field<A>::divBy (type& u, const type& v) const
{
   if (v.is0())  throwDivisionByZero();
   if (u.is0())  return;

   if (u.numCoeff() < v.numCoeff())  throwDivisionByZero();

   u.getC().erase (polyReduceS (this->a, u, v), u.toA0());
}
#endif


/**
 *  rem()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::Rem,A,T> x)
{
   const P* v = x.p2;
   if (v->is0())  throwDivisionByZero();

   const P* u = x.p1;
   if (u->numCoeff() < v->numCoeff())
   {
      c = u->c;
      return;
   }

   P uu (*u);

   const DownI begin
      = Private::findLc (*x.a, uu, Private::polyReduce (*x.a, uu, *v));
   const DownI end   = uu.toA0();
   
   reserve (end - begin);
   std::copy (begin, end, back_inserter(c));
}


/**
 *  reduce()
 */

template<typename A>
void
HIntLib::Private::PRBA_Field<A>::reduce (type& u, const type& v) const
{
   if (v.is0())  throwDivisionByZero();
   if (u.numCoeff() < v.numCoeff())  return;

   u.getC().erase (u.fromLc(), findLc (this->a, u, polyReduce (this->a, u, v)));
}


/**
 *  isDivisor()
 */

template<class A>
bool
HIntLib::Private::PRBA_Field<A>::isDivisor (const type& u, const type& v) const
{
   if (v.is0())  throwDivisionByZero();
   if (&u == &v || u.is0())  return true;
   if (u.numCoeff() < v.numCoeff())  return false;

   // we could get rid of one iteration here, but its probably not worth it

   type uu (u);

   typename type::DownI rend = uu.toA0();
      
   for (typename type::DownI i = polyReduce (this->a, uu, v); i != rend; ++i)
   {
      if (! this->a.is0(*i))  return false;
   }

   return true;
}

template<class A>
bool
HIntLib::Private::PRBA_Field<A>::
isDivisor (const type& u, const type& v, type& result) const
{
   if (v.is0())  throwDivisionByZero();
   if (&u == &v || u.is0())
   {
      result.makeZero();
      return true;
   }
   if (u.numCoeff() < v.numCoeff())  return false;

   if (&result == &u)
   {
      typename type::DownI qend = polyReduceS (this->a, result, v);
      typename type::DownI rend = result.toA0();
      
      for (typename type::DownI i = qend; i != rend; ++i)
      {
         if (! this->a.is0(*i))  return false;
      }

      result.getC().erase (qend, rend);
      return true;
   } 
   else
   {
      type uu (u);

      if (&result == &v)  // watch out for aliasing
      {
         type temp;
         temp.reserve (u.numCoeff() - v.numCoeff() + 1);
         bool ok =   findLc (this->a, uu, polyReduce (this->a, uu, v, temp))
                  == uu.toA0();
         if (ok)  destructiveAssign (result, temp);
         return ok;
      }
      else
      {
         result.makeZero();
         result.reserve (u.numCoeff() - v.numCoeff() + 1);

         return findLc (this->a, uu, polyReduce (this->a, uu, v, result))
             == uu.toA0();
      }
   }
}


/**
 *  make Canonical ()
 */

template<typename A>
typename HIntLib::Private::PRBA_Field<A>::unit_type
HIntLib::Private::PRBA_Field<A>::makeCanonical (type &p) const
{
   if (is0(p) || isCanonical(p))  return this->a.one();

   unit_type l = p.lc();
   unit_type il = this->a.recip (l);

   p.lc() = this->a.one();

   const typename type::DownI end = p.toA0();
   for (typename type::DownI i = p.fromLc() + 1; i != end; ++i)
   {
      this->a.mulBy (*i, il);
   }

   return l;
}


/**
 *  is Associate ()
 *
 *  Determines if  p1  is associate to  p2, i.e. if there exists a unit  m
 *  such that p2 * m = p1
 */

template<typename A>
bool
HIntLib::Private::PRBA_Field<A>::
isAssociate (const type &p1, const type &p2, unit_type &u) const
{
   if (p1.degree() != p2.degree())  return false;
   if (is0(p1))
   {
      u = this->a.one();
      return true;
   }

   u = this->a.div (p1.lc(), p2.lc());

   const typename type::CDownI end1 = p1.toA0();
   for (typename type::CDownI i1 = p1.fromLc() + 1,
                              i2 = p2.fromLc() + 1; i1 != end1; ++i1, ++i2)
   {
      if (! (this->a.mul (*i2, u) == *i1))  return false;
   }

   return true;
}


/**
 *  charRoot()
 *
 *  Determins the char-th root from an element in the coefficient field.
 */

namespace HIntLib { namespace Private
{

template<typename CA>
inline
bool
charRoot (const CA& a, typename CA::type& c, gf_tag)
{
   c = a.invFrobenius (c);
   return true;
}

template<typename CA>
inline
bool
charRoot (const CA& a, typename CA::type& c, ratfunfield_tag)
{
   if (a.is0(c) || a.is1(c))  return true;

   typedef typename CA::base_algebra PA;
   const PA& pa = a.getBaseAlgebra();
   typedef typename PA::Factorization F;

   const unsigned ch = a.characteristic();

   // factor numerator into squarefree polynomials and check result

   F facNum;
   typename PA::unit_type uNum =
      pa.squarefreeFactor (facNum, c.numerator()); 

   for (typename F::const_iterator i = facNum.begin(); i != facNum.end(); ++i)
   {
      if (i->second % ch != 0)  return false;
   }

   // factor denominator into squarefree polynomials and check result

   F facDen;
   pa.squarefreeFactor (facDen, c.denominator()); 

   for (typename F::const_iterator i = facDen.begin(); i != facDen.end(); ++i)
   {
      if (i->second % ch != 0)  return false;
   }

   // build numerator

   typename PA::type num = pa.mul (pa.one(), uNum);

   for (typename F::const_iterator i = facNum.begin(); i != facNum.end(); ++i)
   {
      pa.mulBy (num, pa.power (i->first, i->second / ch));
   }

   // build denominator

   typename PA::type den = pa.one();

   for (typename F::const_iterator i = facDen.begin(); i != facDen.end(); ++i)
   {
      pa.mulBy (den, pa.power (i->first, i->second / ch));
   }

   // Replace  c  by  c ^ (1/char)

   c = a.makeElement (num, den);

   return true;
}

} }  // namespace Private & HIntLib


/**
 *  isSquarefree ()
 *
 *  Determines if  p  is squarefree, i.e. if it does not contain two identical
 *  prime factors.
 */

namespace HIntLib { namespace Private {

// char = prime, finite field

template<typename A>
inline
bool
isSquarefreeImp (const A& a, const typename A::type& p, char_prime, gf_tag)
{
   if (p.numCoeff() <= 2)  return true;

   const typename A::type der = a.derivative(p);

   return ! der.is0() && genIsCoprime (a, der, p);
}

// char = prime, rational function field

template<typename A>
inline
bool
isSquarefreeImp (const A& a, const typename A::type& p,
                 char_prime, ratfunfield_tag)
{
   if (p.numCoeff() <= 2)  return true;

   const typename A::type der = a.derivative(p);

   if (der.is0())
   {
      const unsigned c = a.characteristic();

      for (typename A::type::CDownI i = p.fromLc(); i < p.toA0(); i += c)
      {
         typename A::coeff_type coeff = *i;

         if (! charRoot (a.getCoeffAlgebra(), coeff, ratfunfield_tag()))
         {
            return true;
         }
      }

      return false;
   }
   else
   {
      return genIsCoprime (a, der, p);
   }
}

// char = 0

template<typename A>
inline
bool
isSquarefreeImp (const A& a, const typename A::type& p, char_zero, field_tag)
{
   return p.numCoeff() <= 2 || genIsCoprime (a, a.derivative(p), p);
}

} }  // namespace Private & HIntLib

// dispatcher

template<typename A>
bool
HIntLib::Private::PRBA_Field<A>::isSquarefree (const type& p) const
{
   return isSquarefreeImp (
      *static_cast<const PolynomialRing<A>*> (this),
      p, typename A::char_category(), typename A::algebra_category());
}


/**
 *  squarefreeFactor()
 *
 *  Factors a polynomial into squarefree and coprime factors.
 */

namespace HIntLib { namespace Private {

/*
 *  For char = prime
 *
 *  This is essentially the algorithm from Cohen, CANT, Algorithm 3.4.2
 *
 *  The main difference is that we have to be much more carefull when reducing
 *  p(x^c) to q(x), since we have to draw c-th roots of the coefficients.
 */

template<typename A>
inline
typename A::unit_type
squarefreeFactorImp (const A& a, typename A::Factorization& f,
                     const typename A::type& p, char_prime)
{
   typedef typename A::type type;

   // Step 1  --  Initialize

   typedef typename A::coeff_algebra CA;

   unsigned e = 1;
   unsigned k = 0;
   const unsigned c = a.characteristic();

   // v ... each factor appears here once (module char problem)
   // t ... the remaining factors
   // v * t ... the polynomial to be factored

   type t0 (p);
   const typename A::unit_type unit = a.makeCanonical (t0);
   type t = genGcd (a, a.derivative (t0), t0); a.makeCanonical (t);
   type v = a.div (t0, t);

   for(;;)
   {
      if (v.isConstant())
      {
         if (t.isConstant())  break;

         t0.makeZero();
         for (typename A::type::CDownI i = t.fromLc(); i < t.toA0(); i += c)
         {
            typename type::coeff_type coeff = *i;

            if (! charRoot (a.getCoeffAlgebra(), coeff,
                            typename CA::algebra_category()))
            {
               f.push_back (std::make_pair (t, e));
               return unit;
            }
            t0.mulAndAdd (coeff);
         }

         t = genGcd (a, a.derivative (t0), t0); a.makeCanonical (t);
         v = a.div (t0, t);
         e *= c;
         k = 0;
      }
      else
      {
         // Step 4  -- Special case

         ++k;
         if (k % c == 0)  // Special treatment for  k  a multiple of char
         {
            a.divBy (t, v);
            ++k;
         }

         // Step 5  --  Compute factor

         // w = product of factors apearing more than once

         type w = genGcd (a, t, v); a.makeCanonical (w);

         if (v.numCoeff() > w.numCoeff())
         {
            f.push_back (std::make_pair (a.div (v,w), e * k));
         }
         destructiveAssign (v, w);
         a.divBy (t, v);
      }
   }

   return unit;
}

/*
 *  char = 0
 *
 *  Same algorithm as above, with all code for treating c|k and p'=0 removed.
 */

template<typename A>
inline
typename A::unit_type
squarefreeFactorImp (const A& a, typename A::Factorization& f,
                     const typename A::type& p, char_zero)
{
   typedef typename A::type type;

   // Step 1  --  Initialize

   unsigned k = 0;

   // v ... each factor appears here once
   // t ... the remaining factors
   // v * t ... the polynomial to be factored

   type t0 (p);
   const typename A::unit_type unit = a.makeCanonical (t0);
   type t = genGcd (a, a.derivative (t0), t0);
   a.makeCanonical (t);
   type v = a.div (t0, t);

   while (! v.isConstant())
   {
      ++k;

      // w = product of factors apearing more than once

      type w = genGcd (a, t, v);
      a.makeCanonical (w);

      if (v.numCoeff() > w.numCoeff())
      {
         f.push_back (std::make_pair (a.div(v,w), k));
      }

      destructiveAssign (v, w);
      a.divBy (t, v);
   }

   return unit;
}

} }  // namespace Private & HIntLib

// dispatcher

template<typename A>
typename HIntLib::Private::PRBA_Field<A>::unit_type
HIntLib::Private::PRBA_Field<A>::
squarefreeFactor (Factorization& f, const type& p) const
{
   if (p.is0())  throwDivisionByZero();
   if (p.isConstant())  return p.lc();
   if (p.degree() == 1)
   {
      f.push_back (std::make_pair (p, 1));
      return makeCanonical (f.back().first);
   }

   return squarefreeFactorImp (
      *static_cast<const PolynomialRing<A>*> (this),
      f, p, typename A::char_category());
}


/*********************  Instantiations  **************************************/


// Instantiate PolynomialRingBase<field_tag> and subclasses

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_FIELD(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_BB(X) \
   template Polynomial<X::type>::Polynomial(Private::PG<Private::Rem,X >);\
   template Polynomial<X::type>::Polynomial(Private::PG<Private::Quot,X >);\
   namespace Private { \
   template unsigned PRBA_Field<X >::order (const type&) const; \
   template void PRBA_Field<X >::div \
      (const type&, const type&, type&, type&) const; \
   template void PRBA_Field<X >::reduce   (type&, const type&) const; \
   template void PRBA_Field<X >::quotient (type&, const type&) const; \
   template Polynomial<X::type>::DownI polyReduce \
      (const X&, Polynomial<X::type>&, const Polynomial<X::type>&); \
   template Polynomial<X::type>::DownI polyReduce \
      (const X&, Polynomial<X::type>&, const Polynomial<X::type>&, \
                 Polynomial<X::type>&); \
   template bool PRBA_Field<X >::isAssociate \
      (const type&, const type&, unit_type&) const; \
   template bool PRBA_Field<X >::isDivisor (const type&, const type&) const; \
   template bool PRBA_Field<X >::isDivisor \
      (const type&, const type&, type&) const; \
   template PRBA_Field<X >::unit_type \
            PRBA_Field<X >::makeCanonical (type&) const; \
   template bool PRBA_Field<X >::isSquarefree (const type&) const; \
   template PRBA_Field<X >::unit_type \
            PRBA_Field<X >::squarefreeFactor \
      (Factorization&, const type&) const; \
   } \
   HINTLIB_INSTANTIATE_GENGCD(PolynomialRing<X >)

#endif

