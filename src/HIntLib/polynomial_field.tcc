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

#ifndef HINTLIB_POLYNOMIAL_FIELD_TCC
#define HINTLIB_POLYNOMIAL_FIELD_TCC 1

#ifndef HINTLIB_POLYNOMIAL_TCC
#include <HIntLib/polynomial.tcc>
#endif

#include <HIntLib/gcd.tcc>
#include <HIntLib/counter.h>
#include <HIntLib/linearalgebra.h>
#include <HIntLib/prime.h>


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

   return a.order (p.lc());   // no x, no polynomial
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
      r.getC().erase (r.fromLc(), findLc (a, r, polyReduce (a, r, v, q)));
   }
   else if (&u == &q)
   {
      if (num <= 0)
      {
          destructiveAssign (r, q);
          q.makeZero();
      }

      const typename type::DownI qend   = polyReduceS (a, q, v);
      const typename type::DownI rbegin = findLc(a, q, qend);
      const typename type::DownI rend   = q.toA0();

      r.makeZero();
      r.reserve (rend - rbegin);
      r.getC().assign (rbegin, rend);

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
            = findLc (a, uu, polyReduce (a, uu, v, q));
      const typename type::CDownI end   = uu.toA0();

      r.makeZero();
      r.reserve (end - start);
      r.getC().assign (start, end);
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

   u.getC().erase (polyReduceS (a, u, v), u.toA0());
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

   u.getC().erase (polyReduceS (a, u, v), u.toA0());
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

   const P::DownI begin
      = Private::findLc (*x.a, uu, Private::polyReduce (*x.a, uu, *v));
   const P::DownI end   = uu.toA0();
   
   reserve (end - begin);
   c.assign (begin, end);
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

   u.getC().erase (u.fromLc(), findLc (a, u, polyReduce (a, u, v)));
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
      
   for (typename type::DownI i = polyReduce (a, uu, v); i != rend; ++i)
   {
      if (! a.is0(*i))  return false;
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
      typename type::DownI qend = polyReduceS (a, result, v);
      typename type::DownI rend = result.toA0();
      
      for (typename type::DownI i = qend; i != rend; ++i)
      {
         if (! a.is0(*i))  return false;
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
         bool ok = findLc (a, uu, polyReduce (a, uu, v, temp)) == uu.toA0();
         if (ok)  destructiveAssign (result, temp);
         return ok;
      }
      else
      {
         result.makeZero();
         result.reserve (u.numCoeff() - v.numCoeff() + 1);

         return findLc (a, uu, polyReduce (a, uu, v, result)) == uu.toA0();
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
   if (is0(p) || isCanonical(p)) return a.one();

   unit_type l = p.lc();
   unit_type il = a.recip (l);

   p.lc() = a.one();

   const typename type::DownI end = p.toA0();
   for (typename type::DownI i = p.fromLc() + 1; i != end; ++i)
   {
      a.mulBy (*i, il);
   }

   return l;
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
   if (p.is0())  throwDivisionByZero();
   if (p.isConstant())  return p.lc();

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

   if (p.is0())  throwDivisionByZero();
   if (p.isConstant())  return p.lc();

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
   return squarefreeFactorImp (
      *static_cast<const PolynomialRing<A>*> (this),
      f, p, typename A::char_category());
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
HIntLib::Private::PRBA_Rational<A>::isPrime (const type& p) const
{
   const int degree = p.degree();

   if (degree <= 0)  return false;
   if (degree == 1)  return true;  //  bx+a  is always prime

   if (a.is0 (p.ct()))  return false; // ...+ cx^2 + bx + 0  is never prime

   // Make sure p is square-free

   if (! isSquarefree (p))  return false;
   
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
   
   int_type gcdNum = p.ct().numerator();
   int_type lcmDen = p.ct().denominator();

   for (int i = 1; i <= degree; ++i)
   {
      if (! a.is0 (p[i]))
      {
         gcdNum = genGcd (intAlg, gcdNum, p[i].numerator());
         lcmDen = genLcm (intAlg, lcmDen, p[i].denominator());
      }
   }
   
   intPoly_type intp; intp.reserve (degree + 1);

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
      pp.mulAndAdd (a.makeElement (intAlg.one(), den));

      // construct numerator

      for (int j = 0; j <= s; ++j)
      {
         if (i == j)  continue;

         type copy (pp);
         pp.mulByX();
         mulBy (copy, a.makeElement (intAlg.neg (nodes[j])));
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
   unsigned* n = x.n;
   
   for (;;)
   {
      P p = a->elementMonic ((*n)++);
      if (p.degree() >= 6)  throw Overflow();

      if (a->isPrime (p))
      {
         swap (p);
         return;
      }
   }
}


/********************  Polynomial Ring <real_tag>  ***************************/


/**
 *  isPrime()
 */

template<class A>
bool
HIntLib::Private::PRBA_Real<A>::isPrime (const type& p) const
{
   const int degree = p.degree();

   if (degree == 1)  return true;  //  bx+a  is always prime

   // 0 is zero
   // constants are units
   // polynomials of degree >= 3 as well as  ax^2+bx  can always be factored
   //    over the reals

   if (degree <= 0 || degree >= 3 || a.is0(p.ct()))  return false;

   // Calculate discriminant to determine if there is a (real) solution

   const typename A::type disc
      = a.sub (a.sqr (p[1]), a.mul (coeff_type(4), a.mul (p[2], p.ct())));
   
   return disc < 0.0 && ! a.is0 (disc);
} 


/**
 *  next()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::NextReal,A,T> x)
{
   const A* alg = x.a;
   unsigned n = (*(x.n))++;
   
   if (n & 1)
   {
      // create polynomial  x^2 + bx + c
      // with arbitrary  b, and  c  such that  b^2 < 4c

      reserve (3);
      unsigned ib, ic;
      unthread (n / 2, ib, ic);

      real b = alg->element (ib);  // arbitrary real number
      real c = HIntLib::sqr(b) / 4.0 + alg->element (ic * 2 + 1);
                   // (x > 0) + b^2 / 4

      mulAndAdd (alg->one()).mulAndAdd (b).mulAndAdd (c);
   }
   else
   {
      // create polynomial  x + c

      reserve (2);
      mulAndAdd (alg->one()).mulAndAdd (alg->element (n / 2));
   }
}


/********************  Polynomial Ring <complex_tag>  ************************/


/**
 *  next()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::NextComplex,A,T> x)
{
   const A* a = x.a;
   
   reserve (2);
   mulAndAdd (a->one()).mulAndAdd (a->element ((*(x.n))++));
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
bool
HIntLib::Private::PRBA_GF<A>::isPrimitive (const type& p) const
{
   const int deg = p.degree();

   if (deg <= 0 || ! isCanonical (p))  return false;

   // calculate  a0 = (-1)^deg * p[0]

   typename A::type a0 = p.ct();
   if (odd (deg))  a.negate (a0);

   // i)  a0  must be a primitive element

   if (! a.isPrimitiveElement (a0))  return false;
   
   // ii) x^r mod p  must be equal to a0

   const unsigned r = (powInt (a.size(), deg) - 1) / (a.size() - 1);
   const type polyx = x();
   const type xPowR
      = powerMod (*static_cast<const PolynomialRing<A>*>(this), polyx, r, p);

   if (! isUnit (xPowR) || xPowR.ct() != a0)  return false;

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
 *  See Knuth, TACP, vol 2, 4.6.2,
 *      Lidl/Niederreiter, Finite Fields, 4.1, and
 *      Cohen, CANT, Algorithm 3.4.10.
 */

template<class A>
bool
HIntLib::Private::PRBA_GF<A>::isPrime (const type& p) const
{
   const int degree = p.degree();
   const unsigned s = a.size();

   // check for the trivial cases

   if (degree == 1)  return true;       //  bx+a  is always prime
   if (degree <= 0 || a.is0(p.ct()))  return false;  // zero, units or ...+0

   // There is a large chance for a linear divisor.
   // So try linear factors first is  s  is not too large.

   if (s < 50 || (degree >= 5 && s < 80))
   {
      type q = x();

      for (unsigned i = 1; i < s; ++i)   // x+1, x+2,...
      {
         q.ct() = a.element (i);
         if (isDivisor (p, q))  return false;
      }

      if (degree <= 3)  return true;

      // if  deg  is large compared to  s, try all quadratic divisors

      if (s < 16 && degree > squareBeatsLinear[s])
      {
         type q = x(2);

         for (unsigned i = 1; i < s; ++i)   // x+1, x+2,...
         {
            q.ct() = a.element (i);

            for (unsigned j = 0; j < s; ++j)
            {
               q[1] = a.element (j);
               if (isDivisor (p, q))  return false;
            }
         }

         if (degree <= 5)  return true;
      }
   }

   // Make sure p is square-free

   if (! isSquarefree (p))  return false;

   // Berlekamp's Algorithm

   type factor =
      (int (s) < degree) ? x(s) :
      powerMod (*static_cast<const PolynomialRing<A>*>(this), x(), s, p);

   // The frist row contains only zeros, so we do not include it in the matrix

   Array<typename A::type> matrix ((degree - 1) * degree);
   type q = one();

   // Calculate the matrix  B - I

   for (int row = 0; row < degree - 1; ++row)
   {
      mulBy (q, factor);
      reduce (q, p);

      for (int col = 0; col <= q.degree(); ++col)
      {
         matrix [row * degree + col] = q[col];
      }

      for (int col = q.degree() + 1; col < degree; ++col)
      {
         matrix [row * degree + col] = coeff_type();
      }

      a.subFrom (matrix [row * degree + row + 1], a.one());
   }

   // degree - rank (of the untruncated matrix) gives the number of factors

   return isLinearlyIndependent (a, matrix.begin(), degree, degree - 1);
}


/**
 *  next()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::NextGF,A,T> x)
{
   const Private::PRBA_GF<A>* a = x.a;
   unsigned* n = x.n;
   
   for (;;)
   {
      P p = a->elementMonic ((*n)++);

      if (a->isPrime (p))
      {
         swap (p);
         return;
      }
   }
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
   template bool PRBA_Field<X >::isDivisor \
      (const type&, const type&) const; \
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

// Instantiate PolynomialRingBase<rational_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_RATIONAL(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_FIELD(X) \
   template Polynomial<X::type>::Polynomial \
      (Private::PG<Private::NextRational,X >);\
   namespace Private { \
   template bool PRBA_Rational<X >::isPrime (const type&) const; \
   }

// Instantiate PolynomialRingBase<real_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_REAL(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_FIELD(X) \
   template Polynomial<X::type>::Polynomial \
      (Private::PG<Private::NextReal,X >); \
   namespace Private { \
   template bool PRBA_Real<X >::isPrime (const type&) const; \
   }

// Instantiate PolynomialRingBase<complex_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_COMPLEX(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_FIELD(X) \
   template Polynomial<X::type>::Polynomial \
      (Private::PG<Private::NextComplex,X >);

// Instantiate PolynomialRingBase<gf_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_GF(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_FIELD(X) \
   template Polynomial<X::type>::Polynomial (Private::PG<Private::NextGF,X >);\
   namespace Private { \
   template bool PRBA_GF<X >::isPrimitive (const type&) const;\
   template bool PRBA_GF<X >::isPrime (const type&) const; \
   }

#endif

