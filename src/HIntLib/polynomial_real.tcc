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

#ifndef HINTLIB_POLYNOMIAL_REAL_TCC
#define HINTLIB_POLYNOMIAL_REAL_TCC 1

#ifndef HINTLIB_POLYNOMIAL_FIELD_TCC
#include <HIntLib/polynomial_field.tcc>
#endif


/********************  Polynomial Ring <complex_tag>  ************************/


/**
 *  next()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::NextComplex,A,T> x)
{
   const A* a = x.a;
   
   reserve (2);
   mulAndAdd (a->one()).mulAndAdd (a->element ((*(x.np))++));
}


/**
 *  roots()
 *
 *  Factor a squarefree complex polynomial.
 *
 *  If cpx is false, the polynomial is assumed to be a real polynomial.
 *
 *  See Cohen, CANT, Alg. 3.6.6
 */

template<typename A>
void
HIntLib::Private::PRBA_Complex<A>::roots (
      typename PRBA_Field<A>::Factorization& r,
      const type& p, unsigned multiplicity, bool cpx) const
{
   typedef typename A::type coeff_type;

   const A& aa = this->a;

   // Shortcut the trivial case deg = 1

   if (p.degree() == 1)
   {
      r.push_back (std::make_pair (p, multiplicity));
      return;
   }

   typedef typename A::real_field RF;
   typedef typename RF::type real_type;

   const RF rf (aa.getRealField());

   // Step 1: Initialization

   const type pDeriv (derivative (p));
   type q (p);

   while (! q.isConstant())
   {
      const type qDeriv (derivative (q));

      // Step 2: Initialize root finding

      coeff_type x (1.3, 0.31415926);
      coeff_type y = evaluate (q, x);

      // Newton iteration

      for (;;)
      {
         unsigned tries = 0;

         // Determine step

         coeff_type dx = aa.div (y, evaluate (qDeriv, x));

         if (aa.is0 (dx))  break;

         for (;;)
         {
            // Calculate next approximation

            coeff_type xNew = aa.sub (x, dx);
            coeff_type yNew = evaluate (q, xNew);

            // if the new value is an improvement, use it

            if (abs(yNew) < abs(y))
            {
               x = xNew;
               y = yNew;
               break;
            }

            // if not, reduce step size by a factor of 4 (up to 20 times)

            if (++tries > 20)  throw 1;

            aa.mulBy (dx, coeff_type (.25));
         }
      }

      // Step 5: Polish root
      // It is important to use the original polynomial!

      aa.subFrom (x, aa.div (evaluate (p, x), evaluate (pDeriv, x)));
      aa.subFrom (x, aa.div (evaluate (p, x), evaluate (pDeriv, x)));

      // Step 6: Divide

      // If x is almost real, make it real

      if (rf.is0 (aa.im(x)))  x = coeff_type (aa.re(x));

      if (cpx || (rf.is0 (aa.im(x))))
      {
         divByLinearFactor (q, x);

         type root;
         root.mulAndAdd (aa.one()).mulAndAdd (aa.neg (x));

         r.push_back (std::make_pair (root, multiplicity));
      }
      else
      {
         // if  a + bi  is a complex root of a real polynomial p,
         // then  x^2 - 2ax + (a^2 + b^2)  is a factor of p.

         type root;
         root.reserve (3);
         root.mulAndAdd (aa.one())
             .mulAndAdd (coeff_type (rf.neg (rf.dbl (aa.re (x)))))
             .mulAndAdd (coeff_type (rf.add (rf.sqr (aa.re (x)),
                                             rf.sqr (aa.im (x)))));

         divBy (q, root);

         r.push_back (std::make_pair (root, multiplicity));
      }
   }
}


/**
 *  factor()
 *
 *  Factors a complex polynomial into linear factors.
 */

template<typename A>
typename HIntLib::Private::PRBA_Field<A>::unit_type
HIntLib::Private::PRBA_Complex<A>::factor (
      typename PRBA_Field<A>::Factorization& f, const type& p) const
{
   typedef typename PRBA_Field<A>::Factorization Factorization;

   // Perform squarefree factorization

   Factorization sff;
   typename A::type unit = squarefreeFactor (sff, p);

   // Final splitting of all squarefree factors

   for (typename Factorization::iterator i = sff.begin(); i != sff.end(); ++i)
   {
      roots (f, i->first, i->second, true);
   }

   return unit;
}


/********************  Polynomial Ring <real_tag>  ***************************/


/**
 *  isPrime()
 */

template<class A>
bool
HIntLib::Private::PRBA_Real<A>::isPrime (const type& p) const
{
   typedef typename A::type coeff_type;

   const int degree = p.degree();

   if (degree == 1)  return true;  //  bx+a  is always prime

   // 0 is zero
   // constants are units
   // polynomials of degree >= 3 as well as  ax^2+bx  can always be factored
   //    over the reals

   if (degree <= 0 || degree >= 3 || this->a.is0(p.ct()))  return false;

   // Calculate discriminant to determine if there is a (real) solution

   const coeff_type disc
      = this->a.sub (this->a.sqr (p[1]),
                     this->a.mul (coeff_type(4), this->a.mul (p[2], p.ct())));
   
   // a degree 2 polynomial is irreducible if disc is strictly negative

   return disc < coeff_type();
} 


/**
 *  next()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::NextReal,A,T> x)
{
   const A* alg = x.a;
   unsigned n = (*(x.np))++;
   
   if (n & 1)
   {
      // create polynomial  x^2 + bx + c
      // with arbitrary  b, and  c  such that  b^2 < 4c

      unsigned ib, ic;
      unthread (n / 2, ib, ic);

      const coeff_type b = alg->element (ib);  // arbitrary real number
      const coeff_type c = alg->add (alg->div (alg->sqr(b), coeff_type (4.)),
                                    alg->element (ic * 2 + 1));
                   // (x > 0) + b^2 / 4

      reserve (3);
      mulAndAdd (alg->one()).mulAndAdd (b).mulAndAdd (c);
   }
   else
   {
      // create polynomial  x + c

      reserve (2);
      mulAndAdd (alg->one()).mulAndAdd (alg->element (n / 2));
   }
}


/**
 *  factor()
 *
 *  Factors a polynomial over the reals using PRBA_Complex<>::roots()
 */

template<typename A>
typename HIntLib::Private::PRBA_Field<A>::unit_type
HIntLib::Private::PRBA_Real<A>::factor (
      typename PRBA_Field<A>::Factorization& f, const type& p) const
{
   typedef typename PRBA_Field<A>::Factorization Factorization;

   // Perform squarefree factorization

   Factorization sff;
   typename A::type unit = squarefreeFactor (sff, p);

   // Construct complex extension field and complex polynomial ring

   typedef typename A::complex_field CF;
   typedef typename CF::type complex_type;
   typedef PolynomialRing<CF> CP;
   typedef typename CP::type complex_poly_type;

   const CF cField (this->a.getComplexField());
   const CP cPolys (cField);

   // Final splitting of all squarefree factors

   for (typename Factorization::iterator sffi = sff.begin();
        sffi != sff.end(); ++sffi)
   {
      const type& p = sffi->first;

      // if p is prime, don't search for roots

      if (isPrime (p))
      {
         f.push_back (*sffi);
      }
      else  // degree > 1
      {
         const unsigned multiplicity = sffi->second;

         // make a complex copy of p

         complex_poly_type cp;
         for (typename type::CDownI i = p.fromLc(); i != p.toA0(); ++i)
         {
            cp.mulAndAdd (complex_type (*i));
         }

         // factor using complex Newton iteration

         typename CP::Factorization cf;
         cPolys.roots (cf, cp, multiplicity, false);

         // process all irreducible factors

         for (typename CP::Factorization::iterator it = cf.begin();
              it != cf.end(); ++it)
         {
            // convert back to a real polynomial

            type rp;
            for (typename complex_poly_type::CDownI i = it->first.fromLc();
                 i != it->first.toA0(); ++i)
            {
               rp.mulAndAdd (cField.re (*i));
            }

            // store real irreducible factor

            f.push_back (std::make_pair (rp, multiplicity));
         }
      }
   }

   return unit;
}


/*********************  Instantiations  **************************************/


// Instantiate PolynomialRingBase<real_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_REAL(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_FIELD(X) \
   template Polynomial<X::type>::Polynomial \
      (Private::PG<Private::NextReal,X >); \
   namespace Private { \
   template bool PRBA_Real<X >::isPrime (const type&) const; \
   template PRBA_Real<X >::unit_type \
            PRBA_Real<X >::factor (Factorization&, const type&) const; \
   }

// Instantiate PolynomialRingBase<complex_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_COMPLEX(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_FIELD(X) \
   template Polynomial<X::type>::Polynomial \
      (Private::PG<Private::NextComplex,X >); \
   namespace Private { \
   template PRBA_Complex<X >::unit_type \
            PRBA_Complex<X >::factor (Factorization&, const type&) const; \
   }

#endif

