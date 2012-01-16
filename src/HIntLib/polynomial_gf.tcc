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

#ifndef HINTLIB_POLYNOMIAL_GF_TCC
#define HINTLIB_POLYNOMIAL_GF_TCC 1

#ifndef HINTLIB_POLYNOMIAL_FIELD_TCC
#include <HIntLib/polynomial_field.tcc>
#endif

#include <HIntLib/array.h>
#include <HIntLib/linearalgebragen.h>
#include <HIntLib/prime.h>
#include <HIntLib/hlmath.h>


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
   if (deg & 1)  this->a.negate (a0);

   // i)  a0  must be a primitive element

   if (! this->a.isPrimitiveElement (a0))  return false;
   
   // ii) x^r mod p  must be equal to a0

   const unsigned r = (powInt (this->a.size(), deg) - 1) / (this->a.size() - 1);
   const type polyx = this->x();
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
   typedef typename A::type coeff_type;

   const int degree = p.degree();
   const unsigned s = this->a.size();

   // check for the trivial cases

   if (degree == 1)  return true;       //  bx+a  is always prime
   if (degree <= 0 || this->a.is0(p.ct()))  return false;  // zero, units or .+0

   // There is a large chance for a linear divisor.
   // So try linear factors first is  s  is not too large.

   if (s < 50 || (degree >= 5 && s < 80))
   {
      for (unsigned i = 1; i < s; ++i)   // x+1, x+2,...
      {
         if (this->a.is0 (evaluate (p, this->a.element(i))))  return false;
      }

      if (degree <= 3)  return true;

      // if  deg  is large compared to  s, try all quadratic divisors

      if (s < 16 && degree > this->squareBeatsLinear[s])
      {
         type q = this->x(2);

         for (unsigned i = 1; i < s; ++i)   // x+1, x+2,...
         {
            q.ct() = this->a.element (i);

            for (unsigned j = 0; j < s; ++j)
            {
               q[1] = this->a.element (j);
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
      (int (s) < degree) ? this->x(s) :
      powerMod (*static_cast<const PolynomialRing<A>*>(this), this->x(), s, p);

   // The frist row contains only zeros, so we do not include it in the matrix

   Array<coeff_type> matrix ((degree - 1) * degree);
   type q = this->one();

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

      this->a.subFrom (matrix [row * degree + row + 1], this->a.one());
   }

   // degree - rank (of the untruncated matrix) gives the number of factors

   return isLinearlyIndependent (this->a, matrix.begin(), degree - 1, degree);
}


/**
 *  factor()
 *
 *  Berlekamp's algorithm for small p
 *
 *  See Cohen, CANT, Alg. 3.4.10 and Knuth, TACP, 4.6.2, Alg B
 */

template<typename A>
typename HIntLib::Private::PRBA_Field<A>::unit_type
HIntLib::Private::PRBA_GF<A>::factor (
      typename PRBA_Field<A>::Factorization& f, const type& p) const
{
   typedef typename A::type coeff_type;
   typedef typename PRBA_Field<A>::Factorization Factorization;

   // Perform squarefree factorization

   Factorization sff;
   coeff_type unit = squarefreeFactor (sff, p);

   const unsigned s = this->a.size();

   // Final splitting of all squarefree factors

   for (typename Factorization::iterator it = sff.begin();
        it != sff.end(); ++it)
   {
      type& p = it->first;

      // Check for trivial cases

      if (p.degree() == 1)   // degree 1 is always prime
      {
         f.push_back (std::make_pair (p, it->second));
         continue;
      }

      // There is a large chance for a linear factor.
      // So try linear factors first.

      for (unsigned i = 0; i < s; ++i)   // x+1, x+2,...
      {
         const coeff_type u = this->a.element (i);

         if (this->a.is0 (evaluate (p, u)))
         {
            f.push_back (std::make_pair (
                     type().mulAndAdd (this->a.one())
                           .mulAndAdd (this->a.neg (u)),
                     it->second));
            divByLinearFactor (p, u);
            if (p.degree() == 1)  break;
         }
      }

      // Everything with degree less than 4 must be prime now

      const int degree = p.degree();

      if (degree < 4)
      {
         f.push_back (std::make_pair (p, it->second));
         continue;
      }

      // Berlekamp's Algorithm

      type factor =
         (int (s) < degree) ? this->x(s) :
         powerMod (*static_cast<const PolynomialRing<A>*>(this),
                   this->x(), s, p);

      Array<coeff_type> matrix (2 * HIntLib::sqr (degree));
      type q = this->one();

      // Calculate the matrix  B - I

      for (int col = 0; col < degree; ++col)
      {
         for (int row = 0; row <= q.degree(); ++row)
         {
            matrix [row * degree + col] = q[row];
         }

         for (int row = q.degree() + 1; row < degree; ++row)
         {
            matrix [row * degree + col] = coeff_type();
         }

         this->a.subFrom (matrix [col * (degree + 1)], this->a.one());

         mulBy (q, factor);
         reduce (q, p);
      }

      // Calculate null space

      coeff_type* ns = &matrix [HIntLib::sqr(degree)];

      unsigned numFactors
         = nullSpace (this->a, matrix.begin(), degree, degree, ns);

      // If  p  does not factor, we are done

      if (numFactors == 1)
      {
         f.push_back (std::make_pair (p, it->second));
         continue;
      }

      // Split  p  until we have  numFactors  factors

      std::vector<type> polys;
      polys.reserve (numFactors + 1);
      polys.push_back (p);

      for (;;)
      {
         ns += degree;

         coeff_type* leadingCoeff = ns + degree - 1;
         while (this->a.is0(*leadingCoeff))  --leadingCoeff;
         type splitter (ns, leadingCoeff + 1);

         unsigned numPolys = polys.size();
         for (unsigned i = 0; i < numPolys; ++i)
         {
            type& currentPoly = polys[i];
            if (currentPoly.degree() < 4)  continue;

            bool success = false;

            for (unsigned element = 0; element < s; ++element)
            {
               splitter.ct() = this->a.element(element);

               const type split =
                  genGcd (*static_cast<const PolynomialRing<A>*>(this),
                          currentPoly, splitter);

               if (split.degree() >= 1 && split.degree() < currentPoly.degree())
               {
                  success = true;
                  polys.push_back (split);
               }
            }

            if (success)
            {
               destructiveAssign (currentPoly, polys.back());
               polys.pop_back();
               if (polys.size() == numFactors)  goto end;
            }
         }
      }
end:

      // copy irreducible factors to output Factorization

      for (typename std::vector<type>::iterator i = polys.begin();
           i != polys.end(); ++i)
      {
         // use makeCanonicla(), because the result of genGcd() used above is
         //    not guaranteed to have this property
         makeCanonical (*i);
         f.push_back (std::make_pair (*i, it->second));
      }
   }

   // return the unit (calculated by squarefreeFactor() )

   return unit;
}


/**
 *  next()
 */

template<typename T> template<typename A>
HIntLib::Polynomial<T>::Polynomial (Private::PG<Private::NextGF,A,T> x)
{
   const Private::PRBA_GF<A>* a = x.a;
   unsigned* np = x.np;
   
   for (;;)
   {
      P p = a->elementMonic ((*np)++);

      if (a->isPrime (p))
      {
         swap (p);
         return;
      }
   }
}


/*********************  Instantiations  **************************************/

// Instantiate PolynomialRingBase<gf_tag>

#define HINTLIB_INSTANTIATE_POLYNOMIALRING_GF(X) \
   HINTLIB_INSTANTIATE_POLYNOMIALRING_FIELD(X) \
   template Polynomial<X::type>::Polynomial (Private::PG<Private::NextGF,X >);\
   namespace Private { \
   template bool PRBA_GF<X >::isPrimitive (const type&) const;\
   template bool PRBA_GF<X >::isPrime (const type&) const; \
   template PRBA_GF<X >::unit_type \
            PRBA_GF<X >::factor (Factorization&, const type&) const; \
   }

#endif

