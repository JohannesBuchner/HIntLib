/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/galoisfield.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/prime.h>
#include <HIntLib/exception.h>
#include <HIntLib/array.h>
#include <HIntLib/linearalgebragen.h>
#include <HIntLib/hlmath.h>


namespace L = HIntLib;


/**
 *  constructor
 */

template<typename B>
L::GaloisField<B>::GaloisField (unsigned base, unsigned exponent, bool xPrim)
   : FactorField<PolynomialRing<ModularArithmeticField<B> > >
      (Poly (Field (base)),
       findPoly (base, exponent, xPrim))
{}

template<typename B>
L::GaloisField<B>::GaloisField (unsigned size, bool xPrim)
   : FactorField<PolynomialRing<ModularArithmeticField<B> > >
      (Poly (Field (Prime::factorPrimePowerPrime (size))),
       findPoly (size, xPrim))
{}


/**
 *  find Poly ()
 *
 *  We find a primitive polynomials of degree _deg_+1, such that the polynomial
 *  p = x  becomes a primive root in the resulting modular arithmetic.
 *
 *  The algorithm is based on the what is done in the  Mathematica Package
 *  Algebra`FiniteFields` (IrreduciblePolynomial[] and TransformIrreducible[]).
 */

template<typename B>
typename L::GaloisField<B>::T
L::GaloisField<B>::findPoly (unsigned base, unsigned deg, bool xPrim)
{
   if (deg == 0)  throw GaloisFieldExponent ();

   Field field (base);
   Poly poly (field);

   if (deg == 1)  return poly.x();

   // Find an irreducible polynomial

   typename Poly::PrimeGenerator ig (poly, deg);
   T p = ig.next();

   if (! xPrim)  return p;

   // Construct the extension field

   ExtensionField ef (poly, p);

   // search for a primitive root

   T prim;

   for (unsigned n = base; ; ++n)
   {
      prim = poly.element (n);

      if (ef.isPrimitiveElement (prim))
      {
         // We want the polynomial x to be a primitive root of the extension
         // field.  If this is the case already, we are done

         if (n == base)  return p;

         break;
      }
   }

   // prim is a primitive element, therefore  1, prim, prim^2, prim^3,...  is a
   // base of the field-vector space.
   // We construct the matrix for the linear transformation mapping
   //   1, x, x^2, x^3,...,x^(deg-1)
   // onto
   //   1, prim, prim^2, prim^3,...,prim^(deg-1)

   Array<B> matrix (HIntLib::sqr (deg), B());

   T x = ef.one();

   for (unsigned col = 0; col < deg; ++col)
   {
      for (int row = 0; row <= x.degree(); ++row)
      {
         matrix[row * deg + col] = x[row];
      }

      ef.mulBy (x, prim);
   }

   // we need the inverse mapping, so we invert the matrix

   matrixInverse (field, matrix.begin(), deg);

   // our primitive polynomial is  x^deg - mapping^(-1)(prim^deg).
   // prim^deg is exactly what x contains right now.

   p[deg] = field.one();

   for (unsigned row = 0; row < deg; ++row)
   {
      B y = B();
      for (unsigned col = 0; col < deg; ++ col)
      {
         field.addTo (y, field.mul (matrix[row * deg + col], x[col]));
      }
      p[row] = field.neg (y);
   }

   return p;
}

template<typename B>
typename L::GaloisField<B>::T
L::GaloisField<B>::findPoly (unsigned size, bool xPrim)
{
   unsigned base, exponent;
   Prime::factorPrimePower (size, base, exponent);
   return findPoly (base, exponent, xPrim);
}


namespace HIntLib
{
   template class GaloisField<unsigned char>;
   template class GaloisField<unsigned short>;
}

