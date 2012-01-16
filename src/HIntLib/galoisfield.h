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

#ifndef HINTLIB_GALOIS_FIELD_H
#define HINTLIB_GALOIS_FIELD_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/polynomial.h>
#include <HIntLib/factorring.h>


namespace HIntLib
{

/**
 *  Galois Field
 */

template <typename B>
class GaloisField : public FactorField<PolynomialRing<FactorField<B> > >
{
public:
   GaloisField (unsigned base, unsigned exponent, bool xPrim = false);
   GaloisField (unsigned size, bool xPrim = false);

private:
   typedef FactorField<B> Field;
   typedef PolynomialRing<Field> Poly;
   typedef FactorField<Poly> ExtensionField;
   typedef typename Poly::type T;

   T findPoly (unsigned base, unsigned exponent, bool xPrim);
   T findPoly (unsigned size, bool xPrim);
};

}  // namespace HIntLib

#endif

