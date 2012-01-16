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

#include <HIntLib/galoisfield.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;


/**
 *  constructor
 */

template<class B>
L::GaloisField<B>::GaloisField (unsigned base, unsigned exponent)
   : ModularArithField<PolynomialRing<ModularArithField<B> > >
      (Poly (Field (base)),
       findPoly (base, exponent))
{}

template<class B>
typename L::GaloisField<B>::T
L::GaloisField<B>::findPoly (unsigned base, unsigned exponent)
{
   if (exponent == 0)  throw GaloisFieldExponent ();

   Field field (base);
   Poly poly (field);
   T p (exponent + 1);

   p.mulAndAdd (1);

   for (unsigned i = 0; i < exponent; ++i)  p.mulAndAdd (0);

   for (unsigned nn = 0; ; ++nn) 
   {
      unsigned n = nn;
      unsigned k = 0;

      while (n)
      {
         p[k++] = n % base;
         n /= base;
      }

      if (poly.isPrime (p))  return p;
   }
}


namespace HIntLib
{
   template class GaloisField<unsigned char>;
   template class GaloisField<unsigned short>;
}

