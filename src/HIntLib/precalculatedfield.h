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

#ifndef HINTLIB_PRECALCULATED_RING_H
#define HINTLIB_PRECALCULATED_RING_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <iosfwd>

#include <HIntLib/mymath.h>


namespace HIntLib
{

/**
 *  Precalculated Field Base
 */

class PrecalculatedFieldBase
{
public:
   unsigned size() const  { return s; }
   unsigned getCharacteristic() const  { return characteristic; }

protected:
   PrecalculatedFieldBase (unsigned _s) : s(_s), characteristic(_s) {}
   PrecalculatedFieldBase (const PrecalculatedFieldBase &f)
      : s (f.s), characteristic (f.characteristic) {}

   const unsigned s;
         unsigned characteristic;
};

std::ostream& operator<< (std::ostream &, const PrecalculatedFieldBase &);


/**
 *  Precalculated Field
 */

template <class T>
class PrecalculatedField : public PrecalculatedFieldBase
{
protected:
   T* shared;
   T* addTable;
   T* mulTable;
   T* negTable;
   T* recTable;

private:
   void copy();
   void destroy();
   void setPointers();

public:

   typedef T type;

   PrecalculatedField (unsigned _s);
   PrecalculatedField (const PrecalculatedField<T> &);

   ~PrecalculatedField()  { destroy(); }

   PrecalculatedField & operator= (const PrecalculatedField<T> &);

   bool operator== (const PrecalculatedField<T> &)  const;
   bool operator!= (const PrecalculatedField<T> &r) const
      { return ! operator==(r); }

   T zero() const  { return T(0); }
   T one()  const  { return T(1); }

   bool is0 (T a) const  { return a == T(0); }
   bool is1 (T a) const  { return a == T(1); }

   T element(unsigned i) const  { return T(i); }
   unsigned index (T x) const   { return unsigned (x); }

   T  add   (T  a, T b) const  { return addTable [a + s * b]; }
   T& addTo (T& a, T b) const  { return a = add (a,b); }

   T  neg    (T  a) const  { return negTable [a]; }
   T& negate (T& a) const  { return a = neg(a); }

   T  sub     (T  a, T b) const  { return add   (a, neg (b)); }
   T& subFrom (T& a, T b) const  { return addTo (a, neg (b)); }

   T  mul   (T  a, T b) const  { return mulTable [a + s * b]; }
   T& mulBy (T& a, T b) const  { return a = mul (a,b); }

   T  recip (T  a) const  { return recTable [a]; }

   T  div   (T  a, T b) const  { return mul   (a, recip (b)); }
   T& divBy (T& a, T b) const  { return mulBy (a, recip (b)); }

   T times (T a, unsigned k) const;
   T power (T a, unsigned k) const { return a ? powInt (a, k % (s-1)) : 0; }

   void dump (std::ostream &) const;
   
   // Change tables

   void setAdd (T, T, T);
   void setMul (T, T, T);
   void setNeg (T, T);
   void setRecip (T, T);
};


/**
 *  Produces a precalculated Galois Field
 */

template<class T>
void makeGaloisField (PrecalculatedField<T> &, unsigned base, unsigned exp);

template<class T>
class PrecalculatedGaloisField : public PrecalculatedField<T>
{
public:
   PrecalculatedGaloisField (unsigned prime, unsigned power);
   PrecalculatedGaloisField (unsigned size);
};


/**
 *  copy()
 *
 *  Initialize a Precalculated Field using some Arithmetic A
 */

template<class T, class A>
void copy (PrecalculatedField<T> & dest, const A src);

}  // namespace HIntLib

#endif

