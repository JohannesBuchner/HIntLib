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

#ifndef HINTLIB_MODULAR_ARITHMETIC_H
#define HINTLIB_MODULAR_ARITHMETIC_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <iosfwd>

#include <HIntLib/algebra.h>
#include <HIntLib/prime.h>
#include <HIntLib/hlmath.h>

namespace HIntLib
{

namespace Private
{

/**
 *  Modular Arithmetic Base
 */

class ModularArithmeticBase
{
public:
   unsigned size() const     { return m; }
   unsigned modulus() const  { return m; }

   void printSuffix (std::ostream &) const;

protected:
   ModularArithmeticBase (unsigned modulus)
      : m(modulus), nilradical (calcNilradical(modulus)) {}

   static void invalidType();

   // some methods for fields

   void checkField (unsigned) const;
   unsigned recipImp (unsigned) const;
   unsigned orderImp (unsigned) const;

   // some methods for rings

   static unsigned calcNilradical (unsigned);

   void checkRing (unsigned) const;
   unsigned additiveOrderImp (unsigned) const;
   bool isUnitImp (unsigned) const;
   unsigned unitElementImp(unsigned) const;
   unsigned unitIndexImp(unsigned) const;
   unsigned unitRecipImp(unsigned) const;

   // I/O

   static void prnShort  (std::ostream &, unsigned);
          void prn       (std::ostream &, unsigned) const;

   const unsigned m;
   const unsigned nilradical;
};

std::ostream& operator<< (std::ostream &, const ModularArithmeticBase &);


/**
 * Modular Arithmetic
 */

template <class T>
class ModularArithmetic : public Private::ModularArithmeticBase
{
protected:
   ModularArithmetic (unsigned modulus)
      : Private::ModularArithmeticBase (modulus)
   {
      if (   std::numeric_limits<T>::is_signed
          || ! std::numeric_limits<T>::is_integer
          || 2 * std::numeric_limits<T>::digits
               > std::numeric_limits<unsigned>::digits)
      {
         invalidType();
      }
   } 

public:
   typedef T type;
   typedef nopolynomial_tag polynomial_category;
   typedef finite_tag size_category;

   T one()  const  { return T(1); }

   bool is0 (T a) const  { return a == 0; }
   bool is1 (T a) const  { return a == 1; }

   T element(unsigned i) const  { return T(i); }
   unsigned index (T x) const   { return x; }

   T add (const T& a, const T& b) const
      { unsigned t = unsigned(a) + unsigned(b); return     (t>=m) ? t - m : t; }
   void addTo (T& a,    const T& b) const
      { unsigned t = unsigned(a) + unsigned(b); a = (t>=m) ? t - m : t; }

   T neg (const T& a) const  { return a != 0  ?  m - a  :  0; }
   void negate (T& a) const    { if (a != 0) a = m - a; }

   T sub (const T& a, const T& b) const
      { int t = int(a) - int(b); return     (t<0) ? t+m : t; }
   void subFrom (T& a,  const T& b) const
      { int t = int(a) - int(b); a = (t<0) ? t+m : t; }

   T dbl (const T& a) const
      { unsigned t = unsigned (a) << 1; return t >= m ? t - m : t; }
   void times2 (T& a) const    { a = dbl(a); }
   T times (const T& a, unsigned k) const
      { if (k >= m) k %= m;  return (a * k) % m; }

   T mul (const T& a, const T& b) const
      { return (unsigned(a) * unsigned(b)) % m; }
   void mulBy (T& a,    const T& b) const
      { a = (unsigned(a) * unsigned(b)) % m; }

   T sqr (const T& a) const  { return (a * a) % m; }
   void square (T& a) const    { a = (a * a) % m; }

   void printShort (
         std::ostream &o, const T& a, PrintShortFlag = PrintShortFlag()) const
      { prnShort (o, a); }
   void print (std::ostream &o, const T& a) const  { prn (o, a); }
};

}  // namespace Private


/**
 *  Modular Arithmetic Ring
 *
 *  Modular arithmetic modulo an arbitrary number gives a ring which in general
 *  is no domain.
 */

template<typename T>
class ModularArithmeticRing : public Private::ModularArithmetic<T>
{
public:
   typedef ring_tag algebra_category;
   typedef char_non char_category;
   typedef zerodivisor_tag zerodivisor_category;
   typedef T type;
   typedef T unit_type;

   ModularArithmeticRing (unsigned modulus)
      : Private::ModularArithmetic<T> (modulus)
      { checkRing (std::numeric_limits<T>::max()); }

   bool isUnit (const T& a) const  { return isUnitImp (a); }
   bool isZerodivisor (const T& a) const { return ! isUnitImp (a); }
   bool isNilpotent (const T& a) const  { return a % nilradical == 0; }

   unsigned numUnits()      const  { return eulerPhi (m); }
   unsigned numNilpotents() const  { return m / nilradical; }

   unsigned unitIndex (const T& a) const  { return unitIndexImp (a); }
   T unitElement (unsigned n)      const  { return unitElementImp (n); }

   static T toUnit   (const T& a)  { return a; }
   static T fromUnit (const T& a)  { return a; }

   T unitRecip (const T& a) const  { return unitRecipImp(a); }
   T mulUnit (const T& a, const T& b) const  { return mul(a,b); }
   void mulByUnit (T& a,    const T& b) const  { mulBy(a,b); }

   unsigned additiveOrder (T a) const  { return additiveOrderImp (a); }

   T power (const T& a, unsigned k) const
      { return powerMod (unsigned(a), k, unsigned(m)); }
};


/**
 *  Modular Arithmetic Field
 *
 *  Modular arithmetic modulo a prime number gives a cyclic finite field.
 */

template<typename T>
class ModularArithmeticField : public Private::ModularArithmetic<T>
{
public:
   typedef cyclic_tag algebra_category;
   typedef char_prime char_category;
   typedef nozerodivisor_tag zerodivisor_category;
   typedef T type;

   ModularArithmeticField (unsigned modulus)
      : Private::ModularArithmetic<T> (modulus)
      { checkField(std::numeric_limits<T>::max()); }

   unsigned characteristic() const  { return m; }

   T recip (const T& x) const { return T (recipImp (x)); }
   void reciprocal (T& x) const { x = recipImp (x); }

   T div  (const T& a, const T& b) const  { return mul (a, recip(b)); }
   void divBy (T& a,  const T& b) const  { mulBy (a, recip(b)); }

   unsigned order (T a) const    { return orderImp (a); }
   bool isPrimitiveElement (T a) const  { return isPrimitiveRoot (a, m); }
   T power (const T& a, unsigned k) const
      { return powerModReduce (unsigned(a), k, unsigned(m)); }

   HINTLIB_TRIVIAL_CYCLIC_MEMBERS
};

}  // namespace HIntLib

#endif

