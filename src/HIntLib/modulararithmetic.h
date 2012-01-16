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

#ifndef HINTLIB_MODULAR_ARITHMETIC_H
#define HINTLIB_MODULAR_ARITHMETIC_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
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
#ifdef HINTLIB_BUILD_WCHAR
   void printSuffix (std::wostream &) const;
#endif

protected:
   ModularArithmeticBase (unsigned modulus, unsigned max, bool field);
   ModularArithmeticBase (const ModularArithmeticBase& mab)
      : m(mab.m), nilradical(mab.nilradical), fac (mab.fac) {}

   static void invalidType();

   // some methods for fields

   unsigned recipImp (unsigned) const;
   unsigned orderImp (unsigned) const;
   bool isPrimitiveImp (unsigned) const;

   // some methods for rings

   unsigned additiveOrderImp (unsigned) const;
   bool isUnitImp (unsigned) const;
   unsigned unitElementImp(unsigned) const;
   unsigned unitIndexImp(unsigned) const;
   unsigned unitRecipImp(unsigned) const;

   // I/O

   static void prnShort  (std::ostream &, unsigned);
          void prn       (std::ostream &, unsigned) const;
#ifdef HINTLIB_BUILD_WCHAR
   static void prnShort  (std::wostream &, unsigned);
          void prn       (std::wostream &, unsigned) const;
#endif

   // state

   const unsigned m;
   unsigned nilradical;
   Prime::Factorization fac;

   typedef Prime::Factorization::const_iterator FacI;
};

std::ostream& operator<< (std::ostream &, const ModularArithmeticBase &);
#ifdef HINTLIB_BUILD_WCHAR
std::wostream& operator<< (std::wostream &, const ModularArithmeticBase &);
#endif


/**
 * Modular Arithmetic
 */

template <class T>
class ModularArithmetic : public Private::ModularArithmeticBase
{
protected:
   ModularArithmetic (unsigned modulus, bool field)
      : Private::ModularArithmeticBase (
            modulus, std::numeric_limits<T>::max(), field)
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
#ifdef HINTLIB_BUILD_WCHAR
   void printShort (
         std::wostream &o, const T& a, PrintShortFlag = PrintShortFlag()) const
      { prnShort (o, a); }
   void print (std::wostream &o, const T& a) const  { prn (o, a); }
#endif
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
   typedef char_none char_category;
   typedef zerodivisor_tag zerodivisor_category;
   typedef T type;
   typedef T unit_type;

   ModularArithmeticRing (unsigned modulus)
      : Private::ModularArithmetic<T> (modulus, false) {}

   bool isUnit (const T& a) const  { return isUnitImp (a); }
   bool isZerodivisor (const T& a) const { return ! isUnitImp (a); }
   bool isNilpotent (const T& a) const  { return a % this->nilradical == 0; }

   unsigned numUnits()      const  { return Prime::eulerPhi (this->m); }
   unsigned numNilpotents() const  { return this->m / this->nilradical; }

   unsigned unitIndex (const T& a) const  { return unitIndexImp (a); }
   T unitElement (unsigned n)      const  { return this->unitElementImp (n); }

   static T toUnit   (const T& a)  { return a; }
   static T fromUnit (const T& a)  { return a; }

   T unitRecip (const T& a) const  { return unitRecipImp(a); }
   T mulUnit (const T& a, const T& b) const  { return mul(a,b); }
   void mulByUnit (T& a,    const T& b) const  { mulBy(a,b); }

   unsigned additiveOrder (T a) const  { return additiveOrderImp (a); }

   T power (const T& a, unsigned k) const
      { return powerMod (unsigned(a), k, unsigned(this->m)); }
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
      : Private::ModularArithmetic<T> (modulus, true) {}

   unsigned characteristic() const  { return this->m; }

   T recip (const T& x) const { return T (recipImp (x)); }
   void reciprocal (T& x) const { x = recipImp (x); }

   T div  (const T& a, const T& b) const  { return mul (a, recip(b)); }
   void divBy (T& a,  const T& b) const  { mulBy (a, recip(b)); }

   unsigned order (const T& a) const    { return orderImp (a); }
   bool isPrimitiveElement (const T& a) const  { return isPrimitiveImp (a); }
   T power (const T& a, unsigned k) const
      { return powerModReduce (unsigned(a), k, unsigned(this->m)); }

   HINTLIB_TRIVIAL_CYCLIC_MEMBERS
};

}  // namespace HIntLib

#endif

