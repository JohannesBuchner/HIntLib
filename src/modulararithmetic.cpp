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

#include <HIntLib/modulararithmetic.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/output.h>

namespace P = HIntLib::Private;


/**
 *  Constructor
 */

P::ModularArithmeticBase::ModularArithmeticBase
      (unsigned modulus, unsigned max, bool field)
   : m (modulus)
{
   if (m > max + 1 || m < 2)  throw InvalidModularFieldSize (m);

   if (field)
   {
      if (! Prime::test (m))  throw InvalidModularFieldSize (m);

      // calculate factorization of m-1 and set nilradical

      Prime::factor (fac, m - 1);
      nilradical = m;
   }
   else // ring
   {
      // calculate nilradical

      PrimeDivisors pd (m);

      nilradical = 1;
      while (unsigned p = pd.next())  nilradical *= p;
   }
}


/**
 *  I/O
 */

std::ostream&
P::operator<< (std::ostream &o, const ModularArithmeticBase &a)
{
   Printer ss (o);
#if HINTLIB_CHARACTER_SET == 4 && defined (HINTLIB_UTF8_SELECT)
   HINTLIB_UTF8_SELECT(ss.utf8(),
   ss << "\xe2\x84\xa4",  // DOUBLE-STRUCK CAPITAL Z
   ss << 'Z')
#else
   ss << 'Z';
#endif
   ss.subscript(a.modulus());
   return o;
}
#ifdef HINTLIB_BUILD_WCHAR
std::wostream&
P::operator<< (std::wostream &o, const ModularArithmeticBase &a)
{
   WPrinter ss (o);
#if HINTLIB_CHARACTER_SET == 4
   ss << L'\x2124';  // DOUBLE-STRUCK CAPITAL Z
#else
   ss << L'Z';
#endif
   ss.subscript(a.modulus());
   return o;
}
#endif

void
P::ModularArithmeticBase::printSuffix (std::ostream &o) const
{
   Printer ss (o);
   ss << '(' << m << ')';
}
#ifdef HINTLIB_BUILD_WCHAR
void
P::ModularArithmeticBase::printSuffix (std::wostream &o) const
{
   WPrinter ss (o);
   ss << L'(' << m << L')';
}
#endif

void
P::ModularArithmeticBase::prnShort (std::ostream &o, unsigned a)
{
   o << int (a);  // convert to int to make showpos work
}
#ifdef HINTLIB_BUILD_WCHAR
void
P::ModularArithmeticBase::prnShort (std::wostream &o, unsigned a)
{
   o << int (a);  // convert to int to make showpos work
}
#endif

void
P::ModularArithmeticBase::prn (std::ostream &o, unsigned a) const
{
   Printer ss (o);
   ss << int (a) << " (" << m << ')';
}
#ifdef HINTLIB_BUILD_WCHAR
void
P::ModularArithmeticBase::prn (std::wostream &o, unsigned a) const
{
   WPrinter ss (o);
   ss << int (a) << L" (" << m << L')';
}
#endif


/**
 *  invalidType()
 */

void
P::ModularArithmeticBase::invalidType ()
{
   throw InvalidType ("ModularArithmeticRing/Field");
}


/**
 * recip()  (for fields)
 */

unsigned
P::ModularArithmeticBase::recipImp (unsigned x) const
{
   throwDivisionByZero (x);

   if (x == 1 || x == m - 1)  return x;

   if (m < 8)
   {
      for (unsigned i = 2; ; ++i)
      {
         if ((i * x) % m == 1)  return i;
      }
   }
   else
   {
      int a;
      gcd (int (x), int(m), a);
      return a < 0 ? a + m : a;
   }
}


/**
 *  unitRecip()  (for rings)
 */

unsigned
P::ModularArithmeticBase::unitRecipImp (unsigned x) const
{
   throwDivisionByZero (x);

   if (x == 1 || x == m - 1)  return x;

   int a;
   int g = gcd (int(x), int(m), a);

   if (g != 1)  throw InvalidModularFieldSize (m);
   return a < 0 ? a + m : a;
}


/**
 *  isUnit()  (for rings)
 */

bool
P::ModularArithmeticBase::isUnitImp (unsigned x) const
{
   if (x == 0)  return false;
   if (x == 1 || x == m - 1)  return true;

   return isCoprime(x, m);
}


/**
 *  additiveOrder()  (for rings)
 */

unsigned
P::ModularArithmeticBase::additiveOrderImp (unsigned a) const
{
   return m / gcd (a, m);
}


/**
 *  order()  (for fields)
 *
 *  Order of an element in a finite group.
 *  See Algorithm 1.4.3 in H.Kohen, CANT
 */

unsigned
P::ModularArithmeticBase::orderImp (unsigned a) const
{
   throwDivisionByZero (a);

   unsigned e = m - 1;

   for (FacI i = fac.begin(); i != fac.end(); ++i)
   {
      const unsigned& prime = i-> first;

      e /= powInt (prime, i->second);
      unsigned g1 = powerMod (a, e, m);

      while (g1 != 1)
      {
         g1 = powerMod (g1, prime, m);
         e *= prime;
      }
   }

   return e;
}


/**
 *  isPrimitiveElement()  (for fields)
 *
 *  Order of an element in a finite group.
 *  See Algorithm 1.4.3 in H.Kohen, CANT
 */

bool
P::ModularArithmeticBase::isPrimitiveImp (unsigned a) const
{
   if (a == 0) return false;

   const unsigned q = m - 1;

   for (FacI i = fac.begin(); i != fac.end(); ++i)
   {
      if (powerMod(a, q / i->first, m) == 1)  return false;
   }

   return true;
}


/**
 *  unitElement()  (for rings)
 */

unsigned
P::ModularArithmeticBase::unitElementImp (unsigned n) const
{
   unsigned x = 1;
   for (unsigned i = 0; i < n; ++i)  while (! isUnitImp(++x)) ;
   return x;
}


/**
 *  unitIndex()  (for rings)
 */

unsigned
P::ModularArithmeticBase::unitIndexImp (unsigned x) const
{
   unsigned count = 0;
   for (unsigned i = 1; i <= x; ++i)  if (isUnitImp(i))  ++count;
   return count;
}

