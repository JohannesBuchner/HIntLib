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

/**
 *  This Class models the polynomiols over the field {0,1} = GF_2
 *
 *  An integer type has to be supplied. The bits of this type are used to 
 *  store and manage the coefficients of the polynom, i.e.
 *  unsigned long int can store 32 coeffizients.
 *
 */

#ifndef HINTLIB_POLYNOMIAL2_H
#define HINTLIB_POLYNOMIAL2_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <iosfwd>

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_LIMITS
  #include <limits>
#else
  #include <HIntLib/fallback_limits.h>
#endif

#include <HIntLib/bitop.h>
#include <HIntLib/prime.h>
#include <HIntLib/exception.h>
#include <HIntLib/algebra.h>


namespace HIntLib
{

/**
 *  GF2
 */

class GF2 : public BitOpBasedAddition<unsigned char>
{
public:
   GF2 () {}

   typedef unsigned char type;
   typedef cyclic_tag algebra_category;
   typedef nopolynomial_tag polynomial_category;

   static unsigned size()  { return 2; }
   static unsigned modulus()  { return 2; }
   static unsigned characteristic()  { return 2; }
   static unsigned extensionDegree()  { return 1; }

   static type one()  { return type(1); }

   static bool is1 (type a)  { return a; }

   static type element  (unsigned i)  { return type(i); }
   static unsigned index  (type x)  { return x; }

   static type  mul   (const type& a, const type& b)  { return a & b; }
   static type& mulBy (      type& a, const type& b)  { return a &= b; }

   static type power (const type& a, unsigned)  { return a; }

   static type recip     (const type&)  { return 1; }

   static type  div   (const type& a, const type&)  { return a; }
   static type& divBy (      type& a, const type&)  { return a; }

   static unsigned additiveOrder (const type& a)  { return a + 1; }
   static unsigned order (const type&)  { return 1; }
   static bool isPrimitiveElement (const type&)  { return true; }
   
   static void print (std::ostream &, type);
   static void printShort (std::ostream &, type);
   static void printSuffix (std::ostream &);
};


/**
 *  GF2VectorSpace
 */

template<typename T>
class BitRef
{
public:
   BitRef (T* _ptr, unsigned bit) : ptr (_ptr), mask (T(1) << bit) {}

   operator unsigned char () const { return (*ptr & mask) != 0; }
   BitRef<T>&  operator= (unsigned char x)
      { if (x) *ptr |= mask; else *ptr &= ~mask; return *this; }
private:
   T* ptr;
   T  mask;
};

class GF2VectorSpaceBase
{
public:
   GF2 getScalarAlgebra() const  { return GF2(); }

   unsigned size() const  { return 1u << dim; }
   unsigned dimension() const  { return dim; }

   void printSuffix (std::ostream &o) const { GF2().printSuffix (o); }

protected:
   GF2VectorSpaceBase (unsigned _dim) : dim(_dim) {}

   const unsigned dim;
};

template<typename T>
class GF2VectorSpace : public GF2VectorSpaceBase,
                       public BitOpBasedAddition<T>
{
public:
   typedef vectorspace_tag algebra_category;
   typedef nopolynomial_tag polynomial_category;

   typedef GF2 scalar_algebra;
   typedef GF2::type scalar_type;
   typedef T type;
   typedef BitRef<T> scalar_reference;

   GF2VectorSpace (unsigned _dim) : GF2VectorSpaceBase (_dim)
   {
      if (dim < 1 || dim > unsigned (std::numeric_limits<T>::digits))
      {
         throw FIXME (__FILE__, __LINE__);
      }
   } 

   GF2VectorSpace ()
      : GF2VectorSpaceBase (std::numeric_limits<T>::digits) {}

   template<typename I> void toCoord (type a, I p) const
   {
      for (unsigned i = 0; i != dim; ++i,a>>=1)  *p++ = scalar_type (a & 1);
   }
   template<typename I> void fromCoord (type& a, I p) const
   {
      a = 0;
      p += dim;
      for (unsigned i = 0; i != dim; ++i) a = (a << 1) | (*--p);
   }

   scalar_type coord (const type& a, unsigned k) const
      { return (a >> k) & 1; }
   scalar_reference coord (type& a, unsigned k) const
      { return BitRef<T> (&a, k); }
   
   type element(unsigned i) const  { return type(i); }
   unsigned index (type x) const   { return x; }

   type  mul   (const type& a, scalar_type l) const  { return     l ? a : 0; }
   type& scale (      type& a, scalar_type l) const  { return a = l ? a : 0; }

   void print (std::ostream &, const type&) const;
   void printShort (std::ostream &, const type&) const;

   HINTLIB_TRIVIAL_DOMAIN_MEMBERS
};


/**
 *  Polynomial 2
 */

template<typename T>
class Polynomial2
{
private:
   typedef Polynomial2<T> P;

   static const T pattern0101 = T (0x5555555555555555ull);

   bool msb() const  { return d & (T(1) << std::numeric_limits<T>::digits-1); }

public:
   typedef unsigned char coeff_type;
   typedef BitRef<T>     coeff_reference;

   // Constructors

   Polynomial2 ()             : d (T()) {}   // Default is f(x) = 0
   Polynomial2 (const P& p)   : d (p.d) {}   // Copy constructor
   explicit Polynomial2 (T d) : d (d)   {}   // Conversion from int
   template<typename I> Polynomial2 (I a0, I an) : d (0)
   {
      while (an != a0)  d = (d << 1) | (*--an);
   }

   // Assignment

   P& operator= (const P& p)  { d = p.d; return *this; }

   // You can convert a polynomial to an int to get the coefficient bitmap

   operator T () { return d; }

   // Returns a certain coefficient (0 or 1).

   unsigned char operator[] (unsigned i) const { return (d >> i) & 1; }
   BitRef<T>     operator[] (unsigned i)  { return BitRef<T> (&d, i); }

   unsigned char lc() const  { return 1; }
   BitRef<T>     lc() { return BitRef<T> (&d, ms1(d)); }

   template<typename I>
   void toCoeff (I p) const
   {
      T x = d;
      while (x)  *p++ = (x & 1), x >>= 1;
   }

   // Degree of a polynomial. degree(f(x)=0) = -1

   int degree() const { return ms1 (d); }

   // Addition and Substraction

   P operator+ (const P& p) const  { return P (d ^ p.d); }
   P operator- (const P& p) const  { return P (d ^ p.d); }
   P operator+ () const   { return *this; }   // Unary plus
   P operator- () const   { return *this; }   // Unary minus
   P& operator+= (const P& p)  { d ^= p.d; return *this; }
   P& operator-= (const P& p)  { d ^= p.d; return *this; }

   // Multiplikation and Division

   static void div (const P u, const P v, P &q, P &r);

   P  operator*  (const P  p) const;
   P  operator/  (const P& p) const { P q, r; div (*this, p, q, r); return q; }
   P  operator%  (const P& p) const { P q, r; div (*this, p, q, r); return r; }
   P& operator*= (const P& p) { return *this = *this * p; }
   P& operator/= (const P& p) { P r; div (*this, p, *this, r); return *this; }
   P& operator%= (const P& p) { P q; div (*this, p, q, *this); return *this; }

   P& divByX () { d >>= 1; return *this; }
   P& mulByX () { if (msb()) throw Overflow(); d <<= 1; return *this; }
   P& divByX (unsigned k)  { d >>= k; return *this; }
   P& mulByX (unsigned k)  { d <<= k; return *this; }

   bool is0() const  { return d == 0; }
   bool is1() const  { return d == 1; }
   bool isX() const  { return d == 2; }

   P derivative () const { return P((d >> 1) & pattern0101); }

   // Constants

   static P one ()  { return P (1); }
   static P x ()    { return P (2); }
   static P xPow (unsigned k) { return P (T(1) << k); }

   // Test for certain classs

   bool isPrimitive() const;
   bool isPrime() const;
   bool isComposit() const { return d > 3 && ! isPrime(); }

   unsigned char evaluate (unsigned char x) const;

   void printShort (std::ostream&, char) const;

   template<typename TT>
   friend
   bool operator== (const Polynomial2<TT> p1, const Polynomial2<TT> p2);

private:
   T d;      // coefficient bitmap
};

// Comparation

template<typename T>
bool operator== (const Polynomial2<T> p1, const Polynomial2<T> p2)
{
   return p1.d == p2.d;
}


// Define other operations

template<typename T>
inline
Polynomial2<T> operator* (const Polynomial2<T> p, int k)
{
   return (k & 1) ? p : Polynomial2<T>();
}

template<typename T>
inline
Polynomial2<T> operator* (int k, const Polynomial2<T> p)
{
   return (k & 1) ? p : Polynomial2<T>();
}

template<typename T>
std::ostream& operator<< (std::ostream &o, const Polynomial2<T> p)
{
   p.printShort (o, 'x'); return o;
}


/**
 *  Polynomial 2 Ring
 */

class Polynomial2RingBase
{
public:
   struct unit_type {};

   static unsigned size()  { return 0; }
   static unsigned characteristic()  { return 2; }

   static unsigned numUnits()  { return 1; }
   static unit_type unitRecip (unit_type)  { return unit_type(); }
   static unit_type unitElement (unsigned) { return unit_type(); }
   static unsigned unitIndex (unit_type)  { return 0; }

   static void printSuffix (std::ostream &);

   char getVar() const  { return var; }

protected:
   Polynomial2RingBase (char _var) : var (_var) {}
   const char var;
};

inline
bool operator== (Polynomial2RingBase::unit_type, Polynomial2RingBase::unit_type)
{
   return true;
}

template <typename T>
class Polynomial2Ring : public Polynomial2RingBase
{
public:
   typedef Polynomial2<T> type;
   typedef euclidean_tag algebra_category;
   typedef polynomial_tag polynomial_category;
   typedef GF2 coeff_algebra;
   typedef GF2::type coeff_type;
   typedef BitRef<T> coeff_reference;

   Polynomial2Ring (char _var = 'x') : Polynomial2RingBase (_var) {}

   static coeff_algebra getCoeffAlgebra()  { return GF2(); }

   static type one()  { return type(1); }
   static type x (unsigned k = 1)  { return type (T(1) << k); }

   static bool is0 (type p)  { return p.is0(); }
   static bool is1 (type p)  { return p.is1(); }

   static type element  (unsigned i) { return type(i); }
   static unsigned index  (type x) { return T(x); }

   static type add (const type& a, const type& b)  { return a +  b; }
   static type& addTo (type& a,    const type& b)  { return a += b; }

   static type neg (const type& a)  { return a; }
   static type& negate (type& a)    { return a; }

   static type sub (const type& a, const type& b)  { return a -  b; }
   static type& subFrom (type& a,  const type& b)  { return a -= b; }

   static type mul (const type& a, const type& b)  { return a *  b; }
   static type& mulBy (type& a,    const type& b)  { return a *= b; }

   static unit_type  mulUnit   (unit_type,    unit_type) { return unit_type(); }
   static unit_type& mulByUnit (unit_type &a, unit_type) { return a; }
   static type  mulUnit   (const type& a, unit_type)  { return a; }
   static type& mulByUnit (      type& a, unit_type)  { return a; }

   static type quot (const type& a, const type& b)  { return a / b; }
   static type rem  (const type& a, const type& b)  { return a % b; }

   static void div (const type& a, const type& b, type& q, type& r)
      { type::div (a, b, q, r); }

   static type times (const type& p, unsigned k) { return odd(k) ? p : type(); }
   static type power (const type& p, unsigned k) { return powInt (p, k); }

   static bool isUnit (const type& p)        { return p.is1(); }
   static bool isPrime (const type& p)       { return p.isPrime(); }
   static bool isComposit (const type& p)    { return p.isComposit(); }
   static bool isCanonical (const type&)     { return true; }
   static bool isPrimitive (const type& p)   { return p.isPrimitive(); }

   static unit_type makeCanonical (type &)  { return unit_type(); }

   static type fromUnit (unit_type)  { return one(); }
   static unit_type toUnit (const type&)  { return unit_type(); }
   
   static unsigned numOfRemainders(const type& p)  { return 1u << p.degree(); }

   static unsigned norm (const type& p)  { return p.degree() + 1; }

   static type derivative (const type& p)  { return p.derivative(); }
   static coeff_type evaluate (const type& p, const coeff_type& x)
      { return p.evaluate (x); }

   void print (std::ostream &, const type&) const;
   void printShort (std::ostream &o, const type& p) const
      { p.printShort (o, var); }

   static unsigned order (const type& p)  { return is1(p) ? 1 : 0; }

   class PrimeGenerator
   {
   public:
      PrimeGenerator(const Polynomial2Ring<T> &) : n(2) {}
      PrimeGenerator(const Polynomial2Ring<T>& alg, const type& p) : n (p) {}
      type next();
   private:
      T n;
      PrimeGenerator ();
      PrimeGenerator& operator= (const PrimeGenerator&);
   };

   HINTLIB_TRIVIAL_DOMAIN_MEMBERS
}; 


std::ostream& operator<< (std::ostream &, const GF2 &);
std::ostream& operator<< (std::ostream &, const GF2VectorSpaceBase &);
std::ostream& operator<< (std::ostream &, const Polynomial2RingBase &);


} // namespace HIntLib

#endif

