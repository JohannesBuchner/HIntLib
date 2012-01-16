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
 *  A Polynomial over an arbitrary ring or field.
 */

#ifndef HINTLIB_POLYNOMIAL_H
#define HINTLIB_POLYNOMIAL_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <vector>
#include <utility>
#include <iosfwd>

#include <HIntLib/algebra.h>
#include <HIntLib/hlmath.h>

namespace HIntLib
{

/*
 *  Forward declarations
 */

template<typename T> class Polynomial;
namespace Private
{
   class PRB;
   template<class A> class PRBA;
   template<class A> class PRBA_Ring;
   template<class A> class PRBA_Domain;
   template<class A> class PRBA_UFD;
   template<class A> class PRBA_Field;
   template<class A> class PRBA_Rational;
   template<class A> class PRBA_Real;
   template<class A> class PRBA_Complex;
   template<class A> class PRBA_GF;

/*
 * Stub objects
 */

#define HINTLIB_NEQ template<typename Gen2> bool operator!= (const PG<Gen2,A,T>& g2) const  { return ! (*this == g2); }

#define HINTLIB_PG_A(Gen) \
   struct Gen {}; \
   template<typename A, typename T> struct PG<Gen,A,T> \
   { \
      HINTLIB_NEQ \
      const A* a; \
      PG (const A& _a) : a(&_a) {} \
   };
#define HINTLIB_PG_A_N(Gen) \
   struct Gen {}; \
   template<typename A, typename T> struct PG<Gen,A,T> \
   { \
      HINTLIB_NEQ \
      const A* a; \
      unsigned n; \
      PG (const A& _a, unsigned _n) : a(&_a), n(_n) {} \
   };
#define HINTLIB_PG_A_P(Gen) \
   struct Gen {}; \
   template<typename A, typename T> struct PG<Gen,A,T> \
   { \
      HINTLIB_NEQ \
      const A* a; \
      const Polynomial<T>* p; \
      PG (const A& _a, const Polynomial<T>& _p) : a(&_a), p(&_p) {} \
   };
#define HINTLIB_PG_A_P_P(Gen) \
   struct Gen {}; \
   template<typename A, typename T> struct PG<Gen,A,T> \
   { \
      HINTLIB_NEQ \
      const A* a; \
      const Polynomial<T>* p1; \
      const Polynomial<T>* p2; \
      PG (const A& _a, const Polynomial<T>& _p1, const Polynomial<T>& _p2) \
         : a(&_a), p1(&_p1), p2(&_p2) {} \
   };
#define HINTLIB_PG_A_P_N(Gen) \
   struct Gen {}; \
   template<typename A, typename T> struct PG<Gen,A,T> \
   { \
      HINTLIB_NEQ \
      const A* a; \
      const Polynomial<T>* p; \
      unsigned n; \
      PG (const A& _a, const Polynomial<T>& _p, unsigned _n) \
         : a(&_a), p(&_p), n(_n) {} \
   };

   /**
    *  PG  -  Polynomial generator
    *
    *  Various specializations (tages with the classes defined above) are
    *  returned by the methods of PolynomialRing wich are supposed to return
    *  type.
    *
    *  The constructors and operator= of Polynomial<> accepts these PGs as
    *  arguments and constructs the proper result in the target.
    */
      
   template<typename Gen, typename A, typename T = typename A::type>
   struct PG {};

   // Zero, One, Copy

   struct Zero {};
   template<typename A, typename T> struct PG<Zero,A,T> {};
   
   HINTLIB_PG_A (One)
   
   struct Copy {};
   template<typename A, typename T> struct PG<Copy,A,T>
   {
      HINTLIB_NEQ
      const Polynomial<T>* p;
      PG (const Polynomial<T>& _p) : p(&_p) {}
   };

   // X, Xn

   HINTLIB_PG_A (X)
   HINTLIB_PG_A_N (Xn)

   // Element, ElementMonic

   template<typename X> struct Element {};
   template<typename A, typename T, typename X> struct PG<Element<X>,A,T>
   {
      HINTLIB_NEQ
      const A* a;
      unsigned n;
      PG (const A& _a, unsigned _n) : a(&_a), n(_n) {}
   };

   template<typename X> struct ElementMonic {};
   template<typename A, typename T,typename X> struct PG<ElementMonic<X>,A,T>
   {
      HINTLIB_NEQ
      const A* a;
      unsigned n;
      PG (const A& _a, unsigned _n) : a(&_a), n(_n) {}
   };

   // Neg, Dbl, Times

   template<typename X> struct Neg {};
   template<typename A, typename T, typename X> struct PG<Neg<X>,A,T>
   {
      HINTLIB_NEQ
      const A* a;
      const Polynomial<T>* p;
      PG (const A& _a, const Polynomial<T>& _p) : a(&_a), p(&_p) {}
   };

   template<typename A, typename T> struct PG<Neg<char_two>,A,T>
      : public PG<Copy,A,T>
   {
      HINTLIB_NEQ
      PG (const A&, const Polynomial<T>& _p) : PG<Copy,A,T> (_p) {}
   };

   template<typename X> struct Dbl{};
   template<typename A, typename T, typename X> struct PG<Dbl<X>,A,T>
   {
      HINTLIB_NEQ
      const A* a;
      const Polynomial<T>* p;
      PG (const A& _a, const Polynomial<T>& _p) : a(&_a), p(&_p) {}
   };

   template<typename A, typename T> struct PG<Dbl<char_two>,A,T>
      : public PG<Zero,A,T>
   {
      PG (const A&, const Polynomial<T>&) {}
   };

   template<typename X> struct Times{};
   template<typename A, typename T, typename X> struct PG<Times<X>,A,T>
   {
      HINTLIB_NEQ
      const A* a;
      const Polynomial<T>* p;
      unsigned n;
      PG (const A& _a, const Polynomial<T>& _p, unsigned _n)
         : a(&_a), p(&_p), n(_n) {}
   };

   template<typename A, typename T> struct PG<Times<char_two>,A,T>
   {
      HINTLIB_NEQ
      const Polynomial<T>* p;
      unsigned n;
      PG (const A&, const Polynomial<T>& _p, unsigned _n)
         : p(&_p), n(_n) {}
   };

   // Add

   HINTLIB_PG_A_P_P (Add)
   HINTLIB_PG_A_P_P (Sub)

   // sqr(), square(), power()

   HINTLIB_PG_A_P (Derivative)

   template<typename X, typename Y> struct Sqr {};
   template<typename A, typename T, typename X, typename Y>
   struct PG<Sqr<X,Y>,A,T>
   {
      HINTLIB_NEQ
      const A* a;
      const Polynomial<T>* p;
      PG (const A& _a, const Polynomial<T>& _p) : a(&_a), p(&_p) {}
   };

   HINTLIB_PG_A_P_N (Power)

   // Mul

   template<typename X> struct Mul {};
   template<typename A, typename T, typename X> struct PG<Mul<X>,A,T>
   {
      HINTLIB_NEQ
      const A* a;
      const Polynomial<T>* p1;
      const Polynomial<T>* p2;
      PG (const A& _a, const Polynomial<T>& _p1, const Polynomial<T>& _p2)
         : a(&_a), p1(&_p1), p2(&_p2) {}
   };

   struct MulCoeff {};
   template<typename A, typename T> struct PG<MulCoeff,A,T>
   {
      HINTLIB_NEQ
      const A* a;
      const Polynomial<T>* p;
      const T* u;
      PG (const A& _a, const Polynomial<T>& _p, const T& _u)
         : a(&_a), p(&_p), u(&_u) {}
   };

   // Ring stuff

   struct UnitRecipR {};
   template<typename A, typename T> struct PG<UnitRecipR,A,T>
   {
      HINTLIB_NEQ
      const PRBA_Ring<A>* a;
      const Polynomial<T>* p;
      PG (const PRBA_Ring<A>* _a, const Polynomial<T>& _p)
         : a(_a), p(&_p) {}
   };

   // Domain stuff

   HINTLIB_PG_A_P_P (DivDomain)

   struct MulUnitD {};
   template<typename A, typename T> struct PG<MulUnitD,A,T>
   {
      HINTLIB_NEQ
      const A* a;
      const Polynomial<T>* p;
      const typename A::unit_type* u;
      PG (const A& _a, const Polynomial<T>& _p, const typename A::unit_type& _u)
         : a(&_a), p(&_p), u(&_u) {}
   };

   struct FromUnitD {};
   template<typename A, typename T> struct PG<FromUnitD,A,T>
   {
      HINTLIB_NEQ
      const A* a;
      const typename A::unit_type* u;
      PG (const A& _a, const typename A::unit_type& _u) : a(&_a), u(&_u) {}
   };

   // Field stuff

   HINTLIB_PG_A_P_P (Quot)
   HINTLIB_PG_A_P_P (Rem)

   struct FromCoeff {};
   template<typename A, typename T> struct PG<FromCoeff,A,T>
   {
      HINTLIB_NEQ
      const T* u;
      PG (const T& _u) : u(&_u) {}
   };

   // Next

   struct NextGF {};
   template<typename A, typename T> struct PG<NextGF,A,T>
   {
      HINTLIB_NEQ
      const PRBA_GF<A>* a;
      unsigned* n;
      PG (const PRBA_GF<A>& _a, unsigned& _n)
         : a(&_a), n(&_n) {}
   };

   struct NextRational {};
   template<typename A, typename T> struct PG<NextRational,A,T>
   {
      HINTLIB_NEQ
      const PRBA_Rational<A>* a;
      unsigned* n;
      PG (const PRBA_Rational<A>& _a, unsigned& _n)
         : a(&_a), n(&_n) {}
   };

   struct NextReal {};
   template<typename A, typename T> struct PG<NextReal,A,T>
   {
      HINTLIB_NEQ
      const A* a;
      unsigned* n;
      PG (const A& _a, unsigned& _n) : a(&_a), n(&_n) {}
   };

   struct NextComplex {};
   template<typename A, typename T> struct PG<NextComplex,A,T>
   {
      HINTLIB_NEQ
      const A* a;
      unsigned* n;
      PG (const A& _a, unsigned& _n) : a(&_a), n(&_n) {}
   };

#undef HINTLIB_PG_A
#undef HINTLIB_PG_A_N
#undef HINTLIB_PG_A_P
#undef HINTLIB_PG_A_P_N
#undef HINTLIB_PG_A_P_P
#undef HINTLIB_NEQ

}  // namespace Private


/***********  Polynomial *****************************************************/


/**
 *  Polynomial
 *
 *  Polynomial with an arbitrary number of terms, coefficients of type T
 */

template <typename T>
class Polynomial
{
private:
   typedef Polynomial<T> P;
   typedef std::vector<T> V;

   V c;

public:
   // the following stuff should be private

   P& mulAndAdd (const T& x)  { c.push_back (x); return *this; }
   P& reserve (unsigned num)  { c.reserve (num); return *this; }

// private:
   typedef typename V::iterator                DownI;
   typedef typename V::const_iterator         CDownI;

   CDownI fromLc() const { return c.begin(); }
    DownI fromLc()       { return c.begin(); }
   CDownI toA0()   const { return c.end(); }
    DownI toA0()         { return c.end(); }

   unsigned numCoeff () const  { return c.size(); }
   void makeZero ()  { c.clear(); }
   void erase (unsigned num)  { if (num)  c.erase(c.begin(), c.begin() + num); }
   V& getC()  { return c; }
   const V& getC() const  { return c; }

public:
   // typedefs

   typedef T coeff_type;
   typedef typename V::reference coeff_reference;
   typedef typename V::const_reference coeff_const_reference;

   // Constructors

   Polynomial () {}
   Polynomial (const P&);
   template<typename I> Polynomial (I a0, I an)
      : c (std::reverse_iterator<I> (an), std::reverse_iterator<I> (a0)) {}
   ~Polynomial();

   /*
    * Assignment and relatives
    */

   void swap (P& p)  { c.swap (p.c); }

   // Polynomials<> are assigned by using assignemt of vector<>

   P& operator= (const P&);

   // By default, PGs are assigned by constructing a temporary and swaping its
   // content.

   template<typename Gen, typename A>
   P& operator= (const Private::PG<Gen,A,T>& x)
      { P p (x); swap(p); return *this; }

   // For some PGs, an optimized version can be provided

   template<typename A>
   P& operator= (const Private::PG<Private::Zero,A,T>&)
      { makeZero(); return *this; }

   template<typename A>
   P& operator= (const Private::PG<Private::Copy,A,T>& x)
      { return operator= (*(x.p)); }

   /*
    * Constructors with creators
    */

   template<typename A> Polynomial (Private::PG<Private::One,A,T> x)
   {
      mulAndAdd (x.a->one());
   }
   template<typename A> Polynomial (Private::PG<Private::X,A,T>);
   template<typename A> Polynomial (Private::PG<Private::Xn,A,T>);
   template<typename A> Polynomial
      (Private::PG<Private::Element<finite_tag>,A,T>);
   template<typename A> Polynomial
      (Private::PG<Private::Element<infinite_tag>,A,T>);
   template<typename A> Polynomial
      (Private::PG<Private::ElementMonic<finite_tag>,A,T>);
   template<typename A> Polynomial
      (Private::PG<Private::ElementMonic<infinite_tag>,A,T>);
   // Neg
   template<typename A, typename X>
                        Polynomial (Private::PG<Private::Neg<X>,A,T>);
   template<typename A> Polynomial (Private::PG<Private::Neg<char_two>,A,T>);
   // Dbl
   template<typename A> Polynomial (Private::PG<Private::Dbl<char_non>,A,T>);
   template<typename A> Polynomial (Private::PG<Private::Dbl<char_zero>,A,T>);
   template<typename A> Polynomial (Private::PG<Private::Dbl<char_prime>,A,T>);
   template<typename A> Polynomial (Private::PG<Private::Dbl<char_two>,A,T>) {}
   // Times
   template<typename A> Polynomial (Private::PG<Private::Times<char_non>,A,T>);
   template<typename A> Polynomial (Private::PG<Private::Times<char_zero>,A,T>);
   template<typename A> Polynomial (Private::PG<Private::Times<char_prime>,A,T>);
   template<typename A> Polynomial (Private::PG<Private::Times<char_two>,A,T>);
   // Add
   template<typename A> Polynomial (Private::PG<Private::Add,A,T>);
   template<typename A> Polynomial (Private::PG<Private::Sub,A,T>);
   template<typename A> Polynomial (Private::PG<Private::Derivative,A,T>);
   // Sqr and Power
   template<typename A, typename Y> Polynomial
      (Private::PG<Private::Sqr<zerodivisor_tag,Y>,A,T>);
   template<typename A, typename Y> Polynomial
      (Private::PG<Private::Sqr<nozerodivisor_tag,Y>,A,T>);
   template<typename A> Polynomial
      (Private::PG<Private::Sqr<nozerodivisor_tag,char_two>,A,T>);

   template<typename A> Polynomial (Private::PG<Private::Power,A,T>);
   // Mul
   template<typename A> Polynomial
      (Private::PG<Private::Mul<nozerodivisor_tag>,A,T>);
   template<typename A> Polynomial
      (Private::PG<Private::Mul<zerodivisor_tag>,A,T>);
   template<typename A> Polynomial (Private::PG<Private::MulCoeff,A,T>);

   // Ring stuff
   template<typename A> Polynomial (Private::PG<Private::UnitRecipR,A,T>);

   // Domain stuff
   template<typename A> Polynomial (Private::PG<Private::DivDomain,A,T>);
   template<typename A> Polynomial (Private::PG<Private::MulUnitD,A,T>);
   template<typename A> Polynomial (Private::PG<Private::FromUnitD,A,T> x)
   {
      mulAndAdd (x.a->fromUnit (*x.u));
   }

   // Field stuff
   template<typename A> Polynomial (Private::PG<Private::Quot,A,T>);
   template<typename A> Polynomial (Private::PG<Private::Rem,A,T>);
   template<typename A> Polynomial (Private::PG<Private::FromCoeff,A,T> x)
   {
      mulAndAdd (*(x.u));
   }
   template<typename A> Polynomial (Private::PG<Private::NextGF,A,T>);
   template<typename A> Polynomial (Private::PG<Private::NextRational,A,T>);
   template<typename A> Polynomial (Private::PG<Private::NextComplex,A,T>);
   template<typename A> Polynomial (Private::PG<Private::NextReal,A,T>);

   /*
    * Normal member functions
    */

   int degree() const  { return int (c.size()) - 1; }
   bool is0() const  { return c.empty(); }
   bool isConstant() const  { return numCoeff() <= 1; }

   coeff_const_reference operator[] (unsigned i) const { return c[degree()-i]; }
   coeff_reference       operator[] (unsigned i)       { return c[degree()-i]; }
       
   coeff_const_reference lc() const  { return c.front(); }
   coeff_reference       lc()        { return c.front(); }

   coeff_const_reference ct() const  { return c.back(); }
   coeff_reference       ct()        { return c.back(); }

   template<typename I>
   void toCoeff (I a0) const  { std::copy (c.rbegin(), c.rend(), a0); }

   // Comparation

   P& divByX () { if (degree() >= 0)  c.pop_back(); return *this; }
   P& mulByX () { if (degree() >= 0)  c.push_back(coeff_type()); return *this; }
   P& divByX (unsigned);
   P& mulByX (unsigned);
};


/**
 *  operator==()  for Polynomial<> and PGs
 */

/*
 * P == P
 *
 * In general, we let vector<T> decide if Ps are equal. Only for Real<T> and
 * Complex<T> a different version is provided.
 *
 * There is no P != P, because that can easily be generated by std::rel_ops
 *
 * The STL of the SGI compiler has a bug: equal() uses operator!= (which is not
 * guaranteed to be available) instead of operator==().
 */

#ifdef HINTLIB_EQUAL_BUG
   template<typename T>
   bool operator== (const Polynomial<T>&, const Polynomial<T>&);
#else  // no equal bug
   template<typename T>
   inline
   bool operator== (const Polynomial<T>& p1, const Polynomial<T>& p2)
   {
      return p1.getC() == p2.getC();
   }
#endif


/*
 *  P == G
 *  G == P
 *  P != G
 *  G != P
 */

//  Case 2, 3 and 4 are reduced to case 1.

template<typename A, typename Gen>
inline
bool operator== (const Private::PG<Gen,A>& g,
                 const Polynomial<typename A::type>& p)
{
   return p == g;
}

template<typename A, typename Gen>
inline
bool operator!= (const Polynomial<typename A::type>& p,
                 const Private::PG<Gen,A>& g)
{
   return ! (p == g);
}

template<typename A, typename Gen>
inline
bool operator!= (const Private::PG<Gen,A>& g,
                 const Polynomial<typename A::type>& p)
{
   return ! (p == g);
}


// In case 1, by default, we create the polynomial and compare.
// Specializations can be provided for Zero and One

template<typename A, typename Gen>
inline
bool operator== (const Polynomial<typename A::type>& p,
                 const Private::PG<Gen,A>& g)
{
   return Polynomial<typename A::type>(g) == p;
}

template<typename A>
inline
bool operator== (const Polynomial<typename A::type>& p,
                 const Private::PG<Private::Zero,A>&)
{
   return p.is0();
}

template<typename A>
inline
bool operator== (const Polynomial<typename A::type>& p,
                 const Private::PG<Private::One,A>& x)
{
   return x.a->is1 (p);
}


/*
 *  G == G
 *
 *  If both arguments are Gs, we create one of them and are back at the case
 *  above.
 */

template<typename A, typename Gen1, typename Gen2>
inline
bool operator== (const Private::PG<Gen1,A>& g1, const Private::PG<Gen2,A>& g2)
{
   return Polynomial<typename A::type>(g1) == g2;
}

}  // namespace HIntLib


/**
 *  swap()  for Polynomial<>
 *
 *  We dump this into namespace std.
 *  This is not really standard conforming, but the only way that works.
 *  See various discussions on various C++ mailing lists.
 */

namespace std
{
   template<typename T>
   inline
   void swap (::HIntLib::Polynomial<T>& p1, ::HIntLib::Polynomial<T>& p2)
   {
      p1.swap (p2);
   }
}  // namespcae std

namespace HIntLib
{

   /**
    *  destructiveAssign()
    *
    *  The cheapest way of returning a polynomial is to swap src with dst.
    */

   template<typename T>
   inline
   void destructiveAssign (Polynomial<T>& dst, Polynomial<T>& src)
   {
      src.swap (dst);
   }


/***********  Polynomial Ring  ***********************************************/


namespace Private
{

/**
 *  Polynomial Ring Base  --  no template
 *
 *  Non-template base class of polynomial rings
 */

class PRB
{
public:
   static unsigned size()  { return 0; }

   void printVariable (std::ostream &) const;
   void printVariableWithBrackets (std::ostream &) const;
   void printVariablePow (std::ostream&, unsigned) const;

protected:
   explicit PRB (char _var) : var (_var) {}

   static const int squareBeatsLinear [];

private:
   const char var;
};


/**
 *  Polynomial Ring Base  --  A
 *
 *  All members which use only ring_base-properties of A
 */

template <class A>
class PRBA : public PRB
{
private:
   typedef typename A::size_category SC;
   typedef typename A::char_category CC;
   typedef typename A::zerodivisor_category ZC;

public:
   typedef polynomial_tag polynomial_category;
   typedef infinite_tag size_category;
   typedef CC char_category;
   typedef ZC zerodivisor_category;
   
   typedef A coeff_algebra;
   typedef typename A::type coeff_type;
   typedef Polynomial<coeff_type> type;
   typedef typename type::coeff_reference coeff_reference;

   const coeff_algebra& getCoeffAlgebra() const  { return a; }
   
   // construction

   PG<Xn, A> x(unsigned k) const  { return PG<Xn, A> (a, k + 1); }
   PG<X,  A> x()           const  { return PG<X,  A> (a); }
   PG<One,A> one()         const  { return PG<One,A> (a); }

   PG<Element<SC>,A> element (unsigned n) const
   {
      return PG<Element<SC>,A> (a, n);
   }
   PG<ElementMonic<SC>,A> elementMonic (unsigned n) const
   {
      return PG<ElementMonic<SC>,A> (a, n);
   }

private:
   unsigned indexImp (const type&, finite_tag) const;
   unsigned indexImp (const type&, infinite_tag) const;
public:
   unsigned index (const type& p) const  { return indexImp (p, SC()); }

   // basic properties

   bool is0 (const type& p) const  { return p.is0(); }
   bool is1 (const type&)  const;
   bool isMonic (const type&) const;

   // Negation, doubeling and times

private:
   void negateImp (type&, char_any) const;
   void negateImp (type&, char_two) const  {}

   void times2Imp (type&, char_non) const;
   void times2Imp (type&, char_prime) const;
   void times2Imp (type&, char_zero) const;
   void times2Imp (type& p, char_two) const  { p.makeZero(); }

public:
   void negate (type& p) const  { negateImp (p, CC()); }
   void times2 (type& p) const  { times2Imp (p, CC()); }

   PG<Neg<CC>,A> neg (const type& p) const { return PG<Neg<CC>,A> (a, p); }
   PG<Dbl<CC>,A> dbl (const type& p) const { return PG<Dbl<CC>,A> (a, p); }

   PG<Times<CC>,A> times (const type& p, unsigned n) const
      { return PG<Times<CC>,A> (a, p, n); }

   // Addition

   void addTo (type&, const type&) const;
   PG<Add,A> add (const type& p1, const type& p2) const
      { return PG<Add,A> (a, p1, p2); }

   void subFrom (type&, const type&) const;
   PG<Sub,A> sub (const type& p1, const type& p2) const
      { return PG<Sub,A> (a, p1, p2); }

   // sqr(), square(), power()

   PG<Sqr<ZC,CC>,A> sqr (const type& p) const
      { return PG<Sqr<ZC,CC>,A> (a, p); }
   PG<Power,A> power (const type& p, unsigned k) const
      { return PG<Power,A> (a, p, k); }
   void square (type& p) const  { p = sqr(p); }

   // mul()

   PG<Mul<ZC>,A> mul (const type& p1, const type& p2) const
      { return PG<Mul<ZC>,A> (a, p1, p2); }
   void mulBy (type& p1, const type& p2) const  { p1 = mul (p1, p2); }

   PG<MulCoeff,A> mul (const type& p, const coeff_type& c) const
      { return PG<MulCoeff,A> (a, p, c); }
   void mulBy (type&, const coeff_type&) const;

   // other stuff

   PG<Derivative,A> derivative (const type& p) const
      { return PG<Derivative,A> (a, p); }

   coeff_type evaluate (const type&, const coeff_type&) const;

   void print (std::ostream &, const type&) const;
   void printShort (
      std::ostream &, const type&, PrintShortFlag = PrintShortFlag()) const;
   void printSuffix (std::ostream &o) const  { a.printSuffix(o); }

protected:
   PRBA (const A& _a, char _var) : PRB (_var), a(_a) {}

   const A a;
};

template<class A>
std::ostream& operator<< (std::ostream &, const PRBA<A> &);


/**
 *  Polynomial Ring Base
 */

// ring_tag

template<class A>
class PRBA_Ring : public PRBA<A>
{
protected:
   PRBA_Ring (const A& _a, char _var) : PRBA<A> (_a, _var) {}

public:
   typedef typename PRBA<A>::type type;
   typedef type unit_type;
   typedef ring_tag algebra_category;
   
   unsigned numNilpotents() const
      { return a.numNilpotents() == 1 ? 1 : 0; }
   unsigned numUnits() const
      { return a.numNilpotents() == 1 ? a.numUnits() : 0; }

   const type&    fromUnit (const unit_type& p) const  { return p; }
   const unit_type& toUnit (const      type& p) const  { return p; }

   bool isNilpotent   (const type&) const;
   bool isUnit        (const type&) const;
   bool isZerodivisor (const type&) const;

   PG<UnitRecipR,A> unitRecip (const unit_type& u) const
      { return PG<UnitRecipR,A> (this, u); }

   PG<Mul<zerodivisor_tag>,A> mulUnit (const type& p1, const type& p2) const
      { return PG<Mul<zerodivisor_tag>,A> (a, p1, p2); }
   void mulByUnit (type& p1, const type& p2) const  { mulBy (p1, p2); }

   unsigned additiveOrder (const type&) const;
};


// domain_tag

template<class A>
class PRBA_Domain : public PRBA<A>
{
protected:
   PRBA_Domain (const A& _a, char _var) : PRBA<A> (_a, _var) {}

public:
   typedef typename PRBA<A>::type type;
   typedef typename A::unit_type unit_type;
   
   typedef domain_tag algebra_category;
   typedef noprimedetection_tag primedetection_category;

   bool isUnit (const type& p) const
      { return p.numCoeff() == 1 && a.isUnit (p.lc()); }

   unsigned numUnits() const  { return a.numUnits(); }
   unsigned unitIndex (const unit_type& u) const  { return a.unitIndex (u); }
   unit_type unitElement (unsigned i) const  { return a.unitElement (i); }

   PG<FromUnitD,A> fromUnit (const unit_type& u) const
      { return PG<FromUnitD,A> (a, u); }
   unit_type toUnit (const type& p) const  { return a.toUnit (p.lc()); }

   unit_type unitRecip (const unit_type& u) const  { return a.unitRecip (u); }

   PG<MulUnitD,A> mulUnit (const type& p, const unit_type& u) const
      { return PG<MulUnitD,A> (a, p, u); }
   void mulByUnit (type&, const unit_type&) const;

   unit_type mulUnit (const unit_type& u1, const unit_type& u2) const
      { return a.mulUnit (u1, u2); }
   void mulByUnit    (      unit_type& u1, const unit_type& u2) const
      { a.mulByUnit (u1, u2); }

   unsigned order(const type&) const;
   unsigned characteristic() const  { return a.characteristic(); }

   void divBy (type& u, const type& v) const;
   PG<DivDomain,A> div (const type& p1, const type& p2) const
      { return PG<DivDomain,A> (a, p1, p2); }

   bool isDivisor (const type&, const type&) const;
   bool isDivisor (const type&, const type&, type&) const;

   HINTLIB_TRIVIAL_DOMAIN_MEMBERS
};


// ufd_tag

template<class A>
class PRBA_UFD : public PRBA_Domain<A>
{
protected:
   PRBA_UFD (const A& _a, char _var) : PRBA_Domain<A> (_a, _var) {}

public:
   typedef typename PRBA<A>::type type;
   typedef typename PRBA_Domain<A>::unit_type unit_type;

   typedef ufd_tag algebra_category;

   bool isCanonical (const type&) const;
   unit_type makeCanonical (type&) const;
};


// field_tag

template<class A>
class PRBA_Field : public PRBA<A>
{
protected:
   PRBA_Field (const A& _a, char _var) : PRBA<A> (_a, _var) {}

private:
   typedef typename A::size_category SC;

public:
   typedef typename PRBA<A>::type type;
   typedef typename PRBA<A>::coeff_type unit_type;
   typedef std::vector<std::pair<type,unsigned> > Factorization;

   typedef euclidean_tag algebra_category;
   typedef noprimedetection_tag primedetection_category;

   unit_type unitRecip (const unit_type& u) const  { return a.recip (u); }
   void div (const type&, const type&, type&, type&) const;
   void reduce   (type&, const type&) const;
   void quotient (type&, const type&) const;
   void divBy    (type& u, const type& v) const  { quotient (u, v); }
   PG<Rem,A> rem (const type& p1, const type& p2) const
      { return PG<Rem,A> (a, p1, p2); }
   PG<Quot,A> quot (const type& p1, const type& p2) const
      { return PG<Quot,A> (a, p1, p2); }
   PG<Quot,A> div (const type& p1, const type& p2) const
      { return PG<Quot,A> (a, p1, p2); }

   bool isDivisor (const type&, const type&) const;
   bool isDivisor (const type&, const type&, type&) const;

   unsigned numOfRemainders (const type& p) const
      { return SC::finite ? powInt (a.size(), p.degree()) : 0; }

   unsigned norm (const type& p) const  { return p.degree() + 1; }

   bool isSquarefree (const type&) const;
   unit_type squarefreeFactor (Factorization&, const type&) const;
   bool isCanonical (const type& p) const  { return isMonic(p); }

   // units

   bool isUnit  (const type&p) const  { return p.numCoeff() == 1;}
   unsigned numUnits() const  { return SC::finite ? a.size() - 1 : 0; }

   unit_type unitElement (unsigned i) const  { return a.element(i + 1); }
   unsigned unitIndex (const unit_type& u) const  { return a.index (u) - 1; }
   unit_type makeCanonical (type &) const;

   PG<FromCoeff,A> fromUnit (const unit_type& u) const
      { return PG<FromCoeff,A> (u); }
   unit_type toUnit (const type& p) const  { return p.lc(); }

   PG<MulCoeff,A> mulUnit (const type& p, const unit_type& u) const
      { return PG<MulCoeff,A> (a, p, u); }
   void mulByUnit (type& p, const unit_type& u) const  { mulBy (p, u); }

   unit_type mulUnit (const unit_type& u1, const unit_type& u2) const
      { return a.mul (u1, u2); }
   void mulByUnit    (      unit_type& u1, const unit_type& u2) const
      { a.mulBy (u1, u2); }

   unsigned order(const type&) const;
   unsigned characteristic() const  { return a.characteristic(); }

   HINTLIB_TRIVIAL_DOMAIN_MEMBERS
};


// rational_tag

template<class A>
class PRBA_Rational : public PRBA_Field<A>
{
protected:
   PRBA_Rational (const A& _a, char _var) : PRBA_Field<A> (_a, _var) {}

public:
   typedef typename PRBA<A>::type type;
   typedef primedetection_tag primedetection_category;

   bool isPrime    (const type& p) const;
   bool isComposit (const type& p) const
      { return p.degree() > 1 && ! isPrime (p); }

   class PrimeGenerator;
};

template<class A>
class PRBA_Rational<A>::PrimeGenerator
{
public:
   PrimeGenerator(const PRBA_Rational<A> &_alg)
      : alg (_alg), n(2) {}

   PG<NextRational,A> next()  { return PG<NextRational,A> (alg, n); }

private:
   PrimeGenerator ();
   PrimeGenerator& operator= (const PrimeGenerator&);
   const PRBA_Rational<A> alg;
   unsigned n;
};


// real_tag

template<class A>
class PRBA_Real : public PRBA_Field<A>
{
protected:
   PRBA_Real (const A& _a, char _var) : PRBA_Field<A> (_a, _var) {}

public:
   typedef typename PRBA<A>::type type;
   typedef primedetection_tag primedetection_category;

   bool isPrime    (const type& p) const;
   bool isComposit (const type& p) const
      { return p.degree() > 1 && ! isPrime (p); }

   class PrimeGenerator;
};

template<class A>
class PRBA_Real<A>::PrimeGenerator
{
public:
   PrimeGenerator(const PRBA_Real<A> &_alg)
      : alg (_alg.getCoeffAlgebra()), n(0) {}

   PG<NextReal,A> next()  { return PG<NextReal,A> (alg, n); }

private:
   PrimeGenerator ();
   PrimeGenerator& operator= (const PrimeGenerator&);
   const A alg;
   unsigned n;
};


// complex_tag

template<class A>
class PRBA_Complex : public PRBA_Field<A>
{
protected:
   PRBA_Complex (const A& _a, char _var) : PRBA_Field<A> (_a, _var) {}

public:
   typedef typename PRBA<A>::type type;
   typedef primedetection_tag primedetection_category;

   bool isPrime    (const type& p) const  { return p.degree() == 1; };
   bool isComposit (const type& p) const  { return p.degree() >  1; };

   class PrimeGenerator;
};

template<class A>
class PRBA_Complex<A>::PrimeGenerator
{
public:
   PrimeGenerator(const PRBA_Complex<A> &_alg)
      : alg (_alg.getCoeffAlgebra()), n(0) {}

   PG<NextComplex,A> next()  { return PG<NextComplex,A> (alg, n); }

private:
   PrimeGenerator ();
   PrimeGenerator& operator= (const PrimeGenerator&);
   const A alg;
   unsigned n;
};


// gf_tag

unsigned funnySum (int, unsigned);

template<class A>
class PRBA_GF : public PRBA_Field<A>
{
protected:
   PRBA_GF (const A& _a, char _var) : PRBA_Field<A> (_a, _var) {}

public:
   typedef typename PRBA<A>::type type;
   typedef primedetection_tag primedetection_category;

   bool isPrimitive (const type&) const;
   bool isPrime (const type&) const;
   bool isComposit (const type& p) const
      { return p.degree() > 1 && ! isPrime (p); }

   class PrimeGenerator;
};

template<class A>
class PRBA_GF<A>::PrimeGenerator
{
public:
   PrimeGenerator(const PRBA_GF<A> &_alg)
      : alg (_alg), n(1) {}
   PrimeGenerator(const PRBA_GF<A> &_alg, int deg)
      : alg (_alg), n(funnySum (deg, alg.getCoeffAlgebra().size()))
      {}

   PG<NextGF,A> next()  { return PG<NextGF,A> (alg, n); }

private:
   const PRBA_GF<A> alg;
   unsigned n;
   PrimeGenerator ();
   PrimeGenerator& operator= (const PrimeGenerator&);
};


/**
 *  RingId
 *
 *  Determines the optimal implementation based on the algebra_category of the
 *  coefficient ring.
 */

template<typename A,typename TAG> struct RingId;

template<typename A> struct RingId<A,ring_tag>
   { typedef PRBA_Ring<A> base; };
template<typename A> struct RingId<A,domain_tag>
   { typedef PRBA_Domain<A> base; };
template<typename A> struct RingId<A,ufd_tag>
   { typedef PRBA_UFD<A> base; };
template<typename A> struct RingId<A,field_tag>
   { typedef PRBA_Field<A> base; };
template<typename A> struct RingId<A,rational_tag>
   { typedef PRBA_Rational<A> base; };
template<typename A> struct RingId<A,real_tag>
   { typedef PRBA_Real<A> base; };
template<typename A> struct RingId<A,complex_tag>
   { typedef PRBA_Complex<A> base; };
template<typename A> struct RingId<A,gf_tag>
   { typedef PRBA_GF<A> base; };

}  // namespace Private


/**
 *  Polynomial Ring
 *
 *  - Use Private::RingId to identify the proper implementation
 *  - Provide a default for the constructor
 */

template<class A>
class PolynomialRing
   : public Private::RingId<A,typename A::algebra_category::P>::base
{
public:
   explicit PolynomialRing (const A& _a, char _var = 'x')
      : Private::RingId<A,typename A::algebra_category::P>::base (_a, _var)
   {}
};

} // namespace HIntLib

#endif

