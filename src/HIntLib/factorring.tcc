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

#include <HIntLib/factorring.h>

#include <HIntLib/exception.h>
#include <HIntLib/output.h>
#include <HIntLib/array.h>
#include <HIntLib/prime.h>


/***************  Factor Base  ***********************************************/


/**
 *  I/O
 */

template<typename A>
std::ostream&
HIntLib::Private::operator<< (std::ostream &o, const FactorBB<A> &a)
{
   Private::Printer ss (o);

   ss << a.arithmetic() << "/(";
   a.arithmetic().printShort (ss, a.modulus());
   ss << ")";

   return o;
}

template<typename A>
void
HIntLib::Private::FactorBB<A>::print (std::ostream &o, const type& x) const
{
   Private::Printer ss (o);

   A::printShort (ss, x);
   ss << " (";
   A::printShort (ss, m);
   ss << ')';
   A::printSuffix (ss);
}

template<typename A>
void
HIntLib::Private::FactorBB<A>::printSuffix (std::ostream &o) const
{
   Private::Printer ss (o);

   ss << '(';
   A::printShort (ss, m);
   ss << ')';
   A::printSuffix (ss);
}


/**
 *  throwIfZero()
 */

template<class A>
void
HIntLib::Private::FactorBB<A>::throwIfZero (const type& x) const
{
   if (is0 (x)) throw DivisionByZero();
}


/**
 *  element()
 */

template<class A>
typename HIntLib::Private::FactorPolyB<A,HIntLib::infinite_tag>::type
         HIntLib::Private::FactorPolyB<A,HIntLib::infinite_tag>::
element (unsigned index) const
{
   unsigned num = m.degree();

   Array<unsigned> indices (num);

   unthreadn (index, indices.begin(), num);

   while (num && indices[num - 1] == 0)  --num;

   if (! num)  return type();

   type p = A::x (num - 1);

   for (unsigned i = 0; i < num; ++i)
   {
      p[i] = A::getCoeffAlgebra().element (indices[i]); 
   }

   return p;
}


/**
 *  index()
 */

template<class A>
unsigned
HIntLib::Private::FactorPolyB<A,HIntLib::infinite_tag>::
index (const type& p) const
{
   if (is0 (p))  return 0;

   const unsigned total = m.degree();
   const unsigned num = p.degree() + 1;

   Array<unsigned> indices (total);

   for (unsigned i = 0; i < num; ++i)
   {
      indices[i] = A::getCoeffAlgebra().index(p[i]);
   }

   for (unsigned i = num; i < total; ++i)  indices[i] = 0;

   return threadn (indices.begin(), total);
}


/**
 *  recip()
 */

template<class A>
void
HIntLib::Private::FactorB<A,HIntLib::euclidean_tag,HIntLib::polynomial_tag>::
recipImp (const type& x, type& r) const
{
   throwIfZero (x);
   type g = genGcd (arithmetic(), x, m, r);
   if (! A::isUnit (g))  throw InvalidModularFieldSize (0);
   A::mulByUnit (r, A::unitRecip(A::toUnit (g)));
}

template<class A>
void
HIntLib::Private::FactorB<A,HIntLib::integer_tag,HIntLib::nopolynomial_tag>::
recipImp (const type& x, type& r) const
{
   throwIfZero (x);
   type g = genGcd (arithmetic(), x, m, r);
   if (! A::is1 (g))  throw InvalidModularFieldSize (int (m));
   if (r < 0)  A::addTo (r, m);
}


/*****************  Factor Ring  *********************************************/

/**
 *  Constructor
 */

template<typename A>
HIntLib::Private::FactorRingB<A,HIntLib::integer_tag,HIntLib::nopolynomial_tag>::
FactorRingB (const A& a, const typename A::type& modulus)
   : FactorB<A,integer_tag,nopolynomial_tag> (a, modulus)
{
   typedef typename A::Factorization F;
   F f;
   factor (f, modulus);

   nilradical = one();

   for (typename F::const_iterator i = f.begin(); i != f.end(); ++i)
   {
      A::mulBy (nilradical, i->first);
   }
}

template<typename A>
HIntLib::Private::FactorRingB<A,HIntLib::euclidean_tag,HIntLib::polynomial_tag>::
FactorRingB (const A& a, const typename A::type& modulus)
   : FactorB<A,euclidean_tag,polynomial_tag> (a, modulus)
{
   typedef typename A::Factorization F;
   F f;
   squarefreeFactor (f, modulus);

   nilradical = one();

   for (typename F::const_iterator i = f.begin(); i != f.end(); ++i)
   {
      A::mulBy (nilradical, i->first);
   }
}


/**
 *  numNilpotents()  --  for polynomials
 */

template<typename A>
unsigned
HIntLib::Private::FactorRingB<A,HIntLib::euclidean_tag,HIntLib::polynomial_tag>::
numNilpotents () const
{
   unsigned degDiff = m.degree() - nilradical.degree();

   if (degDiff == 0)  return 1;

   typedef typename A::coeff_algebra::size_category SC;
   return SC::finite ? powInt (getCoeffAlgebra().size(), degDiff) : 0;
}


/**
 *  additiveOrder()  --  for integers
 */

template<typename A>
unsigned
HIntLib::Private::FactorRingB<A,HIntLib::integer_tag,HIntLib::nopolynomial_tag>::
additiveOrder (const typename A::type& u) const
{
   return int (A::div (m, genGcd (arithmetic(), u, m)));
}


/***************  Factor Field  **********************************************/


/**
 *  order()
 */

template<class A>
unsigned HIntLib::FactorField<A>::orderImp (const type& u, gf_tag) const
{
   throwIfZero (u);

   //  Order of an element in a finite group.
   //  See Algorithm 1.4.3 in H.Kohen, CANT

   unsigned e = size() - 1;
   PrimeDivisors pd (e);
   unsigned exponent;

   while (unsigned prime = pd.next (exponent))
   {
      e /= powInt (prime, exponent);
      type g1 = power (u, e);

      while (! is1 (g1))
      {
         g1 = power (g1, prime);
         e *= prime;
      }
   }

   return e;
}

template<class A>
unsigned HIntLib::FactorField<A>::orderImp (const type& u, field_tag) const
{
   throwIfZero (u);

   // XXX This is a complete hack!
   // We dream up some upper bound and check if we get to 1 in less than this
   // number of multiplications. If not, we declare the order to be infinite.

   const unsigned max = 1u << m.degree();
   type x (u);
   unsigned n = 1;

   while (! is1 (x))
   {
      mulBy (x, u);
      ++n;
      if (n > max) return 0;
   }
   return n;
}


/**
 *  isPrimitiveElement()
 */

template<class A>
bool HIntLib::FactorField<A>::isPrimitiveElement (const type& u) const
{
   throwIfZero (u);

   if (! size())  throw InternalError (__FILE__, __LINE__);

   const unsigned q = size() - 1;
   if (Prime::test(q))  return ! is1(u);

   PrimeDivisors pd (q);

   while (unsigned prime = pd.next())
   {
      if (is1 (power(u, q / prime)))  return false;
   }
   return true;
}


// Instantiations

#define HINTLIB_INSTANTIATE_FACTORRING(X) \
   namespace Private { \
   template std::ostream& operator<< (std::ostream&, const FactorBB<X >&); \
   template void FactorBB<X >::print (std::ostream&, const type&) const; \
   template void FactorBB<X >::printSuffix (std::ostream&) const; \
   template void FactorBB<X >::throwIfZero(const type&) const; \
   template void FactorB<X,X::algebra_category,X::polynomial_category>:: \
            recipImp (const type&, type&) const; \
   } \
   template unsigned FactorField<X >::order(const type&) const; \
   template bool FactorField<X >::isPrimitiveElement(const type&) const;
   
#define HINTLIB_INSTANTIATE_FACTORRING_INT(X) \
   HINTLIB_INSTANTIATE_FACTORRING(X) \
   namespace Private { \
   template unsigned FactorRingB<X,integer_tag,nopolynomial_tag>:: \
      additiveOrder (const type&) const; \
   template FactorRingB<X,integer_tag,nopolynomial_tag>::FactorRingB \
      (const X&, const type&); \
   }

#define HINTLIB_INSTANTIATE_FACTORRING_PB(X) \
   HINTLIB_INSTANTIATE_FACTORRING(X) \
   namespace Private { \
   template FactorRingB<X,euclidean_tag,polynomial_tag>::FactorRingB \
      (const X&, const type&); \
   template unsigned FactorRingB<X,euclidean_tag,polynomial_tag>:: \
      numNilpotents() const; \
   }

#define HINTLIB_INSTANTIATE_FACTORRING_POL(X) \
   HINTLIB_INSTANTIATE_FACTORRING_PB(X)

#define HINTLIB_INSTANTIATE_FACTORRING_POLI(X) \
   HINTLIB_INSTANTIATE_FACTORRING_PB(X) \
   namespace Private { \
   template FactorPolyB<X,infinite_tag>::type \
            FactorPolyB<X,infinite_tag>::element (unsigned) const; \
   template unsigned \
            FactorPolyB<X,infinite_tag>::index (const type&) const; \
   }

