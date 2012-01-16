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
#include <HIntLib/hlmath.h>


/***************  Factor Base  ***********************************************/


/**
 *  ~FactorBB()
 *
 *  If this is inline, GCC 3.3 chocks.
 */

template<typename A>
HIntLib::Private::FactorBB<A>::~FactorBB() {}


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
         HIntLib::Private::FactorPolyB<A,HIntLib::infinite_tag>::element
      (unsigned index) const
{
   unsigned num = this->m.degree();

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

   const unsigned total = this->m.degree();
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
HIntLib::Private::FactorPoly<A>::recipImp (const type& x, type& r) const
{
   throwIfZero (x);
   type g = genGcd (this->arithmetic(), x, this->m, r);
   if (! A::isUnit (g))  throw InvalidModularFieldSize (0);
   A::mulByUnit (r, A::unitRecip(A::toUnit (g)));
}

template<class A>
void
HIntLib::Private::FactorInteger<A>::recipImp (const type& x, type& r) const
{
   throwIfZero (x);
   type g = genGcd (this->arithmetic(), x, this->m, r);
   if (! A::is1 (g))  throw InvalidModularFieldSize (int (this->m));
   if (r < 0)  A::addTo (r, this->m);
}


/*****************  Factor Ring  *********************************************/

/**
 *  Constructor
 */

// for integers

template<typename A>
HIntLib::Private::FactorRingB<A,HIntLib::integer_tag>::FactorRingB
      (const A& a, const typename A::type& modulus)
   : FactorInteger<A> (a, modulus)
{
   typedef typename A::Factorization F;
   F f;
   factor (f, modulus);

   nilradical = this->one();
   numberOfUnits = 1;

   for (typename F::const_iterator i = f.begin(); i != f.end(); ++i)
   {
      A::mulBy (nilradical, i->first);
      numberOfUnits = numberOfUnits * (int (i->first) - 1)
                                    * powInt (int (i->first), i->second - 1);
   }
}

// for polynomials over finite fields

template<typename A>
void
HIntLib::Private::FactorRingB<A,HIntLib::polyoverfield_tag>::init (finite_tag)
{
   typedef typename A::Factorization F;
   F f;
   factor (f, this->m);

   nilradical = this->one();
   numberOfUnits = 1;

   for (typename F::const_iterator i = f.begin(); i != f.end(); ++i)
   {
      A::mulBy (nilradical, i->first);
      unsigned nor = numOfRemainders (i->first);
      numberOfUnits = numberOfUnits * (nor - 1) * powInt (nor, i->second - 1);
   }
}

// for polynomials over infinite fields

template<typename A>
void
HIntLib::Private::FactorRingB<A,HIntLib::polyoverfield_tag>::init (infinite_tag)
{
   typedef typename A::Factorization F;
   F f;
   squarefreeFactor (f, this->m);  // squarefree factor sufficient for nilrad.

   nilradical = this->one();

   for (typename F::const_iterator i = f.begin(); i != f.end(); ++i)
   {
      A::mulBy (nilradical, i->first);
   }

   numberOfUnits = 0;
}


/**
 *  numNilpotents()  --  for polynomials
 */

template<typename A>
unsigned
HIntLib::Private::FactorRingB<A,HIntLib::polyoverfield_tag>
      ::numNilpotents () const
{
   unsigned degDiff = this->m.degree() - nilradical.degree();

   if (degDiff == 0)  return 1;

   typedef typename A::coeff_algebra::size_category SC;
   return SC::finite ? powInt (this->getCoeffAlgebra().size(), degDiff) : 0;
}


/**
 *  additiveOrder()  --  for integers
 */

template<typename A>
unsigned
HIntLib::Private::FactorRingB<A,HIntLib::integer_tag>::additiveOrder
      (const typename A::type& u) const
{
   return int (A::div (this->m, genGcd (this->arithmetic(), u, this->m)));
}


/***************  Factor Field  **********************************************/


/**
 *  order()
 */

template<class A>
unsigned
HIntLib::Private::FactorFieldPolyB<A,HIntLib::finite_tag>::
order (const type& u) const
{
   throwIfZero (u);

   if (facNMinus1.empty())  Prime::factor (facNMinus1, this->size() -1);

   //  Order of an element in a finite group.
   //  See Algorithm 1.4.3 in H.Kohen, CANT

   unsigned e = this->size() - 1;
   for (FacI i = facNMinus1.begin(); i != facNMinus1.end(); ++i)
   {
      const unsigned& prime = i->first;

      e /= powInt (prime, i->second);
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
unsigned
HIntLib::Private::FactorFieldB<A,HIntLib::integer_tag>::
order(const type& u) const
{
   throwIfZero (u);

   if (facNMinus1.empty())  A::factor (facNMinus1, A::sub (this->m, A::one()));

   //  Order of an element in a finite group.
   //  See Algorithm 1.4.3 in H.Kohen, CANT

   unsigned e = this->size() - 1;
   for (FacI i = facNMinus1.begin(); i != facNMinus1.end(); ++i)
   {
      const unsigned prime = i->first;

      e /= powInt (prime, i->second);
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
unsigned
HIntLib::Private::FactorFieldPolyB<A,HIntLib::infinite_tag>::
order(const type& u) const
{
   throwIfZero (u);

   // XXX This is a complete hack!
   // We dream up some upper bound and check if we get to 1 in less than this
   // number of multiplications. If not, we declare the order to be infinite.

   const unsigned max = 1u << this->m.degree();
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
 *
 *  This works only for finite fields!
 */

template<class A>
bool
HIntLib::Private::FactorFieldPolyB<A,HIntLib::finite_tag>::
isPrimitiveElement (const type& u) const
{
   if (is0 (u)) return false;

   const unsigned q = this->size() - 1;
   if (facNMinus1.empty())  Prime::factor (facNMinus1, q);

   for (FacI i = facNMinus1.begin(); i != facNMinus1.end(); ++i)
   {
      if (is1 (power(u, q / i->first)))  return false;
   }

   return true;
}

template<class A>
bool
HIntLib::Private::FactorFieldB<A,HIntLib::integer_tag>::
isPrimitiveElement (const type& u) const
{
   if (is0 (u)) return false;

   if (facNMinus1.empty())  A::factor (facNMinus1, A::sub (this->m, A::one()));

   const unsigned q = this->size() - 1;

   for (FacI i = facNMinus1.begin(); i != facNMinus1.end(); ++i)
   {
      if (is1 (power(u, q / i->first)))  return false;
   }

   return true;
}


// Instantiations


#define HINTLIB_INSTANTIATE_FACTORRING(X) \
   namespace Private { \
   template FactorBB<X >::~FactorBB(); \
   template std::ostream& operator<< (std::ostream&, const FactorBB<X >&); \
   template void FactorBB<X >::print (std::ostream&, const type&) const; \
   template void FactorBB<X >::printSuffix (std::ostream&) const; \
   template void FactorBB<X >::throwIfZero(const type&) const; \
   }
   
#define HINTLIB_INSTANTIATE_FACTORRING_INT(X) \
   HINTLIB_INSTANTIATE_FACTORRING(X) \
   namespace Private { \
   template void FactorInteger<X >::recipImp (const type&, type&) const; \
   template unsigned FactorRingB<X,integer_tag>:: \
      additiveOrder (const type&) const; \
   template FactorRingB<X,integer_tag>::FactorRingB \
      (const X&, const type&); \
   template unsigned FactorFieldB<X,integer_tag>:: \
      order(const type&) const; \
   template bool FactorFieldB<X,integer_tag>:: \
      isPrimitiveElement(const type&) const; \
   }

#define HINTLIB_INSTANTIATE_FACTORRING_PB(X) \
   HINTLIB_INSTANTIATE_FACTORRING(X) \
   namespace Private { \
   template void FactorPoly<X >::recipImp (const type&, type&) const; \
   template unsigned FactorRingB<X,polyoverfield_tag>:: \
      numNilpotents() const; \
   }

#define HINTLIB_INSTANTIATE_FACTORRING_POL(X) \
   HINTLIB_INSTANTIATE_FACTORRING_PB(X) \
   namespace Private { \
   template void FactorRingB<X,polyoverfield_tag>::init(finite_tag);\
   template unsigned FactorFieldPolyB<X,finite_tag>::order(const type&) const; \
   template bool FactorFieldPolyB<X,finite_tag>:: \
      isPrimitiveElement(const type&) const; \
   }

#define HINTLIB_INSTANTIATE_FACTORRING_POLI(X) \
   HINTLIB_INSTANTIATE_FACTORRING_PB(X) \
   namespace Private { \
   template void \
            FactorRingB<X,polyoverfield_tag>::init(infinite_tag);\
   template FactorPolyB<X,infinite_tag>::type \
            FactorPolyB<X,infinite_tag>::element (unsigned) const; \
   template unsigned \
            FactorPolyB<X,infinite_tag>::index (const type&) const; \
   template unsigned \
            FactorFieldPolyB<X,infinite_tag>::order(const type&) const; \
   }


