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

#include <algorithm>
#include <iomanip>

#include <HIntLib/lookupfield.h>

#include <HIntLib/hlmath.h>
#include <HIntLib/array.h>
#include <HIntLib/output.h>
#include <HIntLib/exception.h>
#include <HIntLib/bitop.h>
#include <HIntLib/prime.h>

namespace L = HIntLib;

using std::setw;

/******************  Ref Counting Algebra  ***********************************/

inline unsigned* allocateMemory (size_t memory)
{
   unsigned* p = new unsigned [(memory - 1) / sizeof(unsigned) + 2];
   *p = 1;
   return p;
}

/**
 *  Constructor
 */

L::RefCountingAlgebra::RefCountingAlgebra (size_t memory)
   : refCount (allocateMemory (memory)),
     size (memory)
{}

L::RefCountingAlgebra::RefCountingAlgebra (size_t memory, unsigned* data)
   : refCount (data ? data : allocateMemory (memory)),
     size (memory)
{}


/**
 *  operator=
 */

L::RefCountingAlgebra &
L::RefCountingAlgebra::operator= (const RefCountingAlgebra &r)
{
   if (refCount != r.refCount)
   {
      destroy();
      refCount = r.refCount;
      size = r.size;
      setPointers();
      copy();
   }
   return *this;
}


/**
 *  write()
 */

void L::RefCountingAlgebra::write ()
{
   if (*refCount == 1)  return;    // if nobody else uses it, we are fine

   if (*refCount >= 2)  --*refCount;

   const char* oldData = charPtr();
   refCount = allocateMemory (size);
   setPointers();

   std::copy (oldData, oldData + size, charPtr());
}


/**
 *  clear()
 */

void L::RefCountingAlgebra::clear()
{
   write();
   std::fill (charPtr(), charPtr() + size, char (0));
}


/**
 *  destroy()
 */

void L::RefCountingAlgebra::destroy ()
{
   if (*refCount == 0)  return;
   if (! --*refCount)  delete[] refCount;
}


/**
 *  operator==
 */

bool L::RefCountingAlgebra::operator== (const RefCountingAlgebra &r) const
{
   return refCount == r.refCount
       || (   size == r.size
           && std::equal (charPtr(), charPtr() + size, r.charPtr()));
}


/******************  Printing  ***********************************************/


void
L::Private::printVectorSpaceName (std::ostream &o, unsigned b, unsigned d)
{
   Printer ss (o);
   ss << "GF" << b << '^' << d;
}


void
L::Private::printFieldName (std::ostream &o, unsigned b)
{
   Printer ss (o);
   ss << "GF" << b;
}

void
L::Private::printSuffix (std::ostream &o, unsigned b)
{
   Printer ss (o);
   ss << "(GF" << b << ')';
}

void
L::Private::printNumberSuffix (std::ostream &o, int x, unsigned b)
{
   Printer ss (o);
   ss << x << " (GF" << b << ')';
}

void
L::Private::printNumber (std::ostream &o, int x)
{
   o << x;
}


/******************  Lookup Field Base Base  *********************************/


/**
 *  operator=
 */

L::LookupFieldBB &
L::LookupFieldBB::operator= (const LookupFieldBB &r)
{
   RefCountingAlgebra::operator=(r);
   s = r.s;
   return *this;
}

/**
 *  operator==
 */

bool L::LookupFieldBB::operator== (const LookupFieldBB &r) const
{
   return size() == r.size() && RefCountingAlgebra::operator== (r);
}

/**
 *  setCharacteristic()
 */

void L::LookupFieldBB::setCharacteristic (unsigned c, unsigned deg)
{
   write();
   static_cast<unsigned*> (getDataPtr())[0] = c;
   static_cast<unsigned*> (getDataPtr())[1] = deg;
}


/******************  Lookup Field Base  **************************************/


namespace HIntLib
{
   namespace Private
   {
      extern
      const unsigned* lookupFields [HINTLIB_PRECALCULATED_FIELD_MAX_SIZE + 1];
   }
}

namespace
{
   template<typename T>
   inline unsigned* getPrecalculatedFieldData (unsigned size)
   {
      return 0;
   }

   template<>
   inline unsigned* getPrecalculatedFieldData<unsigned char> (unsigned size)
   {
      return size <= HINTLIB_PRECALCULATED_FIELD_MAX_SIZE
          ? const_cast<unsigned*> (HIntLib::Private::lookupFields[size]) : 0;
   }
}


/**
 *  privSetPointers()
 */

template<typename T>
void L::LookupFieldBase<T>::privSetPointers ()
{
   mulTable = reinterpret_cast<T*> (static_cast<unsigned*> (getDataPtr()) + 2);
   recTable     = mulTable   + HIntLib::sqr(size());
   sqrTable     = recTable   + size();
   orderTable   = sqrTable   + size();
   frobTable    = orderTable + size();
   invFrobTable = frobTable  + size();
}


/**
 *  Constructor
 */

template<typename T>
L::LookupFieldBase<T>::LookupFieldBase (unsigned _s, unsigned numData)
   : LookupFieldBB (_s, numData * sizeof (T) + 2 * sizeof (unsigned),
                    getPrecalculatedFieldData<T> (_s))
{
   if (     std::numeric_limits<T>::is_signed
       || ! std::numeric_limits<T>::is_integer)
   {
      throw InvalidType ("LookupField");
   }

   privSetPointers();
} 


/**
 *  dump()
 */

template<typename T>
void L::LookupFieldBase<T>::dump (std::ostream &o) const
{
   o << '\n' << setw (3) << "*";
   for (unsigned j = 0; j < size(); ++j)  o << setw(3) << j;
   o << "   shared: " << getRefCount() << '\n'
     << setw (3) << "^-1" << setw(3) << " ";
   for (unsigned j = 1; j < size(); ++j)  o << setw(3) << unsigned (recip (j));
   o << '\n';

   for (unsigned i = 0; i < size(); ++i)
   {
      o << setw(3) << i;

      for (unsigned j = 0; j < size(); ++j)
      {
         o << setw(3) << unsigned (mul (i,j));
      }
      o << '\n';
   }
   o << '\n';
}
   

/**
 *  power()
 */

template<typename T>
T  L::LookupFieldBase<T>::power (T x, unsigned k) const
{
   if (! x)  return T();

   const unsigned o = order (x);
   if (k >= o)  k %= o;

   T result = one();
   for(;;)
   {
      if (k & 1)  mulBy (result, x);
      if ((k >>= 1) == 0)  return result;
      square (x);
   }
}


/**
 *  setMul ()
 */

template<typename T>
void L::LookupFieldBase<T>::setMul (T a, T b, T prod)
{
   if (a    >= size())  throw LookupFieldSet (a,    size());
   if (b    >= size())  throw LookupFieldSet (b,    size());
   if (prod >= size())  throw LookupFieldSet (prod, size());

   write();
   mulTable [a + size() * b] = prod;
   mulTable [b + size() * a] = prod;

   if (a == b)  sqrTable[a] = prod;

   if (prod == 1)
   {
      recTable [a] = b;
      recTable [b] = a;
   }
}


/**
 *  setOrder()
 */

template<typename T>
void L::LookupFieldBase<T>::setOrder (T a, unsigned o)
{
   if (a >= size())  throw LookupFieldSet (a, size());
   if (o >= size() || o == 0)  throw LookupFieldSet (o, size());
   write();
   orderTable [a] = o;
}


/**
 *  setFrobenius()
 */

template<typename T>
void L::LookupFieldBase<T>::setFrobenius (T a, T b)
{
   if (a >= size())  throw LookupFieldSet (a, size());
   if (b >= size())  throw LookupFieldSet (b, size());
   write();
      frobTable [a] = b;
   invFrobTable [b] = a;
}


/******************  Lookup Field /Pow2 /Prime *******************************/


/**
 *  privSetPointers()
 */

template<typename T>
void L::LookupField<T>::privSetPointers ()
{
   addTable = invFrobTable + size();
   negTable = addTable + HIntLib::sqr(size());
   dblTable = negTable + size();
}

template<typename T>
void L::LookupField<T>::setPointers ()
{
   LookupFieldBase<T>::setPointers();
   privSetPointers();
}


/**
 *  Constructor
 */

template<typename T>
L::LookupField<T>::LookupField (unsigned _s)
   : LookupFieldBase<T> (_s, 2 *_s*_s + 7 *_s)
{
   privSetPointers();
}


/**
 *  times()
 */

template<typename T>
T  L::LookupField<T>::times (T x, unsigned k) const
{
   const unsigned c = characteristic();
   if (k >= c)  k %= c;

   T result = T();

   for(;;)
   {
      if (k & 1)  addTo (result, x);
      if ((k >>= 1) == 0)  return result;
      times2 (x);
   }
}


/**
 *  dump()
 */

template<typename T>
void L::LookupField<T>::dump (std::ostream &o) const
{
   LookupFieldBase<T>::dump (o);

   o << '\n' << setw (3) << "+";
   for (unsigned j = 0; j < size(); ++j)  o << setw(3) << j;
   o << '\n' << setw (3) << "-";
   for (unsigned j = 0; j < size(); ++j)  o << setw(3) << unsigned (neg (j));
   o << '\n';

   for (unsigned i = 0; i < size(); ++i)
   {
      o << setw(3) << i;

      for (unsigned j = 0; j < size(); ++j)
      {
         o << setw(3) << unsigned (add (i,j));
      }
      o << '\n';
   }
}
   

/**
 *  setAdd ()
 */

template<typename T>
void L::LookupField<T>::setAdd (T a, T b, T sum)
{
   if (a   >= size())  throw LookupFieldSet (a,   size());
   if (b   >= size())  throw LookupFieldSet (b,   size());
   if (sum >= size())  throw LookupFieldSet (sum, size());

   write();
   addTable [a + size() * b] = sum;
   addTable [b + size() * a] = sum;

   if (a == b)  dblTable [a] = sum;

   if (sum == 0)
   {
      negTable [a] = b;
      negTable [b] = a;
   }
}


/******************  Lookup Vector Space Base Base  **************************/


/**
 *  operator==
 */

bool L::LookupVectorSpaceBB::operator== (const LookupVectorSpaceBB &r) const
{
   return size() == r.size() && dimension() == r.dimension()
       && RefCountingAlgebra::operator== (r);
}


/**
 *  operator=
 */

L::LookupVectorSpaceBB &
L::LookupVectorSpaceBB::operator= (const LookupVectorSpaceBB &r)
{
   RefCountingAlgebra::operator=(r);
   s = r.s;
   dim = r.dim;
   return *this;
}


/******************  Lookup Vector Space Base  *******************************/

/**
 *  Constructor
 */

template<typename T, typename C>
L::LookupVectorSpaceBase<T,C>::LookupVectorSpaceBase
   (unsigned aSize, unsigned _dim, unsigned numData)
   : LookupVectorSpaceBB (powInt (aSize, _dim), _dim, numData * sizeof (T))
{
   if (     std::numeric_limits<T>::is_signed
       || ! std::numeric_limits<T>::is_integer
       ||   std::numeric_limits<T>::digits < std::numeric_limits<C>::digits)
   {
      throw InvalidType ("LookupVectorSpaceBB");
   }

   if (dimension() == 0)  throw FIXME(__FILE__, __LINE__);

   if (logInt (unsigned (T(-1)) + 1, aSize) < int (_dim))
   {
      throw FIXME(__FILE__, __LINE__);
   }
} 


/**
 *  dump()
 */

template<typename T, typename C>
void L::LookupVectorSpaceBase<T,C>::dump (std::ostream &o, unsigned as) const
{
   o << '\n' << setw (3) << "*";
   for (unsigned j = 0; j < as; ++j)  o << setw(3) << j;
   o << "   shared: " << getRefCount() << '\n';

   for (unsigned i = 0; i < size(); ++i)
   {
      o << setw(3) << i;

      for (unsigned j = 0; j < as; ++j)
      {
         o << setw(3) << unsigned (mul (i,j));
      }
      o << '\n';
   }
   o << '\n';
}
   

/******************  Lookup Vector Space /Pow2  ******************************/


/**
 *  setPointers()
 */

template<typename T, typename C>
void L::LookupVectorSpace<T,C>::setPointers ()
{
   mulTable = static_cast<T*> (getDataPtr());
   addTable = mulTable + size() * algebra.size();
   negTable = mulTable + size() * algebra.index(algebra.neg(algebra.one()));
}

template<typename T, typename C>
void L::LookupVectorSpacePow2<T,C>::setPointers ()
{
   mulTable = static_cast<T*> (getDataPtr());
}


/**
 *  Constructor
 */

template<typename T, typename C>
L::LookupVectorSpace<T,C>::LookupVectorSpace
   (const scalar_algebra& a, unsigned _dim)
   : LookupVectorSpaceBase<T,C> (
         a.size(),
         _dim,
         powInt (a.size(), _dim) * a.size()    // mul table
           + powInt (a.size(), 2 * _dim)),        // add table
     algebra (a),
     shift (algebra.size())
{
   setPointers();

   Array<C> as (dimension());
   Array<C> bs (dimension());
   Array<C> cs (dimension());

   for (unsigned i = 0; i < size(); ++i)
   {
      toCoord (element(i), &as[0]);

      // fill mulTable

      for (unsigned j = 0; j < algebra.size(); ++j)
      {
         C s = algebra.element (j);
        
         for (unsigned k = 0; k < dimension(); ++k)
         {
            cs[k] = algebra.mul (as[k], s);
         }

         fromCoord (mulTable [i + size() * j], &cs[0]);
      }

      // Fill addTable

      for (unsigned j = 0; j < size(); ++j)
      {
         toCoord (element(j), &bs[0]);

         for (unsigned k = 0; k < dimension(); ++k)
         {
            cs[k] = algebra.add (as[k], bs[k]);
         }

         fromCoord (addTable [i + size() * j], &cs[0]);
      }
   }
}

template<typename T, typename C>
L::LookupVectorSpacePow2<T,C>::LookupVectorSpacePow2
   (const scalar_algebra& a, unsigned _dim)
   : LookupVectorSpaceBase<T,C> (
         a.size(),
         _dim,
         powInt (a.size(), _dim) * a.size()),  // mul table
     algebra (a),
     shift (ms1 (algebra.size())),
     mask (algebra.size() - 1)
{
   setPointers();

   Array<C> as (dimension());
   Array<C> cs (dimension());

   for (unsigned i = 0; i < size(); ++i)
   {
      toCoord (element(i), &as[0]);

      for (unsigned j = 0; j < algebra.size(); ++j)
      {
         C s = element (j);
        
         for (unsigned k = 0; k < dimension(); ++k)
         {
            cs[k] = algebra.mul (as[k], s);
         }

         fromCoord (mulTable [i + size() * j], &cs[0]);
      }
   }
}


/**
 *  Copy Constructor
 */

template<typename T, typename C>
L::LookupVectorSpace<T,C>::LookupVectorSpace (const LookupVectorSpace<T,C> &r)
   : LookupVectorSpaceBase<T,C> (r),
     algebra (r.algebra),
     shift (r.shift)
{
   setPointers();
}

template<typename T, typename C>
L::LookupVectorSpacePow2<T,C>::LookupVectorSpacePow2
   (const LookupVectorSpacePow2<T,C> &r)
   : LookupVectorSpaceBase<T,C>(r),
     BitOpBasedAddition<T> (),
     algebra (r.algebra),
     shift (r.shift),
     mask (r.mask)
{
   setPointers();
}


/**
 *  operator=
 */

template<typename T, typename C>
L::LookupVectorSpace<T,C> &
L::LookupVectorSpace<T,C>::operator= (const LookupVectorSpace<T,C> &r)
{
   LookupVectorSpaceBB::operator=(r);
   algebra  = r.algebra;
   shift    = r.shift;
   return *this;
}

template<typename T, typename C>
L::LookupVectorSpacePow2<T,C> &
L::LookupVectorSpacePow2<T,C>::operator= (const LookupVectorSpacePow2<T,C> &r)
{
   LookupVectorSpaceBB::operator=(r);
   algebra = r.algebra;
   shift   = r.shift;
   mask    = r.mask;
   return *this;
}


/**
 *  times()
 */

template<typename T, typename C>
T  L::LookupVectorSpace<T,C>::times (T x, unsigned k) const
{
   T result = T();
   if (k >= algebra.characteristic())  k %= algebra.characteristic();

   for(;;)
   {
      if (k & 1)  addTo (result, x);
      if ((k >>= 1) == 0)  return result;
      times2 (x);
   }
}


/**
 *  dump()
 */

template<typename T, typename C>
void L::LookupVectorSpace<T,C>::dump (std::ostream &o) const
{
   LookupVectorSpaceBase<T,C>::dump (o, algebra.size());

   o << '\n' << setw (3) << "+";
   for (unsigned j = 0; j < size(); ++j)  o << setw(3) << j;
   o << '\n' << setw (3) << "-";
   for (unsigned j = 0; j < size(); ++j)  o << setw(3) << unsigned (neg (j));
   o << '\n';

   for (unsigned i = 0; i < size(); ++i)
   {
      o << setw(3) << i;

      for (unsigned j = 0; j < size(); ++j)
      {
         o << setw(3) << unsigned (add (i,j));
      }
      o << '\n';
   }
}


/**
 *  print Short ()
 */

template<typename T, typename C>
void
L::LookupVectorSpace<T,C>::printShort (std::ostream &o, T x) const
{
   Private::Printer ss (o);

   ss << '(';
   for (unsigned i = 0; i < dimension(); ++i)
   {
      if (i > 0)  ss << ',';
      algebra.printShort (ss, coord (x, i));
   }
   ss << ')';
}

template<typename T, typename C>
void
L::LookupVectorSpacePow2<T,C>::printShort (std::ostream &o, T x) const
{
   Private::Printer ss (o);

   ss << '(';
   for (unsigned i = 0; i < dimension(); ++i)
   {
      if (i > 0)  ss << ',';
      algebra.printShort (ss, coord (x, i));
   }
   ss << ')';
}


/**
 *  print ()
 */

template<typename T, typename C>
void
L::LookupVectorSpace<T,C>::print (std::ostream &o, T x) const
{
   Private::Printer ss (o);

   printShort (ss, x);
   ss << ' ';
   printSuffix (ss);
}

template<typename T, typename C>
void
L::LookupVectorSpacePow2<T,C>::print (std::ostream &o, T x) const
{
   Private::Printer ss (o);

   printShort (ss, x);
   ss << ' ';
   printSuffix (ss);
}


/******************  Vector Space Pow 2  *************************************/


/**
 *  Constructor
 */

template<typename T>
L::VectorSpacePow2<T>::VectorSpacePow2
   (const scalar_algebra& a, unsigned _dim)
   : algebra (a),
     vecalg (a, logInt (unsigned (static_cast<unsigned char>(-1))+1, a.size())),
     baseBits (ms1 (algebra.size())),
     vecBits  (baseBits * vecalg.dimension()),
     baseMask (~T(0) >> (std::numeric_limits<T>::digits - baseBits)),
     vecMask  (~T(0) >> (std::numeric_limits<T>::digits - vecBits )),
     dim (_dim),
     vecDim ((dim - 1) / (vecalg.dimension()) + 1)
{
   if (baseBits * dim > unsigned (std::numeric_limits<T>::digits))
   {
      throw FIXME (__FILE__, __LINE__);
   }
}


/**
 *  size ()
 */

template<typename T>
unsigned
L::VectorSpacePow2<T>::size () const
{
   return unsigned (std::numeric_limits<unsigned>::digits) > dim * baseBits
        ? 1u << (dim * baseBits)
#if 0
        : std::numeric_limits<unsigned>::max();
#else
        : 0;
#endif
}


/**
 *  mul ()
 */

template<typename T>
T
L::VectorSpacePow2<T>::mul (const T& a, scalar_type l) const
{
   T x (0);

   for (unsigned i = vecDim * vecBits; i > 0; )
   {
      x = (x << vecBits) | vecalg.mul ((a >> (i -= vecBits)) & vecMask, l);
   }

   return x;
}


/**
 *  print Short ()
 */

template<typename T>
void
L::VectorSpacePow2<T>::printShort (std::ostream &o, const T& x) const
{
   Private::Printer ss (o);

   ss << '(';
   for (unsigned i = 0; i < dimension(); ++i)
   {
      if (i > 0)  ss << ',';
      algebra.printShort (ss, coord (x, i));
   }
   ss << ')';
}


/**
 *  print ()
 */

template<typename T>
void
L::VectorSpacePow2<T>::print (std::ostream &o, const T& x) const
{
   Private::Printer ss (o);

   printShort (ss, x);
   ss << ' ';
   printSuffix (ss);
}


namespace HIntLib
{
#define HINTLIB_INSTANTIATE(X) \
   template LookupFieldBase<X>::LookupFieldBase (unsigned, unsigned); \
   template void LookupFieldBase<X>::dump (std::ostream &) const; \
   template X    LookupFieldBase<X>::power (X, unsigned) const; \
   template void LookupFieldBase<X>::setMul (X, X, X); \
   template void LookupFieldBase<X>::setOrder (X, unsigned); \
   template void LookupFieldBase<X>::setFrobenius (X, X); \
   template void LookupFieldBase<X>::privSetPointers (); \
   template void LookupField<X>::setPointers (); \
   template void LookupField<X>::privSetPointers (); \
   template LookupField<X>::LookupField (unsigned); \
   template X    LookupField<X>::times (X, unsigned) const; \
   template void LookupField<X>::dump (std::ostream &) const; \
   template void LookupField<X>::setAdd (X, X, X);

   HINTLIB_INSTANTIATE (unsigned char)
#undef HINTLIB_INSTANTIATE

#define HINTLIB_INSTANTIATE(X,Y) \
   template LookupVectorSpaceBase<X,Y>::LookupVectorSpaceBase \
               (unsigned, unsigned, unsigned); \
   template void LookupVectorSpaceBase<X,Y>::dump \
               (std::ostream &,unsigned) const; \
   template void LookupVectorSpace<X,Y>::dump (std::ostream &) const; \
   template void LookupVectorSpace<X,Y>::setPointers (); \
   template void LookupVectorSpacePow2<X,Y>::setPointers (); \
   template LookupVectorSpace<X,Y>::LookupVectorSpace \
               (const scalar_algebra&, unsigned); \
   template LookupVectorSpacePow2<X,Y>::LookupVectorSpacePow2 \
               (const scalar_algebra&, unsigned); \
   template LookupVectorSpace<X,Y>::LookupVectorSpace \
               (const LookupVectorSpace<X,Y> &); \
   template LookupVectorSpacePow2<X,Y>::LookupVectorSpacePow2 \
               (const LookupVectorSpacePow2<X,Y> &); \
   template LookupVectorSpace<X,Y> & \
            LookupVectorSpace<X,Y>::operator= \
               (const LookupVectorSpace<X,Y>&); \
   template LookupVectorSpacePow2<X,Y> & \
            LookupVectorSpacePow2<X,Y>::operator= \
               (const LookupVectorSpacePow2<X,Y>&);\
   template X LookupVectorSpace<X,Y>::times (X, unsigned) const; \
   template void LookupVectorSpace<X,Y>::print (std::ostream&,X) const; \
   template void LookupVectorSpacePow2<X,Y>::print (std::ostream&,X) const; \
   template void LookupVectorSpace<X,Y>::printShort (std::ostream&,X) const; \
   template void LookupVectorSpacePow2<X,Y>::printShort (std::ostream&,X) const;

   HINTLIB_INSTANTIATE (unsigned char, unsigned char)
#undef HINTLIB_INSTANTIATE

#define HINTLIB_INSTANTIATE(X) \
   template VectorSpacePow2<X>::VectorSpacePow2 \
                (const scalar_algebra&, unsigned); \
   template unsigned VectorSpacePow2<X>::size() const; \
   template X VectorSpacePow2<X>::mul (const X& a, scalar_type l) const; \
   template void VectorSpacePow2<X>::print (std::ostream&,const X&) const; \
   template void VectorSpacePow2<X>::printShort (std::ostream&, const X&) const;

   HINTLIB_INSTANTIATE (u32)
#ifdef HINTLIB_U32_NOT_EQUAL_U64
   HINTLIB_INSTANTIATE (u64)
#endif

#undef HINTLIB_INSTANTIATE
}

