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

#include <algorithm>
#include <iomanip>

#include <HIntLib/lookupfield.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/hlmath.h>
#include <HIntLib/array.h>
#include <HIntLib/output.h>
#include <HIntLib/exception.h>
#include <HIntLib/bitop.h>
#include <HIntLib/prime.h>

namespace L = HIntLib;
namespace P = L::Private;

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
P::printVectorSpaceName (std::ostream &o, unsigned b, unsigned d)
{
   Printer ss (o);
   ss << "(GF" << b << ')';;
   ss.power (d);
}
#ifdef HINTLIB_BUILD_WCHAR
void
P::printVectorSpaceName (std::wostream &o, unsigned b, unsigned d)
{
   WPrinter ss (o);
   ss << L"(GF" << b << L')';
   ss.power (d);
}
#endif

void
P::printFieldName (std::ostream &o, unsigned b)
{
   Printer ss (o);
   ss << "GF" << b;
}
#ifdef HINTLIB_BUILD_WCHAR
void
P::printFieldName (std::wostream &o, unsigned b)
{
   WPrinter ss (o);
   ss << L"GF" << b;
}
#endif

void
P::printSuff (std::ostream &o, unsigned b)
{
   Printer ss (o);
   ss << "(GF" << b << ')';
}
#ifdef HINTLIB_BUILD_WCHAR
void
P::printSuff (std::wostream &o, unsigned b)
{
   WPrinter ss (o);
   ss << L"(GF" << b << L')';
}
#endif

void
P::printNumberSuffix (std::ostream &o, int x, unsigned b)
{
   Printer ss (o);
   ss << x << " (GF" << b << ')';
}
#ifdef HINTLIB_BUILD_WCHAR
void
P::printNumberSuffix (std::wostream &o, int x, unsigned b)
{
   WPrinter ss (o);
   ss << x << L" (GF" << b << L')';
}
#endif

void
P::printNumber (std::ostream &o, int x)
{
   o << x;
}
#ifdef HINTLIB_BUILD_WCHAR
void
P::printNumber (std::wostream &o, int x)
{
   o << x;
}
#endif


/******************  Lookup Field Base Base  *********************************/


/**
 *  operator=
 */

P::LookupFieldBB &
P::LookupFieldBB::operator= (const LookupFieldBB &r)
{
   RefCountingAlgebra::operator=(r);
   s = r.s;
   return *this;
}

/**
 *  operator==
 */

bool
P::LookupFieldBB::operator== (const LookupFieldBB &r) const
{
   return size() == r.size() && RefCountingAlgebra::operator== (r);
}

/**
 *  setCharacteristic()
 */

void
P::LookupFieldBB::setCharacteristic (unsigned c, unsigned deg)
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
void
P::LookupFieldBase<T>::privSetPointers ()
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
P::LookupFieldBase<T>::LookupFieldBase (unsigned _s, unsigned numData)
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
void
P::LookupFieldBase<T>::dump (std::ostream &o) const
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
#ifdef HINTLIB_BUILD_WCHAR
template<typename T>
void
P::LookupFieldBase<T>::dump (std::wostream &o) const
{
   o << L'\n' << setw (3) << L"*";
   for (unsigned j = 0; j < size(); ++j)  o << setw(3) << j;
   o << L"   shared: " << getRefCount() << L'\n'
     << setw (3) << L"^-1" << setw(3) << L" ";
   for (unsigned j = 1; j < size(); ++j)  o << setw(3) << unsigned (recip (j));
   o << L'\n';

   for (unsigned i = 0; i < size(); ++i)
   {
      o << setw(3) << i;

      for (unsigned j = 0; j < size(); ++j)
      {
         o << setw(3) << unsigned (mul (i,j));
      }
      o << L'\n';
   }
   o << L'\n';
}
#endif


/**
 *  power()
 */

template<typename T>
T
P::LookupFieldBase<T>::power (T x, unsigned k) const
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
void
P::LookupFieldBase<T>::setMul (T a, T b, T prod)
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
void
P::LookupFieldBase<T>::setOrder (T a, unsigned o)
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
void
P::LookupFieldBase<T>::setFrobenius (T a, T b)
{
   if (a >= size())  throw LookupFieldSet (a, size());
   if (b >= size())  throw LookupFieldSet (b, size());
   write();
      frobTable [a] = b;
   invFrobTable [b] = a;
}


/**
 *  frobenius()
 */

template<typename T>
T
P::LookupFieldBase<T>::frobenius (T a, unsigned k) const
{
   while (k-- > 0)  a = frobenius(a);
   return a;
}


/******************  Lookup Field /Pow2 /Prime *******************************/


/**
 *  privSetPointers()
 */

template<typename T>
void
L::LookupField<T>::privSetPointers ()
{
   const unsigned S = this->size();

   addTable = this->invFrobTable + S;
   negTable = addTable + HIntLib::sqr(S);
   dblTable = negTable + S;
}

template<typename T>
void
L::LookupField<T>::setPointers ()
{
   Private::LookupFieldBase<T>::setPointers();
   privSetPointers();
}


/**
 *  Constructor
 */

template<typename T>
L::LookupField<T>::LookupField (unsigned _s)
   : Private::LookupFieldBase<T> (_s, 2 *_s*_s + 7 *_s)
{
   privSetPointers();
}


/**
 *  times()
 */

template<typename T>
T
L::LookupField<T>::times (T x, unsigned k) const
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
void
L::LookupField<T>::dump (std::ostream &o) const
{
   Private::LookupFieldBase<T>::dump (o);

   const unsigned S = this->size();

   o << '\n' << setw (3) << "+";
   for (unsigned j = 0; j != S; ++j)  o << setw(3) << j;
   o << '\n' << setw (3) << "-";
   for (unsigned j = 0; j != S; ++j)  o << setw(3) << unsigned (neg (j));
   o << '\n';

   for (unsigned i = 0; i != S; ++i)
   {
      o << setw(3) << i;
      for (unsigned j = 0; j != S; ++j)  o << setw(3) << unsigned (add (i,j));
      o << '\n';
   }
}
#ifdef HINTLIB_BUILD_WCHAR
template<typename T>
void
L::LookupField<T>::dump (std::wostream &o) const
{
   Private::LookupFieldBase<T>::dump (o);

   const unsigned S = this->size();

   o << L'\n' << setw (3) << L"+";
   for (unsigned j = 0; j != S; ++j)  o << setw(3) << j;
   o << L'\n' << setw (3);
#if HINTLIB_CHARACTER_SET >= 3
   o << L"\x2212";   // MINUS SIGN
#else
   o << L"-";
#endif
   for (unsigned j = 0; j != S; ++j)  o << setw(3) << unsigned (neg (j));
   o << L'\n';

   for (unsigned i = 0; i != S; ++i)
   {
      o << setw(3) << i;
      for (unsigned j = 0; j != S; ++j)  o << setw(3) << unsigned (add (i,j));
      o << L'\n';
   }
}
#endif


/**
 *  setAdd ()
 */

template<typename T>
void
L::LookupField<T>::setAdd (T a, T b, T sum)
{
   const unsigned S = this->size();

   if (a   >= S)  throw LookupFieldSet (a,   S);
   if (b   >= S)  throw LookupFieldSet (b,   S);
   if (sum >= S)  throw LookupFieldSet (sum, S);

   this->write();
   addTable [a + S * b] = sum;
   addTable [b + S * a] = sum;

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

bool
P::LookupVectorSpaceBB::operator== (const LookupVectorSpaceBB &r) const
{
   return size() == r.size() && dimension() == r.dimension()
       && RefCountingAlgebra::operator== (r);
}


/**
 *  operator=
 */

P::LookupVectorSpaceBB &
P::LookupVectorSpaceBB::operator= (const LookupVectorSpaceBB &r)
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
P::LookupVectorSpaceBase<T,C>::LookupVectorSpaceBase
   (unsigned aSize, int _dim, unsigned numData)
   : LookupVectorSpaceBB (powInt (aSize, _dim), _dim, numData * sizeof (T))
{
   if (     std::numeric_limits<T>::is_signed
       || ! std::numeric_limits<T>::is_integer
       ||   std::numeric_limits<T>::digits < std::numeric_limits<C>::digits)
   {
      throw InvalidType ("LookupVectorSpaceBB");
   }

   if (dimension() <= 0)  throw FIXME(__FILE__, __LINE__);

   if (logInt (unsigned(T(-1)) + 1, aSize) < _dim)
   {
      throw FIXME(__FILE__, __LINE__);
   }
}


/**
 *  dump()
 */

template<typename T, typename C>
void
P::LookupVectorSpaceBase<T,C>::dump (std::ostream &o, unsigned as) const
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
#ifdef HINTLIB_BUILD_WCHAR
template<typename T, typename C>
void
P::LookupVectorSpaceBase<T,C>::dump (std::wostream &o, unsigned as) const
{
   o << L'\n' << setw (3) << "*";
   for (unsigned j = 0; j < as; ++j)  o << setw(3) << j;
   o << L"   shared: " << getRefCount() << L'\n';

   for (unsigned i = 0; i < size(); ++i)
   {
      o << setw(3) << i;

      for (unsigned j = 0; j < as; ++j)
      {
         o << setw(3) << unsigned (mul (i,j));
      }
      o << L'\n';
   }
   o << L'\n';
}
#endif


/******************  Lookup Vector Space /Pow2  ******************************/


/**
 *  setPointers()
 */

template<typename T, typename C>
void L::LookupVectorSpace<T,C>::setPointers ()
{
   const unsigned S = this->size();

   this->mulTable = static_cast<T*> (this->getDataPtr());
   addTable = this->mulTable + S * algebra.size();
   negTable = this->mulTable + S * algebra.index(algebra.neg(algebra.one()));
}

template<typename T, typename C>
void L::LookupVectorSpacePow2<T,C>::setPointers ()
{
   this->mulTable = static_cast<T*> (this->getDataPtr());
}


/**
 *  Constructor
 */

template<typename T, typename C>
L::LookupVectorSpace<T,C>::LookupVectorSpace
   (const scalar_algebra& a, unsigned _dim)
   : P::LookupVectorSpaceBase<T,C> (
         a.size(),
         _dim,
         powInt (a.size(), _dim) * a.size()    // mul table
           + powInt (a.size(), 2 * _dim)),     // add table
     algebra (a),
     shift (algebra.size())
{
   setPointers();

   const unsigned DIM = this->dimension();
   const unsigned S = this->size();

   Array<C> as (DIM);
   Array<C> bs (DIM);
   Array<C> cs (DIM);

   for (unsigned i = 0; i != S; ++i)
   {
      toCoord (this->element(i), &as[0]);

      // fill mulTable

      for (unsigned j = 0; j != algebra.size(); ++j)
      {
         C s = algebra.element (j);

         for (unsigned k = 0; k != DIM; ++k)
         {
            cs[k] = algebra.mul (as[k], s);
         }

         fromCoord (this->mulTable [i + S * j], &cs[0]);
      }

      // Fill addTable

      for (unsigned j = 0; j != S; ++j)
      {
         toCoord (this->element(j), &bs[0]);

         for (unsigned k = 0; k != DIM; ++k)
         {
            cs[k] = algebra.add (as[k], bs[k]);
         }

         fromCoord (addTable [i + S * j], &cs[0]);
      }
   }
}

template<typename T, typename C>
L::LookupVectorSpacePow2<T,C>::LookupVectorSpacePow2
   (const scalar_algebra& a, unsigned _dim)
   : P::LookupVectorSpaceBase<T,C> (
         a.size(),
         _dim,
         powInt (a.size(), _dim) * a.size()),  // mul table
     algebra (a),
     shift (ms1 (algebra.size())),
     mask (algebra.size() - 1)
{
   setPointers();

   const unsigned DIM = this->dimension();
   const unsigned S = this->size();

   Array<C> as (DIM);
   Array<C> cs (DIM);

   for (unsigned i = 0; i != S; ++i)
   {
      toCoord (this->element(i), &as[0]);

      for (unsigned j = 0; j < algebra.size(); ++j)
      {
         C s = this->element (j);

         for (unsigned k = 0; k != DIM; ++k)
         {
            cs[k] = algebra.mul (as[k], s);
         }

         fromCoord (this->mulTable [i + S * j], &cs[0]);
      }
   }
}


/**
 *  Copy Constructor
 */

template<typename T, typename C>
L::LookupVectorSpace<T,C>::LookupVectorSpace (const LookupVectorSpace<T,C> &r)
   : P::LookupVectorSpaceBase<T,C> (r),
     algebra (r.algebra),
     shift (r.shift)
{
   setPointers();
}

template<typename T, typename C>
L::LookupVectorSpacePow2<T,C>::LookupVectorSpacePow2
   (const LookupVectorSpacePow2<T,C> &r)
   : P::LookupVectorSpaceBase<T,C>(r),
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
   P::LookupVectorSpaceBB::operator=(r);
   algebra  = r.algebra;
   shift    = r.shift;
   return *this;
}

template<typename T, typename C>
L::LookupVectorSpacePow2<T,C> &
L::LookupVectorSpacePow2<T,C>::operator= (const LookupVectorSpacePow2<T,C> &r)
{
   P::LookupVectorSpaceBB::operator=(r);
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
   P::LookupVectorSpaceBase<T,C>::dump (o, algebra.size());

   const unsigned S = this->size();

   o << '\n' << setw (3) << "+";
   for (unsigned j = 0; j != S; ++j)  o << setw(3) << j;
   o << '\n' << setw (3) << "-";
   for (unsigned j = 0; j != S; ++j)  o << setw(3) << unsigned (neg (j));
   o << '\n';

   for (unsigned i = 0; i != S; ++i)
   {
      o << setw(3) << i;

      for (unsigned j = 0; j != S; ++j)
      {
         o << setw(3) << unsigned (add (i,j));
      }
      o << '\n';
   }
}
#ifdef HINTLIB_BUILD_WCHAR
template<typename T, typename C>
void L::LookupVectorSpace<T,C>::dump (std::wostream &o) const
{
   P::LookupVectorSpaceBase<T,C>::dump (o, algebra.size());

   const unsigned S = this->size();

   o << L'\n' << setw (3) << L"+";
   for (unsigned j = 0; j != S; ++j)  o << setw(3) << j;
   o << L'\n' << setw (3);
#if HINTLIB_CHARACTER_SET >= 3
   o << L"\x2212";   // MINUS SIGN
#else
   o << L"-";
#endif
   for (unsigned j = 0; j != S; ++j)  o << setw(3) << unsigned (neg (j));
   o << L'\n';

   for (unsigned i = 0; i != S; ++i)
   {
      o << setw(3) << i;

      for (unsigned j = 0; j != S; ++j)
      {
         o << setw(3) << unsigned (add (i,j));
      }
      o << L'\n';
   }
}
#endif


/**
 *  printVector()
 */

#ifdef HINTLIB_OSTREAM_IS_BASIC_OSTREAM
template<typename A, typename C>
void
P::vectorPrintShort (
      const A& a, const typename A::type& x, std::basic_ostream<C>& o)
{
#else
template<typename A>
void
P::vectorPrintShort (
      const A& a, const typename A::type& x, std::ostream& o)
{
   typedef char C;
#endif
   typename PrinterSelector<C>::printer ss (o);

   const int DIM = a.dimension();

   ss << C('(');
   for (int i = 0; i != DIM; ++i)
   {
      if (i > 0)  ss << C(',');
      a.getScalarAlgebra().printShort (ss, a.coord (x, i));
   }
   ss << C(')');

}

#ifdef HINTLIB_OSTREAM_IS_BASIC_OSTREAM
template<typename A, typename C>
void
P::vectorPrint (
      const A& a, const typename A::type& x, std::basic_ostream<C>& o)
{
#else
template<typename A>
void
P::vectorPrint (
      const A& a, const typename A::type& x, std::ostream& o)
{
   typedef char C;
#endif
   typename PrinterSelector<C>::printer ss (o);

   vectorPrintShort (a, x, ss);
   ss << C(' ');
   a.printSuffix (ss);
}

#ifdef HINTLIB_OSTREAM_IS_BASIC_OSTREAM

#define HINTLIB_INSTANTIATE_VECTOR_PRINT_1(X,C) \
   namespace Private { \
   template void vectorPrintShort \
      (const X&, const X::type& x, std::basic_ostream<C>&); \
   template void vectorPrint \
      (const X&, const X::type& x, std::basic_ostream<C>&); \
   }

#ifdef HINTLIB_BUILD_WCHAR
#define HINTLIB_INSTANTIATE_VECTOR_PRINT(X) \
   HINTLIB_INSTANTIATE_VECTOR_PRINT_1(X,char) \
   HINTLIB_INSTANTIATE_VECTOR_PRINT_1(X,wchar_t)
#else
#define HINTLIB_INSTANTIATE_VECTOR_PRINT(X) \
   HINTLIB_INSTANTIATE_VECTOR_PRINT_1(X,char)
#endif

#else

#define HINTLIB_INSTANTIATE_VECTOR_PRINT(X) \
   namespace Private { \
   template void vectorPrintShort (const X&, const X::type& x, std::ostream&); \
   template void vectorPrint      (const X&, const X::type& x, std::ostream&); \
   }

#endif


/******************  Vector Space Pow 2  *************************************/


/**
 *  Constructor
 */

template<typename T>
L::VectorSpacePow2<T>::VectorSpacePow2
   (const scalar_algebra& a, int _dim)
   : algebra (a),
     vecalg (a, logInt (unsigned (static_cast<unsigned char>(-1))+1, a.size())),
     baseBits (ms1 (algebra.size())),
     vecBits  (baseBits * vecalg.dimension()),
     baseMask (~T(0) >> (std::numeric_limits<T>::digits - baseBits)),
     vecMask  (~T(0) >> (std::numeric_limits<T>::digits - vecBits )),
     dim (_dim),
     vecDim ((dim - 1) / (vecalg.dimension()) + 1)
{
   if (baseBits * dim > std::numeric_limits<T>::digits)
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
   return std::numeric_limits<unsigned>::digits > dim * baseBits
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
      x = (x << vecBits)
        | vecalg.mul (scalar_type((a >> (i -= vecBits)) & vecMask), l);
   }

   return x;
}


namespace HIntLib
{
   // LookupField, LookupFieldPow2, LookupFieldPrime, and base classes

#define HINTLIB_INSTANTIATE(X) \
   namespace Private { \
   template LookupFieldBase<X >::LookupFieldBase (unsigned, unsigned); \
   template void LookupFieldBase<X >::dump (std::ostream &) const; \
   template X    LookupFieldBase<X >::power (X, unsigned) const; \
   template X    LookupFieldBase<X >::frobenius (X, unsigned) const; \
   template void LookupFieldBase<X >::setMul (X, X, X); \
   template void LookupFieldBase<X >::setOrder (X, unsigned); \
   template void LookupFieldBase<X >::setFrobenius (X, X); \
   template void LookupFieldBase<X >::privSetPointers (); \
   } \
   template void LookupField<X >::setPointers (); \
   template void LookupField<X >::privSetPointers (); \
   template LookupField<X >::LookupField (unsigned); \
   template X    LookupField<X >::times (X, unsigned) const; \
   template void LookupField<X >::dump (std::ostream &) const; \
   template void LookupField<X >::setAdd (X, X, X);

   HINTLIB_INSTANTIATE (unsigned char)
#undef HINTLIB_INSTANTIATE

   // LookupVectorSpace, LookupVectorSpacePow2, and base classes

#define HINTLIB_INSTANTIATE(X,Y,NAME) \
   namespace Private { \
   template LookupVectorSpaceBase<X,Y >::LookupVectorSpaceBase \
               (unsigned, int, unsigned); \
   template void LookupVectorSpaceBase<X,Y >::dump \
               (std::ostream &,unsigned) const; \
   } \
   template void LookupVectorSpace<X,Y >::dump (std::ostream &) const; \
   template void LookupVectorSpace<X,Y >::setPointers (); \
   template void LookupVectorSpacePow2<X,Y >::setPointers (); \
   template LookupVectorSpace<X,Y >::LookupVectorSpace \
               (const scalar_algebra&, unsigned); \
   template LookupVectorSpacePow2<X,Y >::LookupVectorSpacePow2 \
               (const scalar_algebra&, unsigned); \
   template LookupVectorSpace<X,Y >::LookupVectorSpace \
               (const LookupVectorSpace<X,Y > &); \
   template LookupVectorSpacePow2<X,Y >::LookupVectorSpacePow2 \
               (const LookupVectorSpacePow2<X,Y > &); \
   template LookupVectorSpace<X,Y > & \
            LookupVectorSpace<X,Y >::operator= \
               (const LookupVectorSpace<X,Y >&); \
   template LookupVectorSpacePow2<X,Y > & \
            LookupVectorSpacePow2<X,Y >::operator= \
               (const LookupVectorSpacePow2<X,Y >&);\
   template X LookupVectorSpace<X,Y >::times (X, unsigned) const; \
   typedef LookupVectorSpace    <unsigned char,unsigned char > NAME ## type1; \
   typedef LookupVectorSpacePow2<unsigned char,unsigned char > NAME ## type2; \
   HINTLIB_INSTANTIATE_VECTOR_PRINT(NAME ## type1) \
   HINTLIB_INSTANTIATE_VECTOR_PRINT(NAME ## type2)

   HINTLIB_INSTANTIATE (unsigned char,unsigned char,char)
#undef HINTLIB_INSTANTIATE

   // VectorSpacePow2

#define HINTLIB_INSTANTIATE(X) \
   template VectorSpacePow2<X >::VectorSpacePow2 \
                (const scalar_algebra&, int); \
   template unsigned VectorSpacePow2<X >::size() const; \
   template X VectorSpacePow2<X >::mul (const X& a, scalar_type l) const; \
   HINTLIB_INSTANTIATE_VECTOR_PRINT(VectorSpacePow2<X >)

   HINTLIB_INSTANTIATE (u32)
#ifdef HINTLIB_U32_NOT_EQUAL_U64
   HINTLIB_INSTANTIATE (u64)
#endif

#undef HINTLIB_INSTANTIATE
}

