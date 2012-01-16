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

#ifndef HINTLIB_LOOKUP_FIELD_H
#define HINTLIB_LOOKUP_FIELD_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <iosfwd>

#include <HIntLib/mymath.h>
#include <HIntLib/algebra.h>


namespace HIntLib
{

/**
 *  Class Diagram
 *
 *  RefCountingAlgebra
 *  |
 *  \--LookupFieldBB
 *     |
 *     \--LookupFieldBase<T>
 *        |
 *        |--LookupField<T>
 *        |
 *        \--LookupFieldMulOnly<T>
 *           |
 *           |--LookupFieldPow2
 *           \--LookupFieldPrime
 */

/**
 *  Ref Counting Algebra
 */

class RefCountingAlgebra
{
private:
   unsigned* refCount;
   size_t size;

   void copy()  { if (*refCount)  ++*refCount; }
   void destroy();
   char* charPtr() const  { return reinterpret_cast<char*> (refCount + 1); }

protected:
   explicit RefCountingAlgebra (size_t memory);
            RefCountingAlgebra (size_t memory, unsigned* data);
   RefCountingAlgebra (const RefCountingAlgebra &f)
      : refCount (f.refCount), size (f.size)  { copy(); }

   RefCountingAlgebra & operator= (const RefCountingAlgebra &);

   bool operator== (const RefCountingAlgebra &)  const;

   void write();
   virtual void setPointers() = 0;

   unsigned getRefCount() const  { return *refCount; }
   void* getDataPtr() const  { return reinterpret_cast<void*> (refCount + 1); }

public:
   virtual ~RefCountingAlgebra()  { destroy(); }

   void clear ();
};


/*****************************************************************************/

/**
 *  Printing in namespace Priv
 */

namespace Priv
{
   void printVectorSpaceName (std::ostream &, unsigned, unsigned);
   void printFieldName (std::ostream &, unsigned);
   void printNumber (std::ostream&, unsigned);
   void printNumberSuffix (std::ostream&, unsigned, unsigned);
   void printSuffix (std::ostream&, unsigned);
}

/**
 *  Lookup Field Base Base
 */

class LookupFieldBB : public RefCountingAlgebra
{
public:
   unsigned size() const  { return s; }
   void setCharacteristic(unsigned, unsigned);

   void printSuffix (std::ostream &o) const  { Priv::printSuffix (o, s); }

protected:
   LookupFieldBB (unsigned _s, size_t memory)
      : RefCountingAlgebra (memory), s (_s)  { setCharacteristic (s, 1); }
   LookupFieldBB (unsigned _s, size_t memory, unsigned* data)
      : RefCountingAlgebra (memory, data), s (_s)  {}
   LookupFieldBB (const LookupFieldBB &f)
      : RefCountingAlgebra (f), s (f.s)  {}

   LookupFieldBB& operator=(const LookupFieldBB&);
   
   bool operator== (const LookupFieldBB &)  const;

   unsigned getCharacteristic() const
      { return static_cast<unsigned*> (getDataPtr())[0]; }
   unsigned getExtensionDegree() const
      { return static_cast<unsigned*> (getDataPtr())[1]; }

private:
   unsigned s;
};

inline std::ostream& operator<< (std::ostream &o, const LookupFieldBB &f)
   { Priv::printFieldName (o, f.size()); return o; }


/**
 *  Lookup Field Base
 */

template <typename T>
class LookupFieldBase : public LookupFieldBB
{
private:
   void privSetPointers();

protected:
   T* mulTable;
   T* recTable;
   T* orderTable;

   LookupFieldBase (unsigned _s, unsigned numData);
   LookupFieldBase (const LookupFieldBase<T> &r) : LookupFieldBB (r)
      { privSetPointers(); }
   LookupFieldBase & operator= (const LookupFieldBase<T> &r)
      { LookupFieldBB::operator=(r); privSetPointers(); return *this; }

   void setPointers()  { privSetPointers(); }

public:
   typedef nopolynomial_tag polynomial_category;
   typedef T type;

   T one()  const  { return 1; }

   bool is1 (T a) const  { return a == 1; }

   T element(unsigned i) const  { return T(i); }
   unsigned index (T x) const   { return unsigned (x); }

   T  mul   (T  a, T b) const  { return mulTable [a + size() * b]; }
   T& mulBy (T& a, T b) const  { return a = mul (a,b); }

   T  recip (T  a) const  { return recTable [a]; }

   T  div   (T  a, T b) const  { return mul   (a, recip (b)); }
   T& divBy (T& a, T b) const  { return mulBy (a, recip (b)); }

   T power (T a, unsigned k) const;

   unsigned order (T a) const  { return orderTable [a]; }
   bool isPrimitiveElement (T a) const  { return order (a) == size() - 1; }

   void print (std::ostream &o, const type& x) const
      { Priv::printNumberSuffix (o, x, size()); }
   void printShort (std::ostream &o, const type& x) const 
      { Priv::printNumber (o, x); }

   void dump (std::ostream &) const;
   
   // Change tables

   void setMul (T, T, T);
   void setRecip (T, T);
   void setOrder (T, unsigned);
};


/**
 *  Lookup Field
 */

template <typename T>
class LookupField : public LookupFieldBase<T>
{
private:
   T* addTable;
   T* negTable;

   void privSetPointers();

protected:
   void setPointers();

public:
   typedef gf_tag algebra_category;
   typedef T type;

   explicit LookupField (unsigned);
   LookupField (const LookupField<T> &r) : LookupFieldBase<T> (r)
      { privSetPointers(); }
   LookupField & operator= (const LookupField<T> &r)
      { LookupFieldBase<T>::operator=(r); privSetPointers(); return *this; }

   bool is0 (T a) const  { return a == T(0); }

   T  add   (T  a, T b) const  { return addTable [a + size() * b]; }
   T& addTo (T& a, T b) const  { return a = add (a,b); }

   T  neg    (T  a) const  { return negTable [a]; }
   T& negate (T& a) const  { return a = neg(a); }

   T  sub     (T  a, T b) const  { return add   (a, neg (b)); }
   T& subFrom (T& a, T b) const  { return addTo (a, neg (b)); }

   T times (T a, unsigned k) const;

   void dump (std::ostream &) const;
   
   // Change tables

   void setAdd (T, T, T);
   void setNeg (T, T);

   unsigned characteristic() const  { return getCharacteristic(); }
   unsigned extensionDegree() const  { return getExtensionDegree(); }

   HINTLIB_TRIVIAL_DOMAIN_MEMBERS
};


/**
 *  Lookup Field Mul Only
 */

template <typename T>
class LookupFieldMulOnly : public LookupFieldBase<T>
{
protected:
   explicit LookupFieldMulOnly (unsigned _s)
      : LookupFieldBase<T> (_s,  (_s*_s) + 2*_s) {}
   LookupFieldMulOnly (const LookupFieldMulOnly<T> &r)
      : LookupFieldBase<T> (r) {}
};


/**
 *  Lookup Field Pow 2
 */

template <typename T>
class LookupFieldPow2 : public LookupFieldMulOnly<T>,
                        public BitOpBasedAddition<T>
{
public:
   typedef gf_tag algebra_category;
   typedef T type;

   explicit LookupFieldPow2 (unsigned _s) : LookupFieldMulOnly<T> (_s) {}
   LookupFieldPow2 (const LookupFieldPow2<T> &r)
      : LookupFieldMulOnly<T>(r),
        BitOpBasedAddition<T>() {}

   LookupFieldPow2 & operator= (const LookupFieldPow2<T> &r)
      { LookupFieldBase<T>::operator=(r); return *this; }

   unsigned extensionDegree() const  { return getExtensionDegree(); }

   HINTLIB_TRIVIAL_DOMAIN_MEMBERS
};


/**
 *  Lookup Field Prime
 */

template <typename T>
class LookupFieldPrime : public LookupFieldMulOnly<T>
{
public:
   typedef cyclic_tag algebra_category;
   typedef T type;

   explicit LookupFieldPrime (unsigned _s) : LookupFieldMulOnly<T> (_s) {}
   LookupFieldPrime (const LookupFieldPow2<T> &r) : LookupFieldMulOnly<T> (r) {}

   LookupFieldPrime & operator= (const LookupFieldPrime<T> &r)
      { LookupFieldBase<T>::operator=(r); return *this; }

   bool is0(const T& a) const  { return !a; }

   T add (const T& a, const T& b) const
      { unsigned t = unsigned(a) + unsigned(b);
        return     (t>=size()) ? t - size() : t; }
   T& addTo (T& a,    const T& b) const  { return a = add (a, b); }

   T neg (const T& a) const  { return a != 0  ?  size() - a  :  0; }
   T& negate (T& a) const    { if (a != 0) a = size() - a;  return a; }

   T sub (const T& a, const T& b) const
      { int t = int(a) - int(b); return     (t<0) ? t+size() : t; }
   T& subFrom (T& a,  const T& b) const  { return a = sub (a, b); }

   T times (const T& a, unsigned k) const
      { if (k >= size()) k %= size();  return (a * k) % size(); }

   unsigned characteristic () const  { return size(); }

   HINTLIB_TRIVIAL_CYCLIC_MEMBERS
};


/**
 *  Produces a precalculated Galois Field
 */

template<class T>
void makeGaloisField (T &, unsigned base, unsigned exp);

template<typename T>
class LookupGaloisField : public LookupField<T>
{
public:
   LookupGaloisField (unsigned prime, unsigned power);
   explicit LookupGaloisField (unsigned size);
};

template<typename T>
class LookupGaloisFieldPow2 : public LookupFieldPow2<T>
{
public:
   LookupGaloisFieldPow2 (unsigned prime, unsigned power);
   explicit LookupGaloisFieldPow2 (unsigned size);
};

template<typename T>
class LookupGaloisFieldPrime : public LookupFieldPrime<T>
{
public:
   LookupGaloisFieldPrime (unsigned prime, unsigned power);
   explicit LookupGaloisFieldPrime (unsigned size);
};


/**
 *  copy()
 *
 *  Initialize a Lookup Field using some Arithmetic A
 */

template<typename T, class A>
void copy (LookupFieldMulOnly<T> & dest, const A src);
template<typename T, class A>
void copy (LookupField<T> & dest, const A src);


/*****************************************************************************/

/**
 *  Lookup Vector Space Base Base
 */

class LookupVectorSpaceBB : public RefCountingAlgebra
{
public:
   typedef vectorspace_tag algebra_category;
   typedef nopolynomial_tag polynomial_category;

   unsigned size() const  { return s; }
   unsigned dimension() const  { return dim; }

protected:
   LookupVectorSpaceBB (unsigned _s, unsigned _dim, unsigned memory)
      : RefCountingAlgebra (memory), s(_s), dim(_dim) {}
   LookupVectorSpaceBB (const LookupVectorSpaceBB &v)
      : RefCountingAlgebra (v), s(v.s), dim(v.dim) {}
   LookupVectorSpaceBB& operator= (const LookupVectorSpaceBB& v);

   bool operator== (const LookupVectorSpaceBB &)  const;
   void prnSuffix (std::ostream &o, unsigned) const
      { Priv::printSuffix (o, s); }

private:
   unsigned s;
   unsigned dim;
};


/**
 *  Lookup Vector Space Base
 */

template <typename T, typename C>
class LookupVectorSpaceBase : public LookupVectorSpaceBB
{
protected:
   T* mulTable;

   LookupVectorSpaceBase (unsigned aSize, unsigned _dim, unsigned numData);
   LookupVectorSpaceBase (const LookupVectorSpaceBase<T,C> &r)
      : LookupVectorSpaceBB (r) {}

   void dump (std::ostream &, unsigned) const;
   
public:

   typedef T type;
   typedef C scalar_type;

   T element(unsigned i) const  { return T(i); }
   unsigned index (T x) const   { return unsigned (x); }

   T  mul   (T  a, scalar_type l) const  { return mulTable [a + size() * l]; }
   T& scale (T& a, scalar_type l) const  { return a = mul(a,l); }
};


/**
 *  Lookup Vector Space
 */

template <typename T, typename C>
class LookupVectorSpace : public LookupVectorSpaceBase<T,C>
{
public:
   typedef LookupField<C> scalar_algebra;

   class scalar_reference
   {
   public:
      operator C () const { return (*ptr % u) / l; }
      scalar_reference&  operator= (C x)
         {  *ptr = (*ptr % l) + (*ptr / u) * u + x * l; return *this; }
   private:
      scalar_reference (T* _ptr, unsigned _l, unsigned shift)
         : ptr (_ptr), l (_l), u (_l * shift) {}

      T* ptr;
      unsigned l;
      unsigned u;

      friend class LookupVectorSpace<T,C>;
   };

private:
   T* addTable;
   T* negTable;
   scalar_algebra algebra;
   unsigned shift;

   void setPointers();

public:
   LookupVectorSpace (const scalar_algebra&, unsigned);
   LookupVectorSpace (const LookupVectorSpace<T,C> &);
   LookupVectorSpace & operator= (const LookupVectorSpace<T,C> &);

   scalar_algebra getScalarAlgebra() const  { return algebra; }

   template<typename I> void toCoord (T a, I p) const
   {
      for (unsigned i = 0; i != dimension(); ++i, a /= shift)
         *p++ = C (a % shift);
   }
   template<typename I> void fromCoord (T& a, I p) const
   {
      a = 0;
      p += dimension();
      for (unsigned i = 0; i != dimension(); ++i) a = (a * shift) + (*--p);
   }

   C coord (const T& a, unsigned k) const
      { return (a / (powInt(unsigned (shift), k))) % shift; }
   scalar_reference coord (T& a, unsigned k) const
      { return scalar_reference (&a, powInt (unsigned (shift), k), shift); }

   bool is0 (T a) const  { return !a; }

   T  add   (T  a, T b) const  { return addTable [a + size() * b]; }
   T& addTo (T& a, T b) const  { return a = add (a,b); }

   T  neg    (T  a) const  { return negTable [a]; }
   T& negate (T& a) const  { return a = neg(a); }

   T  sub     (T  a, T b) const  { return add   (a, neg (b)); }
   T& subFrom (T& a, T b) const  { return addTo (a, neg (b)); }

   T times (T a, unsigned k) const;

   void dump (std::ostream &) const;
   void print (std::ostream &, T) const;
   void printShort (std::ostream &, T) const;
   void printSuffix (std::ostream &o) const
      { Priv::printSuffix (o, algebra.size()); }

   unsigned additiveOrder (const T& a) const
      { return a ? algebra.characteristic() : 1; }
};

template<typename T, typename C>
inline
std::ostream& operator<< (std::ostream &o, const LookupVectorSpace<T,C> &v)
{
   Priv::printVectorSpaceName (o, v.getScalarAlgebra().size(), v.dimension());
   return o;
}


/**
 *  Lookup Vector Space Pow 2
 */

template <typename T, typename C>
class LookupVectorSpacePow2 : public LookupVectorSpaceBase<T,C>,
                              public BitOpBasedAddition<T>
{
public:
   typedef LookupFieldPow2<C> scalar_algebra;

   class scalar_reference
   {
   public:
      operator C () const { return (*ptr & mask) >> shift; }
      scalar_reference&  operator= (C x)
         {  *ptr = (*ptr & ~mask) | (x << shift); return *this; }
   private:
      scalar_reference (T* _ptr, T _mask, unsigned _shift)
         : ptr (_ptr), mask (_mask), shift (_shift) {}

      T* ptr;
      T  mask;
      unsigned shift;

      friend class LookupVectorSpacePow2<T,C>;
   };

private:
   void setPointers();

   scalar_algebra algebra;
   unsigned shift;
   T mask;

public:

   scalar_algebra getScalarAlgebra() const  { return algebra; }

   LookupVectorSpacePow2 (const scalar_algebra&, unsigned);
   LookupVectorSpacePow2 (const LookupVectorSpacePow2<T,C> &);
   LookupVectorSpacePow2 & operator= (const LookupVectorSpacePow2<T,C> &);

   template<typename I> void toCoord (T a, I p) const
   {
      for (unsigned i = 0; i != dimension(); ++i, a>>=shift)
         *p++ = C (a & mask);
   }
   template<typename I> void fromCoord (T& a, I p) const
   {
      a = 0;
      p += dimension();
      for (unsigned i = 0; i != dimension(); ++i) a = (a << shift) | (*--p);
   }

   C coord (const T& a, unsigned k) const
      { return (a >> (k * shift)) & mask; }
   scalar_reference coord (T& a, unsigned k) const
      { return scalar_reference (&a, mask << (k * shift), k * shift); }

   void dump (std::ostream &o) const
      { LookupVectorSpaceBase<T,C>::dump (o, algebra.size()); }
   void print (std::ostream &, T) const;
   void printShort (std::ostream &, T) const;
   void printSuffix (std::ostream &o) const
      { Priv::printSuffix (o, algebra.size()); }
};

template<typename T, typename C>
inline
std::ostream& operator<< (std::ostream &o, const LookupVectorSpacePow2<T,C> &v)
{
   Priv::printVectorSpaceName (o, v.getScalarAlgebra().size(), v.dimension());
   return o;
}


/**
 *  Vector Space Pow 2
 */

template <typename T>
class VectorSpacePow2 : public BitOpBasedAddition<T>
{
public:
   typedef vectorspace_tag algebra_category;
   typedef nopolynomial_tag polynomial_category;

   typedef LookupFieldPow2<unsigned char> scalar_algebra;
   typedef T type;
   typedef unsigned char scalar_type;

   class scalar_reference
   {
   public:
      operator scalar_type () const { return (*ptr & mask) >> shift; }
      scalar_reference&  operator= (scalar_type x)
         {  *ptr = (*ptr & ~mask) | (x << shift); return *this; }
   private:
      scalar_reference (T* _ptr, T _mask, unsigned _shift)
         : ptr (_ptr), mask (_mask), shift (_shift) {}

      T* ptr;
      T  mask;
      unsigned shift;

      friend class VectorSpacePow2<T>;
   };

private:

   scalar_algebra algebra;
   LookupVectorSpacePow2<scalar_type,scalar_type> vecalg;
   unsigned baseBits;
   unsigned vecBits;
   T baseMask;
   T vecMask;
   unsigned dim;
   unsigned vecDim;

public:

   scalar_algebra getScalarAlgebra() const  { return algebra; }

   VectorSpacePow2 (const scalar_algebra&, unsigned);
   // default copy and assignment just all right

   unsigned size() const;
   unsigned dimension() const  { return dim; }

   template<typename I> void toCoord (T a, I p) const
   {
      for (unsigned i = 0; i != dimension(); ++i, a>>=baseBits)
         *p++ = scalar_type(a & baseMask);
   }
   template<typename I> void fromCoord (T& a, I p) const
   {
      a = 0;
      p += dimension();
      for (unsigned i = 0; i != dimension(); ++i) a = (a << baseBits) | (*--p);
   }

   scalar_type coord (const T& a, unsigned k) const
      { return (a >> (k * baseBits)) & baseMask; }
   scalar_reference coord (T& a, unsigned k) const
      { return scalar_reference (&a, baseMask << (k*baseBits), k*baseBits); }

   T element(unsigned i) const  { return T(i); }
   unsigned index (T x) const   { return unsigned (x); }

   T  mul   (const T& a,       scalar_type l) const;
   T& scale (      T& a, const scalar_type& l) const  { return a = mul(a,l); }

   void print (std::ostream &, const T&) const;
   void printShort (std::ostream &, const T&) const;
   void printSuffix (std::ostream &o) const
      { Priv::printSuffix (o, algebra.size()); }

   HINTLIB_TRIVIAL_DOMAIN_MEMBERS
};

template<typename T>
inline
std::ostream& operator<< (std::ostream &o, const VectorSpacePow2<T> &v)
{
   Priv::printVectorSpaceName (o, v.getScalarAlgebra().size(), v.dimension());
   return o;
}

}  // namespace HIntLib

#endif

