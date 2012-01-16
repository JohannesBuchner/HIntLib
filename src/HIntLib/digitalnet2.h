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
 *  Digital Net 2
 *  Digital Net 2 <T>
 *
 *  Constracts digital nets (and sequences) in base 2 based on GeneratorMatrix2
 */

#ifndef DIGITAL_NET_2_H
#define DIGITAL_NET_2_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/generatormatrix.h>
#include <HIntLib/qrnsequencebase.h>
#include <HIntLib/shiftscale.h>
#include <HIntLib/pointset.h>
#include <HIntLib/qmcintegration.h>

namespace HIntLib
{
#ifdef HINTLIB_IEEE_MAGIC_WORKS
   template<class X> struct FloatType {};
   template<> struct FloatType<u32> {  typedef float  floatType; };
   #ifdef HINTLIB_U32_NOT_EQUAL_U64
   template<> struct FloatType<u64> {  typedef double floatType; };
   #endif
#endif


/**
 *  Digital Net 2
 *
 *  A Digital Net in base 2
 *
 *  T .. unsigned integer type at least  _prec_  bits
 *         usually either "u32" or "u64"
 */

template<class T>
class DigitalNet2: public QRNSequenceBase, public DigitalNet
{
public:
   const GeneratorMatrix2<T> & getGeneratorMatrix() const  { return c; }
   Index getOptimalNumber(Index max) const;

   void setCube (const Hypercube &);
   void randomize (PRNG &);

protected:
   const unsigned prec;
   GeneratorMatrix2Copy<T> c;
   Array<T> x;      // current vector (size dim)
   Array<T> xStart; // Inital values for x (size dim)
   const Truncation trunc;
   const real trivialScale;
   ShiftScale ss;
#ifdef HINTLIB_IEEE_MAGIC_WORKS
   ShiftScale ssMagic;
   enum MODE {STANDARD, DIRECT} mode;
   typedef typename FloatType<T>::floatType floatType;
   static const unsigned floatBits = std::numeric_limits<floatType>::digits - 1;
#endif

   DigitalNet2
      (const GeneratorMatrix2<T> &, const Hypercube &,
       unsigned m, Index index, bool equi, Truncation);

   void copyXtoP          (real*);
   void copyXtoPDontScale (real*);

   void resetX          (Index n, real*p) { resetX(n); copyXtoP         (p); }
   void resetXDontScale (Index n, real*p) { resetX(n); copyXtoPDontScale(p); }

private:
   void resetX (Index);
};


/**
 *  Digital Net 2 Naive
 *
 *  Traverses the net in normal order - naive implementation
 */

template<class T>
class DigitalNet2Naive : public DigitalNet2<T>
{
public:
   DigitalNet2Naive
      (const GeneratorMatrix2<T> &gm, const Hypercube &_h,
       unsigned m, Index i, bool equi, DigitalNet::Truncation t)
      : DigitalNet2<T> (gm, _h, m, i, equi, t) {}

   DigitalNet2Naive
      (const GeneratorMatrix2<T> &gm, const Hypercube &_h)
      : DigitalNet2<T> (gm, _h, gm.getM(), 0, false, FULL) {}

   void first          (real *p, Index _n = 0)  { resetX          (n = _n, p); }
   void firstDontScale (real *p, Index _n = 0)  { resetXDontScale (n = _n, p); }

   void next           (real *p)  { resetX          (++n, p); }
   void nextDontScale  (real *p)  { resetXDontScale (++n, p); }
};


/**
 *  Digital Net 2 Gray
 *
 *  Traverses the net in gray-code order
 */

template<class T>
class DigitalNet2Gray : public DigitalNet2<T>
{
public:
   DigitalNet2Gray
      (const GeneratorMatrix2<T> &, const Hypercube &,
       unsigned m, Index i, bool equi, DigitalNet::Truncation,
       bool correct = true);

   DigitalNet2Gray
      (const GeneratorMatrix2<T> &, const Hypercube &, bool correct = true);

   void first          (real*p, Index _n = 0)
      { resetX          (grayCode(n = _n), p); }
   void firstDontScale (real*p, Index _n = 0)
      { resetXDontScale (grayCode(n = _n), p); }

   void next          (real*);
   void nextDontScale (real*);
};



/**
 *  Digital Net 2 Point Set
 */

// Non-template base

class DigitalNet2PointSetBase : public PartitionablePointSet
{
public:
#ifdef HINTLIB_IEEE_MAGIC_WORKS
   typedef u64 BaseType;
#else
   typedef Index BaseType;
#endif
   typedef DigitalNet::Truncation Truncation;

protected:
   const GeneratorMatrix2Copy<BaseType> gm;
   bool copy;
   bool equi;
   Truncation trunc;
   Index index;
   const Hypercube *h;

   int calculateM (Index) const;

   DigitalNet2PointSetBase (
      const GeneratorMatrix2<u32>& _gm, bool _equi, Truncation t, Index i)
      : gm(_gm), equi(_equi), trunc(t), index (i), h (0)  {}
#ifdef HINTLIB_U32_NOT_EQUAL_U64
   DigitalNet2PointSetBase (
      const GeneratorMatrix2<u64>& _gm, bool _equi, Truncation t, Index i)
      : gm(_gm), equi(_equi), trunc(t), index (i), h (0)  {}
#endif

public:
   void setCube (const Hypercube*);

   bool doJobRep (real *, ReportingJob &, Index);
   void doJobPartition (real *, Job &, Index, Index, Index);

   Index getOptimalNumber (Index max, const Hypercube &);

private:
   DigitalNet2PointSetBase (const DigitalNet2PointSetBase &);
   const DigitalNet2PointSetBase & operator= (const DigitalNet2PointSetBase &);
};

// Template

template<class Sum>
class DigitalNet2PointSet: public DigitalNet2PointSetBase
{
public:
   DigitalNet2PointSet (
      const GeneratorMatrix2<u32>& _gm,
      bool e = true, DigitalNet::Truncation t = DigitalNet::TRUNCATE,
      Index i = 0)
   : DigitalNet2PointSetBase(_gm, e, t, i) {}
#ifdef HINTLIB_U32_NOT_EQUAL_U64
   DigitalNet2PointSet (
      const GeneratorMatrix2<u64>& _gm,
      bool e = true, DigitalNet::Truncation t = DigitalNet::TRUNCATE,
      Index i = 0)
   : DigitalNet2PointSetBase(_gm, e, t, i) {}
#endif

   void integratePartition (real *, Function &, Index, Index, Index, Stat&);
   void integratePartition (real *, Function &, Index, Index, Index, StatVar&);
};

// Specializatioin for Sum=real

template<>
class DigitalNet2PointSet<real>: public DigitalNet2PointSetBase
{
public:
   DigitalNet2PointSet<real> (
      const GeneratorMatrix2<u32>& gm,
      bool e = true, DigitalNet::Truncation t = DigitalNet::TRUNCATE,
      Index i = 0)
   : DigitalNet2PointSetBase(gm, e, t, i) {}
#ifdef HINTLIB_U32_NOT_EQUAL_U64
   DigitalNet2PointSet<real> (
      const GeneratorMatrix2<u64>& gm,
      bool e = true, DigitalNet::Truncation t = DigitalNet::TRUNCATE,
      Index i = 0)
   : DigitalNet2PointSetBase(gm, e, t, i) {}
#endif

   void integratePartition (real *, Function &, Index, Index, Index, Stat&);
   void integratePartition (real *, Function &, Index, Index, Index, StatVar&);
};


/**
 *  Digital Sequence 2 Point Set
 */

// Non-template base

class DigitalSeq2PointSetBase : public PartitionablePointSet
{
public:
#ifdef HINTLIB_IEEE_MAGIC_WORKS
   typedef u64 BaseType;
#else
   typedef Index BaseType;
#endif
   void init()  { delete net; net = 0; offset = 0; }

protected:
   GeneratorMatrix2Copy<BaseType> gm;
   DigitalNet2Gray<BaseType>* net;
   Index offset;
   const bool reset;

   void checkSize (const DigitalNet2<BaseType> &, Index) const;

   DigitalSeq2PointSetBase (const GeneratorMatrix2<u32> &_gm, bool _reset)
      : gm (_gm), net (0), offset (0), reset (_reset) {}
#ifdef HINTLIB_U32_NOT_EQUAL_U64
   DigitalSeq2PointSetBase (const GeneratorMatrix2<u64> &_gm, bool _reset)
      : gm (_gm), net (0), offset (0), reset (_reset) {}
   ~DigitalSeq2PointSetBase ();
#endif

public:
   Index getOptimalNumber (Index max, const Hypercube &);
   void setCube (const Hypercube *);

   bool doJobRep       (real *, ReportingJob &, Index);
   void doJobPartition (real *, Job &, Index, Index, Index);

private:
   DigitalSeq2PointSetBase (const DigitalNet2PointSetBase &);
   const DigitalSeq2PointSetBase & operator= (const DigitalNet2PointSetBase &);
};

// Template

template<class Sum>
class DigitalSeq2PointSet: public DigitalSeq2PointSetBase
{
public:
   DigitalSeq2PointSet (const GeneratorMatrix2<u32>& gm, bool reset)
      : DigitalSeq2PointSetBase(gm, reset) {}
#ifdef HINTLIB_U32_NOT_EQUAL_U64
   DigitalSeq2PointSet (const GeneratorMatrix2<u64>& gm, bool reset)
      : DigitalSeq2PointSetBase(gm, reset) {}
#endif

   void integratePartition (real *, Function &, Index, Index, Index, Stat&);
   void integratePartition (real *, Function &, Index, Index, Index, StatVar&);
};

// Specializatioin for Sum=real

template<>
class DigitalSeq2PointSet<real>: public DigitalSeq2PointSetBase
{
public:
   DigitalSeq2PointSet<real> (const GeneratorMatrix2<u32>& gm, bool reset)
      : DigitalSeq2PointSetBase(gm, reset) {}
#ifdef HINTLIB_U32_NOT_EQUAL_U64
   DigitalSeq2PointSet<real> (const GeneratorMatrix2<u64>& gm, bool reset)
      : DigitalSeq2PointSetBase(gm, reset) {}
#endif

   void integratePartition (real *, Function &, Index, Index, Index, Stat&);
   void integratePartition (real *, Function &, Index, Index, Index, StatVar&);
};

} // namespace HIntLib


/*****  Implementation  **************/


/**
 *  copy X to P
 */

template<class T>
inline
void HIntLib::DigitalNet2<T>::copyXtoP (real* point)
{
#ifdef HINTLIB_IEEE_MAGIC_WORKS
   if (      std::numeric_limits<floatType>::digits
          >= std::numeric_limits<real>::digits
       || mode == DIRECT)
   {
      floatType *p = reinterpret_cast<floatType*>(&x[0]);
      for (unsigned d = 0; d < getDimension(); ++d)
      {
         point[d] = ssMagic[d] (p[d] + -1.0);
      }
   }
   else
#endif
   for (unsigned d = 0; d < getDimension(); ++d)  point[d] = ss[d] (x[d]);
}


/**
 *  copy X to P Dont Scale ()
 */

template<class T>
inline
void HIntLib::DigitalNet2<T>::copyXtoPDontScale (real* point)
{
#ifdef HINTLIB_IEEE_MAGIC_WORKS
   if (      std::numeric_limits<floatType>::digits
          >= std::numeric_limits<real>::digits
       || mode == DIRECT)
   {
      floatType *p = reinterpret_cast<floatType*>(&x[0]);
      for (unsigned d = 0; d < getDimension(); ++d)  point[d] = p[d] + -1.0;
   }
   else
#endif
   for (unsigned d = 0; d < getDimension(); ++d)
   {
      point[d] = x[d] * trivialScale;
   }
}


/**
 *  Gray :: next()
 *
 *  The next method updates the array x depending on n and the vectors in v
 *  Once x is updated, these values are converted to real numbers and returned.
 */

template<class T>
inline
void HIntLib::DigitalNet2Gray<T>::next (real* point)
{
   const T *vp = c (ls0 (n++));  // determin bit that has to be changed

   for (unsigned d = 0; d < getDimension(); ++d)  x[d] ^= vp[d];

   copyXtoP (point);
} 

template<class T>
inline
void HIntLib::DigitalNet2Gray<T>::nextDontScale (real* point)
{
   const T *vp = c (ls0 (n++));  // determin bit that has to be changed

   for (unsigned d = 0; d < getDimension(); ++d)  x[d] ^= vp[d];

   copyXtoPDontScale (point);
} 


/**
 *  DigitalNet2PointSet :: integrate()
 */

template<class Sum>
void HIntLib::DigitalNet2PointSet<Sum>::integratePartition
   (real* point, Function &f, Index num, Index begin, Index end, Stat& stat)
{
   DigitalNet2Gray<BaseType> net
      (gm, *h, calculateM (num), index, equi, trunc);

   Statistic<real, Sum> s;
   qmcIntegration (point, net, f, begin, end, s);
   stat = s;
}

template<class Sum>
void HIntLib::DigitalNet2PointSet<Sum>::integratePartition
   (real* point, Function &f, Index num, Index begin, Index end, StatVar& stat)
{
   DigitalNet2Gray<BaseType> net
      (gm, *h, calculateM (num), index, equi, trunc);

   StatisticVar<real, Sum> s;
   qmcIntegration (point, net, f, begin, end, s);
   stat = s;
}


/**
 *  DigitalSeq2PointSet :: integrate()
 */

template<class Sum>
void HIntLib::DigitalSeq2PointSet<Sum>::integratePartition
   (real* point, Function &f, Index num, Index begin, Index end, Stat& stat)
{
   if (! reset)
   {
      begin  += offset;
      end    += offset;
      offset += num;
   }
   checkSize (*net, end);

   Statistic<real, Sum> s;
   qmcIntegration (point, *net, f, begin, end, s);
   stat = s;
}

template<class Sum>
void HIntLib::DigitalSeq2PointSet<Sum>::integratePartition
   (real* point, Function &f, Index num, Index begin, Index end, StatVar& stat)
{
   if (! reset)
   {
      begin  += offset;
      end    += offset;
      offset += num;
   }
   checkSize (*net, end);

   StatisticVar<real, Sum> s;
   qmcIntegration (point, *net, f, begin, end, s);
   stat = s;
}

#endif

