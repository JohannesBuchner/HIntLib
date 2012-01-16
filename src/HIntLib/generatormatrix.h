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


#ifndef HINTLIB_GENERATOR_MATRIX_H
#define HINTLIB_GENERATOR_MATRIX_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <iosfwd>

#include <HIntLib/defaults.h>


namespace HIntLib
{

/**
 *  Generator Matrix
 *
 *  Base class for all Generator Matrices
 *
 *  A Generator Matrix contains the set of matrices used for the generation of
 *  a digital net.
 *
 *  A Generator Matrix contains _dim_ matrices (for dimension 0,..,dim-1).
 *  Each matrix hold _m_*_prec_ entries from a ring with _base_ elements.
 *  The maximum number of points that can be generated from the matrix is given
 *  by  _base_ ^ _m_.
 *  Each point has _prec_ significant base-_base_ digits.
 */

class GeneratorMatrix
{
public:

   // no public constructor!
   virtual ~GeneratorMatrix() {};

   unsigned getBase() const          { return base; }
   unsigned getDimension() const     { return dim; }
   unsigned getM() const             { return m; }

   unsigned getVectorization() const { return vec; }
   unsigned getPrec() const          { return prec; }      // # vecBase digits
   unsigned getTotalPrec() const     { return totalPrec; } // # base digits
   unsigned getNumOfLeadingDigits() const
      { return totalPrec - (vec * (prec - 1)); }
   unsigned getNumOfMissingDigits() const { return prec * vec - totalPrec; }

   virtual void     setDigit  (unsigned d, unsigned r, unsigned b, unsigned x);
   virtual void     setVector (unsigned d, unsigned r, unsigned b, u64 x);
   virtual unsigned getDigit  (unsigned d, unsigned r, unsigned b) const = 0;
   virtual u64      getVector (unsigned d, unsigned r, unsigned b) const = 0;

   void dump       (std::ostream &) const;
   void vectorDump (std::ostream &) const;
   void libSeqDump (std::ostream &) const;
   void libSeqDump (const char *) const;
   void binaryDump (std::ostream &) const;
   void binaryDump (const char *) const;

   // checkXXX(). Only used by assign().
   
   static void checkCopy (const GeneratorMatrix &, const GeneratorMatrix &);
   static void checkCopyDim (
         const GeneratorMatrix &, const GeneratorMatrix &, int offset = 0);

protected:
   GeneratorMatrix (unsigned _base, unsigned _vec,
                    unsigned _dim, unsigned _m, unsigned _totalPrec)
      : base (_base),
        vec (_vec),
        prec ((_totalPrec-1) / _vec + 1),
        dim (_dim),
        m   (_m),
        totalPrec (_totalPrec) {}
   GeneratorMatrix (const GeneratorMatrix &);

   const unsigned base;
   const unsigned vec;
   const unsigned prec;
   const unsigned dim;
   const unsigned m;
         unsigned totalPrec;
private:
   GeneratorMatrix& operator= (const GeneratorMatrix&);
};

bool operator== (const GeneratorMatrix &, const GeneratorMatrix &);


/**
 *  assign()
 */

void assign (const GeneratorMatrix &, unsigned, GeneratorMatrix &, unsigned);
void assign (const GeneratorMatrix &, GeneratorMatrix &);


/**
 *  GM Copy
 */

class GMCopy
{
public:
   GMCopy ();

   GMCopy& dim (int x)          { dimValue = x; return *this; }
   GMCopy& maxDim (int x)       { maxDimValue = x; return *this; }
   GMCopy& m (int x)            { mValue = x; return *this; }
   GMCopy& maxM (int x)         { maxMValue = x; return *this; }
   GMCopy& totalPrec (int x)    { totalPrecValue = x; return *this; }
   GMCopy& maxTotalPrec (int x) { maxTotalPrecValue = x; return *this; }
   GMCopy& vec (int x)   { vecValue = x; return *this; }
   GMCopy& noVec ()      { vecValue = 1; return *this; }
   GMCopy& keepVec()     { vecValue = -2; return *this; }
   // GMCopy& optimalVec()  { vecValue = -3; return *this; }
   GMCopy& equi()        { equiValue = true; return *this; }
   GMCopy& equi(bool e)  { equiValue = e; return *this; }

   unsigned getDimension (const GeneratorMatrix &) const;
   unsigned getM         (const GeneratorMatrix &) const;
   unsigned getTotalPrec (const GeneratorMatrix &, unsigned max = 0) const;
   unsigned getVectorization (const GeneratorMatrix &, unsigned bits) const;
   bool getEqui () const  { return equiValue; }
   void checkNoVec () const;

private:
   int dimValue;
   int maxDimValue;
   int mValue;
   int maxMValue;
   int totalPrecValue;
   int maxTotalPrecValue;
   int vecValue;
   bool equiValue;
};

}  // namespace HIntLib

#endif

