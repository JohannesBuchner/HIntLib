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

/**
 *  exception.h
 *
 *  Defines the exceptions used by HIntLib
 */

#ifndef HINTLIB_EXCEPTION_H
#define HINTLIB_EXCEPTION_H 1
 
#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <exception>

namespace HIntLib
{
   // This function is defend nowhere.
   // It is called for invalid template arguments, causing a link-time error

   void invalidArgument (const char*);

   // Every exception defined here is derived from Exception

   class Exception : public std::exception
   {
   public:
      ~Exception() throw()  { delete[] s; }
      virtual const char* what() const throw();
   protected:
      Exception() : s(0) {}
      Exception(const Exception &);
      void setStringCopy (const char*) const;
      void setString (char*) const;
   private:
      Exception& operator=(const Exception &);    // no assignment
      virtual void makeString() const = 0;
      mutable char* s;
   };

   // Invalid Dimension

   class RangeException : public Exception {};

   class InvalidDimension : public RangeException
   {
      virtual void makeString() const;
   public:
      InvalidDimension (int d) : dim(d) {}
      int getDimension() const { return dim; }
   private:
      const int dim;
   };

   class DimensionZero : public InvalidDimension
   {
      virtual void makeString() const;
   public:
      DimensionZero() : InvalidDimension(0) {}
   };

   void checkDimensionNotZero (int dim);

   class DimensionTooHigh : public InvalidDimension
   {
      virtual void makeString() const;
   public:
      DimensionTooHigh(int d, int m) : InvalidDimension(d), max (m) {}
      int getMaximum() const  {  return max; }
   private:
      const int max;
   };

   template <int max>
   inline
   void checkDimensionLeq (int dim)
   {
      if (dim > max)  throw DimensionTooHigh (dim, max);
   }

   class DimensionTooLow : public InvalidDimension
   {
      virtual void makeString() const;
   public:
      DimensionTooLow(int d, int m) : InvalidDimension(d), min (m) {}
      int getMinimum() const  {  return min; }
   private:
      const int min;
   };

   template <int min>
   inline
   void checkDimensionGeq (int dim)
   {
      if (dim < min)  throw DimensionTooLow (dim, min);
   }

   // Dimension Mismatch

   class DimensionMismatch : public Exception
   {
      virtual void makeString() const;
   public:
      DimensionMismatch (int d1, int d2) : dim1(d1), dim2(d2) {}
      int getDimension1() const  { return dim1; }
      int getDimensino2() const  { return dim2; }
   private:
      const int dim1;
      const int dim2;
   };

   void throwDimensionMismatch(int dim1, int dim2) HINTLIB_GNU_NORETURN;

   inline
   void checkDimensionEqual (int dim1, int dim2)
   {
      if (dim1 != dim2)  throwDimensionMismatch (dim1, dim2);
   }

   // Integrator Exceptions

   class IntegratorException : public Exception {};

   class NoEvaluationsPossible : public IntegratorException
   {
      virtual void makeString() const;
      const Index n;
   public:
      NoEvaluationsPossible (Index nn) : n(nn) {}
      Index getN() const  { return n; }
   };

   class RequestedErrorNegative : public IntegratorException
   {
      virtual void makeString() const;
      const real error;
   public:
      RequestedErrorNegative (real e) : error(e) {}
      real getError() const  { return error; }
   };

   class MaxEvaluationsRequired : public IntegratorException
   {
      virtual void makeString() const;
   };

   class TerminationCriterionMissing : public IntegratorException
   {
      virtual void makeString() const;
   };

   // Generator Matrix Exceptions

   class GeneratorMatrixException : public Exception {};

   class GM_PrecTooHigh : public GeneratorMatrixException
   {
      virtual void makeString() const;
      const int b;
      const int max;
      const int prec;
   public:
      GM_PrecTooHigh (int _b, int _max, int _prec)
         : b(_b), max(_max), prec(_prec)  {}
   };

   class GM_BaseTooLarge : public GeneratorMatrixException
   {
      virtual void makeString() const;
      const int base;
      const int max;
   public:
      GM_BaseTooLarge (int _base, int _max)
         : base(_base), max(_max)  {}
   };

   class GM_CopyPrec : public GeneratorMatrixException
   {
      virtual void makeString() const;
      const int n;
      const int o;
   public:
      GM_CopyPrec (int _n, int _o) : n(_n), o(_o)  {}
   };

   class GM_CopyM : public GeneratorMatrixException
   {
      virtual void makeString() const;
      const int n;
      const int o;
   public:
      GM_CopyM (int _n, int _o) : n(_n), o(_o)  {}
   };

   class GM_CopyDim : public GeneratorMatrixException
   {
      virtual void makeString() const;
      const int n;
      const int o;
   public:
      GM_CopyDim (int _n, int _o) : n(_n), o(_o)  {}
   };

   class GM_CopyBase: public GeneratorMatrixException
   {
      virtual void makeString() const;
      const int n;
      const int o;
   public:
      GM_CopyBase (int _n, int _o) : n(_n), o(_o)  {}
   };

   // Sequence Exceptions

   class SequenceException : public Exception {};

   class DigitalNetTooLarge : public SequenceException
   {
      virtual void makeString() const;
      const int b;
      const int max;
      const int m;
   public:
      DigitalNetTooLarge (int _b, int _max, int _m)
         : b(_b), max(_max), m(_m)  {}
   };

   class NumbersExhausted : public SequenceException
   {
      virtual void makeString() const;
   };

   class NetIndexTooHigh : public SequenceException
   {
      const Index requestedIndex;
      const int originalM;
      const int newM;
      virtual void makeString() const;
   public:
      NetIndexTooHigh (Index i, int _originalM, int _newM)
         : requestedIndex (i), originalM(_originalM), newM(_newM)  {}
   };

   // Arithmetic

   class ArithmeticException : public Exception {};

   class Overflow : public ArithmeticException
   {
      virtual void makeString() const;
   };

   void throwOverflow();

   class DivisionByZero : public ArithmeticException
   {
      virtual void makeString() const;
   };

   void throwDivisionByZero ();

   template<typename T>
   void throwDivisionByZero (const T& x)
   {
      if (! x)  throwDivisionByZero();
   }

   class InvalidLogBase : public ArithmeticException
   {
      virtual void makeString() const;
   public:
      InvalidLogBase (real _n) : n (_n)  {}
   private:
      real n;
   };

   class InvalidModularFieldSize : public ArithmeticException
   {
      virtual void makeString() const;
   public:
      InvalidModularFieldSize (unsigned _n) : n (_n)  {}
   private:
      unsigned n;
   };

   class NotAPrimePower : public ArithmeticException
   {
      virtual void makeString() const;
   public:
      NotAPrimePower (unsigned _n) : n (_n)  {}
   private:
      unsigned n;
   };

   class GaloisFieldExponent : public ArithmeticException
   {
      virtual void makeString() const;
   };

   class LookupFieldException : public ArithmeticException {};

   class LookupFieldCopy : public LookupFieldException
   {
      virtual void makeString() const;
   };

   class LookupFieldSet : public LookupFieldException
   {
      virtual void makeString() const;
   public:
      LookupFieldSet (unsigned _index, unsigned _size)
         : index (_index), size (_size)  {}
   private:
      const unsigned index;
      const unsigned size;
   };

   // PrimeNumberExceptions

   class PrimeNumberNth : public RangeException 
   {
      virtual void makeString() const;
   public:
      PrimeNumberNth (unsigned _n) : n (_n)  {}
      unsigned getN() const  { return n; }
   private:
      const unsigned n;
   };

   // Not supported Functionality

   class NotSupportedException : public Exception {};

   class SetIndexNotSupported : public NotSupportedException
   {
      virtual void makeString() const;
   }; 

   class GetExactResultNotSupported : public NotSupportedException
   {
      virtual void makeString() const;
   }; 

   class DerivativeNotSupported : public NotSupportedException
   {
      virtual void makeString() const;
   };

   class MCRoutinesCreateNotSupported : public NotSupportedException
   {
      virtual void makeString() const;
   };

   // Other exceptions

   class OtherException : public Exception
   {
      virtual void makeString() const;
   public:
      OtherException (const char* _msg) : msg(_msg) {}
   private:
      const char* msg;
   };

   void checkPercentageRange (double v);

   // Internal Error

   class InternalError : public Exception
   {
      virtual void makeString() const;
   public:
      InternalError (const char *fn, unsigned l) : fname (fn), line (l) {}
   private:
      const char* fname;
      const unsigned line;
   };

   class FIXME : public InternalError
   {
   public:
      FIXME (const char*s, unsigned line) : InternalError (s, line) {}
   };

   class InvalidType : public Exception
   {
      virtual void makeString() const;
   public:
      InvalidType (const char* _s) : className (_s)  {}
   private:
      const char* className;
   };

#ifndef HINTLIB_UNREACHABLE_CALLS_REMOVED
   class InvalidArgument : public Exception
   {
      virtual void makeString() const;
   public:
      InvalidArgument (const char* _s) : name (_s)  {}
   private:
      const char* name;
   };
#endif

   // PRNG

   class BuiltInPRNGUsedTwice : public Exception
   {
      virtual void makeString() const;
   public:
      BuiltInPRNGUsedTwice () {}
   };

   void throwInvalidLCGParameters () HINTLIB_GNU_NORETURN;

   class InvalidLCGParameters : public Exception
   {
      virtual void makeString () const;
   public:
      InvalidLCGParameters () {}
   };

   // Line Reader Exception

   class LineReaderException : public Exception
   {
      virtual void makeString () const;
   public:
      LineReaderException (unsigned, unsigned, const char*, const char* = 0);
      LineReaderException (const LineReaderException &);
      ~LineReaderException() throw();
   private:
      unsigned ln;
      unsigned pos;
      char* line;
      char* msg;
   };

}  // namespace HIntLib

#endif

