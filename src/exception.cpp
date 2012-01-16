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
 *  exception.h
 *
 *  Defines the exceptions used by HIntLib
 */

#ifdef __GNUG__
#pragma implementation
#endif

#include <string>

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_CSTRING
  #include <cstring>
  #define HINTLIB_SSN std::
#else
  #include <string.h>
  #define HINTLIB_SSN
#endif

#ifdef HINTLIB_HAVE_SSTREAM
  #include <sstream>
#else
  #include <HIntLib/fallback_sstream.h>
#endif

#include <HIntLib/exception.h>

#include <HIntLib/hypercube.h>

using std::ostringstream;

namespace L = HIntLib;

namespace
{
   char* ms (ostringstream &ss)
   {
      std::string s = ss.str();

      char* p = new char [s.length() + 1];
      HINTLIB_SSN strcpy (p, s.c_str());

      return p;
   }
}

/**
 *  Copy Constructor
 */

L::Exception::Exception(const Exception &e)
   : std::exception(e), s(0)
{
   if (e.s)  setStringCopy (e.s);
}


/**
 *  setStringCopy()
 *
 *  Sets s. New memory is allocated an the string is copied.
 */

void L::Exception::setStringCopy (const char*ss) const
{
   delete[] s;
   s = new char [strlen(ss) + 1];
   HINTLIB_SSN strcpy (s, ss);
}


/**
 *  setString()
 *
 *  Sets the string s to the given pointer.
 *
 *  ss has to be pointed to a char[] on the free store which will eventually
 *  be deallcoated by the destructor of Exception.
 */

inline void L::Exception::setString(char*ss) const
{
   s = ss;
}

const char* L::Exception::what() const throw()
{
   if (!s)  makeString();

   return s;
}

// Invalid Dimension

void L::InvalidDimension::makeString() const
{
   ostringstream ss;
   ss << "Dimension " << dim << " is invalid!";
   setString (ms(ss));
}

void L::DimensionZero::makeString() const
{
   setStringCopy ("Dimension must not be zero!");
}

void L::checkDimensionNotZero (unsigned dim)
{
   if (dim == 0)  throw DimensionZero();
}

void L::DimensionTooHigh::makeString() const
{
   ostringstream ss;
   ss << "Dimension " << getDimension() << " is too high!  Maximum is "
      << max << ".";
   setString (ms(ss));
}

void L::DimensionTooLow::makeString() const
{
   ostringstream ss;
   ss << "Dimension " << getDimension() << " is too low!  Minimum is "
      << min << ".";
   setString (ms(ss));
}

// Dimension Mismatch

void L::DimensionMismatch::makeString() const
{
   ostringstream ss;
   ss << "Dimension of objects does not match (" << dim1 << "!=" << dim2 << "!";
   setString (ms(ss));
}

void L::throwDimensionMismatch (unsigned dim1, unsigned dim2)
{
   throw DimensionMismatch(dim1, dim2);
}

// Integrator Exceptions

void L::NoEvaluationsPossible::makeString() const
{
   ostringstream ss;
   ss << "Cannot apply this Integrator with only " << n << " abscissas!";
   setString (ms(ss));
}

void L::MaxEvaluationsRequired::makeString() const
{
   setStringCopy ("Integrator needs upper bound on the number of abscissas!");
}

// Generator Matrix Exceptions

void L::GM_PrecTooHigh::makeString() const
{
   ostringstream ss;
   ss << "Cannot create Generator Matrix with " << prec << " base-" << b
      << " digits precision! Maximum are " << max << " digits.";
   setString (ms(ss));
}

void L::GM_BaseTooLarge::makeString() const
{
   ostringstream ss;
   ss << "Cannot create Generator Matrix with base " << base
      << ". Maxium is " << max << "!";
   setString (ms(ss));
}

void L::GM_CopyPrec::makeString() const
{
   ostringstream ss;
   ss << "Trying to create Generator Matrix with precision " << n
      << " from matrix with precision " << o << '!';
   setString (ms(ss));
}

void L::GM_CopyM::makeString() const
{
   ostringstream ss;
   ss << "Trying to create Generator Matrix with M=" << n
      << " from matrix with m=" << o << '!';
   setString (ms(ss));
}

void L::GM_CopyDim::makeString() const
{
   ostringstream ss;
   ss << "Trying to create Generator Matrix for dimension " << n
      << " from a matrix for dimension " << o << '!';
   setString (ms(ss));
}

void L::GM_CopyBase::makeString() const
{
   ostringstream ss;
   ss << "Trying to create Generator Matrix for base " << n
      << " from a matrix for base " << o << '!';
   setString (ms(ss));
}


// Sequence Exceptions

void L::DigitalNetTooLarge::makeString() const
{
   ostringstream ss;
   ss << "Cannot create Digital Net with " << b << '^' << m << " points! " 
         "Maximum is " << b << '^' << max << '.';
   setString (ms(ss));
}

void L::NumbersExhausted::makeString() const
{
   setStringCopy ("All numbers of a Sequence are exhausted!");
}

void L::NetIndexTooHigh::makeString() const
{
   ostringstream ss;
   ss << "Trying to get Digital Net #" << requestedIndex
      << " with m = " << newM << " from a net with size m = " << originalM
      << '!';
   setString (ms(ss));
}

// Arithmetic

void L::Overflow::makeString() const
{
   setStringCopy ("Arithmetic Overflow!");
}

void L::throwOverflow()
{
   throw Overflow();
}

void L::DivisionByZero::makeString() const
{
   setStringCopy ("Division by Zero!");
}

void L::throwDivisionByZero()
{
   throw DivisionByZero();
}

void L::InvalidLogBase::makeString() const
{
   ostringstream ss;
   ss << "Cannot determine base " << n << " logarithm!";
   setString (ms(ss));
}

void L::GaloisFieldExponent::makeString() const
{
   setStringCopy ("There is no GaloisField with size p^0!");
}

void L::InvalidModularFieldSize::makeString() const
{
   ostringstream ss;
   ss << "Modulus " << n << " is invalid for ModularIntegerField. Must be prime!";
   setString (ms(ss));
}

void L::NotAPrimePower::makeString() const
{
   ostringstream ss;
   ss << n << " is not a prime power!";
   setString (ms(ss));
}

void L::LookupFieldCopy::makeString() const
{
   setStringCopy ("Error while copying an arithmetic into a PrecalculatedField!");
}

void L::LookupFieldSet::makeString() const
{
   ostringstream ss;
   ss << "Index " << index << " is too large for field of size " << size << '!';
   setString (ms(ss));
}

// Prime Number exceptions

void L::PrimeNumberNth::makeString() const
{
   ostringstream ss;
   ss << "Cannot determine the " << getN() << "th prime number!";
   setString (ms(ss));
}

// Not supported functionality

void L::SetIndexNotSupported::makeString() const
{
   setStringCopy ("setIndex() not supported!");
}

void L::GetExactResultNotSupported::makeString() const
{
   setStringCopy ("getExactResult() not supported!");
}

void L::DerivativeNotSupported::makeString() const
{
   setStringCopy ("derivative() not implemented for this Integrand!");
}

void L::MCRoutinesCreateNotSupported::makeString() const
{
   setStringCopy ("MCRoutinesCreate does not support this operation!");
}

// Internal error

void L::InternalError::makeString() const
{
   ostringstream ss;
   ss << "Internal error in file \'" << fname << "\', line " << line << "!";
   setString (ms(ss));
}

void L::InvalidType::makeString() const
{
   ostringstream ss;
   ss << "Invalid type to instantiate template " << className << "<>!";
   setString (ms(ss));
}


/**
 *  invalid Argument ()
 */

#ifndef HINTLIB_UNREACHABLE_CALLS_REMOVED
void L::InvalidArgument::makeString() const
{
   ostringstream ss;
   ss << "Invalid argument used in " << name << " (this should have been detected at link-time)!";
   setString (ms(ss));
}

void L::invalidArgument (const char* s)
{
   throw InvalidArgument (s);
}
#endif


// RRNGs

void L::BuiltInPRNGUsedTwice::makeString() const
{
   setStringCopy ("BuiltInPRNG constructed twice!");
}

void L::throwInvalidLCGParameters ()
{
   throw InvalidLCGParameters();
}

void L::InvalidLCGParameters::makeString() const
{
   setStringCopy ("Parameters for LCG are invalid!");
}

// Line Reader Exception

L::LineReaderException::LineReaderException
   (unsigned _ln, unsigned _pos, const char* _line, const char* _msg)
   : ln (_ln), pos(_pos),
     line (new char [strlen(_line) + 1]),
     msg (new char [strlen(_msg ? _msg : "Error") + 1])
{
   HINTLIB_SSN strcpy (line, _line);
   HINTLIB_SSN strcpy (msg, _msg ? _msg : "Error");
}

L::LineReaderException::LineReaderException (const LineReaderException &lr)
   : Exception (lr),
     ln (lr.ln), pos(lr.pos),
     line (new char [strlen(lr.line) + 1]),
     msg  (new char [strlen(lr.msg)  + 1])
{
   HINTLIB_SSN strcpy (line, lr.line);
   HINTLIB_SSN strcpy (msg, lr.msg);
}

L::LineReaderException::~LineReaderException () throw ()
{
   delete (line);
   delete (msg);
}

void L::LineReaderException::makeString() const
{
   ostringstream ss;
   ss << msg << " in line " << ln << ", " << pos+1 << ": " << line;
   setString (ms(ss));
}

#undef HINTLIB_SSN


