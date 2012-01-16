/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schürer <rudolf.schuerer@sbg.ac.at>
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
 *  Make
 */

#ifdef __GNUG__
#pragma implementation
#endif

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_SSTREAM
  #include <sstream>
#else
  #include <HIntLib/fallback_sstream.h>
#endif

#include <HIntLib/make.h>

namespace L = HIntLib;


/**
 *   GeneratorMatrixDoesNotExist Exception
 */

void L::Make::GeneratorMatrixDoesNotExist::makeString() const
{
   std::ostringstream ss;
   ss << "GeneratorMatrix(2/Gen) #" << number << " does not exist!";
   setStringCopy (ss.str().c_str());
}


/**
 *   QRNSequenceDoesNotExist Exception
 */

void L::Make::QRNSequenceDoesNotExist::makeString() const
{
   std::ostringstream ss;
   ss << "QRNSequence/QRNNet #" << number << " does not exist!";
   setStringCopy (ss.str().c_str());
}


/**
 *   CubatureRuleDoesNotExist Exception
 */

void L::Make::CubatureRuleDoesNotExist::makeString() const
{
   std::ostringstream ss;
   ss << "CubatureRule/Embeddedrule #" << number << " does not exist!";
   setStringCopy (ss.str().c_str());
}


