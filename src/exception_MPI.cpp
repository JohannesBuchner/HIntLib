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

#include <HIntLib/hlmpi.h>

#include <HIntLib/exception_MPI.h>

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


void L::TooFewNodes::makeString() const
{
   setStringCopy ("Too few processing nodes for this Integrator!");
} 

void L::MPIError::makeString() const
{
   ostringstream ss;
   ss << "MPI error number " << error << " encountered!";
   setString (ms(ss));
}

#undef HINTLIB_SSN



