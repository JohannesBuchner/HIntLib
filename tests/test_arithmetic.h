/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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

#ifndef HINTLIB_TEST_ARITHMETIC_H
#define HINTLIB_TEST_ARITHMETIC_H 1

#include <iostream>
#include <string>

#include "test.h"

#ifdef HINTLIB_HAVE_SSTREAM
  #include <sstream>
#else
  #include <HIntLib/fallback_sstream.h>
#endif

#ifdef HINTLIB_BUILD_WCHAR
#undef USE_WCHAR
#else
#undef USE_WCHAR
#endif

#ifdef USE_WCHAR

typedef wchar_t CHAR;
#define COUT std::wcout
typedef std::wstring STRING;
typedef std::basic_stringstream<CHAR> STRINGSTREAM;

#else

typedef char CHAR;
#define COUT std::cout
typedef std::string STRING;
typedef std::stringstream STRINGSTREAM;
#endif

// Userdefined constants

extern unsigned SIZE;
extern bool FLUSH;
extern unsigned W;

// counters

extern unsigned nilpotentsCounter;
extern unsigned unitsCounter;
extern unsigned primitivesCounter;
extern unsigned lastNorm;
extern int lastDegree;

// function prototypes

bool performTest (
      const STRING& cat, const STRING& name, const std::string& type);
void printNumberOrInf (unsigned n);
void printInfinity ();
void checkCounter (unsigned expected, unsigned actual, bool all, const char* s);
void checkIndex (unsigned size, unsigned index, const char* s);
void checkInfiniteInFinite (unsigned size, unsigned n, const char* s);

#endif

