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

/**
 *  test.cpp
 *  test.h
 *
 *  Provides a common main() routine and some utility functions shared by all
 *  test programs.
 */

#ifndef HINTLIB_TEST_H
#define HINTLIB_TEST_H

#include <iosfwd>

#include <HIntLib/defaults.h>

namespace L = HIntLib;

extern int verbose;

/*
 *  Define 
 */

#ifdef HINTLIB_ENCODING_LOCALE
extern bool utf8;
# if HINTLIB_CHARACTER_SET >= 4
#  define UnicodeAscii(u,a) (utf8 ? u : a)
# endif
# if HINTLIB_CHARACTER_SET >= 3
#  define Wgl4Ascii(u,a) (utf8 ? u : a)
# endif
# if HINTLIB_CHARACTER_SET >= 2
#  define Latin1Ascii(u,l,a) (utf8 ? u : l)
# endif
#endif

#ifdef HINTLIB_ENCODING_LATIN1
#define utf8 false
# if HINTLIB_CHARACTER_SET >= 2
#  define Latin1Ascii(u,l,a) (l)
# endif
#endif

#ifdef HINTLIB_ENCODING_UTF8
#define utf8 true
# if HINTLIB_CHARACTER_SET >= 4
#  define UnicodeAscii(u,a) (u)
# endif
# if HINTLIB_CHARACTER_SET >= 3
#  define Wgl4Ascii(u,a) (u)
# endif
# if HINTLIB_CHARACTER_SET >= 2
#  define Latin1Ascii(u,l,a) (u)
# endif
#endif

#ifndef UnicodeAscii
# define UnicodeAscii(u,a) (a)
#endif

#ifndef Wgl4Ascii
# define Wgl4Ascii(u,a) (a)
#endif

#ifndef Latin1Ascii
# define Latin1Ascii(u,l,a) (a)
#endif


extern const char options[];
extern const char option_msg[];
extern const char testProgramParameters[];
extern const char testProgramUsage[];
extern const char testProgramName[];
extern const int  testProgramCopyright;

void error();
void error(const char*);
void usage(const char*);
int parseInt (const char*);
void parseRange (const char*, int, int&, int&);
void printHeader (std::ostream&);
void doubleQuote (std::ostream&, const char*);
#ifdef HINTLIB_BUILD_WCHAR
void printHeader (std::wostream&);
void doubleQuote (std::wostream&, const char*);
void doubleQuote (std::wostream&, const wchar_t*);
#endif

void test (int argc, char** argv);
bool opt (int, const char*);

#define SILENT if (verbose <= 0)
#define NORMAL if (verbose >= 1)
#define DEB1   if (verbose >= 2)
#define DEB2   if (verbose >= 3)
#define DEB3   if (verbose >= 4)
#define DEB4   if (verbose >= 5)

#endif

