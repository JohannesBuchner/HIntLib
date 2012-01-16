/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2006  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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

#include <iostream>
#include <stdexcept>

#include "test.h"
#include "test_charset.h"

using std::cout;

const char testProgramName[] = "test_charset";

void displayCharacter (unsigned long c)
{
#ifdef HINTLIB_UTF8_SELECT
   HINTLIB_UTF8_SELECT (utf8,
      char buffer [5];
      encode_utf8 (c, buffer);
      cout << "'" << buffer << "'",
      if (c <= 255)  cout << "'" << char(c) << "'";
      else cout << "???";
   )
#else
   cout << "???";
#endif

}


void test(int argc, char**)
{
   if (argc != 0)  usage("Too many arguments!");

   printHeader (cout);

   cout << "Current locale: ";
#ifdef HINTLIB_STREAMS_SUPPORT_LOCALE
   try
   {
      cout << std::locale("").name() << '\n';
   }
   catch (std::runtime_error&)
   {
      cout << "Invalid!  std::locale(\"\") failed!!!\n";
   }
#else
   cout << "Streams do not support locale\n";
#endif
   cout << "Locale of cout: "
#ifdef HINTLIB_STREAMS_SUPPORT_LOCALE
        << cout.getloc().name() << '\n'
#else
        << "Streams do not support locale\n"
#endif
        << "Character set: " << Type(HINTLIB_CHARACTER_SET) << '\n'
        << "Encoding for char: "
#if defined(HINTLIB_ENCODING_LATIN1)
      "Latin-1\n"
#elif defined(HINTLIB_ENCODING_UTF8)
      "UTF-8\n"
#elif defined(HINTLIB_ENCODING_LOCALE)
      << (utf8 ? "UTF-8" : "Latin-1") << " (auto detected from locale)\n"
#else
      "ASCII\n"
#endif
      << '\n';


   unsigned long lastC = characters[0].code - 1;
   
   for (const Data* p = characters; p != characters + numChars; ++p)
   {
      printLine (cout, *p, lastC);
   }
}

