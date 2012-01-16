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

#ifndef HINTLIB_TEST_CHARSET_H
#define HINTLIB_TEST_CHARSET_H

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_SSTREAM
  #include <sstream>
#else
  #include <HIntLib/fallback_sstream.h>
#endif

enum Type { ASCII = 1, LATIN1 = 2, WGL4 = 3, UNICODE = 4 };

std::ostream& operator<< (std::ostream&, Type);

extern const struct Data
{
   unsigned long code;
   Type          type;
   const char*   name;
}
characters [];

extern const unsigned numChars;

void displayCharacter (unsigned long c);
void encode_utf8 (unsigned long, char*);
std::string beginOfLine (const Data& data);


template<class OSTREAM>
void printLine(OSTREAM& COUT, const Data& data, unsigned long lastC)
{
   unsigned long c = data.code;
   if (c > lastC + 1)  COUT << '\n';

   std::string line = beginOfLine (data);

   COUT << line.c_str();

   if (data.type <= HINTLIB_CHARACTER_SET)
   {
      displayCharacter (c);
   }
   else
   {
      COUT << "---";
   }
   
   COUT << "  " << data.name << '\n';

   lastC = c;
}


#endif

