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

#include "test.h"
#include "test_charset.h"

const char testProgramName[] = "test_charsetw";

#ifdef HINTLIB_BUILD_WCHAR
void displayCharacter (unsigned long c)
{
   std::wcout << L'\'' << wchar_t (c) << L'\'';
}


void test(int argc, char**)
{
   if (argc != 0)  usage("Too many arguments!");

   printHeader (std::wcout);

   unsigned long lastC = characters[0].code - 1;
   
   for (const Data* p = characters; p != characters + numChars; ++p)
   {
      printLine (std::wcout, *p, lastC);
   }
}
#else
void test(int argc, char**)
{
   usage ("Type wchar_t is not supported!");
}
#endif

