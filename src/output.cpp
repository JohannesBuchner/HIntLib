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

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/defaults.h>

#ifdef HINTLIB_STREAMS_SUPPORT_LOCALE
#include <string>
#endif

#include <HIntLib/output.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

namespace P = HIntLib::Private;

using std::string;

/**
 *  utf8Support()
 *
 *  Determines whether UTF-8 symbols should be produced.
 */

#ifdef HINTLIB_STREAMS_SUPPORT_LOCALE
bool
P::utf8Support (const std::ostream& o)
{
   string name = o.getloc().name();
   return (name.find ("UTF-8") != string::npos
        || name.find ("utf-8") != string::npos);
}
#endif


/**
 *  Constructor
 */

P::Printer::Printer (std::ostream &_o)
   : o (_o)
#ifdef HINTLIB_ENCODING_LOCALE
   , haveUtf8Support (UNKNOWN)
#endif
{
   flags (o.flags());
   precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCALE
   imbue (o.getloc());
#endif
}

#ifdef HINTLIB_BUILD_WCHAR
P::WPrinter::WPrinter (std::wostream &_o)
   : o (_o)
{
   flags (o.flags());
   precision (o.precision());
#ifdef HINTLIB_STREAMS_SUPPORT_LOCALE
   imbue (o.getloc());
#endif
}
#endif


/**
 *  Destructor
 */

P::Printer::~Printer ()
{
   // o << str(); // string<> does not handle width properly on GNU C++ 2.95.x
   o << str().c_str();
}

#ifdef HINTLIB_BUILD_WCHAR
P::WPrinter::~WPrinter ()
{
   // o << str(); // string<> does not handle width properly on GNU C++ 2.95.x
   o << str().c_str();
}
#endif


/**
 *  utf8()
 */

#ifdef HINTLIB_ENCODING_LOCALE
bool
P::Printer::utf8 ()
{
   if (haveUtf8Support == UNKNOWN)
   {
      haveUtf8Support = utf8Support (o) ? YES : NO;
   }

   return haveUtf8Support == YES;
}
#endif


/**
 *  power ()
 */

void
P::Printer::power (unsigned i)
{
#if HINTLIB_CHARACTER_SET == 4 && (defined(HINTLIB_ENCODING_UTF8) || defined(HINTLIB_ENCODING_LOCALE))
   if (i <= 9)
   {
      HINTLIB_UTF8_SELECT (utf8(),
         switch (i)
         {
         case 0:  *this << "\xe2\x81\xb0"; return;
         case 1:  *this << "\xc2\xb9";     return;
         case 2:  *this << "\xc2\xb2";     return;
         case 3:  *this << "\xc2\xb3";     return;
         case 4:  *this << "\xe2\x81\xb4"; return;
         case 5:  *this << "\xe2\x81\xb5"; return;
         case 6:  *this << "\xe2\x81\xb6"; return;
         case 7:  *this << "\xe2\x81\xb7"; return;
         case 8:  *this << "\xe2\x81\xb8"; return;
         case 9:  *this << "\xe2\x81\xb9"; return;
         }
      ,
         switch (i)
         {
         case 1:  *this << '\xb9'; return;
         case 2:  *this << '\xb2'; return;
         case 3:  *this << '\xb3'; return;
         }
      )
   }
#elif HINTLIB_CHARACTER_SET >= 2 && defined(HINTLIB_UTF8_SELECT)
   if (i >= 1 && i <= 3)
   {
      HINTLIB_UTF8_SELECT (utf8(),
         switch (i)
         {
         case 1:  *this << "\xc2\xb9"; return;
         case 2:  *this << "\xc2\xb2"; return;
         case 3:  *this << "\xc2\xb3"; return;
         }
      ,
         switch (i)
         {
         case 1:  *this << '\xb9'; return;
         case 2:  *this << '\xb2'; return;
         case 3:  *this << '\xb3'; return;
         }
      )
   }
#endif

   *this << '^' << i;
}

#ifdef HINTLIB_BUILD_WCHAR
void
P::WPrinter::power (unsigned i)
{
#if HINTLIB_CHARACTER_SET == 4
   static const wchar_t powers [] =
   {
      L'\x2070',  // ^0
      L'\x00b9',  // ^1
      L'\x00b2',  // ^2
      L'\x00b3',  // ^3
      L'\x2074',  // ^4
      L'\x2075',  // ^5
      L'\x2076',  // ^6
      L'\x2077',  // ^7
      L'\x2078',  // ^8
      L'\x2079'   // ^9
   };

   if (i <= 9)
   {
      *this << powers [i];
      return;
   }
#elif HINTLIB_CHARACTER_SET >= 2
   static const wchar_t powers [] =
   {
      L'\x00b9',  // ^1
      L'\x00b2',  // ^2
      L'\x00b3',  // ^3
   };
   if (i >= 1 && i <= 3)
   {
      *this << powers [i - 1];
      return;
   }
#endif
   *this << L'^' << i;
}
#endif


/**
 *  subscript()
 */

void
P::Printer::subscript (unsigned i)
{
#if HINTLIB_CHARACTER_SET == 4
#if defined(HINTLIB_ENCODING_UTF8) || defined(HINTLIB_ENCODING_LOCALE)
   if (i <= 9
#ifdef HINTLIB_ENCODING_LOCALE
              && utf8()
#endif
      )
   {
      *this << "\xe2\x82" << '\x80' + char(i);
      return;
   }
#endif
#endif

   *this << '_' << i;
}

#ifdef HINTLIB_BUILD_WCHAR
void
P::WPrinter::subscript (unsigned i)
{
#if HINTLIB_CHARACTER_SET == 4
   if (i <= 9)
   {
      *this << L'\x2080' + wchar_t(i);
      return;
   }
#endif

   *this << L'_' << i;
}
#endif


/**
 *  minusSign()
 */

void
P::Printer::minusSign()
{
#if HINTLIB_CHARACTER_SET >= 3 && defined(HINTLIB_UTF8_SELECT)
   HINTLIB_UTF8_SELECT (utf8(), *this << "\xe2\x88\x92", *this << '-')
#else
   *this << '-';
#endif
}

#ifdef HINTLIB_BUILD_WCHAR
void
P::WPrinter::minusSign()
{
#if HINTLIB_CHARACTER_SET >= 3
   *this << L'\x2212';   // MINUS SIGN
#else
   *this << L'-';
#endif
}
#endif

