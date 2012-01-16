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

#ifndef HINTLIB_OUTPUT_H
#define HINTLIB_OUTPUT_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#ifdef HINTLIB_HAVE_SSTREAM
  #include <sstream>
#else
  #include <HIntLib/fallback_sstream.h>
#endif

namespace HIntLib
{
namespace Private
{

#ifdef HINTLIB_STREAMS_SUPPORT_LOCALE
bool utf8Support (const std::ostream&);
#endif

template<typename S>
class PrinterSelector
{};

class Printer : public std::ostringstream
{
public:
   Printer (std::ostream &);
   ~Printer ();

#ifdef HINTLIB_ENCODING_LOCALE
   bool utf8 ();
#endif
   void power (unsigned i);
   void subscript (unsigned i);
   void minusSign();

private:
   std::ostream& o;
#ifdef HINTLIB_ENCODING_LOCALE
   enum { UNKNOWN, YES, NO } haveUtf8Support;
#endif
};

template<>
class PrinterSelector<char>
{
public:
   typedef Printer printer;
};

#ifdef HINTLIB_BUILD_WCHAR
class WPrinter : public std::wostringstream
{
public:
   WPrinter (std::wostream &);
   ~WPrinter ();

   void power (unsigned i);
   void subscript (unsigned i);
   void minusSign();

private:
   std::wostream& o;
};

template<>
class PrinterSelector<wchar_t>
{
public:
   typedef WPrinter printer;
};
#endif

   
}  // namespace Private
}  // namespace HIntLib

#endif
