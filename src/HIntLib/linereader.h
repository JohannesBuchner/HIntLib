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

#ifndef HINTLIB_LINE_READER_H
#define HINTLIB_LINE_READER_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma interface
#endif

#include <iosfwd>
#include <string>


namespace HIntLib
{

/**
 *  Line Reader
 */

class LineReader
{
public:
   LineReader (std::istream &_str)
      : str (_str), line (""), pos (0), ln (0), nl (true) {}

   int get();
   void putBack ();
   void putBackLine() { nl = false; pos = 0; }
   void ignoreLine()  { nl = false; pos = line.length(); }

   unsigned getUnsigned ();
   const char* getName ();

   int getLineNumber () const  { return ln; }
   int getPosition () const  { return pos + 1; }
   const char* getLine () const  { return line.c_str(); }

   void throwException (const char*) const;
   void throwEOF () const;

private:
   std::istream &str;
   std::string line;
   unsigned pos;
   int ln;
   bool nl;
   std::string name;
};


/**
 *  Tokenizer
 */

class Tokenizer
{
public:

   Tokenizer (std::istream &_str)
      : lr(_str), nextToken (NONE), lastToken (NONE) {}

   enum Token { END, NUMBER, NAME, NONE };
   Token next ();
   void putBack ();

   int getNumber ()  { return number; }
   const char* getName()  { return s; }
   void expectName (const char*);
   void expectName ();
   void expectNumber ();
   void ignoreLine ();
   void throwException (const char* msg) const  { lr.throwException (msg); }

private:
   LineReader lr;
   Token nextToken;
   Token lastToken;
   int number;
   const char* s;
};

} // namespace HIntLib

#endif

