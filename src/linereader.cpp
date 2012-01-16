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

#ifdef HINTLIB_HAVE_SSTREAM
  #include <sstream>
#else
  #include <HIntLib/fallback_sstream.h>
#endif

#include <HIntLib/linereader.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

#include <HIntLib/exception.h>


namespace L = HIntLib;


/**
 *  Line Reader
 */

int
L::LineReader::get ()
{
   if (nl)
   {
      nl = false;

      std::getline (str, line);
      ++ln;
      pos = 0;

      if (!str)  return EOF;
   }

   if (pos == line.length())
   {
      nl = true;
      return '\n';
   }

   return line[pos++];
}

void
L::LineReader::putBack ()
{
   if (nl)  nl = false;
   else if (pos != 0)  --pos;
   else throw InternalError (__FILE__, __LINE__);
}

unsigned
L::LineReader::getUnsigned ()
{
   int c = get();
   if (c == EOF)  throwEOF ();
   if (! isdigit(c))  throwException ("Number expected");
   putBack();
   unsigned x = 0;

   while ((c = get()) != EOF && isdigit (c))  x = 10 * x + (c - '0');

   if (c != EOF)  putBack();

   return x;
}

const char*
L::LineReader::getName ()
{
   int c = get();
   if (c == EOF)  throwEOF ();
   if (! isalpha(c))  throwException ("Invalid first character for name");
   putBack();

   name = "";

   while ((c = get()) != EOF && (isalnum (c) || c == '_'))  name += c;

   if (c != EOF)  putBack();

   return name.c_str();
}

void
L::LineReader::throwException (const char* msg) const
{
   throw L::LineReaderException (ln, pos, line.c_str(), msg);
}
void
L::LineReader::throwEOF () const
{
   throw L::LineReaderException (ln, pos, line.c_str(), "End of file");
}


/**
 *  Tokenizer
 */

void
L::Tokenizer::putBack ()
{
   if (nextToken != NONE || lastToken == NONE)
   {
      throw L::InternalError (__FILE__, __LINE__);
   }
   nextToken = lastToken;
}

L::Tokenizer::Token
L::Tokenizer::next()
{
   if (nextToken != NONE)
   {
      lastToken = nextToken;
      nextToken = NONE;
      return lastToken;
   }

   int c;

   for (;;)
   {
      // eat white space

      while ((c = lr.get()) != EOF)
      {
         if (!isspace(c))  break;
      }

      if (c == EOF)  return lastToken = END;

      if (c != '#')  break;

      ignoreLine();
   }

   if (isdigit (c))
   {
      lr.putBack ();
      number = lr.getUnsigned();
      return lastToken = NUMBER;
   }

   lr.putBack ();
   s = lr.getName();
   return lastToken = NAME;
}

void
L::Tokenizer::ignoreLine()
{
   lr.ignoreLine();
}

void
L::Tokenizer::expectNumber ()
{
   if (next() != NUMBER)  lr.throwException ("Expected a number");
}

void
L::Tokenizer::expectName ()
{
   if (next() != NAME)  lr.throwException ("Expected a name");
}

void
L::Tokenizer::expectName (const char* name)
{
   if (next() != NAME || strcmp (getName(), name) != 0)
   {
      std::ostringstream ss;
      ss << "Expected name \"" << name << "\" not found";
      lr.throwException (ss.str().c_str());
   }
}


