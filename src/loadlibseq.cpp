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

#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include <HIntLib/defaults.h>

#ifdef HAVE_SSTREAM
  #include <sstream>
#else
  #include <HIntLib/hintlib_sstream.h>
#endif

#include <HIntLib/generatormatrix.h>
#include <HIntLib/prime.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;

using std::numeric_limits;


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
   void ignoreLine()  { pos = line.length(); nl = false; }

   unsigned getUnsigned ();
   const char* getName ();

   unsigned getLineNumber () const  { return ln; }
   unsigned getPosition () const  { return pos + 1; }
   const char* getLine () const  { return line.c_str(); }

   void throwException (const char*) const;
   void throwEOF () const;

private:
   std::istream &str;
   std::string line;
   unsigned pos;
   unsigned ln;
   bool nl;
   std::string name;
};

int LineReader::get ()
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

void LineReader::putBack ()
{
   if (nl)  nl = false;
   else if (pos != 0)  --pos;
   else throw L::InternalError (__FILE__, __LINE__);
}

unsigned LineReader::getUnsigned ()
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

const char* LineReader::getName ()
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

void LineReader::throwException (const char* msg) const
{
   throw L::LineReaderException (ln, pos, line.c_str(), msg);
}
void LineReader::throwEOF () const
{
   throw L::LineReaderException (ln, pos, line.c_str(), "End of file");
}


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

   unsigned getNumber ()  { return number; }
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
   unsigned number;
   const char* s;
};

void Tokenizer::putBack ()
{
   if (nextToken != NONE || lastToken == NONE)
   {
      throw L::InternalError (__FILE__, __LINE__);
   }
   nextToken = lastToken;
}

Tokenizer::Token Tokenizer::next()
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

void Tokenizer::ignoreLine()
{
   lr.ignoreLine();
}

void Tokenizer::expectNumber ()
{
   if (next() != NUMBER)  lr.throwException ("Expected a number");
}

void Tokenizer::expectName ()
{
   if (next() != NAME)  lr.throwException ("Expected a name");
}

void Tokenizer::expectName (const char* name)
{
   if (next() != NAME || strcmp (getName(), name) != 0)
   {
      std::ostringstream ss;
      ss << "Expected name \"" << name << "\" not found";
      lr.throwException (ss.str().c_str());
   }
}


/**
 *  load Lib Seq ()
 */

namespace
{
   const char keyBeginMatrix [] = "C_MATRIX";
   const char keyEndMatrix []   = "END_C_MATRIX";
   const char keyName []        = "C_MATRIX_NAME";
   const char keyRingName []    = "C_MATRIX_RING_NAME";
   const char keyBase [    ]    = "C_MATRIX_RING_CARD";
   const char keyPrecision []   = "DIGIT_ACCURACY";
   const char keyDimension []   = "MAX_DIMENSION";
   const char keyBeginDim []    = "DIMENSION";
   const char keyEndDim []      = "END_DIMENSION";
   const char keyEndLine []     = "END_LINE";
}


L::MutableGeneratorMatrixGen<unsigned char>*
L::loadLibSeq (std::istream &str)
{
   Tokenizer t (str);

   t.expectName (keyBeginMatrix);

   unsigned base = 0, dim = 0, prec = 0;

   for (;;)
   {
      t.expectName ();
      const char* token = t.getName();

      if (strcmp (token, keyName) == 0) 
      {
         t.ignoreLine();
      }
      else if (strcmp (token, keyRingName) == 0)
      {
         t.ignoreLine();
      }
      else if (strcmp (token, keyBase) == 0)
      {
         t.expectNumber ();
         if (base != 0)  t.throwException ("C_MATRIX_RING_CARD defined twice");
         base = t.getNumber();
      }
      else if (strcmp (token, keyPrecision) == 0)
      {
         t.expectNumber ();
         if (prec != 0)  t.throwException ("DIGIT_ACCURACY defined twice");
         prec = t.getNumber();
      }
      else if (strcmp (token, keyDimension) == 0)
      {
         t.expectNumber ();
         if (dim != 0)  t.throwException ("MAX_DIMENSION defined twice");
         dim = t.getNumber();
      }
      else if (strcmp (token, keyBeginDim) == 0)
      {
         t.putBack ();
         break;
      }
      else  t.ignoreLine();
   }

   if (base == 0)  t.throwException ("C_MATRIX_RING_CARD missing");
   if (dim == 0)   t.throwException ("MAX_DIMENSION missing");
   if (prec == 0)  t.throwException ("DIGIT_ACCURACY missing");

   unsigned m = 0;

   HeapAllocatedGeneratorMatrixGen<unsigned char>* gm = 0;

   try
   {
      for (unsigned d = 0; d < dim; ++d)
      {
         t.expectName (keyBeginDim);
         t.expectNumber ();
         if (t.getNumber() != d)
         {
            t.throwException ("Another dimension expected");
         }

         for (unsigned b = 0; b < prec; ++b)
         {
            if (d == 0 && b == 0)
            {
               std::vector<char> v;

               while (t.next() == Tokenizer::NUMBER)
               {
                  unsigned x = t.getNumber();
                  if (x >= base)  t.throwException ("Entry larger than base");
                  v.push_back (x);
               }

               t.putBack();
               t.expectName (keyEndLine);

               m = v.size();
               gm = new HeapAllocatedGeneratorMatrixGen<unsigned char>
                              (base, dim, m, prec);

               for (unsigned r = 0; r < m; ++r)  gm->set (d, r, b, v[r]);
            }
            else
            {
               for (unsigned r = 0; r < m; ++r)
               {
                  t.expectNumber ();
                  unsigned x = t.getNumber();
                  if (x >= base)  t.throwException ("Entry larger than base");
                  gm->set (d, r, b, x);
               }
               t.expectName (keyEndLine);
            }
         }
         t.expectName (keyEndDim);
      }
      t.expectName (keyEndMatrix);
   }
   catch (...)
   {
      delete gm;
      throw;
   }

   return gm;
}

L::MutableGeneratorMatrixGen<unsigned char>*
L::loadLibSeq (const char* s)
{
   std::ifstream str (s);
   return loadLibSeq (str);
}


/**
 *  dump Lib Seq ()
 */


void L::GeneratorMatrix::dumpLibSeq (std::ostream &o) const
{
   o <<
      keyBeginMatrix << '\n' <<
      keyName << " HIntLibMatrix\n" <<
      keyRingName <<
      (Prime::test (getBase()) ? " Z/qZ\n" : " GF(p,n)\n") <<
      keyBase       << ' ' << getBase() << '\n' <<
      keyPrecision  << ' ' << getPrecision() << '\n' <<
      keyDimension  << ' ' << getDimension() << "\n\n";

   for (unsigned d = 0; d < getDimension(); ++d)
   {
      o << keyBeginDim << ' ' << d << '\n';

      for (unsigned b = 0; b < getPrecision(); ++b)
      {
         for (unsigned r = 0; r < getM(); ++r)
         {
            o << getDigit (d, r, b) << ' ';
         }

         o << keyEndLine << '\n';
      }

      o << keyEndDim << "\n\n";
   }

   o << keyEndMatrix << '\n';
}

void L::GeneratorMatrix::dumpLibSeq (const char* fName) const
{
   std::ofstream f (fName);
   dumpLibSeq (f);
}



/**
 *  load Binary ()
 */

namespace
{
   const char binaryMagic [] = "HIntLib GeneratorMatrix\n";
}


L::MutableGeneratorMatrixGen<unsigned char>*
L::loadBinary (std::istream &str)
{
   char s [sizeof (binaryMagic) - 1]; 
   str.read (s, sizeof (binaryMagic) - 1);

   if (! std::equal (binaryMagic, binaryMagic + sizeof(binaryMagic), s))
   {
      throw FIXME (__FILE__, __LINE__);
   }

   unsigned base = str.get();
   unsigned dim  = str.get();
   dim |= str.get() << numeric_limits<char>::digits;
   unsigned m    = str.get();
   unsigned prec = str.get();

   if (! str)  throw FIXME (__FILE__, __LINE__);

   HeapAllocatedGeneratorMatrixGen<unsigned char>* gm =
      new HeapAllocatedGeneratorMatrixGen<unsigned char> (base, dim, m, prec);

   char check = 0;
   
   for (unsigned d = 0; d < gm->getDimension(); ++d)
   {
      for (unsigned b = 0; b < gm->getPrecision(); ++b)
      {
         for (unsigned r = 0; r < gm->getM(); ++r)
         {
            char digit = str.get ();
            gm->set (d, r, b, digit);
            check ^= digit;
         }
      }
   }

   if (! str || check != str.get() || ! str || str.get() != EOF)
   {
      throw FIXME (__FILE__, __LINE__);
   }

   return gm;
}

L::MutableGeneratorMatrixGen<unsigned char>*
L::loadBinary (const char* s)
{
   std::ifstream str (s);
   return loadBinary (str);
}


/**
 *  dump Binary ()
 */


void L::GeneratorMatrix::dumpBinary (std::ostream &o) const
{
   if (   getBase ()      > numeric_limits<unsigned char>::max()
       || getDimension()  > (1u << 2 * numeric_limits<char>::digits)
       || getM ()         > numeric_limits<unsigned char>::max()
       || getPrecision () > numeric_limits<unsigned char>::max())
   {
      throw FIXME (__FILE__, __LINE__);
   }

   o << binaryMagic;
   o.put (getBase());
   o.put (getDimension() & char (-1));
   o.put (getDimension() >> numeric_limits<char>::digits);
   o.put (getM());
   o.put (getPrecision());

   char check = 0;

   for (unsigned d = 0; d < getDimension(); ++d)
   {
      for (unsigned b = 0; b < getPrecision(); ++b)
      {
         for (unsigned r = 0; r < getM(); ++r)
         {
            char digit = getDigit (d, r, b);
            o.put (digit);
            check ^= digit;
         }
      }
   }

   o.put (check);
}

void L::GeneratorMatrix::dumpBinary (const char* fName) const
{
   std::ofstream f (fName);
   dumpBinary (f);
}


/**
 *  load Niederreiter Xing ()
 */

L::MutableGeneratorMatrixGen<unsigned char>*
L::loadNiederreiterXing (unsigned dim)
{
   std::ostringstream ss;
   ss << DATADIR << "/fkmat"
      << std::setw(2) << std::setfill ('0') << dim << ".bin";
   try
   {
      MutableGeneratorMatrixGen<unsigned char>* m =
         loadBinary (ss.str().c_str());

      return m;
   }
   catch (...)
   {
      throw InvalidDimension (dim);
   }
}

