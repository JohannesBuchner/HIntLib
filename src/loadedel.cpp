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

#include <vector>

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_ISTREAM
  #include <istream>
#else
  #include <iostream>
#endif

#include <HIntLib/linereader.h>
#include <HIntLib/hlalgorithm.h>
#include <HIntLib/generatormatrixgen.h>


namespace L = HIntLib;


/**
 *  load Edel()
 */

L::GeneratorMatrixGen<unsigned char>*
L::loadEdel (std::istream &str, unsigned base)
{
   LineReader lr (str);

   // Ignore leading space or newline

   {
      int c;

      do
      {
         c = lr.get();
      }
      while (c == '\n' || c == ' ');

      lr.putBack();
   }

   // Read the first column vector in order to determin  prec

   unsigned prec = 0;
   while (isdigit (lr.get()))  ++prec;

   if (prec == 0)  lr.throwException ("expecting digits");
   
   lr.putBack();

   // Read first line in order to determin  s

   unsigned dim = 1;

   while (lr.get() == '*')
   {
      for (unsigned b = 0; b < prec; ++b)
      {
         if (! isdigit (lr.get()))  lr.throwException ("digit expected");
      }

      ++dim;
   }

   lr.putBackLine();

   // read all lines
   
   unsigned size = prec * dim;
   std::vector<unsigned char*> lines;
   unsigned char* line = 0;

   try
   {
      for (;;)
      {
         // go to next line

         {
            int c = lr.get();
            if (c == EOF)  break;
            if (c == '\n' || c == ' ')  continue;
            lr.putBack();
         }

         // read line into buffer

         line = new unsigned char [size];
         unsigned char* p = line;

         for (unsigned d = 0; d < dim; ++d)
         {
            for (unsigned b = 0; b < prec; ++b)
            {
               int c = lr.get();
               unsigned digit = c - '0';
               if (digit >= base)  lr.throwException ("digit larger than base");
               *p++ = digit;
            }

            if (d + 1 < dim && lr.get() != '*')
            {
               lr.throwException ("'*' expected");
            }
         }

         // store buffer

         lr.ignoreLine();
         lines.push_back (line); line = 0;
      }
   }
   catch (...)
   {
      delete[] line;
      purge (lines.begin(), lines.end());

      throw;
   }

   // create matrix of appropriate size

   unsigned m = lines.size();

   GeneratorMatrixGen<unsigned char>* gm =
      new GeneratorMatrixGen<unsigned char> (base, dim, m, prec);

   // copy content into matrix

   for (unsigned r = 0; r < m; ++r)
   {
      unsigned char* currentLine = lines[r];

      for (unsigned d = 0; d < dim; ++d)
      {
         for (unsigned b = 0; b < prec; ++b)
         {
            gm->setd (d, r, b, *currentLine++);
         }
      }

      delete[] lines[r];
   }

   return gm;
}

