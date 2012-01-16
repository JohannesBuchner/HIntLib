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

#if defined(HINTLIB_HAVE_OSTREAM) && defined(HINTLIB_HAVE_ISTREAM)
#  include <ostream>
#  include <istream>
#else
#  include <iostream>
#endif

#ifdef HINTLIB_HAVE_CSTRING
#  include <cstring>
#  define HINTLIB_SSN std::
#else
#  include <string.h>
#  define HINTLIB_SSN
#endif

#include <iomanip>

#include <HIntLib/linereader.h>
#include <HIntLib/generatormatrixgen.h>
#include <HIntLib/prime.h>


namespace L = HIntLib;


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


L::GeneratorMatrixGen<unsigned char>*
L::loadLibSeq (std::istream &str)
{
   Tokenizer t (str);

   t.expectName (keyBeginMatrix);

   int base = 0, dim = 0, prec = 0;

   for (;;)
   {
      t.expectName ();
      const char* token = t.getName();

      if (HINTLIB_SSN strcmp (token, keyName) == 0)
      {
         t.ignoreLine();
      }
      else if (HINTLIB_SSN strcmp (token, keyRingName) == 0)
      {
         t.ignoreLine();
      }
      else if (HINTLIB_SSN strcmp (token, keyBase) == 0)
      {
         t.expectNumber ();
         if (base != 0)  t.throwException ("C_MATRIX_RING_CARD defined twice");
         base = t.getNumber();
      }
      else if (HINTLIB_SSN strcmp (token, keyPrecision) == 0)
      {
         t.expectNumber ();
         if (prec != 0)  t.throwException ("DIGIT_ACCURACY defined twice");
         prec = t.getNumber();
      }
      else if (HINTLIB_SSN strcmp (token, keyDimension) == 0)
      {
         t.expectNumber ();
         if (dim != 0)  t.throwException ("MAX_DIMENSION defined twice");
         dim = t.getNumber();
      }
      else if (HINTLIB_SSN strcmp (token, keyBeginDim) == 0)
      {
         t.putBack ();
         break;
      }
      else  t.ignoreLine();
   }

   if (base == 0)  t.throwException ("C_MATRIX_RING_CARD missing");
   if (dim == 0)   t.throwException ("MAX_DIMENSION missing");
   if (prec == 0)  t.throwException ("DIGIT_ACCURACY missing");

   int m = 0;

   GeneratorMatrixGen<unsigned char>* gm = 0;

   try
   {
      for (int d = 0; d < dim; ++d)
      {
         t.expectName (keyBeginDim);
         t.expectNumber ();
         if (t.getNumber() != d)
         {
            t.throwException ("Another dimension expected");
         }

         for (int b = 0; b < prec; ++b)
         {
            if (d == 0 && b == 0)
            {
               std::vector<char> v;

               while (t.next() == Tokenizer::NUMBER)
               {
                  int x = t.getNumber();
                  if (x >= base)  t.throwException ("Entry larger than base");
                  v.push_back (x);
               }

               t.putBack();
               t.expectName (keyEndLine);

               m = v.size();
               gm = new GeneratorMatrixGen<unsigned char> (base, dim, m, prec);

               for (int r = 0; r < m; ++r)  gm->setd (d, r, b, v[r]);
            }
            else
            {
               for (int r = 0; r < m; ++r)
               {
                  t.expectNumber ();
                  int x = t.getNumber();
                  if (x >= base)  t.throwException ("Entry larger than base");
                  gm->setd (d, r, b, x);
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


/**
 *  Lib Seq Export ()
 */


void L::GeneratorMatrix::libSeqExport (std::ostream &o) const
{
   o <<
      keyBeginMatrix << '\n' <<
      keyName << " HIntLibMatrix\n" <<
      keyRingName <<
      (Prime::test (getBase()) ? " Z/qZ\n" : " GF(p,n)\n") <<
      keyBase       << ' ' << getBase() << '\n' <<
      keyPrecision  << ' ' << getPrec() << '\n' <<
      keyDimension  << ' ' << getDimension() << "\n\n";

   for (int d = 0; d < getDimension(); ++d)
   {
      o << keyBeginDim << ' ' << d << '\n';

      for (int b = 0; b < getPrec(); ++b)
      {
         for (int r = 0; r < getM(); ++r)
         {
            o << getDigit (d, r, b) << ' ';
         }

         o << keyEndLine << '\n';
      }

      o << keyEndDim << "\n\n";
   }

   o << keyEndMatrix << '\n';
}


