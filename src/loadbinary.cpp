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

#if defined(HINTLIB_HAVE_OSTREAM) && defined(HINTLIB_HAVE_ISTREAM)
#  include <ostream>
#  include <istream>
#else
#  include <iostream>
#endif

#include <fstream>
#include <iomanip>
#include <algorithm>

#ifdef HINTLIB_HAVE_SSTREAM
#  include <sstream>
#else
#  include <HIntLib/fallback_sstream.h>
#endif

#include <HIntLib/generatormatrixgen.h>
#include <HIntLib/exception.h>


namespace L = HIntLib;

using std::numeric_limits;


/**
 *  load Binary ()
 */

namespace
{
   const char binaryMagic [] = "HIntLib GeneratorMatrix\n";
}


L::GeneratorMatrixGen<unsigned char>* L::loadBinary (std::istream &str)
{
   char s [sizeof (binaryMagic) - 1];
   str.read (s, sizeof (binaryMagic) - 1);
   if (! str)  throw FIXME (__FILE__, __LINE__);

   if (! std::equal (s, s + sizeof (s), binaryMagic))
   {
      throw FIXME (__FILE__, __LINE__);
   }

   int base = str.get();
   int dim  = str.get();
   dim |= str.get() << numeric_limits<char>::digits;
   int m    = str.get();
   int prec = str.get();

   if (! str)  throw FIXME (__FILE__, __LINE__);

   GeneratorMatrixGen<unsigned char>* gm =
      new GeneratorMatrixGen<unsigned char>(base, dim, m, prec);

   char check = 0;

   for (int d = 0; d < gm->getDimension(); ++d)
   {
      for (int b = 0; b < gm->getPrec(); ++b)
      {
         for (int r = 0; r < gm->getM(); ++r)
         {
            char digit = str.get ();
            gm->setd (d, r, b, digit);
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


/**
 *  Binary Export ()
 */


void L::GeneratorMatrix::binaryExport (std::ostream &o) const
{
   if (   unsigned(getBase()) > numeric_limits<unsigned char>::max()
       || getDimension()      > (1 << 2 * numeric_limits<char>::digits)
       || unsigned(getM ())   > numeric_limits<unsigned char>::max()
       || unsigned(getPrec()) > numeric_limits<unsigned char>::max())
   {
      throw FIXME (__FILE__, __LINE__);
   }

   o << binaryMagic;
   o.put (getBase());
   o.put (getDimension() & ~char(0));
   o.put (getDimension() >> numeric_limits<char>::digits);
   o.put (getM());
   o.put (getPrec());

   char check = 0;

   for (int d = 0; d < getDimension(); ++d)
   {
      for (int b = 0; b < getPrec(); ++b)
      {
         for (int r = 0; r < getM(); ++r)
         {
            char digit = getDigit (d, r, b);
            o.put (digit);
            check ^= digit;
         }
      }
   }

   o.put (check);
}


/**
 *  load Niederreiter Xing ()
 */

namespace
{
   const char* nxFileNames [] =
   {
      HINTLIB_DATADIR "/fkmat",
      HINTLIB_SRCDATADIR "/fkmat"
   };
}

L::GeneratorMatrixGen<unsigned char>* L::loadNiederreiterXing (int dim)
{
   for (int i = 0; i < 2; ++i)
   {
      std::ostringstream ss;
      ss << nxFileNames [i]
         << std::setw(2) << std::setfill ('0') << dim << ".bin";

      try
      {
         std::ifstream ifs (ss.str().c_str(), std::ios::binary);
         if (! ifs)  continue; 

         return loadBinary (ifs);
      }
      catch (...)  { continue; }
   }

   throw InvalidDimension (dim);
}

