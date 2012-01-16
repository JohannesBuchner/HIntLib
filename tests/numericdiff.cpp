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

/**
 *  numeric diff
 *
 *  Compares the numbers listed in two files.
 */

#include <fstream>
#include <iostream>
#include <HIntLib/hlmath.h>

using std::cerr;
using std::endl;

int main (int argc, char** argv)
{
   if (argc != 3)
   {
      cerr << "Usage: numericdiff file1 file2\n\n";

      return 1;
   }

   std::ifstream f1 (argv[1]);
   if (! f1)
   {
      cerr << "Cannot open '" << argv[1] << "'!\n\n";
      return 2;
   }
   std::ifstream f2 (argv[2]);
   if (! f2)
   {
      cerr << "Cannot open '" << argv[2] << "'!\n\n";
      return 2;
   }

   for (unsigned n = 1; ; ++n)
   {
      double x1, x2;

      f1 >> x1;
      f2 >> x2;

      if (f1.eof() && f2.eof())  return 0;

      if (f1.eof())
      {
         cerr << "File 1 is shorter!\n\n";
         return 3;
      }

      if (f2.eof())
      {
         cerr << "File 2 is shorter!\n\n";
         return 3;
      }

      if (HIntLib::abs (x1 - x2) > 0.000001)
      {
         cerr << "Files differ! Position: " << n << ", "
              << x1 << " != " << x2 << endl;
         return 4;
      }
   }
}


