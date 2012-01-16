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

#include <memory>
#include <iostream>
#include <fstream>
#include <exception>

#include <HIntLib/generatormatrixgen.h>

using std::cerr;
using std::endl;

int main (int argc, char** argv)
{
   if (argc != 3)
   {
      cerr << "Usage: convert libSeqMatrix BinaryMatirx\n\n";
      exit (1);
   }

   try
   {
      // Read matrix from input file

      std::ifstream ifs (argv[1]);

      std::auto_ptr<HIntLib::GeneratorMatrix> m (HIntLib::loadLibSeq (ifs));

      if (! ifs)
      {
         cerr << "Error reading input file '" << argv[1] << "'!\n\n";
         return 1;
      }

      // Write matrix to output file

      std::ofstream ofs (argv[2]);

      m->binaryExport (ofs);

      if (! ofs)
      {
         cerr << "Error writing output file '" << argv[2] << "'!\n\n";
         return 1;
      }
   }
   catch (std::exception &e)
   {
      cerr << "Exception: " << e.what() << endl;
   }
   catch (...)
   {
      cerr << "Unknown exception!" << endl;
   }

   return 0;
}

