/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration
 *
 *  Copyright (C) 2002,03,04,05  Rudolf Schürer <rudolf.schuerer@sbg.ac.at>
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


#include <iostream>
#include <iomanip>

// This is not really a library object. However, since library files are
// linked static, we do not require DLL-export names.
#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/old_bitop.h>

using std::cout;
using std::setw;

using namespace HIntLib::NormalBitOp;

template<class F>
inline
void create_data (F* f)
{
   cout << "{\n";

   for (unsigned i = 0; i < 256; ++i)
   {
      if (i % 16 == 0)  cout << "  ";

      unsigned char c = i;

      cout << setw (2) << f (c) << ",";

      if (i % 16 == 15)  cout << '\n';
   }

   cout << "};\n\n";
}

int main()
{
   cout << "/*******************************************************/\n"
           "/***   This file is program-generated!               ***/\n"
           "/***                                                 ***/\n"
           "/***   Do not change!!!                              ***/\n"
           "/***                                                 ***/\n"
           "/***   Update " __FILE__ " to update this file.  ***/\n"
           "/*******************************************************/\n"
           "\n"
           "#ifdef __GNUG__\n"
           "#pragma implementation\n"
           "#endif\n"
           "\n"
           "#include <HIntLib/bitop.h>\n"
           "\n"
           "namespace L = HIntLib;\n"
           "\n"
           "const int L::Private::ms1_data [] =\n";

   typedef int F1 (unsigned char);

   create_data (static_cast<F1*> (ms1));

   cout << "const int L::Private::ls0_data [] =\n";

   typedef unsigned int F2 (unsigned char);

   create_data (static_cast<F2*> (ls0));

   return 0;
}

