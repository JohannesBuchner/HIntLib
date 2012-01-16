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

#include <iostream>
#include <iomanip>

#include <HIntLib/galoisfield.h>
#include <HIntLib/prime.h>

using std::cout;
using std::setw;

namespace L = HIntLib;

void indent (unsigned in)
{
   for (unsigned i = 0; i < in; ++i)  cout << ' ';
}

void CArrayDump (const L::GaloisField<unsigned char> f, unsigned in)
{
   typedef L::GaloisField<unsigned char>::type T;

   indent (in);
   cout << "// Multiplication\n";

   for (unsigned i = 0; i < f.size(); ++i)
   {
      indent (in);
      T a = f.element (i);

      for (unsigned j = 0; j < f.size(); ++j)
      {
         cout << setw(2)
              << unsigned (f.index (f.mul (a, f.element (j)))) << ',';
      }
      cout << '\n';
   }

   indent (in);
   cout << "// Reciprocals\n";
   indent (in);
   cout << " 0,";

   for (unsigned i = 1; i < f.size(); ++i)
   {
      cout << setw(2) << unsigned (f.index (f.recip (f.element (i)))) << ',';
   }
   cout << '\n';

   indent (in);
   cout << "// Squares\n";
   indent (in);

   for (unsigned i = 0; i < f.size(); ++i)
   {
      cout << setw(2) << unsigned (f.index (f.sqr (f.element (i)))) << ',';
   }
   cout << '\n';

   indent (in);
   cout << "// Multiplicative order\n";
   indent (in);
   cout << " 0,";

   for (unsigned i = 1; i < f.size(); ++i)
   {
      cout << setw(2) << f.order (f.element (i)) << ',';
   }
   cout << '\n';

   indent (in);
   cout << "// Frobenius\n";
   indent (in);

   for (unsigned i = 0; i < f.size(); ++i)
   {
      cout << setw(2)
           << unsigned (f.index (f.frobenius (f.element (i)))) << ',';
   }
   cout << '\n';

   indent (in);
   cout << "// Inverse Frobenius\n";
   indent (in);

   for (unsigned i = 0; i < f.size(); ++i)
   {
      cout << setw(2)
           << unsigned (f.index (f.invFrobenius (f.element (i)))) << ',';
   }
   cout << '\n';

   indent (in);
   cout << "// Addition\n";

   for (unsigned i = 0; i < f.size(); ++i)
   {
      indent (in);
      T a = f.element (i);

      for (unsigned j = 0; j < f.size(); ++j)
      {
         cout << setw(2)
              << unsigned (f.index (f.add (a, f.element (j)))) << ',';
      }
      cout << '\n';
   }

   indent (in);
   cout << "// Negatives\n";
   indent (in);

   for (unsigned i = 0; i < f.size(); ++i)
   {
      cout << setw(2) << unsigned (f.index (f.neg (f.element (i)))) << ',';
   }
   cout << "\n";

   indent (in);
   cout << "// Doubles\n";
   indent (in);

   for (unsigned i = 0; i < f.size(); ++i)
   {
      cout << setw(2) << unsigned (f.index (f.dbl (f.element (i)))) << ',';
   }
   cout << "\n";
}

unsigned arraySize (unsigned size)
{
   return 2 * size * size + 7 * size;
}

int main()
{
   cout <<
      "/************************************************************/\n"
      "/***   This file is program-generated!                    ***/\n"
      "/***                                                      ***/\n"
      "/***   Do not change!!!                                   ***/\n"
      "/***                                                      ***/\n"
      "/***   Update " __FILE__ " to update this file. ***/\n"
      "/************************************************************/\n"
      "\n"
      "#include <HIntLib/defaults.h>\n"
      "\n"
      "namespace\n"
      "{\n"
      "\n"
      "template<int size>\n"
      "struct Data\n"
      "{\n"
      "   unsigned count;\n"
      "   unsigned characteristic;\n"
      "   unsigned degree;\n"
      "   unsigned char data [size];\n"
      "};\n"
      "\n";

   for (unsigned size = 2; size <= HINTLIB_PRECALCULATED_FIELD_MAX_SIZE; ++size)
   {
      unsigned prime, power;
      if (L::Prime::isPrimePower(size, prime, power))
      {
         L::GaloisField<unsigned char> gf (prime, power);
         cout <<
            "const Data<" << arraySize(size) << ">"
               " precalculatedLookupField" << gf.size() << " = \n"
            "{ 0, " << gf.characteristic() << ", "
                    << gf.extensionDegree() << ",\n"
            "   {\n";
         CArrayDump (gf, 6);
         cout <<
            "   }\n"
            "};\n\n";
      }
   }

   cout <<
      "}  // anonymous namespace\n"
      "\n"
      "namespace HIntLib\n"
      "{\n"
      "   namespace Private\n"
      "   {\n"
      "      const unsigned* lookupFields"
               " [HINTLIB_PRECALCULATED_FIELD_MAX_SIZE + 1] =\n"
      "      {\n";

   for (unsigned size = 0; size <= HINTLIB_PRECALCULATED_FIELD_MAX_SIZE; ++size)
   {
      if (L::Prime::isPrimePower(size))
      {
         cout << "         reinterpret_cast<const unsigned*>"
                 " (&precalculatedLookupField" << size << "),\n";
      }
      else
      {
         cout << "         0,\n";
      }
   }

   cout <<
      "      };\n"
      "   }  // namespace Private\n"
      "}  // namespace HIntLib\n\n";

   return 0;
}

