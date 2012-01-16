/**
 *  numeric diff
 *
 *  Compares the numbers listed in two files.
 */

#include <fstream>
#include <iostream>
#include <math.h>

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
      cerr << "Can not open '" << argv[1] << "'!\n\n";
      return 2;
   }
   std::ifstream f2 (argv[2]);
   if (! f2)
   {
      cerr << "Can not open '" << argv[2] << "'!\n\n";
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

      if (fabs (x1 - x2) > 0.000001)
      {
         cerr << "Files differ! Position: " << n << ", "
              << x1 << " != " << x2 << endl;
         return 4;
      }
   }
}


