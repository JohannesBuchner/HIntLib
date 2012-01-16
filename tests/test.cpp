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

/**
 *  test.cpp
 *  test.h
 *
 *  Provides a common main() routine and some utility functions shared by all
 *  test programs.
 */

#include <iostream>
#include <unistd.h>
#include <string>

#include <HIntLib/exception.h>

#include "test.h"

using std::cerr;
using std::cout;
using std::endl;

int verbose = 1;
int status  = 0;

/**
 *  error()
 *
 *  Increments the error counter.
 *  If a total of 10 errors is reached, the program aborts.
 */

void error()
{
  if (++status == 10)
  {
     cout << "\nToo many errors encountered! Terminating test!\n" << endl;

     exit (status);
  }
}


/**
 *  Displays an error message before calling error().
 */

void error(const char* s)
{
  cout << "\nError: " << s << "!!!\n" << endl;
  error();
}

const char option_msg [] =
  "  -s     Silent mode (verbosity level 0)\n"
  "  -q     Quiet. Same as -s\n"
  "  -v     Increases verbosity level by one (can be applied multiple times)\n"
  "  -h     Display usage information\n"
  "  -?     Same as -h\n";


int main (int argc, char** argv)
{
   opterr = 0;  // do not print error messages

   std::string allOptions = "sqvh?";
   allOptions += options;

   for(;;)
   {
      int c = getopt (argc, argv, allOptions.c_str());

      if (c == -1)  break;

      switch (c)
      {
      case 'v':  ++ verbose; break;
      case 's':
      case 'q':  verbose = (verbose > 0) ? 0 : verbose - 1; break;
      case 'h':
      case '?':  usage(); break;
      default:   if (opt(c, optarg))  break;
                 usage(); break;
      }
   }

   try
   {
      test(argc - optind, &argv[optind]);

      exit(status);
   }
   catch (L::Exception &e)
   {
      cerr << "\n\nException caught: " << e.what() << "\n\n";

      exit (20);
   }
   catch (std::exception &e)
   {
      cerr << "\n\nStandard exception caught: " << e.what() << "\n\n";

      exit (20);
   }
   catch (...)
   {
      cerr << "\n\nUnknown exception caught!\n\n";

      exit (20);
   }
}




