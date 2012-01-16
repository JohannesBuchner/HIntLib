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
 *  test.cpp
 *  test.h
 *
 *  Provides a common main() routine and some utility functions shared by all
 *  test programs.
 */

#include <iostream>
#include <string>
#include <stdexcept>
#include <unistd.h>

#include <HIntLib/exception.h>
#include <HIntLib/output.h>

#include "test.h"

using std::cerr;
using std::cout;
using std::endl;

enum { CURRENT_YEAR = 2006 };

int verbose = 1;
int status  = 0;
#ifdef HINTLIB_ENCODING_LOCALE
bool utf8 = false;
#endif

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
     cout << "\nToo many errors encountered. Terminating test!\n" << endl;

     exit (status);
  }
}

//  Displays an error message before calling error().

void error(const char* s)
{
   cout << "\nError: " << s << "!!!\n" << endl;
   error();
}


/**
 *  printHeader()
 */

void
printHeader (std::ostream& o)
{
   o << HINTLIB_PACKAGE_STRING " " << Wgl4Ascii ("\xe2\x80\x94", "--")
     << ' ' << testProgramName
     << "\nCopyright "
     << Latin1Ascii ("\xc2\xa9 ", "\xa9 ", "(C) ");

   for (int i = testProgramCopyright; ; ++i)
   {
      o << i;
      if (i == CURRENT_YEAR)  break;
      o << ',';
   }

   o << Latin1Ascii (" by Rudolf Sch\xc3\xbcrer\n\n",
                     " by Rudolf Sch\xfcrer\n\n",
                     " by Rudolf Schuerer\n\n");
}

#ifdef HINTLIB_BUILD_WCHAR
void
printHeader (std::wostream& o)
{
#if HINTLIB_CHARACTER_SET >= 3
   o << HINTLIB_PACKAGE_STRING << L" \x2014 " << testProgramName
     << L"\nCopyright \x00a9 ";
   
   for (int i = testProgramCopyright; ; ++i)
   {
      o << i;
      if (i == CURRENT_YEAR)  break;
      o << L',';
   }

   o << L" by Rudolf Sch\x00fcrer\n\n"; 
#elif HINTLIB_CHARACTER_SET >= 1
   o << HINTLIB_PACKAGE_STRING << L" -- " << testProgramName
     << L"\nCopyright \x00a9 ";
   
   for (int i = testProgramCopyright; ; ++i)
   {
      o << i;
      if (i == CURRENT_YEAR)  break;
      o << L',';
   }

   o << L" by Rudolf Sch\x00fcrer\n\n"; 
#else
   o << HINTLIB_PACKAGE_STRING << L" -- " << testProgramName
     << L"\nCopyright (C) ";
   
   for (int i = testProgramCopyright; ; ++i)
   {
      o << i;
      if (i == CURRENT_YEAR)  break;
      o << L',';
   }

   o << L" by Rudolf Schuerer\n\n"; 
#endif
}
#endif


/**
 *  doubleQuote()
 */

void doubleQuote (std::ostream& o, const char* s)
{
   o << Wgl4Ascii ("\xe2\x80\x9c", "\"")  // LEFT DOUBLE QUOTATION MARK
     << s
     << Wgl4Ascii ("\xe2\x80\x9d", "\""); // RIGHT DOUBLE QUOTATION MARK
}

#ifdef HINTLIB_BUILD_WCHAR
void doubleQuote (std::wostream& o, const char* s)
{
#if HINTLIB_CHARACTER_SET >= 3
   o << L'\x201c'  // LEFT DOUBLE QUOTATION MARK
     << s
     << L'\x201d'; // RIGHT DOUBLE QUOTATION MARK
#else
   o << L'"' << s << L'"';
#endif
}

void doubleQuote (std::wostream& o, const wchar_t* s)
{
#if HINTLIB_CHARACTER_SET >= 3
   o << L'\x201c'  // LEFT DOUBLE QUOTATION MARK
     << s
     << L'\x201d'; // RIGHT DOUBLE QUOTATION MARK
#else
   o << L'"' << s << L'"';
#endif
}
#endif


/**
 *  usage()
 */

void usage(const char* msg)
{
   printHeader (cerr);

   if (msg)  cerr << msg << "\n\n";

   cerr << "Usage: " << testProgramName << ' '
        << testProgramParameters << "\n\n"
        << testProgramUsage
        << option_msg
        << "  -s     Silent mode (verbosity level 0)\n"
           "  -q     Quiet. Same as -s\n"
           "  -v     Increases verbosity level by one "
                    "(can be applied multiple times)\n"
           "  -h     Display usage information\n\n";

   exit (1);
}


/**
 *  main()
 */

int main (int argc, char** argv)
{
   // Initialize output stream

#ifdef HINTLIB_STREAMS_SUPPORT_LOCALE
   if (! setlocale (LC_CTYPE, ""))
   {
      cerr << ("Locale is invalid (setlocale() returned null)!\n");
   }

   try
   {
      cout.imbue (std::locale(std::locale::classic(),
                              std::locale(""),
                              std::locale::ctype));
   }
   catch (std::runtime_error&)
   {
      cerr << ("Locale is invalid (std::locale() threw exeception)!\n");
   }
#endif

#ifdef HINTLIB_ENCODING_LOCALE
   utf8 = HIntLib::Private::utf8Support(cout);
#endif
   
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
      case 'h':  usage(0); break;
      default:
         {
            if (opt(c, optarg))  break;
            std::ostringstream ss;
            ss << "Invalid option '" << char(optopt) << "'!";
            usage (ss.str().c_str());
            break;
         }
      }
   }

   try
   {
      test(argc - optind, &argv[optind]);

      exit(status);
   }
   catch (L::Exception &e)
   {
      cerr << "\n\nHIntLib exception caught: " << e.what() << "\n\n";

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

