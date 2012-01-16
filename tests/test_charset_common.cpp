/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2006  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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

#include <iomanip>

#include "test.h"
#include "test_charset.h"

using std::setw;
using std::setfill;
using std::ostream;

const char options[] = "";
const char option_msg[] = "";
const char testProgramParameters[] = "[OPTION]...";
const char testProgramUsage[] =
   "Lists all non-ASCII characters used by HIntLib "
   "and the HIntLib test programs.\n\n";
const int  testProgramCopyright = 2006;

bool opt(int, const char*)
{
   return false;
}


/**
 *  operator<< for Type
 */

ostream& operator<< (ostream& o, Type t)
{
   switch (t)
   {
      case ASCII:   o << "ASCII";   break;
      case LATIN1:  o << "LATIN1";  break;
      case WGL4:    o << "WGL4";    break;
      case UNICODE: o << "UNICODE"; break;
   }

   return o;
}


/**
 *  Character Data
 *
 *  The following non-ASCII characters are currently used by HIntLib:
 */

const Data characters [] =
{
//   UNICODE Subset   Name
   { 0x00a9, LATIN1,  "COPYRIGHT SIGN" },                     // test_*

   { 0x00B1, LATIN1,  "PLUS-MINUS SIGN" },
   { 0x00B2, LATIN1,  "SUPERSCRIPT TWO" },
   { 0x00B3, LATIN1,  "SUPERSCRIPT THREE" },

   { 0x00B7, LATIN1,  "MIDDLE DOT" },                         // test_arithmetic

   { 0x00B9, LATIN1,  "SUPERSCRIPT ONE" },

   { 0x00BC, LATIN1,  "VULGAR FRACTION ONE QUARTER" },
   { 0x00BD, LATIN1,  "VULGAR FRACTION ONE HALF" },
   { 0x00BE, LATIN1,  "VULGAR FRACTION THREE QUARTERS" },

   { 0x00D7, LATIN1,  "MULTIPLICATION SIGN" },

   { 0x00FC, LATIN1,  "LATIN SMALL LETTER U WITH DIAERESIS" },// test_*

   { 0x0394, WGL4,    "GREEK CAPITAL LETTER DELTA" },         // test_rule

   { 0x03B1, WGL4,    "GREEK SMALL LETTER ALPHA" },           // test_arithmetic

   { 0x03B5, WGL4,    "GREEK SMALL LETTER EPSILON" },         // test_rule

   { 0x2014, WGL4,    "EM DASH" },                            // test_*

   { 0x2016, UNICODE, "DOUBLE VERTICAL LINE" },               // test_arithmetic

   { 0x201C, WGL4,    "LEFT DOUBLE QUOTATION MARK" },         // test_*
   { 0x201D, WGL4,    "RIGHT DOUBLE QUOTATION MARK" },        // test_*

   { 0x2070, UNICODE, "SUPERSCRIPT ZERO" },

   { 0x2074, UNICODE, "SUPERSCRIPT FOUR" },
   { 0x2075, UNICODE, "SUPERSCRIPT FIVE" },
   { 0x2076, UNICODE, "SUPERSCRIPT SIX" },
   { 0x2077, UNICODE, "SUPERSCRIPT SEVEN" },
   { 0x2078, UNICODE, "SUPERSCRIPT EIGHT" },
   { 0x2079, UNICODE, "SUPERSCRIPT NINE" },

   { 0x207F, WGL4,    "SUPERSCRIPT LATIN SMALL LETTER N" },   // test_arithmetic

   { 0x2080, UNICODE, "SUPSCRIPT ZERO" },
   { 0x2081, UNICODE, "SUPSCRIPT ONE" },
   { 0x2082, UNICODE, "SUPSCRIPT TWO" },
   { 0x2083, UNICODE, "SUPSCRIPT THREE" },
   { 0x2084, UNICODE, "SUPSCRIPT FOUR" },
   { 0x2085, UNICODE, "SUPSCRIPT FIVE" },
   { 0x2086, UNICODE, "SUPSCRIPT SIX" },
   { 0x2087, UNICODE, "SUPSCRIPT SEVEN" },
   { 0x2088, UNICODE, "SUPSCRIPT EIGHT" },
   { 0x2089, UNICODE, "SUPSCRIPT NINE" },
   
   { 0x2102, UNICODE, "DOUBLE-STRUCK CAPITAL C" },
   { 0x211A, UNICODE, "DOUBLE-STRUCK CAPITAL Q" },
   { 0x211D, UNICODE, "DOUBLE-STRUCK CAPITAL R" },
   { 0x2124, UNICODE, "DOUBLE-STRUCK CAPITAL Z" },
   
   { 0x2153, UNICODE, "VULGAR FRACTION ONE THIRD" },
   { 0x2154, UNICODE, "VULGAR FRACTION TWO THIRDS" },
   { 0x2155, UNICODE, "VULGAR FRACTION ONE FIFTH" },
   { 0x2156, UNICODE, "VULGAR FRACTION TWO FIFTHS" },
   { 0x2157, UNICODE, "VULGAR FRACTION THREE FIFTHS" },
   { 0x2158, UNICODE, "VULGAR FRACTION FOUR FIFTHS" },
   { 0x2159, UNICODE, "VULGAR FRACTION ONE SIXTH" },
   { 0x215A, UNICODE, "VULGAR FRACTION FIVE SIXTHS" },
   { 0x215B, WGL4,    "VULGAR FRACTION ONE EIGHTH" },
   { 0x215C, WGL4,    "VULGAR FRACTION THREE EIGHTHS" },
   { 0x215D, WGL4,    "VULGAR FRACTION FIVE EIGHTHS" },
   { 0x215E, WGL4,    "VULGAR FRACTION SEVEN EIGHTHS" },

   { 0x2212, WGL4,    "MINUS SIGN" },

   { 0x2215, WGL4,    "DIVISION SLASH" },

   { 0x221E, WGL4,    "INFINITY" },                           // test_arithmetic

   { 0x2223, UNICODE, "DIVIDES" },                            // test_arithmetic

   { 0x223C, UNICODE, "TILDE OPERATOR" },                     // test_arithmetic

   { 0x22C5, UNICODE, "DOT OPERATOR" },                       // test_arithmetic

   { 0x2329, UNICODE, "LEFT-POINTING ANGLE BRACKET" },        // test_arithmetic
   { 0x232A, UNICODE, "RIGHT-POINTING ANGLE BRACKET" },       // test_arithmetic

   { 0x25CB, WGL4,    "WHITE CIRCLE" },                       // test_combin

   { 0x25CF, WGL4,    "BLACK CIRCLE" },                       // test_combin

#if 0
   { 0x27E8, UNICODE, "MATHEMATICAL LEFT ANGLE BRACKET" },    // test_arithmetic
   { 0x27E9, UNICODE, "MATHEMATICAL RIGHT ANGLE BRACKET" },   // test_arithmetic
#endif
};

const unsigned numChars = sizeof(characters) / sizeof(Data);


/**
 *  encode_utf8()
 */

void encode_utf8 (unsigned long c, char* buffer)
{
   if (c < 0x80)  // 0 - 127 => 1 byte
   {
      *buffer++ = c;
   }
   else if (c < 0x800)  // 128 - 2047 => 2 bytes
   {
      *buffer++ =  (c >> 6)       | 192;
      *buffer++ = ( c       & 63) | 128;
   }
   else if (c < 0x10000)  // 2048 - 66535 => 3 bytes
   {
      *buffer++ = ( c >> 12)       | 224;
      *buffer++ = ((c >>  6) & 63) | 128;
      *buffer++ = ( c        & 63) | 128;
   }
   else  // everything else:  4 bytes
   {
      *buffer++ = ( c >> 18)       | 240;
      *buffer++ = ((c >> 12) & 63) | 128;
      *buffer++ = ((c >>  6) & 63) | 128;
      *buffer++ = ( c        & 63) | 128;
   }

   *buffer = '\0';
}


/**
 *  beginOfLine()
 */

std::string beginOfLine (const Data& data)
{
   unsigned long c = data.code;
   char buffer [5];
   encode_utf8 (c, buffer);

   std::ostringstream utf8Formatted;
   utf8Formatted.setf (ostream::uppercase);
   utf8Formatted.setf (ostream::hex, ostream::basefield);
   utf8Formatted.setf (ostream::right, ostream::adjustfield);
   utf8Formatted << setfill ('0');
   for (char* i = buffer; *i != '\0'; ++i)
   {
      unsigned char utf8Char = *i;
      utf8Formatted << setw(2) << unsigned(utf8Char) << ' ';
   }

   std::ostringstream line;

   line.setf (ostream::uppercase);
   line.setf (ostream::hex, ostream::basefield);
   line.setf (ostream::right, ostream::adjustfield);
   line << "U+" << setfill('0') << setw(4) << data.code;
   line.setf (ostream::left, ostream::adjustfield);
   line << "  " << setfill(' ') << setw(10) << utf8Formatted.str().c_str()
        << setw(9) << data.type;

   return line.str();
}


