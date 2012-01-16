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

#ifndef HINTLIB_POLYNOMIAL_BASE_H
#define HINTLIB_POLYNOMIAL_BASE_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <iosfwd>

namespace HIntLib
{
namespace Private
{

/**
 *  Polynomial Ring Base  --  no template
 *
 *  Non-template base class of polynomial rings
 */

class PRB
{
public:
   static unsigned size()  { return 0; }

   void printVariable (std::ostream &) const;
   void printVariableWithBrackets (std::ostream &) const;
   void printVariablePow (std::ostream&, unsigned) const;

protected:
   explicit PRB (char _var) : var (_var) {}

   static const int squareBeatsLinear [];

private:
   const char var;
};

unsigned funnySum (int, unsigned);

} // namespace Private
} // namespace HIntLib

#endif

