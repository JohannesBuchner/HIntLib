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

#ifndef HINTLIB_POLYNOMIAL_2_BASE_H
#define HINTLIB_POLYNOMIAL_2_BASE_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <iosfwd>


namespace HIntLib
{
namespace Private
{

/**
 *  Polynomial 2 Ring Base
 */

class Polynomial2RingBase
{
public:
   struct unit_type {};

   static unsigned size()  { return 0; }
   static unsigned characteristic()  { return 2; }

   static unsigned numUnits()  { return 1; }
   static unit_type unitRecip (unit_type)  { return unit_type(); }
   static unit_type unitElement (unsigned) { return unit_type(); }
   static unsigned unitIndex (unit_type)  { return 0; }

   static void printSuffix (std::ostream &);
   void printVariable (std::ostream &) const;

protected:
   Polynomial2RingBase (char _var) : var (_var) {}
   const char var;
};

inline
bool
operator== (Polynomial2RingBase::unit_type, Polynomial2RingBase::unit_type)
{
   return true;
}

std::ostream& operator<< (std::ostream &, const Polynomial2RingBase &);


} // namespace Private
} // namespace HIntLib

#endif

