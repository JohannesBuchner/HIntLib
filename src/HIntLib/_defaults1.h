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

/**
 *  defaults.h
 *
 *  Defines a number of global data types, most noteworthy the type used for
 *  normal floating-point arithmetic (real) and for counting loop iterations
 *  for long-running loops like integrand evaluations, sequence-numbers,...
 */

#ifndef HINTLIB_DEFAULTS_H
#define HINTLIB_DEFAULTS_H 1

#ifdef __GNUG__
#pragma interface
#endif

namespace HIntLib
{

