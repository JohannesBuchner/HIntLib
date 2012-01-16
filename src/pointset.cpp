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
 *  PRNG
 *
 *  Abstract base class for Pseudo Random Number Generators
 */


#ifdef __GNUG__
#pragma implementation
#endif

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/pointset.h>

#include <HIntLib/mersennetwister.h>
#include <HIntLib/mcpointset.h>

namespace L = HIntLib;

void L::PartitionablePointSet::integrate (
   real *point, Integrand &f, Index num, Stat &stat)
{
   integratePartition (point, f, num, 0, num, stat);
}

void L::PartitionablePointSet::integrate (
   real *point, Integrand &f, Index num, StatVar &stat)
{
   integratePartition (point, f, num, 0, num, stat);
}

void L::PartitionablePointSet::doJob (real *point, Job &job, Index num)
{
   doJobPartition (point, job, num, 0, num);
}

