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

#ifndef TEST_H
#define TEST_H

#include <HIntLib/defaults.h>

namespace L = HIntLib;

extern int verbose;
extern const char* options;

extern const char option_msg [];

void error();
void error(const char*);
void usage();

void test (int argc, char** argv);
bool opt (int, const char*);

#define SILENT if (verbose <= 0)
#define NORMAL if (verbose >= 1)
#define DEB1   if (verbose >= 2)
#define DEB2   if (verbose >= 3)
#define DEB3   if (verbose >= 4)
#define DEB4   if (verbose >= 5)

#endif

