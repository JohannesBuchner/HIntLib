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

#include <HIntLib/defaults.h>

#ifdef HINTLIB_HAVE_CSTRING
#  include <cstring>
#  define HINTLIB_SSN std::
#else
#  include <string.h>
#  define HINTLIB_SSN
#endif

#define HINTLIB_LIBRARY_OBJECT

#include <HIntLib/make.h>

#include <HIntLib/pseudoembeddedrule.h>

#include <HIntLib/rule1midpoint.h>
#include <HIntLib/rule1trapezoidal.h>
#include <HIntLib/rule2simplex.h>
#include <HIntLib/rule2thacher.h>
#include <HIntLib/rule2ionescu.h>
#include <HIntLib/rule3octahedron.h>
#include <HIntLib/rule3cross.h>
#include <HIntLib/rule3tyler.h>
#include <HIntLib/rule3gauss.h>
#include <HIntLib/rule3ewing.h>
#include <HIntLib/rule3simpson.h>
#include <HIntLib/rule5hammer.h>
#include <HIntLib/rule5mustardlynessblatt.h>
#include <HIntLib/rule5stroud.h>
#include <HIntLib/rule5stroud2.h>
#include <HIntLib/rule5gauss.h>
#include <HIntLib/rule7phillips.h>
#include <HIntLib/rule75genzmalik.h>
#include <HIntLib/rule9stenger.h>
#include <HIntLib/rulegauss.h>

namespace L = HIntLib;

using L::CubatureRuleFactory;
using L::EmbeddedRuleFactory;

namespace
{

CubatureRuleFactory* makeGauss1() { return L::RuleGauss::getFactory(1); }
CubatureRuleFactory* makeGauss2() { return L::RuleGauss::getFactory(2); }
CubatureRuleFactory* makeGauss3() { return L::RuleGauss::getFactory(3); }
CubatureRuleFactory* makeGauss4() { return L::RuleGauss::getFactory(4); }
CubatureRuleFactory* makeGauss5() { return L::RuleGauss::getFactory(5); }
CubatureRuleFactory* makeGauss6() { return L::RuleGauss::getFactory(6); }
CubatureRuleFactory* makeGauss7() { return L::RuleGauss::getFactory(7); }
CubatureRuleFactory* makeGauss8() { return L::RuleGauss::getFactory(8); }
CubatureRuleFactory* makeGauss9() { return L::RuleGauss::getFactory(9); }

struct RuleRecord
{
   int n;
   const char* name;
   CubatureRuleFactory* ((*function)());
};

const RuleRecord factories [] =
{
   {  11, "1-Midpoint",           L::Rule1Midpoint::getFactory },
   {  12, "1-Trapezoidal",        L::Rule1Trapezoidal::getFactory },
   {  13, "1-Gauss1D",            makeGauss1 },
   {  21, "2-Simplex",            L::Rule2Simplex::getFactory },
   {  22, "2-Thacher",            L::Rule2Thacher::getFactory },
   {  23, "2-Ionescu",            L::Rule2Ionescu::getFactory },
   {  31, "3-Octahedron",         L::Rule3Octahedron::getFactory },
   {  32, "3-Cross",              L::Rule3Cross::getFactory },
   {  33, "3-Tyler",              L::Rule3Tyler::getFactory },
   {  34, "3-Gauss",              L::Rule3Gauss::getFactory },
   {  35, "3-Ewing",              L::Rule3Ewing::getFactory },
   {  36, "3-Simpson",            L::Rule3Simpson::getFactory },
   {  37, "3-Gauss1D",            makeGauss2 },
   {  51, "5-Hammer",             L::Rule5Hammer::getFactory },
   {  52, "5-Stroud",             L::Rule5Stroud::getFactory },
   {  53, "5-Gauss",              L::Rule5Gauss::getFactory },
   {  54, "5-Gauss1D",            makeGauss3 },
   {  55, "5-Stroud2",            L::Rule5Stroud2::getFactory },
   {  56, "5-MustardLynessBlatt", L::Rule5MustardLynessBlatt::getFactory },
   {  71, "7-GenzMalik",          (CubatureRuleFactory* (*)())
                                        L::Rule75GenzMalik::getFactory },
   {  72, "7-Phillips",           L::Rule7Phillips::getFactory },
   {  73, "7-Gauss1D",            makeGauss4 },
   {  91, "9-Stenger",            L::Rule9Stenger::getFactory },
   {  92, "9-Gauss1D",            makeGauss5 },
   { 111, "11-Gauss1D",           makeGauss6 },
   { 131, "13-Gauss1D",           makeGauss7 },
   { 151, "15-Gauss1D",           makeGauss8 },
   { 171, "17-Gauss1D",           makeGauss9 },
};

inline bool operator< (const RuleRecord& r1, const RuleRecord& r2)
{
   return r1.n < r2.n;
}

const RuleRecord* lookup (int n)
{
   RuleRecord r = {n, 0, 0};
   const RuleRecord* factoriesEnd =
      factories + sizeof (factories) / sizeof (RuleRecord);

   const RuleRecord* p = std::lower_bound (factories, factoriesEnd, r);

   if (p == factoriesEnd || p->n != n)
   {
      throw L::Make::CubatureRuleDoesNotExist (n);
   }

   return p;
}

}  // anonymous namespace

#ifdef HINTLIB_INSTANTIATE_STL
#ifdef HINTLIB_SGI
   template const RuleRecord* std::__lower_bound<> (
      const RuleRecord*, const RuleRecord*, const RuleRecord&, ptrdiff_t*);
#endif
#endif

CubatureRuleFactory* L::Make::cubatureRuleFactory (int rule)
{
   return lookup (rule)->function();
}

const char* L::Make::getCubatureRuleFactoryName (int rule)
{
   return lookup (rule)->name;
}

namespace
{
   int rulePairs [][2] =
   {  {31, 11},  // Octahedron / Midpoint
      {51, 31},  // Hammer / Octahedron
      {72, 51},  // Phillips / Hammer
      {72, 52},  // Phillips / Stroud
      {91, 72},  // Stenger / Phillips
   };


EmbeddedRuleFactory* makePERF (int rule1, int rule2)
{
   CubatureRuleFactory *p1 = 0, *p2 = 0;

   try
   {
      p1 = L::Make::cubatureRuleFactory (rule1);
      p2 = L::Make::cubatureRuleFactory (rule2);
      return new L::PseudoEmbeddedRuleFactory (p1, p2);
   }
   catch (...)
   {
      delete p1;
      delete p2;
      throw;
   }
}

EmbeddedRuleFactory* makePDERF (int rule1, int rule2, int rule3)
{
   CubatureRuleFactory *p1 = 0, *p2 = 0, *p3 = 0;

   try
   {
      p1 = L::Make::cubatureRuleFactory (rule1);
      p2 = L::Make::cubatureRuleFactory (rule2);
      p3 = L::Make::cubatureRuleFactory (rule3);
      return new L::PseudoDoubleEmbeddedRuleFactory (p1, p2, p3);
   }
   catch (...)
   {
      delete p1;
      delete p2;
      delete p3;
      throw;
   }
}

const char* namePERF (int rule1, int rule2)
{
   static char s [100];  // FIXME   Use string/sstream

   const char* s1 = L::Make::getCubatureRuleFactoryName (rule1);
   const char* s2 = L::Make::getCubatureRuleFactoryName (rule2);

   // FIXME this fails if the degree is >= 10

   s[0] = s1[0];
   s[1] = '/';
   s[2] = s2[0];
   s[3] = '-';
   s[4] = '\0';

   HINTLIB_SSN strcat (s, s1+2);
   HINTLIB_SSN strcat (s, "/");
   HINTLIB_SSN strcat (s, s2+2);

   return s;
}

const char* namePDERF (int rule1, int rule2, int rule3)
{
   static char s [150];  // FIXME   Use string/sstream

   const char* s1 = L::Make::getCubatureRuleFactoryName (rule1);
   const char* s2 = L::Make::getCubatureRuleFactoryName (rule2);
   const char* s3 = L::Make::getCubatureRuleFactoryName (rule3);

   // FIXME this fails if the degree is >= 10

   s[0] = s1[0];
   s[1] = '/';
   s[2] = s2[0];
   s[3] = '/';
   s[4] = s3[0];
   s[5] = '-';
   s[6] = '\0';

   HINTLIB_SSN strcat (s, s1+2);
   HINTLIB_SSN strcat (s, "/");
   HINTLIB_SSN strcat (s, s2+2);
   HINTLIB_SSN strcat (s, "/");
   HINTLIB_SSN strcat (s, s3+2);

   return s;
}

}  // anonymous namespace

EmbeddedRuleFactory* L::Make::embeddedRuleFactory (int rule)
{
   if (rule >= 1 && rule <= 5)
   {
      return makePERF (rulePairs [rule][0], rulePairs [rule][1]);
   }

   switch (rule)
   {
   case 11: return makePDERF (9, 6, 7);    // Phillips / Hammer / Stroud

   case 21: return Rule75GenzMalik::getFactory();

   default: throw CubatureRuleDoesNotExist (rule);
   }
}

const char* L::Make::getEmbeddedRuleFactoryName (int rule)
{
   if (rule >= 1 && rule <= 5)
   {
      return namePERF (rulePairs [rule][0], rulePairs [rule][1]);
   }

   switch (rule)
   {
   case 11: return namePDERF (9, 6, 7);    // Phillips / Hammer / Stroud

   case 21: return "7/5-GenzMalik";

   default: throw CubatureRuleDoesNotExist (rule);
   }
}

