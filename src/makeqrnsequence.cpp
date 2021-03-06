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

#define HINTLIB_LIBRARY_OBJECT

/**
 *  makeQRNSequence()
 *  makeQRNNet()
 *  getQRNSequenceName()
 *  getQRNNetName()
 *
 *  Creates instances of all implemented QRN Sequences and nets
 */

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation "qrnsequence.h"
#endif

#ifdef HINTLIB_HAVE_SSTREAM
#  include <sstream>
#else
#  include <HIntLib/fallback_sstream.h>
#endif

#ifdef HINTLIB_HAVE_CSTRING
#  include <cstring>
#  define HINTLIB_SSN std::
#else
#  include <string.h>
#  define HINTLIB_SSN
#endif

#include <memory>

#include <HIntLib/make.h>

#include <HIntLib/tparameter.h>
#include <HIntLib/qrnsequence.h>
#include <HIntLib/sobolmatrix.h>
#include <HIntLib/niederreitermatrix.h>
#include <HIntLib/niederreitermatrixgen.h>
#include <HIntLib/faure.h>
#include <HIntLib/halton.h>
#include <HIntLib/shiftnet.h>
#include <HIntLib/digitalnet2.h>
#include <HIntLib/digitalnetgen.h>
#include <HIntLib/generatormatrixgen.h>
#include <HIntLib/generatormatrixgenrow.h>
#include <HIntLib/generatormatrixvirtual.h>
#include <HIntLib/modulararithmetic.h>
#include <HIntLib/lookupfield.h>
#include <HIntLib/gf2.h>
#include <HIntLib/onedimvectorspace.h>
#include <HIntLib/prime.h>
#include <HIntLib/mersennetwister.h>
#include <HIntLib/mcpointset.h>

namespace L = HIntLib;


/*****************************************************************************/
/***            Generator Matrix 2                                         ***/
/*****************************************************************************/

namespace
{
   using namespace L;

   typedef GeneratorMatrix2<u64> GM2;
   typedef GeneratorMatrixGen<unsigned char> GMGen;

   const GM2* gm64 [2];
   PRNG *mt;

   /**
    *  init()
    */

   void init()
   {
      static int initDone = 0;

      if (initDone == 0xaffe)  return;

      initDone = 0xaffe;

      gm64 [0] = new SobolMatrix();
      gm64 [1] = new NiederreiterMatrix();
      mt = new PRNGImp<MersenneTwister> ();
   }

   const char* generatorNames [] = {
      "Sobol",
      "Niederreiter",
      "NiederreiterXing",
      "Random",
      "ShiftNet",
      "Zero"};

}  // end anomymous namespace


GM2* L::Make::generatorMatrix2 (int n, int dim)
{
   init();

   /*
    *  Generator Matrix in Base 2
    *
    *  0 - Sobol
    *  1 - Niederreiter
    *  2 - Niederreiter / Xing
    *  3 - Random
    *  4 - ShiftNet
    *  5 - Zero
    */

   switch (n)
   {
   case 0:  // Sobol
   case 1:  // Niederreiter
      {
         if (dim > gm64[n]->getDimension())  throw InvalidDimension(dim);
         return new GM2 (DiscardDimensions (dim, *gm64 [n]));
      }
   case 2:  // NX
      {
         if (dim == 0)  return new GM2 (0);

         std::auto_ptr<GMGen> m (loadNiederreiterXing (dim < 4 ? 4 : dim));

         return new GM2 (DiscardDimensions (dim, *m));
      }
   case 3:  // Random
      {
         GM2* m = new GM2 (dim);

         for (int d = 0; d < m->getDimension(); ++d)
         {
            for (int r = 0; r < m->getM(); ++r)
               for (int b = 0; b < m->getPrec(); ++b)
                  m->setd(d, r, b, mt->equidist (2));
         }
         return m;
      }
   case 4:  // ShiftNet
      {
         if (dim > 39)  throw InvalidDimension(dim);
         GM2* m = new GM2 (dim, dim);
         initShiftNet (*m);
         return m;
      }
   case 5:  // Zero
      return new GM2 (dim);
   }

   throw GeneratorMatrixDoesNotExist (n);
}

const char* L::Make::getGeneratorMatrix2Name (int n)
{
   init ();

   /*
    *  Generator Matrix in Base 2
    *
    *  0 - Sobol
    *  1 - Niederreiter
    *  2 - Niederreiter / Xing
    *  3 - Random
    *  4 - ShiftNet
    *  5 - Zero
    */

   if (n >= 0 && n < int (sizeof (generatorNames) / sizeof (generatorNames[0])))
   {
      return generatorNames [n];
   }

   throw GeneratorMatrixDoesNotExist (n);
}


/*****************************************************************************/
/***            Generator Matrix Gen                                       ***/
/*****************************************************************************/


GMGen* L::Make::generatorMatrixGen (int n, int dim, int m)
{
   init();

   /*
    *  Generator Matrix Gen
    *
    *  xx
    *  |
    *  \____ 00 Faure (prime base >= dim)
    *        01 Sobol (base = 2)
    *        06 Niederreiter / Xing
    *        10 ShiftNet
    *        12 Faure (prime power base >= dim)
    *        14 Random
    *        15 Zero
    *        xx base of Niederreiter Matrix
    */

   switch (n)
   {
   case 0:  // Faure
      {
         const int b = Prime::next (dim);
         GMGen* p = (m == -1 ) ? new GMGen (b, dim) : new GMGen (b, dim, m);
         initFaure (*p);
         return p;
      }
   case 1:  // Sobol
      {
         if (dim > gm64[0]->getDimension())  throw InvalidDimension(dim);
         if (m   > int(gm64[0]->getM()))     throw GM_CopyM(m, gm64[0]->getM());

         return new GMGen (DiscardDimensions (dim, AdjustM (m, *gm64 [0])));
      }
   case 6:  // NX
      {
         if (dim == 0) return new GMGen (2, 0);

         std::auto_ptr<GMGen> p (loadNiederreiterXing (dim < 4 ? 4 : dim));
         if (m > int(p->getM()))  throw GM_CopyM(m, p->getM());

         GMGen* pp = new GMGen (DiscardDimensions (dim, AdjustM (m, *p)));

         return pp;
      }
   case 10:  // shift net
      {
         std::auto_ptr<GM2> p (generatorMatrix2 (4, dim));
         return new GMGen (AdjustM (m, *p));
      }
   case 12:  // generalized faure
      {
         int b = dim;
         while (! Prime::isPrimePower (b))  ++b;
         GMGen* p = (m == -1 ) ? new GMGen (b, dim) : new GMGen (b, dim, m);
         initFaure (*p);
         return p;
      }
   case 14:  // Random, base 2
      {
         GMGen* p = (m == -1 ) ? new GMGen (2, dim) : new GMGen (2, dim, m);

         for (int d = 0; d < p->getDimension(); ++d)
         {
            for (int r = 0; r < p->getM(); ++r)
               for (int b = 0; b < p->getPrec(); ++b)
                  p->setd(d, r, b, mt->equidist (2));
         }
         return p;
      }
   case 15:
      return (m == -1 ) ? new GMGen (2, dim) : new GMGen (2, dim, m);
   default:
      if (n >= 1 && n < 100 && Prime::isPrimePower (n))  // Niederreiter
      {
         GMGen* p = (m == -1 ) ? new GMGen (n, dim) : new GMGen (n, dim, m, m);
         initNiederreiter (*p);
         return p;
      }
   }

   throw GeneratorMatrixDoesNotExist (n);
}

const char* L::Make::getGeneratorMatrixGenName (int n)
{
   /*
    *  Generator Matrix Gen
    *
    *  xx
    *  |
    *  \____ 00 Faure (base >= dim)
    *        01 Sobol (base = 2)
    *        06 Niederreiter / Xing
    *        10 ShiftNet
    *        12 Faure (prime power base >= dim)
    *        14 Random
    *        15 Zero
    *        xx base of Niederreiter Matrix
    */

   switch (n)
   {
   case  0:  return "Faure";
   case  1:  return "Sobol";
   case  6:  return "NiederreiterXing";
   case 10:  return "ShiftNet";
   case 12:  return "Faure (prime power base)";
   case 14:  return "Random (base 2)";
   case 15:  return "Zero (base 2)";
   }

   if (n >= 1 && n < 100)
   {
      if (Prime::isPrimePower (n))
      {
         static char s[100];

         std::ostringstream ss;
         ss << generatorNames [1] << "_base" << n;

         std::string str = ss.str();
         HINTLIB_SSN strcpy (s, str.c_str());
         return s;
      }
   }

   throw GeneratorMatrixDoesNotExist (n);
}


/*****************************************************************************/
/***            QRN Sequence                                               ***/
/*****************************************************************************/


/**
 *  makeQRNSequence()
 *
 *  Creates Generator # n.
 *
 *  The Generator is created on the heap.  The user is responsible for freeing
 *  it when it is not used any more.
 */

L::QRNSequence* L::Make::qrnSequence (int n, const Hypercube &h)
{
   /*
    *  1agg   Digital Sequences in Base 2
    *   | |
    *   | \__ Generator Matrix 2
    *   |
    *   \____ 1 - Normal
    *         2 - Gray-Code
    */

   if (n >= 1000 && n < 2000)
   {
      int nn = n;
      int gen = nn % 100; nn /= 100;

      GM2* gm;
      try
      {
         gm = generatorMatrix2 (gen, h.getDimension());
      }
      catch (GeneratorMatrixDoesNotExist &)
      {
         throw QRNSequenceDoesNotExist (n);
      }

      int type = nn % 10; nn /= 10;
      if (type >= 3)  throw QRNSequenceDoesNotExist (n);

      typedef DigitalNet2Gray<u64> DN2G;

      QRNSequence* seq;

      switch (type)
      {
         case 1:  seq = new QRNSequenceP<DN2G> (new DN2G (*gm, h)); break;
         case 2:  seq = new QRNSequenceP<DN2G> (new DN2G (*gm, h, false));
                  break;

         default: throw InternalError (__FILE__, __LINE__);
      }
      delete gm;
      return seq;
   }

   /*
    *  2agg   Digital Sequences in general base
    *   | |
    *   | \__ Generator Matrix Gen
    *   |
    *   \____ 0 - Naive
    *         1 - Normal
    *         2 - Gray-Code
    */

   if (n >= 2000 && n < 3000)
   {
      int nn = n;
      int gen = nn % 100; nn /= 100;

      GMGen* gm;
      try
      {
         gm = generatorMatrixGen (gen, h.getDimension());
      }
      catch (GeneratorMatrixDoesNotExist &)
      {
         throw QRNSequenceDoesNotExist (n);
      }

      int type = nn % 10; nn /= 10;
      if (type >= 3)  throw QRNSequenceDoesNotExist (n);

      int base = gm->getBase();
      bool prime = Prime::test (base);

      int vec = digitsRepresentable(static_cast<unsigned char>(base));

      QRNSequence* seq;

      if (vec == 1)
      // if (1)
      {
         // Field types

         typedef LookupField           <unsigned char> Fgen;
         typedef LookupFieldPow2       <unsigned char> Fpow2;
         typedef ModularArithmeticField<unsigned char> Fprime;

         // Field types used for Vector Space construction

         typedef LookupGaloisField    <unsigned char> GFgen;
         typedef LookupGaloisFieldPow2<unsigned char> GFpow2;

         // Vector Space types

         typedef OneDimVectorSpace<Fgen>   VSgen;
         typedef OneDimVectorSpace<Fpow2>  VSpow2;
         typedef OneDimVectorSpace<Fprime> VSprime;
         typedef OneDimVectorSpace<GF2>    VS2;

         // Digital Net types

         typedef DigitalNetGenNormal    <VSgen,   Index> NormalGen;
         typedef DigitalNetGenNaive     <VSgen,   Index> NaiveGen;
         typedef DigitalNetGenGray      <VSgen,   Index> GrayGen;
         typedef DigitalNetGenGray      <VSpow2,  Index> GrayPow2;
         typedef DigitalNetGenCyclicGray<VSprime, Index> CyclicPrime;
         typedef DigitalNetGenCyclicGray<VS2,     Index> Cyclic2;

         switch (type)
         {
            case 0:  seq = new QRNSequenceP<NaiveGen> (new NaiveGen (
                              VSgen(GFgen(base)), *gm, h));
                     break;
            case 1:  seq = new QRNSequenceP<NormalGen> (new NormalGen (
                              VSgen(GFgen(base)), *gm, h));
                     break;
            case 2:
            {
               if (base == 2)
                  seq = new QRNSequenceP<Cyclic2> (new Cyclic2 (
                           VS2(GF2()), *gm, h));
               else if (prime)
                  seq = new QRNSequenceP<CyclicPrime> (new CyclicPrime (
                           VSprime(Fprime(base)), *gm, h));
               else if (base % 2 == 0)
                  seq = new QRNSequenceP<GrayPow2> (new GrayPow2 (
                           VSpow2(GFpow2(base)), *gm, h));
               else
                  seq = new QRNSequenceP<GrayGen> (new GrayGen (
                           VSgen(GFgen(base)), *gm, h));
               break;
            }
            default: throw InternalError (__FILE__, __LINE__);
         }
      }
      else   // we can do vectorization
      {
         // Field types used for Vector Space construction

         typedef LookupGaloisField    <unsigned char> GFgen;
         typedef LookupGaloisFieldPow2<unsigned char> GFpow2;

         // Vector Space types

         typedef LookupVectorSpace    <unsigned char,unsigned char> VSgen;
         typedef LookupVectorSpacePow2<unsigned char,unsigned char> VSpow2;

         typedef VectorSpacePow2<u64> VSpow2full;

         // Digital Net types

         typedef DigitalNetGenNormal    <VSgen,   Index> NormalGen;
         typedef DigitalNetGenNaive     <VSgen,   Index> NaiveGen;
         typedef DigitalNetGenGray      <VSgen,   Index> GrayGen;
         typedef DigitalNetGenGray      <VSpow2,  Index> GrayPow2;
         typedef DigitalNetGenCyclicGray<VSgen,   Index> CyclicGen;
         typedef DigitalNetGenCyclicGray<VSpow2,  Index> CyclicPow2;

         VSgen x = VSgen(GFgen(base), vec);

         switch (type)
         {
            case 0:  seq = new QRNSequenceP<NaiveGen> (new NaiveGen (
                              VSgen(GFgen(base), vec), *gm, h));
                     break;
            case 1:  seq = new QRNSequenceP<NormalGen> (new NormalGen (
                              VSgen(GFgen(base), vec), *gm, h));
                     break;
            case 2:
            {
               if (base == 2)
                  seq = new QRNSequenceP<CyclicPow2> (new CyclicPow2 (
                           VSpow2(GFpow2(base), vec), *gm, h));
               else if (base % 2 == 0)
               {
                  seq = new QRNSequenceP<GrayPow2> (new GrayPow2 (
                           VSpow2(GFpow2(base), vec), *gm, h));
               }
               else if (prime)
                  seq = new QRNSequenceP<CyclicGen> (new CyclicGen (
                           VSgen(GFgen(base), vec), *gm, h));
               else
                  seq = new QRNSequenceP<GrayGen> (new GrayGen (
                           VSgen(GFgen(base), vec), *gm, h));
               break;
            }
            default: throw InternalError (__FILE__, __LINE__);
         }
      }

      delete gm;
      return seq;
   }

   /*
    *  9xxx    Special Cases
    *    |
    *    \____ 0 - Halton
    */

   switch (n)
   {
   case 9000: return new QRNSequenceT<Halton> (h);

   default: throw QRNSequenceDoesNotExist (n);
   }
}


/**
 *  get QRN Sequence Name()
 *
 *  Returns a short name for QRN Sequence # n.
 */

namespace
{
   const char* typeNames3 [3] = { "Naive", "Normal", "Gray-Code" };
}

const char* L::Make::getQrnSequenceName (int n)
{
   static char s [100];

   /*
    *  1agg   Digital Sequences in Base 2
    *   | |
    *   | \__ Generator Matrix 2
    *   |
    *   \____ 1 - Normal
    *         2 - Gray-Code
    */

   if (n >= 1000 && n < 2000)
   {
      int nn = n;
      int gen = nn % 100; nn /= 100;

      const char* genName;
      try
      {
         genName = getGeneratorMatrix2Name (gen);
      }
      catch (GeneratorMatrixDoesNotExist &)
      {
         throw QRNSequenceDoesNotExist (n);
      }

      int type = nn % 10; nn /= 10;
      if (type == 0 || type >= 3)  throw QRNSequenceDoesNotExist (n);

      std::ostringstream ss;
      ss << genName << '_' << typeNames3 [type];

      std::string str = ss.str();
      HINTLIB_SSN strcpy (s, str.c_str());
      return s;
   }

   /*
    *  2agg   Digital Sequences in general base
    *   | |
    *   | \__ Generator Matrix Gen
    *   |
    *   \____ 0 - Naive
    *         1 - Normal
    *         2 - Gray-Code
    */

   if (n >= 2000 && n < 3000)
   {
      int nn = n;
      int gen = nn % 100; nn /= 100;

      const char* genName;
      try
      {
         genName = getGeneratorMatrixGenName (gen);
      }
      catch (GeneratorMatrixDoesNotExist &)
      {
         throw QRNSequenceDoesNotExist (n);
      }

      int type = nn % 10; nn /= 10;
      if (type >= 3)  throw QRNSequenceDoesNotExist (n);

      std::ostringstream ss;
      ss << genName << '_' << typeNames3 [type];

      std::string str = ss.str();
      HINTLIB_SSN strcpy (s, str.c_str());
      return s;
   }

   /*
    *  9xxx    Special Cases
    *    |
    *    \____ 0 - Halton
    */

   switch (n)
   {
   case 9000: return "Halton";

   default: throw QRNSequenceDoesNotExist (n);
   }
}


/*****************************************************************************/
/***            QRN Net                                                    ***/
/*****************************************************************************/


namespace
{
   const int indices [10] = { 0, 1, 2, 16, 100, 0, 0, 0, 0, 0 };

   const DigitalNet::Truncation truncs [3]
      = { DigitalNet::TRUNCATE,
          DigitalNet::CENTER,
          DigitalNet::FULL };
}

/**
 *  makeQRNNet()
 *
 *  Creates Generator # n.
 *
 *  The Generator is created on the heap.  The user is responsible for freeing
 *  when it is not used any more.
 */

L::QRNSequence* L::Make::qrnNet (int n, const Hypercube &h, Index size)
{
   /*
    *  1rtg   Digital Sequences in Base 2
    *   |||
    *   ||\__ Generator Matrix 2
    *   ||
    *   |\___ 0 - Truncate
    *   |     1 - Center
    *   |     2 - Full
    *   |
    *   |     5 - Add equidistributed coordinate
    *   |
    *   \____ 0 - Index 0
    *         1 - Index 1
    *         2 - Index 2
    *         3 - Index 16
    *         4 - Index 100
    *
    *         9 - Randomize
    */

   if (n >= 1000 && n < 2000)
   {
      int nn = n;
      int gen = nn % 10; nn /= 10;

      GM2* gm;
      try
      {
         gm = generatorMatrix2 (gen, h.getDimension());
      }
      catch (GeneratorMatrixDoesNotExist &)
      {
         throw QRNSequenceDoesNotExist (n);
      }

      if (nn % 5 >= 3)  throw QRNSequenceDoesNotExist (n);
      DigitalNet::Truncation trunc = truncs [nn % 5]; nn /= 5;

      bool equi = nn % 2; nn /= 2;
      int index = nn % 10; nn /= 10;
      if (index > 4 && index < 9)  throw QRNSequenceDoesNotExist (n);

      typedef DigitalNet2Gray<u64> DN2G;

      DN2G* dn = new DN2G (*gm, h, ms1(size), indices[index], equi, trunc);
      if (index == 9)  dn->randomize(*mt);

      delete gm;

      return new QRNSequenceP<DN2G> (dn);
   }

#if 0
   /*
    *  2agg   Digital Sequences in general base
    *   | |
    *   | \__ Generator Matrix Gen
    *   |
    *   \____ 0 - Naive
    *         1 - Normal
    *         2 - Gray-Code
    */

   if (n >= 2000 && n < 3000)
   {
      int nn = n;
      int gen = nn % 100; nn /= 100;

      GMGen* gm;
      try
      {
         gm = generatorMatrixGen (gen);
      }
      catch (GeneratorMatrixDoesNotExist &)
      {
         throw QRNSequenceDoesNotExist (n);
      }

      int type = nn % 10; nn /= 10;
      if (type >= 3)  throw QRNSequenceDoesNotExist (n);

      int base = gm->getBase();
      QRNSequence* seq;

      typedef LookupField<unsigned char> Ring;
      LookupGaloisField<unsigned char> r (base);

      typedef DigitalNetGenNormal    <Ring, Index> Normal;
      typedef DigitalNetGenNaive     <Ring, Index> Naive;
      typedef DigitalNetGenGray      <Ring, Index> Gray;
      typedef DigitalNetGenCyclicGray<Ring, Index> Cyclic;

      switch (type)
      {
         case 0:  seq = new QRNSequenceP<Naive> (new Naive (*gm, r, h)); break;
         case 1:  seq = new QRNSequenceP<Normal> (new Normal (*gm, r, h));
                  break;
         case 2:
         {
            if (Prime::test (base))
               seq = new QRNSequenceP<Cyclic> (new Cyclic (*gm, r, h));
            else
               seq = new QRNSequenceP<Gray> (new Gray (*gm, r, h));

            break;
         }
         default: throw InternalError (__FILE__, __LINE__);
      }

      delete gm;
      return seq;
   }
#endif
   /*
    *  2agg   Digital Sequences in general base
    *   | |
    *   | \__ Generator Matrix Gen
    *   |
    *   \____ 1 - Fix one-dim projection
    *         2 - Fix two-dim projection
    */

   if (n >= 2000 && n < 3000)
   {
      int nn = n;
      int gen = nn % 100; nn /= 100;

      GMGen* gm;
      try
      {
         gm = generatorMatrixGen (gen, h.getDimension());
      }
      catch (GeneratorMatrixDoesNotExist &)
      {
         throw QRNSequenceDoesNotExist (n);
      }
      int base = gm->getBase();

      int m = logInt (size, Index (base));
      GMGen gm2 (AdjustM (m, AdjustPrec (m, *gm)));

      int fix = nn % 10; nn /= 10;
      if (fix > 2)  throw QRNSequenceDoesNotExist (n);

      if (fix)
      {
         GeneratorMatrixGenRow<unsigned char> gm3 (gm2);
         if (fix == 1)  fixOneDimensionalProjections (gm3);
         else           fixTwoDimensionalProjections (gm3);
         assign (gm3, gm2);
      }

      QRNSequence* seq;

      typedef OneDimVectorSpace<LookupField<unsigned char> > VS;
      LookupGaloisField<unsigned char> r (base);
      VS vs (r);

      typedef DigitalNetGenGray      <VS, Index> Gray;
      typedef DigitalNetGenCyclicGray<VS, Index> Cyclic;

      if (Prime::test (base))
         seq = new QRNSequenceP<Cyclic> (new Cyclic (vs, gm2, h));
      else
         seq = new QRNSequenceP<Gray> (new Gray (vs, gm2, h));

      delete gm;
      return seq;
   }

   /*
    *  9xxx    Special Cases
    *    |
    *    \____ 0 - Halton
    */

   switch (n)
   {
   case 9000: return new QRNSequenceT<Halton> (h);

   default: throw QRNSequenceDoesNotExist (n);
   }
}


/**
 *  getQRNSequenceName()
 *
 *  Returns a short name for Integrator # n.
 */

namespace
{
   const char* truncationNames [3] = { "Truncate", "Center", "Full" };
}

const char* L::Make::getQrnNetName (int n)
{
   static char s [100];

   /*
    *  1rtg   Digital Sequences in Base 2
    *   |||
    *   ||\__ Generator Matrix 2
    *   ||
    *   |\___ 0 - Truncate
    *   |     1 - Center
    *   |     2 - Full
    *   |
    *   |     5 - Add equidistributed coordinate
    *   |
    *   \____ 0 - Index 0
    *         1 - Index 1
    *         2 - Index 2
    *         3 - Index 16
    *         4 - Index 100
    *
    *         9 - Randomize
    */

   if (n >= 1000 && n < 2000)
   {
      int nn = n;
      int gen = nn % 10; nn /= 10;

      const char* matrixName;
      try
      {
         matrixName = getGeneratorMatrix2Name (gen);
      }
      catch (GeneratorMatrixDoesNotExist &)
      {
         throw QRNSequenceDoesNotExist (n);
      }

      int trunc = nn % 5; nn /= 5;
      if (trunc >= 3)  throw QRNSequenceDoesNotExist (n);

      bool equi = nn % 2; nn /= 2;
      int index = nn % 10; nn /= 10;
      if (index > 4 && index < 9)  throw QRNSequenceDoesNotExist (n);

      std::ostringstream ss;
      ss << matrixName << '_' << truncationNames [trunc];

      if (index > 0 && index < 9)  ss << '_' << indices[index];
      if (index == 9)  ss << "_Rand";
      if (equi)        ss << "_Equi";

      std::string str = ss.str();
      HINTLIB_SSN strcpy (s, str.c_str());
      return s;
   }

#if 0
   /*
    *  2agg   Digital Sequences in general base
    *   | |
    *   | \__ Generator Matrix Gen
    *   |
    *   \____ 0 - Naive
    *         1 - Normal
    *         2 - Gray-Code
    */

   if (n >= 2000 && n < 3000)
   {
      int nn = n;
      int gen = nn % 100; nn /= 100;

      const char* matrixName;
      try
      {
         matrixName = getGeneratorMatrixName (gen);
      }
      catch (GeneratorMatrixDoesNotExist &)
      {
         throw QRNSequenceDoesNotExist (n);
      }

      int trunc = nn % 10; nn /= 10;
      if (trunc >= 3)  throw QRNSequenceDoesNotExist (n);

      bool equi = nn % 2; nn /= 2;
      int index = nn % 10; nn /= 10;
      if (index > 4 && index < 9)  throw QRNSequenceDoesNotExist (n);

      std::ostringstream ss;
      ss << matrixName << '_' << truncationNames [trunc];

      if (index > 0 && index < 9)  ss << '_' << indices[index];
      if (index == 9)  ss << "_Rand";
      if (equi)        ss << "_Equi";

      std::string str = ss.str();
      HINTLIB_SSN strcpy (s, str.c_str());
      return s;
   }
#endif
   /*
    *  2agg   Digital Sequences in general base
    *   | |
    *   | \__ Generator Matrix Gen
    *   |
    *   \____ 1 - Fix one-dim projection
    *         2 - Fix two-dim projection
    */

   if (n >= 2000 && n < 3000)
   {
      int nn = n;
      int gen = nn % 100; nn /= 100;

      const char* matrixName;
      try
      {
         matrixName = getGeneratorMatrixGenName (gen);
      }
      catch (GeneratorMatrixDoesNotExist &)
      {
         throw QRNSequenceDoesNotExist (n);
      }
      int fix = nn % 10; nn /= 10;
      if (fix > 2)  throw QRNSequenceDoesNotExist (n);

      std::ostringstream ss;
      ss << matrixName;
      if      (fix == 1)  ss << "_fix1dim";
      else if (fix == 2)  ss << "_fix2dim";

      std::string str = ss.str();
      HINTLIB_SSN strcpy (s, str.c_str());
      return s;
   }


   /*
    *  9xxx    Special Cases
    *    |
    *    \____ 0 - Halton
    */

   switch (n)
   {
   case 9000: return "Halton";

   default: throw QRNSequenceDoesNotExist (n);
   }
}

