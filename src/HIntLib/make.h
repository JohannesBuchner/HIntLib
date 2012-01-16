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
 *  make.h
 */

#ifndef MAKE_H
#define MAKE_H 1

#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/exception.h>

namespace HIntLib
{
   class CubatureRuleFactory;
   class EmbeddedRuleFactory;
   template<class T> class HeapAllocatedGeneratorMatrix2;
   template<class T> class HeapAllocatedGeneratorMatrixGen;
   class QRNSequence;
   class Hypercube;

   class Make
   {
   public:

      /**
       *  Exceptions
       */
      
      class GeneratorMatrixDoesNotExist : public Exception
      {
         virtual void makeString() const;
         int number;
      public:
         GeneratorMatrixDoesNotExist (int n) : number(n) {}
         int getNumber() const  { return number; }
      };

      class QRNSequenceDoesNotExist : public Exception
      {
         virtual void makeString() const;
         int number;
      public:
         QRNSequenceDoesNotExist (int n) : number(n) {}
         int getNumber() const  { return number; }
      };

      class CubatureRuleDoesNotExist : public Exception
      {
         virtual void makeString() const;
         int number;
      public:
         CubatureRuleDoesNotExist(int n) : number(n) {}
         int getNumber() const  { return number; }
      };

      /**
       *  cubature Rule Factory()
       *  embedded Rule Factory()
       *
       *  get CubatureRule Name()
       *  get EmbeddedRule Name()
       *
       *  Creats factories for CubatureRules and EmbeddedRules
       */

      static CubatureRuleFactory* cubatureRuleFactory (int);
      static EmbeddedRuleFactory* embeddedRuleFactory (int);

      static const char* getCubatureRuleFactoryName (int);
      static const char* getEmbeddedRuleFactoryName (int);

      /**
       *  generator Matrix 2 ()
       *  generator Matrix Gen ()
       *
       *  get Generator Matrix 2 Name ()
       *  get Generator Matrix Gen Name ()
       */

      static HeapAllocatedGeneratorMatrix2<u64>*
         generatorMatrix2 (int, unsigned dim);
      static HeapAllocatedGeneratorMatrixGen<unsigned char>*
         generatorMatrixGen (int, unsigned dim);

      static const char* getGeneratorMatrix2Name (int);
      static const char* getGeneratorMatrixGenName (int);
      
      /**
       *  qrn Sequence ()
       *  qrn Net ()
       *
       *  get QRNSequence Name()
       *  get QRNNet Name()
       *
       *  Creates QRN Sequence/net # n.
       */

      static QRNSequence* qrnSequence (int, const Hypercube &);
      static QRNSequence* qrnNet (int, const Hypercube &, Index);

      static const char* getQrnSequenceName (int);
      static const char* getQrnNetName (int);
   };

}  // namespace HIntLib

#endif

