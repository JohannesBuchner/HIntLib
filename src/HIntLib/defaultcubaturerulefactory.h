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
 *  DefaultCubatureRuleFactory
 *  DefaultEmbeddedRuleFactory
 *
 *  Default implementation of CubatureRuleFactory and EmbeddedRuleFactory
 *
 *  This file should only be included in the *.cpp-file of cubature rules
 */
 
#ifndef DEFAULT_CUBATURERULE_FACTORY_H
#define DEFAULT_CUBATURERULE_FACTORY_H 1
 
#ifdef __GNUG__
#pragma interface
#endif

#include <HIntLib/embeddedrule.h>
 
namespace HIntLib
{
   template<class T>
   class DefaultCubatureRuleFactory : public CubatureRuleFactory
   {
   public:
      DefaultCubatureRuleFactory<T> () {}

      virtual T* create (unsigned dim) { return new T(dim); }
      virtual DefaultCubatureRuleFactory<T>* clone() const
         {  return new DefaultCubatureRuleFactory<T>; }
   };
 
   template<class T>
   class DefaultEmbeddedRuleFactory : public EmbeddedRuleFactory
   {
   public:
      DefaultEmbeddedRuleFactory<T> () {}

      virtual T* create (unsigned dim) { return new T(dim); }
      virtual DefaultEmbeddedRuleFactory<T>* clone() const
         {  return new DefaultEmbeddedRuleFactory<T>; }
   };
 
}  // namespace HIntLib
 
#endif

