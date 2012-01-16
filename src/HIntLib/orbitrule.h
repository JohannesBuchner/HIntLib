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

/**
 *  Orbit Rule
 *
 *  Baseclass for all Cubature Rules using some kind of orbits as abscissa sets
 */

#ifndef HINTLIB_ORBIT_RULE_H
#define HINTLIB_ORBIT_RULE_H 1

#include <HIntLib/defaults.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
// #pragma interface
#endif

#include <HIntLib/hlmath.h>
#include <HIntLib/integrand.h>
#include <HIntLib/array.h>
#include <HIntLib/minmaxfinder.h>


namespace HIntLib
{

class OrbitRule
{
protected:
   OrbitRule (int dim) : dim(dim), p(dim) {}

   // No copy constructor or assignment, because these operations are private
   // for Point

   class Scaler
   {
   public:
      Scaler (const real* w, real s) : scale(s), width(w) {}
      real operator[] (int i) const { return width[i] * scale; }
   private:
      const real scale;
      const real *width;
   };

   // Sets p to the given center.
   // This method must be called before any of the following evalXXX() methods
   // are used.

   void setCenter (const real* center);

   // define some methods that are used for evaluating the integrand on
   // certain (fully) symmetric point sets.

   real eval0_0      (Integrand &);
   template<class T>
   real evalR0_0fs   (Integrand &, const real*, T r);
   real evalR0_0fs   (Integrand &, const real*, Point& a);
   template<class T1, class T2>
   int evalR0_0fs4d (Integrand &, const real*,
      real &sum0, T1 r1, real &sum1, T2 r2, real &sum2);
   real evalRR0_0fs  (Integrand &, const real*, const real*);
   real evalRS0_0fs  (Integrand &, const real*, const real*, const real*);
   real evalRRR0_0fs (Integrand &, const real*, const real*);
   real evalRRRR0_0fs(Integrand &, const real*, const real*);
   real evalRR0_0s   (Integrand &, const real*, const real*);
   real evalRS0_0s   (Integrand &, const real*, const real*, const real*);
   real evalR_Rfs    (Integrand &, const real*, const real*);
   real eval3powS    (Integrand &, const real*, const real*, real w0, real w1);

   // The following functions calculate the number of sampling points used
   // in each of these cubature formulas.

   // The formulas are structured such that none of the intermediate results
   // can overflow.

   Index num0_0       () const { return Index(1); }
   Index numR0_0fs    () const { return Index(2) * dim; }
   Index numRR0_0fs   () const { return Index(2) * dim * (dim-1); }
   Index numRS0_0fs   () const { return Index(4) * dim * (dim-1); }
   Index numRRR0_0fs  () const
      { return Index(dim) * (dim-1) * (dim-2) / 3 * 4; }
   Index numRRRR0_0fs () const
      { return Index(dim) * (dim-1) * (dim-2) / 3 * (dim-3) * 2; }
   Index numRR0_0s    () const
      { return (dim & 1) ? Index(dim) * ((dim-1) / 2) 
                         : Index(dim / 2) * (dim-1); }
   Index numRS0_0s    () const { return Index(dim) * (dim-1); }
   Index numR_Rfs     () const { return Index(1) << dim; }
   Index num3powS     () const { return powInt (Index(3), dim); }

protected:

   const int dim;
   mutable Point p;
};

}  // namespace HIntLib


/*******************  Implementation ****************/


inline
void
HIntLib::OrbitRule::setCenter (const real* center)
{
   for (int i = 0; i != dim; ++i)  p[i] = center[i];
}

inline
HIntLib::real
HIntLib::OrbitRule::eval0_0 (Integrand &f)
{
   return f(p);
}

template<class T>
inline
HIntLib::real
HIntLib::OrbitRule::evalR0_0fs (Integrand &f, const real* c, const T r)
{
   real sum = 0;

   for (int i = 0; i != dim; ++i)
   {
      p[i] = c[i] + r[i]; sum += f(p);
      p[i] = c[i] - r[i]; sum += f(p);
      p[i] = c[i];
   }

   return sum;
}

inline
HIntLib::real
HIntLib::OrbitRule::evalR0_0fs (Integrand &f, const real* c, Point& a)
{
   return evalR0_0fs (f, c, static_cast<real*> (a));
}

template<class T1, class T2>
inline
int
HIntLib::OrbitRule::evalR0_0fs4d (
   Integrand &f, const real* c,
   real &sum0_, T1 r1, real &sum1, T2 r2, real &sum2)
{
   MaxFinder<real> mf;
   int dimDiffMax = 0;

   // Keep a local sum0
   // This is required in the case that sum0 and sum1/sum2 reference the same
   // object. In this case sum0 would be overwritten before it is used to
   // calculate the 4th difference

   real sum0 = f(p);

   real ratio = r1[0] / r2[0];
   ratio *= ratio;

   for (int i = 0; i < dim; i++)
   {
      real f1a, f1b, f2a, f2b;
 
      p[i] = c[i] - r1[i]; sum1 += (f1a = f(p));
      p[i] = c[i] + r1[i]; sum1 += (f1b = f(p));
      p[i] = c[i] - r2[i]; sum2 += (f2a = f(p));
      p[i] = c[i] + r2[i]; sum2 += (f2b = f(p));
      p[i] = c[i];
 
      real diff = abs (f1a + f1b - 2 * sum0 - ratio * (f2a + f2b - 2 * sum0));

      if (mf << diff)  dimDiffMax = i;
   }

   sum0_ += sum0;

   return dimDiffMax;
}

inline
HIntLib::real
HIntLib::OrbitRule::evalRR0_0fs (Integrand &f, const real*c, const real* r)
{
   real sum = 0;

   for (int i = 0; i != dim - 1; ++i)
   {
      p[i] = c[i] - r[i];

      for (int j = i + 1; j != dim; ++j)
      {
         p[j] = c[j] - r[j]; sum += f(p);
         p[i] = c[i] + r[i]; sum += f(p);
         p[j] = c[j] + r[j]; sum += f(p);
         p[i] = c[i] - r[i]; sum += f(p);

         p[j] = c[j];  // Done with j  ->  Restore p [j]
      }
      p[i] = c[i];     // Done with i  ->  Restore p [i]
   }

   return sum;
}

inline
HIntLib::real
HIntLib::OrbitRule::evalRS0_0fs (
   Integrand &f, const real*c, const real* r, const real *s)
{
   real sum = 0;
 
   for (int i = 0; i != dim; ++i)
   {
      p[i] = c[i] - r[i];

      for (int j = 0; j != dim; ++j)
      {
         if (i != j)
         {
            p[j] = c[j] - s[j]; sum += f(p);
            p[i] = c[i] + r[i]; sum += f(p);
            p[j] = c[j] + s[j]; sum += f(p);
            p[i] = c[i] - r[i]; sum += f(p);

            p[j] = c[j];  // Done with j ->  Restore p [j]
         }
      }
      p [i] = c [i];   // Done with i  ->  Restore p [i]
   }
 
   return sum;
}


inline
HIntLib::real
HIntLib::OrbitRule::evalRR0_0s (Integrand &f, const real*c, const real* r)
{
   real sum = 0;

   for (int i = 0; i != dim - 1; ++i)
   {
      p[i] = c[i] + r[i];

      for (int j = i + 1; j != dim; ++j)
      {
         p[j] = c[j] + r[j]; sum += f(p);
         p[j] = c[j];
      }
      p[i] = c[i];
   }

   return sum;
}

inline
HIntLib::real
HIntLib::OrbitRule::evalRS0_0s (
      Integrand &f, const real*c, const real* r, const real* s)
{
   real sum = 0;

   for (int i = 0; i != dim; ++i)
   {
      p[i] = c[i] + r[i];

      for (int j = 0; j != dim; ++j)
      {
         if (i != j)
         {
            p[j] = c[j] + s[j]; sum += f(p);
            p[j] = c[j];
         }
      }
      p[i] = c[i];
   }

   return sum;
}

#endif


