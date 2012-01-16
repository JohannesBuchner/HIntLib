/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration 
 *
 *  Copyright (C) 2002  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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

#include <HIntLib/gcd.h>

/**
 * genGcd() - both multiplicators
 */

template<class A>
typename A::type HIntLib::genGcd (
   const A& a,
   const typename A::type&  v, const typename A::type&  u,
         typename A::type& mv,       typename A::type& mu)
{
   typedef typename A::type T;
   
   T q;

   if (a.is0 (v))
   {
      destructiveAssign (mv, q);  // zero

      if (&mu == &u)
      {
         T temp (mu);
         mu = a.one();
         return temp;
      }
      else
      {
         mu = a.one();
         return u;
      }
   }

   T u3;
   T u2;
   a.div (u, v, u2, u3);
   a.negate (u2);

   if (a.is0 (u3))
   {
      destructiveAssign (mu, q);  // zero

      if (&mv == &v)
      {
         T temp (mv);
         mv = a.one();
         return temp;
      }
      else
      {
         mv = a.one();
         return v;
      }
   }

   T v3;
   T v1;
   a.div (v, u3, v1, v3);
   a.negate (v1);
   T v2 (a.add (a.one(), a.mul (v1, u2)));

   T u1 = a.one();

   for (;;)
   {
      if (a.is0 (v3))
      {
         destructiveAssign (mu, u1);
         destructiveAssign (mv, u2);

         return u3;
      }

      a.div (u3, v3, q, u3);
      a.subFrom (u1, a.mul (q, v1));
      a.subFrom (u2, a.mul (q, v2));

      if (a.is0 (u3))
      {
         destructiveAssign (mu, v1);
         destructiveAssign (mv, v2);

         return v3;
      }

      a.div (v3, u3, q, v3);
      a.subFrom (v1, a.mul (q, u1));
      a.subFrom (v2, a.mul (q, u2));
   }
}


/**
 * genGcd() - one multiplicator
 */

template<class A>
typename A::type HIntLib::genGcd (
   const A& a,
   const typename A::type&  v, const typename A::type&  u,
         typename A::type& mv)
{
   typedef typename A::type T;
   
   T q;

   if (a.is0 (v))
   {
      destructiveAssign (mv, q);  // zero
      return u;
   }

   T u3;
   T u2;
   a.div (u, v, u2, u3);
   a.negate (u2);

   if (a.is0 (u3))
   {
      if (&v == &mv)
      {
         T temp (mv);
         mv = a.one();
         return temp;
      }
      else
      {
         mv = a.one();
         return v;
      }
   }

   T v3;
   a.div (v, u3, q, v3);
   T v2 (a.sub (a.one(), a.mul (q, u2)));

   for (;;)
   {
      if (a.is0 (v3))
      {
         destructiveAssign (mv, u2);
         return u3;
      }

      a.div (u3, v3, q, u3);
      a.subFrom (u2, a.mul (q, v2));

      if (a.is0 (u3))
      {
         destructiveAssign (mv, v2);
         return v3;
      }

      a.div (v3, u3, q, v3);
      a.subFrom (v2, a.mul (q, u2));
   }
}


/**
 * genGcd() - no multiplicator
 */

template<class A>
typename A::type
HIntLib::genGcd (const A& a, const typename A::type& u, const typename A::type& v)
{
   typedef typename A::type T;
   
   if (a.is0 (u))  return v;
   T v3 (a.rem (v, u));

   if (a.is0 (v3))  return u;
   T u3 (a.rem (u, v3));

   for (;;)
   {
      if (a.is0 (u3)) return v3;
      a.reduce (v3, u3);

      if (a.is0 (v3)) return u3;
      a.reduce (u3, v3);
   }
}


/**
 * genIsCoprime()
 */

template<class A>
bool
HIntLib::genIsCoprime (
      const A& a, const typename A::type& u, const typename A::type& v)
{
   typedef typename A::type T;
   
   if (a.is0 (u))  return a.isUnit (v);
   T v3 (a.rem (v, u));

   if (a.is0 (v3))  return a.isUnit (u);
   T u3 (a.rem (u, v3));

   const T* p;

   for (;;)
   {
      if (a.is0 (u3))  { p = &v3; break; }
      a.reduce (v3, u3);

      if (a.is0 (v3))  { p = &u3; break; }
      a.reduce (u3, v3);
   }

   return a.isUnit (*p);
}


/**
 *  powerMod()
 */

template<class A>
typename A::type
HIntLib::powerMod (const A &a, const typename A::type& x, unsigned exponent,
                   const typename A::type &m)
{
   typename A::type result (a.one());
   typename A::type xx (x);

   for (;;)
   {
      if (exponent & 1)
      {
         a.mulBy (result, xx);
         a.reduce (result, m);
      }

      if ((exponent >>= 1) == 0)  return result;

      a.square (xx);
      a.reduce (xx, m);
   }
}

template<class A>
typename A::type
HIntLib::powerMod (
      const A &a, const typename A::type& x, unsigned exp1, unsigned exp2,
                  const typename A::type& m)
{
   unsigned i = 0;
   typename A::type xx (x);

   for (;;)
   {
      typename A::type result (a.one());
      unsigned e = exp1;

      for (;;)
      {
         if (e & 1)
         {
            a.mulBy (result, xx);
            a.reduce (result, m);
         }

         if ((e >>= 1) == 0)  break;

         a.square (xx);
         a.reduce (xx, m);
      }

      if (++i == exp2)  return result;

      destructiveAssign (xx, result);
   }
}

#define HINTLIB_INSTANTIATE_GENGCD(X) \
   template bool genIsCoprime (const X&, const X::type&, const X::type&); \
   template X::type genGcd (const X&, const X::type&, const X::type&); \
   template X::type genGcd \
      (const X&, const X::type&, const X::type&, X::type &); \
   template X::type genGcd \
      (const X&, const X::type&, const X::type&, X::type &, X::type &); \
   template X::type powerMod \
      (const X&, const X::type&, unsigned, const X::type&); \
   template X::type powerMod \
      (const X&, const X::type&, unsigned, unsigned, const X::type&);

