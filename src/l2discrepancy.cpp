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

#include <HIntLib/discrepancy.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation "discrepancy.h"
#endif

#include <HIntLib/array.h>
#include <HIntLib/hlmath.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;

using std::max;

using L::real;
using L::sqr;
using L::Array;

// #define HINTLIB_CHECK_CORRECTNESS

#ifdef HINTLIB_CHECK_CORRECTNESS
#include <iostream>
#include <iomanip>
using namespace std;
#endif


/**
 *  l2Discrepancy()
 *
 *  Calculate the L_2-discrepancy of a point set based on Warnock's forumla.
 *
 *  Runtime is O(n^2).
 */

real
L::l2DiscrepancyWarnock (const real* points, unsigned dim, Index n)
{
   checkDimensionNotZero (dim);

   // Term 2 and diagonal part of term 3

   real term2 = 0;
   real term3Diag = 0;

   const real* const lastPoint = points + n * dim;

   for (const real* point = points; point != lastPoint; )
   {
      const real* const lastCoord = point + dim;

      real product2 = 1 - sqr (*point);
      real product3 = 1 - *point;

      while (++point != lastCoord)
      {
         product2 *= 1 - sqr (*point);
         product3 *= 1 - *point;
      }

      // cout << "product2 = " << product2 << '\n';

      term2     += product2;
      term3Diag += product3;
   }

   // One triangle of term 3

   real term3Triangle = 0;

   for (const real* point1 = points + dim; point1 != lastPoint; point1 += dim)
   {
      const real* const lastCoord1 = point1 + dim;

      for (const real* point2 = points;       point2 != point1; )
      {
         real product = 1 - max (*point1, *point2++);

         for (const real* p1 = point1 + 1; p1 != lastCoord1; ++p1, ++point2)
         {
            product *= 1 - max (*p1, *point2);
         }

         // cout << "product3 = " << product << '\n';

         term3Triangle += product;
      }
   }

   // return result

   return HINTLIB_MN sqrt(
         1 / HINTLIB_MN pow (real(3), int(dim))
        - term2 / (HINTLIB_MN pow (real(2), int(dim) - 1) * n)
        + (term3Diag + 2 * term3Triangle) / sqr(real(n)));
}

namespace
{
   struct L2Data
   {
      const real* point;  // pointer to the array of coordinates
      real weight;        // weight (initially 1/n)
   };


   /**
    *  partition()
    *
    *  Swaps elements in [begin,end)  such that
    *    - all elements in [begin,center) are < split, and
    *    - all elements in [center,end) are >= split.
    *
    *  center is returned.
    */

   inline
   L2Data* partition (L2Data* begin, L2Data* end, int dim, real split)
   {
// #define HINTLIB_DEBUG

#if defined HINTLIB_DEBUG || defined HINTLIB_CHECK_CORRECTNESS
      L2Data *begin2 = begin, *end2 = end;
#endif
#ifdef HINTLIB_DEBUG
      cout << "------------- Before --------------\n";
      for (L2Data* p = begin2; p != end2; ++p)
      {
         cout << setw(16) << p->point[dim] << setw(16) << p->weight << '\n';
      }
#endif

      L2Data* center = begin + (end - begin) / 2;
      --end;

      for (;;)
      {
         // skip elements that are in the proper part of the list

         while (begin->point[dim] <  split) if (++begin == end)  goto done;
         while (end  ->point[dim] >= split) if (--end == begin)  goto done;

         // exchange offending elements

         std::swap (*begin, *end);

         // move beyond them

         if (++begin >= --end)  break;
      }
   done:

      // Determe to which part the last record belongs

      if (begin->point[dim] < split)  ++begin;

      // If there are equal values, try to move towards the center

      while (begin > center && (begin - 1)->point[dim] >= split)  --begin;
      while (begin < center && (begin + 1)->point[dim] <  split)  ++begin;

#ifdef HINTLIB_DEBUG
      cout << "------------- After ---------------\n";
      for (L2Data* p = begin2; p != end2; ++p)
      {
         if (begin == p)   cout << "        ---- " << split << '\n';
         if (center == p)  cout << "        -\n";
         cout << setw(16) << p->point[dim] << setw(16) << p->weight << '\n';
      }
      if (begin == end2)  cout << "        ---- " << split << '\n';
      cout << "-----------------------------------\n";
#endif

#ifdef HINTLIB_CHECK_CORRECTNESS
      for (L2Data* p = begin2; p != center; ++p)
      {
         if (p->point[dim] >= split)  {  cout << "Error!\n\n"; exit(1); }
      }
      for (L2Data* p = center; p != end2; ++p)
      {
         if (p->point[dim] < split)  {  cout << "Error!\n\n"; exit(1); }
      }
#endif
#undef HINTLIB_DEBUG

      return begin;
   }


   /**
    *  directCalculation()
    */

   inline
   real
   directCalculation (
      L2Data* beginA, L2Data* endA, L2Data* beginB, L2Data* endB, int dim)
   {
      real outerSum = 0;

      for ( ; beginA != endA; ++beginA)
      {
         real innerSum = 0;
         
         for (L2Data* pb = beginB; pb != endB; ++pb)
         {
            real product = pb->weight;

            for (int d = 0; d <= dim; ++d)
            {
               product *= 1 - max (beginA->point[d], pb->point[d]);
            }
            
            innerSum += product;
         }

         outerSum += innerSum * beginA->weight;
      }

      return outerSum;
   }


   /**
    *  directCalculationSymmetric()
    */

   inline
   real
   directCalculationSymmetric (L2Data* begin, L2Data* end, int dim)
   {
      real sum = 0;

      for (L2Data* p1 = begin ; p1 != end; ++p1)
      {
         for (L2Data* p2 = p1 + 1; p2 != end; ++p2)
         {
            real product = 1;

            for (int d = 0; d <= dim; ++d)
            {
               product *= 1 - max (p1->point[d], p2->point[d]);
            }
            
            sum += product;
         }
      }

      return sum;
   }

   const unsigned THRESHOLD = 12;


   /**
    * printList()
    */

#ifdef HINTLIB_CHECK_CORRECTNESS
   void
   printList (const L2Data* begin, const L2Data* end, int dim)
   {
      cout << "Dim = " << dim << '\n';
      for (const L2Data* p = begin; p != end; ++p)
      {
         cout << p->weight << '\t';
         for (int d = 0; d <= dim; ++d)
         {
            cout << p->point[d] << ' ';
         }
         cout << '\n';
      }
   }
#endif


   /**
    *  heinrichRecursion()
    */

   real
   heinrichRecursion (
         L2Data* beginA, L2Data* endA, L2Data* beginB, L2Data* endB,
         int dim, real lb, real ub)
   {
      // const int DEB = 1;

      const unsigned na = endA - beginA;
      const unsigned nb = endB - beginB;

#ifdef HINTLIB_CHECK_CORRECTNESS
      if (na == 0 || nb == 0)
      {
         cerr << "na or nb are 0!!!\n\n";
         exit(1);
      }
#endif

#if 0
      if (DEB > 0)
      {
         cout << "na = " << na << ", nb = " << nb << ", dim = " << dim
            << ", values in [" << lb << ',' << ub << "].  ";
      }
#endif

      // Case 2: There are no coordinates/dimensions left
      // This is the improved basic step that is responsible for the better
      // asymptotic performance

      if (dim == -1)
      {
#ifdef HINTLIB_CHECK_CORRECTNESS
         const real correctResult = directCalculation
            (beginA, endA, beginB, endB, dim);
#endif
         real sumA = 0;
         while (beginA != endA)  sumA += (beginA++)->weight;
         real sumB = 0;
         while (beginB != endB)  sumB += (beginB++)->weight;

         const real result = sumA * sumB;

#ifdef HINTLIB_CHECK_CORRECTNESS
         if (result != correctResult)
         {
            cout << "Error calculating base optimization!\n"
                 << "Correct: " << correctResult << ", but " << result << "\n";
            printList (beginA, endA, dim);
            printList (beginB, endB, dim);
            exit(1);
         }
#endif
         // if (DEB > 0)  cout << "Dim = 0.  Returning with " << result << '\n';
         return result;
      }

      // Case 3a: List A contains only a single point

#if 1
      if (na == 1)
      {
#ifdef HINTLIB_CHECK_CORRECTNESS
         const real correctSum= directCalculation
            (beginA, endA, beginB, endB, dim);
#endif

         // cout << "A\n";
         real sum = 0;
         for ( ; beginB != endB; ++beginB)
         {
            real product = beginB->weight;

            for (int d = 0; d <= dim; ++d)
            {
               product *= 1 - max (beginA->point[d], beginB->point[d]);
            }

            sum += product;
         }

         sum *= beginA->weight;

#if 0
         if (DEB > 0)
         {
            cout << "List A contains only one element. " << sum << '\n';
         }
#endif

#ifdef HINTLIB_CHECK_CORRECTNESS
         if (sum != correctSum)
         {
            cout << "Error calculating list-A-trivial term!\n"
                 << "Correct: " << correctSum << ", but " << sum << "\n";
            printList (beginA, endA, dim);
            printList (beginB, endB, dim);
            exit(1);
         }
#endif

         return sum;
      }

      // Case 3b: List B contains only a single point

      if (nb == 1)
      {
#ifdef HINTLIB_CHECK_CORRECTNESS
         const real correctSum= directCalculation
            (beginA, endA, beginB, endB, dim);
#endif

         // cout << "B\n";
         real sum = 0;
         for ( ; beginA != endA; ++beginA)
         {
            real product = beginA->weight;

            for (int d = 0; d <= dim; ++d)
            {
               product *= 1 - max (beginA->point[d], beginB->point[d]);
            }

            sum += product;
         }

         sum *= beginB->weight;

#if 0
         if (DEB > 0)
         {
            cout << "List B contains only one element. " << sum << '\n';
         }
#endif
#ifdef HINTLIB_CHECK_CORRECTNESS
         if (sum != correctSum)
         {
            cout << "Error calculating list-B-trivial term!\n"
                 << "Correct: " << correctSum << ", but " << sum << "\n";
            printList (beginA, endA, dim);
            printList (beginB, endB, dim);
            exit(1);
         }
#endif

         return sum;
      }
#endif

      // If there are too little elements left, calculate the product directly

#if 1
      if (std::min (na, nb) <= (THRESHOLD << dim))
      {
         real result
            = directCalculation (beginA, endA, beginB, endB, dim);

         // if (DEB > 0) cout << "Small list. " << result << '\n';

         return result;
      }
#endif

      // Case 4: Recursion
      // We can assume that na >=2, nb >= 2, and at least one coordinate left

      // if (DEB > 0)  cout << "Recursion..." << '\n';
      const real mid = (ub + lb) * 0.5;

      // Sort and split both lists

      L2Data* splitA = partition (beginA, endA, dim, mid); 
      L2Data* splitB = partition (beginB, endB, dim, mid); 
      const unsigned restA = endA - splitA;
      const unsigned restB = endB - splitB;

      // Get results for the four sub-problems
    
#ifdef HINTLIB_CHECK_CORRECTNESS
      const real correctTermLL = directCalculation
         (beginA, splitA, beginB, splitB, dim);
#endif
      const real termLL = (splitA == beginA || splitB == beginB) ? 0 :
         heinrichRecursion (beginA, splitA, beginB, splitB, dim, lb, mid);
#ifdef HINTLIB_CHECK_CORRECTNESS
      const real correctTermLL2 = directCalculation
         (beginA, splitA, beginB, splitB, dim);
      if (termLL != correctTermLL || termLL != correctTermLL2)
      {
         cout << "Error calculating LL term!\n"
              << "Correct: " << correctTermLL << ", but " << termLL
              << ", afterwards: " << correctTermLL2 << "\n";
         printList (beginA, splitA, dim);
         printList (beginA, splitB, dim);
         exit(1);
      }
#endif

      if (restA == 0 && restB == 0)  return termLL;

#ifdef HINTLIB_CHECK_CORRECTNESS
      const real correctTermHH = directCalculation
         (splitA, endA, splitB, endB, dim);
#endif
      const real termHH = (restA == 0 || restB == 0) ? 0 :
         heinrichRecursion (splitA, endA, splitB, endB, dim, mid, ub);
#ifdef HINTLIB_CHECK_CORRECTNESS
      const real correctTermHH2 = directCalculation
         (splitA, endA, splitB, endB, dim);
      if (termHH != correctTermHH || termHH != correctTermHH2)
      {
         cout << "Error calculating HH term!\n"
              << "Correct: " << correctTermHH << ", but " << termHH
              << ", afterwards: " << correctTermHH2 << "\n";
         printList (splitA, endA, dim);
         printList (splitB, endB, dim);
         exit(1);
      }
#endif

      if (splitA == beginA && splitB == beginB)  return termHH;

      Array<L2Data> newList (max (restA, restB));

      real termHL = 0, termLH = 0;

      if (restA != 0 && splitB != beginB)
      {
         for (L2Data* src = splitA, *dest = newList.begin(); src != endA;
              ++src, ++dest)
         {
            dest->point  = src->point;
            dest->weight = src->weight * (1 - src->point[dim]);
         }
#ifdef HINTLIB_CHECK_CORRECTNESS
         const real correctTermHL = directCalculation
            (newList.begin(), newList.begin() + restA, beginB, splitB, dim - 1);
#endif
         termHL = heinrichRecursion (
               newList.begin(), newList.begin() + restA, beginB, splitB,
               dim - 1, 0.0, 1.0);
#ifdef HINTLIB_CHECK_CORRECTNESS
         const real correctTermHL2 = directCalculation
            (newList.begin(), newList.begin() + restA, beginB, splitB, dim - 1);
         if (termHL != correctTermHL || termHL != correctTermHL2)
         {
            cout << "Error calculating HL term!\n"
                 << "Correct: " << correctTermHL << ", but " << termHL
                 << ", afterwards: " << correctTermHL2 << "\n";
            printList (newList.begin(), newList.begin() + restA, dim);
            printList (beginB, splitB, dim);
            exit(1);
         }
#endif
      }

      if (restB != 0 && splitA != beginA)
      {
         for (L2Data* src = splitB, *dest = newList.begin(); src != endB;
              ++src, ++dest)
         {
            dest->point  = src->point;
            dest->weight = src->weight * (1 - src->point[dim]);
         }

#ifdef HINTLIB_CHECK_CORRECTNESS
         const real correctTermLH = directCalculation
            (beginA, splitA, newList.begin(), newList.begin() + restB, dim - 1);
#endif
         termLH = heinrichRecursion (
               beginA, splitA, newList.begin(), newList.begin() + restB,
               dim - 1, 0.0, 1.0);
#ifdef HINTLIB_CHECK_CORRECTNESS
         const real correctTermLH2 = directCalculation
            (beginA, splitA, newList.begin(), newList.begin() + restB, dim - 1);
         if (termLH != correctTermLH || termLH != correctTermLH2)
         {
            cout << "Error calculating LH term!\n"
                 << "Correct: " << correctTermLH << ", but " << termLH
                 << ", afterwards: " << correctTermLH2 << "\n";
            printList (beginA, splitA, dim);
            printList (newList.begin(), newList.begin() + restB, dim);
            exit(1);
         }
#endif
      }

#if 0
      if (DEB)
      {
         cout << "Returning with " << termLL << " + " << termHH << " + "
              << termLH << " + " << termHL << " = "
              << termLL + termHH + termLH + termHL << "\n";
      }
#endif

      return termLL + termHH + termLH + termHL;
   }


   /**
    *  heinrichRecursionSymmetric()
    */

   real
   heinrichRecursionSymmetric (
         L2Data* begin, L2Data* end, int dim, real lb, real ub)
   {
      // const int DEB = 0;

      const unsigned n = end - begin;

#ifdef HINTLIB_CHECK_CORRECTNESS
      if (n <= 1)
      {
         cout << "n is 0 or 1!!!\n\n";
         exit(1);
      }
#endif

#if 0
      if (DEB > 0)
      {
         cout << "n = " << n << ", dim = " << dim
            << ", values in [" << lb << ',' << ub << "].  ";
      }
#endif

      // If there are too little elements left, calculate the product directly

      if (n <= (THRESHOLD << dim))
      {
         real result = directCalculationSymmetric (begin, end, dim);

         // if (DEB > 0) cout << "Small list. " << result << '\n';

         return result;
      }

      // Case 4: Recursion
      // We can assume that na >=2, nb >= 2, and at least one coordinate left

      // if (DEB > 0)  cout << "Recursion..." << '\n';
      const real mid = (ub + lb) * 0.5;

      // Sort and split both lists

      L2Data* split = partition (begin, end, dim, mid); 
      const unsigned n1 = split - begin;
      const unsigned n2 = end - split;

      // Get results for the four sub-problems
    
      const real termLL = (n1 <= 1) ? 0 :
         heinrichRecursionSymmetric (begin, split, dim, lb, mid);

      if (n2 == 0)  return termLL;

      const real termHH = (n2 <= 1) ? 0 :
         heinrichRecursionSymmetric (split, end, dim, mid, ub);

      if (n1 == 0)  return termHH;

      Array<L2Data> newList (n2);

      for (L2Data* src = split, *dest = newList.begin();
           src != end; ++src, ++dest)
      {
         dest->point  = src->point;
         dest->weight = src->weight * (1 - src->point[dim]);
      }

      const real termHL = heinrichRecursion (
            newList.begin(), newList.begin() + n2, begin, split, dim - 1, 0, 1);

      return termLL + termHH + termHL;
   }

} // anonimous namespace


/**
 *  l2DiscrepancyHeinrich()
 */

real
L::l2DiscrepancyHeinrich (const real* points, unsigned dim, Index n)
{
   checkDimensionNotZero (dim);

   // Term 2 and diagonal part of term 3

   real term2 = 0;
   real term3Diag = 0;

   const real* const lastPoint = points + n * dim;

   for (const real* point = points; point != lastPoint; )
   {
      const real* const lastCoord = point + dim;

      real product2 = 1 - sqr (*point);
      real product3 = 1 - *point;

      while (++point != lastCoord)
      {
         product2 *= 1 - sqr (*point);
         product3 *= 1 - *point;
      }

      // cout << "product2 = " << product2 << '\n';

      term2     += product2;
      term3Diag += product3;
   }

   // Calculate one triangle of term 3 using Heinrich's recursive method

   Array<L2Data> list (n);

   L2Data* dest = list.begin();
   for (const real* src = points; src != lastPoint; src += dim, ++dest)
   {
      dest->point  = src;
      dest->weight = 1;
   }

   const real term3
      = heinrichRecursionSymmetric (list.begin(), dest, dim - 1, 0.0, 1.0);

   // return result

   return HINTLIB_MN sqrt (
         1 / HINTLIB_MN pow (real(3), int(dim))
        - term2 / (HINTLIB_MN pow (real(2), int(dim) - 1) * n)
        + (2 * term3 + term3Diag) / sqr(real(n)));
}


/**
 *  l2Discrepancy()
 *
 *  Calculate the L_2-discrepancy of a point set.
 *
 *  If the number of points is small (compared to 2^dim), then we use Warnock's
 *  forumla. Otherwise, Heinrich's algorithm is used.
 */

real
L::l2Discrepancy (const real* points, unsigned dim, Index n)
{
   if (dim + 6 >= 32 || n <= 80 || (1u << (dim + 2)) >= n)
   {
      return l2DiscrepancyWarnock (points, dim, n);
   }
   else
   {
      return l2DiscrepancyHeinrich (points, dim, n);
   }
}


