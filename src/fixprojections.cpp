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

#include <iostream>
#include <iomanip>
#include <memory>

#include <HIntLib/tparameter.h>
#include <HIntLib/generatormatrixvirtual.h>
#include <HIntLib/generatormatrixgenrow.h>
#include <HIntLib/minmaxfinder.h>
#include <HIntLib/array.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;
typedef L::GeneratorMatrixGenRow<unsigned char> GM;

using std::cout;
using std::endl;


/**
 *  makeRegular()
 */

void
L::makeRegular (GM& gm, unsigned d)
{
   const unsigned M = gm.getM();

   if (gm.getPrec() < M)
   {
      throw OtherException (
            "makeRegular(): Precision must be at least as large as m!");
   }

   if (d >= gm.getDimension())  throw InvalidDimension (d);

   gm.la().basisSupplement (gm(d), M, M);
}


/**
 *  fixOneDimensionalProjections()
 */

void
L::fixOneDimensionalProjections (GM& gm)
{
   const unsigned M = gm.getM();
   const unsigned DIM = gm.getDimension();

   if (gm.getPrec() < M)
   {
      throw OtherException (
            "fixOneDimensionalProjections(): "
            "Precision must be at least as large as m!");
   }


   for (unsigned d = 0; d < DIM; ++d)  gm.la().basisSupplement (gm(d), M, M);
}


/**
 *  withIdentityMatrix ()
 */

void
L::withIdentityMatrix (GeneratorMatrixGenRow<unsigned char>& gm, unsigned d)
{
   const unsigned DIM = gm.getDimension();
   if (d >= DIM)  throw InvalidDimension (d);

   LinearAlgebra& la = gm.la();

   const unsigned m  = gm.getM();
   const unsigned m2 = m * m;
   const unsigned numRows  = gm.getPrec();
   const unsigned numRows2 = std::min (m, numRows);

   Array<unsigned char> scratch (m * (m + numRows));
   unsigned char* scratch1 = scratch.begin();
   unsigned char* scratch2 = scratch.begin() + m2;

   // Copy first matrix to scrap. Add 0-rows in order to ensure that it is a
   // square matrix

   std::copy (gm(d), gm(d) + m * numRows2, scratch1);
   std::fill (scratch1 + m * numRows2, scratch1 + m2, 0);

   // replace trailing lin.depending vectors by basis vectors, then invert

   la.basisSupplement (scratch1, m, m);
   la.matrixInverse   (scratch1, m);

   // multiply with all matrices

   for (unsigned dd = 0; dd < gm.getDimension(); ++dd)
   {
      std::copy (gm(dd), gm(dd) + m * numRows, scratch2);
      la.matrixMul (scratch2, scratch1, numRows, m, m, gm(dd));
   }
}


/**
 *  withIdentityMatrix2()
 */

void
L::withIdentityMatrix2 (GM& gm, unsigned d)
{
   const unsigned DIM = gm.getDimension();
   if (d >= DIM)  throw InvalidDimension (d);

   const unsigned M = gm.getM();
   if (M != gm.getPrec())  throw FIXME (__FILE__, __LINE__);

   const unsigned M2 = M * M;
   LinearAlgebra& la = gm.la();

   // Construct inverse transformation, tranforming original matrix into E

   if (! la.matrixInverse (gm(d), M))
   {
      throw FIXME (__FILE__, __LINE__);
   }

#if 0
   cout << "Transformation: \n";
   printMatrix (gm(d), M, M);
#endif

   // Transform the other matrices

   Array<unsigned char> temp (M2);

   for (unsigned dd = 0; dd < DIM; ++dd)
   {
      if (dd == d)  continue;

      la.matrixMul (gm(dd), gm(d), M, temp.begin());
      std::copy (temp.begin(), temp.begin() + M2, gm(dd));
   }

#if 0
   cout << "Transformed matrices: \n";
   printMatrix (tm, M * DIM, M);
#endif

   // Replace transformation matrix with identity matrix

   gm.makeIdentityMatrix (d);
}


/**
 *  With Identity Matrix
 */

void
L::WithIdentityMatrix::init (unsigned d, const GeneratorMatrix& src)
{
   GeneratorMatrixGenRow<unsigned char>& gm =
      *new GeneratorMatrixGenRow<unsigned char> (src);
   setMatrix (&gm);

   withIdentityMatrix (gm, d);
}


namespace
{

void
printMatrix (std::ostream& o,
             const unsigned char* m, unsigned numRows, unsigned numCols)
{
   o << '\n';

   for (unsigned row = 0; row < numRows; ++row)
   {
      for (unsigned col = 0; col < numCols; ++col)
      {
         o << std::setw (2) << int (*m++);
      }
      o << '\n';
   }
   o << '\n';
}

void
printVector (std::ostream& o, const unsigned char* i, const unsigned char* end)
{
   o << '(';
   while (i != end)  o << int (*i++);
   o << ')';
}


/**
 *  copyToMatrix()
 *
 *  Copies the first num (or M) vectors of a generator matrix to the given
 *  array.
 */

inline
unsigned char*
copyToMatrix (const GM& gm, unsigned dim, unsigned num, unsigned char* m)
{
   const unsigned char* base = gm(dim);
   return std::copy (base, base + gm.getM() * num, m);
}


/**
 *  AchievableThickness
 *
 *  Records an upper bound on the maximum achievable thickness
 *    - between two coordinates,
 *    - from a given coordinate to all others, and
 *    - for the whole net.
 */

class AchievableThickness
{
public:
   AchievableThickness (const GM&);

   int get() const  { return data0; }
   int get (unsigned d) const  { return data1 [d]; }
   int get (unsigned d1, unsigned d2) const
      { return data2 [(d1 < d2) ? (d1 * DIM + d2) : (d2 * DIM + d1)]; }

   int getMax (unsigned d) const;

   bool set (unsigned d1, unsigned d2, int x);

   void init   (const GM&, unsigned);
   void update (const GM&, int);

   void print() const;

private:
   const unsigned DIM;
   const unsigned M;

   L::Array<unsigned char> m;

   L::Array<int> data2;
   L::Array<int> data1;
   int data0;
};

AchievableThickness::AchievableThickness (const GM& gm)
   : DIM (gm.getDimension()), M (gm.getM()),
     m (M * M),
     data2 (DIM * DIM, M), data1 (DIM, M), data0 (M)
{}

int AchievableThickness::getMax (unsigned d) const
{
   L::MaxFinder<int> mf;
   for (unsigned i = 0; i < DIM; ++i)  if (i != d)  mf << get (i, d);
   return mf.getMaximum();
}

bool AchievableThickness::set (unsigned d1, unsigned d2, int x)
{
   int& r = data2 [(d1 < d2) ? (d1 * DIM + d2) : (d2 * DIM + d1)];

   if (x < r)
   {
      r = x;

      if (x < data1[d1])  data1[d1] = x;
      if (x < data1[d2])  data1[d2] = x;
      if (x < data0)  data0 = x;

      return true;
   }
   else  return false;
}

/**
 *  init()
 *
 *  Determines if there are any thickness-bounds resulting from the first
 *  _numRows_ rows.
 */

void AchievableThickness::init (const GM& gm, unsigned numRows)
{
   const unsigned MAX_THICKNESS = std::min (2u * numRows, M);

   for (unsigned d1 = 1; d1 < DIM; ++d1)
   for (unsigned d2 = 0; d2 < d1;  ++d2)
   {
      for (unsigned k = 1; k <= MAX_THICKNESS; ++k)
      {
         for (unsigned k1 = 0; k1 <= k; ++k1)
         {
            if (k1 > numRows || k - k1 > numRows)  continue;

            unsigned char* p = m.begin();
            p = copyToMatrix (gm, d1, k1,     p);
                copyToMatrix (gm, d2, k - k1, p);

            if (! gm.la().isLinearlyIndependent(m.begin(), k, M))
            {
               set (d1, d2, k - 1);
               goto failed;
            }
         }
      }

      failed: ;
   }
}


/**
 *  update()
 *
 *  Check if there are any new thickness-bounds resulting from the last one of 
 *  _numRows_ rows.
 */

void AchievableThickness::update (const GM& gm, int numRows)
{
   for (unsigned d1 = 1; d1 < DIM; ++d1)
   for (unsigned d2 = 0; d2 < d1;  ++d2)
   {
      {
         if (get(d1,d2) < numRows) continue;

         unsigned char* p = m.begin();
         p = copyToMatrix (gm, d1, numRows,              p);
             copyToMatrix (gm, d2, get(d1,d2) - numRows, p);

         int num =
            gm.la().numLinearlyIndependentVectors (m.begin(), get(d1,d2), M);
         if (num < 2 * numRows)  set (d1, d2, num);
      }

      {
         if (get(d1,d2) < numRows) continue;

         unsigned char* p = m.begin();
         p = copyToMatrix (gm, d2, numRows,              p);
             copyToMatrix (gm, d1, get(d1,d2) - numRows, p);

         int num =
            gm.la().numLinearlyIndependentVectors (m.begin(), get(d1,d2), M);
         if (num < 2 * numRows)  set (d1, d2, num);
      }
   }
}

void AchievableThickness::print() const
{
   cout << "AchievableThickness\n" << get() << '\n';
   cout << "    ";
   for (unsigned i = 0; i < DIM; ++i)  cout << std::setw(3) << i;
   cout << "\n    ";
   for (unsigned i = 0; i < DIM; ++i)  cout << std::setw(3) << get(i);
   cout << "\n    ";
   for (unsigned i = 0; i < DIM; ++i)  cout << std::setw(3) << getMax(i);
   cout << "\n\n";
   for (unsigned i = 0; i < DIM - 1; ++i)
   {
      cout << std::setw (3) << i << ' ';
      for (unsigned j = 0; j <= i; ++j)  cout << "   ";
      for (unsigned j = i + 1; j < DIM; ++j)  cout << std::setw(3) << get(i,j);
      cout << '\n';
   }
}


/**
 *  increment()
 *
 *  Iterate through all linearly independent vectors.
 *
 *  The least significant non-zero digit is always one.
 */

inline
bool increment (unsigned char* begin, unsigned char* end, unsigned base)
{
   if (*begin == 0)  // a trailing zero can always be incremented
   {
      *begin = 1;
   }
   else   // trailing one ...
   {
      unsigned char* i = begin;

      *i = 0;  // ... is reset to zero and triggers a carry
      if (++i == end)  return false;

      // ... carry is propagated by all digits  base-1, which are reset to 0

      while (*i == base - 1)
      {
         *i = 0;
         if (++i == end)  return false;
      }

      // Finaly, we reach a digit which can be incremented.
      // However, this is now the least significant non-zero digits. If it not
      // equal to one, we increment again by one, which resets the least
      // significant digit.

      if (++ *i > 1)  *begin = 1;
   }

   return true;
}


/**
 *  fillFirstRow()
 */

inline
void fillFirstRow (GM& gm)
{
   const unsigned B = gm.getBase();
   const unsigned M = gm.getM();
   const unsigned DIM = gm.getDimension();

   unsigned char* p = gm(0);
   std::fill (p, p + M, 0);
   increment (p, p + M, B);

   for (unsigned d = 1; d < DIM; ++d)
   {
      std::copy (p, p + M, gm(d));
      p = gm(d);
      increment (p, p + M, B);
   }
}


/**
 *  fullCorrect()
 *
 *  Correct all 2-dimensional projections simultaniously
 */

bool
fullCorrect (
      GM& gm, unsigned d, int out, unsigned char* work,
      bool specialAlgo, AchievableThickness& at)
{
   const unsigned DEB = 0;

   // Geometry of the matrix

   const unsigned M   = gm.getM();
   const unsigned M2  = M * M;
   const unsigned DIM = gm.getDimension();
   const unsigned B   = gm.getBase();
   L::LinearAlgebra& la = gm.la();

   if (DEB > 1) cout << endl;
   if (DEB > 0)
   {
      cout << "Doing dimension " << d << "/" << DIM << ", out = " << out;
      if (specialAlgo)  cout << ", special algo";
      cout << endl;
   }

#if 0  // no, we need to make sure that a new l.i. vector is appended!
   if (maxBreadth == 0)
   {
      if (DEB > 1) cout << "breadth > 0 cannot be achieved!" << endl;
      return 0;
   }
#endif
 
   // we need  m^2 (3 + dim) + 2 M  bytes

   unsigned char* const trans1     = work;          // m^2
   unsigned char* const trans2     = trans1 + M2;   // m^2
   unsigned char* const tm         = trans2 + M2;   // m^2 * dim
   unsigned char* const m          = tm + M2 * DIM; // m^2
   unsigned char* const candidate  = m + M2;        // m
   unsigned char* const optimalVec = candidate + M; // m

   // determine breadth resulting from original vector

   int originalBreadth;

   if (specialAlgo)
   {
      // special algo ignores the original vector
      originalBreadth = -1;
   }
   else
   {
      originalBreadth = M;

      for (unsigned d2 = 0; d2 < DIM; ++d2)
      {
         if (d == d2)  continue;

         unsigned char* p = m;
         p = copyToMatrix (gm, d,  out + 1,       p);
             copyToMatrix (gm, d2, M - (out + 1), p);

         int num = la.numLinearlyIndependentVectors (m, M, M);
         if (num < out + 1)
         {
            originalBreadth = -1;
            break;
         }

         int localBreadth = num - (out + 1);
         
         if (localBreadth < originalBreadth)  originalBreadth = localBreadth;
      }

      if (DEB > 1)  cout << "Original vector has breadth " << originalBreadth
                         << endl;

      // If the original value results in a sufficient breadth, we return
      // immediately

      if (originalBreadth >= out + 1)
      {
         if (DEB > 1)  cout << "Inital breadth sufficient!" << endl;
         return false;
      }
   }

   int optimalBreadth = originalBreadth;

   // Build transformation matrix, tranforming E into original matrix

   copyToMatrix (gm, d, out, trans2);
   if (int (la.basisSupplement (trans2, out, M)) < out)
   {
      throw L::InternalError (__FILE__, __LINE__);
   }

   if (DEB > 2)
   {
      cout << "Inverse Transformation: \n";
      printMatrix (cout, trans2, M, M);
   }

   // Construct inverse transformation, tranforming original matrix into E

   std::copy (trans2, trans2 + M2, trans1);
   if (! la.matrixInverse (trans1, M))
   {
      throw L::InternalError (__FILE__, __LINE__);
   }

   if (DEB > 2)
   {
      cout << "Transformation: \n";
      printMatrix (cout, trans1, M, M);
   }

   // Transform all generator matrices

   for (unsigned dd = 0; dd < DIM; ++dd)
   {
      copyToMatrix (gm, dd, M, m);
      la.matrixMul (m, trans1, M, tm + dd * M2);
   }

   if (DEB > 3)
   {
      cout << "Transformed matrices: \n";
      printMatrix (cout, tm, M * DIM, M);
   }

   const int MM = M - out;

   if (specialAlgo)
   {
      // Entries 0,..,M-1 contain num dims for truncated breadth
      // Entires M,..,2M-1 contain num dims for untracted breadth
      L::Array<unsigned> optimalNumDimensions (2 * M, 0u);
      L::Array<unsigned> numDimensions (2 * M);

      L::Array<int> breadths (DIM);
      L::Array<int> optimalBreadths (DIM);

      // try all possible row vectors

      std::fill (candidate, candidate + MM, 0);
      while (increment(candidate, candidate + MM, B))
      {
         // determine maximum mumber of l.i. vectors with each dimension
            
         std::fill (numDimensions.begin(), numDimensions.begin() + 2 * M, 0);

         for (unsigned d2 = 0; d2 < DIM; ++d2)
         {
            if (d == d2)  continue;

            // copy test vector and transformed matrix

            unsigned char* p = m;
            p = std::copy (candidate, candidate + MM, p);
            for (int b = 0; b < MM - 1; ++b)
            {
               const unsigned char* in = tm + d2 * M2 + b * M;
               p = std::copy (in + out, in + M, p);
            }

            // count linearly independent vectors for this dimension

            int localBreadth
               = la.numLinearlyIndependentVectors (m, MM, MM) - 1;
            breadths[d2] = localBreadth;

            // Discard extra precision

            int localBreadthTrunc =
               std::min (localBreadth, std::max (at.get (d, d2) - out - 1, 0));

            // update combined result

            for (int i = 1; i <= localBreadthTrunc; ++i) ++numDimensions[i];
            for (int i = 1; i <= localBreadth; ++i)      ++numDimensions[M + i];
         }

         // is the new vector superior to previous results?

         for (unsigned i = 1; i < 2 * M; ++i)
         {
            if (numDimensions[i] > optimalNumDimensions[i])  break;
            if (numDimensions[i] < optimalNumDimensions[i])  goto vectorFailed;
         }

//         optimalBreadth = breadth;

         // we got a new optimal breadth!

         std::copy (candidate, candidate + MM, optimalVec);
         std::copy (numDimensions.begin(), numDimensions.begin() + 2 * M,
                    optimalNumDimensions.begin());
         std::copy (breadths.begin(), breadths.begin() + DIM,
                    optimalBreadths.begin());

         if (DEB > 1)
         {
            cout << "New optimal num dims:";
            for (unsigned i = 1; i < 2 * M; ++i)
            {
               cout << ' ' << optimalNumDimensions[i];
            }
            cout << ", breadth = " << optimalBreadth << endl;
         }

   vectorFailed: ;
      }

      if (DEB > 1)  cout << "Breadths:";
      for (unsigned i = 0; i < DIM; ++i)
      {
         if (i == d)  continue;
         if (DEB > 1)  cout << ' ' << optimalBreadths[i];
         at.set (d, i, optimalBreadths[i] + out + 1);
      }
      if (DEB > 1) cout << endl;

      optimalBreadth = out;
   }
   else   // normal algo
   {
      for (int breadth = optimalBreadth + 1; breadth <= out + 1; ++breadth)
      {
         if (DEB > 1)
         {
            cout << " Trying breadth " << breadth << " ... " << std::flush;
         }

         // try all possible row vectors

         std::fill (candidate, candidate + MM, 0);
         while (increment(candidate, candidate + MM, B))
         {
            // Print candidate vector for debugging

            if (DEB > 3)
            {
               cout <<"     checking ";
               printVector (cout, candidate, candidate + MM);
               cout << "  ";
            }

            // check linear independence with each matrix
            
            for (unsigned d2 = 0; d2 < DIM; ++d2)
            {
               if (d == d2)  continue;

               unsigned char* p = m;
               p = std::copy (candidate, candidate + MM, p);
               for (int b = 0; b < breadth; ++b)
               {
                  const unsigned char* in = tm + d2 * M2 + b * M;
                  p = std::copy (in + out, in + M, p);
               }

               if (! la.isLinearlyIndependent(m, breadth + 1, MM))
               {
                  if (DEB > 3)  cout << 'x' << endl;
                  goto vectorFailed1;
               }
               if (DEB > 3)  cout << '.';
            }

            // vector passed all tests

            if (DEB > 3)  cout << "  ok, done!" << endl;

            std::copy (candidate, candidate + MM, optimalVec);
            optimalBreadth = breadth;
            break;

            vectorFailed1: ;
         }

         // if we have not found a vector, we quit the search

         if (optimalBreadth != breadth)
         {
            if (DEB > 1) cout << "failed" << endl;
            break;
         }
         else
         {
            if (DEB > 1) cout << "success" << endl;
         }
      }
   }

   // Do we have to use the original vector?

   if (DEB > 1)  cout << "Final breadth: " << optimalBreadth << endl;

   if (optimalBreadth <= originalBreadth)
   {
      if (DEB > 1) cout << "Original vector is all right!" << endl;

      return false;
   }
   else
   {
      if (DEB > 1)
      {
         cout << "Transformed output vector " << out << " is ";
         printVector (cout, &optimalVec[0], &optimalVec[MM]);
         cout << endl;
      }

      // transform back vector

      la.matrixMul (optimalVec, trans2 + out * M, 1, MM, M, m);

      // copy new vector back to generator matrix

      if (DEB > 1)
      {
         cout << "Output vector " << out << " is ";
         printVector (cout, &m[0], &m[M]);
         cout << endl;
      }

      for (unsigned i = 0; i < M; ++i)  gm.setd (d, i, out, m[i]);

      return true;
   }
}

}  // anonimous namespace


/**
 *  fixTwoDimensionalProjections()
 */

void L::fixTwoDimensionalProjections (GM& gm)
{
   const unsigned DEB = 0;

   const unsigned M = gm.getM();
   const unsigned DIM = gm.getDimension();

   if (gm.getPrec() < M)
   {
      throw OtherException (
            "fixTwoDimensionalProjections(): "
            "Precision must be at least as large as m!");
   }

   Array<unsigned char> work ((3 + DIM) * M * M + 2 * M);

   int initialOut = std::max (
      std::max (M - tParameterMax3DimProjection(gm), 2u) - 2,  // what is taboo
      (M - tParameterMax2DimProjection(gm)) / 2     // what is alright already
   );

   AchievableThickness at (gm);
   at.init (gm, initialOut);
   if (DEB > 2)  at.print();

   Array<bool> dimensionDone (DIM);

   if (initialOut == 0)
   {
      fillFirstRow (gm);
      initialOut = 1;
   }

   for (int out = initialOut; out < int (M); ++out)
   {
      if (DEB > 1)  cout << "\nout = " << out << endl;

      fixOneDimensionalProjections (gm);

      std::fill (dimensionDone.begin(), dimensionDone.begin() + DIM, false);
      bool changes;

      do
      {
         changes = false;

         for (unsigned d = 0; d < DIM; ++d)
         {
            if (dimensionDone [d])  continue;
            if (DEB > 1)  cout << d << ": ";

            // Check if the coordinate is independent

            if (at.getMax (d) <= 2 * out + 1)
            {
               if (DEB > 1)  cout << '^';
               fullCorrect (gm, d, out, work.begin(), true, at);
               dimensionDone [d] = true;
            }
            else   // non-independent coordinate
            {
               bool c = fullCorrect (gm, d, out, work.begin(), false, at);
               changes |= c;
               if (DEB > 1)  cout << (c ? '*' : '-');
            }
            if (DEB > 1)  cout << endl;

            if (DEB > 3)  gm.print (cout);
         }
      }
      while (changes);

      at.update (gm, out + 1);
      if (DEB > 2)  at.print();
   }

   fixOneDimensionalProjections (gm);
}



