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

#ifdef __GNUG__
#pragma implementation
#endif

#include <HIntLib/shiftnet.h>

#include <HIntLib/generatormatrix.h>
#include <HIntLib/exception.h>

namespace L = HIntLib;

using L::u64;

namespace
{
   typedef const u64 cu64;

   struct SNData
   {
      unsigned length;
      const u64* data;
   };

   /**
    * Lists of row vectors for the first matrix.
    *
    * All other matrices are created by appropriately shifting the first matrix
    * (therefore ShiftNet).
    */

   /*
    * base   b = 2
    *
    * The table for b = 2 was originally calculated by Wolfgang Schmid, see
    *
    * [1] Wolfgnag Ch. Schmid. Shift-Nets: A New Class of Binary Digital
    *     (t,m,s)-Nets. In Harald Niederreiter et al, editors, Monte Carlo and
    *     Quasi Monte Carlo Methods 1996, volume 127 of Lecture Notes in
    *     Statistics. Springer, 1997.
    *
    * We use the numbers appearing in this paper with the following exceptions:
    * 
    * - For m=21, the second-to-last number is 6 instead of 62. This results in
    *   the same t value (t=11), but is the lexicographically first solution
    *   with this t value (as are all other entries).
    *
    * - For m >= 23, superior results were calculated by Rudolf Schürer.
    *
    * Optimality is prooven up to  m = 21.
    */

   cu64 sn2_1 [] = {1};
   cu64 sn2_2 [] = {1,    2};
   cu64 sn2_3 [] = {1,    6,      2};
   cu64 sn2_4 [] = {1,    6,      2};
   cu64 sn2_5 [] = {1,   14,     18,     2};
   cu64 sn2_6 [] = {1,   14,     18,     2};
   cu64 sn2_7 [] = {1,   46,     26,     6,    2};
   cu64 sn2_8 [] = {1,   46,     26,     6,    2};
   cu64 sn2_9 [] = {1,   94,    178,    38,   10,    2};
   cu64 sn2_10[] = {1,   94,    170,    38,   10,    2};
   cu64 sn2_11[] = {1,  492,   1638,    86,   14,   18,   2};
   cu64 sn2_12[] = {1,  378,    590,   158,   42,    6,   2};
   cu64 sn2_13[] = {1,  366,    572,   154,   14,   18,   2};
   cu64 sn2_14[] = {1,  734,  10652,   372,   46,   22,  10,   2};
   cu64 sn2_15[] = {1,  734,   3308,   372,   46,   22,  10,   2};
   cu64 sn2_16[] = {1,  758,   3132,   414,   46,   22,  10,   2};
   cu64 sn2_17[] = {1, 3562,  29276,   870,   62,  198,  74,   6,  2};
   cu64 sn2_18[] = {1, 3450,   9814,  2020,  458,   54,  26,   6,  2};
   cu64 sn2_19[] = {1, 3562,   8566,  1884,  236,   30,  38,  10,  2};
   cu64 sn2_20[] = {1,10094, 301178, 17650, 2300,  314,  78,  22, 10, 2};
   cu64 sn2_21[] = {1,10106,  54364,  3974, 4214,  158, 294,  42,  6, 2};
   cu64 sn2_22[] = {1,10198,  54476, 18296, 2674,  124, 150,  26,  6, 2};
   cu64 sn2_23[] = {1,40302, 468876, 69080, 1852,  222,2378,  46, 22,10, 2};
   cu64 sn2_24[] = {1,40302, 213628, 87254, 4978,  238,1640,  30, 38,10, 2};
   cu64 sn2_25[] = {1,40314, 158318, 74218, 5302,  606, 314, 156, 14,20, 2};
   cu64 sn2_26[] = {1,85470,7399828,177330,16084,  702,1388,2126,114,22,10,2};
   cu64 sn2_27[] = {1,81366,1969322,164826, 5950,16630,2382, 604, 46,22,10,2};
   cu64 sn2_28[] = {1,81366,1350324,548382,17334, 1756,2222, 286,102,42, 6,2};
   cu64 sn2_29[] = {1,80630,1443834,449758,11678,16762, 686,1078, 78,26, 6,2};

   // We use the following vector (with  m - t = 10) for m = 29,...,39

   cu64 sn2_x [] = {1,10094,  51322, 22822,  190,  626, 270,  22, 10, 2};

   const SNData snets2 [] =
   {
      { 0, 0 },
      { sizeof (sn2_1)  / sizeof (u64), sn2_1  },
      { sizeof (sn2_2)  / sizeof (u64), sn2_2  },
      { sizeof (sn2_3)  / sizeof (u64), sn2_3  },
      { sizeof (sn2_4)  / sizeof (u64), sn2_4  },
      { sizeof (sn2_5)  / sizeof (u64), sn2_5  },
      { sizeof (sn2_6)  / sizeof (u64), sn2_6  },
      { sizeof (sn2_7)  / sizeof (u64), sn2_7  },
      { sizeof (sn2_8)  / sizeof (u64), sn2_8  },
      { sizeof (sn2_9)  / sizeof (u64), sn2_9  },
      { sizeof (sn2_10) / sizeof (u64), sn2_10 },
      { sizeof (sn2_11) / sizeof (u64), sn2_11 },
      { sizeof (sn2_12) / sizeof (u64), sn2_12 },
      { sizeof (sn2_13) / sizeof (u64), sn2_13 },
      { sizeof (sn2_14) / sizeof (u64), sn2_14 },
      { sizeof (sn2_15) / sizeof (u64), sn2_15 },
      { sizeof (sn2_16) / sizeof (u64), sn2_16 },
      { sizeof (sn2_17) / sizeof (u64), sn2_17 },
      { sizeof (sn2_18) / sizeof (u64), sn2_18 },
      { sizeof (sn2_19) / sizeof (u64), sn2_19 },
      { sizeof (sn2_20) / sizeof (u64), sn2_20 },
      { sizeof (sn2_21) / sizeof (u64), sn2_21 },
      { sizeof (sn2_22) / sizeof (u64), sn2_22 },
      { sizeof (sn2_23) / sizeof (u64), sn2_23 },
      { sizeof (sn2_24) / sizeof (u64), sn2_24 },
      { sizeof (sn2_25) / sizeof (u64), sn2_25 },
      { sizeof (sn2_26) / sizeof (u64), sn2_26 },
      { sizeof (sn2_27) / sizeof (u64), sn2_27 },
      { sizeof (sn2_28) / sizeof (u64), sn2_28 },
      { sizeof (sn2_29) / sizeof (u64), sn2_29 },
      { sizeof (sn2_x)  / sizeof (u64), sn2_x  },  // m = 30
      { sizeof (sn2_x)  / sizeof (u64), sn2_x  },  // m = 31
      { sizeof (sn2_x)  / sizeof (u64), sn2_x  },  // m = 32
      { sizeof (sn2_x)  / sizeof (u64), sn2_x  },  // m = 33
      { sizeof (sn2_x)  / sizeof (u64), sn2_x  },  // m = 34
      { sizeof (sn2_x)  / sizeof (u64), sn2_x  },  // m = 35
      { sizeof (sn2_x)  / sizeof (u64), sn2_x  },  // m = 36
      { sizeof (sn2_x)  / sizeof (u64), sn2_x  },  // m = 37
      { sizeof (sn2_x)  / sizeof (u64), sn2_x  },  // m = 38
      { sizeof (sn2_x)  / sizeof (u64), sn2_x  },  // m = 39
   };


   /*
    *  bases   b > 2
    *
    *  For all bases b > 2, the parameters for optimal shift nets were
    *  calculated by Rudolf Schürer based on work by Wolfgang Schmid.
    *
    *  See an upcoming paper by W.Ch.Schmid and R.Schürer.
    */

   // b = 3
   //
   // Optimality prooven up to  m = 16

   cu64 sn3_1 [] = {1};
   cu64 sn3_2 [] = {1,      3};
   cu64 sn3_3 [] = {1,     12,       3};
   cu64 sn3_4 [] = {1,     12,       3};
   cu64 sn3_5 [] = {1,     39,      21,       3};
   cu64 sn3_6 [] = {1,    129,     264,      12,     3};
   cu64 sn3_7 [] = {1,    480,    1602,     111,    21,    3};
   cu64 sn3_8 [] = {1,    372,     822,      39,    21,    3};
   cu64 sn3_9 [] = {1,   3357,    7824,     120,   255,   21,    3};
   cu64 sn3_10[] = {1,   2577,    2061,     201,   273,   21,    3};
   cu64 sn3_11[] = {1,   8391,   23034,    2856,   120,  255,   21,   3};
   cu64 sn3_12[] = {1,   7680,   24660,    2388,   426,   66,   21,   3};
   cu64 sn3_13[] = {1,  32205,  657297,   62004,   885,  282,   93,  21,  3};
   cu64 sn3_14[] = {1,  25509,  134274,  177600,  1197,  120,  309,  12,  3};
   cu64 sn3_15[] = {1, 280578, 5344227,  557100,  1452, 2388,  309, 102, 12, 3};
   cu64 sn3_16[] = {1, 204762,  659883,   37101,  6222, 7509,  468,  48, 12, 3};
   cu64 sn3_17[] = {1, 204591,  562188,   19236, 61677,  372,  768,  93, 21, 3};
   cu64 sn3_18[] = {1, 615666, 6805764, 1651467, 15222,24186, 6300, 282,102,12, 3};
   cu64 sn3_19[] = {1, 614289, 5110500, 1668576, 37533, 3036,  399,6600, 93,21, 3};
   cu64 sn3_20[] = {1, 614289, 2122158, 4796073, 28965, 1416, 2655, 354, 93,21, 3};
   cu64 sn3_21[] = {1,1842582,24946203,43316163,114276,27444,10074, 696,831,39,21,3};

   const SNData snets3 [] =
   {
      { 0, 0 },
      { sizeof (sn3_1)  / sizeof (u64), sn3_1  },
      { sizeof (sn3_2)  / sizeof (u64), sn3_2  },
      { sizeof (sn3_3)  / sizeof (u64), sn3_3  },
      { sizeof (sn3_4)  / sizeof (u64), sn3_4  },
      { sizeof (sn3_5)  / sizeof (u64), sn3_5  },
      { sizeof (sn3_6)  / sizeof (u64), sn3_6  },
      { sizeof (sn3_7)  / sizeof (u64), sn3_7  },
      { sizeof (sn3_8)  / sizeof (u64), sn3_8  },
      { sizeof (sn3_9)  / sizeof (u64), sn3_9  },
      { sizeof (sn3_10) / sizeof (u64), sn3_10 },
      { sizeof (sn3_11) / sizeof (u64), sn3_11 },
      { sizeof (sn3_12) / sizeof (u64), sn3_12 },
      { sizeof (sn3_13) / sizeof (u64), sn3_13 },
      { sizeof (sn3_14) / sizeof (u64), sn3_14 },
      { sizeof (sn3_15) / sizeof (u64), sn3_15 },
      { sizeof (sn3_16) / sizeof (u64), sn3_16 },
      { sizeof (sn3_17) / sizeof (u64), sn3_17 },
      { sizeof (sn3_18) / sizeof (u64), sn3_18 },
      { sizeof (sn3_19) / sizeof (u64), sn3_19 },
      { sizeof (sn3_20) / sizeof (u64), sn3_20 },
      { sizeof (sn3_21) / sizeof (u64), sn3_21 },
   };

   // b = 4
   //
   // Optimality prooven up to  m = 14

   cu64 sn4_1 [] = {1};
   cu64 sn4_2 [] = {1,      4};
   cu64 sn4_3 [] = {1,     20,        4};
   cu64 sn4_4 [] = {1,    148,       36,       4};
   cu64 sn4_5 [] = {1,    420,      116,      36,     4};
   cu64 sn4_6 [] = {1,    356,      148,      36,     4};
   cu64 sn4_7 [] = {1,   1380,     4532,     148,    36,     4};
   cu64 sn4_8 [] = {1,  10644,    52048,    1604,    84,    52,   4};
   cu64 sn4_9 [] = {1,   5732,    17108,    1172,   100,    20,   4};
   cu64 sn4_10[] = {1,  23028,   997012,    4068,  4372,   116,  20,   4};
   cu64 sn4_11[] = {1,  22132,    77332,    2004,  4644,   100,  20,   4};
   cu64 sn4_12[] = {1, 289396,  1307728,    5972, 17636,   404, 100,  20,  4};
   cu64 sn4_13[] = {1, 284084,    93972, 1050228,  4756,  1300,  84,  36,  4};
   cu64 sn4_14[] = {1,1141140, 30078544,  278164, 19444,  5428,1108, 228, 20, 4};
   cu64 sn4_15[] = {1,1137572,  4591956,  298612, 43988,  6532, 340, 180, 20, 4};
   cu64 sn4_16[] = {1,4560612,119772276, 3019604, 58708, 72244,2916,4260,116,36,4};
   cu64 sn4_17[] = {1,4560612, 19068292,67139444,118100,263668,5204, 660,244,20,4};

   const SNData snets4 [] =
   {
      { 0, 0 },
      { sizeof (sn4_1)  / sizeof (u64), sn4_1  },
      { sizeof (sn4_2)  / sizeof (u64), sn4_2  },
      { sizeof (sn4_3)  / sizeof (u64), sn4_3  },
      { sizeof (sn4_4)  / sizeof (u64), sn4_4  },
      { sizeof (sn4_5)  / sizeof (u64), sn4_5  },
      { sizeof (sn4_6)  / sizeof (u64), sn4_6  },
      { sizeof (sn4_7)  / sizeof (u64), sn4_7  },
      { sizeof (sn4_8)  / sizeof (u64), sn4_8  },
      { sizeof (sn4_9)  / sizeof (u64), sn4_9  },
      { sizeof (sn4_10) / sizeof (u64), sn4_10 },
      { sizeof (sn4_11) / sizeof (u64), sn4_11 },
      { sizeof (sn4_12) / sizeof (u64), sn4_12 },
      { sizeof (sn4_13) / sizeof (u64), sn4_13 },
      { sizeof (sn4_14) / sizeof (u64), sn4_14 },
      { sizeof (sn4_15) / sizeof (u64), sn4_15 },
      { sizeof (sn4_16) / sizeof (u64), sn4_16 },
      { sizeof (sn4_17) / sizeof (u64), sn4_17 },
   };

   // b = 5
   //
   // Optimality prooven up to  m = 13

   cu64 sn5_1 [] = {1};
   cu64 sn5_2 [] = {1,       5};
   cu64 sn5_3 [] = {1,      30,        5};
   cu64 sn5_4 [] = {1,     280,       55,       5};
   cu64 sn5_5 [] = {1,     930,      580,      55,      5};
   cu64 sn5_6 [] = {1,     805,      280,      55,      5};
   cu64 sn5_7 [] = {1,    3930,     1780,     355,     80,     5};
   cu64 sn5_8 [] = {1,   19730,   240830,    1530,    180,    30,    5};
   cu64 sn5_9 [] = {1,  149855,   802400,   41755,   1405,   330,   30,    5};
   cu64 sn5_10[] = {1, 2074705,   911180,  274525,  19280, 14930,  155,   55,  5};
   cu64 sn5_11[] = {1, 2920055, 30041005,  111605,   4030,  1705,  555,  105,  5};
   cu64 sn5_12[] = {1,  504805,  2026530,  132380,  12305,  2455,  230,   80,  5};
   cu64 sn5_13[] = {1, 2496455,186973250,10577780,  22730, 79730, 1080,  230, 30,  5};
   cu64 sn5_14[] = {1, 2458455, 39774780,48849830,  94905,  4555, 2155,  430, 80,  5};
   cu64 sn5_15[] = {1,51301705,502048480,14355055, 279780,411905, 4605,15805,280, 55,5};
   cu64 sn5_16[] = {1,51274805,245305530,10582980,2027855, 83530,17030,  980,205,105,5};

   const SNData snets5 [] =
   {
      { 0, 0 },
      { sizeof (sn5_1)  / sizeof (u64), sn5_1  },
      { sizeof (sn5_2)  / sizeof (u64), sn5_2  },
      { sizeof (sn5_3)  / sizeof (u64), sn5_3  },
      { sizeof (sn5_4)  / sizeof (u64), sn5_4  },
      { sizeof (sn5_5)  / sizeof (u64), sn5_5  },
      { sizeof (sn5_6)  / sizeof (u64), sn5_6  },
      { sizeof (sn5_7)  / sizeof (u64), sn5_7  },
      { sizeof (sn5_8)  / sizeof (u64), sn5_8  },
      { sizeof (sn5_9)  / sizeof (u64), sn5_9  },
      { sizeof (sn5_10) / sizeof (u64), sn5_10 },
      { sizeof (sn5_11) / sizeof (u64), sn5_11 },
      { sizeof (sn5_12) / sizeof (u64), sn5_12 },
      { sizeof (sn5_13) / sizeof (u64), sn5_13 },
      { sizeof (sn5_14) / sizeof (u64), sn5_14 },
      { sizeof (sn5_15) / sizeof (u64), sn5_15 },
      { sizeof (sn5_16) / sizeof (u64), sn5_16 },
   };

   // b = 7
   //
   // Optimality prooven for all m

   cu64 sn7_1 [] = {1};
   cu64 sn7_2 [] = {1,         7};
   cu64 sn7_3 [] = {1,        56,              7};
   cu64 sn7_4 [] = {1,       742,            105,        7};
   cu64 sn7_5 [] = {1,      3192,            497,       56,        7};
   cu64 sn7_6 [] = {1,     22547,          10052,      987,       56,       7};
   cu64 sn7_7 [] = {1,    763035,         108542,    14168,      987,     154,      7};
   cu64 sn7_8 [] = {1,    137795,          44303,   826294,      399,     105,      7};
   cu64 sn7_9 [] = {1,    964621,       17595221,    21567,     6034,     448,    203,     7};
   cu64 sn7_10[] = {1,  48858789,       34126302,   227563,  2491902,    6818,    840,   105,     7};
   cu64 sn7_11[] = {1,  61764262,      579706407, 31566297,  1188257,  229768,   8484,   399,   105,  7};
   cu64 sn7_12[] = {1,2487999850u,     957498325, 11239676,246221913, 1479905, 267498, 26565,   791,154, 7};
   cu64 sn7_13[] = {1,3307531955u,69918946524ull,775266443,159430082,31916108,3958612,235011,106379,448,56,7};

   const SNData snets7 [] =
   {
      { 0, 0 },
      { sizeof (sn7_1)  / sizeof (u64), sn7_1  },
      { sizeof (sn7_2)  / sizeof (u64), sn7_2  },
      { sizeof (sn7_3)  / sizeof (u64), sn7_3  },
      { sizeof (sn7_4)  / sizeof (u64), sn7_4  },
      { sizeof (sn7_5)  / sizeof (u64), sn7_5  },
      { sizeof (sn7_6)  / sizeof (u64), sn7_6  },
      { sizeof (sn7_7)  / sizeof (u64), sn7_7  },
      { sizeof (sn7_8)  / sizeof (u64), sn7_8  },
      { sizeof (sn7_9)  / sizeof (u64), sn7_9  },
      { sizeof (sn7_10) / sizeof (u64), sn7_10 },
      { sizeof (sn7_11) / sizeof (u64), sn7_11 },
      { sizeof (sn7_12) / sizeof (u64), sn7_12 },
      { sizeof (sn7_13) / sizeof (u64), sn7_13 },
   };

   // b = 8
   //
   // Optimality prooven up to m = 9

   cu64 sn8_1 [] = {1};
   cu64 sn8_2 [] = {1,       8};
   cu64 sn8_3 [] = {1,      72,       8};
   cu64 sn8_4 [] = {1,    1096,     136,        8};
   cu64 sn8_5 [] = {1,    5256,     712,       72,    8};
   cu64 sn8_6 [] = {1,   42440,   26696,     1736,   72,    8};
   cu64 sn8_7 [] = {1,  341640,  121928,     5000, 2312,   72,   8};
   cu64 sn8_8 [] = {1,  300296,   78664,     6600,  776,   72,   8};
   cu64 sn8_9 [] = {1, 2402184,17925896,    77640, 4808, 1160,  72,  8};
   cu64 sn8_10[] = {1, 2400968,  628808,    65160,11528,  648, 264,  8};
   cu64 sn8_11[] = {1,19208072, 5030664,134274760,75592,13256,1416,136,8};
   cu64 sn8_12[] = {1,19207368, 5033608,134255880,78728, 6600,1224,136,8};
   cu64 sn8_13[] = {1,19207368, 5033608,134255368,80200, 6088, 712,264,8};
   cu64 sn8_14[] = {1,19207368, 5033608,134255368,80200, 6088, 712,264,8};

   const SNData snets8 [] =
   {
      { 0, 0 },
      { sizeof (sn8_1)  / sizeof (u64), sn8_1  },
      { sizeof (sn8_2)  / sizeof (u64), sn8_2  },
      { sizeof (sn8_3)  / sizeof (u64), sn8_3  },
      { sizeof (sn8_4)  / sizeof (u64), sn8_4  },
      { sizeof (sn8_5)  / sizeof (u64), sn8_5  },
      { sizeof (sn8_6)  / sizeof (u64), sn8_6  },
      { sizeof (sn8_7)  / sizeof (u64), sn8_7  },
      { sizeof (sn8_8)  / sizeof (u64), sn8_8  },
      { sizeof (sn8_9)  / sizeof (u64), sn8_9  },
      { sizeof (sn8_10) / sizeof (u64), sn8_10 },
      { sizeof (sn8_11) / sizeof (u64), sn8_11 },
      { sizeof (sn8_12) / sizeof (u64), sn8_12 },
      { sizeof (sn8_13) / sizeof (u64), sn8_13 },
      { sizeof (sn8_14) / sizeof (u64), sn8_14 },
   };

   // b = 9
   //
   // Optimality prooven up to m = 9

   cu64 sn9_1 [] = {1};
   cu64 sn9_2 [] = {1,       9};
   cu64 sn9_3 [] = {1,      90,       9};
   cu64 sn9_4 [] = {1,    2277,     171,        9};
   cu64 sn9_5 [] = {1,    8352,     819,      171,     9};
   cu64 sn9_6 [] = {1,   81576,   43101,     1872,    90,    9};
   cu64 sn9_7 [] = {1,  607023,  151722,    15885,  4140,   90,   9};
   cu64 sn9_8 [] = {1,20173383, 2293605,   369936, 25929, 3006, 252,  9};
   cu64 sn9_9 [] = {1, 5389263, 1412001,   106524, 17424, 1791,  90,  9};
   cu64 sn9_10[] = {1, 5387562, 1262961,   311211, 25767, 1305, 414,  9};
   cu64 sn9_11[] = {1,48488310,11448711,387495504,139329,10539, 900,252,9};
   cu64 sn9_12[] = {1,48488310,11423520,  4032594,290394,36621,1953,252,9};

   const SNData snets9 [] =
   {
      { 0, 0 },
      { sizeof (sn9_1)  / sizeof (u64), sn9_1  },
      { sizeof (sn9_2)  / sizeof (u64), sn9_2  },
      { sizeof (sn9_3)  / sizeof (u64), sn9_3  },
      { sizeof (sn9_4)  / sizeof (u64), sn9_4  },
      { sizeof (sn9_5)  / sizeof (u64), sn9_5  },
      { sizeof (sn9_6)  / sizeof (u64), sn9_6  },
      { sizeof (sn9_7)  / sizeof (u64), sn9_7  },
      { sizeof (sn9_8)  / sizeof (u64), sn9_8  },
      { sizeof (sn9_9)  / sizeof (u64), sn9_9  },
      { sizeof (sn9_10) / sizeof (u64), sn9_10 },
      { sizeof (sn9_11) / sizeof (u64), sn9_11 },
      { sizeof (sn9_12) / sizeof (u64), sn9_12 },
   };

   // The "list of lists"
   
   struct SNLists
   {
      unsigned num;
      const SNData* lists;
   };

   const SNLists snets [] =
   {
      { 0, 0 },
      { 0, 0 },
      { sizeof (snets2) / sizeof (SNData), snets2 },
      { sizeof (snets3) / sizeof (SNData), snets3 },
      { sizeof (snets4) / sizeof (SNData), snets4 },
      { sizeof (snets5) / sizeof (SNData), snets5 },
      { 0, 0 },
      { sizeof (snets7) / sizeof (SNData), snets7 },
      { sizeof (snets8) / sizeof (SNData), snets8 },
      { sizeof (snets9) / sizeof (SNData), snets9 }
   };
}


/**
 *  max Base For Shift Net ()
 */

unsigned L::maxBaseForShiftNet ()
{
   return sizeof (snets) / sizeof (SNLists) - 1;
}


/**
 *  max M For Shift Net ()
 */

unsigned L::maxMForShiftNet (unsigned base)
{
   if (base > maxBaseForShiftNet())  throw InternalError (__FILE__, __LINE__);

   return snets[base].num - 1;
}


/**
 *  optimal Shift Net T ()
 */

unsigned L::optimalShiftNetT (unsigned base, unsigned m)
{
   if (m > maxMForShiftNet (base))  throw InternalError (__FILE__, __LINE__);

   return m - snets[base].lists[m].length;
}


/**
 *  init Shift Net ()
 */

void L::Private::initShiftNetFirstDim (GeneratorMatrix &gm)
{
   unsigned b = gm.getBase();
   unsigned m = gm.getM();

   if (gm.getDimension() > m)  throw InternalError (__FILE__, __LINE__);

   if (m > maxMForShiftNet (b))  throw InternalError (__FILE__, __LINE__);

   const SNData& r = snets[b].lists[m];
   const u64* data = r.data;

   for (unsigned i = 0; i < r.length; ++i)
   {
      gm.vSetPackedRowVector (0, i, data[i]);
   }
}


