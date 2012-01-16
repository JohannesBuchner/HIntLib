/*
 *  HIntLib  -  Library for High-dimensional Numerical Integration
 *
 *  Copyright (C) 2002,03,04,05,06  Rudolf Schuerer <rudolf.schuerer@sbg.ac.at>
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

#include <HIntLib/shiftnet.h>

#ifdef HINTLIB_USE_INTERFACE_IMPLEMENTATION
#pragma implementation
#endif

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
    * The orgiginal table, which is still used for  m <= 21, was calculated by
    * Wolfgang Schmid and is published in
    *
    * [1] Wolfgnag Ch. Schmid. Shift-Nets: A New Class of Binary Digital
    *     (t,m,s)-Nets. In Harald Niederreiter et al, editors, Monte Carlo and
    *     Quasi Monte Carlo Methods 1996, volume 127 of Lecture Notes in
    *     Statistics. Springer, 1997.
    *
    * Additional results have been calcuated and are published in
    *
    * [2] Wolfgang Schmid and Rudolf Schuerer. Shift-Nets and Salzburg Tables:
    *     Power Computing in Number-Theoretical Numerics.  In H. Efinger and
    *     A. Uhl, editors, Scientific Computing in Salzburg--Festschrift on
    *     the Occasion of Peter Zinterhofs 60th Birthday,
    *     volume 189 of books@ocg.at, pages 175-184.
    *     Oesterreichische Computer Gesellschaft, 2005.
    *
    * - For m=21, the second-to-last number is 6 instead of 62. This results in
    *   the same t value (t=11), but is the lexicographically first solution
    *   with this t value (as are all other entries).
    *
    * - For m >= 22, superior results were calculated by Rudolf Schuerer.
    *
    * - The results for m >= 32 were calculated by a randomized algorithm.
    *   Therefore they are not the lexicographically first solution.
    *
    * Optimality is proven for  m <= 27,  m = 29  and  m = 32.
    *
    * The values have been confirmed by test_shiftnet up to m <= 35.
    */

   //     b m          k  t
   cu64 sn2_1 [] = /*  1  0 */ {1};
   cu64 sn2_2 [] = /*  2  0 */ {1,     2};
   cu64 sn2_3 [] = /*  3  0 */ {1,     6,       2};
   cu64 sn2_4 [] = /*  3  1 */ {1,     6,       2};
   cu64 sn2_5 [] = /*  4  1 */ {1,    14,      18,      2};
   cu64 sn2_6 [] = /*  4  2 */ {1,    14,      18,      2};
   cu64 sn2_7 [] = /*  5  2 */ {1,    46,      26,      6,     2};
   cu64 sn2_8 [] = /*  5  3 */ {1,    46,      26,      6,     2};
   cu64 sn2_9 [] = /*  6  3 */ {1,    94,     178,     38,    10,    2};
   cu64 sn2_10[] = /*  6  4 */ {1,    94,     170,     38,    10,    2};
   cu64 sn2_11[] = /*  7  4 */ {1,   492,    1638,     86,    14,   18,    2};
   cu64 sn2_12[] = /*  7  5 */ {1,   378,     590,    158,    42,    6,    2};
   cu64 sn2_13[] = /*  7  6 */ {1,   366,     572,    154,    14,   18,    2};
   cu64 sn2_14[] = /*  8  6 */ {1,   734,   10652,    372,    46,   22,   10,   2};
   cu64 sn2_15[] = /*  8  7 */ {1,   734,    3308,    372,    46,   22,   10,   2};
   cu64 sn2_16[] = /*  8  8 */ {1,   758,    3132,    414,    46,   22,   10,   2};
   cu64 sn2_17[] = /*  9  8 */ {1,  3562,   29276,    870,    62,  198,   74,   6,   2};
   cu64 sn2_18[] = /*  9  9 */ {1,  3450,    9814,   2020,   458,   54,   26,   6,   2};
   cu64 sn2_19[] = /*  9 10 */ {1,  3562,    8566,   1884,   236,   30,   38,  10,   2};
   cu64 sn2_20[] = /* 10 10 */ {1, 10094,  301178,  17650,  2300,  314,   78,  22,  10,  2};
   cu64 sn2_21[] = /* 10 11 */ {1, 10106,   54364,   3974,  4214,  158,  294,  42,   6,  2};
   cu64 sn2_22[] = /* 11 11 */ {1,248292, 2954238,  30118,  1598,  462, 2134, 106,  14, 18, 2};
   cu64 sn2_23[] = /* 11 12 */ {1, 40302,  468876,  69080,  1852,  222, 2378,  46,  22, 10, 2};
   cu64 sn2_24[] = /* 11 13 */ {1, 40302,  213628,  87254,  4978,  238, 1640,  30,  38, 10, 2};
   cu64 sn2_25[] = /* 12 13 */ {1,953816,19043324,1386616,  3454, 4766, 8422, 426,  54, 14,18, 2};
   cu64 sn2_26[] = /* 12 14 */ {1, 85470, 7399828, 177330, 16084,  702, 1388,2126, 114, 22,10, 2};
   cu64 sn2_27[] = /* 12 15 */ {1, 81366, 1969322, 164826,  5950,16630, 2382, 604,  46, 22,10, 2};
   cu64 sn2_28[] = /* 12 16 */ {1, 81366, 1350324, 548382, 17334, 1756, 2222, 286, 102, 42, 6, 2};
   cu64 sn2_29[] = /* 13 16 */ {1,441530,80331402,3209780, 16210, 3502,16734,1590, 186,102,14,18,2};
   cu64 sn2_30[] = /* 13 17 */ {1,380762,12737598,1623494, 27322, 4958, 1462,2158, 572,142,22,10,2};
   cu64 sn2_31[] = /* 13 18 */ {1,380762, 5422020,2209652,195832, 5606,  926,2166,1196, 58,14,18,2};
   // From here on the results are from a random search
   cu64 sn2_32[] = /* 14 18 */ {1ull,3400052811ull,3597241718ull,1836987300ull,1594329648ull,2363152018ull,3663481202ull,254ull,2944933244ull,1288542503ull,85100844ull,14ull,2241408841ull,2};
   cu64 sn2_33[] = /* 14 19 */ {1ull,2759196615ull,3834602671ull,7921089133ull,1280288929ull,3625398249ull,510ull,7374213406ull,2196151764ull,5437991504ull,30ull,2962542481ull,6ull,2};
   cu64 sn2_34[] = /* 14 20 */ {1ull,12826654228ull,2707469858ull,16374587625ull,8194448799ull,2572780966ull,510ull,4943951277ull,16040614009ull,15695348889ull,30ull,3238678136ull,6ull,2};
   cu64 sn2_35[] = /* 14 21 */ {1ull,24577979269ull,21243745729ull,12004716326ull,3409890995ull,33037036302ull,510ull,7615552080ull,11617703845ull,28469892938ull,30ull,10842457427ull,6ull,2};
   cu64 sn2_36[] = /* 15 21 */ {1ull,53537520163ull,13671743026ull,6702687059ull,66747765040ull,38244040695ull,45645222695ull,510ull,2239436683ull,3771625885ull,26405587062ull,30ull,53777595871ull,6ull,2};
   cu64 sn2_37[] = /* 15 22 */ {1ull,30011467176ull,113605748158ull,36293300906ull,12640341894ull,7563717446ull,2163147254ull,510ull,78309912ull,554048202ull,455581864ull,30ull,17842890ull,6ull,2};
   cu64 sn2_38[] = /* 15 23 */ {1ull,148492682932ull,11518364980ull,124590427048ull,22615365708ull,7189328024ull,2581307720ull,510ull,35219344476ull,1875510500ull,392236660ull,30ull,1052906ull,6ull,2};
   cu64 sn2_39[] = /* 15 24 */ {1ull,82570019104ull,434076037826ull,174728916490ull,33771737014ull,38815839758ull,5875249768ull,510ull,9883641034ull,157263040ull,2006917276ull,30ull,52431948ull,6ull,2};
   cu64 sn2_40[] = /* 16 24 */ {1ull,1035402270526ull,390985785214ull,99558255746ull,158217747614ull,14599991424ull,24851978862ull,1022ull,37741046940ull,2082600178ull,6808527264ull,3210584282ull,30ull,335298856ull,6ull,2};
   cu64 sn2_41[] = /* 16 25 */ {1ull,1277142714912ull,466582673580ull,113476261492ull,709435844222ull,23901600938ull,172012960994ull,1022ull,47229464928ull,5615167710ull,676679792ull,3222222122ull,30ull,319398920ull,6ull,2};
   cu64 sn2_42[] = /* 16 26 */ {1ull,1433371124986ull,737221789938ull,421658708568ull,2234522916674ull,50693493886ull,213240453194ull,1022ull,84672829724ull,20566587638ull,1122472306ull,6915589358ull,30ull,842651072ull,6ull,2};
   cu64 sn2_43[] = /* 16 27 */ {1ull,293209544920ull,8401701700494ull,2795908687792ull,212968921370ull,1697992550714ull,7235861092ull,1022ull,72516661730ull,17657192562ull,560680050ull,595333875798ull,30ull,8962439406ull,6ull,2};
   cu64 sn2_44[] = /* 17 27 */ {1ull,298315930204ull,8311495190206ull,2420256721058ull,1258653779276ull,185153503386ull,45416911232ull,8893150996346ull,1022ull,9108962308ull,4728650218ull,551289606540ull,89155065228ull,30ull,126506020ull,6ull,2};
   cu64 sn2_45[] = /* 17 28 */ {1ull,26014575068178ull,10708696090380ull,3578856221492ull,718133188606ull,4448982373512ull,1532191869804ull,28044619920ull,1022ull,345180565764ull,12246413340ull,176967824756ull,76968517744ull,30ull,5290591396ull,6ull,2};
   cu64 sn2_46[] = /* 17 29 */ {1ull,46533585702876ull,12025878902106ull,23100588854114ull,7350078185148ull,4220226572076ull,216366703432ull,136567547600ull,1022ull,1668962684324ull,617321042286ull,320902339706ull,55428896000ull,30ull,4793209192ull,6ull,2};
   cu64 sn2_47[] = /* 17 30 */ {1ull,133673232793888ull,64819726304770ull,12472881511556ull,2957890168044ull,6303984716318ull,434759128228ull,217035929718ull,1022ull,82789985520ull,19281268114922ull,58444043330ull,1106694571308ull,30ull,564457325806ull,6ull,2};
   cu64 sn2_48[] = /* 18 30 */ {1ull,184144246828030ull,102437169982012ull,41586268975600ull,21931174314804ull,17217612111798ull,3803617259032ull,5350474515212ull,2046ull,1700606148868ull,1077250535570ull,6380229310ull,173947840546ull,62ull,71127994630ull,26986195210ull,6ull,2};
   cu64 sn2_49[] = /* 18 31 */ {1ull,507681286982174ull,9990737124620ull,21292769649496ull,41537653652992ull,5379760118938ull,141727606880992ull,2860800581240ull,2046ull,71097656099700ull,197257452014ull,1760205335410ull,326904410570ull,62ull,626417222164ull,35385409674ull,6ull,2};
   cu64 sn2_50[] = /* 18 32 */ {1ull,1070631240403486ull,211253128345358ull,337998534917822ull,77868296958120ull,9584267301886ull,7962547328126ull,37935937179502ull,2046ull,776575435034ull,18099853053976ull,3734654243798ull,12502182878ull,62ull,1122590091358ull,140433269702ull,6ull,2};
   cu64 sn2_51[] = /* 18 33 */ {1ull,2133475928475316ull,1056227505680260ull,500692499039652ull,29161687783876ull,43717129559390ull,73384511578028ull,157640003446884ull,9157819674750ull,1022ull,2195859117384ull,2471851981150ull,5186684162452ull,609197967832ull,30ull,51949723106ull,6ull,2};

#define HINTLIB_SNBM(b,m) { sizeof(sn##b##_##m)  / sizeof(u64), sn##b##_##m },
#define HINTLIB_SNM(m) HINTLIB_SNBM(2,m)

   const SNData snets2 [] =
   {
      { 0, 0 },
      HINTLIB_SNM(1)
      HINTLIB_SNM(2)
      HINTLIB_SNM(3)
      HINTLIB_SNM(4)
      HINTLIB_SNM(5)
      HINTLIB_SNM(6)
      HINTLIB_SNM(7)
      HINTLIB_SNM(8)
      HINTLIB_SNM(9)
      HINTLIB_SNM(10)
      HINTLIB_SNM(11)
      HINTLIB_SNM(12)
      HINTLIB_SNM(13)
      HINTLIB_SNM(14)
      HINTLIB_SNM(15)
      HINTLIB_SNM(16)
      HINTLIB_SNM(17)
      HINTLIB_SNM(18)
      HINTLIB_SNM(19)
      HINTLIB_SNM(20)
      HINTLIB_SNM(21)
      HINTLIB_SNM(22)
      HINTLIB_SNM(23)
      HINTLIB_SNM(24)
      HINTLIB_SNM(25)
      HINTLIB_SNM(26)
      HINTLIB_SNM(27)
      HINTLIB_SNM(28)
      HINTLIB_SNM(29)
      HINTLIB_SNM(30)
      HINTLIB_SNM(31)
      HINTLIB_SNM(32)
      HINTLIB_SNM(33)
      HINTLIB_SNM(34)
      HINTLIB_SNM(35)
      HINTLIB_SNM(36)
      HINTLIB_SNM(37)
      HINTLIB_SNM(38)
      HINTLIB_SNM(39)
      HINTLIB_SNM(40)
      HINTLIB_SNM(41)
      HINTLIB_SNM(42)
      HINTLIB_SNM(43)
      HINTLIB_SNM(44)
      HINTLIB_SNM(45)
      HINTLIB_SNM(46)
      HINTLIB_SNM(47)
      HINTLIB_SNM(48)
      HINTLIB_SNM(49)
      HINTLIB_SNM(50)
      HINTLIB_SNM(51)
   };


   /*
    *  bases   b > 2
    *
    *  For all bases b > 2, the parameters for optimal shift nets were
    *  calculated by Rudolf Schuerer based on work by Wolfgang Schmid.
    *
    *  See [2]
    */

   /*  b = 3
    *
    *  Optimality proven up to  m = 18
    *
    *  The values are confirmed by test_shiftnet up to m <= 28.
    */

   //     b m          t
   cu64 sn3_1 [] = /*  0 */ {1};
   cu64 sn3_2 [] = /*  0 */ {1,       3};
   cu64 sn3_3 [] = /*  0 */ {1,      12,         3};
   cu64 sn3_4 [] = /*  1 */ {1,      12,         3};
   cu64 sn3_5 [] = /*  1 */ {1,      39,        21,           3};
   cu64 sn3_6 [] = /*  1 */ {1,     129,       264,          12,      3};
   cu64 sn3_7 [] = /*  1 */ {1,     480,      1602,         111,     21,     3};
   cu64 sn3_8 [] = /*  2 */ {1,     372,       822,          39,     21,     3};
   cu64 sn3_9 [] = /*  2 */ {1,    3357,      7824,         120,    255,    21,     3};
   cu64 sn3_10[] = /*  3 */ {1,    2577,      2061,         201,    273,    21,     3};
   cu64 sn3_11[] = /*  3 */ {1,    8391,     23034,        2856,    120,   255,    21,     3};
   cu64 sn3_12[] = /*  4 */ {1,    7680,     24660,        2388,    426,    66,    21,     3};
   cu64 sn3_13[] = /*  4 */ {1,   32205,    657297,       62004,    885,   282,    93,    21,     3};
   cu64 sn3_14[] = /*  5 */ {1,   25509,    134274,      177600,   1197,   120,   309,    12,     3};
   cu64 sn3_15[] = /*  5 */ {1,  280578,   5344227,      557100,   1452,  2388,   309,   102,    12,  3};
   cu64 sn3_16[] = /*  6 */ {1,  204762,    659883,       37101,   6222,  7509,   468,    48,    12,  3};
   cu64 sn3_17[] = /*  7 */ {1,  204591,    562188,       19236,  61677,   372,   768,    93,    21,  3};
   cu64 sn3_18[] = /*  7 */ {1,  615666,   6805764,     1651467,  15222, 24186,  6300,   282,   102, 12,  3};
   cu64 sn3_19[] = /*  8 */ {1,  614289,   5110500,     1668576,  37533,  3036,   399,  6600,    93, 21,  3};
   cu64 sn3_20[] = /*  8 */ {1, 1886214, 546400605,    10398675,  36912, 62517,  7743,  3249,   282, 93, 21, 3};
   cu64 sn3_21[] = /*  9 */ {1, 1842582,  24946203,    43316163, 114276, 27444, 10074,   696,   831, 39, 21, 3};
   cu64 sn3_22[] = /* 10 */ {1, 1842582,  19423098,    43328289, 321699,  6528, 27489,  6708,   282,102, 12, 3};
   cu64 sn3_23[] = /* 10 */ {1,16562037, 224237595,    93675468, 345963, 69636, 23313,532137,  1083,156,282,12, 3};
   cu64 sn3_24[] = /* 11 */ {1,16562037, 173212947,    48652140, 194637, 97992, 46326,532371,  1317,237,279,21, 3};
   cu64 sn3_25[] = /* 11 */ {1,49871991,5575327680ull,658310898,3000666,194304,102270, 29439,532848,372,867,66,21,3};
   // From here on the results are from a random search
   cu64 sn3_26[] = /* 12 */ {1ull,831545819163ull,956521391988ull,260688658275ull,48429921477ull,9960011049ull,9840ull,538743594ull,11958528354ull,2445441666ull,120ull,237579225ull,12ull,3ull};
   cu64 sn3_27[] = /* 13 */ {1ull,4349953663557ull,748079138850ull,989655954963ull,121657319115ull,4769969253ull,9840ull,55260402624ull,3369516996ull,11100742107ull,120ull,11201406ull,12ull,3ull};
   cu64 sn3_28[] = /* 13 */ {1ull,15388703864799ull,6893307523626ull,2537674332930ull,335970385143ull,191482050813ull,86038792785ull,9840ull,3669504225ull,37261977ull,10812030951ull,120ull,1085657724ull,12ull,3ull};
   cu64 sn3_29[] = /* 14 */ {1ull,56060805292983ull,12806710597830ull,376327175565ull,4460437265097ull,4491938274ull,1122732918360ull,9840ull,734257374ull,211629491199ull,1183013814ull,120ull,21260239686ull,12ull,3ull};
   cu64 sn3_30[] = /* 15 */ {1ull,187364293303245ull,55127765350518ull,14263812008178ull,6978000419871ull,1039256092212ull,153472697580ull,655379448273ull,92368032138ull,26765484834ull,3062203995ull,7275052251ull,150293694ull,439979991ull,3754257ull};
   cu64 sn3_31[] = /* 15 */ {1ull,345273610454289ull,66611023781451ull,76793907820770ull,10560503714442ull,7493607109287ull,2193213354558ull,241407949449ull,594358940127ull,39184112973ull,31227510828ull,4254110004ull,128311578ull,2601080433ull,1040573004ull,168147957ull};
   cu64 sn3_32[] = /* 16 */ {1ull,1063338288959436ull,484282411052565ull,33931505674572ull,139967914791717ull,14664454659066ull,6061894363398ull,1222059368691ull,31471203669ull,774934985424ull,119852701026ull,1258882191ull,11476832352ull,7837272021ull,954380469ull,382935660ull};
   cu64 sn3_33[] = /* 17 */ {1ull,4922860331480172ull,167374415042151ull,1064833886562846ull,264485254558635ull,68491359297993ull,20682536574753ull,1200620318337ull,5672741795577ull,221685130824ull,648838299303ull,2856365043ull,49412434383ull,650979576ull,28232975934ull,175406580ull};
   cu64 sn3_34[] = /* 17 */ {1ull,15719548369416267ull,1817150935833450ull,2275042912304676ull,33717623166240ull,285529903397313ull,78345305936289ull,18513079974066ull,3134510931600ull,1293456551319ull,773124176034ull,60838617066ull,94361698467ull,24348670053ull,10271091828ull,1238366325ull,968429841ull};
   cu64 sn3_35[] = /* 18 */ {1ull,27159102597702264ull,8468303868011523ull,4704232689980814ull,1650033922554870ull,179060628675741ull,239869327759818ull,57867230294124ull,11000451936555ull,7198961575410ull,1931228761743ull,512519296860ull,257836881960ull,19783770516ull,69986925537ull,5291662836ull,2299590537ull};
   cu64 sn3_36[] = /* 19 */ {1ull,77190647696701971ull,185050975369170ull,43074129199209300ull,11364751046404659ull,4349529151503501ull,232911687384336ull,1249410141409935ull,53363735943030ull,20091120417720ull,4369845304806ull,253731928644ull,1132454999388ull,648838299303ull,2856365043ull,49412434383ull,650979576ull};
   cu64 sn3_37[] = /* 19 */ {1ull,77190647696701971ull,300470932508463474ull,10428196941054006ull,21476162021395626ull,4866190026058824ull,507007675151793ull,680673399647676ull,37017416643204ull,1825883309205ull,145705350912588ull,18323395842990ull,3140159970765ull,678969036978ull,62887105554ull,197633321403ull,6136823187ull,13854941820ull};
   cu64 sn3_38[] = /* 20 */ {1ull,527474553587699334ull,28700862730990161ull,254383508257041063ull,114075225489970767ull,2836290903770058ull,428393248093257ull,6228421646530467ull,1379458309371078ull,27492545070948ull,71732411393340ull,4867944869949ull,16258589757183ull,2206878490593ull,673210185009ull,125271702402ull,21299805339ull,5737812291ull};

#undef  HINTLIB_SNM
#define HINTLIB_SNM(m) HINTLIB_SNBM(3,m)
   const SNData snets3 [] =
   {
      { 0, 0 },
      HINTLIB_SNM(1)
      HINTLIB_SNM(2)
      HINTLIB_SNM(3)
      HINTLIB_SNM(4)
      HINTLIB_SNM(5)
      HINTLIB_SNM(6)
      HINTLIB_SNM(7)
      HINTLIB_SNM(8)
      HINTLIB_SNM(9)
      HINTLIB_SNM(10)
      HINTLIB_SNM(11)
      HINTLIB_SNM(12)
      HINTLIB_SNM(13)
      HINTLIB_SNM(14)
      HINTLIB_SNM(15)
      HINTLIB_SNM(16)
      HINTLIB_SNM(17)
      HINTLIB_SNM(18)
      HINTLIB_SNM(19)
      HINTLIB_SNM(20)
      HINTLIB_SNM(21)
      HINTLIB_SNM(22)
      HINTLIB_SNM(23)
      HINTLIB_SNM(24)
      HINTLIB_SNM(25)
      HINTLIB_SNM(26)
      HINTLIB_SNM(27)
      HINTLIB_SNM(28)
      HINTLIB_SNM(29)
      HINTLIB_SNM(30)
      HINTLIB_SNM(31)
      HINTLIB_SNM(32)
      HINTLIB_SNM(33)
      HINTLIB_SNM(34)
      HINTLIB_SNM(35)
      HINTLIB_SNM(36)
      HINTLIB_SNM(37)
      HINTLIB_SNM(38)
   };

   /*
    *  b = 4
    *
    *  Optimality proven up to  m = 16
    *
    *  The values have been confirmed by test_shiftnet up to m <= 27.
    */

   //     b m          t
   cu64 sn4_1 [] = /*  0 */ {1};
   cu64 sn4_2 [] = /*  0 */ {1,       4};
   cu64 sn4_3 [] = /*  0 */ {1,      20,            4};
   cu64 sn4_4 [] = /*  0 */ {1,     148,           36,       4};
   cu64 sn4_5 [] = /*  0 */ {1,     420,          116,      36,     4};
   cu64 sn4_6 [] = /*  1 */ {1,     356,          148,      36,     4};
   cu64 sn4_7 [] = /*  1 */ {1,    1380,         4532,     148,    36,      4};
   cu64 sn4_8 [] = /*  1 */ {1,   10644,        52048,    1604,    84,     52,      4};
   cu64 sn4_9 [] = /*  2 */ {1,    5732,        17108,    1172,   100,     20,      4};
   cu64 sn4_10[] = /*  2 */ {1,   23028,       997012,    4068,  4372,    116,     20,    4};
   cu64 sn4_11[] = /*  3 */ {1,   22132,        77332,    2004,  4644,    100,     20,    4};
   cu64 sn4_12[] = /*  3 */ {1,  289396,      1307728,    5972, 17636,    404,    100,   20,   4};
   cu64 sn4_13[] = /*  4 */ {1,  284084,        93972, 1050228,  4756,   1300,     84,   36,   4};
   cu64 sn4_14[] = /*  4 */ {1, 1141140,     30078544,  278164, 19444,   5428,   1108,  228,  20,  4};
   cu64 sn4_15[] = /*  5 */ {1, 1137572,      4591956,  298612, 43988,   6532,    340,  180,  20,  4};
   cu64 sn4_16[] = /*  5 */ {1, 4560612,    119772276, 3019604, 58708,  72244,   2916, 4260, 116, 36, 4};
   cu64 sn4_17[] = /*  6 */ {1, 4560612,     19068292,67139444,118100, 263668,   5204,  660, 244, 20, 4};
   cu64 sn4_18[] = /*  6 */ {1,18503396,2487730820ull,75383828,701332,  83412,1051508, 6096, 340,180,20,4};
   cu64 sn4_19[] = /*  7 */ {1,18242452,    337414756,76492596,775380,1144548,   5748,16788,1236,164,20,4};
   // From here on the results are from a random search
   cu64 sn4_20[] = /*  7 */ {1ull,937784895952ull,249417171952ull,56210677404ull,8272765088ull,2435667156ull,132622044ull,591894324ull,57573840ull,3595752ull,493624ull,210912ull,4248364ull};
   cu64 sn4_21[] = /*  8 */ {1ull,1772846415736ull,440855311344ull,152052235328ull,51366103596ull,5122122452ull,4129501316ull,444074212ull,60811204ull,69699204ull,11047956ull,3468968ull,482832ull};
   cu64 sn4_22[] = /*  9 */ {1ull,3623992660952ull,865508287776ull,9007680111560ull,248796364948ull,17723518476ull,7842066276ull,3005546384ull,28125192ull,683956688ull,74385128ull,215964ull,15759096ull};
   cu64 sn4_23[] = /*  9 */ {1ull,70266251015452ull,8935502466700ull,3261158107940ull,989694124476ull,79510624512ull,45658599452ull,10220742120ull,970824916ull,1250289064ull,23006492ull,4384612ull,204887164ull,388916ull};
   cu64 sn4_24[] = /* 10 */ {1ull,218711403506560ull,29504483308196ull,4168179186844ull,4732450184080ull,1001563897780ull,31859542148ull,214869027772ull,1320736844ull,5079493156ull,650474120ull,26223236ull,144858760ull,16677708ull};
   cu64 sn4_25[] = /* 10 */ {1ull,537461395837324ull,252680203920264ull,40566860476324ull,4652501602108ull,4066766875196ull,656175525612ull,128208610056ull,66532140064ull,13244247392ull,223398544ull,574207732ull,58909132ull,13085464ull,1077144084ull};
   cu64 sn4_26[] = /* 11 */ {1ull,3742325116467864ull,766707040973232ull,22098817439256ull,2996474399896ull,70559079903472ull,13943533391908ull,435142368696ull,267436481304ull,56506379464ull,4671730180ull,1458597512ull,977326820ull,43138164ull,1183848ull};
   cu64 sn4_27[] = /* 12 */ {1ull,12749524371208856ull,2733838284066092ull,178998454738648ull,6331963249076ull,38156812383560ull,847426424823488ull,4010925652704ull,828122160260ull,150461396412ull,1742054028ull,35188679480ull,216269972ull,5142566648ull,557900976ull};
   cu64 sn4_28[] = /* 12 */ {1ull,35783389814156876ull,12772709039164100ull,2898193552600520ull,762581604812552ull,131499109788860ull,8826849762584ull,21558167562352ull,2557206086728ull,904137223032ull,259185783852ull,19013578696ull,11990942168ull,1821326444ull,38813840ull,820262168ull};
   cu64 sn4_29[] = /* 13 */ {1ull,30602972298046328ull,78402032901426348ull,6878131701808304ull,1066578001440552ull,2449752510451940ull,35628600179180ull,71593728865876ull,5085974439164ull,3042343897112ull,994649578400ull,64923537460ull,4849046692ull,141243824632ull,3726286364ull,669795700ull};
   cu64 sn4_30[] = /* 13 */ {1ull,1017647077718317700ull,103217403691044304ull,1394208590519988ull,27941234897645624ull,4780250131328924ull,45838144836496ull,1069862669066844ull,86023953463120ull,6610353398308ull,1869627471588ull,150152766716ull,829350111428ull,3968765668ull,21765128780ull,13665334280ull,265808296ull};
   cu64 sn4_31[] = /* 14 */ {1ull,4264147572608957080ull,746995863390780136ull,48635266115903652ull,226500436901935528ull,13083478898102448ull,1959694928829312ull,823299057928012ull,98766130172544ull,33358749797040ull,9834869621576ull,413529122088ull,126402468764ull,23453992612ull,2215456548844ull,5988248684ull,1804951684ull};
   cu64 sn4_32[] = /* 14 */ {1ull,18308419720600247156ull,1890250419963433164ull,328597764567677192ull,169519399395178992ull,58965199545906304ull,2253011623428256ull,907302055028648ull,9205707942388344ull,259630366807196ull,2608829125220ull,62586635076100ull,13790225728004ull,460151186612ull,252055709852ull,25053864028ull,1766846672ull,13161730292ull};
   // Results for larger m are available, but cannot be stored in this format

#undef  HINTLIB_SNM
#define HINTLIB_SNM(m) HINTLIB_SNBM(4,m)
   const SNData snets4 [] =
   {
      { 0, 0 },
      HINTLIB_SNM(1)
      HINTLIB_SNM(2)
      HINTLIB_SNM(3)
      HINTLIB_SNM(4)
      HINTLIB_SNM(5)
      HINTLIB_SNM(6)
      HINTLIB_SNM(7)
      HINTLIB_SNM(8)
      HINTLIB_SNM(9)
      HINTLIB_SNM(10)
      HINTLIB_SNM(11)
      HINTLIB_SNM(12)
      HINTLIB_SNM(13)
      HINTLIB_SNM(14)
      HINTLIB_SNM(15)
      HINTLIB_SNM(16)
      HINTLIB_SNM(17)
      HINTLIB_SNM(18)
      HINTLIB_SNM(19)
      HINTLIB_SNM(20)
      HINTLIB_SNM(21)
      HINTLIB_SNM(22)
      HINTLIB_SNM(23)
      HINTLIB_SNM(24)
      HINTLIB_SNM(25)
      HINTLIB_SNM(26)
      HINTLIB_SNM(27)
      HINTLIB_SNM(28)
      HINTLIB_SNM(29)
      HINTLIB_SNM(30)
      HINTLIB_SNM(31)
      HINTLIB_SNM(32)
   };

   /*
    *  b = 5
    *
    *  Optimality proven for  m <= 15  and for  m = 17
    *
    *  The values have been confirmed by test_shiftnet up to m <= 25,
    *  i.e., for all m.
    */

   //     b m          t
   cu64 sn5_1 [] = /*  0 */ {1};
   cu64 sn5_2 [] = /*  0 */ {1,        5};
   cu64 sn5_3 [] = /*  0 */ {1,       30,         5};
   cu64 sn5_4 [] = /*  0 */ {1,      280,        55,          5};
   cu64 sn5_5 [] = /*  0 */ {1,      930,       580,         55,      5};
   cu64 sn5_6 [] = /*  1 */ {1,      805,       280,         55,      5};
   cu64 sn5_7 [] = /*  1 */ {1,     3930,      1780,        355,     80,      5};
   cu64 sn5_8 [] = /*  1 */ {1,    19730,    240830,       1530,    180,     30,    5};
   cu64 sn5_9 [] = /*  1 */ {1,   149855,    802400,      41755,   1405,    330,   30,    5};
   cu64 sn5_10[] = /*  1 */ {1,  2074705,    911180,     274525,  19280,  14930,  155,   55,   5};
   cu64 sn5_11[] = /*  2 */ {1,  2920055,  30041005,     111605,   4030,   1705,  555,  105,   5};
   cu64 sn5_12[] = /*  3 */ {1,   504805,   2026530,     132380,  12305,   2455,  230,   80,   5};
   cu64 sn5_13[] = /*  3 */ {1,  2496455, 186973250,   10577780,  22730,  79730, 1080,  230,  30,  5};
   cu64 sn5_14[] = /*  4 */ {1,  2458455,  39774780,   48849830,  94905,   4555, 2155,  430,  80,  5};
   cu64 sn5_15[] = /*  4 */ {1, 51301705, 502048480,   14355055, 279780, 411905, 4605,15805, 280, 55, 5};
   cu64 sn5_16[] = /*  5 */ {1, 51274805, 245305530,   10582980,2027855,  83530,17030,  980, 205,105, 5};
   cu64 sn5_17[] = /*  5 */ {1,256435955,7356237775ull,83847525, 490230,2064055,22930,79355,1530,455,30,5};
   // From here on the results are from a random search
   cu64 sn5_18[] = /*  6 */ {1ull,3163037327215ull,58283536500ull,10186171490ull,4243419415ull,458856974740ull,924839705ull,194292875ull,17551815ull,870745ull,2311250ull,351900ull};
   cu64 sn5_19[] = /*  6 */ {1ull,18188911134265ull,3200562885465ull,423137816015ull,104741497570ull,23890293485ull,4109631270ull,992577600ull,116734170ull,38865935ull,7200360ull,1245415ull,369405ull};
   cu64 sn5_20[] = /*  7 */ {1ull,46958726785690ull,2228843224785ull,15344033282910ull,459162120820ull,18763043215ull,66712744095ull,215612395ull,5147262630ull,528682095ull,40914380ull,6450955ull,1810625ull};
   cu64 sn5_21[] = /*  7 */ {1ull,103407331808035ull,73906432359645ull,14197366504185ull,3789788783115ull,303588414810ull,20072670785ull,91701111035ull,3920814355ull,956761030ull,193859450ull,29109725ull,1815395ull,6069095ull};
   cu64 sn5_22[] = /*  8 */ {1ull,2325028983837050ull,339609525637745ull,54686045818620ull,10191810138315ull,390237228110ull,3162350456955ull,6810987690ull,35238246310ull,4088834195ull,1120789640ull,153790185ull,518435ull,15781435ull};
   cu64 sn5_23[] = /*  8 */ {1ull,5978899564948140ull,645555181735885ull,146616157931560ull,16010720207070ull,59818870960225ull,1442997770420ull,352079040540ull,112070355925ull,29050108975ull,1876119220ull,1071138555ull,11275455ull,8521300ull,98719515ull};
   cu64 sn5_24[] = /*  9 */ {1ull,18094905600096760ull,742298092644455ull,7311812839032355ull,60105264403240ull,194558147301940ull,8487871366095ull,578605692970ull,1618880704895ull,104547451740ull,4395543140ull,24887819600ull,80104200ull,763076640ull,43732110ull};
   cu64 sn5_25[] = /* 10 */ {1ull,84751088663821040ull,26801030544636345ull,9836245401792775ull,1231645253023925ull,346246883670360ull,61416730607525ull,2515375409670ull,15631717761340ull,521906283245ull,120571492815ull,4975750870ull,13232226170ull,1218162080ull,115215025ull};

#undef  HINTLIB_SNM
#define HINTLIB_SNM(m) HINTLIB_SNBM(5,m)
   const SNData snets5 [] =
   {
      { 0, 0 },
      HINTLIB_SNM(1)
      HINTLIB_SNM(2)
      HINTLIB_SNM(3)
      HINTLIB_SNM(4)
      HINTLIB_SNM(5)
      HINTLIB_SNM(6)
      HINTLIB_SNM(7)
      HINTLIB_SNM(8)
      HINTLIB_SNM(9)
      HINTLIB_SNM(10)
      HINTLIB_SNM(11)
      HINTLIB_SNM(12)
      HINTLIB_SNM(13)
      HINTLIB_SNM(14)
      HINTLIB_SNM(15)
      HINTLIB_SNM(16)
      HINTLIB_SNM(17)
      HINTLIB_SNM(18)
      HINTLIB_SNM(19)
      HINTLIB_SNM(20)
      HINTLIB_SNM(21)
      HINTLIB_SNM(22)
      HINTLIB_SNM(23)
      HINTLIB_SNM(24)
      HINTLIB_SNM(25)
   };


   /*
    *  b = 7
    *
    *  Optimality proven for m <= 14
    *
    *  The values have been confirmed by test_shiftnet up to m <= 20,
    *  i.e., for all m.
    */

   //     b m          t
   cu64 sn7_1 [] = /*  0 */ {1};
   cu64 sn7_2 [] = /*  0 */ {1,         7};
   cu64 sn7_3 [] = /*  0 */ {1,        56,           7};
   cu64 sn7_4 [] = /*  0 */ {1,       742,         105,           7};
   cu64 sn7_5 [] = /*  0 */ {1,      3192,         497,          56,        7};
   cu64 sn7_6 [] = /*  0 */ {1,     22547,       10052,         987,       56,       7};
   cu64 sn7_7 [] = /*  0 */ {1,    763035,      108542,       14168,      987,     154,      7};
   cu64 sn7_8 [] = /*  1 */ {1,    137795,       44303,      826294,      399,     105,      7};
   cu64 sn7_9 [] = /*  1 */ {1,    964621,    17595221,       21567,     6034,     448,    203,     7};
   cu64 sn7_10[] = /*  1 */ {1,  48858789,    34126302,      227563,  2491902,    6818,    840,   105,     7};
   cu64 sn7_11[] = /*  1 */ {1,  61764262,   579706407,    31566297,  1188257,  229768,   8484,   399,   105,  7};
   cu64 sn7_12[] = /*  1 */ {1,2487999850u,  957498325,    11239676,246221913, 1479905, 267498, 26565,   791,154, 7};
   cu64 sn7_13[] = /*  1 */ {1,3307531955u,69918946524ull,775266443,159430082,31916108,3958612,235011,106379,448,56,7};
   cu64 sn7_14[] = /*  3 */ {1, 330413426, 23833203058ull, 17092579, 40547115,  274995, 827813,  5201,  1477,203, 7};
   cu64 sn7_15[] = /*  4 */ {1, 330405684,  2019198230,   146524217, 10032414, 1858080,  33474,  3486,   791, 56, 7};
   // From here on the results are from a random search
   cu64 sn7_16[] = /*  4 */ {1ull,20070618966839ull,3695363667261ull,585482814588ull,31105501385ull,8459248287ull,1717296196ull,139782342ull,289506ull,16537206ull,4994801ull,20258ull};
   cu64 sn7_17[] = /*  5 */ {1ull,75901594765066ull,16544166645139ull,2365334504638ull,578712107960ull,85470377062ull,11049706535ull,566100542ull,101488457ull,20049813ull,797181ull,1745387ull};
   cu64 sn7_18[] = /*  5 */ {1ull,1178502960295066ull,89700433808783ull,11601266684875ull,2316368614959ull,489500959339ull,90881318864ull,11035685248ull,326855585ull,152531421ull,39545688ull,3633133ull,627746ull};
   cu64 sn7_19[] = /*  6 */ {1ull,4613553869755621ull,1104802219166989ull,174363747574030ull,26325503718570ull,514044682074ull,15321275900ull,6926772279ull,2034958454725ull,185542581ull,574426034ull,1735909ull,23420663ull};
   cu64 sn7_20[] = /*  6 */ {1ull,61781329968928934ull,5922621996449130ull,7520584469531ull,568175782001800ull,166380597908238ull,4445232361766ull,57461780299ull,490406160706ull,3255398867ull,1964489331ull,211033438ull,6881812ull,2226686ull};

#undef  HINTLIB_SNM
#define HINTLIB_SNM(m) HINTLIB_SNBM(7,m)
   const SNData snets7 [] =
   {
      { 0, 0 },
      HINTLIB_SNM(1)
      HINTLIB_SNM(2)
      HINTLIB_SNM(3)
      HINTLIB_SNM(4)
      HINTLIB_SNM(5)
      HINTLIB_SNM(6)
      HINTLIB_SNM(7)
      HINTLIB_SNM(8)
      HINTLIB_SNM(9)
      HINTLIB_SNM(10)
      HINTLIB_SNM(11)
      HINTLIB_SNM(12)
      HINTLIB_SNM(13)
      HINTLIB_SNM(14)
      HINTLIB_SNM(15)
      HINTLIB_SNM(16)
      HINTLIB_SNM(17)
      HINTLIB_SNM(18)
      HINTLIB_SNM(19)
      HINTLIB_SNM(20)
   };


   /*
    *  b = 8
    *
    *  Optimality proven for m <= 14
    *
    *  The values have been confirmed by test_shiftnet up to m <= 17,
    *  i.e., for all m.
    */

   //     b m          t
   cu64 sn8_1 [] = /*  0 */ {1};
   cu64 sn8_2 [] = /*  0 */ {1,           8};
   cu64 sn8_3 [] = /*  0 */ {1,          72,              8};
   cu64 sn8_4 [] = /*  0 */ {1,        1096,            136,             8};
   cu64 sn8_5 [] = /*  0 */ {1,        5256,            712,            72,            8};
   cu64 sn8_6 [] = /*  0 */ {1,       42440,          26696,          1736,           72,           8};
   cu64 sn8_7 [] = /*  0 */ {1,      341640,         121928,          5000,         2312,          72,      8};
   cu64 sn8_8 [] = /*  1 */ {1,      300296,          78664,          6600,          776,          72,      8};
   cu64 sn8_9 [] = /*  0 */ {1,    24792264,        9042376,       1217736,       108360,       31304,    904,     200,      8};
   cu64 sn8_10[] = /*  1 */ {1,    43906888,      946942528,      12046856,       189832,        4744,   1224,      72,      8};
   cu64 sn8_11[] = /*  2 */ {1,    19208072,        5030664,     134274760,        75592,       13256,   1416,     136,      8};
   cu64 sn8_12[] = /*  1 */ {1,  1447352136,    60741297216ull, 1009705480,     38821832,    11915848,1095176,    7432,    712,    72,  8};
   cu64 sn8_13[] = /*  1 */ {1, 34135165896ull,208126646848ull, 5750870536ull, 767327432,    22280328,2433600,  715528,  20296,   712, 72, 8};
   cu64 sn8_14[] = /*  1 */ {1,744054312072ull,104965919944ull,43410671240ull,5808415560ull,222732552,7904008,17877568,1172360,106696,904,72,8};
   // From here on the results are from a random search
   cu64 sn8_15[] = /*  4 */ {1ull,26717256822192ull,2068210416208ull,470087866968ull,48582330488ull,1434428648ull,449901200ull,128717504ull,6905320ull,1440584ull,140080ull};
   cu64 sn8_16[] = /*  4 */ {1ull,105514649966584ull,3666470442136ull,4746695639816ull,175917574688ull,44628974880ull,2344228848ull,766610032ull,81677456ull,8596200ull,133752ull,2616ull};
   cu64 sn8_17[] = /*  5 */ {1ull,1927941266758312ull,236138896714080ull,9053124369216ull,1848132093728ull,296765829640ull,41277655816ull,2144286408ull,588405888ull,102881008ull,14089080ull,1986376ull};

#undef  HINTLIB_SNM
#define HINTLIB_SNM(m) HINTLIB_SNBM(8,m)
   const SNData snets8 [] =
   {
      { 0, 0 },
      HINTLIB_SNM(1)
      HINTLIB_SNM(2)
      HINTLIB_SNM(3)
      HINTLIB_SNM(4)
      HINTLIB_SNM(5)
      HINTLIB_SNM(6)
      HINTLIB_SNM(7)
      HINTLIB_SNM(8)
      HINTLIB_SNM(9)
      HINTLIB_SNM(10)
      HINTLIB_SNM(11)
      HINTLIB_SNM(12)
      HINTLIB_SNM(13)
      HINTLIB_SNM(14)
      HINTLIB_SNM(15)
      HINTLIB_SNM(16)
      HINTLIB_SNM(17)
   };


   /*
    *  b = 9
    *
    *  Optimality proven for m <= 15  (i.e., for all m)
    *
    *  The values have been confirmed by test_shiftnet up to m <= 15,
    *  i.e., for all m.
    */

   //     b m          t
   cu64 sn9_1 [] = /*  0 */ {1};
   cu64 sn9_2 [] = /*  0 */ {1,            9};
   cu64 sn9_3 [] = /*  0 */ {1,           90,                 9};
   cu64 sn9_4 [] = /*  0 */ {1,         2277,               171,              9};
   cu64 sn9_5 [] = /*  0 */ {1,         8352,               819,            171,              9};
   cu64 sn9_6 [] = /*  0 */ {1,        81576,             43101,           1872,             90,             9};
   cu64 sn9_7 [] = /*  0 */ {1,       607023,            151722,          15885,           4140,            90,            9};
   cu64 sn9_8 [] = /*  0 */ {1,     20173383,           2293605,         369936,          25929,          3006,          252,        9};
   cu64 sn9_9 [] = /*  1 */ {1,      5389263,           1412001,         106524,          17424,          1791,           90,        9};
   cu64 sn9_10[] = /*  1 */ {1,     69705693,         794334528,       10206009,         128799,         13941,         1710,       90,       9};
   cu64 sn9_11[] = /*  1 */ {1,    458275734,       24488544753ull,   292560102,       31392126,        189468,        38727,     1791,     252,      9};
   cu64 sn9_12[] = /*  1 */ {1, 100592277012ull,    11724155469ull,  2886056496ull,    92972781,      26907399,      4223187,     8433,     819,    171,     9};
   cu64 sn9_13[] = /*  1 */ {1, 100628135631ull,  1998701214435ull, 29646653436ull,   839167299,     207580410,      7002126,  3993390,   11997,    819,   252,  9};
   cu64 sn9_14[] = /*  1 */ {1,5130930877983ull,  1186123550958ull,182687248449ull, 12725655183ull,  781952220,    350898408, 42038199, 3424203,  61893,   819,252, 9};
   cu64 sn9_15[] = /*  1 */ {1,9592049562741ull,115475165576424ull,682311915261ull,194082825798ull,28149470496ull,1795925691,139260474,41596668,1110357,167922,900,90,9};

#undef  HINTLIB_SNM
#define HINTLIB_SNM(m) HINTLIB_SNBM(9,m)
   const SNData snets9 [] =
   {
      { 0, 0 },
      HINTLIB_SNM(1)
      HINTLIB_SNM(2)
      HINTLIB_SNM(3)
      HINTLIB_SNM(4)
      HINTLIB_SNM(5)
      HINTLIB_SNM(6)
      HINTLIB_SNM(7)
      HINTLIB_SNM(8)
      HINTLIB_SNM(9)
      HINTLIB_SNM(10)
      HINTLIB_SNM(11)
      HINTLIB_SNM(12)
      HINTLIB_SNM(13)
      HINTLIB_SNM(14)
      HINTLIB_SNM(15)
   };

#undef HINTLIB_SNM
#undef HINTLIB_SNBM

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
      { sizeof(snets2) / sizeof(SNData), snets2 },
      { sizeof(snets3) / sizeof(SNData), snets3 },
      { sizeof(snets4) / sizeof(SNData), snets4 },
      { sizeof(snets5) / sizeof(SNData), snets5 },
      { 0, 0 },
      { sizeof(snets7) / sizeof(SNData), snets7 },
      { sizeof(snets8) / sizeof(SNData), snets8 },
      { sizeof(snets9) / sizeof(SNData), snets9 }
   };
}


/**
 *  max Base For Shift Net ()
 */

int L::maxBaseForShiftNet ()
{
   return sizeof (snets) / sizeof (SNLists) - 1;
}


/**
 *  max M For Shift Net ()
 */

int L::maxMForShiftNet (int base)
{
   if (base > maxBaseForShiftNet())  throw InternalError (__FILE__, __LINE__);

   return snets[base].num - 1;
}


/**
 *  optimal Shift Net T ()
 */

int L::optimalShiftNetT (int base, int m)
{
   if (m > maxMForShiftNet (base))  throw InternalError (__FILE__, __LINE__);

   return m - snets[base].lists[m].length;
}


/**
 *  init Shift Net ()
 */

void L::Private::initShiftNetFirstDim (GeneratorMatrix &gm)
{
   int b = gm.getBase();
   int m = gm.getM();

   if (gm.getDimension() > m)  throw InternalError (__FILE__, __LINE__);

   if (m > maxMForShiftNet (b))  throw InternalError (__FILE__, __LINE__);

   const SNData& r = snets[b].lists[m];
   const u64* data = r.data;

   for (unsigned i = 0; i < r.length; ++i)
   {
      gm.vSetPackedRowVector (0, i, data[i]);
   }
}


