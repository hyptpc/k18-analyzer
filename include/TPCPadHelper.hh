// -*- C++ -*-


/*
 * Update note:
 *  - Pad / Layer / Row indices are unified to start from 0
 *
 * Author: Haein Lee
 * Date  : 2025-12-26
 */

#ifndef TPC_PAD_HELPER_HH
#define TPC_PAD_HELPER_HH

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <TDirectory.h>
#include <TH2Poly.h>
#include <TMath.h>
#include <TVector3.h>

#include "DetectorID.hh"
#include "Exception.hh"

namespace tpc
{
const Double_t Z_TARGET = -143.; // Target from center

// E72 cylindrical target (axis along y) with a surrounding holder.
//   xz plane: raw target radius = 40 mm, +10 mm holder shell, +10 mm clearance
//   y  axis : raw half-length   = 50 mm, +10 mm clearance
// Use TARGET_RADIUS for radial (cylindrical) cuts in the xz plane, and
// TARGET_HALF_X/Y/Z for axis-aligned (box) cuts.
const Double_t TARGET_RADIUS = 40. + 10. + 10.;  // = 60 mm  (raw + holder + clearance)
const Double_t TARGET_HALF_X = TARGET_RADIUS;    // xz plane is the cylindrical cross section
const Double_t TARGET_HALF_Y = 50. + 10.;        // = 60 mm  (raw + clearance)
const Double_t TARGET_HALF_Z = TARGET_RADIUS;    // xz plane is the cylindrical cross section

enum EPadParameter
{
  kLayerID,
  kNumOfPad,
  kRadius,
  kNumOfDivision,
  kDummy,
  kLength,
  NPadParameter
};

//for kinematic fitting
//E42
static const Double_t POSITION_SCALE = 1.;

//for clustering
//E42
static const Int_t MAX_ROW_DIF_TPC = 2;

// unit is [m/ns], used for p = q B r in the helix momentum formula and so on
static const Double_t C_LIGHT = TMath::C() * 1.e-9;
  
//_____________________________________________________________________________
//#OfPad #division #radius padLength
static const Double_t padParameter[NumOfLayersTPC][NPadParameter] =
{{0, 48,    14.75, 48, 0,  9.},
 {1, 48,    24.25, 48, 0,  9.},
 {2, 72,    33.75, 72, 0,  9.},
 {3, 96,    43.25, 96, 0,  9.},
 {4, 120,    52.75,120,0,   9.},
 {5, 144,    62.25,144,0,   9.},
 {6, 168,    71.75,168,0,   9.},
 {7, 192,    81.25,192,0,   9.},
 {8, 216,    90.75,216,0,   9.},
 {9, 240,    100.25,240,0,  9.},
 {10,208,    111.5,241, 0,  12.5},
 {11,218,    124.5,271, 0,  12.5},
 {12,230,    137.5,300, 0,  12.5},
 {13,214,    150.5,330, 0,  12.5},
 {14,212,    163.5,360, 0,  12.5},
 {15,214,    176.5,390, 0,  12.5},
 {16,220,    189.5,420, 0,  12.5},
 {17,224,    202.5,449, 0,  12.5},
 {18,232,    215.5,479, 0,  12.5},
 {19,238,    228.5,509, 0,  12.5},
 {20,244,    241.5,539, 0,  12.5},
 {21,232,    254.5,569, 0,  12.5},
 {22,218,    267.5,599, 0,  12.5},
 {23,210,    280.5,628, 0,  12.5},
 {24,206,    293.5,658, 0,  12.5},
 {25,202,    306.5,688, 0,  12.5},
 {26,200,    319.5,718, 0,  12.5},
 {27,196,    332.5,748, 0,  12.5},
 {28,178,    345.5,777, 0,  12.5},
 {29,130,    358.5,807, 0,  12.5},
 {30,108,    371.5,837, 0,  12.5},
 {31,90,     384.5,867, 0, 12.5}};

//_____________________________________________________________________________
static const Int_t padOnCenterFrame[] =
{

//Pads on the frame
964,965,966,967,968,969,970,971,972,973,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1176,1177,1178,1179,1180,1181,1182,1183,1184,1185,1186,1187,1188,1189,1190,1191,1192,1193,1194,1195,1196,1197,1198,1199,1200,1201,1202,1203,1204,1205,1206,1207,1208,1209,1210,1211,1212,1236,1237,1238,1239,1240,1241,1242,1243,1244,1245,1246,1247,1248,1249,1250,1251,1252,1253,1254,1255,1256,1257,1258,1259,1260,1261,1262,1263,1264,1265,1266,1267,1268,1269,1270,1271,1394,1395,1396,1397,1398,1399,1400,1404,1405,1406,1407,1408,1409,1410,1411,1412,1413,1414,1415,1416,1417,1418,1419,1420,1421,1422,1423,1424,1425,1426,1427,1428,1429,1430,1431,1435,1436,1437,1438,1439,1440,1441,1454,1455,1456,1457,1458,1459,1460,1464,1465,1466,1467,1468,1469,1470,1471,1472,1473,1474,1475,1476,1477,1478,1479,1480,1481,1482,1483,1484,1485,1486,1487,1488,1489,1490,1491,1495,1496,1497,1498,1499,1500,1501,1594,1595,1596,1597,1603,1604,1605,1606,1607,1646,1647,1648,1649,1650,1651,1657,1658,1659,1662,1663,1664,1670,1671,1672,1673,1674,1675,1714,1715,1716,1717,1718,1724,1725,1726,1727,1806,1807,1808,1809,1810,1811,1812,1813,1814,1815,1816,1817,1878,1879,1880,1881,1888,1889,1890,1891,1952,1953,1954,1955,1956,1957,1958,1959,1960,1961,1962,1963,2017,2018,2019,2024,2025,2026,2027,2099,2100,2101,2102,2103,2110,2111,2112,2113,2114,2186,2187,2188,2189,2194,2195,2196,2218,2219,2220,2221,2222,2223,2224,2225,2226,2227,2308,2309,2310,2316,2317,2322,2323,2329,2330,2331,2412,2413,2414,2415,2416,2417,2418,2419,2420,2421,2426,2427,2428,2429,2517,2518,2519,2520,2521,2522,2523,2524,2525,2526,2539,2540,2541,2542,2543,2544,2545,2546,2547,2548,2636,2637,2638,2639,2731,2732,2738,2739,2760,2761,2767,2768,2949,2950,2951,2952,2953,2954,2955,2956,2957,2986,2987,2988,2989,2990,2991,2992,2993,2994,3173,3174,3175,3176,3177,3178,3179,3180,3181,3218,3219,3220,3221,3222,3223,3224,3225,3226,3404,3405,3406,3407,3408,3409,3410,3411,3412,3457,3458,3459,3460,3461,3462,3463,3464,3465,3641,3642,3643,3644,3645,3646,3647,3648,3649,3702,3703,3704,3705,3706,3707,3708,3709,3710,3876,3877,3878,3879,3880,3881,3882,3883,3944,3945,3946,3947,3948,3949,3950,3951,4097,4098,4104,4173,4179,4180,4307,4308,4309,4310,4311,4312,4313,4314,4315,4390,4391,4392,4393,4394,4395,4396,4397,4398,4511,4512,4519,4602,4609,4610,4712,4713,4714,4715,4716,4717,4718,4719,4810,4811,4812,4813,4814,4815,4816,4817,4909,4910,4911,4912,4913,4914,4915,4916,5015,5016,5017,5018,5019,5020,5021,5022,5103,5104,5107,5110,5111,5216,5217,5220,5223,5224,5286,5287,5288,5289,5290,5291,5292,5293,5294,5407,5408,5409,5410,5411,5412,5413,5414,5415,5440,5441,5442,5443,5444,5565,5566,5567,5568,5569,

//Empty Pads
1016,1017,1393,1401,1402,1403,1432,1433,1434,1461,1462,1463,1492,1493,1494,1502,1598,1599,1600,1601,1602,1603,1652,1653,1654,1655,1656,1665,1666,1667,1668,1669,1719,1720,1721,1722,1723,1876,1877,1882,1883,1884,1885,1886,1887,1892,1893,2020,2021,2022,2023,2104,2105,2106,2107,2108,2109,2190,2191,2192,2193,2311,2312,2313,2314,2315,2324,2325,2326,2327,2328,2733,2734,2735,2736,2737,2762,2763,2764,2765,2766,4099,4100,4101,4102,4103,4174,4175,4176,4177,4178,4513,4514,4515,4516,4517,4518,4603,4604,4605,4606,4607,4608,5105,5106,5108,5109,5218,5219,5221,5222

};

//Pads on the gem supporting frame (No active area)
//_____________________________________________________________________________
static const Int_t deadChannel[] =
{
  72,134,135,222,332,408,467,626,809,59,110,185,284,555,726,921,1140,1363,1564,1565,1778,1818,1819,2036,2037,2245,2246,2247,2455,2456,2668,2669,2887,3111,3112,3342,3343,3579,3580,3814,4035,1629,1850,2068,2277,2486,2487,2699,2700,2918,3142,3374,3611,3845,4066,4277,4481,4681,4879,5073
};

//============================================================
// Pads affected by section-frame-related noise
// Noise pad information taken from E42 data
// Used for the 2026 April beamtime LSB threshold setup
//============================================================
//_____________________________________________________________________________
static const Int_t padOnSectionFrame_E42[] =
{
//Gem section 1 frame 0
1503 ,1504, 1505 ,1506 ,1507 ,1508 ,1509 ,1510 ,1511 ,1512 ,1513 ,1514 ,1515 ,1728 ,1729 ,1730,

//Gem section 1 frame 1
72 ,134 ,135 ,222 ,332 , 408,  467 ,626 ,809 ,1016, 59 ,110 ,185 ,284,555 ,726, 921 ,1140 ,1363, 1564 ,1565 ,1778,

//Gem section 2 frame 0, frame 1 + adjacent pads
1628,1629,1630,1817,1818,1819,1820,1849,1850,1851,2035,2036,2037,2038,2067,2068,2069,2244,2245,2246,2247,2248,2276,2277,2278,2454,2455,2456,2457,2485,2486,2487,2488,2667,2668,2669,2670,2698,2699,2700,2701,2886,2887,2888,2917,2918,2919,3110,3111,3112,3113,3141,3142,3143,3341,3342,3343,3344,3373,3374,3375,3578,3579,3580,3581,3610,3611,3612,3813,3814,3815,3844,3845,3846,4034,4035,4036,4065,4066,4067,4276,4277,4278,4480,4481,4482,4680,4681,4682,4878,4879,4880,5072,5073,5074,

//Gem section 3 frame 0
4720 ,4721 ,4722 ,4723 ,4724 ,4725 ,4726 ,4727 ,4728 ,4929 ,4930 ,4931 ,4932 ,4933 ,4934 ,4935 ,4936 ,5134 ,5135 ,5136 ,5137 ,5138 ,5139 ,5140 ,5328 ,5329 ,5330 ,5331 ,5332 ,5333 ,5487 ,5488 ,5489 ,5490 ,5491 ,5492 ,5612 ,5613 ,5614 ,5615 ,5616 ,5716 ,5717 ,5718 ,5719 ,5720, 5532 ,5533, 5612 ,5613 ,5614 ,5615 ,5616 ,5654 ,5655 ,5716 ,5717 ,5718 ,5719 ,5720 ,5757 ,5758,

//Gem section 3 frame 1
3413 ,3414 ,3415 ,3416 ,3417 ,3418 ,3419 ,3660 ,3661 ,3662 ,3663 ,3664 ,3665 ,3904 ,3905 ,3906 ,3907 ,3908 ,4134 ,4135 ,4136 ,4137 ,4353 ,4354 ,4355 ,4356 ,4565 ,4566 ,4567 ,4568 ,4774 ,4775 ,4776 ,4981 ,4982 ,5183 ,5184 ,5372 ,5373 ,5374 ,5530 ,5531 ,5532 ,5533 ,5654 ,5655 ,5757 ,5758,

//Gem section 4 frame 0
4409 ,4410 ,4411 ,4412 ,4413 ,4414 ,4415 ,4416 ,4417 ,4418 ,4419 ,4420 ,4421 ,4422 ,4423 ,4424 ,4425 ,4426 ,4427 ,4428 ,4429 ,4430 ,4431 ,4432 ,4433 ,4434 ,4435 ,4436 ,4437 ,4438 ,4439 ,4440 ,4441 ,4442 ,4443 ,4444 ,4445 ,4446 ,4447 ,4448 ,4449 ,4450 ,4451 ,4452 ,4611 ,4612 ,4613 ,4614 ,4615 ,4616 ,4617 ,4618 ,4619,

//Gem section 4 frame 1
2794 ,2795 ,2796 ,2797 ,2798 ,2799 ,2800 ,2801 ,2802 ,2803 ,2804 ,2805 ,2806 ,2807 ,2808 ,2809 ,3001 ,3002 ,3003 ,3004 ,3005 ,3006 ,3007 ,3008 ,3009 ,3010 ,3011 ,3012 ,3013 ,3014 ,3015 ,3016 ,3017, 3018, 3037, 3038, 3039 ,3040 ,3041 ,3042 ,3043 ,3044 ,3045 ,3046 ,3047 ,3048 ,3049 ,3050 ,3051 ,3052 ,3053 ,3054 ,3227 ,3228 ,3229 ,3230 ,3231 ,3289 ,3290 ,3291 ,3292 ,3293 ,3294 ,3295 ,3296 ,3297 ,3540 ,3541 ,3542 ,3543 ,3544 ,3545 ,3546
};

//=======================================================
// Pads affected by section-frame-related noise
// Selected based on hit-pattern studies using unbiased data
// run03396 (HS On, clk trigger, Threshold = 0)
//=======================================================
//_____________________________________________________________________________
static const Int_t padOnSectionFrame_E72[] =
{
//section1 left
1508,1509,1510,1511,1512,
1287,1288,1289,1290,
1057,1058,1059,1060,1061,1062,1063,
854,855,856,857,858,859,860,861,862,863,864,865,866,867,868,869,870,871,872,873,874,
1090,1091,1092,1093,1094,1095,1096,
1336,1337,1338,1339,

//section1 right
812,
629,
470,
335,
224,
137,
74,
23,24,
11,
57,
108,109,
183,184,
283,
406,
553,
724,
919,
1138,
1361,
1562,1563,
1775,1776,

//section1 right edge
1450,1451,1452,1453,
1442,1443,1444,1445,1446,1447,
1660,

//section2 left
1818,
2034,2035,2036,2037,
2244,2245,2246,
2453,2454,2455,
2667,2668,
2885,2886,2887,
3110,3111,
3341,3342,
3578,3579,
3812,3813,3814,
4034,4035,

//section2 right
1630,
1850,1851,
2069,
2278,
2487,
2700,2701,
2919,
3143,
3374,
3611,3612,
3846,
4067,
4277,4278,
4482,
4682,
4879,
5073,5074,

//section3 left
3414,3415,3416,3417,3418,3419,3420,3421,
3661,3662,3663,3664,3665,3666,3667,
3905,3906,3907,3908,3909,3910,
4135,4136,4137,4138,4139,
4354,4355,4356,4357,4358,
4567,4568,4569,4570,
4775,4776,4777,4778,
4981,4982,4983,
5183,5184,5185,
5374,5375,5376,
5532,5533,5534,
5655,5656,5657,
5758,5759,5760,

//section3 right
4720,4721,4722,4723,4724,4725,4726,
4927,4928,4929,4930,4931,4932,4933,4934,4935,
5132,5133,5134,5135,5136,5137,5138,5139,
5326,5327,5328,5329,5330,5331,5332,
5486,5487,5488,5489,5490,5491,
5611,5612,5613,5614,5615,
5715,5716,5717,5718,5719,

//section3 bottom edge
4105,4106,
5112,5113,

//section3 top edge
4171,4172,
5214,5215,
5404,5405,5406,
5560,5561,5562,5563,5564,

//section4 left
3229,3230,
2999,3000,3001,3002,3003,3004,3005,3006,3007,3008,3009,3010,3011,3012,
2789,2790,2791,2792,2793,2794,2795,2796,2797,2798,2799,2800,2801,2802,2803,2804,2805,2806,2807,2808,2809,2810,2811,2812,2813,2814,2815,
3042,3043,3044,3045,3046,3047,3048,3049,3050,3051,3052,3053,3054,3055,3056,
3292,3293,3294,3295,3296,3297,3298,3299,3300,3301,
3543,3544,3545,3546,3547,3548,3549,
3796,3797,

//section4 right
4611,4612,4613,4614,4615,4616,4617,4618,4619,4620,4621,4622,4623,4624,4625,
4414,4415,4416,4417,4418,4419,4420,4421,4422,4423,4424,4425,4426,4427,4428,4429,4430,4431,4432,4433,4434,4435,4436,4437,4438,4439,4440,4441,4442,4443,4444,4445,4446,4447,4448,
4661,4662,4663
};

//=======================================================
// Pads showing abnormal waveform shapes
// Identified from waveform inspection
// run2278 (HS On, K- 735 MeV/c, physics trigger)
//=======================================================
//_____________________________________________________________________________

static const Int_t padAbnormalWaveform_E72[] =
{
//section1 left -> No


//section1 right
812,
629,
470,
335,
224,
137,
74,
57,58,
108,109,
183,184,
283,
406,
553,
724,
919,
1138,
1361,
1562,1563,
1776,

//section1 right edge
1450,1451,1452,1453,

//section2 left
1818,
2034,2035,2036,
2244,2245,2246,
2453,2454,2455,
2667,2668,
2885,2886,2887,
3110,3111,
3341,3342,
3578,3579,
3813,
4034,4035,

//section2 right
1630,
1850,1851,
2069,
2278,
2487,
2700,2701,
2919,
3143,
3374,
3612,
3846,
4067,
4277,4278,
4482,
4682,
4879,
5073,5074,

//section3 left -> No

//section3 right -> No

//section3 bottom edge
4105,4106,
5112,5113,

//section3 top edge
4171,4172,
5214,5215

//section4 left -> No

//section4 right -> No
};


//=======================================================
// Pads affected by section-frame-related noise
// Selected from Tohoku 2024 Jan cosmic test data
// Used for the 2025 Nov beamtime LSB threshold setup
// (run2556 - run3044)
//=======================================================
//_____________________________________________________________________________
static const Int_t padOnSectionFrame_Tohoku[] =
{
56,57,58,68,72,73,74,108,109,137,183,224,282,335,405,470,552,553,629,723,812,854,855,856,857,858,859,860,861,862,863,864,918,919,1042,1058,1059,1060,1061,1062,1063,1094,1137,1138,1288,1289,1290,1291,1361,1384,1386,1387,1388,1390,1509,1510,1511,1512,1513,1562,1563,1586,1587,1588,1589,1590,1592,1608,1630,1706,1729,1776,1850,2034,2035,2037,2069,2217,2242,2243,2244,2245,2246,2278,2453,2454,2455,2467,2469,2470,2483,2487,2667,2701,2790,2792,2794,2795,2796,2797,2798,2799,2800,2801,2802,2804,2805,2806,2807,2808,2809,2810,2811,2812,2813,2814,2885,2887,2919,2999,3000,3001,3002,3003,3004,3005,3006,3007,3008,3009,3010,3011,3012,3013,3045,3050,3051,3052,3053,3054,3055,3056,3086,3110,3111,3143,3227,3228,3291,3292,3293,3294,3295,3296,3297,3298,3299,3341,3342,3344,3374,3375,3414,3415,3416,3417,3418,3419,3420,3449,3451,3452,3467,3473,3475,3478,3479,3480,3481,3482,3525,3526,3527,3528,3541,3542,3543,3544,3545,3546,3547,3565,3566,3568,3569,3570,3578,3580,3612,3661,3662,3663,3664,3665,3666,3667,3668,3796,3797,3814,3846,3905,3906,3907,3908,3909,4034,4035,4067,4105,4106,4135,4136,4137,4138,4139,4171,4172,4277,4278,4354,4355,4356,4357,4385,4386,4387,4415,4416,4417,4418,4419,4420,4421,4422,4423,4424,4425,4426,4427,4428,4429,4430,4431,4432,4433,4434,4435,4436,4437,4438,4439,4440,4441,4442,4443,4444,4445,4446,4447,4482,4611,4612,4613,4614,4615,4616,4617,4618,4619,4620,4621,4622,4623,4624,4682,4720,4721,4722,4723,4724,4725,4726,4775,4776,4777,4778,4879,4927,4928,4929,4930,4931,4932,4934,4935,4945,4981,4982,4983,5073,5074,5102,5112,5113,5132,5133,5134,5135,5136,5137,5138,5139,5183,5184,5185,5214,5215,5297,5375,5445,5448,5449,5488,5533,5562,5572,5611,5612,5613,5614,5655,5656,5657,5715,5716,5717,5718,5740,5741,5758,5759,5760
};


//_____________________________________________________________________________
static const Int_t frameLowEdge[NumOfLayersTPC][5] =
{
  //dummy value : -100
  //Layer0~7 : no pad on the frame
  //Layer10~31(Outer layers) : the last pad treated as a pad on the frame

  {-100,-100,-100,-100,-100}, //0
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100}, //5
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {75,127,-100,-100,-100},
  {71,131,-100,-100,-100},
  {48,109,207,-100,-100}, //10
  {41,93,109,161,217},
  {35,105,181,229,-100},
  {16,98,185,213,-100},
  {3,93,107,197,211},
  {90,112,209,-100,-100}, //15
  {90,119,219,-100,-100},
  {88,125,223,-100,-100},
  {88,133,231,-100,-100},
  {87,140,237,-100,-100},
  {86,147,243,-100,-100}, //20
  {77,145,231,-100,-100},
  {66,142,217,-100,-100},
  {58,141,209,-100,-100},
  {52,143,205,-100,-100},
  {47,145,201,-100,-100}, //25
  {42,148,199,-100,-100},
  {36,149,195,-100,-100},
  {23,144,177,-100,-100},
  {124,-100,-100,-100,-100},
  {107,-100,-100,-100,-100}, //30
  {89,-100,-100,-100,-100}
};

//_____________________________________________________________________________
static const Int_t frameHighEdge[NumOfLayersTPC][5] =
{
  //dummy value : -100
  //Layer0~7 : no pad on the frame
  //Layer10~31(Outer layers) : the first pad treated as a pad on the frame

  {-100,-100,-100,-100,-100}, //0
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100}, //5
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {86,140,-100,-100,-100},
  {108,168,-100,-100,-100},
  {0,98,159,-100,-100}, //10
  {0,56,108,124,176},
  {0,48,124,194,-100},
  {0,28,115,197,-100},
  {0,14,104,118,208},
  {4,101,123,-100,-100}, //15
  {0,100,129,-100,-100},
  {0,98,135,-100,-100},
  {0,98,143,-100,-100},
  {0,97,150,-100,-100},
  {0,96,157,-100,-100}, //20
  {0,86,154,-100,-100},
  {0,75,151,-100,-100},
  {0,68,151,-100,-100},
  {0,62,153,-100,-100},
  {0,56,154,-100,-100}, //25
  {0,51,157,-100,-100},
  {0,46,159,-100,-100},
  {0,33,154,-100,-100},
  {5,-100,-100,-100,-100},
  {0,-100,-100,-100,-100}, //30
  {0,-100,-100,-100,-100}
};

//_____________________________________________________________________________
//E42
static const Double_t clusterSizeInner[2][10] ={
  {0.12912,0.12912,0.129121,0.164862,0.211372,0.256798,0.33125,0.315186,0.27439,0.303371},//proton
  {0.274143,0.550639,0.683333,0.683333,0.683333,0.683333,0.683333,0.683333,0.683333,0.683333}//pion
  //0.0 - 0.1 , 0.1 - 0.2, ... 0.9 - 1.0 GeV/c
};

//_____________________________________________________________________________
//E42
static const Double_t clusterSizeOuter[2][10] ={
  {0.285714,0.285714,0.326683,0.280962,0.243775,0.234132,0.217313,0.256798,0.293478,0.285319},//proton
  {0.257669,0.348881,0.308642,0.308642,0.308642,0.308642,0.308642,0.308642,0.308642,0.308642}//pion
  //0.0 - 0.1 , 0.1 - 0.2, ... 0.9 - 1.0 GeV/c
};

//_____________________________________________________________________________
// Check validity of Layer ID
inline void ValidateLayer(Int_t layer, const char* func_name)
{
  if (layer < 0 || NumOfLayersTPC <= layer) {
    throw Exception(Form("[tpc::%s] Invalid layerID %d (Max: %d)", 
                         func_name, layer, NumOfLayersTPC - 1));
  }
}

//_____________________________________________________________________________
// Check validity of Row ID for a specific Layer
inline void ValidateRow(Int_t layer, Int_t row, const char* func_name)
{
  // Ensure layer is valid before accessing array to avoid segfault
  ValidateLayer(layer, func_name);

  Int_t max_row = static_cast<Int_t>(padParameter[layer][kNumOfPad]);
  if (row < 0 || max_row <= row) {
    throw Exception(Form("[tpc::%s] Invalid rowID %d for layer %d (Max: %d)", 
                         func_name, row, layer, max_row - 1));
  }
}

//_____________________________________________________________________________
inline Bool_t IsCircularLayer(Int_t layer)
{
#ifdef PAD_HELPER_DEBUG
  ValidateLayer(layer, __func__);
#endif
  // 0-origin layer index: 0-9 are full-ring layers.
  return (layer < 10);
}

//_____________________________________________________________________________
inline Bool_t TryResolveMRow(Int_t layer, Double_t m_row, Double_t& resolved_m_row)
{
#ifdef PAD_HELPER_DEBUG
  ValidateLayer(layer, __func__);
#endif
  if (!std::isfinite(m_row)) return false;

  const Double_t n_pad = padParameter[layer][kNumOfPad];
  const Double_t n_div = padParameter[layer][kNumOfDivision];
  constexpr Double_t kRoundHalfWidth = 0.5;
  constexpr Double_t kMRowTolerance = 1.e-9;

  if (IsCircularLayer(layer)) {
    Double_t wrapped = std::fmod(m_row, n_pad);
    if (wrapped < 0.) wrapped += n_pad;
    resolved_m_row = wrapped;
    return true;
  }

  // Sector layers: try equivalent one-turn-shifted representations, then
  // accept only values that can round to a physical row index.
  const Double_t min_m_row = -kRoundHalfWidth - kMRowTolerance;
  const Double_t max_m_row = n_pad - kRoundHalfWidth + kMRowTolerance;
  const Double_t candidates[3] = {m_row, m_row + n_div, m_row - n_div};
  for (Double_t candidate : candidates) {
    if (candidate < min_m_row || max_m_row <= candidate) continue;
    if (candidate < 0.) {
      resolved_m_row = 0.;
    } else if (candidate > n_pad - 1.) {
      resolved_m_row = n_pad - 1.;
    } else {
      resolved_m_row = candidate;
    }
    return true;
  }
  return false;
}

//_____________________________________________________________________________
inline Bool_t TryResolveRow(Int_t layer, Double_t m_row, Int_t& resolved_row)
{
  Double_t resolved_m_row = TMath::QuietNaN();
  if (!TryResolveMRow(layer, m_row, resolved_m_row)) return false;

  const Int_t n_pad = static_cast<Int_t>(padParameter[layer][kNumOfPad]);
  Int_t row = TMath::Nint(resolved_m_row);
  if (row < 0 || n_pad <= row) {
    if (!IsCircularLayer(layer)) return false;
    row %= n_pad;
    if (row < 0) row += n_pad;
  }

  resolved_row = row;
  return true;
}

//_____________________________________________________________________________
// Check validity of Pad ID (Global ID)
inline void ValidatePadID(Int_t pad_id, const char* func_name)
{
  if (pad_id < 0 || NumOfPadTPC <= pad_id) {
    throw Exception(Form("[tpc::%s] Invalid padID %d (Max: %d)", 
                         func_name, pad_id, NumOfPadTPC - 1));
  }
}

//_____________________________________________________________________________
inline Int_t GetAGETId(Int_t asad, Int_t layer, Int_t row)
{
#ifdef PAD_HELPER_DEBUG
  ValidateRow(layer, row, __func__);
#endif
  Int_t flag=-1;
  switch(asad){
  case 0:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==1 && row==22-ch*2)      flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==46-(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==28-(ch-23)*2) flag=0;
    if(layer==1 && row==0) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==1 && row==23-ch*2)      flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==47-(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==29-(ch-23)*2) flag=1;
    if(layer==1 && row==1) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==0 && row==22-ch*2)      flag=2;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==34-(ch-13)*2) flag=2;
    for(Int_t ch=23; ch<=31; ++ch) if(layer==2 && row==16-(ch-23)*2) flag=2;
    if(layer==0 && row==0) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==0 && row==23-ch*2)      flag=3;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==35-(ch-13)*2) flag=3;
    for(Int_t ch=23; ch<=31; ++ch) if(layer==2 && row==17-(ch-23)*2) flag=3;
    if(layer==0 && row==1) flag=3;
    return flag;
  case 1:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==1 && row==47-ch*2)      flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==95-(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==77-(ch-23)*2) flag=0;
    if(layer==1 && row==25) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==1 && row==46-ch*2)      flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==94-(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==76-(ch-23)*2) flag=1;
    if(layer==1 && row==24) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==0 && row==47-ch*2)      flag=2;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==71-(ch-13)*2) flag=2;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==2 && row==53-(ch-23)*2) flag=2;
    if(layer==0 && row==25) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==0 && row==46-ch*2)      flag=3;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==70-(ch-13)*2) flag=3;
    for(Int_t ch=23; ch<=31; ++ch) if(layer==2 && row==52-(ch-23)*2) flag=3;
    if(layer==0 && row==24) flag=3;
    return flag;
  case 2:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==5 && row==22-ch*2)       flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==5 && row==143-(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==125-(ch-23)*2) flag=0;
    if(layer==5 && row==0) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==5 && row==23-ch*2)       flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==5 && row==142-(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==124-(ch-23)*2) flag=1;
    if(layer==5 && row==1) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<= 9; ++ch) if(layer==4 && row==18-ch*2)       flag=2;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==117-(ch-12)*2) flag=2;
    if(layer==4 && row==119) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<= 9; ++ch) if(layer==4 && row==19-ch*2)       flag=3;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==116-(ch-12)*2) flag=3;
    if(layer==4 && row==118) flag=3;
    return flag;
  case 3:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==5 && row==24+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==46+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==66+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==5 && row==25+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==47+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==67+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==4 && row==20+ch*2)      flag=2;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==42+(ch-12)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==4 && row==21+ch*2)      flag=3;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==43+(ch-12)*2) flag=3;
    return flag;
  case 4:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==5 && row==73+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==95+(ch-12)*2)  flag=0;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==115+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==5 && row==72+ch*2)       flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==94+(ch-12)*2)  flag=1;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==114+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==4 && row==61+ch*2)      flag=2;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==83+(ch-12)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==4 && row==60+ch*2)      flag=3;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==82+(ch-12)*2) flag=3;
    return flag;
  case 5:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==7 && row==30-ch*2)       flag=0;
    for(Int_t ch=12; ch<=16; ++ch) if(layer==7 && row==8-(ch-12)*2)   flag=0;
    for(Int_t ch=17; ch<=21; ++ch) if(layer==7 && row==191-(ch-17)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==181-(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==7 && row==31-ch*2)       flag=1;
    for(Int_t ch=12; ch<=16; ++ch) if(layer==7 && row==9-(ch-12)*2)   flag=1;
    for(Int_t ch=17; ch<=21; ++ch) if(layer==7 && row==190-(ch-17)*2) flag=1;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==180-(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==6 && row==26-ch*2)       flag=2;
    for(Int_t ch=12; ch<=14; ++ch) if(layer==6 && row==4-(ch-12)*2)   flag=2;
    for(Int_t ch=15; ch<=21; ++ch) if(layer==6 && row==167-(ch-15)*2) flag=2;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==153-(ch-23)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==6 && row==27-ch*2)       flag=3;
    for(Int_t ch=12; ch<=14; ++ch) if(layer==6 && row==5-(ch-12)*2)   flag=3;
    for(Int_t ch=15; ch<=21; ++ch) if(layer==6 && row==166-(ch-15)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==152-(ch-23)*2) flag=3;
    return flag;
  case 6:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==7 && row==32+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==54+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==74+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==7 && row==33+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==55+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==75+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==6 && row==28+ch*2)      flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==50+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==70+(ch-23)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==6 && row==29+ch*2)      flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==51+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==71+(ch-23)*2) flag=3;
    return flag;
  case 7:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==7 && row==97+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==119+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==139+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==7 && row==96+ch*2)       flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==118+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==138+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==6 && row==85+ch*2)       flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==107+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==127+(ch-23)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==6 && row==84+ch*2)       flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==106+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==126+(ch-23)*2) flag=3;
    return flag;
  case 8:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==9 && row==ch*2)          flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==22+(ch-12)*2)  flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==42+(ch-23)*2)  flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==86+(ch-46)*2)  flag=0;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==106+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==9 && row==1+ch*2)        flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==23+(ch-12)*2)  flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==43+(ch-23)*2)  flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==87+(ch-46)*2)  flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==107+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==8 && row==ch*2)         flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==86+(ch-46)*2) flag=2;
    if(layer==8 && row==106) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==8 && row==1+ch*2)       flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==87+(ch-46)*2) flag=3;
    if(layer==8 && row==107) flag=3;
    return flag;
  case 9:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==9 && row==121+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==143+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==163+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==207+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==227+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==9 && row==120+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==142+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==162+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==206+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==226+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==8 && row==109+ch*2)      flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==131+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==151+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==195+(ch-46)*2) flag=2;
    if(layer==8 && row==215) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==8 && row==108+ch*2)      flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==130+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==150+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==194+(ch-46)*2) flag=3;
    if(layer==8 && row==214) flag=3;
    return flag;
  case 10:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==11 && row==1+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==11 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==11 && row==87+(ch-46)*2) flag=0;
    if(layer==11 && row==107) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==11 && row==ch*2)          flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==22+(ch-12)*2)  flag=1;
    for(Int_t ch=23; ch<=24; ++ch) if(layer==11 && row==42+(ch-23)*2)  flag=1;
    for(Int_t ch=28; ch<=44; ++ch) if(layer==11 && row==52+(ch-28)*2)  flag=1;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==11 && row==86+(ch-46)*2)  flag=1;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==11 && row==106+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==10 && row==ch*2)         flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=30; ++ch) if(layer==10 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==10 && row==60+(ch-32)*2) flag=2;
    for(Int_t ch=49; ch<=54; ++ch) if(layer==10 && row==92+(ch-49)*2) flag=2;
    if(layer==10 && row==86) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==10 && row==1+ch*2)       flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==10 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=27; ch<=29; ++ch) if(layer==10 && row==51+(ch-27)*2) flag=3;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==10 && row==61+(ch-32)*2) flag=3;
    for(Int_t ch=48; ch<=54; ++ch) if(layer==10 && row==91+(ch-48)*2) flag=3;
    if(layer==10 && row==87) flag=3;
    return flag;
  case 11:
    //AGET-0
    for(Int_t ch= 0; ch<= 1; ++ch) if(layer==11 && row==110+ch*2)      flag=0;
    for(Int_t ch= 4; ch<=10; ++ch) if(layer==11 && row==118+(ch-4)*2)  flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==132+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=30; ++ch) if(layer==11 && row==152+(ch-23)*2) flag=0;
    for(Int_t ch=33; ch<=44; ++ch) if(layer==11 && row==172+(ch-33)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==11 && row==196+(ch-46)*2) flag=0;
    if(layer==11 && row==216) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<= 1; ++ch) if(layer==11 && row==109+ch*2)      flag=1;
    for(Int_t ch= 5; ch<=10; ++ch) if(layer==11 && row==119+(ch-5)*2)  flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==131+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=30; ++ch) if(layer==11 && row==151+(ch-23)*2) flag=1;
    for(Int_t ch=34; ch<=44; ++ch) if(layer==11 && row==173+(ch-34)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==11 && row==195+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==11 && row==215+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<= 5; ++ch) if(layer==10 && row==105+ch*2)      flag=2;
    for(Int_t ch= 8; ch<=10; ++ch) if(layer==10 && row==121+(ch-8)*2)  flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==127+(ch-12)*2) flag=2;
    for(Int_t ch=25; ch<=44; ++ch) if(layer==10 && row==151+(ch-25)*2) flag=2;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==10 && row==191+(ch-46)*2) flag=2;
    if(layer==10 && row==147) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<= 6; ++ch) if(layer==10 && row==104+ch*2)      flag=3;
    for(Int_t ch= 8; ch<=10; ++ch) if(layer==10 && row==120+(ch-8)*2)  flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==126+(ch-12)*2) flag=3;
    for(Int_t ch=26; ch<=28; ++ch) if(layer==10 && row==152+(ch-26)*2) flag=3;
    for(Int_t ch=30; ch<=44; ++ch) if(layer==10 && row==160+(ch-30)*2) flag=3;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==10 && row==190+(ch-46)*2) flag=3;
    if(layer==10 && row==146) flag=3;
    return flag;
  case 12:
    //AGET-0
    for(Int_t ch= 0; ch<= 9; ++ch) if(layer==13 && row==1+ch*2)       flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==13 && row==25+(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==13 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==13 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<= 9; ++ch) if(layer==13 && row==ch*2)         flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==13 && row==24+(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==13 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==13 && row==86+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==12 && row==1+ch*2)        flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==23+(ch-12)*2)  flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==43+(ch-23)*2)  flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==87+(ch-46)*2)  flag=2;
    for(Int_t ch=58; ch<=59; ++ch) if(layer==12 && row==109+(ch-58)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==12 && row==ch*2)          flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==22+(ch-12)*2)  flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==42+(ch-23)*2)  flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==86+(ch-46)*2)  flag=3;
    for(Int_t ch=58; ch<=59; ++ch) if(layer==12 && row==108+(ch-58)*2) flag=3;
    return flag;
  case 13:
    //AGET-0
    for(Int_t ch= 1; ch<=10; ++ch) if(layer==13 && row==110+(ch-1)*2)  flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==13 && row==130+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=42; ++ch) if(layer==13 && row==150+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==13 && row==194+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 2; ch<=10; ++ch) if(layer==13 && row==111+(ch-2)*2)  flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==13 && row==129+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=43; ++ch) if(layer==13 && row==149+(ch-23)*2) flag=1;
    for(Int_t ch=47; ch<=55; ++ch) if(layer==13 && row==195+(ch-47)*2) flag=1;
    if(layer==13 && row==213) flag=1;
    //AGET-2
    for(Int_t ch= 4; ch<=10; ++ch) if(layer==12 && row==124+(ch-4)*2)  flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==138+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==158+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==202+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=60; ++ch) if(layer==12 && row==222+(ch-57)*2) flag=2;
    if(layer==12 && row==118) flag=2;
    if(layer==12 && row==120) flag=2;
    //AGET-3
    for(Int_t ch= 5; ch<=10; ++ch) if(layer==12 && row==125+(ch-5)*2)  flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==137+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==157+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==201+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==12 && row==221+(ch-57)*2) flag=3;
    if(layer==12 && row==119) flag=3;
    if(layer==12 && row==121) flag=3;
    return flag;
  case 14:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==15 && row==1+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==15 && row==ch*2)         flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==86+(ch-46)*2) flag=1;
    if(layer==15 && row==106) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==14 && row==ch*2)          flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==22+(ch-12)*2)  flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==42+(ch-23)*2)  flag=2;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==14 && row==86+(ch-46)*2)  flag=2;
    for(Int_t ch=54; ch<=55; ++ch) if(layer==14 && row==102+(ch-54)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==14 && row==1+ch*2)        flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==23+(ch-12)*2)  flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==43+(ch-23)*2)  flag=3;
    for(Int_t ch=46; ch<=50; ++ch) if(layer==14 && row==87+(ch-46)*2)  flag=3;
    for(Int_t ch=54; ch<=55; ++ch) if(layer==14 && row==103+(ch-54)*2) flag=3;
    return flag;
  case 15:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==15 && row==108+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==130+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==150+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==194+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==15 && row==107+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==129+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==149+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==193+(ch-46)*2) flag=1;
    if(layer==15 && row==213) flag=1;
    //AGET-2
    if(layer==14 && row==107) flag=2;
    if(layer==14 && row==109) flag=2;
    for(Int_t ch= 4; ch<=10; ++ch) if(layer==14 && row==115+(ch-4)*2)  flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==129+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==149+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==14 && row==193+(ch-46)*2) flag=2;
    //AGET-3
    if(layer==14 && row==106) flag=3;
    if(layer==14 && row==108) flag=3;
    for(Int_t ch= 5; ch<=10; ++ch) if(layer==14 && row==116+(ch-5)*2)  flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==128+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==148+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==14 && row==192+(ch-46)*2) flag=3;
    return flag;
  case 16:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==17 && row==ch*2)          flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==22+(ch-12)*2)  flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==42+(ch-23)*2)  flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==86+(ch-46)*2)  flag=0;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==106+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==17 && row==1+ch*2)        flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==23+(ch-12)*2)  flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==43+(ch-23)*2)  flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==87+(ch-46)*2)  flag=1;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==107+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==16 && row==ch*2)          flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==22+(ch-12)*2)  flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==42+(ch-23)*2)  flag=2;
    for(Int_t ch=46; ch<=49; ++ch) if(layer==16 && row==86+(ch-46)*2)  flag=2;
    for(Int_t ch=52; ch<=55; ++ch) if(layer==16 && row==98+(ch-52)*2)  flag=2;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==16 && row==106+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==16 && row==1+ch*2)       flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=48; ++ch) if(layer==16 && row==87+(ch-46)*2) flag=3;
    for(Int_t ch=52; ch<=55; ++ch) if(layer==16 && row==99+(ch-52)*2) flag=3;
    if(layer==16 && row== 10) flag=3;
    if(layer==16 && row==109) flag=3;
    return flag;
  case 17:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==17 && row==113+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==135+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==155+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==199+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==219+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==17 && row==112+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==134+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==154+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==198+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==218+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<= 5; ++ch) if(layer==16 && row==111+ch*2)      flag=2;
    for(Int_t ch= 8; ch<=10; ++ch) if(layer==16 && row==127+(ch-8)*2)  flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==133+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==153+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==16 && row==197+(ch-46)*2) flag=2;
    if(layer==16 && row==217) flag=2;
    if(layer==16 && row==219) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<= 5; ++ch) if(layer==16 && row==110+ch*2)      flag=3;
    for(Int_t ch= 9; ch<=10; ++ch) if(layer==16 && row==128+(ch-9)*2)  flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==132+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==152+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==16 && row==196+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==16 && row==216+(ch-57)*2) flag=3;
    return flag;
  case 18:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==19 && row==1+ch*2)        flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==23+(ch-12)*2)  flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==43+(ch-23)*2)  flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==87+(ch-46)*2)  flag=0;
    for(Int_t ch=57; ch<=62; ++ch) if(layer==19 && row==107+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==19 && row==ch*2)          flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==22+(ch-12)*2)  flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==42+(ch-23)*2)  flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==86+(ch-46)*2)  flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==19 && row==106+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==18 && row==ch*2)          flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==22+(ch-12)*2)  flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==42+(ch-23)*2)  flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==86+(ch-46)*2)  flag=2;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==106+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==18 && row==1+ch*2)        flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==23+(ch-12)*2)  flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==43+(ch-23)*2)  flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==87+(ch-46)*2)  flag=3;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==107+(ch-57)*2) flag=3;
    return flag;
  case 19:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==19 && row==120+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==142+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==162+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==206+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=62; ++ch) if(layer==19 && row==226+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==19 && row==119+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==141+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==161+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==205+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==19 && row==225+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==18 && row==117+ch*2)      flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==139+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==159+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==203+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==223+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==18 && row==116+ch*2)      flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==138+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==158+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==202+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==222+(ch-57)*2) flag=3;
    return flag;
  case 20:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==21 && row==ch*2)          flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==22+(ch-12)*2)  flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==42+(ch-23)*2)  flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==86+(ch-46)*2)  flag=0;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==106+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==21 && row==1+ch*2)        flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==23+(ch-12)*2)  flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==43+(ch-23)*2)  flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==87+(ch-46)*2)  flag=1;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==107+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==20 && row==ch*2)          flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==22+(ch-12)*2)  flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==42+(ch-23)*2)  flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==86+(ch-46)*2)  flag=2;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==106+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==20 && row==1+ch*2)        flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==23+(ch-12)*2)  flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==43+(ch-23)*2)  flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==87+(ch-46)*2)  flag=3;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==107+(ch-57)*2) flag=3;
    return flag;
  case 21:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==21 && row==117+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==139+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==159+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==203+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==223+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==21 && row==116+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==138+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==158+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==202+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==222+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==20 && row==123+ch*2)      flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==145+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==165+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==209+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==229+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==20 && row==122+ch*2)      flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==144+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==164+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==208+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==228+(ch-57)*2) flag=3;
    return flag;
  case 22:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==23 && row==1+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==23 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==23 && row==ch*2)         flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==23 && row==86+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==22 && row==1+ch*2)       flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==22 && row==23+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=35; ++ch) if(layer==22 && row==43+(ch-23)*2) flag=2;
    for(Int_t ch=39; ch<=44; ++ch) if(layer==22 && row==75+(ch-39)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==87+(ch-46)*2) flag=2;
    if(layer==22 && row==107) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==22 && row==ch*2)          flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==22 && row==22+(ch-12)*2)  flag=3;
    for(Int_t ch=23; ch<=36; ++ch) if(layer==22 && row==42+(ch-23)*2)  flag=3;
    for(Int_t ch=39; ch<=44; ++ch) if(layer==22 && row==74+(ch-39)*2)  flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==86+(ch-46)*2)  flag=3;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==22 && row==106+(ch-57)*2) flag=3;
    return flag;
  case 23:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==23 && row==106+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==128+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==148+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==23 && row==192+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==23 && row==105+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==127+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==147+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==23 && row==191+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==22 && row==110+ch*2)      flag=2;
    for(Int_t ch=12; ch<=17; ++ch) if(layer==22 && row==132+(ch-12)*2) flag=2;
    if(layer==22 && row==150) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==22 && row==152+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==196+(ch-46)*2) flag=2;
    if(layer==22 && row==216) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==22 && row==109+ch*2)      flag=3;
    for(Int_t ch=12; ch<=18; ++ch) if(layer==22 && row==131+(ch-12)*2) flag=3;
    if(layer==22 && row==149) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==22 && row==151+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==195+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==22 && row==215+(ch-57)*2) flag=3;
    return flag;
  case 24:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==25 && row==1+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==25 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==25 && row==ch*2)         flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==25 && row==86+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==24 && row==1+ch*2)       flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==23+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=28; ++ch) if(layer==24 && row==43+(ch-23)*2) flag=2;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==24 && row==61+(ch-32)*2) flag=2;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==24 && row==87+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==24 && row==ch*2)         flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==22+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==24 && row==42+(ch-23)*2) flag=3;
    for(Int_t ch=33; ch<=44; ++ch) if(layer==24 && row==62+(ch-33)*2) flag=3;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==24 && row==86+(ch-46)*2) flag=3;
    return flag;
  case 25:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==25 && row==102+ch*2)      flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==124+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==144+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==25 && row==188+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==25 && row==101+ch*2)      flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==123+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==143+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==25 && row==187+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==24 && row==104+ch*2)      flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==126+(ch-12)*2) flag=2;
    for(Int_t ch=26; ch<=44; ++ch) if(layer==24 && row==152+(ch-26)*2) flag=2;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==24 && row==190+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==24 && row==103+ch*2)      flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==125+(ch-12)*2) flag=3;
    for(Int_t ch=26; ch<=44; ++ch) if(layer==24 && row==151+(ch-26)*2) flag=3;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==24 && row==189+(ch-46)*2) flag=3;
    return flag;
  case 26:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==27 && row==ch*2)         flag=0;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==27 && row==22+(ch-12)*2) flag=0;
    for(Int_t ch=24; ch<=44; ++ch) if(layer==27 && row==44+(ch-24)*2) flag=0;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==86+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==27 && row==1+ch*2)       flag=1;
    for(Int_t ch=12; ch<=19; ++ch) if(layer==27 && row==23+(ch-12)*2) flag=1;
    if(layer==27 && row==41) flag=1;
    for(Int_t ch=24; ch<=44; ++ch) if(layer==27 && row==45+(ch-24)*2) flag=1;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==87+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==26 && row==ch*2)         flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==86+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==26 && row==1+ch*2)       flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==87+(ch-46)*2) flag=3;
    return flag;
  case 27:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==27 && row==99+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==27 && row==121+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=28; ++ch) if(layer==27 && row==141+(ch-23)*2) flag=0;
    for(Int_t ch=31; ch<=44; ++ch) if(layer==27 && row==157+(ch-31)*2) flag=0;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==185+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==27 && row==98+ch*2)       flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==27 && row==120+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=28; ++ch) if(layer==27 && row==140+(ch-23)*2) flag=1;
    if(layer==27 && row==154) flag=1;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==27 && row==158+(ch-32)*2) flag=1;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==184+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==26 && row==101+ch*2)      flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==123+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==143+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==187+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==26 && row==100+ch*2)      flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==122+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==142+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==186+(ch-46)*2) flag=3;
    return flag;
  case 28:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==29 && row==1+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==29 && row==43+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==29 && row==ch*2)         flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=34; ++ch) if(layer==29 && row==42+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==28 && row==1+ch*2)       flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==23+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==43+(ch-23)*2) flag=2;
    if(layer==28 && row==87) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==28 && row==ch*2)         flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==22+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==42+(ch-23)*2) flag=3;
    if(layer==28 && row==86) flag=3;
    if(layer==28 && row==88) flag=3;
    return flag;
  case 29:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==29 && row==66+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==88+(ch-12)*2)  flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==29 && row==108+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==29 && row==65+ch*2)       flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==87+(ch-12)*2)  flag=1;
    for(Int_t ch=23; ch<=34; ++ch) if(layer==29 && row==107+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==28 && row==90+ch*2)       flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==112+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==132+(ch-23)*2) flag=2;
    if(layer==28 && row==176) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==28 && row==89+ch*2)       flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==111+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==131+(ch-23)*2) flag=3;
    if(layer==28 && row==175) flag=3;
    if(layer==28 && row==177) flag=3;
    return flag;
  case 30:
    //AGET-0
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==31 && row==1+ch*2)       flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==31 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==31 && row==43+(ch-23)*2) flag=0;
    if(layer==31 && row==87) flag=0;
    if(layer==31 && row==89) flag=0;
    //AGET-1
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==31 && row==ch*2)         flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==31 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==31 && row==42+(ch-23)*2) flag=1;
    if(layer==31 && row==86) flag=1;
    if(layer==31 && row==88) flag=1;
    //AGET-2
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==30 && row==ch*2)         flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==30 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==30 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==30 && row==86+(ch-46)*2) flag=2;
    if(layer==30 && row==106) flag=2;
    //AGET-3
    for(Int_t ch= 0; ch<=10; ++ch) if(layer==30 && row==1+ch*2)       flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==30 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==30 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==30 && row==87+(ch-46)*2) flag=3;
    if(layer==30 && row==107) flag=3;
    return flag;
  default:
    return -1;
  }
}

//_____________________________________________________________________________
inline Int_t GetASADId(Int_t layer, Int_t row) //0~30
{
#ifdef PAD_HELPER_DEBUG
  ValidateRow(layer, row, __func__);
#endif
  Int_t flag = layer/4;
  Int_t section;
  if(flag==0) section=0; //layer 0~3
  else if(flag==1) section=1; //layer 4~7
  else if(layer==30||layer==31) section=3; //layer 30~31
  else section=2; //layer 8~29

  Int_t half = padParameter[layer][1]/2;
  Int_t division1 = padParameter[layer][1]/6;
  Int_t division2 = padParameter[layer][1]*5/6;

  switch(section){
  case 0:
    if(row<half)
      return 0;
    else
      return 1;
  case 1:
    Int_t dummy;
    if(layer%4<2) dummy=2;
    if(2<=layer%4) dummy=5;
    if(row<division1||division2<=row)
      return dummy;
    else if(division1<=row&&row<half)
      return dummy+1;
    else
      return dummy+2;
  case 3:
    return 30;
  default:
    if(row<half)
      return layer-layer%2;
    else
      return layer-layer%2+1;
  }
}

//_____________________________________________________________________________
inline Int_t GetCoBoId(Int_t layer, Int_t row)
{
#ifdef PAD_HELPER_DEBUG
  ValidateRow(layer, row, __func__);
#endif
  switch(layer){
  case 4:
    if(60<=row && row<=99)
      return 1;
    else
      return 0;
  case 5:
    if(72<=row && row<=119)
      return 1;
    else
      return 0;
  default:
    return layer/4;
  }
}

//_____________________________________________________________________________
inline Int_t GetPadId(Int_t layer_id, Int_t row_id)
{
  // Return -1 (invalid) for negative row_id to prevent crashes.
  // This is primarily to handle the INT_MIN(-2147483648) value resulting from (Int_t)NaN
  // when initializing TPCHit with NaN in TPCCluster as below;
  //  m_mean_hit(new TPCHit(layer, TMath::QuietNaN())) <- this
  if (row_id < 0) return -1;
#ifdef PAD_HELPER_DEBUG
  ValidateRow(layer_id, row_id, __func__);
#endif

  Int_t pad_id = 0;
  for (Int_t layer = 0; layer < layer_id; layer++)
    pad_id += static_cast<Int_t>(padParameter[layer][kNumOfPad]);
  pad_id += row_id;
  return pad_id;
}

//_____________________________________________________________________________
inline Int_t getLayerID(Int_t pad_id)
{
#ifdef PAD_HELPER_DEBUG
  ValidatePadID(pad_id, __func__);
#endif

  Int_t sum = 0;
  for (Int_t layer = 0; layer < NumOfLayersTPC; layer++) {
    sum += static_cast<Int_t>(padParameter[layer][kNumOfPad]);
    if (pad_id < sum) return layer;
  }
  
  throw Exception(Form("[tpc::%s] Invalid padID %d (Max: %d)", 
                       __func__, pad_id, NumOfPadTPC - 1));
}

//_____________________________________________________________________________
inline Int_t getRowID(Int_t pad_id)
{
#ifdef PAD_HELPER_DEBUG
  ValidatePadID(pad_id, __func__);
#endif

  Int_t sum = 0;
  for (Int_t layer = 0; layer < NumOfLayersTPC; layer++) {
    Int_t n_pad = static_cast<Int_t>(padParameter[layer][kNumOfPad]);
    if (pad_id < sum + n_pad) return pad_id - sum;
    sum += n_pad;
  }

  // Theoretically unreachable. Safeguard for parameter inconsistencies.
  throw Exception(Form("[tpc::%s] Logic Error: padID %d not found in parameter table.", __func__, pad_id));
}

//_____________________________________________________________________________
inline Double_t GetTheta(Int_t pad_id)
{
#ifdef PAD_HELPER_DEBUG
  ValidatePadID(pad_id, __func__);
#endif

  Int_t sum = 0;
  for (Int_t layer = 0; layer < NumOfLayersTPC; layer++) {
    Int_t n_pad = static_cast<Int_t>(padParameter[layer][kNumOfPad]);

    // Check if pad_id belongs to this layer
    if (pad_id < sum + n_pad) {
      Int_t row = pad_id - sum;
      Double_t n_div = padParameter[layer][kNumOfDivision];

      // Calculate theta
      Double_t s_theta = 180. - (360. / n_div) * n_pad / 2.;
      Double_t theta  = s_theta + (row + 0.5) * 360. / n_div - 180;

      return theta;
    }
    sum += n_pad;
  }

  // Theoretically unreachable. Safeguard for parameter inconsistencies.
  throw Exception(Form("[tpc::%s] Logic Error: padID %d not found.", __func__, pad_id));
}

//_____________________________________________________________________________
inline Double_t GetTheta(Int_t layer, Double_t m_row)
{
  Double_t resolved_m_row = TMath::QuietNaN();
  if (!TryResolveMRow(layer, m_row, resolved_m_row)) {
    Double_t max_row = padParameter[layer][kNumOfPad];
    throw Exception(Form("[tpc::%s] Invalid m_row %f for layer %d (Limit: < %f)", 
                         __func__, m_row, layer, max_row));
  }

  Int_t    n_pad   = static_cast<Int_t>(padParameter[layer][kNumOfPad]);
  Double_t n_div   = padParameter[layer][kNumOfDivision];
  Double_t s_theta = 180. - (360. / n_div) * n_pad / 2.;
  Double_t theta  = s_theta + (resolved_m_row + 0.5) * 360. / n_div - 180;

  return theta;
}

//_____________________________________________________________________________
inline Double_t GetMrow(Int_t layer, Double_t m_phi)
{
#ifdef PAD_HELPER_DEBUG
  ValidateLayer(layer, __func__);
#endif
  constexpr Double_t wrap_epsilon = 1.e-4;
  Double_t n_pad = padParameter[layer][kNumOfPad];
  Double_t n_div = padParameter[layer][kNumOfDivision];
  Double_t mrow = 0.5*(n_pad-1.) + (90.-m_phi)*n_div/360.;
  if (mrow < -wrap_epsilon) {
    mrow = 0.5*(n_pad-1.) + (450.-m_phi)*n_div/360.;
  }
  return mrow;
}

//_____________________________________________________________________________
inline Double_t GetRadius(Int_t layer)
{
#ifdef PAD_HELPER_DEBUG
  ValidateLayer(layer, __func__);
#endif
  return padParameter[layer][kRadius];
}

//_____________________________________________________________________________
inline Double_t GetR(Int_t pad_id)
{
#ifdef PAD_HELPER_DEBUG
  ValidatePadID(pad_id, __func__);
#endif

  Int_t sum = 0;
  for (Int_t layer = 0; layer < NumOfLayersTPC; layer++) {
    Int_t n_pad = static_cast<Int_t>(padParameter[layer][kNumOfPad]);

    // Check if pad_id belongs to this layer
    if (pad_id < sum + n_pad) {
      return padParameter[layer][kRadius];
    }
    sum += n_pad;
  }

  // Theoretically unreachable. Safeguard for parameter inconsistencies.
  throw Exception(Form("[tpc::%s] Logic Error: padID %d not found.", __func__, pad_id));
}

//_____________________________________________________________________________
inline TVector3 GetPosition(Int_t pad_id)
{
#ifdef PAD_HELPER_DEBUG
  ValidatePadID(pad_id, __func__);
#endif

  Int_t sum = 0;
  for (Int_t layer = 0; layer < NumOfLayersTPC; layer++) {
    Int_t n_pad = static_cast<Int_t>(padParameter[layer][kNumOfPad]);

    // Check if pad_id belongs to this layer
    if (pad_id < sum + n_pad) {
      Int_t row = pad_id - sum;
      Double_t theta  = GetTheta(layer, static_cast<Double_t>(row))*TMath::DegToRad();
      Double_t radius = padParameter[layer][kRadius];
      Double_t x = radius * std::sin(theta);
      Double_t z = radius * std::cos(theta) + Z_TARGET;
      return TVector3(x, 0., z);
    }
    sum += n_pad;
  }

  // return TVector3(0., -1., 0.);
  return TVector3(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN());
}

//_____________________________________________________________________________
inline TVector3 GetPosition(Int_t layer, Double_t m_row)
{
#ifdef PAD_HELPER_DEBUG
  ValidateLayer(layer, __func__);
#endif

  Double_t resolved_m_row = TMath::QuietNaN();
  if (!TryResolveMRow(layer, m_row, resolved_m_row)) {
    return TVector3(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN());
  }

  Double_t theta  = GetTheta(layer, resolved_m_row) * TMath::DegToRad();
  Double_t radius = padParameter[layer][kRadius];
  Double_t x = radius * std::sin(theta);
  Double_t z = radius * std::cos(theta) + Z_TARGET;
  return TVector3(x, 0., z);
}

//_____________________________________________________________________________
// TPC GEM section (1-4) in the xz plane, split by the z = +/- x diagonals.
// Returns -1 on the sector boundary where |x| == |z|.
inline Int_t GetSection(Double_t x, Double_t z)
{
  Int_t section = -1;
  if (std::abs(x) == std::abs(z)) return section;
  if (std::abs(x) <  std::abs(z)) {
    if (z < 0) section = 1;
    else       section = 3;
  } else {
    if (x < 0) section = 2;
    else       section = 4;
  }
  return section;
}

//_____________________________________________________________________________
inline Int_t GetSection(Int_t pad_id)
{
  TVector3 pad_pos = GetPosition(pad_id);
  return GetSection(pad_pos.x(), pad_pos.z());
}

//_____________________________________________________________________________
inline Int_t GetSection(Int_t layer, Double_t m_row)
{
  const TVector3 pad_pos = GetPosition(layer, m_row);
  return GetSection(pad_pos.x(), pad_pos.z());
}

//_____________________________________________________________________________
// Find PadID from global position (z, x)
// Returns:
//    0 or positive : Valid PadID
//    -layer        : Hit inside the gap between layer and layer-1
//    -1000         : Not found (outside detector volume)
inline Int_t FindPadID(Double_t z, Double_t x)
{
  // 0 <= angle < 360
  Double_t radius = std::hypot(x, z-Z_TARGET);
  Double_t angle  = 180.0 + std::atan2(x, z-Z_TARGET) * TMath::RadToDeg();
  if (angle >= 360.0) angle -= 360.0;
  if (angle <    0.0) angle += 360.0;

  Int_t hit_layer = -1;
  for (Int_t layer = 0; layer < NumOfLayersTPC; layer++) {
    Double_t r_pad = padParameter[layer][kRadius];
    Double_t l_pad = padParameter[layer][kLength];
    Double_t r_in  = r_pad - l_pad * 0.5;
    Double_t r_out = r_pad + l_pad * 0.5;
    if (r_in <= radius && radius <= r_out) {
      hit_layer = layer;
      break;
    }

    // Gap check: If the hit falls into the gap between the previous and current layers
    if (layer > 0) {
      Double_t r_prev_out = padParameter[layer-1][kRadius] + padParameter[layer-1][kLength] * 0.5;
      if (r_prev_out < radius && radius < r_in) {
        return -layer; 
      }
    }
  }
  if (hit_layer == -1) return -1000;

  // Calc row  
  Double_t n_pad = padParameter[hit_layer][kNumOfPad];
  Double_t n_div = padParameter[hit_layer][kNumOfDivision];
  Double_t s_theta = 180. - (360. / n_div) * n_pad / 2.;
  Double_t d_theta = 360. / n_div;

  Double_t diff = angle - s_theta;
  if (std::isnan(diff) || diff < 0) return -1000;
  Int_t row = static_cast<Int_t>(diff / d_theta);
  if (row < 0 || static_cast<Int_t>(n_pad) <= row) return -1000;

  return GetPadId(hit_layer, row);
}

//_____________________________________________________________________________
inline Double_t ArcLength(Int_t layer, Double_t row1, Double_t row2)
{
#ifdef PAD_HELPER_DEBUG
  ValidateLayer(layer, __func__);
#endif

  Double_t radius = padParameter[layer][kRadius];

  // Unit: Degree
  Double_t theta1 = GetTheta(layer, row1);
  Double_t theta2 = GetTheta(layer, row2);
  Double_t diff   = std::abs(theta1 - theta2);

  // diff >= 0, then std::fmod(diff, 360.0) >= 0
  diff = std::fmod(diff, 360.0);

  // There are two paths between two points on a circle: 
  // the minor arc (short) and the major arc (long).
  // We always use the minor arc (shortest distance).
  if (diff > 180.0) diff = 360.0 - diff;

  return radius * diff * TMath::DegToRad();
}

//_____________________________________________________________________________
inline void 
InitializeHistograms(const char* name)
{
  auto h1 = gDirectory->Get<TH2Poly>(name);
  if (!h1) {
    std::cerr << "[tpc::InitializeHistograms] Error: TH2Poly '" << name << "' not found in gDirectory" << std::endl;
    return;
  }

  Double_t X[5];
  Double_t Y[5];
  for (Int_t layer = 0; layer < NumOfLayersTPC; ++layer) {    
    Int_t n_pad      = static_cast<Int_t>(padParameter[layer][kNumOfPad]);
    Double_t n_div   = padParameter[layer][kNumOfDivision];
    Double_t radius = padParameter[layer][kRadius];
    Double_t length = padParameter[layer][kLength];
    Double_t r_min  = radius - length / 2.0;
    Double_t r_max  = radius + length / 2.0;
    Double_t d_theta = 360.0 / n_div;
    Double_t s_theta = - d_theta * n_pad / 2.0;

    for (Int_t j = 0; j < n_pad; ++j) {      
      Double_t theta1 = (s_theta + j*d_theta)     * TMath::DegToRad();
      Double_t theta2 = (s_theta + (j+1)*d_theta) * TMath::DegToRad();

      X[0] = r_max * std::cos(theta1);
      X[1] = r_max * std::cos(theta2);
      X[2] = r_min * std::cos(theta2);
      X[3] = r_min * std::cos(theta1);
      X[4] = X[0];
      for (Int_t k = 0; k < 5; ++k) X[k] += Z_TARGET; 
            
      Y[0] = r_max * std::sin(theta1);
      Y[1] = r_max * std::sin(theta2);
      Y[2] = r_min * std::sin(theta2);
      Y[3] = r_min * std::sin(theta1);
      Y[4] = Y[0];

      h1->AddBin(5, X, Y);
    }
  }
}

//_____________________________________________________________________________
inline Bool_t IsClusterable(Int_t layer, Int_t row_a, Int_t row_b)
{
#ifdef PAD_HELPER_DEBUG
  ValidateRow(layer, row_a, __func__);
  ValidateRow(layer, row_b, __func__);
#endif

  // Ensure min < max
  Int_t min_row = std::min(row_a, row_b);
  Int_t max_row = std::max(row_a, row_b);
  Int_t dist_linear = max_row - min_row;

  // ---------------------------------------------------------
  // Case A: Sector Layers (Layer 10+) - No wrapping
  // ---------------------------------------------------------
  if (layer >= 10) {
    Int_t dead_pads = 0;
    for (Int_t r = min_row + 1; r < max_row; ++r) {
      Int_t pad_id = GetPadId(layer, r);
      if (std::find(std::begin(deadChannel), std::end(deadChannel), pad_id) != std::end(deadChannel)) {
        dead_pads++;
      }
    }
    return (dist_linear <= MAX_ROW_DIF_TPC + dead_pads);
  }

  // ---------------------------------------------------------
  // Case B: Ring Layers (Layer 0-9) - Check Shortest Path
  // ---------------------------------------------------------
  Int_t n_pad = static_cast<Int_t>(padParameter[layer][kNumOfPad]);
  Int_t dist_wrap = n_pad - dist_linear; // Distance across the boundary (0)

  // Determine the shortest path
  // If Linear distance is shorter (or equal), check Linear path
  if (dist_linear <= dist_wrap) {
    Int_t dead_pads = 0;
    // Count dead pads between min_row and max_row
    for (Int_t r = min_row + 1; r < max_row; ++r) {
      Int_t pad_id = GetPadId(layer, r);
      if (std::find(std::begin(deadChannel), std::end(deadChannel), pad_id) != std::end(deadChannel)) {
        dead_pads++;
      }
    }
    return (dist_linear <= MAX_ROW_DIF_TPC + dead_pads);
  }
  // If Wrap-around distance is shorter, check Wrap path
  else {
    Int_t dead_pads = 0;
    // Count dead pads in the wrap-around path:
    // 1. [max_row+1, n_pad-1] (End of array)
    for (Int_t r = max_row + 1; r < n_pad; ++r) {
      Int_t pad_id = GetPadId(layer, r);
      if (std::find(std::begin(deadChannel), std::end(deadChannel), pad_id) != std::end(deadChannel)) {
        dead_pads++;
      }
    }
    // 2. [0, min_row-1] (Start of array)
    for (Int_t r = 0; r < min_row; ++r) {
      Int_t pad_id = GetPadId(layer, r);
      if (std::find(std::begin(deadChannel), std::end(deadChannel), pad_id) != std::end(deadChannel)) {
        dead_pads++;
      }
    }
    return (dist_wrap <= MAX_ROW_DIF_TPC + dead_pads);
  }
}

//_____________________________________________________________________________
// Check if the pad is considered "Dead" (broken channel or on the structural frame)
inline Bool_t IsDead(Int_t pad_id)
{
#ifdef PAD_HELPER_DEBUG
  ValidatePadID(pad_id, __func__);
#endif

  Bool_t on_centerframe = std::find(std::begin(padOnCenterFrame), std::end(padOnCenterFrame), pad_id) 
                          != std::end(padOnCenterFrame);
  if (on_centerframe) return true;

  Bool_t is_dead_channel = std::find(std::begin(deadChannel), std::end(deadChannel), pad_id) 
                           != std::end(deadChannel);
  return is_dead_channel ? true : false;
}

//_____________________________________________________________________________
inline Bool_t IsDead(Int_t layer, Int_t row){
  Int_t pad_id = GetPadId(layer, row);
  return IsDead(pad_id);
}

//_____________________________________________________________________________
inline Bool_t Noise(Int_t pad_id){
#ifdef PAD_HELPER_DEBUG
  ValidatePadID(pad_id, __func__);
#endif

  Bool_t noise = std::find(std::begin(padOnSectionFrame_E72), std::end(padOnSectionFrame_E72), pad_id) != std::end(padOnSectionFrame_E72);
  if(noise) return true;
  else return false;
  
}

//_____________________________________________________________________________
inline Bool_t Noise(Int_t layer, Int_t row){

  Int_t pad_id = GetPadId(layer, row);
  return Noise(pad_id);
}


//Functions for G4 simulation
//_____________________________________________________________________________
//E42
inline Double_t GetClSize1Prob(Double_t mom, Int_t pid, Int_t layer){
  Bool_t inner = false;
  if(layer < 10) inner = true;
  Int_t pid_flag = 0;
  if(abs(pid)<300) pid_flag = 0;
  else pid_flag = 1;
  Int_t mom_flag = (Int_t)(mom/0.1);
  Double_t prob = 0;
  if(mom_flag > 9){
    mom_flag = 9;
    if(inner){
      prob = clusterSizeInner[pid_flag][mom_flag];
    }
    else{
      prob = clusterSizeOuter[pid_flag][mom_flag];
    }
  }
  else{
    Double_t resi = (mom )/0.1 - mom_flag;
    if(inner){
      prob = clusterSizeInner[pid_flag][mom_flag]*(1-resi) + clusterSizeInner[pid_flag][mom_flag+1]*resi; 
    }
    else{
      prob = clusterSizeOuter[pid_flag][mom_flag]*(1-resi) + clusterSizeOuter[pid_flag][mom_flag+1]*resi; 
    }
  }
  return prob; 
}

//_____________________________________________________________________________
//E42
inline Double_t GetDetectionEfficiency(TVector3 pos, Int_t pid, TVector3 mom, Double_t de){
  Int_t pad = FindPadID(pos.z(), pos.x());
  Int_t layer = getLayerID(pad);
  bool inner = false;
  if(layer < 10) inner = true;
  int pid_flag = 0;// 0 -> pion, 1 -> kaon, 2 -> proton or any other baryons.
  if(abs(pid)== 211) pid_flag = 0;
  else if(abs(pid)== 321) pid_flag = 1;
  else if(abs(pid)== 2212 or abs(pid)> 3000) pid_flag = 2;
  double eff = 0.;
  if(inner){
    if(pid_flag == 0 or pid_flag == 1){
      eff = 1;
      if(layer>4){
        eff = 1 - 0.02*(layer-4);
      }
    }
    else if (pid_flag == 2){
      eff = 1;
    }
  }
  else {
    if(pid_flag == 0 or pid_flag == 1){
      eff = 1;
      if(layer > 15 and layer < 29){
        eff = 1 - 0.1*(layer-15)/13;
      }
      else if (layer > 28 and layer < 31){
        eff = 0.9 - 0.05*(layer-28);
      }
      else if(layer > 30){
        eff = 0.8 + 0.1 *(layer - 30);
      }
    }
    else if (pid_flag == 2){
      eff = 1;
    }
  }
  return eff;
}
}
#endif
