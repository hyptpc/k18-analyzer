// -*- C++ -*-

#ifndef TPC_PAD_HELPER_HH
#define TPC_PAD_HELPER_HH

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>

#include <TDirectory.h>
#include <TH2Poly.h>
#include <TMath.h>
#include <TVector3.h>

#include "DetectorID.hh"

namespace tpc
{
const Double_t ZTarget = -143.; // Target from center

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
inline Int_t GetAGETId(Int_t asad, Int_t layer, Int_t row)
{
  Int_t flag=-1;
  switch(asad){
  case 0:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==1 && row==22-ch*2) flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==46-(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==28-(ch-23)*2) flag=0;
    if(layer==1 && row==0) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==1 && row==23-ch*2) flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==47-(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==29-(ch-23)*2) flag=1;
    if(layer==1 && row==1) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==0 && row==22-ch*2) flag=2;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==34-(ch-13)*2) flag=2;
    for(Int_t ch=23; ch<=31; ++ch) if(layer==2 && row==16-(ch-23)*2) flag=2;
    if(layer==0 && row==0) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==0 && row==23-ch*2) flag=3;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==35-(ch-13)*2) flag=3;
    for(Int_t ch=23; ch<=31; ++ch) if(layer==2 && row==17-(ch-23)*2) flag=3;
    if(layer==0 && row==1) flag=3;
    return flag;
  case 1:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==1 && row==47-ch*2) flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==95-(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==77-(ch-23)*2) flag=0;
    if(layer==1 && row==25) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==1 && row==46-ch*2) flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==94-(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==76-(ch-23)*2) flag=1;
    if(layer==1 && row==24) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==0 && row==47-ch*2) flag=2;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==71-(ch-13)*2) flag=2;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==2 && row==53-(ch-23)*2) flag=2;
    if(layer==0 && row==25) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==0 && row==46-ch*2) flag=3;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==70-(ch-13)*2) flag=3;
    for(Int_t ch=23; ch<=31; ++ch) if(layer==2 && row==52-(ch-23)*2) flag=3;
    if(layer==0 && row==24) flag=3;
    return flag;
  case 2:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==5 && row==22-ch*2) flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==5 && row==143-(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==125-(ch-23)*2) flag=0;
    if(layer==5 && row==0) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==5 && row==23-ch*2) flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==5 && row==142-(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==124-(ch-23)*2) flag=1;
    if(layer==5 && row==1) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=9; ++ch) if(layer==4 && row==18-ch*2) flag=2;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==117-(ch-12)*2) flag=2;
    if(layer==4 && row==119) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=9; ++ch) if(layer==4 && row==19-ch*2) flag=3;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==116-(ch-12)*2) flag=3;
    if(layer==4 && row==118) flag=3;
    return flag;
  case 3:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==5 && row==24+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==46+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==66+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==5 && row==25+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==47+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==67+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==4 && row==20+ch*2) flag=2;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==42+(ch-12)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==4 && row==21+ch*2) flag=3;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==43+(ch-12)*2) flag=3;
    return flag;
  case 4:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==5 && row==73+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==95+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==115+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==5 && row==72+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==94+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==114+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==4 && row==61+ch*2) flag=2;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==83+(ch-12)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==4 && row==60+ch*2) flag=3;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==82+(ch-12)*2) flag=3;
    return flag;
  case 5:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==7 && row==30-ch*2) flag=0;
    for(Int_t ch=12; ch<=16; ++ch) if(layer==7 && row==8-(ch-12)*2) flag=0;
    for(Int_t ch=17; ch<=21; ++ch) if(layer==7 && row==191-(ch-17)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==181-(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==7 && row==31-ch*2) flag=1;
    for(Int_t ch=12; ch<=16; ++ch) if(layer==7 && row==9-(ch-12)*2) flag=1;
    for(Int_t ch=17; ch<=21; ++ch) if(layer==7 && row==190-(ch-17)*2) flag=1;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==180-(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==6 && row==26-ch*2) flag=2;
    for(Int_t ch=12; ch<=14; ++ch) if(layer==6 && row==4-(ch-12)*2) flag=2;
    for(Int_t ch=15; ch<=21; ++ch) if(layer==6 && row==167-(ch-15)*2) flag=2;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==153-(ch-23)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==6 && row==27-ch*2) flag=3;
    for(Int_t ch=12; ch<=14; ++ch) if(layer==6 && row==5-(ch-12)*2) flag=3;
    for(Int_t ch=15; ch<=21; ++ch) if(layer==6 && row==166-(ch-15)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==152-(ch-23)*2) flag=3;
    return flag;
  case 6:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==7 && row==32+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==54+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==74+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==7 && row==33+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==55+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==75+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==6 && row==28+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==50+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==70+(ch-23)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==6 && row==29+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==51+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==71+(ch-23)*2) flag=3;
    return flag;
  case 7:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==7 && row==97+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==119+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==139+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==7 && row==96+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==118+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==138+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==6 && row==85+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==107+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==127+(ch-23)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==6 && row==84+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==106+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==126+(ch-23)*2) flag=3;
    return flag;
  case 8:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==9 && row==ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==22+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==42+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==86+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==106+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==9 && row==1+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==23+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==43+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==87+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==107+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==8 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==86+(ch-46)*2) flag=2;
    if(layer==8 && row==106) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==8 && row==1+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==87+(ch-46)*2) flag=3;
    if(layer==8 && row==107) flag=3;
    return flag;
  case 9:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==9 && row==121+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==143+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==163+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==207+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==227+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==9 && row==120+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==142+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==162+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==206+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==226+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==8 && row==109+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==131+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==151+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==195+(ch-46)*2) flag=2;
    if(layer==8 && row==215) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==8 && row==108+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==130+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==150+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==194+(ch-46)*2) flag=3;
    if(layer==8 && row==214) flag=3;
    return flag;
  case 10:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==11 && row==1+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==11 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==11 && row==87+(ch-46)*2) flag=0;
    if(layer==11 && row==107) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==11 && row==ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=24; ++ch) if(layer==11 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=28; ch<=44; ++ch) if(layer==11 && row==52+(ch-28)*2) flag=1;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==11 && row==86+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==11 && row==106+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==10 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=30; ++ch) if(layer==10 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==10 && row==60+(ch-32)*2) flag=2;
    for(Int_t ch=49; ch<=54; ++ch) if(layer==10 && row==92+(ch-49)*2) flag=2;
    if(layer==10 && row==86) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==10 && row==1+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==10 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=27; ch<=29; ++ch) if(layer==10 && row==51+(ch-27)*2) flag=3;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==10 && row==61+(ch-32)*2) flag=3;
    for(Int_t ch=48; ch<=54; ++ch) if(layer==10 && row==91+(ch-48)*2) flag=3;
    if(layer==10 && row==87) flag=3;
    return flag;
  case 11:
    //AGET-0
    for(Int_t ch=0; ch<=1; ++ch) if(layer==11 && row==110+ch*2) flag=0;
    for(Int_t ch=4; ch<=10; ++ch) if(layer==11 && row==118+(ch-4)*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==132+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=30; ++ch) if(layer==11 && row==152+(ch-23)*2) flag=0;
    for(Int_t ch=33; ch<=44; ++ch) if(layer==11 && row==172+(ch-33)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==11 && row==196+(ch-46)*2) flag=0;
    if(layer==11 && row==216) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=1; ++ch) if(layer==11 && row==109+ch*2) flag=1;
    for(Int_t ch=5; ch<=10; ++ch) if(layer==11 && row==119+(ch-5)*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==131+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=30; ++ch) if(layer==11 && row==151+(ch-23)*2) flag=1;
    for(Int_t ch=34; ch<=44; ++ch) if(layer==11 && row==173+(ch-34)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==11 && row==195+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==11 && row==215+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=5; ++ch) if(layer==10 && row==105+ch*2) flag=2;
    for(Int_t ch=8; ch<=10; ++ch) if(layer==10 && row==121+(ch-8)*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==127+(ch-12)*2) flag=2;
    for(Int_t ch=25; ch<=44; ++ch) if(layer==10 && row==151+(ch-25)*2) flag=2;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==10 && row==191+(ch-46)*2) flag=2;
    if(layer==10 && row==147) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=6; ++ch) if(layer==10 && row==104+ch*2) flag=3;
    for(Int_t ch=8; ch<=10; ++ch) if(layer==10 && row==120+(ch-8)*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==126+(ch-12)*2) flag=3;
    for(Int_t ch=26; ch<=28; ++ch) if(layer==10 && row==152+(ch-26)*2) flag=3;
    for(Int_t ch=30; ch<=44; ++ch) if(layer==10 && row==160+(ch-30)*2) flag=3;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==10 && row==190+(ch-46)*2) flag=3;
    if(layer==10 && row==146) flag=3;
    return flag;
  case 12:
    //AGET-0
    for(Int_t ch=0; ch<=9; ++ch) if(layer==13 && row==1+ch*2) flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==13 && row==25+(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==13 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==13 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=9; ++ch) if(layer==13 && row==ch*2) flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==13 && row==24+(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==13 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==13 && row==86+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==12 && row==1+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==23+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==43+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==87+(ch-46)*2) flag=2;
    for(Int_t ch=58; ch<=59; ++ch) if(layer==12 && row==109+(ch-58)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==12 && row==ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==22+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==42+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==86+(ch-46)*2) flag=3;
    for(Int_t ch=58; ch<=59; ++ch) if(layer==12 && row==108+(ch-58)*2) flag=3;
    return flag;
  case 13:
    //AGET-0
    for(Int_t ch=1; ch<=10; ++ch) if(layer==13 && row==110+(ch-1)*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==13 && row==130+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=42; ++ch) if(layer==13 && row==150+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==13 && row==194+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=2; ch<=10; ++ch) if(layer==13 && row==111+(ch-2)*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==13 && row==129+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=43; ++ch) if(layer==13 && row==149+(ch-23)*2) flag=1;
    for(Int_t ch=47; ch<=55; ++ch) if(layer==13 && row==195+(ch-47)*2) flag=1;
    if(layer==13 && row==213) flag=1;
    //AGET-2
    for(Int_t ch=4; ch<=10; ++ch) if(layer==12 && row==124+(ch-4)*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==138+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==158+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==202+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=60; ++ch) if(layer==12 && row==222+(ch-57)*2) flag=2;
    if(layer==12 && row==118) flag=2;
    if(layer==12 && row==120) flag=2;
    //AGET-3
    for(Int_t ch=5; ch<=10; ++ch) if(layer==12 && row==125+(ch-5)*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==137+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==157+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==201+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==12 && row==221+(ch-57)*2) flag=3;
    if(layer==12 && row==119) flag=3;
    if(layer==12 && row==121) flag=3;
    return flag;
  case 14:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==15 && row==1+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==15 && row==ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==86+(ch-46)*2) flag=1;
    if(layer==15 && row==106) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==14 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==14 && row==86+(ch-46)*2) flag=2;
    for(Int_t ch=54; ch<=55; ++ch) if(layer==14 && row==102+(ch-54)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==14 && row==1+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=50; ++ch) if(layer==14 && row==87+(ch-46)*2) flag=3;
    for(Int_t ch=54; ch<=55; ++ch) if(layer==14 && row==103+(ch-54)*2) flag=3;
    return flag;
  case 15:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==15 && row==108+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==130+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==150+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==194+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==15 && row==107+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==129+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==149+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==193+(ch-46)*2) flag=1;
    if(layer==15 && row==213) flag=1;
    //AGET-2
    if(layer==14 && row==107) flag=2;
    if(layer==14 && row==109) flag=2;
    for(Int_t ch=4; ch<=10; ++ch) if(layer==14 && row==115+(ch-4)*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==129+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==149+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==14 && row==193+(ch-46)*2) flag=2;
    //AGET-3
    if(layer==14 && row==106) flag=3;
    if(layer==14 && row==108) flag=3;
    for(Int_t ch=5; ch<=10; ++ch) if(layer==14 && row==116+(ch-5)*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==128+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==148+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==14 && row==192+(ch-46)*2) flag=3;
    return flag;
  case 16:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==17 && row==ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==22+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==42+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==86+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==106+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==17 && row==1+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==23+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==43+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==87+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==107+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==16 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=49; ++ch) if(layer==16 && row==86+(ch-46)*2) flag=2;
    for(Int_t ch=52; ch<=55; ++ch) if(layer==16 && row==98+(ch-52)*2) flag=2;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==16 && row==106+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==16 && row==1+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=48; ++ch) if(layer==16 && row==87+(ch-46)*2) flag=3;
    for(Int_t ch=52; ch<=55; ++ch) if(layer==16 && row==99+(ch-52)*2) flag=3;
    if(layer==16 && row==10) flag=3;
    if(layer==16 && row==109) flag=3;
    return flag;
  case 17:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==17 && row==113+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==135+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==155+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==199+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==219+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==17 && row==112+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==134+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==154+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==198+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==218+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=5; ++ch) if(layer==16 && row==111+ch*2) flag=2;
    for(Int_t ch=8; ch<=10; ++ch) if(layer==16 && row==127+(ch-8)*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==133+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==153+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==16 && row==197+(ch-46)*2) flag=2;
    if(layer==16 && row==217) flag=2;
    if(layer==16 && row==219) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=5; ++ch) if(layer==16 && row==110+ch*2) flag=3;
    for(Int_t ch=9; ch<=10; ++ch) if(layer==16 && row==128+(ch-9)*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==132+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==152+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==16 && row==196+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==16 && row==216+(ch-57)*2) flag=3;
    return flag;
  case 18:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==19 && row==1+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==87+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=62; ++ch) if(layer==19 && row==107+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==19 && row==ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==86+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==19 && row==106+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==18 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==86+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==106+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==18 && row==1+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==87+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==107+(ch-57)*2) flag=3;
    return flag;
  case 19:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==19 && row==120+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==142+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==162+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==206+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=62; ++ch) if(layer==19 && row==226+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==19 && row==119+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==141+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==161+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==205+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==19 && row==225+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==18 && row==117+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==139+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==159+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==203+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==223+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==18 && row==116+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==138+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==158+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==202+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==222+(ch-57)*2) flag=3;
    return flag;
  case 20:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==21 && row==ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==22+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==42+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==86+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==106+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==21 && row==1+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==23+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==43+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==87+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==107+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==20 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==86+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==106+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==20 && row==1+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==87+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==107+(ch-57)*2) flag=3;
    return flag;
  case 21:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==21 && row==117+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==139+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==159+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==203+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==223+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==21 && row==116+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==138+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==158+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==202+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==222+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==20 && row==123+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==145+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==165+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==209+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==229+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==20 && row==122+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==144+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==164+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==208+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==228+(ch-57)*2) flag=3;
    return flag;
  case 22:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==23 && row==1+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==23 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==23 && row==ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==23 && row==86+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==22 && row==1+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==22 && row==23+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=35; ++ch) if(layer==22 && row==43+(ch-23)*2) flag=2;
    for(Int_t ch=39; ch<=44; ++ch) if(layer==22 && row==75+(ch-39)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==87+(ch-46)*2) flag=2;
    if(layer==22 && row==107) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==22 && row==ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==22 && row==22+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=36; ++ch) if(layer==22 && row==42+(ch-23)*2) flag=3;
    for(Int_t ch=39; ch<=44; ++ch) if(layer==22 && row==74+(ch-39)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==86+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==22 && row==106+(ch-57)*2) flag=3;
    return flag;
  case 23:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==23 && row==106+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==128+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==148+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==23 && row==192+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==23 && row==105+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==127+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==147+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==23 && row==191+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==22 && row==110+ch*2) flag=2;
    for(Int_t ch=12; ch<=17; ++ch) if(layer==22 && row==132+(ch-12)*2) flag=2;
    if(layer==22 && row==150) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==22 && row==152+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==196+(ch-46)*2) flag=2;
    if(layer==22 && row==216) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==22 && row==109+ch*2) flag=3;
    for(Int_t ch=12; ch<=18; ++ch) if(layer==22 && row==131+(ch-12)*2) flag=3;
    if(layer==22 && row==149) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==22 && row==151+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==195+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==22 && row==215+(ch-57)*2) flag=3;
    return flag;
  case 24:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==25 && row==1+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==25 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==25 && row==ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==25 && row==86+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==24 && row==1+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==23+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=28; ++ch) if(layer==24 && row==43+(ch-23)*2) flag=2;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==24 && row==61+(ch-32)*2) flag=2;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==24 && row==87+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==24 && row==ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==22+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==24 && row==42+(ch-23)*2) flag=3;
    for(Int_t ch=33; ch<=44; ++ch) if(layer==24 && row==62+(ch-33)*2) flag=3;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==24 && row==86+(ch-46)*2) flag=3;
    return flag;
  case 25:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==25 && row==102+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==124+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==144+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==25 && row==188+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==25 && row==101+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==123+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==143+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==25 && row==187+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==24 && row==104+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==126+(ch-12)*2) flag=2;
    for(Int_t ch=26; ch<=44; ++ch) if(layer==24 && row==152+(ch-26)*2) flag=2;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==24 && row==190+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==24 && row==103+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==125+(ch-12)*2) flag=3;
    for(Int_t ch=26; ch<=44; ++ch) if(layer==24 && row==151+(ch-26)*2) flag=3;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==24 && row==189+(ch-46)*2) flag=3;
    return flag;
  case 26:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==27 && row==ch*2) flag=0;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==27 && row==22+(ch-12)*2) flag=0;
    for(Int_t ch=24; ch<=44; ++ch) if(layer==27 && row==44+(ch-24)*2) flag=0;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==86+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==27 && row==1+ch*2) flag=1;
    for(Int_t ch=12; ch<=19; ++ch) if(layer==27 && row==23+(ch-12)*2) flag=1;
    if(layer==27 && row==41) flag=1;
    for(Int_t ch=24; ch<=44; ++ch) if(layer==27 && row==45+(ch-24)*2) flag=1;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==87+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==26 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==86+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==26 && row==1+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==87+(ch-46)*2) flag=3;
    return flag;
  case 27:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==27 && row==99+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==27 && row==121+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=28; ++ch) if(layer==27 && row==141+(ch-23)*2) flag=0;
    for(Int_t ch=31; ch<=44; ++ch) if(layer==27 && row==157+(ch-31)*2) flag=0;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==185+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==27 && row==98+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==27 && row==120+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=28; ++ch) if(layer==27 && row==140+(ch-23)*2) flag=1;
    if(layer==27 && row==154) flag=1;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==27 && row==158+(ch-32)*2) flag=1;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==184+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==26 && row==101+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==123+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==143+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==187+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==26 && row==100+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==122+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==142+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==186+(ch-46)*2) flag=3;
    return flag;
  case 28:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==29 && row==1+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==29 && row==43+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==29 && row==ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=34; ++ch) if(layer==29 && row==42+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==28 && row==1+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==23+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==43+(ch-23)*2) flag=2;
    if(layer==28 && row==87) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==28 && row==ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==22+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==42+(ch-23)*2) flag=3;
    if(layer==28 && row==86) flag=3;
    if(layer==28 && row==88) flag=3;
    return flag;
  case 29:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==29 && row==66+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==88+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==29 && row==108+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==29 && row==65+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==87+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=34; ++ch) if(layer==29 && row==107+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==28 && row==90+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==112+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==132+(ch-23)*2) flag=2;
    if(layer==28 && row==176) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==28 && row==89+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==111+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==131+(ch-23)*2) flag=3;
    if(layer==28 && row==175) flag=3;
    if(layer==28 && row==177) flag=3;
    return flag;
  case 30:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==31 && row==1+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==31 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==31 && row==43+(ch-23)*2) flag=0;
    if(layer==31 && row==87) flag=0;
    if(layer==31 && row==89) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==31 && row==ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==31 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==31 && row==42+(ch-23)*2) flag=1;
    if(layer==31 && row==86) flag=1;
    if(layer==31 && row==88) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==30 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==30 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==30 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==30 && row==86+(ch-46)*2) flag=2;
    if(layer==30 && row==106) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==30 && row==1+ch*2) flag=3;
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
  Int_t flag=layer/4;
  Int_t section;
  if(flag==0) section=0; //layer 0~3
  else if(flag==1) section=1; //layer 4~7
  else if(layer==30||layer==31) section=3; //layer 30~31
  else section=2; //layer 8~29

  Int_t half=padParameter[layer][1]/2;
  Int_t division1=padParameter[layer][1]/6;
  Int_t division2=padParameter[layer][1]*5/6;

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
inline Int_t GetPadId(Int_t layerID, Int_t rowID)
{
  //Original
  //Int_t padID=0;
  // Check!!!!!!!
  Int_t padID=1;
  for(Int_t layi = 0 ; layi<layerID; layi++) padID += padParameter[layi][1];
  padID+=rowID;
  return padID;

}

inline Int_t getLayerID(Int_t padID)
{
  padID-=1;
  Int_t layer;
  Int_t sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  return layer;
}

inline Int_t getRowID(Int_t padID)
{
  padID-=1;
  Int_t layer, row;
  Int_t sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  row = padID - sum;
  return row;
}
/*
  Double_t getTheta(Int_t layerID, Int_t rowID)
  {
  Double_t sTheta = 180.-(360./padParameter[layerID][3])*padParameter[layerID][1]/2.;
  Double_t theta = sTheta+(rowID+0.5)*(360.-2*sTheta)/padParameter[layerID][1];
  return theta;
  }
*/

inline Double_t getTheta(Int_t padID)
{
  padID-=1;
  Int_t layer, row;
  Int_t sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  row = padID - sum;
  //std::cout<<"layer="<<layer<<", row="<<row<<std::endl;
  Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
  //Double_t theta = sTheta+(row+0.5)*(360.-2*sTheta)/padParameter[layer][1];
  Double_t theta = sTheta+(row+0.5)*360./padParameter[layer][3]-180;
  //std::cout<<"theta="<<theta<<std::endl;

  return theta;
}

inline Double_t getTheta(Int_t layer, Double_t m_row)
{
  Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
  //Double_t theta = sTheta+(row+0.5)*(360.-2*sTheta)/padParameter[layer][1];
  Double_t theta = sTheta+(m_row+0.5)*360./padParameter[layer][3]-180;
  //std::cout<<"theta="<<theta<<std::endl;
  return theta;
}

inline Double_t GetRadius(Int_t layer)
{
  return padParameter[layer][2];
}

inline Double_t getR(Int_t padID)
{
  padID-=1;
  Int_t layer;
  Int_t sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  Double_t R = padParameter[layer][2];
  return R;
}

inline TVector3 getPosition(Int_t padID)
{
  padID-=1;
  Int_t layer, row;
  Int_t sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  row = padID - sum;

  TVector3 result;
  if (row > padParameter[layer][1]){ // out of range
    result.SetX(0);
    result.SetY(-1);
    result.SetZ(0);
  }
  else{
    Double_t x, z;
    //Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;

    //    x = padParameter[layer][2] * -sin((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.);
    //    z = padParameter[layer][2] * -cos((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.) + ZTarget;

    // x = padParameter[layer][2] * sin(getTheta(padID+1)*TMath::Pi()/180.);
    // z = padParameter[layer][2] * cos(getTheta(padID+1)*TMath::Pi()/180.) + ZTarget;
    //std::cout<<"layer="<<layer<<", row"<<row<<std::endl;
    x = padParameter[layer][2] * sin(getTheta(layer,row)*TMath::Pi()/180.);
    z = padParameter[layer][2] * cos(getTheta(layer,row)*TMath::Pi()/180.) + ZTarget;

    // Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
    // Double_t x_ = padParameter[layer][2] * -sin((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.);
    // Double_t z_ = padParameter[layer][2] * -cos((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.) + ZTarget;
    //std::cout<<"x="<<x<<", z="<<z<<", x_="<<x_<<", z_="<<z_<<std::endl;
    result.SetX(x);
    result.SetY(0);
    result.SetZ(z);
  }
  return result;
}

inline TVector3 getPosition(Int_t layer, Double_t m_row)
{
  TVector3 result;
  if(m_row > padParameter[layer][1]){ // out of range
    return TVector3(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN());
  }
  else{
    return TVector3(
      padParameter[layer][2] * sin(getTheta(layer, m_row)*TMath::Pi()/180.),
      0.,
      padParameter[layer][2] * cos(getTheta(layer, m_row)*TMath::Pi()/180.) + ZTarget);
  }
}

inline Int_t findPadID(Double_t z, Double_t x)
{
  z -= ZTarget;
  Double_t radius = sqrt(x*x + z*z);
  Double_t angle;
  if (z == 0)
  {
    if (x > 0)   angle = 1.5*TMath::Pi();
    else if (x < 0)   angle = 0.5*TMath::Pi();
    else return -1000; // no padID if (0,0)
  }
  //  else
  //  {
  //    if (z < 0 && x < 0) angle = atan(x / z);
  //    else if (z > 0 && x < 0) angle = TMath::Pi() - atan(-x / z);
  //    else if (z > 0 && x > 0) angle = TMath::Pi() + atan(x / z);
  //    else if (z < 0 && x > 0) angle = 2*TMath::Pi() - atan(-x / z);
  //  }
//  else if (z < 0) angle = atan(x / z);
//  else  angle = TMath::Pi() + atan(x / z);
  else{
		if (z > 0) angle = TMath::Pi()+atan(x / z);
	  else if( z < 0&&x<0) angle = atan(x / z);
		  else angle = 2*TMath::Pi()+ atan(x / z);//angle of z<0&&x>0 plane should be [1.5Pi,2Pi], not [-0.5Pi , 0].
	}
	//cout << " angle: " << angle*180/TMath::Pi() << endl;

  Int_t layer, row;
  // find layer_num.
  for (layer = 0; !(padParameter[layer][2]+padParameter[layer][5]*0.5 >= radius
                    && padParameter[layer][2]-padParameter[layer][5]*0.5 <= radius); layer++)
  {
    if (layer >= 32) return -1000;
    if (layer != 0)
    {
      if (padParameter[layer][2] - padParameter[layer][5] * 0.5 >= radius &&
          padParameter[layer - 1][2] + padParameter[layer - 1][5] * 0.5 <= radius) return -layer;
    }
  }

  //cout << " layer: " << layer << endl;

  Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;

  // find row_num
  if (angle - (sTheta*TMath::Pi()/180.) < 0) return -1000;

  //Double_t a, b, c;
  row = (int)((angle-(sTheta*TMath::Pi()/180.))/(360./padParameter[layer][3]*TMath::Pi()/180.));
  if (row > padParameter[layer][1]) return -1000;

  //cout << " row: " << row << endl;
  //This is original one
  //return GetPadId(layer, row)+1;
  //Please check
  return GetPadId(layer, row);
}

//_____________________________________________________________________________
inline Double_t
ArcLength(Int_t layer, Double_t row1, Double_t row2)
{
  const Int_t R = padParameter[layer][2];
  Double_t theta = getTheta(layer, row1) - getTheta(layer, row2);
  theta = std::fmod(theta, 2*TMath::Pi());
  if(theta < 0) theta += 2*TMath::Pi();
  theta = TMath::Min(theta, 2*TMath::Pi() - theta);
  return R*theta;
}

//_____________________________________________________________________________
inline void
InitializeHistograms()
{
  std::vector<TH2Poly*> target;
  TList* list = gDirectory->GetList();
  TIter itr(list);
  while(itr.Next()){
    const TString& name((*itr)->GetName());
    const TString& cname((*itr)->ClassName());
    std::cout << " " << std::setw(8) << std::left << name
              << "(" << cname << ")" << std::endl;
    if(cname.EqualTo("TH2Poly")){
      target.push_back(dynamic_cast<TH2Poly*>(*itr));
    }
  }

  Double_t X[5];
  Double_t Y[5];

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    Double_t pLength = tpc::padParameter[layer][5];
    Double_t st = 180.-(360./tpc::padParameter[layer][3])
      * tpc::padParameter[layer][1]/2.;
    Double_t sTheta  = (-1+st/180.)*TMath::Pi();
    Double_t dTheta  = (360./tpc::padParameter[layer][3])/180.*TMath::Pi();
    Double_t cRad    = tpc::padParameter[layer][2];
    Int_t    nPad    = tpc::padParameter[layer][1];
    for(Int_t j=0; j<nPad; ++j){
      X[1] = (cRad+(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
      X[2] = (cRad+(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
      X[3] = (cRad-(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
      X[4] = (cRad-(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
      X[0] = X[4];
      Y[1] = (cRad+(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
      Y[2] = (cRad+(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
      Y[3] = (cRad-(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
      Y[4] = (cRad-(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
      Y[0] = Y[4];
      for(Int_t ii=0; ii<5; ++ii) X[ii] += ZTarget;
      for(auto& h: target) h->AddBin(5, X, Y);
    }
  }
}

//_____________________________________________________________________________
inline Bool_t
IsClusterable(Int_t layer, Int_t row_a, Int_t row_b, Int_t maxdif=1)
{
  if(layer < 10){
    const Int_t npad = padParameter[layer][kNumOfPad];
    return (TMath::Abs(row_a - row_b) <= maxdif
            || TMath::Abs(row_a - row_b) >= npad - maxdif);
  }else{
    return (TMath::Abs(row_a - row_b) <= maxdif);
  }
}
}

#endif
