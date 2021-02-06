// -*- C++ -*-

#ifndef TPC_PAD_HELPER_HH
#define TPC_PAD_HELPER_HH

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <TMath.h>
#include <TVector3.h>

#include "DetectorID.hh"

namespace tpc
{

enum EPadParameter {
  kLayerID,
  kNumOfPad,
  kNumOfDivision,
  kRadius,
  kLength,
  NPadParameter
};

// #PadID is defined as 0 origin
// layerID, #OfPad, #OfDivision, radius, length
static const Double_t padParameter[NumOfLayersTPC][NPadParameter]=
{ {  0,  48,  48,  14.75,  9.0 },
  {  1,  48,  48,  24.25,  9.0 },
  {  2,  72,  72,  33.75,  9.0 },
  {  3,  96,  96,  43.25,  9.0 },
  {  4, 120, 120,  52.75,  9.0 },
  {  5, 144, 144,  62.25,  9.0 },
  {  6, 168, 168,  71.75,  9.0 },
  {  7, 192, 192,  81.25,  9.0 },
  {  8, 216, 216,  90.75,  9.0 },
  {  9, 240, 240, 100.25,  9.0 },
  { 10, 208, 241, 111.50, 12.5 },
  { 11, 218, 271, 124.50, 12.5 },
  { 12, 230, 300, 137.50, 12.5 },
  { 13, 214, 330, 150.50, 12.5 },
  { 14, 212, 360, 163.50, 12.5 },
  { 15, 214, 390, 176.50, 12.5 },
  { 16, 220, 420, 189.50, 12.5 },
  { 17, 224, 449, 202.50, 12.5 },
  { 18, 232, 479, 215.50, 12.5 },
  { 19, 238, 509, 228.50, 12.5 },
  { 20, 244, 539, 241.50, 12.5 },
  { 21, 232, 569, 254.50, 12.5 },
  { 22, 218, 599, 267.50, 12.5 },
  { 23, 210, 628, 280.50, 12.5 },
  { 24, 206, 658, 293.50, 12.5 },
  { 25, 202, 688, 306.50, 12.5 },
  { 26, 200, 718, 319.50, 12.5 },
  { 27, 196, 748, 332.50, 12.5 },
  { 28, 178, 777, 345.50, 12.5 },
  { 29, 130, 807, 358.50, 12.5 },
  { 30, 108, 837, 371.50, 12.5 },
  { 31,  90, 867, 384.50, 12.5 } };

//_____________________________________________________________________________
inline Double_t getDTheta(Int_t layerID)
{
  return (360./padParameter[layerID][2]);
}

//_____________________________________________________________________________
inline Double_t getsTheta(Int_t layerID)
{
  Double_t sTheta = 180.-(360./padParameter[layerID][2])*padParameter[layerID][1]/2.;
  return sTheta;
}

//_____________________________________________________________________________
inline Int_t
GetPadId( Int_t layer, Int_t row )
{
  Int_t pad = 0;
  for( Int_t l=0; l<layer; ++l ){
    pad += padParameter[l][kNumOfPad];
  }
  return pad + row;
}

//_____________________________________________________________________________
inline Int_t getLayerID(Int_t padID)
{
  int layer;
  int sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  return layer;
}

//_____________________________________________________________________________
inline Int_t getRowID(Int_t padID)
{
  int layer, row;
  int sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  row = padID - sum;
  return row;
}

//_____________________________________________________________________________
inline Double_t getR(Int_t padID)
{
  return padParameter[getLayerID(padID)][3];
}

//_____________________________________________________________________________
inline int findPadID(double z, double x)
{
  z += 143;
  double radius = sqrt(x*x + z*z);
  double angle;
  if (z == 0)
  {
    if (x > 0)   angle = 1.5*TMath::Pi();
    else if (x < 0)   angle = 0.5*TMath::Pi();
    else return -1000; // no padID at (0,0)
  }
  else if (z < 0)
    angle = atan(x / z);
  else  // z > 0
    angle = TMath::Pi() + atan(x / z);

  int layer, row;
  // find layerID
  for (layer = 0; !(padParameter[layer][3]+padParameter[layer][4]*0.5 >= radius
                    && padParameter[layer][3]-padParameter[layer][4]*0.5 <= radius); layer++)
  {
    if (layer >= 32) return -1000;
    if (layer != 0)
    {
      if (padParameter[layer][3] - padParameter[layer][4] * 0.5 >= radius &&
          padParameter[layer - 1][3] + padParameter[layer - 1][4] * 0.5 <= radius) return -layer;
    }
  }

  // find rowID
  if (angle - (getsTheta(layer)*TMath::Pi()/180.) < 0) return -2000;

  row = (int)((angle-(getsTheta(layer)*TMath::Pi()/180.))
              /(getDTheta(layer)*TMath::Pi()/180.));
  if (row > padParameter[layer][1]) return -1000;

  return GetPadId( layer, row );
}

//_____________________________________________________________________________
inline Double_t getTheta(Int_t padID)
{
  int layer = getLayerID(padID);
  int row   = getRowID(padID);

  Double_t sTheta = getsTheta(layer);
  Double_t theta = sTheta+(row+0.5)*360./padParameter[layer][3] - 180.;
  return theta; // -180 ~ 180
}

//_____________________________________________________________________________
inline TVector3 getPosition(Int_t padID)
{
  int layer, row;
  int sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  row = padID - sum;

  TVector3 pos;
  if (row >= padParameter[layer][1]){ // out of range
    pos.SetX(0);
    pos.SetY(0);
    pos.SetZ(0);
  }
  else{
    double x, z;
    x = padParameter[layer][3] * sin(getTheta(padID)*TMath::Pi()/180.);
    z = padParameter[layer][3] * cos(getTheta(padID)*TMath::Pi()/180.) - 143.;
    pos.SetX(x);
    pos.SetY(0);
    pos.SetZ(z);
  }
  return pos;
}

}

#endif
