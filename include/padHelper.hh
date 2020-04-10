#ifndef PADHELPER_HH
#define PADHELPER_HH
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <TMath.h>
#include <TVector3.h>

namespace padHelper
{
  //#PadID is defined as 0 origin
  //#OfPad #division #radius padLength
  static const Double_t padParameter[32][6]=
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


  inline Double_t getDTheta(Int_t layerID){
    return (360./padParameter[layerID][3]);
  }

  inline Double_t getsTheta(Int_t layerID)
  {
    Double_t sTheta = 180.-(360./padParameter[layerID][3])*padParameter[layerID][1]/2.;
    return sTheta;
  }

  inline Double_t getRadius(Int_t layerID)
  {
    return padParameter[layerID][2];
  }


  inline Int_t getPadID(Int_t layerID, Int_t rowID)
  {
    Int_t padID=0;
    for(int layi = 0 ; layi<layerID; layi++) padID += padParameter[layi][1];
    padID+=rowID;
    return padID;

  }

  inline Int_t getLayerID(Int_t padID)
  {
    //    padID-=1;
    int layer;
    int sum = 0;

    for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
      {
	sum += padParameter[layer][1];
      }
    return layer;
  }

  inline Int_t getRowID(Int_t padID)
  {
    //    padID-=1;
    int layer, row;
    int sum = 0;

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


  inline int findPadID(double z, double x) 
  { 
    z += 143;
    double radius = sqrt(x*x + z*z);
    double angle;
    if (z == 0)
      {
	if (x > 0)   angle = 1.5*TMath::Pi();
	else if (x < 0)   angle = 0.5*TMath::Pi();
	else return -1000; // no padID if (0,0)
      }
    else
      {
	//angle = tan(-x / z);
	angle = TMath::Pi() - atan(-x / z);
      }

    int layer, row;
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
    
    
    //std::cout<<"padHelper:: layer="<<layer<<", angle="<<angle<<", "<<(getsTheta(layer)*TMath::Pi()/180.)<<std::endl;
    // find row_num
    //  if (angle - (padParameter[layer][4]*TMath::Pi()/180.) < 0) return -1000;
    if (angle - (getsTheta(layer)*TMath::Pi()/180.) < 0) return -2000;

    //    double a, b, c;
    //row = (int)((angle-(padParameter[layer][4]*TMath::Pi()/180.))/(padParameter[layer][3]*TMath::Pi()/180.));
    //    row = (int)((angle-(getsTheta(layer)*TMath::Pi()/180.))/(getDTheta(layer)*TMath::Pi()/180.));
    
    //row = (int)((angle-(getsTheta(layer)*TMath::Pi()/180.))/(getDTheta(layer)*TMath::Pi()/180.))+1;
    row = (int)((angle-(getsTheta(layer)*TMath::Pi()/180.))/(getDTheta(layer)*TMath::Pi()/180.));
    if (row > padParameter[layer][1]) return -1000;

    return getPadID(layer, row);
  }




  inline Double_t getTheta(Int_t padID)
  {
    //    padID-=1;
    int layer, row;
    int sum = 0;

    for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
      {
	sum += padParameter[layer][1];
      }
    row = padID - sum;
    Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
    Double_t theta = sTheta+(row+0.5)*(360.-2*sTheta)/padParameter[layer][1];
    return theta;
  }

  inline Double_t getR(Int_t padID)
  {
    //    padID-=1;
    int layer, row;
    int sum = 0;

    for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
      {
	sum += padParameter[layer][1];
      }
    row = padID - sum;
    Double_t R = padParameter[layer][2];
    return R;
  }

  inline TVector3 getPoint(int padID)
  {
    //padID-=1;
    int layer, row;
    int sum = 0;

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
      double x, z;
      Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
      x = padParameter[layer][2] * -sin((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.);
      z = padParameter[layer][2] * -cos((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.) - 143.0;
      result.SetX(x);
      result.SetY(0);
      result.SetZ(z);
    }
    return result;
  }

}
#endif
