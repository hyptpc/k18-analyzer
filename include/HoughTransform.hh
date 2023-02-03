// -*- C++ -*-

#ifndef HOUGH_TRANSFORM_HH
#define HOUGH_TRANSFORM_HH

#include <vector>

#include <TString.h>
#include <TVector3.h>
#include "DCAnalyzer.hh"

class TPCCluster;

namespace tpc
{

//Hough-Transform functions for track searching
Bool_t HoughTransformLineXZ(std::vector<TVector3> gHitPos,
			    Int_t *MaxBin, Double_t *LinearPar,
			    Int_t MinNumOfHits=8);
Bool_t HoughTransformLineXZ(const std::vector<TPCClusterContainer>& ClCont,
			    Int_t *MaxBin, Double_t *LinearPar,
			    Int_t MinNumOfHits=8);
Bool_t HoughTransformCircleXZ(std::vector<TVector3> gHitPos,
			      Int_t *MaxBin, Double_t *HelixPar,
			      Int_t MinNumOfHits=8);
Bool_t HoughTransformCircleXZ(const std::vector<TPCClusterContainer>& ClCont,
			      Int_t *MaxBin, Double_t *HelixPar,
			      Int_t MinNumOfHits=8);

//Hough-Transform functions for track parameter calculation
void HoughTransformLineYPhi(std::vector<TVector3> gHitPos,
			    Double_t *HelixPar, Double_t MaxHoughWindow);
void HoughTransformLineYPhi(const std::vector<TPCClusterContainer>& ClCont,
			    Double_t *HelixPar, Double_t MaxHoughWindow);

}

#endif
