// -*- C++ -*-

/*
Track searching process
1. Hough-transform on the XZ plane
2. Hough-transform on the vertical plane
    -> Get initial track parameters.
3. Check residual(Houghdist) in the Hough-space and
make a inital track with clusters with small residual.

Track fitting process
1. First fitting with inital track.
2. If fitting is succedeed, check the residual of other cluster with the track.
If residual is acceptable, add cluster into the track.
3. If the track is extended by added clusters, then do fitting again.

Detailed fitting procedures are explained in the TPCLocalTrack/Helix.
*/


#include "TPCTrackSearch.hh"

#include <chrono>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <TH2D.h>
#include <TH3D.h>

#include "DebugTimer.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"
#include "ConfMan.hh"
#include "HoughTransform.hh"
#include "TPCPadHelper.hh"
#include "TPCLTrackHit.hh"
#include "TPCLocalTrack.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCCluster.hh"
#include "TPCBeamRemover.hh"
#include "RootHelper.hh"

#define DebugDisp          0
#define RemainHitTest      0
#define TwoTrackTargetVtx  1 //Distinguish two tracks originating from the target

namespace
{
  const auto qnan = TMath::QuietNaN();
  const double& HS_field_0 = 0.9860;
  const double& valueHSCalc = ConfMan::Get<Double_t>("HSFLDCALC");
  const double& valueHSHall = ConfMan::Get<Double_t>("HSFLDHALL");
  const double Const = 0.299792458;
  const auto& gUser = UserParamMan::GetInstance();

  //const Int_t    MaxNumOfTrackTPC = 5;
  const Int_t    MaxNumOfTrackTPC = 20;
  const Double_t KuramaXZWindow = 30;
  const Double_t KuramaYWindow = 10;
  const Double_t K18XZWindow = 10.5;
  const Double_t K18YWindow = 10;

  const Double_t TargetVtxWindow = 25.;

  // Houghflags
  const Int_t GoodForTracking = 100;
  const Int_t K18Tracks = 200;
  const Int_t KuramaTracks = 300;
  const Int_t AccidentalBeams = 400;
  const Int_t BadForTracking = 500;
  const Int_t Candidate = 1000;

  // Tracks in the Hough-Space
  std::vector<Double_t> XZhough_x;
  std::vector<Double_t> XZhough_y;
  std::vector<Double_t> XZhough_z;
  std::vector<Double_t> Yhough_x;
  std::vector<Double_t> Yhough_y;

  TH1D* hist_y = new TH1D("hist_y", "hist_y", 140, -350, 350);

  //_____________________________________________________________________________
  // Local Functions

  //_____________________________________________________________________________
  template <typename T> void
  CalcTracks(std::vector<T*>& TrackCont)
  {

    for(auto& track: TrackCont) track->Calculate();
  }

  //_____________________________________________________________________________
  template <typename T> void
  ExclusiveTracking(std::vector<T*>& TrackCont)
  {
    for(auto& track: TrackCont){
      track->DoFitExclusive();
      track->CalculateExclusive();
    }
  }

  //_____________________________________________________________________________
  template <typename T> void
  ClearFlags(std::vector<T*>& TrackCont)
  {
    for(const auto& track: TrackCont){
      if(!track) continue;
      Int_t nh = track->GetNHit();
      for(Int_t j=0; j<nh; ++j) track->GetHit(j)->QuitTrack();
    }
  }

  //_____________________________________________________________________________
  //reset houghflag of remain clusters
  void
  ResetHoughFlag(const std::vector<TPCClusterContainer>& ClCont)
  {

    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()==Candidate) hit->SetHoughFlag(0);
      } //ci
    } //layer
  }


  //_____________________________________________________________________________
  template <typename T> void
  GetTrackHitCont(T* track, std::vector<TVector3>& gHitPos)
  {
    Int_t n = track->GetNHit();
    for(Int_t i=0; i<n; ++i){
      TPCLTrackHit *hitp = track->GetHit(i);
      TVector3 pos = hitp->GetLocalHitPos();
      gHitPos.push_back(pos);
    }
  }

  //_____________________________________________________________________________
  Double_t GetMagneticField()
  {
    return HS_field_0*(valueHSHall/valueHSCalc);
  }
}

//_____________________________________________________________________________
namespace tpc
{

//_____________________________________________________________________________
template <typename T> void
FitTrack(T* Track, Int_t Houghflag,
	 const std::vector<TPCClusterContainer>& ClCont,
	 std::vector<T*>& TrackCont,
	 std::vector<T*>& TrackContFailed,
	 Int_t MinNumOfHits)
{

  std::chrono::milliseconds sec;
  auto fit_start = std::chrono::high_resolution_clock::now();

#if DebugDisp
  Track->Print(FUNC_NAME+" Initial track");
  //Track->Print(FUNC_NAME+" Initial track", true);
#endif

  //first fitting
  if(Track->DoFit()){ // No MinNumOfHits cut for the inital Hough track
    auto first_fit = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(first_fit - fit_start);
    Track->SetHoughFlag(Candidate);
    Track->SetFitTime(sec.count());
    Track->SetFitFlag(1);

#if DebugDisp
    Track->Print(FUNC_NAME+" First fitting is succeeded");
    //Track->Print(FUNC_NAME+" First fitting is succeeded", true);
#endif

    T *ExtendedTrack = new T(Track);
    //Residual check with other hits
    bool Add_rescheck = false;
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()>0) continue;
	TVector3 pos = cl->GetPosition();
	TVector3 res = TVector3(cl->ResolutionX(),
				cl->ResolutionY(),
				cl->ResolutionZ());
	double resi=0.;
	if(ExtendedTrack->ResidualCheck(pos, res, resi)){
	  Int_t vtxflag = Track -> GetVtxFlag();
	  Int_t side = Track -> Side(pos);
	  //Vertex inside the target : vtxflag = -1 or 1 / outside vtxflag = 0
	  //if track and new cluster are on the same side and vertex in the target : vtxflag*side = 1
	  //if vertex is outside of the target : vtxflag*side = 0
	  if(vtxflag*side >= 0){
	    ExtendedTrack->AddTPCHit(new TPCLTrackHit(hit));
	    Add_rescheck = true;
	  }
	}
      } //ci
    } //layer

    if(!Add_rescheck){ //No adding and fitting is succeeded
      delete ExtendedTrack;
      if(Track->GetNHit() >= MinNumOfHits){ //track w/ enough clusters
	Track->SetHoughFlag(Houghflag);
	TrackCont.push_back(Track);
      }
      else{
	Track->SetHoughFlag(BadForTracking); //track w/ few clusters
	Track->SetFlag(0);
	TrackContFailed.push_back(Track);
      }
#if DebugDisp
      Track->Print(FUNC_NAME+" No added cluster");
      //Track->Print(FUNC_NAME+" No added cluster", true);
#endif
    }
    else{ //More clusters are added into the track
      auto hitadd = std::chrono::high_resolution_clock::now();
#if DebugDisp
      ExtendedTrack->Print(FUNC_NAME+" After cluster adding");
      //ExtendedTrack->Print(FUNC_NAME+" After cluster adding", true);
#endif
      if(ExtendedTrack->DoFit(MinNumOfHits)){ //2nd fitting success
	delete Track;
	auto second_fit = std::chrono::high_resolution_clock::now();
	sec = std::chrono::duration_cast<std::chrono::milliseconds>(second_fit - hitadd);
	ExtendedTrack->SetFitTime(sec.count());
	ExtendedTrack->SetFitFlag(2);
	ExtendedTrack->SetHoughFlag(Houghflag);
	TrackCont.push_back(ExtendedTrack);

#if DebugDisp
	ExtendedTrack->Print(FUNC_NAME+" Fitting success after adding hits (residual check)");
	//ExtendedTrack->Print(FUNC_NAME+" Fitting success after adding hits (residual check)", true);
#endif
      }
      else{ //2nd fitting failure
	delete ExtendedTrack;
	//push 1st fitting track
	auto hitadd_fail = std::chrono::high_resolution_clock::now();
	sec = std::chrono::duration_cast<std::chrono::milliseconds>(hitadd_fail - hitadd);
	Track->SetFitTime(sec.count());
	Track->SetFitFlag(3);
	if(Track->GetNHit() >= MinNumOfHits){ //track w/ enough clusters
	  Track->SetHoughFlag(Houghflag);
	  TrackCont.push_back(Track);
	}
	else{
	  Track->SetHoughFlag(BadForTracking); //track w/ few clusters
	  TrackContFailed.push_back(Track);
	}

#if DebugDisp
	Track->Print(FUNC_NAME+" 2nd fitting failure!");
	//Track->Print(FUNC_NAME+" 2nd fitting failure!", true);
#endif
      }
    }
  }
  else{ //First fitting is failed. Reset clusters' Hough flag
    auto after_1stfit_fail = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_1stfit_fail - fit_start);
    Track->SetFitTime(sec.count());
    Track->SetFitFlag(4);
    Track->SetHoughFlag(0);
    Track->SetFlag(0);
    Track->FinalizeTrack();
    TrackContFailed.push_back(Track);

#if DebugDisp
    Track->Print(FUNC_NAME+" First fitting is failed");
    //Track->Print(FUNC_NAME+" First fitting is failed", true);
#endif
  }

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

}

//_____________________________________________________________________________
Int_t
LocalTrackSearch(const std::vector<TPCClusterContainer>& ClCont,
		 std::vector<TPCLocalTrack*>& TrackCont,
		 std::vector<TPCLocalTrack*>& TrackContFailed,
		 bool Exclusive,
		 Int_t MinNumOfHits /*=8*/)
{

  const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

  XZhough_x.clear();
  XZhough_y.clear();
  Yhough_x.clear();
  Yhough_y.clear();

  bool prev_add = true;
  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    if(!prev_add) continue;
    prev_add = false;

#if DebugDisp
    std::cout<<FUNC_NAME+" tracki : "<<tracki<<std::endl;
#endif

    std::chrono::milliseconds sec;
    auto before_hough = std::chrono::high_resolution_clock::now();

    //Line Hough-transform on the XZ plane
    Double_t LinearPar[4]; Int_t MaxBinXZ[2];
    if(!tpc::HoughTransformLineXZ(ClCont, MaxBinXZ, LinearPar, MinNumOfHits)){
#if DebugDisp
      std::cout<<FUNC_NAME+" No more track candiate! tracki : "<<tracki<<std::endl;
#endif
      break;
    }

    //Line Hough-transform on the YZ or YX plane
    Int_t MaxBinY[3];
    if(TMath::Abs(LinearPar[2]) < 1) tpc::HoughTransformLineYZ(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);
    else tpc::HoughTransformLineYX(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrack *track = new TPCLocalTrack;
    track->SetHoughParam(LinearPar);
    track->SetParamUsingHoughParam();
    track->CalcClosestDist();

    //If two tracks are merged at the target, seperate them and recalculate params.
    Bool_t vtx_flag = false;
    prev_add = MakeLinearTrack(track, vtx_flag, ClCont, LinearPar, MaxHoughWindowY);
    if(!prev_add) continue;
    if(vtx_flag){
      std::vector<TVector3> gHitPos;
      GetTrackHitCont(track, gHitPos);
      tpc::HoughTransformLineXZ(gHitPos, MaxBinXZ, LinearPar, MinNumOfHits);
      if(TMath::Abs(LinearPar[2]) < 1) tpc::HoughTransformLineYZ(gHitPos, MaxBinY, LinearPar, MaxHoughWindowY);
      else tpc::HoughTransformLineYX(gHitPos, MaxBinY, LinearPar, MaxHoughWindowY);
      track->SetHoughParam(LinearPar);
      track->SetParamUsingHoughParam();
      track->CalcClosestDist();
      track->SetVtxFlag(track->Side(gHitPos[0]));
    }
    if(track -> GetNDF()<1){
      track->SetHoughFlag(Candidate);
      TrackContFailed.push_back(track);
      continue;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<XZhough_x.size(); ++i){
      Int_t bindiffXZ = TMath::Abs(MaxBinXZ[0] - XZhough_x[i]) + TMath::Abs(MaxBinXZ[1] - XZhough_y[i]);
      Int_t bindiffY = TMath::Abs(MaxBinY[0] - Yhough_x[i]) + TMath::Abs(MaxBinY[1] - Yhough_y[i]);
      if(bindiffXZ<=1 && bindiffY<=1){
	hough_flag = false;
#if DebugDisp
	std::cout<<"Previous hough bin on the XZ plane "<<i<<"th x: "
		 <<XZhough_x[i]<<", y: "<<XZhough_y[i]<<" on the vertical plane x: "
		 <<Yhough_x[i]<<", y: "<<Yhough_y[i]<<std::endl;
	std::cout<<"Current hough bin on the XZ plane "<<i<<"th x: "
		 <<MaxBinXZ[0]<<", y: "<<MaxBinXZ[1]<<" on the vertical plane x: "
		 <<MaxBinY[0]<<", y: "<<MaxBinY[1]<<std::endl;
#endif
      }
    }
    XZhough_x.push_back(MaxBinXZ[0]);
    XZhough_y.push_back(MaxBinXZ[1]);
    Yhough_x.push_back(MaxBinY[0]);
    Yhough_y.push_back(MaxBinY[1]);

    if(!hough_flag){
#if DebugDisp
      std::cout<<FUNC_NAME+" The same track is found by Hough-Transform : tracki : "<<tracki<<" hough cont size : "<<XZhough_x.size()<<std::endl;
#endif
      track->SetHoughFlag(Candidate);
      track->SetFitFlag(0);
      TrackContFailed.push_back(track);
      continue;
    }

    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_hough - before_hough);
    track->SetSearchTime(sec.count());

    //Track fitting processes
    FitTrack(track, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  }// tracki
  ResetHoughFlag(ClCont);

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

  CalcTracks(TrackCont);
  CalcTracks(TrackContFailed);
  if(Exclusive) ExclusiveTracking(TrackCont);

  return TrackCont.size();
}

//_____________________________________________________________________________
Bool_t
MakeLinearTrack(TPCLocalTrack *Track, Bool_t &VtxFlag,
		const std::vector<TPCClusterContainer>& ClCont,
		Double_t *LinearPar, Double_t MaxHoughWindowY){

  Bool_t status = false;

  VtxFlag = false; //Vtx in the target, need to check whether two tracks are meged or not
#if TwoTrackTargetVtx
  Double_t closeDist = Track -> GetClosestDist();
  if(closeDist < TargetVtxWindow) VtxFlag = true;
#endif

  //Check Hough-distance and add hits
  int id = 0;
  std::vector<Int_t> side1_hits; std::vector<Int_t> side2_hits;
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(hit->GetHoughFlag()>0) continue;
      TVector3 pos = cl->GetPosition();
      Double_t distXZ = TMath::Abs(LinearPar[2]*(pos.Z() - tpc::ZTarget) - pos.X() +
				   LinearPar[0])/TMath::Sqrt(TMath::Sq(LinearPar[2])+1.);
      Double_t distYZ = TMath::Abs(LinearPar[3]*(pos.Z() - tpc::ZTarget) - pos.Y() +
				   LinearPar[1])/TMath::Sqrt(TMath::Sq(LinearPar[3])+1.);
      if(distXZ < MaxHoughWindowY && distYZ < MaxHoughWindowY){
	hit->SetHoughDist(distXZ);
	hit->SetHoughDistY(distYZ);
	Track->AddTPCHit(new TPCLTrackHit(hit));
	if(VtxFlag){
	  if(Track -> Side(pos) > 0.) side1_hits.push_back(id);
	  else side2_hits.push_back(id);
	}
	id++;
	status = true;
      }
    } //ci
  } //layer

  if(VtxFlag){
    Track -> SetVtxFlag(1);
    if(side1_hits.size() >= side2_hits.size()){
      Track -> EraseHits(side2_hits);
      Track -> SetVtxFlag(1);
    }
    else{
      Track -> EraseHits(side1_hits);
      Track -> SetVtxFlag(-1);
    }
  }

#if DebugDisp
  Track->Print(FUNC_NAME+" Inital track after track searching");
  std::cout<<FUNC_NAME+" Remain clusters : "<<Track -> GetNHit()<<"/"<<id<<std::endl;
#endif

  if(!status) delete Track;
  return status;
}

//_____________________________________________________________________________
Bool_t
MakeHelixTrack(TPCLocalTrackHelix *Track, Bool_t &VtxFlag,
	       const std::vector<TPCClusterContainer>& ClCont,
	       Double_t *HelixPar, Double_t MaxHoughWindow,
	       Double_t MaxHoughWindowY)
{

  Bool_t status = false;

  VtxFlag = false; //Vtx in the target, need to check whether two tracks are meged or not
#if TwoTrackTargetVtx
  Double_t closeDist = Track -> GetClosestDist();
  if(closeDist < TargetVtxWindow) VtxFlag = true;
#endif

  //Check Hough-distance and add hits
  int id = 0;
  std::vector<Int_t> side1_hits; std::vector<Int_t> side2_hits;
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(hit->GetHoughFlag()>0) continue;
      TVector3 pos = cl->GetPosition();
      double tmpx = -pos.x();
      double tmpy = pos.z() - tpc::ZTarget;
      double tmpz = pos.y();
      Double_t r_cal = TMath::Sqrt(pow(tmpx - HelixPar[0], 2) + pow(tmpy - HelixPar[1], 2));
      Double_t dist = TMath::Abs(r_cal - HelixPar[3]);
      if(dist < MaxHoughWindow){
	double tmpt = TMath::ATan2(tmpy - HelixPar[1], tmpx - HelixPar[0]);
	double tmp_xval = HelixPar[3]*tmpt;
	double distY = TMath::Abs(HelixPar[4]*tmp_xval - tmpz + HelixPar[2])/TMath::Sqrt(pow(HelixPar[4], 2) + 1.);
	if(distY < MaxHoughWindowY){
	  hit->SetHoughDist(dist);
	  hit->SetHoughDistY(distY);
	  Track->AddTPCHit(new TPCLTrackHit(hit));
	  if(VtxFlag){
	    if(Track -> Side(pos)==1) side1_hits.push_back(id);
	    else side2_hits.push_back(id);
	  }
	  id++;
	  status = true;
	} //distY
      } //dist
    } //ci
  } //layer

  if(VtxFlag){
    Track -> SetVtxFlag(1);
    if(side1_hits.size() >= side2_hits.size()){
      Track -> EraseHits(side2_hits);
      Track -> SetVtxFlag(1);
    }
    else{
      Track -> EraseHits(side1_hits);
      Track -> SetVtxFlag(-1);
    }
  }

#if DebugDisp
  Track->Print(FUNC_NAME+" Inital track after track searching");
  std::cout<<FUNC_NAME+" Remain clusters : "<<Track -> GetNHit()<<"/"<<id<<std::endl;
#endif

  if(!status) delete Track;
  return status;
}

//_____________________________________________________________________________
void
TestRemainingHits(const std::vector<TPCClusterContainer>& ClCont,
		  std::vector<TPCLocalTrackHelix*>& TrackCont,
		  Int_t MinNumOfHits)
{

  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()==0||hit->GetHoughFlag()==BadForTracking)
	  track->AddTPCHit(new TPCLTrackHit(hit));
      }
    }
    if(track->DoFit(MinNumOfHits)){
      track->SetHoughFlag(GoodForTracking + TrackCont.size());
      TrackCont.push_back(track);
#if DebugDisp
      track->Print(FUNC_NAME+" test with remaining hits");
#endif
    }
    else{
      delete track;
      break;
    }
  }
}

//_____________________________________________________________________________
void
HelixTrackSearch(Int_t Trackflag, Int_t Houghflag,
		 const std::vector<TPCClusterContainer>& ClCont,
		 std::vector<TPCLocalTrackHelix*>& TrackCont,
		 std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		 Int_t MinNumOfHits /*=8*/)
{
  // HoughTransform binning
  const auto MaxHoughWindow = gUser.GetParameter("MaxHoughWindow");
  const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

  bool prev_add = true;
  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    if(!prev_add) continue;
    prev_add = false;

#if DebugDisp
    std::cout<<FUNC_NAME+" tracki : "<<tracki<<std::endl;
#endif

    std::chrono::milliseconds sec;
    auto before_hough = std::chrono::high_resolution_clock::now();

    //Circle Hough-transform
    Double_t HelixPar[5]; Int_t MaxBinXZ[3];
    if(!tpc::HoughTransformCircleXZ(ClCont, MaxBinXZ, HelixPar, MinNumOfHits)){
#if DebugDisp
      std::cout<<FUNC_NAME+" No more circle candiate! tracki : "<<tracki<<std::endl;
#endif
      break;
    }

    //Linear Hough-transform
    Int_t MaxBinY[3];
    tpc::HoughTransformLineYTheta(ClCont, MaxBinY, HelixPar, MaxHoughWindow);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    track->SetHoughParam(HelixPar);
    track->SetParamUsingHoughParam();
    track->CalcClosestDist();
    track->SetFlag(Trackflag);

    //If two tracks are merged at the target, seperate them and recalculate params.
    Bool_t vtx_flag = false;
    prev_add = MakeHelixTrack(track, vtx_flag, ClCont, HelixPar, MaxHoughWindow, MaxHoughWindowY);
    if(!prev_add) continue;
    if(vtx_flag){
      std::vector<TVector3> gHitPos;
      GetTrackHitCont(track, gHitPos);
      tpc::HoughTransformCircleXZ(gHitPos, MaxBinXZ, HelixPar, MinNumOfHits);
      tpc::HoughTransformLineYTheta(gHitPos, MaxBinY, HelixPar, MaxHoughWindow);
      track->SetHoughParam(HelixPar);
      track->SetParamUsingHoughParam();
      track->CalcClosestDist();
      track->SetVtxFlag(track->Side(gHitPos[0]));
    }
    if(track -> GetNDF()<1){
      track->SetHoughFlag(Candidate);
      track->FinalizeTrack();
      TrackContFailed.push_back(track);
      continue;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<XZhough_x.size(); ++i){
      Int_t bindiffXZ = TMath::Abs(MaxBinXZ[0] - XZhough_x[i]) + TMath::Abs(MaxBinXZ[1] - XZhough_y[i]) + TMath::Abs(MaxBinXZ[2] - XZhough_z[i]);
      Int_t bindiffY = TMath::Abs(MaxBinY[0] - Yhough_x[i]) + TMath::Abs(MaxBinY[1] - Yhough_y[i]);
      if(bindiffXZ<=1 && bindiffY<=1){
	hough_flag = false;
#if DebugDisp
	std::cout<<"Previous hough bin on the XZ plane "<<i<<"th x: "
		 <<XZhough_x[i]<<", y: "<<XZhough_y[i]<<", z: "<<XZhough_z[i]<<" on the vertical plane x: "
		 <<Yhough_x[i]<<", y: "<<Yhough_y[i]<<std::endl;
	std::cout<<"Current hough bin on the XZ plane "<<i<<"th x: "
		 <<MaxBinXZ[0]<<", y: "<<MaxBinXZ[1]<<", z: "<<MaxBinXZ[2]<<" on the vertical plane x: "
		 <<MaxBinY[0]<<", y: "<<MaxBinY[1]<<std::endl;
#endif
      }
    }
    XZhough_x.push_back(MaxBinXZ[0]);
    XZhough_y.push_back(MaxBinXZ[1]);
    XZhough_z.push_back(MaxBinXZ[2]);
    Yhough_x.push_back(MaxBinY[0]);
    Yhough_y.push_back(MaxBinY[1]);

    if(!hough_flag){
#if DebugDisp
      std::cout<<FUNC_NAME+" The same track is found by Hough-Transform : tracki : "<<tracki<<" hough cont size : "<<XZhough_x.size()<<std::endl;
#endif
      track->SetHoughFlag(Candidate);
      track->SetFitFlag(0);
      track->FinalizeTrack();
      TrackContFailed.push_back(track);
      continue;
    }

    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_hough - before_hough);
    track->SetSearchTime(sec.count());

    //Track fitting processes
    FitTrack(track, Houghflag, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  }//tracki
  ResetHoughFlag(ClCont);

}

//_____________________________________________________________________________
void
KuramaTrackSearch(std::vector<std::vector<TVector3>> VPs,
		  const std::vector<TPCClusterContainer>& ClCont,
		  std::vector<TPCLocalTrackHelix*>& TrackCont,
		  std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		  Int_t MinNumOfHits /*=8*/)
{

  if(VPs.size()==0) return;
  for(Int_t nt=0; nt<VPs.size(); nt++){
#if DebugDisp
    std::cout<<FUNC_NAME+" Kurama Track "<<nt<<"/"<<VPs.size()<<std::endl;
#endif
    TPCLocalTrackHelix *trackref = new TPCLocalTrackHelix();
    for(Int_t i=0; i<VPs[nt].size(); i++){
      TVector3 pos = VPs[nt][i];
      trackref->AddVPHit(pos);
#if DebugDisp
      std::cout<<FUNC_NAME+" Kurama VP "<<i<<"th pos : "<<VPs[nt][i]<<std::endl;
#endif
    } //i
    trackref->DoVPFit();
    trackref->SetIsKurama();
    TrackContFailed.push_back(trackref);

#if DebugDisp
    std::cout<<FUNC_NAME+" Kurama RK cx : "<<trackref->Getcx()<<" cy : "<<trackref->Getcy()<<" r : "<<trackref->Getr()<<std::endl;
    std::cout<<FUNC_NAME+" Kurama RK p : "<<trackref->Getr()*Const*GetMagneticField()<<std::endl;
#endif

    std::vector<TPCClusterContainer> ClContKurama(NumOfLayersTPC);
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()>0) continue;
	TVector3 pos = cl->GetPosition();
	if(trackref->ResidualCheck(pos, KuramaXZWindow, KuramaYWindow)){
	  ClContKurama[layer].push_back(cl);
	}
      } //ci
    } //layer

    Int_t Trackflag = 1*0 + 2*0 + 4*1 + 8*0; // isBeam, isK18, isKurama, isAccidental
    HelixTrackSearch(Trackflag, KuramaTracks, ClContKurama, TrackCont, TrackContFailed, MinNumOfHits);
    ResetHoughFlag(ClContKurama);
  } //nt

}

//_____________________________________________________________________________
void
K18TrackSearch(std::vector<std::vector<TVector3>> VPs,
	       const std::vector<TPCClusterContainer>& ClCont,
	       std::vector<TPCLocalTrackHelix*>& TrackCont,
	       std::vector<TPCLocalTrackHelix*>& TrackContFailed,
	       Int_t MinNumOfHits /*=8*/)
{

  if(VPs.size()==0) return;
  for(Int_t nt=0; nt<VPs.size(); nt++){
#if DebugDisp
    std::cout<<FUNC_NAME+" K18 Track "<<nt<<"/"<<VPs.size()<<std::endl;
#endif
    TPCLocalTrackHelix *trackref = new TPCLocalTrackHelix();
    for(Int_t i=0; i<VPs[nt].size(); i++){
      TVector3 pos = VPs[nt][i];
      trackref->AddVPHit(pos);
#if DebugDisp
      std::cout<<FUNC_NAME+" K18 RK "<<i<<"th pos : "<<VPs[nt][i]<<std::endl;
#endif
    } //i
    trackref->DoVPFit();
    trackref->SetIsBeam();
    TrackContFailed.push_back(trackref);

#if DebugDisp
    std::cout<<FUNC_NAME+" K18 RK cx : "<<trackref->Getcx()<<" cy : "<<trackref->Getcy()<<" z0 : "<<trackref->Getz0()<<" r : "<<trackref->Getr()<<" dz : "<<trackref->Getdz()<<std::endl;
    std::cout<<FUNC_NAME+" K18 RK p : "<<trackref->Getr()*Const*GetMagneticField()<<std::endl;
#endif

    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();

    Int_t BeforeTGTHits = 0; Int_t Hits = 0;
    std::vector<TPCClusterContainer> ClContK18(NumOfLayersTPC);
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      Double_t minresi = 1000.; Int_t id = -1;
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()>0) continue;
	TVector3 pos = cl->GetPosition();
	if(pos.Z()>tpc::ZTarget) continue;
	Double_t resi;
	if(trackref->ResidualCheck(pos, K18XZWindow, K18YWindow, resi)){
	  ClContK18[layer].push_back(cl);
	  if(minresi > resi){
	    minresi = resi;
	    id = ci;
	  }
	  Hits++;
	}
      } //ci
      if(minresi<100.){
	TPCHit* hit = ClCont[layer][id]->GetMeanHit();
	hit->SetHoughDist(qnan);
	track->AddTPCHit(new TPCLTrackHit(hit));
	BeforeTGTHits++;
      }
    } //layer

#if 0 //undev
    //case 1 : Beam-through
    if(AfterTGTHits>MinNumOfHits){
      Int_t Trackflag = 1*1 + 2*1 + 4*0 + 8*1; // isBeam, isK18, isKurama, isAccidental
      HelixTrackSearch(Trackflag, K18Tracks, ClContK18, TrackCont, TrackContFailed, MinNumOfHits+3);
    }
#endif

#if DebugDisp
    std::cout<<FUNC_NAME+" K18 #Hits in the window : "<<Hits<<" before the target : "<<BeforeTGTHits<<std::endl;
    std::cout<<FUNC_NAME+" Fitting the hits upstream of the target, # hits : "<<BeforeTGTHits<<std::endl;
#endif

    if(BeforeTGTHits<3){
      delete track;
      continue;
    }

    Int_t Trackflag = 1*1 + 2*1 + 4*0 + 8*0; // isBeam, isK18, isKurama, isAccidental
    track->SetFlag(Trackflag);
    track->SetSearchTime(0.);

    //Track fitting processes
    FitTrack(track, K18Tracks, ClContK18, TrackCont, TrackContFailed, 3);
    ResetHoughFlag(ClContK18);
  } //nt

}

//_____________________________________________________________________________
void
AccidentalBeamSearchTemp(const std::vector<TPCClusterContainer>& ClCont,
			 std::vector<TPCLocalTrackHelix*>& TrackCont,
			 std::vector<TPCLocalTrackHelix*>& TrackContFailed,
			 Int_t MinNumOfHits /*=10*/)
{

  Double_t Ywindow = 10.;
  Double_t x_min = -30.; Double_t x_max = 30.;
  Double_t y_min = -50.; Double_t y_max = 50.;

  hist_y->Reset();
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(hit->GetHoughFlag()>0) continue;
      TVector3 pos = cl->GetPosition();
      double x = pos.X();
      double y = pos.Y();
      if(x_min<x && x<x_max && y_min<y && y<y_max) continue;
      else hist_y->Fill(y);
    }
  }

  TSpectrum spec(30);
  double sig=1, th=0.2;
  const int npeaks = spec.Search(hist_y, sig,"goff",th);
  double* peakpos = spec.GetPositionX();
#if DebugEvDisp
  hist_y->Draw("colz");
  gPad->Modified();
  gPad->Update();
  c1.Print(Form("c%d.pdf",cannum));
  getchar();
  cannum++;
#endif

  for(Int_t candidates=0; candidates<npeaks; candidates++){
    std::vector<TPCClusterContainer> ClContBeam(NumOfLayersTPC);
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()>0) continue;
	TVector3 pos = cl->GetPosition();
	if(TMath::Abs(peakpos[candidates] - pos.Y())<Ywindow){
	  ClContBeam[layer].push_back(cl);
	}
      } //ci
    } //layer
    Int_t Trackflag = 1*1 + 2*0 + 4*0 + 8*1; // isBeam, isK18, isKurama, isAccidental
    HelixTrackSearch(Trackflag, AccidentalBeams, ClContBeam, TrackCont, TrackContFailed, MinNumOfHits);
    ResetHoughFlag(ClContBeam);
  } //candidates

}

//_____________________________________________________________________________
Int_t
LocalTrackSearchHelix(const std::vector<TPCClusterContainer>& ClCont,
		      std::vector<TPCLocalTrackHelix*>& TrackCont,
		      std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		      bool Exclusive,
		      Int_t MinNumOfHits /*=8*/)
{

  //Track searching and fitting
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);

#if RemainHitTest
  TestRemainingHits(ClCont, TrackCont, MinNumOfHits);
#endif

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

  CalcTracks(TrackCont);
  CalcTracks(TrackContFailed);
  if(Exclusive) ExclusiveTracking(TrackCont);

  return TrackCont.size();
}

//_____________________________________________________________________________
Int_t
LocalTrackSearchHelix(std::vector<std::vector<TVector3>> K18VPs,
		      std::vector<std::vector<TVector3>> KuramaVPs,
		      const std::vector<TPCClusterContainer>& ClCont,
		      std::vector<TPCLocalTrackHelix*>& TrackCont,
		      std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		      bool Exclusive,
		      Int_t MinNumOfHits /*=8*/)
{

  XZhough_x.clear();
  XZhough_y.clear();
  XZhough_z.clear();
  Yhough_x.clear();
  Yhough_y.clear();

  //Track searching and fitting
  //for Accidental beams and K1.8 & Kurama tracks
  AccidentalBeamSearchTemp(ClCont, TrackCont, TrackContFailed, 12);
  K18TrackSearch(K18VPs, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  KuramaTrackSearch(KuramaVPs, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  //for scattered helix tracks
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);

#if RemainHitTest
  TestRemainingHits(ClCont, TrackCont, MinNumOfHits);
#endif

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

  CalcTracks(TrackCont);
  CalcTracks(TrackContFailed);
  if(Exclusive) ExclusiveTracking(TrackCont);

  return TrackCont.size();
}


//_____________________________________________________________________________
void
HoughTransformTest(const std::vector<TPCClusterContainer>& ClCont,
		   std::vector<TPCLocalTrack*>& TrackCont,
		   Int_t MinNumOfHits /*=8*/)
{

  const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

  XZhough_x.clear();
  XZhough_y.clear();
  Yhough_x.clear();
  Yhough_y.clear();

  bool prev_add = true;
  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    if(!prev_add) continue;
    prev_add = false;

#if DebugDisp
    std::cout<<FUNC_NAME+" tracki : "<<tracki<<std::endl;
#endif

    std::chrono::milliseconds sec;
    auto before_hough = std::chrono::high_resolution_clock::now();

    //Line Hough-transform on the XZ plane
    Double_t LinearPar[4]; Int_t MaxBinXZ[2];
    if(!tpc::HoughTransformLineXZ(ClCont, MaxBinXZ, LinearPar, MinNumOfHits)){
#if DebugDisp
      std::cout<<FUNC_NAME+" No more track candiate! tracki : "<<tracki<<std::endl;
#endif
      break;
    }

    Int_t MaxBinY[3];
    //Line Hough-transform on the YZ or YX plane
    if(TMath::Abs(LinearPar[2]) < 1) tpc::HoughTransformLineYZ(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);
    else tpc::HoughTransformLineYX(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrack *track = new TPCLocalTrack;
    track->SetHoughParam(LinearPar);
    track->SetParamUsingHoughParam();
    track->CalcClosestDist();

    //If two tracks are merged at the target, seperate them and recalculate params.
    Bool_t vtx_flag = false;
    prev_add = MakeLinearTrack(track, vtx_flag, ClCont, LinearPar, MaxHoughWindowY);
    if(!prev_add) continue;
    if(vtx_flag){
      std::vector<TVector3> gHitPos;
      GetTrackHitCont(track, gHitPos);
      tpc::HoughTransformLineXZ(gHitPos, MaxBinXZ, LinearPar, MinNumOfHits);
      if(TMath::Abs(LinearPar[2]) < 1) tpc::HoughTransformLineYZ(gHitPos, MaxBinY, LinearPar, MaxHoughWindowY);
      else tpc::HoughTransformLineYX(gHitPos, MaxBinY, LinearPar, MaxHoughWindowY);
      track->SetHoughParam(LinearPar);
      track->SetParamUsingHoughParam();
      track->CalcClosestDist();
      track->SetVtxFlag(track->Side(gHitPos[0]));
    }
    if(track -> GetNDF()<1){
      track->SetHoughFlag(Candidate);
      delete track;
      continue;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<XZhough_x.size(); ++i){
      Int_t bindiffXZ = TMath::Abs(MaxBinXZ[0] - XZhough_x[i]) + TMath::Abs(MaxBinXZ[1] - XZhough_y[i]);
      Int_t bindiffY = TMath::Abs(MaxBinY[0] - Yhough_x[i]) + TMath::Abs(MaxBinY[1] - Yhough_y[i]);
      if(bindiffXZ<=1 && bindiffY<=1){
	hough_flag = false;
#if DebugDisp
	std::cout<<"Previous hough bin on the XZ plane "<<i<<"th x: "
		 <<XZhough_x[i]<<", y: "<<XZhough_y[i]<<" on the vertical plane x: "
		 <<Yhough_x[i]<<", y: "<<Yhough_y[i]<<std::endl;
	std::cout<<"Current hough bin on the XZ plane "<<i<<"th x: "
		 <<MaxBinXZ[0]<<", y: "<<MaxBinXZ[1]<<" on the vertical plane x: "
		 <<MaxBinY[0]<<", y: "<<MaxBinY[1]<<std::endl;
#endif
      }
    }
    XZhough_x.push_back(MaxBinXZ[0]);
    XZhough_y.push_back(MaxBinXZ[1]);
    Yhough_x.push_back(MaxBinY[0]);
    Yhough_y.push_back(MaxBinY[1]);

    if(!hough_flag){
#if DebugDisp
      std::cout<<FUNC_NAME+" The same track is found by Hough-Transform : tracki : "<<tracki<<" hough cont size : "<<XZhough_x.size()<<std::endl;
#endif
      track->SetHoughFlag(Candidate);
      delete track;
      continue;
    }

    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_hough - before_hough);
    track->SetHoughFlag(tracki + 1);
    track->SetSearchTime(sec.count());
    TrackCont.push_back(track);
  }// tracki

  CalcTracks(TrackCont);
}

//_____________________________________________________________________________
void
HoughTransformTestHelix(const std::vector<TPCClusterContainer>& ClCont,
			std::vector<TPCLocalTrackHelix*>& TrackCont,
			Int_t MinNumOfHits /*=8*/)
{
  // HoughTransform binning
  const auto MaxHoughWindow = gUser.GetParameter("MaxHoughWindow");
  const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

  XZhough_x.clear();
  XZhough_y.clear();
  XZhough_z.clear();
  Yhough_x.clear();
  Yhough_y.clear();

  bool prev_add = true;
  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    if(!prev_add) continue;
    prev_add = false;

#if DebugDisp
    std::cout<<FUNC_NAME+" tracki : "<<tracki<<std::endl;
#endif

    //Circle Hough-transform
    std::chrono::milliseconds sec;
    auto before_hough = std::chrono::high_resolution_clock::now();
    Double_t HelixPar[5]; Int_t MaxBinXZ[3];
    if(!tpc::HoughTransformCircleXZ(ClCont, MaxBinXZ, HelixPar, MinNumOfHits)){
#if DebugDisp
      std::cout<<FUNC_NAME+" No more circle candiate! tracki : "<<tracki<<std::endl;
#endif
      break;
    }

    //Linear Hough-transform
    Int_t MaxBinY[3];
    tpc::HoughTransformLineYTheta(ClCont, MaxBinY, HelixPar, MaxHoughWindow);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    track->SetHoughParam(HelixPar);
    track->SetParamUsingHoughParam();
    track->CalcClosestDist();

    //If two tracks are merged at the target, seperate them and recalculate params.
    Bool_t vtx_flag = false;
    prev_add = MakeHelixTrack(track, vtx_flag, ClCont, HelixPar, MaxHoughWindow, MaxHoughWindowY);
    if(!prev_add) continue;
    if(vtx_flag){
      std::vector<TVector3> gHitPos;
      GetTrackHitCont(track, gHitPos);
      tpc::HoughTransformCircleXZ(gHitPos, MaxBinXZ, HelixPar, MinNumOfHits);
      tpc::HoughTransformLineYTheta(gHitPos, MaxBinY, HelixPar, MaxHoughWindow);
      track->SetHoughParam(HelixPar);
      track->SetParamUsingHoughParam();
      track->CalcClosestDist();
      track->SetVtxFlag(track->Side(gHitPos[0]));
    }
    if(track -> GetNDF()<1){
      track->SetHoughFlag(Candidate);
      delete track;
      continue;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<XZhough_x.size(); ++i){
      Int_t bindiffXZ = TMath::Abs(MaxBinXZ[0] - XZhough_x[i]) + TMath::Abs(MaxBinXZ[1] - XZhough_y[i]) + TMath::Abs(MaxBinXZ[2] - XZhough_z[i]);
      Int_t bindiffY = TMath::Abs(MaxBinY[0] - Yhough_x[i]) + TMath::Abs(MaxBinY[1] - Yhough_y[i]);
      if(bindiffXZ<=1 && bindiffY<=1){
	hough_flag = false;
#if DebugDisp
	std::cout<<"Previous hough bin on the XZ plane "<<i<<"th x: "
		 <<XZhough_x[i]<<", y: "<<XZhough_y[i]<<", z: "<<XZhough_z[i]<<" on the vertical plane x: "
		 <<Yhough_x[i]<<", y: "<<Yhough_y[i]<<std::endl;
	std::cout<<"Current hough bin on the XZ plane "<<i<<"th x: "
		 <<MaxBinXZ[0]<<", y: "<<MaxBinXZ[1]<<", z: "<<MaxBinXZ[2]<<" on the vertical plane x: "
		 <<MaxBinY[0]<<", y: "<<MaxBinY[1]<<std::endl;
#endif
      }
    }
    XZhough_x.push_back(MaxBinXZ[0]);
    XZhough_y.push_back(MaxBinXZ[1]);
    XZhough_z.push_back(MaxBinXZ[1]);
    Yhough_x.push_back(MaxBinY[0]);
    Yhough_y.push_back(MaxBinY[1]);

    if(!hough_flag){
#if DebugDisp
      std::cout<<FUNC_NAME+" The same track is found by Hough-Transform : tracki : "<<tracki<<" hough cont size : "<<XZhough_x.size()<<std::endl;
#endif
      track->SetHoughFlag(Candidate);
      delete track;
      continue;
    }

    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_hough - before_hough);

    track->SetHoughFlag(tracki + 1);
    track->SetSearchTime(sec.count());
    track->FinalizeTrack();
    TrackCont.push_back(track);
  }//tracki

  CalcTracks(TrackCont);
}

} //namespace tpc
