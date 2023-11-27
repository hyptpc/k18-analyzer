// -*- C++ -*-

/*
//Comment by Wooseung

Track finding process
1. Hough-transform on the XZ plane
2. Hough-transform on the vertical plane
    -> Get initial track parameters.
3. Check residual(Houghdist) in the Hough-space and
make a inital track with clusters with small residual.

Track fitting process
1. First fitting with inital track.
2. If fitting is succedeed, check the residual of other cluster with the track.
If residual is acceptable, add cluster into the track.
3. If the track is extended by added clusters, then do fitting again. (recursive)

Especially near the target, two almost parallel tracks can be reconiged as a single track. So, there is a treatment to separate them.
track -> SeparateTracksAtTraget();

LocalTrackSearch/Helix are main functions of TPC tracking.
HoughTransformTest/Helix functions have only track finding algorithm for Hough-transform performace test. (no fitting)

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
#include "RootHelper.hh"

#define DebugDisp 0

namespace
{
  const Double_t Const = 0.299792458;
  const auto qnan = TMath::QuietNaN();
  const auto& gUser = UserParamMan::GetInstance();

  const Int_t    MaxNumOfTrackTPC = 20;
  //const Double_t KuramaXZWindow = 30.;
  //const Double_t KuramaYWindow = 10.;
  const Double_t KuramaXZWindow = 35.;
  const Double_t KuramaYWindow = 20.;
  const Double_t K18XZWindow = 10.5;
  //const Double_t K18YWindow = 10.;
  const Double_t K18YWindow = 15.;

  // Houghflags
  const Int_t GoodForTracking = 100;
  const Int_t K18Tracks = 200;
  //const Int_t KuramaTracks = 300;
  const Int_t BadForTracking = 400;
  const Int_t Candidate = 1000;

  // Tracks in the Hough-Space
  std::vector<Double_t> XZhough_x;
  std::vector<Double_t> XZhough_y;
  std::vector<Double_t> XZhough_z;
  std::vector<Double_t> Yhough_x;
  std::vector<Double_t> Yhough_y;

  //________________________ _____________________________________________________
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
  //reset houghflag of remain clusters
  void
  ResetHoughFlag(const std::vector<TPCClusterContainer>& ClCont, Int_t flagID=Candidate)
  {

    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()==flagID) hit->SetHoughFlag(0);
      } //ci
    } //layer
  }

  //_____________________________________________________________________________
  template <typename T> void
  GetTrackClCont(T* track, std::vector<TVector3>& gHitPos)
  {

    Int_t n = track->GetNHit();
    for(Int_t i=0; i<n; ++i){
      TPCLTrackHit *hitp = track->GetHit(i);
      TVector3 pos = hitp->GetLocalHitPos();
      gHitPos.push_back(pos);
    }
  }

  //_____________________________________________________________________________
  //Convert a straight line track into helix track
  Bool_t
  ConvertTrack(TPCLocalTrack *LinearTrack, TPCLocalTrackHelix *HelixTrack){

    Int_t n = LinearTrack->GetNHit();
    for(Int_t i=0; i<n; ++i){
      TPCLTrackHit *hitp = LinearTrack->GetHit(i);
      HelixTrack->AddTPCHit(hitp);
    }
    Double_t LinearPar[4];
    LinearTrack->GetParam(LinearPar);
    Bool_t status = HelixTrack->ConvertParam(LinearPar);

    delete LinearTrack;
    return status;
  }
}

//_____________________________________________________________________________
namespace tpc
{

//_____________________________________________________________________________
template <typename T> Bool_t
AddClusters(T* Track, const std::vector<TPCClusterContainer>& ClCont)
{

  Bool_t status = false;

  //Residual check with other hits
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(hit->GetHoughFlag()>0) continue;
      Double_t resi=0.;
      if(Track->ResidualCheck(hit, resi)){
	Int_t vtxflag = Track -> GetVtxFlag();
	TVector3 pos = hit -> GetPosition();
	Int_t side = Track -> Side(pos);
	//Vertex inside the target : vtxflag = -1 or 1 / outside vtxflag = 0
	//if track and new cluster are on the same side and vertex in the target : vtxflag*side = 1
	//if vertex is outside of the target : vtxflag*side = 0
	if(vtxflag*side >= 0){
	  Track->AddTPCHit(new TPCLTrackHit(hit));
	  status = true;
	}
      }
    } //ci
  } //layer

  return status;
}

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
  Track->Print(FUNC_NAME+" Initial track after Hough-transform");
  //Track->Print(FUNC_NAME+" Initial track after Hough-transform", true);
#endif

  //first fitting
  if(Track->DoFit()){ // MinNumOfHits cut is not applied in inital tracking
    auto first_fit = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(first_fit - fit_start);
    Track->SetClustersHoughFlag(Candidate);
    Track->SetFitTime(sec.count());
    Track->SetFitFlag(1);

#if DebugDisp
    Track->Print(FUNC_NAME+" First fitting is succeeded");
    //Track->Print(FUNC_NAME+" First fitting is succeeded", true);
#endif

    //1st fitting is succeeded.
    T *ExtendedTrack = new T(Track);
    if(!AddClusters(ExtendedTrack, ClCont)){ //case1. No more cluster to add. -> Move the track into containers
      delete ExtendedTrack;
      if(Track->GetNHit() >= MinNumOfHits){ //track w/ enough clusters
	Track->SetClustersHoughFlag(Houghflag);
	TrackCont.push_back(Track);
      }
      else{
	Track->SetClustersHoughFlag(BadForTracking); //track w/ few clusters
	Track->SetFlag(0);
	TrackContFailed.push_back(Track);
      }
#if DebugDisp
      Track->Print(FUNC_NAME+" No more cluster to add");
      //Track->Print(FUNC_NAME+" No more cluster to add", true);
#endif
    }
    else{ //case2. More clusters are added to the track -> Do track fitting process again.
      auto hitadd = std::chrono::high_resolution_clock::now();
#if DebugDisp
      ExtendedTrack->Print(FUNC_NAME+" After adding clusters");
      //ExtendedTrack->Print(FUNC_NAME+" After adding clusters", true);
#endif
      if(ExtendedTrack->DoFit(MinNumOfHits)){ //2nd fitting success
	delete Track;
	auto second_fit = std::chrono::high_resolution_clock::now();
	sec = std::chrono::duration_cast<std::chrono::milliseconds>(second_fit - hitadd);
	ExtendedTrack->SetFitTime(sec.count());
	ExtendedTrack->SetFitFlag(2);
	ExtendedTrack->SetClustersHoughFlag(Houghflag);
	TrackCont.push_back(ExtendedTrack);

#if DebugDisp
	ExtendedTrack->Print(FUNC_NAME+" Fitting success after adding clusters (residual check)");
	//ExtendedTrack->Print(FUNC_NAME+" Fitting success after adding clusters (residual check)", true);
#endif
      }
      else{ //2nd fitting failure -> Save the 1st fitting result
	delete ExtendedTrack;
	//push 1st fitting track
	auto hitadd_fail = std::chrono::high_resolution_clock::now();
	sec = std::chrono::duration_cast<std::chrono::milliseconds>(hitadd_fail - hitadd);
	Track->SetFitTime(sec.count());
	Track->SetFitFlag(3);
	if(Track->GetNHit() >= MinNumOfHits){ //track w/ enough clusters
	  Track->SetClustersHoughFlag(Houghflag);
	  TrackCont.push_back(Track);
	}
	else{
	  Track->SetClustersHoughFlag(BadForTracking); //track w/ few clusters
	  TrackContFailed.push_back(Track);
	}

#if DebugDisp
	Track->Print(FUNC_NAME+" 2nd fitting failure!");
	//Track->Print(FUNC_NAME+" 2nd fitting failure!", true);
#endif
      }
    }
  }
  else{ //First fitting is failed. (reset clusters' Hough flag)
    auto after_1stfit_fail = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_1stfit_fail - fit_start);
    Track->SetFitTime(sec.count());
    Track->SetFitFlag(4);
    Track->SetClustersHoughFlag(0);
    Track->SetFlag(0);
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
  /*
//_____________________________________________________________________________
//if there is a charge or momentum constraint, they can be applied in tracking.
template <typename T> void
FitTrack(T* Track, Int_t Houghflag,
	 const std::vector<TPCClusterContainer>& ClCont,
	 std::vector<T*>& TrackCont,
	 std::vector<T*>& TrackContFailed,
	 Double_t ChargeConstraint, Double_t RKHelixParam[5],
	 Int_t MinNumOfHits)
{

  std::chrono::milliseconds sec;
  auto fit_start = std::chrono::high_resolution_clock::now();

#if DebugDisp
  Track->Print(FUNC_NAME+" Initial track");
  //Track->Print(FUNC_NAME+" Initial track", true);
#endif

  //first fitting
  if(Track->DoFit(ChargeConstraint, RKHelixParam)){ // No MinNumOfHits cut for the inital Hough track
    auto first_fit = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(first_fit - fit_start);
    Track->SetClustersHoughFlag(Candidate);
    Track->SetFitTime(sec.count());
    Track->SetFitFlag(1);

#if DebugDisp
    Track->Print(FUNC_NAME+" First fitting is succeeded");
    //Track->Print(FUNC_NAME+" First fitting is succeeded", true);
#endif

    //1st fitting is succeeded.
    T *ExtendedTrack = new T(Track);
    if(!AddClusters(ExtendedTrack, ClCont)){ //No more cluster to add
      delete ExtendedTrack;
      if(Track->GetNHit() >= MinNumOfHits){ //track w/ enough clusters
	Track->SetClustersHoughFlag(Houghflag);
	TrackCont.push_back(Track);
      }
      else{
	Track->SetClustersHoughFlag(BadForTracking); //track w/ few clusters
	Track->SetFlag(0);
	TrackContFailed.push_back(Track);
      }
#if DebugDisp
      Track->Print(FUNC_NAME+" No added cluster");
      //Track->Print(FUNC_NAME+" No added cluster", true);
#endif
    }
    else{ //More clusters are added to the track
      auto hitadd = std::chrono::high_resolution_clock::now();
#if DebugDisp
      ExtendedTrack->Print(FUNC_NAME+" After cluster adding");
      //ExtendedTrack->Print(FUNC_NAME+" After cluster adding", true);
#endif
      if(ExtendedTrack->DoFit(ChargeConstraint, RKHelixParam, MinNumOfHits)){ //2nd fitting success
	delete Track;
	auto second_fit = std::chrono::high_resolution_clock::now();
	sec = std::chrono::duration_cast<std::chrono::milliseconds>(second_fit - hitadd);
	ExtendedTrack->SetFitTime(sec.count());
	ExtendedTrack->SetFitFlag(2);
	ExtendedTrack->SetClustersHoughFlag(Houghflag);
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
	  Track->SetClustersHoughFlag(Houghflag);
	  TrackCont.push_back(Track);
	}
	else{
	  Track->SetClustersHoughFlag(BadForTracking); //track w/ few clusters
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
    Track->SetClustersHoughFlag(0);
    Track->SetFlag(0);
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
  */
//_____________________________________________________________________________
Int_t
LocalTrackSearch(const std::vector<TPCClusterContainer>& ClCont,
		 std::vector<TPCLocalTrack*>& TrackCont,
		 std::vector<TPCLocalTrack*>& TrackContFailed,
		 Bool_t Exclusive,
		 Int_t MinNumOfHits)
{

  const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

  XZhough_x.clear();
  XZhough_y.clear();
  Yhough_x.clear();
  Yhough_y.clear();

  Bool_t prev_add = true;
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
    Int_t MaxBinY[2];
    if(TMath::Abs(LinearPar[2]) < 1) tpc::HoughTransformLineYZ(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);
    else tpc::HoughTransformLineYX(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrack *track = new TPCLocalTrack;
    track->SetParam(LinearPar);

    //If two tracks are merged at the target, separate them and recalculate params.
    Bool_t vtx_flag;
    prev_add = MakeLinearTrack(track, vtx_flag, ClCont, LinearPar, MaxHoughWindowY);
    if(!prev_add) continue;
    if(!track -> IsGoodForTracking() || !vtx_flag){
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

  //Check Hough-distance and add hits
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(!hit) continue;
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
	status = true;
      }
    } //ci
  } //layer

  if(status){
    //Vtx in the target, need to check whether two tracks are merged or not
    //VtxFlag = (Track -> SeparateClustersWithGap() && Track -> SeparateTracksAtTarget());
    Track -> SetClustersHoughFlag(Candidate);
    VtxFlag = Track -> SeparateTracksAtTarget();
    if(VtxFlag && AddClusters(Track, ClCont)) Track -> SetClustersHoughFlag(Candidate);
  }

#if DebugDisp
  if(status) Track->Print(FUNC_NAME+" Initial track after track finding");
//if(status) Track->Print(FUNC_NAME+" Initial track after track finding", true);
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

  //Check Hough-distance and add hits
  Int_t id = 0;
  Double_t theta0 = 0; Double_t prev_theta = 0.; Bool_t thetaflip = false;
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(!hit) continue;
      if(hit->GetHoughFlag()>0) continue;
      TVector3 pos = cl->GetPosition();
      Double_t tmpx = -pos.x();
      Double_t tmpy = pos.z() - tpc::ZTarget;
      Double_t tmpz = pos.y();
      Double_t r_cal = TMath::Sqrt(pow(tmpx - HelixPar[0], 2) + pow(tmpy - HelixPar[1], 2));
      Double_t dist = TMath::Abs(r_cal - HelixPar[3]);
      if(dist < MaxHoughWindow){
	Double_t tmpt = TMath::ATan2(tmpy - HelixPar[1], tmpx - HelixPar[0]);
	Double_t tmp_xval = HelixPar[3]*tmpt;
	Double_t distY = TMath::Abs(HelixPar[4]*tmp_xval - tmpz + HelixPar[2])/TMath::Sqrt(pow(HelixPar[4], 2) + 1.);
	if(distY < MaxHoughWindowY){
	  //Check theta flip
	  if(id==0) theta0 = tmpt;
	  if(TMath::Abs(prev_theta - tmpt) > TMath::Pi()) thetaflip = true;
	  prev_theta = tmpt;

	  //Add hit into the track
	  hit->SetHoughDist(dist);
	  hit->SetHoughDistY(distY);
	  Track->AddTPCHit(new TPCLTrackHit(hit));

	  id++;
	  status = true;
	} //distY
      } //dist
    } //ci
  } //layer

  if(thetaflip) Track -> SetIsThetaFlip();
  if(status){
    //Vtx in the target, need to check whether two tracks are merged or not
    Track -> CalcThetaHits();
    Track -> SetClustersHoughFlag(Candidate);
    //VtxFlag = (Track -> SeparateClustersWithGap() && Track -> SeparateTracksAtTarget());
    VtxFlag = Track -> SeparateTracksAtTarget();
    if(VtxFlag && AddClusters(Track, ClCont)){
      Double_t par[5];
      Track -> GetParam(par);
      VtxFlag = Track -> DoPreFit(par);
      Track -> SetClustersHoughFlag(Candidate);
    }
  }

#if DebugDisp
    Track->Print(FUNC_NAME+" Inital track after track finding");
#endif

  if(!status) delete Track;
  return status;
}

//_____________________________________________________________________________
void
HelixTrackSearch(Int_t Trackflag, Int_t Houghflag,
		 const std::vector<TPCClusterContainer>& ClCont,
		 std::vector<TPCLocalTrackHelix*>& TrackCont,
		 std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		 Int_t MinNumOfHits)
{

  // HoughTransform binning
  const auto MaxHoughWindow = gUser.GetParameter("MaxHoughWindow");
  const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

  Bool_t prev_add = true;
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
    Int_t MaxBinY[2];
    tpc::HoughTransformLineYTheta(ClCont, MaxBinY, HelixPar, MaxHoughWindow);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    track->SetParam(HelixPar);
    track->SetFlag(Trackflag);

    //If two tracks are merged at the target, separate them and recalculate params.
    Bool_t vtx_flag;
    prev_add = MakeHelixTrack(track, vtx_flag, ClCont, HelixPar, MaxHoughWindow, MaxHoughWindowY);
    if(!prev_add) continue;
    if(!track -> IsGoodForTracking() || !vtx_flag){
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
      track->SetFitFlag(0);
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
		  std::vector<TPCLocalTrackHelix*>& TrackCont,
		  std::vector<TPCLocalTrackHelix*>& TrackContVP)
{

  if(VPs.size()==0) return;
  for(Int_t nt=0; nt<VPs.size(); nt++){
#if DebugDisp
    std::cout<<FUNC_NAME+" Kurama Track "<<nt<<"/"<<VPs.size()<<std::endl;
#endif

    //Reconstruct the kurama track in the HypTPC
    TPCLocalTrackHelix *trackref = new TPCLocalTrackHelix();
    for(Int_t i=0; i<VPs[nt].size(); i++){
      TVector3 pos = VPs[nt][i];
      trackref->AddVPHit(pos);
#if DebugDisp
      std::cout<<FUNC_NAME+" Kurama VP "<<i<<"th pos : "<<VPs[nt][i]<<std::endl;
#endif
    } //i
    trackref->SetIsKurama();
    trackref->SetTrackID(nt);
    if(!trackref->DoVPFit()){
      delete trackref;
      continue;
    }

    Double_t RKHelixParam[5];
    trackref->GetParam(RKHelixParam);
    TrackContVP.push_back(trackref);

#if DebugDisp
    std::cout<<FUNC_NAME+" Kurama VP Helix cx : "<<trackref->Getcx()<<" cy : "<<trackref->Getcy()<<" z0 : "<<trackref->Getz0()<<" r : "<<trackref->Getr()<<" dz : "<<trackref->Getdz()<<std::endl;
    std::cout<<FUNC_NAME+" Kurama VP Helix p : "<<trackref->Getr()*Const<<std::endl;
#endif

    Int_t Trackflag = 1*0 + 2*0 + 4*1 + 8*0; // isBeam, isK18, isKurama, isAccidental
    for(auto& track: TrackCont){
      if(track){
	//Check whether all track custers within the window along the Kurama track
	if(track -> GetIsK18()==1) continue;
	Int_t nh = track->GetNHit();
	for(Int_t ih=0; ih<nh; ++ih){
	  TPCLTrackHit *hit = track -> GetHit( ih );
	  if( !hit ) continue;
	  const TVector3& hitpos = hit->GetLocalHitPos();
	  if(!trackref->ResidualCheck(hitpos, KuramaXZWindow, KuramaYWindow)) break;
	  //If all clusters are in the window, mark that helix track.
	  if(ih==nh-1){
	    track->AddTrackIDCandidate(nt);
	    track->SetFlag(Trackflag);
	  }
	}
      }
    }
  } //nt

}

//_____________________________________________________________________________
void
K18TrackSearch(std::vector<std::vector<TVector3>> VPs,
	       const std::vector<TPCClusterContainer>& ClCont,
	       std::vector<TPCLocalTrackHelix*>& TrackCont,
	       std::vector<TPCLocalTrackHelix*>& TrackContVP,
	       Int_t MinNumOfHits /*=3*/)
{
  if(VPs.size()==0) return;
  for(Int_t nt=0; nt<VPs.size(); nt++){
#if DebugDisp
    std::cout<<FUNC_NAME+" K18 Track "<<nt<<"/"<<VPs.size()<<std::endl;
#endif

    std::chrono::milliseconds sec;
    auto before_track_search = std::chrono::high_resolution_clock::now();

    TPCLocalTrackHelix *trackref = new TPCLocalTrackHelix();
    for(Int_t i=0; i<VPs[nt].size(); i++){
      TVector3 pos = VPs[nt][i];
      trackref->AddVPHit(pos);
#if DebugDisp
      std::cout<<FUNC_NAME+" K18 RK "<<i<<"th pos : "<<VPs[nt][i]<<std::endl;
#endif
    } //i
    trackref->SetIsK18();
    trackref->SetIsBeam();
    trackref->DoVPFit();

#if DebugDisp
    std::cout<<FUNC_NAME+" K18 VP Helix cx : "<<trackref->Getcx()<<" cy : "<<trackref->Getcy()<<" z0 : "<<trackref->Getz0()<<" r : "<<trackref->Getr()<<" dz : "<<trackref->Getdz()<<std::endl;
    std::cout<<FUNC_NAME+" K18 VP Helix p : "<<trackref->Getr()*Const<<std::endl;
#endif

    Int_t BeforeTGTHits = 0;
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    std::vector<TPCClusterContainer> ClContK18(NumOfLayersTPC);
    for(Int_t layer=0; layer<10; layer++){ //inner layers
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
	}
      } //ci
      if(minresi<100.){
	TPCHit* hit = ClCont[layer][id]->GetMeanHit();
	hit->SetHoughDist(qnan);
	track->AddTPCHit(new TPCLTrackHit(hit));
	BeforeTGTHits++;
      }
    } //layer

    Double_t RKHelixParam[5];
    trackref->GetParam(RKHelixParam);
    trackref->SetNclBeforeTgt(BeforeTGTHits);
    TrackContVP.push_back(trackref);

#if DebugDisp
    std::cout<<FUNC_NAME+" Clusters upstream of the target and within the window : "<<BeforeTGTHits<<std::endl;
    track->Print(FUNC_NAME+ " K18 track candidate");
#endif

    if(BeforeTGTHits<2){
      delete track;
      continue;
    }

    Int_t Trackflag = 1*1 + 2*1 + 4*0 + 8*0; // isBeam, isK18, isKurama, isAccidental
    auto after_track_search = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_track_search - before_track_search);
    track->SetSearchTime(sec.count());
    track->SetFlag(Trackflag);
    track->SetTrackID(nt); //Marking K18 track id

    //Track fitting processes
    auto fit_start = std::chrono::high_resolution_clock::now();
    if(track->DoFit(RKHelixParam, MinNumOfHits)){ //with constraints
      auto first_fit = std::chrono::high_resolution_clock::now();
      sec = std::chrono::duration_cast<std::chrono::milliseconds>(first_fit - fit_start);
      track->SetFitTime(sec.count());
      track->SetFitFlag(1);
      track->SetClustersHoughFlag(Candidate);

#if DebugDisp
      track->Print(FUNC_NAME+" First fitting is succeeded");
      //track->Print(FUNC_NAME+" First fitting is succeeded", true);
#endif

      TPCLocalTrackHelix *ExtendedTrack = new TPCLocalTrackHelix(track);
      //Residual check with other hits
      if(!AddClusters(ExtendedTrack, ClContK18)){ //No more cluster to add
	delete ExtendedTrack;
	track->SetClustersHoughFlag(K18Tracks);
	TrackCont.push_back(track);
#if DebugDisp
	track->Print(FUNC_NAME+" No added cluster");
	//track->Print(FUNC_NAME+" No added cluster", true);
#endif
      }
      else{ //More clusters are added Into the track
	auto hitadd = std::chrono::high_resolution_clock::now();
#if DebugDisp
	ExtendedTrack->Print(FUNC_NAME+" After cluster adding");
	//ExtendedTrack->Print(FUNC_NAME+" After cluster adding", true);
#endif
	if(ExtendedTrack->DoFit(RKHelixParam, MinNumOfHits)){ //2nd fitting success
	  delete track;
	  auto second_fit = std::chrono::high_resolution_clock::now();
	  sec = std::chrono::duration_cast<std::chrono::milliseconds>(second_fit - hitadd);
	  ExtendedTrack->SetFitTime(sec.count());
	  ExtendedTrack->SetFitFlag(2);
	  ExtendedTrack->SetClustersHoughFlag(K18Tracks);
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
	  track->SetFitTime(sec.count());
	  track->SetFitFlag(3);
	  track->SetClustersHoughFlag(K18Tracks);
	  TrackCont.push_back(track);

#if DebugDisp
	  track->Print(FUNC_NAME+" 2nd fitting failure!");
	  //track->Print(FUNC_NAME+" 2nd fitting failure!", true);
#endif
	}
      }
    }
    else{ //First fitting is failed. Reset clusters' Hough flag
#if DebugDisp
      track->Print(FUNC_NAME+" Fitting is failed");
      //track->Print(FUNC_NAME+" Fitting is failed", true);
#endif
      delete track;
    }
    ResetHoughFlag(ClContK18);
  } //nt
}

//_____________________________________________________________________________
Int_t
LocalTrackSearchHelix(const std::vector<TPCClusterContainer>& ClCont,
		      std::vector<TPCLocalTrackHelix*>& TrackCont,
		      std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		      Bool_t Exclusive,
		      Int_t MinNumOfHits)
{

  //Track finding and fitting
  HighMomHelixTrackSearch(ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);

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
		      std::vector<TPCLocalTrackHelix*>& TrackContVP,
		      std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		      Bool_t Exclusive,
		      Int_t MinNumOfHits /*=8*/)
{

  XZhough_x.clear();
  XZhough_y.clear();
  XZhough_z.clear();
  Yhough_x.clear();
  Yhough_y.clear();

  //Track finding and fitting
  //for K1.8 & Kurama tracks
  K18TrackSearch(K18VPs, ClCont, TrackCont, TrackContVP); //default NimNumOfHits = 3
  //for scattered helix tracks
  HighMomHelixTrackSearch(ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  KuramaTrackSearch(KuramaVPs, TrackCont, TrackContVP);

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

  CalcTracks(TrackCont);
  CalcTracks(TrackContVP);
  CalcTracks(TrackContFailed); //Tracking failed cases
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

  Bool_t prev_add = true;
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

    Int_t MaxBinY[2];
    //Line Hough-transform on the YZ or YX plane
    if(TMath::Abs(LinearPar[2]) < 1) tpc::HoughTransformLineYZ(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);
    else tpc::HoughTransformLineYX(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrack *track = new TPCLocalTrack;
    track->SetParam(LinearPar);

    //If two tracks are merged at the target, separate them and recalculate params.
    Bool_t vtx_flag;
    prev_add = MakeLinearTrack(track, vtx_flag, ClCont, LinearPar, MaxHoughWindowY);
    if(!prev_add) continue;
    if(!track -> IsGoodForTracking() || !vtx_flag){
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
      delete track;
      continue;
    }

    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_hough - before_hough);
    track->SetClustersHoughFlag(tracki+1);
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

  Bool_t prev_add = true;
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
    Int_t MaxBinY[2];
    tpc::HoughTransformLineYTheta(ClCont, MaxBinY, HelixPar, MaxHoughWindow);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    track->SetParam(HelixPar);

    //If two tracks are merged at the target, separate them and recalculate params.
    Bool_t vtx_flag;
    prev_add = MakeHelixTrack(track, vtx_flag, ClCont, HelixPar, MaxHoughWindow, MaxHoughWindowY);
    if(!prev_add) continue;
    if(!track -> IsGoodForTracking() || !vtx_flag){
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
      delete track;
      continue;
    }

    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_hough - before_hough);

    track->SetClustersHoughFlag(tracki + 1);
    track->SetSearchTime(sec.count());
    TrackCont.push_back(track);
  }//tracki

  CalcTracks(TrackCont);
}

//_____________________________________________________________________________
void
HighMomHelixTrackSearch(const std::vector<TPCClusterContainer>& ClCont,
			std::vector<TPCLocalTrackHelix*>& TrackCont,
			std::vector<TPCLocalTrackHelix*>& TrackContFailed,
			Int_t MinNumOfHits)
{

  const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

  std::vector<Double_t> tempXZhough_x;
  std::vector<Double_t> tempXZhough_y;
  std::vector<Double_t> tempYhough_x;
  std::vector<Double_t> tempYhough_y;

  Bool_t prev_add = true;
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
    Int_t MaxBinY[2];
    if(TMath::Abs(LinearPar[2]) < 1) tpc::HoughTransformLineYZ(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);
    else tpc::HoughTransformLineYX(ClCont, MaxBinY, LinearPar, MaxHoughWindowY);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrack *trackTemp = new TPCLocalTrack;
    trackTemp->SetParam(LinearPar);

    Bool_t vtx_flag;
    prev_add = MakeLinearTrack(trackTemp, vtx_flag, ClCont, LinearPar, MaxHoughWindowY);
    if(!prev_add) continue;
    if(!trackTemp -> IsGoodForTracking() || !vtx_flag){
      delete trackTemp;
      continue;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<tempXZhough_x.size(); ++i){
      Int_t bindiffXZ = TMath::Abs(MaxBinXZ[0] - tempXZhough_x[i]) + TMath::Abs(MaxBinXZ[1] - tempXZhough_y[i]);
      Int_t bindiffY = TMath::Abs(MaxBinY[0] - tempYhough_x[i]) + TMath::Abs(MaxBinY[1] - tempYhough_y[i]);
      if(bindiffXZ<=1 && bindiffY<=1){
	hough_flag = false;
#if DebugDisp
	std::cout<<"Previous hough bin on the XZ plane "<<i<<"th x: "
		 <<tempXZhough_x[i]<<", y: "<<tempXZhough_y[i]<<" on the vertical plane x: "
		 <<tempYhough_x[i]<<", y: "<<tempYhough_y[i]<<std::endl;
	std::cout<<"Current hough bin on the XZ plane "<<i<<"th x: "
		 <<MaxBinXZ[0]<<", y: "<<MaxBinXZ[1]<<" on the vertical plane x: "
		 <<MaxBinY[0]<<", y: "<<MaxBinY[1]<<std::endl;
#endif
      }
    }
    tempXZhough_x.push_back(MaxBinXZ[0]);
    tempXZhough_y.push_back(MaxBinXZ[1]);
    tempYhough_x.push_back(MaxBinY[0]);
    tempYhough_y.push_back(MaxBinY[1]);

    if(!hough_flag){
#if DebugDisp
      std::cout<<FUNC_NAME+" The same track is found by Hough-Transform : tracki : "<<tracki<<" hough cont size : "<<tempXZhough_x.size()<<std::endl;
#endif
      delete trackTemp;
      continue;
    }

    //Convert the temporary straight line track into helix track
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix;
    //track -> GetNDF();
    if(!ConvertTrack(trackTemp, track)){
#if DebugDisp
      track->Print(FUNC_NAME+ " Track converting is failed");
#endif
      delete track;
      continue;
    }
    else{ //temporary treatment
      if(AddClusters(track, ClCont)){
	//track -> GetNDF();
	Double_t par[5];
	track -> GetParam(par);
	track -> DoPreFit(par);
      }
    }

    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_hough - before_hough);
    track->SetSearchTime(sec.count());

    //Track fitting processes
    FitTrack(track, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  }// tracki
  ResetHoughFlag(ClCont);
  ResetHoughFlag(ClCont, BadForTracking);

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

}

} //namespace tpc
