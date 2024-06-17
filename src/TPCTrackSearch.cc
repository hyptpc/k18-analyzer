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
2. If fitting is succeeded, check the residual of other cluster with the track.
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
#include "TPCVertexHelix.hh"
#include "TPCCluster.hh"
#include "RootHelper.hh"
#include "DebugCounter.hh"
#define DebugDisp 0
#define FragmentedTrackTest 1
#define ReassignClusterTest 1
//#define RemainingClustersTest 1
#define RemainingClustersTest 0

namespace
{
  const auto qnan = TMath::QuietNaN();
  const auto& gUser = UserParamMan::GetInstance();
  const auto& gCounter = debug::ObjectCounter::GetInstance();
  const Int_t    MaxNumOfTrackTPC = 30;
  //const Double_t KuramaXZWindow = 35.; //ref
  const Double_t KuramaXZWindow = 45.;
  //const Double_t KuramaYWindow = 20.; //ref
  const Double_t KuramaYWindow = 30.; //ref
  const Double_t K18XZWindow = 10.5;
  //const Double_t K18YWindow = 10.;
  const Double_t K18YWindow = 15.;

  // Closest distance cut for vertex finding
  const Double_t VertexDistCut = 30.;

  // Maximum number of fitting steps
  //const Int_t MaxFitSteps = 5;
  const Int_t MaxFitSteps = 10;

  // Houghflags
  const Int_t GoodForTracking = 100;
  const Int_t K18Tracks = 200;
  //const Int_t KuramaTracks = 300;
  const Int_t BadHoughTransform = 300;
  const Int_t BadForTracking = 400;
  const Int_t Candidate = 1000;

  // Tracks in the Hough-Space
  std::vector<Double_t> XZhough_x;
  std::vector<Double_t> XZhough_y;
  std::vector<Double_t> XZhough_z;
  std::vector<Double_t> Yhough_x;
  std::vector<Double_t> Yhough_y;

  //_____________________________________________________________________________
  // Local Functions
  //_____________________________________________________________________________
  template <typename T> void
  CalcTracks(std::vector<T*>& TrackCont)
  {
    for(auto& track: TrackCont){
      track->Calculate();
    }
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

  template <typename T> void
  MarkingAccidentalTracks(std::vector<T*>& TrackCont)
  {
    for(auto& track: TrackCont){
      track->CheckIsAccidental();
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
  template <typename T> Bool_t
  IsGood(T* track, Int_t MinNumOfHits)
  {

    Bool_t status = false;
    if(track->IsBackward()) status = true;
    else if(track->GetNHit() >= MinNumOfHits) status = true; //track w/ enough clusters
    return status;
  }

  //_____________________________________________________________________________
  //Convert a straight line track into helix track
  Bool_t
  ConvertTrack(TPCLocalTrack *LinearTrack, TPCLocalTrackHelix *HelixTrack)
  {

    Int_t n = LinearTrack->GetNHit();
    for(Int_t i=0; i<n; ++i){
      TPCLTrackHit *hitp = LinearTrack->GetHit(i);
      HelixTrack->AddTPCHit(hitp);
    }
    Double_t LinearPar[4];
    LinearTrack->GetParam(LinearPar);
    Bool_t status = HelixTrack->ConvertParam(LinearPar);
    delete LinearTrack;
    if(status &&
       HelixTrack-> Getr() < 200.) status = false;
    return status;
  }

  //_____________________________________________________________________________
  template <typename T> Bool_t
  AddClusters(T* Track, const std::vector<TPCClusterContainer>& ClCont, Int_t HoughFlag=0)
  {

    Bool_t status = false;
    //Residual check with other hits
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()!=HoughFlag) continue;
	Double_t resi=0.;
	if(Track->IsGoodHitToAdd(hit, resi)){
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
	//else std::cout<<"resi "<<resi<<std::endl;
      } //ci
    } //layer
    return status;
  }
} //namespace

namespace tpc{
//_____________________________________________________________________________
template <typename T> Bool_t
FitStep(T* Track,
	const std::vector<TPCClusterContainer>& ClCont,
	std::vector<T*>& TrackContFailed,
	Int_t MinNumOfHits)
{

  Bool_t status = true;
  if(Track->DoFit(MinNumOfHits)){ // MinNumOfHits cut is not applied in inital tracking
    Track->SetClustersHoughFlag(Candidate);

#if DebugDisp
    Track->Print(FUNC_NAME+" Fitting is succeeded");
    //Track->Print(FUNC_NAME+" Fitting is succeeded", true);
#endif
  }
  else{ //Fitting is failed. (reset clusters' Hough flag)
    Track->SetClustersHoughFlag(0);
    Track->SetFlag(0);
    TrackContFailed.push_back(Track);

#if DebugDisp
    Track->Print(FUNC_NAME+" Fitting is failed");
#endif
    status = false;
  }

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
  Track->Print(FUNC_NAME+" An initial track for fitting");
  //Track->Print(FUNC_NAME+" An initial track for fitting", true);
#endif

  Int_t nstep = 0;
  while(true){
    T *ExtendedTrack = new T(Track);
    //Int_t thr_ncl = 0.5*MinNumOfHits;
    Int_t thr_ncl = 0;
    if(nstep==0) thr_ncl = 0; //The initial track can be short.
    else if(nstep > MaxFitSteps){ //Tracking is over
      delete ExtendedTrack;
      if(IsGood(Track, MinNumOfHits)){
#if DebugDisp
	Track->Print(FUNC_NAME+" Tracking is over. track is good for tracking");
	//Track->Print(FUNC_NAME+" Tracking is over. track is good for tracking", true);
#endif
	Track->SetClustersHoughFlag(Houghflag);
	TrackCont.push_back(Track);
      }
      else{
#if DebugDisp
	Track->Print(FUNC_NAME+" Tracking is over. track is not good for tracking");
	//Track->Print(FUNC_NAME+" Tracking is over. track is not good for tracking", true);
#endif
	Track->SetClustersHoughFlag(BadForTracking); //track w/ few clusters
	Track->SetFlag(0);
	TrackContFailed.push_back(Track);
      }
      break;
    }
    else if(!AddClusters(ExtendedTrack, ClCont)){
      //No more clusters for addding, then check for clusters of bad tracks.
      Bool_t add_badhough = AddClusters(ExtendedTrack, ClCont, BadHoughTransform);
      Bool_t add_badtracking = AddClusters(ExtendedTrack, ClCont, BadForTracking);
      if(!add_badhough && !add_badtracking){
	delete ExtendedTrack;
	if(IsGood(Track, MinNumOfHits)){
#if DebugDisp
	  Track->Print(FUNC_NAME+" track is good and no more cluster to add");
	  //Track->Print(FUNC_NAME+" track is good and no more cluster to add", true);
#endif
	  Track->SetClustersHoughFlag(Houghflag);
	  TrackCont.push_back(Track);
	}
	else{
#if DebugDisp
	  Track->Print(FUNC_NAME+" Track is not good and no more cluster to add");
	  //Track->Print(FUNC_NAME+" Track is not good and no more cluster to add", true);
#endif
	  Track->SetClustersHoughFlag(BadForTracking); //track w/ few clusters
	  Track->SetFlag(0);
	  TrackContFailed.push_back(Track);
	}
	break; //No more cluster to add
      }
    } // No more cluster to add

    //After adding more clusters, fitting starts.
    if(FitStep(ExtendedTrack, ClCont, TrackContFailed, thr_ncl)){
      delete Track;
      Track = ExtendedTrack; //Updates the track
#if DebugDisp
      std::cout<<"extension"<<std::endl;
#endif
    }
    else{
      if(nstep!=0 && IsGood(Track, MinNumOfHits)){
#if DebugDisp
	Track->Print(FUNC_NAME+" no more clusters to add");
	//Track->Print(FUNC_NAME+" no more clusters to add", true);
#endif
	Track->SetClustersHoughFlag(Houghflag);
	TrackCont.push_back(Track);
      }
      else{
#if DebugDisp
	std::cout<<"delete"<<std::endl;
#endif
	delete Track;
      }
      break; //Extended track's fitting is failed
    }

#if DebugDisp
    std::cout<<FUNC_NAME+" # of fitting steps for tracking: "<<nstep<<std::endl;
    std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<" #failed track : "<<TrackContFailed.size()<<std::endl;
    std::cout<<std::endl;
#endif

    nstep++;
  } //while

  auto fittingtime = std::chrono::high_resolution_clock::now();
  sec = std::chrono::duration_cast<std::chrono::milliseconds>(fittingtime - fit_start);
  if(Track) Track -> SetFitTime(sec.count());
  if(Track) Track -> SetFitFlag(nstep);

}

//_____________________________________________________________________________
Int_t
LocalTrackSearch(const std::vector<TPCClusterContainer>& ClCont,
		 std::vector<TPCLocalTrack*>& TrackCont,
		 std::vector<TPCLocalTrack*>& TrackContFailed,
		 Bool_t Exclusive,
		 Int_t MinNumOfHits)
{

  static const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

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
      track->SetClustersHoughFlag(BadHoughTransform);
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
      track->SetClustersHoughFlag(BadHoughTransform);
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
  ResetHoughFlag(ClCont, BadHoughTransform);

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

  CalcTracks(TrackCont);
  MarkingAccidentalTracks(TrackCont);
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
    Track -> SetClustersHoughFlag(Candidate);
    VtxFlag = (Track -> SeparateClustersWithGap() || Track -> SeparateTracksAtTarget());
    //VtxFlag = Track -> SeparateTracksAtTarget();
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

	//for the inital track, scanning theta within range (-pi, pi)
	Double_t tmpt = TMath::ATan2(tmpy - HelixPar[1], tmpx - HelixPar[0]);
	Double_t tmp_xval = HelixPar[3]*tmpt;
	Double_t distY = TMath::Abs(HelixPar[4]*tmp_xval - tmpz + HelixPar[2])/TMath::Sqrt(pow(HelixPar[4], 2) + 1.);
	Double_t distY_2pi = TMath::Abs(HelixPar[4]*(tmp_xval + 2.*TMath::Pi()*HelixPar[3]) - tmpz + HelixPar[2])/TMath::Sqrt(pow(HelixPar[4], 2) + 1.);
	Double_t distY_m2pi = TMath::Abs(HelixPar[4]*(tmp_xval - 2.*TMath::Pi()*HelixPar[3]) - tmpz + HelixPar[2])/TMath::Sqrt(pow(HelixPar[4], 2) + 1.);
	if(distY < MaxHoughWindowY){
	  //Add hit into the track
	  hit->SetHoughDist(dist);
	  hit->SetHoughDistY(distY);
	  Track->AddTPCHit(new TPCLTrackHit(hit));

	  id++;
	  status = true;
	} //distY
	else if(distY_2pi < MaxHoughWindowY){
	  //Add hit into the track
	  hit->SetHoughDist(dist);
	  hit->SetHoughDistY(distY_2pi);
	  Track->AddTPCHit(new TPCLTrackHit(hit));

	  id++;
	  status = true;
	}
	else if(distY_m2pi < MaxHoughWindowY){
	  //Add hit into the track
	  hit->SetHoughDist(dist);
	  hit->SetHoughDistY(distY_m2pi);
	  Track->AddTPCHit(new TPCLTrackHit(hit));

	  id++;
	  status = true;
	}
      } //dist
    } //ci
  } //layer

  if(status){
    Track -> SetClustersHoughFlag(Candidate);
    Track -> CalcHelixTheta();

    //If track has very large gap, then spilt it.
    if(Track -> SeparateClustersWithGap()) Track -> CalcHelixTheta();

    //If the vtx is in the target, need to check whether two tracks are fragmented tracks or independent tracks.
    VtxFlag = Track -> SeparateTracksAtTarget();
    if(VtxFlag && AddClusters(Track, ClCont)){
      Double_t par[5];
      Track -> GetParam(par);
      VtxFlag = Track -> DoPreFit(par);
      Track -> SetClustersHoughFlag(Candidate);
    }
  }

#if DebugDisp
  if(status) Track->Print(FUNC_NAME+" Initial track after track finding");
  //if(status) Track->Print(FUNC_NAME+" Initial track after track finding", true);
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
  static const auto MaxHoughWindow = gUser.GetParameter("MaxHoughWindow");
  static const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

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
      track->SetClustersHoughFlag(BadHoughTransform);
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
      track->SetClustersHoughFlag(BadHoughTransform);
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
  ResetHoughFlag(ClCont, BadHoughTransform);

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

}

//_____________________________________________________________________________
void
KuramaTrackSearch(std::vector<std::vector<TVector3>> VPs,
		  const std::vector<TPCClusterContainer>& ClCont,
		  std::vector<TPCLocalTrackHelix*>& TrackCont,
		  std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		  std::vector<TPCLocalTrackHelix*>& TrackContVP,
		  std::vector<TPCVertexHelix*>& VertexCont,
		  Bool_t Exclusive,
		  Int_t MinNumOfHits)
{
  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);

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
    std::cout<<FUNC_NAME+" Kurama VP Helix p : "<<trackref->Getr()*0.299792458<<std::endl;
#endif

    Int_t Trackflag = 1*0 + 2*0 + 4*1 + 8*0; // isBeam, isK18, isKurama, isAccidental
    Int_t id = 0;
    std::vector<Int_t> kurama_candidates;
    for(auto& track: TrackCont){
      if(track){
	Int_t ncl_downstream_tgt = 0; //#cluster after the target.

	//Check whether all track custers within the window along the Kurama track
	if(BeamThroughTPC || track -> GetIsK18()!=1){
	  Int_t nh = track->GetNHit();
	  for(Int_t ih=0; ih<nh; ++ih){
	    TPCLTrackHit *hit = track -> GetHit( ih );
	    if( !hit ) continue;
	    const TVector3& hitpos = hit->GetLocalHitPos();
	    if(!trackref->ResidualCheck(hitpos, KuramaXZWindow, KuramaYWindow)) break;
	    if(hitpos.Z()>tpc::ZTarget) ncl_downstream_tgt++;
	    //If all clusters are in the window, mark the track.
	    if(ih==nh-1 && ncl_downstream_tgt>0){
	      track->AddTrackIDCandidate(nt);
	      track->SetFlag(Trackflag);
	      kurama_candidates.push_back(id);
	    }
	  }
	}
	id++;
      }
    }

    //If there is two candidates, checking whether they are a single track or not.
    if(kurama_candidates.size()==2){
      Int_t trackid1 = kurama_candidates.at(0);
      Int_t trackid2 = kurama_candidates.at(1);

      Bool_t order = TrackCont[trackid1] -> GetNHit() >= TrackCont[trackid2] -> GetNHit() ?  true : false;
      if(!order){
	trackid1 = kurama_candidates.at(1);
	trackid2 = kurama_candidates.at(0);
      }

      TPCLocalTrackHelix *track1 = TrackCont[trackid1];
      TPCLocalTrackHelix *track2 = TrackCont[trackid2];
      TPCLocalTrackHelix *MergedTrack = new TPCLocalTrackHelix(track1);
      for(Int_t hit=0;hit<track2 -> GetNHit();hit++)
	MergedTrack -> AddTPCHit(new TPCLTrackHit(track2 -> GetHitInOrder(hit) -> GetHit()));
      //Check whether tracks belong to the same track or not
      if(!MergedTrack -> TestMergedTrack()){
	delete MergedTrack;
	continue;
      }

      //Two candidates are the same track.
      //Calulate the track and erase previous two fragmented tracks.
      std::vector<Int_t> candidates;
      candidates.push_back(trackid1);
      candidates.push_back(trackid2);

      Int_t prev_size = TrackCont.size();
      FitTrack(MergedTrack, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);

      Bool_t fitstatus = true;
      Int_t post_size = TrackCont.size();
      if(prev_size+1 == post_size){
	MergedTrack = TrackCont[post_size-1];
	MergedTrack -> Calculate();
	MergedTrack -> AddTrackIDCandidate(nt);
	MergedTrack -> SetFlag(Trackflag);
	if(Exclusive){
	  MergedTrack -> DoFitExclusive();
	  MergedTrack -> CalculateExclusive();
	}
      }
      else fitstatus = false;

      if(candidates.size()>0){
	std::sort(candidates.begin(), candidates.end());
	if(fitstatus){
	  for(Int_t i=0; i<candidates.size(); ++i){
	    Int_t trackID = candidates[i] - i;
	    TrackCont.erase(TrackCont.begin() + trackID);
	  }

	  //Vertex finding again with new tracks
	  del::ClearContainer(VertexCont);
	  VertexSearch(TrackCont, VertexCont);
	}
	else{
	  for(Int_t i=0; i<candidates.size(); ++i){
	    Int_t trackID = candidates[i] - i;
	    TrackCont[trackID] -> SetClustersHoughFlag(GoodForTracking);
	  }
	}
      }
    } //merging two tracks
  } //nt
}

//_____________________________________________________________________________
void
K18TrackSearch(std::vector<std::vector<TVector3>> VPs,
	       const std::vector<TPCClusterContainer>& ClCont,
	       std::vector<TPCLocalTrackHelix*>& TrackCont,
	       std::vector<TPCLocalTrackHelix*>& TrackContVP,
	       Int_t MinNumOfHits /*=2*/)
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
    std::cout<<FUNC_NAME+" K18 VP Helix p : "<<trackref->Getr()*0.299792458<<std::endl;
#endif

    Int_t BeforeTGTHits = 0;
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    //set min, max theta for checking a closest distance and residual
    track -> SetMint(trackref->GetMint());
    track -> SetMaxt(trackref->GetMaxt());

    trackref->DoVPFit();
    std::vector<TPCClusterContainer> ClContK18(NumOfLayersTPC);
    for(Int_t layer=0; layer<10; layer++){ //inner layers
      Double_t minresi = 1000.; Int_t id = -1;
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()>0) continue;
	TVector3 pos = cl->GetPosition();
	//if(pos.Z()>tpc::ZTarget) continue;
	if(pos.Z()>tpc::ZTarget-10.) continue;
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
      //track->Print(FUNC_NAME+" First fitting is succeeded");
      track->Print(FUNC_NAME+" First fitting is succeeded", true);
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
      else{ //More clusters are added into the track
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
		      std::vector<TPCVertexHelix*>& VertexCont,
		      Bool_t Exclusive,
		      Int_t MinNumOfHits)
{
  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);

  //Scattered helix track searching
  HighMomHelixTrackSearch(ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
#if RemainingClustersTest
  ResetHoughFlag(ClCont, BadForTracking);
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
#endif
  CalcTracks(TrackCont); //before the VertexSearch() calculation should proceed.

  //Vertex finding with tracks in the TrackCont.
  VertexSearch(TrackCont, VertexCont);
#if FragmentedTrackTest
  //Merged fragmented tracks
  RestoreFragmentedTracks(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif
  if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);
#if ReassignClusterTest
  ReassignClustersNearTheTarget(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

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
		      std::vector<TPCVertexHelix*>& VertexCont,
		      Bool_t Exclusive,
		      Int_t MinNumOfHits)
{
  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);

  XZhough_x.clear();
  XZhough_y.clear();
  XZhough_z.clear();
  Yhough_x.clear();
  Yhough_y.clear();

  //Track finding and fitting
  //for K1.8 track searching
  if(!BeamThroughTPC) K18TrackSearch(K18VPs, ClCont, TrackCont, TrackContVP); //default NimNumOfHits = 2
  //Scattered helix track searching
  HighMomHelixTrackSearch(ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
#if RemainingClustersTest
  ResetHoughFlag(ClCont, BadForTracking);
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
#endif
  CalcTracks(TrackCont); //before the VertexSearch() calculation should proceed.

  //Vertex finding with tracks in the TrackCont.
  VertexSearch(TrackCont, VertexCont);
#if FragmentedTrackTest
  //Merged fragmented tracks
  RestoreFragmentedTracks(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif
  if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);
#if ReassignClusterTest
  ReassignClustersNearTheTarget(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

  //Search Kurama track candidates from the TrackCont.
  KuramaTrackSearch(KuramaVPs, ClCont, TrackCont, TrackContFailed, TrackContVP, VertexCont, Exclusive, MinNumOfHits);

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

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
  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);
  static const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

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
      track->SetClustersHoughFlag(BadHoughTransform);
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
      track->SetClustersHoughFlag(BadHoughTransform);
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
  if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);
}

//_____________________________________________________________________________
void
HoughTransformTestHelix(const std::vector<TPCClusterContainer>& ClCont,
			std::vector<TPCLocalTrackHelix*>& TrackCont,
			Int_t MinNumOfHits /*=8*/)
{
  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);

  // HoughTransform binning
  static const auto MaxHoughWindow = gUser.GetParameter("MaxHoughWindow");
  static const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

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
      track->SetClustersHoughFlag(BadHoughTransform);
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
      track->SetClustersHoughFlag(BadHoughTransform);
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
  if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);
}

//_____________________________________________________________________________
void
HighMomHelixTrackSearch(const std::vector<TPCClusterContainer>& ClCont,
			std::vector<TPCLocalTrackHelix*>& TrackCont,
			std::vector<TPCLocalTrackHelix*>& TrackContFailed,
			Int_t MinNumOfHits)
{

  static const auto MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

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
      trackTemp->SetClustersHoughFlag(BadHoughTransform);
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
      trackTemp->SetClustersHoughFlag(BadHoughTransform);
      delete trackTemp;
      continue;
    }

    //Convert the temporary straight line track into helix track
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix;
    if(!ConvertTrack(trackTemp, track)){
      //track -> Print("converting is failed", true);
#if DebugDisp
      track->Print(FUNC_NAME+ " Track converting is failed");
#endif
      track -> SetClustersHoughFlag(BadForTracking); //track w/ few clusters
      delete track;
      continue;
    }

#if DebugDisp
    track->Print(FUNC_NAME+ " Track converting is succeeded");
#endif

    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_hough - before_hough);
    track->SetSearchTime(sec.count());

    //Track fitting processes
    FitTrack(track, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  }// tracki
  ResetHoughFlag(ClCont);
  ResetHoughFlag(ClCont, BadHoughTransform);
  ResetHoughFlag(ClCont, BadForTracking);

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

}

//_____________________________________________________________________________
void
VertexSearch(std::vector<TPCLocalTrackHelix*>& TrackCont,
	     std::vector<TPCVertexHelix*>& VertexCont)
{
  //pair
  for(Int_t trackid1=0; trackid1<TrackCont.size(); trackid1++){
    TPCLocalTrackHelix *track1 = TrackCont[trackid1];
    for(Int_t trackid2=trackid1+1; trackid2<TrackCont.size(); trackid2++){
      TPCLocalTrackHelix *track2 = TrackCont[trackid2];
      TPCVertexHelix *vertex = new TPCVertexHelix(trackid1, trackid2);
      vertex -> Calculate(track1, track2);

      Double_t closedist = vertex -> GetClosestDist();
      if(closedist < VertexDistCut){
	VertexCont.push_back(vertex);

#if DebugDisp
	vertex->Print(FUNC_NAME+" Vertex finding");
	//vertex->Print(FUNC_NAME+" Vertex finding", true);
#endif
      }
      else delete vertex;
    } //trackid2
  } //trackid1

}

//_____________________________________________________________________________
template <typename T> void
RestoreFragmentedTracks(const std::vector<TPCClusterContainer>& ClCont,
			std::vector<T*>& TrackCont,
			std::vector<T*>& TrackContFailed,
			std::vector<TPCVertexHelix*>& VertexCont,
			Bool_t Exclusive,
			Int_t MinNumOfHits)
{

  std::vector<TPCVertexHelix*> candidates_VertexCont;
  std::vector<Int_t> candidates;
  for(auto& vertex: VertexCont){
    //Threshold conditions of a distance and an angle between two tracks.
    TVector3 vtx = vertex -> GetVertex();
    Int_t trackid1 = vertex -> GetTrackId(0);
    Int_t trackid2 = vertex -> GetTrackId(1);
    if(TrackCont[trackid1] -> GetIsK18()==1 ||
       TrackCont[trackid2] -> GetIsK18()==1) continue;
    /*
      std::cout<<"id "<<vertex -> GetTrackId(0)<<" "<<vertex -> GetTrackId(1)<<std::endl;
      std::cout<<"vtx "<<TMath::Abs(vtx.x())<<" "<<TMath::Abs(vtx.y())<<" "<<TMath::Abs(vtx.z())<<std::endl;
      std::cout<<"angle "<<0.1*TMath::Pi()<<" "<<vertex -> GetOpeningAngle()<<" "<<(1. - 0.1)*TMath::Pi()<<std::endl;
    */
    //case1. accidental beam crossing the target is splitted into two tracks
    //case2. merging fragmentations of commom track. for case2. closest point of two parts are not in the track
    Bool_t is_fragmented_accidental = (TrackCont[trackid1] -> GetIsBeam()==1 || TrackCont[trackid2] -> GetIsBeam()==1);

    //two parts are close.
    if(vertex -> GetClosestDist() > 20. ||
       (vertex -> GetOpeningAngle() > 0.1*TMath::Pi() &&
	vertex -> GetOpeningAngle() < (1. - 0.1)*TMath::Pi())) continue;

    //closest point is not in the target.
    //without this, two scattered tracks with opposite direction frequently wrongly merged.
    if(TMath::Abs(vtx.x()) < 15. &&
       TMath::Abs(vtx.y()) < 10. &&
       TMath::Abs(vtx.z() - tpc::ZTarget) < 10.) continue;

#if DebugDisp
    vertex -> Print(FUNC_NAME+" Vertex candidate for merging tracks");
    //vertex -> Print(FUNC_NAME+" Vertex candidate for merging tracks", true);
#endif

    //Avoid duplication
    if(find(candidates.begin(), candidates.end(), trackid1) != candidates.end()) continue;
    if(find(candidates.begin(), candidates.end(), trackid2) != candidates.end()) continue;

    //Longer track(track1) is a reference.
    //Add two tracks.
    Bool_t order = TrackCont[trackid1] -> GetNHit() >= TrackCont[trackid2] -> GetNHit() ?  true : false;
    if(!order){
      trackid1 = vertex -> GetTrackId(0);
      trackid2 = vertex -> GetTrackId(1);
    }

    T *track1 = TrackCont[trackid1];
    T *track2 = TrackCont[trackid2];
    T *MergedTrack = new T(track1);
    for(Int_t hit=0;hit<track2 -> GetNHit();hit++)
      MergedTrack -> AddTPCHit(new TPCLTrackHit(track2 -> GetHitInOrder(hit) -> GetHit()));

    //Check whether tracks belong to the same track or not
    if(!MergedTrack -> TestMergedTrack()){
      candidates_VertexCont.push_back(vertex);
      delete MergedTrack;
      continue;
    }

    candidates.push_back(trackid1);
    candidates.push_back(trackid2);

#if DebugDisp
    MergedTrack -> Print(FUNC_NAME+" Before fitting the merged track");
    //MergedTrack -> Print(FUNC_NAME+" Before fitting the merged track", true);
#endif

    Int_t prev_size = TrackCont.size();
    FitTrack(MergedTrack, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
#if DebugDisp
    MergedTrack -> Print(FUNC_NAME+" After fitting the merged track");
    //MergedTrack -> Print(FUNC_NAME+" After fitting the merged track", true);
#endif

    Int_t post_size = TrackCont.size();
    if(prev_size+1 == post_size){
      MergedTrack = TrackCont[post_size-1];
      MergedTrack -> Calculate();
      if(Exclusive){
	MergedTrack -> DoFitExclusive();
	MergedTrack -> CalculateExclusive();
      }
      if(is_fragmented_accidental){
	MergedTrack -> SetIsAccidental();
	MergedTrack -> SetIsBeam();
      }
    }
    else{
      candidates_VertexCont.push_back(vertex);
      for(Int_t i=0; i<2; ++i){
	Int_t trackID = candidates[candidates.size() - 1];
	TrackCont[trackID] -> SetClustersHoughFlag(GoodForTracking);
	candidates.erase(candidates.begin() + candidates.size() - 1);
      }
    }
  }

  //(Iterative process) For candidates with fitting failed, trying another way for fitting.
  for(auto& vertex: candidates_VertexCont){
    TVector3 vtx = vertex -> GetVertex();
    Int_t trackid1 = vertex -> GetTrackId(0);
    Int_t trackid2 = vertex -> GetTrackId(1);
    Bool_t is_fragmented_accidental = (TrackCont[trackid1] -> GetIsBeam()==1 || TrackCont[trackid2] -> GetIsBeam()==1);

#if DebugDisp
    vertex -> Print(FUNC_NAME+" Vertex candidate for merging tracks");
    //vertex -> Print(FUNC_NAME+" Vertex candidate for merging tracks", true);
#endif

    //Avoid duplication
    if(find(candidates.begin(), candidates.end(), trackid1) != candidates.end()) continue;
    if(find(candidates.begin(), candidates.end(), trackid2) != candidates.end()) continue;

    //Longer track(track1) is a reference.
    //Add two tracks.
    Bool_t order = TrackCont[trackid1] -> GetNHit() >= TrackCont[trackid2] -> GetNHit() ?  true : false;
    if(!order){
      trackid1 = vertex -> GetTrackId(0);
      trackid2 = vertex -> GetTrackId(1);
    }

    T *track1 = TrackCont[trackid1];
    T *track2 = TrackCont[trackid2];
    T *MergedTrack = new T(track1);
    Int_t total_nhit = track1 -> GetNHit() + track2 -> GetNHit();
    std::vector<TPCClusterContainer> clusters_fortest(NumOfLayersTPC);
    for(Int_t hitid=0;hitid<track2 -> GetNHit();hitid++){
      TPCHit *hit = track2 -> GetHitInOrder(hitid) -> GetHit();
      hit -> SetHoughFlag(0);
      Int_t layer = hit -> GetLayer();
      clusters_fortest[layer].push_back(hit -> GetParentCluster());
    }

    Int_t prev_size = TrackCont.size();
    FitTrack(MergedTrack, GoodForTracking, clusters_fortest, TrackCont, TrackContFailed, MinNumOfHits);

    Int_t post_size = TrackCont.size();
    if(prev_size+1 == post_size){
      MergedTrack = TrackCont[post_size-1];

#if DebugDisp
      std::cout<<FUNC_NAME+" #of bad clusters : "<<total_nhit - MergedTrack -> GetNHit()<<std::endl;
#endif

      if(total_nhit - MergedTrack -> GetNHit() > 3){
	delete MergedTrack;
	TrackCont.erase(TrackCont.begin() + post_size - 1);
	track1 -> SetClustersHoughFlag(GoodForTracking);
	track2 -> SetClustersHoughFlag(GoodForTracking);
	continue;
      }
      else{
	MergedTrack -> Calculate();
	if(Exclusive){
	  MergedTrack -> DoFitExclusive();
	  MergedTrack -> CalculateExclusive();
	}
	if(is_fragmented_accidental){
	  MergedTrack -> SetIsAccidental();
	  MergedTrack -> SetIsBeam();
	}
      }
    }
    else{
      track1 -> SetClustersHoughFlag(GoodForTracking);
      track2 -> SetClustersHoughFlag(GoodForTracking);
      continue;
    }

    candidates.push_back(trackid1);
    candidates.push_back(trackid2);
  }
  candidates_VertexCont.clear();

  if(candidates.size()>0){
    std::sort(candidates.begin(), candidates.end());
    for(Int_t i=0; i<candidates.size(); ++i){
      Int_t trackID = candidates[i] - i;
      T *prevtrack = TrackCont[trackID];
      TrackCont.erase(TrackCont.begin() + trackID);
      delete prevtrack;
    }

    //Vertex finding again with new tracks
    del::ClearContainer(VertexCont);
    VertexSearch(TrackCont, VertexCont);
  }
}

//_____________________________________________________________________________
template <typename T> void
ReassignClustersNearTheTarget(const std::vector<TPCClusterContainer>& ClCont,
			      std::vector<T*>& TrackCont,
			      std::vector<T*>& TrackContFailed,
			      std::vector<TPCVertexHelix*>& VertexCont,
			      Bool_t Exclusive,
			      Int_t MinNumOfHits)
{

  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);
  const Int_t MostInnerlayer_Track = 6; //testing layers from 0 to "MostInnerlayer_Cluster".
  const Int_t MostInnerlayer_Cluster = 4; //testing layers from 0 to "MostInnerlayer_Cluster".

  Bool_t status = false;
  Int_t ntracks = TrackCont.size();
  Int_t k18id = -9999;
  std::vector<Int_t> candidates_trackid;
  std::vector<T*> newtracks_fortest;
  std::vector<TPCHit*> clusters_fortest;
  std::vector<Int_t> clusters_original_trackid;
  for(Int_t trackid=0; trackid<ntracks; trackid++){
    T *track = TrackCont[trackid];
    if(track -> GetIsAccidental()==1) continue;
    if(track -> GetIsK18()==1) k18id = trackid;
    if(track -> GetClosestDist() < 50.){ //tracks from the target
      Int_t n = track->GetNHit();
      TPCLTrackHit *hitp0 = track->GetHitInOrder(0); //most inner cluster
      TPCHit *hit0 = hitp0->GetHit();
      Int_t layer0 = hitp0->GetLayer();
      TVector3 pos0 = hitp0->GetLocalHitPos();
      if(layer0 > MostInnerlayer_Track) continue;
      status = true;
      candidates_trackid.push_back(trackid);

      T *CopiedTrack0 = new T(track);
      newtracks_fortest.push_back(CopiedTrack0);

      if(n <= MinNumOfHits || layer0 > MostInnerlayer_Cluster) continue;
      T *CopiedTrack1 = new T(CopiedTrack0);
      Int_t hitorder0 = CopiedTrack0 -> GetOrder(0);
      CopiedTrack1 -> EraseHit(hitorder0);

      if(CopiedTrack1->DoFit(MinNumOfHits)){
	hit0 -> SetHoughFlag(Candidate);
	clusters_fortest.push_back(hit0);
	clusters_original_trackid.push_back(trackid);

	newtracks_fortest[newtracks_fortest.size() - 1] = CopiedTrack1;
	delete CopiedTrack0;

	TPCLTrackHit *hitp1 = track->GetHitInOrder(1); //second most inner cluster
	TPCHit *hit1 = hitp1->GetHit();
	Int_t layer1 = hitp1->GetLayer();
	TVector3 pos1 = hitp1->GetLocalHitPos();
	if(n-1 <= MinNumOfHits || layer1 > MostInnerlayer_Cluster) continue;

	//Exclude most inner cluster and dofit.
	T *CopiedTrack2 = new T(CopiedTrack1);
	Int_t hitorder1 = CopiedTrack2 -> GetOrder(0);
	CopiedTrack2 -> EraseHit(hitorder1);
	if(CopiedTrack2->DoFit(MinNumOfHits)){
	  hit1 -> SetHoughFlag(Candidate);
	  clusters_fortest.push_back(hit1);
	  clusters_original_trackid.push_back(trackid);

	  newtracks_fortest[newtracks_fortest.size() - 1] = CopiedTrack2;
	  delete CopiedTrack1;
	}
	else delete CopiedTrack2;
      }
      else delete CopiedTrack1;
    }
  }
  if(!status) return;

  T *copiedK18track;
  if(k18id!=-9999){
    T *K18track = TrackCont[k18id];
    copiedK18track = new T(K18track);
    candidates_trackid.push_back(k18id);
    newtracks_fortest.push_back(copiedK18track);
  }

  //After excluding clusters near the target, try to add more clusters
  //Because of clusters near the target, frequently tracking quility becomes worse.
  //So without those clusters, it is better to try to add clusters
  for(Int_t i=0; i<newtracks_fortest.size(); i++){
    Int_t trackid = candidates_trackid[i];
    if(trackid == k18id) continue;

    T *track = newtracks_fortest[i];
    T *CopiedTrack = new T(track);

    Int_t prev_size = TrackCont.size();
    FitTrack(CopiedTrack, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
    Int_t post_size = TrackCont.size();
    if(prev_size+1 == post_size){ //Fitting is succeeded
      CopiedTrack = TrackCont[post_size-1];
      if(Exclusive) CopiedTrack -> DoFitExclusive();

      TrackCont.erase(TrackCont.begin() + TrackCont.size() - 1);

      if(CopiedTrack -> GetNHit() >= track -> GetNHit()){
	newtracks_fortest[i] = CopiedTrack;
	delete track;
      }
      else delete CopiedTrack;
    }
  }
  ResetHoughFlag(ClCont); //return Houghflag from "candidate" to "0""

#if DebugDisp
  std::cout<<FUNC_NAME+" Finding the best combination between tracks and cludsters"<<std::endl;
#endif

  Bool_t reassaign = false;
  std::vector<Int_t> reassaign_candidates;
  for(Int_t i=0; i<clusters_fortest.size(); i++){
    TPCHit *hit = clusters_fortest[i]; //cluster for testing
    TVector3 pos = hit -> GetPosition();
    Int_t best_track = -9999;
    Double_t min_residual = 9999;
    for(Int_t j=0; j<newtracks_fortest.size(); j++){
      Double_t residual = 0;
      if(newtracks_fortest[j] -> IsGoodHitToAdd(hit, residual)){
	if(min_residual > residual){
	  min_residual = residual;
	  best_track = j;
	}
      }
    }

    Int_t best_trackid = candidates_trackid[best_track];
    Int_t id = clusters_original_trackid[i]; //cluster's original id
    if(id==best_track){ //cluster is suitable for the original track
      hit -> SetHoughFlag(GoodForTracking);
    }
    else{ //cluster is not suitable for the original track
      reassaign = true;
      if(best_track==-9999){ //cluster is not suitable for all tracks
	reassaign_candidates.push_back(id);
	hit -> SetHoughFlag(0);
      }
      else{ //cluster is suitable for the other track not the original track
	reassaign_candidates.push_back(id);
	reassaign_candidates.push_back(best_trackid);
	newtracks_fortest[best_track] -> AddTPCHit(new TPCLTrackHit(hit));
      }
    }
  }
  if(newtracks_fortest.size()!=candidates_trackid.size()) std::cout<<FUNC_NAME+" FATAL Error #of tracks is not matched"<<std::endl;

  std::sort(reassaign_candidates.begin(), reassaign_candidates.end());
  reassaign_candidates.erase(std::unique(reassaign_candidates.begin(), reassaign_candidates.end()), reassaign_candidates.end());

  if(!reassaign){ //No change
    if(reassaign_candidates.size()!=0) std::cout<<FUNC_NAME+" FATAL Error #of tracks is not matched"<<std::endl;
#if DebugDisp
    std::cout<<FUNC_NAME+" No reassigning happened"<<std::endl;
#endif
    for(Int_t i=0; i<newtracks_fortest.size(); i++) delete newtracks_fortest[i];
  }
  else{ //Reassaigning happens
#if DebugDisp
    std::cout<<FUNC_NAME+" Candidates for reassigning "<<reassaign_candidates.size()<<std::endl;
#endif
    for(Int_t i=0; i<newtracks_fortest.size(); i++){
      Int_t trackid = candidates_trackid[i];
      if(std::find(reassaign_candidates.begin(), reassaign_candidates.end(), trackid)
	 == reassaign_candidates.end()) delete newtracks_fortest[i];
      else{
	Int_t threshold = MinNumOfHits;
	if(trackid == k18id) threshold = 3;

	if(newtracks_fortest[i] -> DoFit(threshold)){
	  newtracks_fortest[i] -> SetClustersHoughFlag(GoodForTracking);
	  newtracks_fortest[i] -> Calculate();
	  if(Exclusive){
	    newtracks_fortest[i] -> DoFitExclusive();
	    newtracks_fortest[i] -> CalculateExclusive();
	  }

	  //update with a new track and delete the previous track
	  T *prev_track = TrackCont[trackid];
	  TrackCont[trackid] = newtracks_fortest[i];
	  if(prev_track -> GetIsBeam()){
	    std::cout<<"isbeam trackid "<<trackid<<newtracks_fortest[i]->GetIsBeam()<<std::endl;
	  }
	  delete prev_track;
	}
	else delete newtracks_fortest[i];
      }
    } //for scattered tracks

    //Vertex finding again with new tracks
    del::ClearContainer(VertexCont);
    VertexSearch(TrackCont, VertexCont);
#if FragmentedTrackTest
    //Merged fragmented tracks
    RestoreFragmentedTracks(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif
    if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);
  } //reassaigning process
}

} //namespace tpc
