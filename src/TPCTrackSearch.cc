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

//Veto process for accidental coincidence event
MarkingAccidentalTracks : find accidental beams
FindAccidentalCoincidenceTracks &&  MarkingClusteredAccidentalTracks: find accidental cpincidence track cluster for the reacted accidental event(not beam-thorugh)

Detailed fitting procedures are explained in the TPCLocalTrack/Helix.
*/

#include "TPCTrackSearch.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "ConfMan.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "HoughTransform.hh"
#include "TPCCluster.hh"
#include "TPCLocalTrack.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCLTrackHit.hh"
#include "TPCPadHelper.hh"
#include "TPCReconstructor.hh"
#include "TPCVertex.hh"
#include "UserParamMan.hh"

#define DebugDisp 0
#define FragmentedTrackTest 1
#define ReassignClusterTest 1
//#define RemainingClustersTest 1
#define RemainingClustersTest 0 //not helpful

namespace
{
  const auto qnan = TMath::QuietNaN();
  const auto& gUser = UserParamMan::GetInstance();
  const auto& gCounter = debug::ObjectCounter::GetInstance();
  const Int_t MaxNumOfTrackTPC = 30;

  const Double_t K18XZWindow = 10.5;
  //const Double_t K18YWindow = 10.;
  const Double_t K18YWindow = 15.;

  // Min helix radius [mm] to accept linear->helix conversion (HighMomHelixTrackSearch).
  // At B = 1 T, 200 mm corresponds to pT ~ 60 MeV/c (pT = 0.3 * B[T] * R[m]).
  const Double_t MIN_HELIX_RADIUS_LINEAR_CONVERT = 200.;

  // Closest distance cut for vertex finding
  const Double_t VertexDistCut = 30.;
  const Double_t ppi_distcut = 10.; //Closest distance for p, pi at the vertex point

  // RestoreFragmentedTracks: max closest distance [mm] for merging fragmented track pairs
  const Double_t FragmentMergeClosestDistMax = 20.;

  // RestoreFragmentedAccidentalTracks: max |dz| for downstream parent-track candidates
  // (same scale as beam-like |dz| cut in TPCLocalTrackHelix.cc)
  const Double_t MaxCandidateAbsDz = 0.05;

  // Maximum number of fitting steps
  const Int_t MaxFitSteps = 10;

  // Houghflags
  const Int_t GoodForTracking = 100;
  const Int_t K18Tracks = 200;

  const Int_t BadHoughTransform = 300;
  const Int_t BadForTracking = 400;
  const Int_t Candidate = 1000;

  // Minimum #hits for Calculate()/GetdEdx() on failed tracks.
  const Int_t MIN_HITS_FOR_FAILED_CALC = 2;

  // Tracks in the Hough-Space
  std::vector<Double_t> XZhough_x;
  std::vector<Double_t> XZhough_y;
  std::vector<Double_t> XZhough_z;
  std::vector<Double_t> Yhough_x;
  std::vector<Double_t> Yhough_y;

  // B-field
  const Double_t& HS_field_0 = ConfMan::Get<Double_t>("HSFLDCALIB");
  const Double_t& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
  const Double_t& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");

  //_____________________________________________________________________________
  // Local Functions
  //_____________________________________________________________________________
  // ReassignClustersNearTheTarget: closest hit in target volume (local coords).
  template <typename T>
  Bool_t FindClosestTargetHit(const T* track, Int_t max_layer, Int_t& hit_index)
  {
    if(!track) return false;
    hit_index = -1;
    Double_t best_r_xy = 1.e10;
    for(Int_t ih = 0, n = track->GetNHit(); ih < n; ++ih){
      TPCLTrackHit* hitp = track->GetHit(ih);
      if(!hitp) continue;
      if(hitp->GetLayer() > max_layer) continue;
      const TVector3& pos = hitp->GetLocalHitPos();
      const Double_t r_xy = TMath::Hypot(pos.x(), pos.y());
      if(r_xy >= tpc::TARGET_RADIUS) continue;
      if(TMath::Abs(pos.z()) >= tpc::TARGET_HALF_Y) continue;
      if(r_xy < best_r_xy){
        best_r_xy = r_xy;
        hit_index = ih;
      }
    }
    return hit_index >= 0;
  }

  template <typename T> void
  CalcTracks(std::vector<T*>& TrackCont)
  {
    for(auto& track: TrackCont){
      track->Calculate();
    }
  }

  template <typename T> void
  DropFailedTracksBelowMinHits(std::vector<T*>& track_cont, Int_t min_hits)
  {
    for(Int_t i = static_cast<Int_t>(track_cont.size()) - 1; i >= 0; --i){
      T* track = track_cont[i];
      if(!track || track->GetNHit() < min_hits){
        delete track;
        track_cont.erase(track_cont.begin() + i);
      }
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
      if(track->GetIsAccidental()!=1) track->CheckIsAccidental();
    }
  }

  //_____________________________________________________________________________
  template <typename T> void
  MarkingBeamTracks(std::vector<T*>& TrackCont)
  {
    static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1.);
    if(BeamThroughTPC) return;

    static const Double_t default_max_abs_dz = 0.05;
    const Double_t max_abs_dz =
      gUser.Has("BeamLikeMaxAbsDzTPC")
        ? gUser.GetParameter("BeamLikeMaxAbsDzTPC")
        : default_max_abs_dz;

    for(auto& track: TrackCont){
      if(!track || track->GetIsBeam()==1 || track->GetIsK18()==1) continue;
      if(track->GetCharge() < 0 &&
         TMath::Abs(track->Getdz()) < max_abs_dz){
        track->SetIsBeam(1);
      }
    }
  }

  template <typename T> void
  MarkingClusteredAccidentalTracks(std::vector<T*>& TrackCont, std::vector<TPCVertex*>& ClusteredVertexCont)
  {
    for(auto& vertex: ClusteredVertexCont){
      vertex->SetIsAccidental();
      Int_t ntracks = vertex->GetNTracks();
      for(Int_t i=0; i<ntracks; i++){
        Int_t id = vertex->GetTrackId(i);
        TrackCont[id]->SetIsAccidental();
      }
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
       HelixTrack->Getr() < MIN_HELIX_RADIUS_LINEAR_CONVERT) status = false;
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
          Int_t vtxflag = Track->GetVtxFlag();
          TVector3 pos  = hit->GetPosition();
          Int_t side    = Track->Side(pos);
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
  if(Track->DoFit(MinNumOfHits)){
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
// FitTrack: Iteratively extend a seed Track by adding clusters and re-fitting.
//
// while(true){
//   ExtendedTrack = copy(Track)            // work on a copy; promote on success
//
//   if(nstep > MaxFitSteps)                // safety stop
//     -> finalize(Track); break
//   else if(nstep > 0 && !AddClusters)     // no good clusters left
//     -> try BadHough / BadTracking pools (last-chance);
//        if still none -> finalize(Track); break
//
//   if(FitStep(ExtendedTrack))             // fit improved -> adopt
//     -> Track = ExtendedTrack
//   else                                   // fit degraded -> stop
//     -> finalize(Track); break
//
//   nstep++
// }
//
// finalize(track):
//   IsGood(track) ? TrackCont (Houghflag) : TrackContFailed (BadForTracking)
//   - hits in TrackContFailed may be reused as seeds by other tracks.
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
    Int_t thr_ncl = 0; // hit-count check is deferred to IsGood() at the end.

    // nstep == 0 is the initial iteration: fit the seed track as-is.
    // MaxFitSteps check is harmless at nstep == 0 (always false),
    // and AddClusters is explicitly guarded by 'nstep > 0' below.
    if(nstep > MaxFitSteps){ //Tracking is over
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
    else if(nstep > 0 && !AddClusters(ExtendedTrack, ClCont)){
      // No more clusters for addding, then check for clusters of bad tracks.
      Bool_t add_badhough    = AddClusters(ExtendedTrack, ClCont, BadHoughTransform);
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
  if(Track) Track->SetFitTime(sec.count());
  if(Track) Track->SetFitFlag(nstep);

}

//_____________________________________________________________________________
Int_t
LocalTrackSearch(const std::vector<TPCClusterContainer>& ClCont,
		 std::vector<TPCLocalTrack*>& TrackCont,
		 std::vector<TPCLocalTrack*>& TrackContFailed,
		 Bool_t exclusive,
		 Int_t MinNumOfHits)
{

  // MaxHoughWindowY: max perp. dist [mm] to XZ line (2nd Hough gate) and to XZ/YZ lines (MakeLinearTrack)
  static const Double_t max_hough_window_y = gUser.GetParameter("MaxHoughWindowY");

  XZhough_x.clear();
  XZhough_y.clear();
  Yhough_x.clear();
  Yhough_y.clear();

  Bool_t prev_add = true;
  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    if(!prev_add) continue;
    prev_add = false;

#if DebugDisp
    std::cout << FUNC_NAME + " tracki : " << tracki << std::endl;
#endif

    std::chrono::milliseconds sec;
    auto before_hough = std::chrono::high_resolution_clock::now();

    //Line Hough-transform on the XZ plane
    Double_t linear_par[4]; Int_t max_bin_xz[2];
    if(!tpc::HoughTransformLineXZ(ClCont, max_bin_xz, linear_par, MinNumOfHits)){
#if DebugDisp
      std::cout<<FUNC_NAME+" No more track candiate! tracki : "<<tracki<<std::endl;
#endif
      break;
    }

    //Line Hough-transform on the YZ or YX plane
    Int_t max_bin_y[2];
    if(TMath::Abs(linear_par[2]) < 1.) tpc::HoughTransformLineYZ(ClCont, max_bin_y, linear_par, max_hough_window_y);
    else tpc::HoughTransformLineYX(ClCont, max_bin_y, linear_par, max_hough_window_y);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrack *track = new TPCLocalTrack;
    track->SetParam(linear_par);

    //If two tracks are merged at the target, separate them and recalculate params.
    Bool_t is_valid_after_sep;
    prev_add = MakeLinearTrack(track, is_valid_after_sep, ClCont, linear_par, max_hough_window_y);
    if(!prev_add) break;  // memo: replaced 'continue' with 'break' since they behave the same here.
    if(!track->IsGoodForTracking() || !is_valid_after_sep){
      track->SetClustersHoughFlag(BadHoughTransform);
      TrackContFailed.push_back(track);
      continue;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<XZhough_x.size(); ++i){
      Int_t bindiff_xz = TMath::Abs(max_bin_xz[0] - XZhough_x[i]) + TMath::Abs(max_bin_xz[1] - XZhough_y[i]);
      Int_t bindiff_y = TMath::Abs(max_bin_y[0] - Yhough_x[i]) + TMath::Abs(max_bin_y[1] - Yhough_y[i]);
      if(bindiff_xz<=1 && bindiff_y<=1){
        hough_flag = false;
#if DebugDisp
        std::cout<<"Previous hough bin on the XZ plane "<<i<<"th x: "
            <<XZhough_x[i]<<", y: "<<XZhough_y[i]<<" on the vertical plane x: "
            <<Yhough_x[i]<<", y: "<<Yhough_y[i]<<std::endl;
        std::cout<<"Current hough bin on the XZ plane "<<i<<"th x: "
            <<max_bin_xz[0]<<", y: "<<max_bin_xz[1]<<" on the vertical plane x: "
            <<max_bin_y[0]<<", y: "<<max_bin_y[1]<<std::endl;
#endif
      }
    }
    XZhough_x.push_back(max_bin_xz[0]);
    XZhough_y.push_back(max_bin_xz[1]);
    Yhough_x.push_back(max_bin_y[0]);
    Yhough_y.push_back(max_bin_y[1]);

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
  CalcTracks(TrackContFailed);
  // MarkingAccidentalTracks(TrackCont); // memo: no-op for DstTPCTracking (IsAccidental is not consumed)
  if(exclusive) ExclusiveTracking(TrackCont);
  return TrackCont.size();
}

//_____________________________________________________________________________
Bool_t
MakeLinearTrack(TPCLocalTrack *track, Bool_t &is_valid_after_sep,
                const std::vector<TPCClusterContainer>& ClCont,
                Double_t *linear_par, Double_t max_hough_window)
{
  Bool_t status = false;

  //Check Hough-distance and add hits
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(!hit) continue;
      if(hit->GetHoughFlag()>0) continue;
      TVector3 pos = cl->GetPosition();
      pos -= TVector3(0., 0., tpc::Z_TARGET);
      Double_t dist_xz = TMath::Abs(linear_par[2]*pos.Z() - pos.X() + linear_par[0])
                         / TMath::Hypot(linear_par[2], 1.);
      Double_t dist_yz = TMath::Abs(linear_par[3]*pos.Z() - pos.Y() + linear_par[1])
                         / TMath::Hypot(linear_par[3], 1.);
      if(dist_xz < max_hough_window && dist_yz < max_hough_window){
        hit->SetHoughDist(dist_xz);
        hit->SetHoughDistY(dist_yz);
        track->AddTPCHit(new TPCLTrackHit(hit));
        status = true;
      }
    } //ci
  } //layer

  if(status){
    // Vtx in the target, need to check whether two tracks are merged or not
    track->SetClustersHoughFlag(Candidate);
    is_valid_after_sep = (track->SeparateClustersWithGap() || track->SeparateTracksAtTarget());
    if(is_valid_after_sep && AddClusters(track, ClCont)) track->SetClustersHoughFlag(Candidate);
  }

#if DebugDisp
  if(status) track->Print(FUNC_NAME+" Initial track after track finding");
#endif

  if(!status) delete track;
  return status;
}

//_____________________________________________________________________________
Bool_t
MakeHelixTrack(
  TPCLocalTrackHelix *Track, Bool_t &is_valid_after_sep,
  const std::vector<TPCClusterContainer>& ClCont,
  Double_t *HelixPar, Double_t max_hough_window,
  Double_t max_hough_window_y
){

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
      Double_t tmpy = pos.z() - tpc::Z_TARGET;
      Double_t tmpz = pos.y();
      Double_t r_cal = TMath::Hypot(tmpx - HelixPar[kHelixCx], tmpy - HelixPar[kHelixCy]);
      Double_t dist = TMath::Abs(r_cal - HelixPar[kHelixR]);
      if(dist < max_hough_window){
        //for the inital track, scanning theta within range (-pi, pi)
        // Y-theta space: x-axis ~ arc length (r*theta), y-axis ~ z
        Double_t theta       = TMath::ATan2(tmpy - HelixPar[kHelixCy], tmpx - HelixPar[kHelixCx]);
        Double_t arc_length  = HelixPar[kHelixR] * theta;            // r * theta
        Double_t denom       = TMath::Hypot(HelixPar[kHelixDz], 1.);  // distance denominator
        Double_t dist_y      = TMath::Abs(HelixPar[kHelixDz]*arc_length - tmpz + HelixPar[kHelixZ0]) / denom;
        Double_t dist_y_p2pi = TMath::Abs(HelixPar[kHelixDz]*(arc_length + 2.*TMath::Pi()*HelixPar[kHelixR]) - tmpz + HelixPar[kHelixZ0]) / denom;
        Double_t dist_y_m2pi = TMath::Abs(HelixPar[kHelixDz]*(arc_length - 2.*TMath::Pi()*HelixPar[kHelixR]) - tmpz + HelixPar[kHelixZ0]) / denom;
        if(dist_y < max_hough_window_y){
          //Add hit into the track
          hit->SetHoughDist(dist);
          hit->SetHoughDistY(dist_y);
          Track->AddTPCHit(new TPCLTrackHit(hit));

          id++;
          status = true;
        } //distY
        else if(dist_y_p2pi < max_hough_window_y){
          //Add hit into the track
          hit->SetHoughDist(dist);
          hit->SetHoughDistY(dist_y_p2pi);
          Track->AddTPCHit(new TPCLTrackHit(hit));

          id++;
          status = true;
        }
        else if(dist_y_m2pi < max_hough_window_y){
          //Add hit into the track
          hit->SetHoughDist(dist);
          hit->SetHoughDistY(dist_y_m2pi);
          Track->AddTPCHit(new TPCLTrackHit(hit));

          id++;
          status = true;
        }
      } //dist
    } //ci
  } //layer

  if(status){
    Track->SetClustersHoughFlag(Candidate);
    Track->CalcHelixTheta();

    //If track has very large gap, then spilt it.
    if(Track->SeparateClustersWithGap()) Track->CalcHelixTheta();

    // Check for merged tracks at target; is_valid_after_sep = separation + prefit ok.
    is_valid_after_sep = Track->SeparateTracksAtTarget();
    if(is_valid_after_sep && AddClusters(Track, ClCont)){
      Double_t par[5];
      Track->GetParam(par);
      is_valid_after_sep = Track->DoPreFit(par);
      Track->SetClustersHoughFlag(Candidate);
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
HelixTrackSearch(
  Int_t Trackflag, Int_t Houghflag,
  const std::vector<TPCClusterContainer>& ClCont,
  std::vector<TPCLocalTrackHelix*>& TrackCont,
  std::vector<TPCLocalTrackHelix*>& TrackContFailed,
  Int_t MinNumOfHits
){

  // MaxHoughWindow: max |r_hit - R_helix| [mm] for Y-theta Hough gate and helix cluster assignment
  static const Double_t max_hough_window  = gUser.GetParameter("MaxHoughWindow");
  // MaxHoughWindowY: max perp. dist [mm] from helix to Y line (MakeHelixTrack distY cut)
  static const Double_t max_hough_window_y = gUser.GetParameter("MaxHoughWindowY");

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
    tpc::HoughTransformLineYTheta(ClCont, MaxBinY, HelixPar, max_hough_window);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    track->SetParam(HelixPar);
    track->SetFlag(Trackflag);

    //If two tracks are merged at the target, separate them and recalculate params.
    Bool_t is_valid_after_sep;
    prev_add = MakeHelixTrack(track, is_valid_after_sep, ClCont, HelixPar, max_hough_window, max_hough_window_y);
    if(!prev_add) continue;
    if(!track->IsGoodForTracking() || !is_valid_after_sep){
      track->SetClustersHoughFlag(BadHoughTransform);
      TrackContFailed.push_back(track);
      continue;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<XZhough_x.size(); ++i){
      Int_t bindiffXZ = TMath::Abs(MaxBinXZ[0] - XZhough_x[i]) + TMath::Abs(MaxBinXZ[1] - XZhough_y[i]) + TMath::Abs(MaxBinXZ[2] - XZhough_z[i]);
      Int_t bindiffY  = TMath::Abs(MaxBinY[0]  - Yhough_x[i])  + TMath::Abs(MaxBinY[1]  - Yhough_y[i]);
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
    std::cout<<FUNC_NAME+" K18 VP Helix p : "<<trackref->Getr()*tpc::C_LIGHT<<std::endl;
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
	//if(pos.Z()>tpc::Z_TARGET) continue;
	if(pos.Z()>tpc::Z_TARGET-10.) continue;
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

    Int_t Trackflag = 1*1 + 2*1 + 4*0 + 8*0; // isBeam, isK18, isKurama(legacy), isAccidental
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
LocalTrackSearchHelix(
  const std::vector<TPCClusterContainer>& ClCont,
  std::vector<TPCLocalTrackHelix*>& TrackCont,
  std::vector<TPCLocalTrackHelix*>& TrackContInvertedCharge,
  std::vector<TPCLocalTrackHelix*>& TrackContFailed,
  std::vector<TPCVertex*>& VertexCont,
  std::vector<TPCVertex*>& ClusteredVertexCont,
  Bool_t Exclusive,
  Int_t MinNumOfHits
){
  // BeamThroughTPC==1: skip accidental-track marking (beam-through mode)
  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1.);

  //Scattered helix track searching
  HighMomHelixTrackSearch(ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);

#if RemainingClustersTest
  ResetHoughFlag(ClCont, BadForTracking);
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
#endif
  CalcTracks(TrackCont); //before the VertexSearch() calculation should proceed.

  if(!BeamThroughTPC) MarkingBeamTracks(TrackCont);
  if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);

  //Vertex finding with tracks in the TrackCont.
  VertexSearch(TrackCont, VertexCont);

#if FragmentedTrackTest
  //Merged fragmented tracks
  RestoreFragmentedTracks(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

#if ReassignClusterTest
  ReassignClustersNearTheTarget(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

  FindAccidentalCoincidenceTracks(TrackCont, VertexCont, ClusteredVertexCont);

#if ReassignClusterTest
  ReassignClustersVertex(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

  TestingCharge(TrackCont, TrackContInvertedCharge, VertexCont, Exclusive);

  RestoreFragmentedAccidentalTracks(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

  DropFailedTracksBelowMinHits(TrackContFailed, MIN_HITS_FOR_FAILED_CALC);
  CalcTracks(TrackContFailed);
  if(Exclusive) ExclusiveTracking(TrackCont);
  return TrackCont.size();
}

//_____________________________________________________________________________
Int_t
LocalTrackSearchHelix(std::vector<std::vector<TVector3>> K18BRVPs,
		      const std::vector<TPCClusterContainer>& ClCont,
		      std::vector<TPCLocalTrackHelix*>& TrackCont,
		      std::vector<TPCLocalTrackHelix*>& TrackContInvertedCharge,
		      std::vector<TPCLocalTrackHelix*>& TrackContVP,
		      std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		      std::vector<TPCVertex*>& VertexCont,
		      std::vector<TPCVertex*>& ClusteredVertexCont,
		      Bool_t Exclusive,
		      Int_t MinNumOfHits)
{
  // BeamThroughTPC==1: skip K18 track search and accidental-track marking
  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1.);

  XZhough_x.clear();
  XZhough_y.clear();
  XZhough_z.clear();
  Yhough_x.clear();
  Yhough_y.clear();

  //Track finding and fitting
  //for K1.8 track searching
  if(!BeamThroughTPC) K18TrackSearch(K18BRVPs, ClCont, TrackCont, TrackContVP); //default NimNumOfHits = 2
  //Scattered helix track searching
  HighMomHelixTrackSearch(ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);

#if RemainingClustersTest
  ResetHoughFlag(ClCont, BadForTracking);
  HelixTrackSearch(0, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
#endif
  CalcTracks(TrackCont); //before the VertexSearch() calculation should proceed.

  if(!BeamThroughTPC) MarkingBeamTracks(TrackCont);
  if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);

  //Vertex finding with tracks in the TrackCont.
  VertexSearch(TrackCont, VertexCont);

#if FragmentedTrackTest
  //Merged fragmented tracks
  RestoreFragmentedTracks(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

#if ReassignClusterTest
  ReassignClustersNearTheTarget(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

  RestoreFragmentedAccidentalTracks(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);

  FindAccidentalCoincidenceTracks(TrackCont, VertexCont, ClusteredVertexCont);

#if ReassignClusterTest
  ReassignClustersVertex(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif

  TestingCharge(TrackCont, TrackContInvertedCharge, VertexCont, Exclusive);

  CalcTracks(TrackContVP);
  DropFailedTracksBelowMinHits(TrackContFailed, MIN_HITS_FOR_FAILED_CALC);
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
  // static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);
  // MaxHoughWindowY: max perp. dist [mm] to XZ line (2nd Hough gate) and to XZ/YZ lines (MakeLinearTrack)
  static const Double_t max_hough_window_y = gUser.GetParameter("MaxHoughWindowY");

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
    Double_t linear_par[4]; Int_t max_bin_xz[2];
    if(!tpc::HoughTransformLineXZ(ClCont, max_bin_xz, linear_par, MinNumOfHits)){
#if DebugDisp
      std::cout<<FUNC_NAME+" No more track candiate! tracki : "<<tracki<<std::endl;
#endif
      break;
    }

    //Line Hough-transform on the YZ or YX plane
    Int_t max_bin_y[2];
    if(TMath::Abs(linear_par[2]) < 1.) tpc::HoughTransformLineYZ(ClCont, max_bin_y, linear_par, max_hough_window_y);
    else tpc::HoughTransformLineYX(ClCont, max_bin_y, linear_par, max_hough_window_y);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrack *track = new TPCLocalTrack;
    track->SetParam(linear_par);

    //If two tracks are merged at the target, separate them and recalculate params.
    Bool_t is_valid_after_sep;
    prev_add = MakeLinearTrack(track, is_valid_after_sep, ClCont, linear_par, max_hough_window_y);
    if(!prev_add) break;  // memo: replaced 'continue' with 'break' since they behave the same here.
    if(!track -> IsGoodForTracking() || !is_valid_after_sep){
      track->SetClustersHoughFlag(BadHoughTransform);
      delete track;
      continue;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<XZhough_x.size(); ++i){
      Int_t bindiff_xz = TMath::Abs(max_bin_xz[0] - XZhough_x[i]) + TMath::Abs(max_bin_xz[1] - XZhough_y[i]);
      Int_t bindiff_y = TMath::Abs(max_bin_y[0] - Yhough_x[i]) + TMath::Abs(max_bin_y[1] - Yhough_y[i]);
      if(bindiff_xz<=1 && bindiff_y<=1){
	hough_flag = false;
#if DebugDisp
	std::cout<<"Previous hough bin on the XZ plane "<<i<<"th x: "
		 <<XZhough_x[i]<<", y: "<<XZhough_y[i]<<" on the vertical plane x: "
		 <<Yhough_x[i]<<", y: "<<Yhough_y[i]<<std::endl;
	std::cout<<"Current hough bin on the XZ plane "<<i<<"th x: "
		 <<max_bin_xz[0]<<", y: "<<max_bin_xz[1]<<" on the vertical plane x: "
		 <<max_bin_y[0]<<", y: "<<max_bin_y[1]<<std::endl;
#endif
      }
    }
    XZhough_x.push_back(max_bin_xz[0]);
    XZhough_y.push_back(max_bin_xz[1]);
    Yhough_x.push_back(max_bin_y[0]);
    Yhough_y.push_back(max_bin_y[1]);

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
  //if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);
}

//_____________________________________________________________________________
void
HoughTransformTestHelix(const std::vector<TPCClusterContainer>& ClCont,
			std::vector<TPCLocalTrackHelix*>& TrackCont,
			Int_t MinNumOfHits /*=8*/)
{
  // static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1);

  // MaxHoughWindow: max |r_hit - R_helix| [mm] for Y-theta Hough gate and helix cluster assignment
  static const Double_t max_hough_window = gUser.GetParameter("MaxHoughWindow");
  // MaxHoughWindowY: max perp. dist [mm] from helix to Y line (MakeHelixTrack distY cut)
  static const Double_t max_hough_window_y = gUser.GetParameter("MaxHoughWindowY");

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
    tpc::HoughTransformLineYTheta(ClCont, MaxBinY, HelixPar, max_hough_window);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    track->SetParam(HelixPar);

    //If two tracks are merged at the target, separate them and recalculate params.
    Bool_t is_valid_after_sep;
    prev_add = MakeHelixTrack(track, is_valid_after_sep, ClCont, HelixPar, max_hough_window, max_hough_window_y);
    if(!prev_add) continue;
    if(!track -> IsGoodForTracking() || !is_valid_after_sep){
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
  //if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);
}

//_____________________________________________________________________________
void
HighMomHelixTrackSearch(
  const std::vector<TPCClusterContainer>& ClCont,
  std::vector<TPCLocalTrackHelix*>& TrackCont,
  std::vector<TPCLocalTrackHelix*>& TrackContFailed,
  Int_t MinNumOfHits
){
  // MaxHoughWindowY: max perp. dist [mm] to XZ line (2nd Hough gate) and to XZ/YZ lines (MakeLinearTrack)
  static const Double_t max_hough_window_y = gUser.GetParameter("MaxHoughWindowY");

  std::vector<Double_t> temp_xz_hough_x;
  std::vector<Double_t> temp_xz_hough_y;
  std::vector<Double_t> temp_y_hough_x;
  std::vector<Double_t> temp_y_hough_y;

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
    Double_t linear_par[4]; Int_t max_bin_xz[2];
    if(!tpc::HoughTransformLineXZ(ClCont, max_bin_xz, linear_par, MinNumOfHits)){
#if DebugDisp
      std::cout<<FUNC_NAME+" No more track candiate! tracki : "<<tracki<<std::endl;
#endif
      break;
    }

    //Line Hough-transform on the YZ or YX plane
    Int_t max_bin_y[2];
    if(TMath::Abs(linear_par[2]) < 1.) tpc::HoughTransformLineYZ(ClCont, max_bin_y, linear_par, max_hough_window_y);
    else tpc::HoughTransformLineYX(ClCont, max_bin_y, linear_par, max_hough_window_y);

    //Make a track(HoughDistCheck)
    //The origin at the target center
    TPCLocalTrack *track_temp = new TPCLocalTrack;
    track_temp->SetParam(linear_par);

    Bool_t is_valid_after_sep;
    prev_add = MakeLinearTrack(track_temp, is_valid_after_sep, ClCont, linear_par, max_hough_window_y);
    if(!prev_add) break;  // memo: replaced 'continue' with 'break' since they behave the same here.
    if(!track_temp -> IsGoodForTracking() || !is_valid_after_sep){
      track_temp->SetClustersHoughFlag(BadHoughTransform);
      delete track_temp;
      continue;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<temp_xz_hough_x.size(); ++i){
      Int_t bindiffXZ = TMath::Abs(max_bin_xz[0] - temp_xz_hough_x[i]) + TMath::Abs(max_bin_xz[1] - temp_xz_hough_y[i]);
      Int_t bindiffY  = TMath::Abs(max_bin_y[0]  - temp_y_hough_x[i])  + TMath::Abs(max_bin_y[1]  - temp_y_hough_y[i]);
      if(bindiffXZ<=1 && bindiffY<=1){
        hough_flag = false;
#if DebugDisp
        std::cout <<"Previous hough bin on the XZ plane "<<i<<"th x: "
                  <<temp_xz_hough_x[i]<<", y: "<<temp_xz_hough_y[i]<<" on the vertical plane x: "
                  <<temp_y_hough_x[i]<<", y: "<<temp_y_hough_y[i]<<std::endl;
        std::cout <<"Current hough bin on the XZ plane "<<i<<"th x: "
                  <<max_bin_xz[0]<<", y: "<<max_bin_xz[1]<<" on the vertical plane x: "
                  <<max_bin_y[0]<<", y: "<<max_bin_y[1]<<std::endl;
#endif
      }
    }
    temp_xz_hough_x.push_back(max_bin_xz[0]);
    temp_xz_hough_y.push_back(max_bin_xz[1]);
    temp_y_hough_x.push_back(max_bin_y[0]);
    temp_y_hough_y.push_back(max_bin_y[1]);

    if(!hough_flag){
#if DebugDisp
      std::cout<<FUNC_NAME+" The same track is found by Hough-Transform : tracki : "<<tracki<<" hough cont size : "<<temp_xz_hough_x.size()<<std::endl;
#endif
      track_temp->SetClustersHoughFlag(BadHoughTransform);
      delete track_temp;
      continue;
    }

    //Convert the temporary straight line track into helix track
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix;
    if(!ConvertTrack(track_temp, track)){
#if DebugDisp
      track->Print(FUNC_NAME+ " Track converting is failed");
#endif
      track->SetClustersHoughFlag(BadForTracking); //track w/ few clusters
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
template <typename T>
void
VertexSearch(std::vector<T*>& TrackCont,
	           std::vector<TPCVertex*>& VertexCont)
{
  //pair
  for(Int_t trackid1=0; trackid1<TrackCont.size(); trackid1++){
    T *track1 = TrackCont[trackid1];
    for(Int_t trackid2=trackid1+1; trackid2<TrackCont.size(); trackid2++){
      T *track2 = TrackCont[trackid2];
      TPCVertex *vertex = new TPCVertex(trackid1, trackid2);
      vertex->Calculate(track1, track2);

      Double_t closest_dist = vertex->GetClosestDist();
      if(closest_dist < VertexDistCut){
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
			std::vector<TPCVertex*>& VertexCont,
			Bool_t Exclusive,
			Int_t MinNumOfHits)
{

  std::vector<TPCVertex*> candidates_VertexCont;
  std::vector<Int_t> candidates;
  for(auto& vertex: VertexCont){
    //Threshold conditions of a distance and an angle between two tracks.
    TVector3 vtx = vertex->GetVertex();
    Int_t trackid1 = vertex->GetTrackId(0);
    Int_t trackid2 = vertex->GetTrackId(1);
    if(TrackCont[trackid1]->GetIsK18()==1 ||
       TrackCont[trackid2]->GetIsK18()==1) continue;

    Bool_t isbeam1 = true;
    for(Int_t ihit=0; ihit<TrackCont[trackid1]->GetNHit(); ihit++){
      TPCHit *hit  = TrackCont[trackid1]->GetHitInOrder(ihit)->GetHit();
      TVector3 pos = hit->GetPosition();
      if(!TrackCont[trackid1]->IsBeamLikeHit(pos)){
        isbeam1 = false;
        break;
      }
    }
    if(isbeam1) continue; //veto beam particle

    Bool_t isbeam2 = true;
    for(Int_t ihit=0; ihit<TrackCont[trackid2]->GetNHit(); ihit++){
      TPCHit *hit  = TrackCont[trackid2]->GetHitInOrder(ihit)->GetHit();
      TVector3 pos = hit->GetPosition();
      if(!TrackCont[trackid2]->IsBeamLikeHit(pos)){
        isbeam2 = false;
        break;
      }
    }
    if(isbeam2) continue; //veto beam particle


    //case1. accidental beam crossing the target is splitted into two tracks
    //case2. merging fragmentations of commom track. for case2. closest point of two parts are not in the track

    // two parts are close.
    // Exclude mid-angle pairs (0.1*pi < opening angle < 0.9*pi):
    // these are unlikely to be fragmented pieces of one track.
    const Double_t angle_window = 0.4*TMath::Pi();
    if(vertex->GetClosestDist() > FragmentMergeClosestDistMax ||
       TMath::Abs(vertex->GetOpeningAngle() - 0.5*TMath::Pi()) < angle_window) continue;

    // closest point is not in the target.
    // without this, two scattered tracks with opposite direction frequently wrongly merged.
    if(TMath::Hypot(vtx.x(), vtx.z() - tpc::Z_TARGET) < tpc::TARGET_RADIUS &&
       TMath::Abs(vtx.y()) < tpc::TARGET_HALF_Y) continue;

#if DebugDisp
    vertex -> Print(FUNC_NAME+" Vertex candidate for merging tracks");
    //vertex -> Print(FUNC_NAME+" Vertex candidate for merging tracks", true);
#endif

    //Avoid duplication
    if(std::find(candidates.begin(), candidates.end(), trackid1) != candidates.end()) continue;
    if(std::find(candidates.begin(), candidates.end(), trackid2) != candidates.end()) continue;

    //Longer track(track1) is a reference.
    //Add two tracks.
    Bool_t order = TrackCont[trackid1]->GetNHit() >= TrackCont[trackid2]->GetNHit() ?  true : false;
    if(!order){
      trackid1 = vertex->GetTrackId(1);
      trackid2 = vertex->GetTrackId(0);
    }

    T *track1 = TrackCont[trackid1];
    T *track2 = TrackCont[trackid2];
    T *MergedTrack = new T(track1);
    for(Int_t ihit=0; ihit<track2->GetNHit(); ++ihit)
      MergedTrack->AddTPCHit(new TPCLTrackHit(track2->GetHitInOrder(ihit)->GetHit()));

    //Check whether tracks belong to the same track or not
    if(!MergedTrack->TestMergedTrack()){
      candidates_VertexCont.push_back(vertex);
      delete MergedTrack;
      continue;
    }
    candidates.push_back(trackid1);
    candidates.push_back(trackid2);

#if DebugDisp
    MergedTrack->Print(FUNC_NAME+" Before fitting the merged track");
    //MergedTrack->Print(FUNC_NAME+" Before fitting the merged track", true);
#endif

    Int_t prev_size = TrackCont.size();
    FitTrack(MergedTrack, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);

    Int_t post_size = TrackCont.size();
    if(prev_size+1 == post_size){ // A merged-track fit succeeded and one new track was added.
      MergedTrack = TrackCont[post_size-1];
      MergedTrack->Calculate();
      if(Exclusive){
        MergedTrack->DoFitExclusive();
        MergedTrack->CalculateExclusive();
      }
      MergedTrack->CheckIsAccidental();

#if DebugDisp
      MergedTrack->Print(FUNC_NAME+" After fitting the merged track");
      //MergedTrack->Print(FUNC_NAME+" After fitting the merged track", true);
#endif
    }
    else{ // A merged-track fit failed (or no new track was added); keep this vertex as a retry candidate.
      candidates_VertexCont.push_back(vertex);
      TrackCont[trackid1]->SetClustersHoughFlag(GoodForTracking);
      TrackCont[trackid2]->SetClustersHoughFlag(GoodForTracking);
      candidates.pop_back();
      candidates.pop_back();
    }
  } // end of for(auto& vertex: VertexCont)

  //(Iterative process) For candidates with fitting failed, trying another way for fitting.
  for(auto& vertex: candidates_VertexCont){
    TVector3 vtx = vertex->GetVertex();
    Int_t trackid1 = vertex->GetTrackId(0);
    Int_t trackid2 = vertex->GetTrackId(1);

#if DebugDisp
    vertex->Print(FUNC_NAME+" Vertex candidate for merging tracks");
    //vertex->Print(FUNC_NAME+" Vertex candidate for merging tracks", true);
#endif

    //Avoid duplication
    if(std::find(candidates.begin(), candidates.end(), trackid1) != candidates.end()) continue;
    if(std::find(candidates.begin(), candidates.end(), trackid2) != candidates.end()) continue;

    //Longer track(track1) is a reference.
    //Add two tracks.
    Bool_t order = TrackCont[trackid1]->GetNHit() >= TrackCont[trackid2]->GetNHit() ?  true : false;
    if(!order){
      trackid1 = vertex->GetTrackId(1);
      trackid2 = vertex->GetTrackId(0);
    }

    T *track1 = TrackCont[trackid1];
    T *track2 = TrackCont[trackid2];
    T *MergedTrack = new T(track1);
    Int_t total_nhit = track1->GetNHit() + track2->GetNHit();
    std::vector<TPCClusterContainer> clusters_fortest(NumOfLayersTPC);
    for(Int_t hitid=0; hitid<track2->GetNHit(); hitid++){
      TPCHit *hit = track2->GetHitInOrder(hitid)->GetHit();
      hit->SetHoughFlag(0);
      Int_t layer = hit->GetLayer();
      clusters_fortest[layer].push_back(hit->GetParentCluster());
    }

    Int_t prev_size = TrackCont.size();
    FitTrack(MergedTrack, GoodForTracking, clusters_fortest, TrackCont, TrackContFailed, MinNumOfHits);

    Int_t post_size = TrackCont.size();
    if(prev_size+1 == post_size){ // A merged-track fit succeeded and one new track was added.
      MergedTrack = TrackCont[post_size-1];

#if DebugDisp
      std::cout<<FUNC_NAME+" #of bad clusters : "<<total_nhit - MergedTrack -> GetNHit()<<std::endl;
#endif

      if(total_nhit - MergedTrack->GetNHit() > 3){ // might be able to optimize the threshold
        delete MergedTrack;
        TrackCont.pop_back();
        track1->SetClustersHoughFlag(GoodForTracking);
        track2->SetClustersHoughFlag(GoodForTracking);
        continue;
      }
      else{
        MergedTrack->Calculate();
        if(Exclusive){
          MergedTrack->DoFitExclusive();
          MergedTrack->CalculateExclusive();
        }
        MergedTrack->CheckIsAccidental();
      }
    }
    else{ // A merged-track fit failed (or no new track was added); keep this vertex as a retry candidate.
      track1->SetClustersHoughFlag(GoodForTracking);
      track2->SetClustersHoughFlag(GoodForTracking);
      continue;
    }

    candidates.push_back(trackid1);
    candidates.push_back(trackid2);
  } // end of for(auto& vertex: candidates_VertexCont)

  if(candidates.size()>0){
    // Sort track indices in descending order so erasing does not shift remaining indices.
    std::sort(candidates.begin(), candidates.end(), std::greater<Int_t>());
    for(Int_t track_id : candidates){
      T *prevtrack = TrackCont[track_id];
      TrackCont.erase(TrackCont.begin() + track_id);
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
			      std::vector<TPCVertex*>& VertexCont,
			      Bool_t Exclusive,
			      Int_t MinNumOfHits)
{

  // BeamThroughTPC==1: skip accidental-track marking after cluster reassignment
  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1.);
  const Int_t MaxTargetLayer = tpc::LAST_TGT_LAYER + 1;
  static const Int_t MaxStrip = 2;

  Bool_t status = false;
  Int_t ntracks = TrackCont.size();
  Int_t k18id = -9999;
  std::vector<Int_t> candidates_trackid;
  std::vector<T*> newtracks_fortest;
  std::vector<TPCHit*> clusters_fortest;
  std::vector<Int_t> clusters_original_trackid;
  for(Int_t trackid=0; trackid<ntracks; trackid++){
    T *track = TrackCont[trackid];
    if(track->GetIsAccidental()==1) continue;
    if(track->GetIsK18()==1) k18id = trackid;
    if(track->GetClosestDist() >= tpc::TARGET_RADIUS) continue;

    Int_t ih = -1;
    if(!FindClosestTargetHit(track, MaxTargetLayer, ih)) continue;
    TPCLTrackHit *hitp = track->GetHit(ih);
    if(!hitp) continue;

    status = true;
    candidates_trackid.push_back(trackid);

    T *trial = new T(track);
    newtracks_fortest.push_back(trial);

    for(Int_t iter=0; iter<MaxStrip; ++iter){
      if(trial->GetNHit() <= MinNumOfHits) break;

      if(iter > 0){
        if(!FindClosestTargetHit(trial, MaxTargetLayer, ih)) break;
        hitp = trial->GetHit(ih);
        if(!hitp) break;
      }

      TPCHit *cl = hitp->GetHit();
      T *strip = new T(trial);
      strip->EraseHit(ih);
      cl->SetHoughFlag(GoodForTracking);
      if(!strip->DoFit(MinNumOfHits)){
        delete strip;
        break;
      }

      clusters_fortest.push_back(cl);
      clusters_original_trackid.push_back(trackid);

      delete trial;
      trial = strip;
      newtracks_fortest.back() = trial;
    }
  }
  if(!status) return;

  if(k18id != -9999){
    T *k18_track = new T(TrackCont[k18id]);
    candidates_trackid.push_back(k18id);
    newtracks_fortest.push_back(k18_track);
  }

  //After excluding clusters near the target, try to add more clusters
  //Because of clusters near the target, frequently tracking quility becomes worse.
  //So without those clusters, it is better to try to add clusters
  Bool_t reassign = false;
  std::vector<Int_t> reassign_candidates;

  for(Int_t i=0; i<(Int_t)newtracks_fortest.size(); i++){
    Int_t trackid = candidates_trackid[i];
    if(trackid == k18id) continue;

    T *trial = newtracks_fortest[i];
    T *extended = new T(trial);

    Int_t n_tracks_before = TrackCont.size();
    FitTrack(extended, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
    Int_t n_tracks_after = TrackCont.size();
    if(n_tracks_before + 1 == n_tracks_after){ //Fitting is succeeded
      extended = TrackCont[n_tracks_after - 1];
      if(Exclusive){
        extended->DoFitExclusive();
        extended->CalculateExclusive();
      }
      TrackCont.pop_back();

      if(extended->GetNHit() > trial->GetNHit()){
        reassign = true;
        reassign_candidates.push_back(trackid);
        newtracks_fortest[i] = extended;
        delete trial;
      }
      else delete extended;
    }
  }

#if DebugDisp
  std::cout<<FUNC_NAME+" Finding the best combination between tracks and clusters"<<std::endl;
#endif

  // Assign each cluster stripped near the target to the best-matching candidate track.
  // Compare residuals among newtracks_fortest; keep on the original track, reassign,
  // or reject the cluster if no candidate accepts it.
  for(Int_t i=0; i<(Int_t)clusters_fortest.size(); i++){
    TPCHit *cl = clusters_fortest[i];
    Int_t orig_trackid = clusters_original_trackid[i];

    Int_t best_idx = -9999;
    Int_t best_trackid = -1;
    Double_t min_residual = 9999.;
    for(Int_t cand_idx=0; cand_idx<(Int_t)newtracks_fortest.size(); cand_idx++){
      Double_t residual = 0.;
      if(newtracks_fortest[cand_idx]->IsGoodHitToAdd(cl, residual)){
        if(min_residual > residual){
          min_residual = residual;
          best_idx = cand_idx;
          best_trackid = candidates_trackid[cand_idx];
        }
      }
    }

    cl->SetHoughFlag(GoodForTracking);
    if(orig_trackid == best_trackid){ // cluster fits the track it was stripped from
      newtracks_fortest[best_idx]->AddTPCHit(new TPCLTrackHit(cl));
    }
    else{ // cluster belongs elsewhere or nowhere
      reassign = true;
      if(best_trackid == -1){ // no candidate accepts this cluster
        reassign_candidates.push_back(orig_trackid);
        cl->SetHoughFlag(0);
      }
      else{ // cluster fits another candidate track
        reassign_candidates.push_back(orig_trackid);
        reassign_candidates.push_back(best_trackid);
        newtracks_fortest[best_idx]->AddTPCHit(new TPCLTrackHit(cl));
      }
    }
  }
  if(newtracks_fortest.size()!=candidates_trackid.size()) 
    std::cout<<FUNC_NAME+" FATAL Error #of tracks is not matched"<<std::endl;

  std::sort(reassign_candidates.begin(), reassign_candidates.end());
  reassign_candidates.erase(std::unique(reassign_candidates.begin(), reassign_candidates.end()), reassign_candidates.end());

  if(!reassign){ //No change
    if(reassign_candidates.size()!=0) std::cout<<FUNC_NAME+" FATAL Error #of tracks is not matched"<<std::endl;
#if DebugDisp
    std::cout<<FUNC_NAME+" No reassigning happened"<<std::endl;
#endif
    del::ClearContainer(newtracks_fortest);
  }
  else{ //Reassigning happens
#if DebugDisp
    std::cout<<FUNC_NAME+" Candidates for reassigning "<<reassign_candidates.size()<<std::endl;
#endif

    // Re-fit and commit tracks in reassign_candidates into TrackCont.
    for(Int_t cand_idx=0; cand_idx<(Int_t)newtracks_fortest.size(); cand_idx++){
      Int_t trackid = candidates_trackid[cand_idx];
      T *trial = newtracks_fortest[cand_idx];

      if(std::find(reassign_candidates.begin(), reassign_candidates.end(), trackid) == reassign_candidates.end()){
        delete trial;
        continue;
      }

      Int_t threshold = MinNumOfHits;
      if(trackid == k18id) threshold = 3;

      if(!trial->DoFit(threshold)){
        delete trial;
        continue;
      }

      trial->SetClustersHoughFlag(GoodForTracking);
      trial->Calculate();
      if(Exclusive){
        trial->DoFitExclusive();
        trial->CalculateExclusive();
      }

      T *prev_track = TrackCont[trackid];
      TrackCont[trackid] = trial;
      delete prev_track;
    }
  }
  ResetHoughFlag(ClCont);

  if(reassign){
    //Vertex finding again with new tracks
    del::ClearContainer(VertexCont);
    VertexSearch(TrackCont, VertexCont);
#if FragmentedTrackTest
    //Merged fragmented tracks
    RestoreFragmentedTracks(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif
    if(!BeamThroughTPC) MarkingBeamTracks(TrackCont);
    if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);
  } //reassigning process
}

//_____________________________________________________________________________
template <typename T> void
ReassignClustersVertex(const std::vector<TPCClusterContainer>& ClCont,
		       std::vector<T*>& TrackCont,
		       std::vector<T*>& TrackContFailed,
		       std::vector<TPCVertex*>& VertexCont,
		       Bool_t Exclusive,
		       Int_t MinNumOfHits)
{

  // BeamThroughTPC==1: skip accidental-track marking after vertex reassignment
  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1.);

  // strip -> assign -> commit -> VertexSearch; one reassignment per outer while pass
  const Double_t vertex_closest_dist_max = 10.;              // max GetClosestDist [mm]
  const Double_t cluster_dist_max = 20.;                     // strip window around GetTrackPos [mm]
  const Double_t arc_margin_mm = 0.5 * tpc::TARGET_RADIUS;   // Mint/Maxt margin for vertex GetTrackTheta [mm]

  std::vector<std::pair<Int_t, Int_t>> candidate_pairs;
  Bool_t reassign_pending = true;  // retry outer loop after each committed reassignment
  while(reassign_pending){
    reassign_pending = false;

    for(auto& vertex : VertexCont){
      const Int_t trackid0 = vertex->GetTrackId(0);
      const Int_t trackid1 = vertex->GetTrackId(1);

      if(std::find(candidate_pairs.begin(), candidate_pairs.end(),
		   std::make_pair(trackid0, trackid1)) != candidate_pairs.end()) continue;

      T* tracks[2] = {TrackCont[trackid0], TrackCont[trackid1]};
      const Int_t trackids[2] = {trackid0, trackid1};

      if(vertex->GetClosestDist() > vertex_closest_dist_max) continue;
      if(tracks[0]->GetIsAccidental() == 1 || tracks[1]->GetIsAccidental() == 1) continue;
      if(tracks[0]->GetIsK18() == 1 || tracks[1]->GetIsK18() == 1) continue;

      // GetTrackTheta: helix parameter at closest approach on each track (not cluster theta)
      Bool_t skip_vertex = false;
      for(Int_t itrack = 0; itrack < 2; ++itrack){
        const Double_t theta = vertex->GetTrackTheta(itrack);
        const Double_t margin_rad = arc_margin_mm / tracks[itrack]->Getr();
        if(theta < tracks[itrack]->GetMint() - margin_rad ||
           theta > tracks[itrack]->GetMaxt() + margin_rad){
          skip_vertex = true;
          break;
        }
      }
      if(skip_vertex) continue;

      // Collect hit orders within cluster_dist_max of GetTrackPos on each track
      Int_t closest_hit_idx[2] = {0, 0};
      std::vector<Int_t> close_hit_orders[2];
      for(Int_t itrack = 0; itrack < 2; ++itrack){
        const TVector3 vtx_pos = vertex->GetTrackPos(itrack);
        Double_t min_dist = 1.e10;
        for(Int_t ih = 0, n = tracks[itrack]->GetNHit(); ih < n; ++ih){
          const Double_t d =
            (tracks[itrack]->GetHitInOrder(ih)->GetLocalCalPosHelix() - vtx_pos).Mag();
          if(d < cluster_dist_max){
            close_hit_orders[itrack].push_back(tracks[itrack]->GetOrder(ih));
          }
          if(d < min_dist){
            closest_hit_idx[itrack] = ih;
            min_dist = d;
          }
        }
      }

      // Skip ambiguous geometry: both tracks' nearest hit is InOrder(0) (track end)
      // and neither has exactly one cluster in cluster_dist_max.
      if(closest_hit_idx[0] == 0 && closest_hit_idx[1] == 0 &&
         close_hit_orders[0].size() != 1 && close_hit_orders[1].size() != 1) continue;

      if(close_hit_orders[0].empty() && close_hit_orders[1].empty()) continue;

      // Strip nearby hits; build trial tracks (not yet in TrackCont)
      std::vector<T*> trial_tracks;
      std::vector<Int_t> candidates_trackid;
      std::vector<TPCHit*> clusters_fortest;
      std::vector<Int_t> clusters_original_trackid;

      Bool_t strip_ok = true;
      for(Int_t itrack = 0; itrack < 2; ++itrack){
        T* trial = new T(tracks[itrack]);
        trial->EraseHits(close_hit_orders[itrack]);
        if(!trial->DoFit(MinNumOfHits)){
          del::ClearContainer(trial_tracks);
          delete trial;
          for(Int_t t = 0; t <= itrack; ++t){
            tracks[t]->SetClustersHoughFlag(GoodForTracking);
          }
          strip_ok = false;
          break;
        }
        trial_tracks.push_back(trial);
        candidates_trackid.push_back(trackids[itrack]);
        for(Int_t ord : close_hit_orders[itrack]){
          TPCHit* hit = tracks[itrack]->GetHit(ord)->GetHit();
          clusters_fortest.push_back(hit);
          clusters_original_trackid.push_back(trackids[itrack]);
        }
      }
      if(!strip_ok) continue;

      // Assign each stripped cluster to the trial with smallest IsGoodHitToAdd residual
      Bool_t reassign = false;
      std::vector<Int_t> reassign_candidates;
      for(Int_t icl = 0, ncl = clusters_fortest.size(); icl < ncl; ++icl){
        TPCHit* hit = clusters_fortest[icl];
        hit->SetHoughFlag(Candidate);

        Int_t best_idx = -9999;
        Int_t best_trackid = -1;
        Double_t min_residual = 1.e10;
        for(Int_t j = 0, nt = trial_tracks.size(); j < nt; ++j){
          Double_t residual = 0.;
          if(trial_tracks[j]->IsGoodHitToAdd(hit, residual) && min_residual > residual){
            min_residual = residual;
            best_idx = j;
            best_trackid = candidates_trackid[j];
          }
        }

        const Int_t orig_trackid = clusters_original_trackid[icl];
        if(orig_trackid == best_trackid){
          hit->SetHoughFlag(GoodForTracking);
          trial_tracks[best_idx]->AddTPCHit(new TPCLTrackHit(hit));
        }
        else{
          reassign = true;
          if(best_trackid == -1){
            reassign_candidates.push_back(orig_trackid);
            hit->SetHoughFlag(0);
          }
          else{
            reassign_candidates.push_back(orig_trackid);
            reassign_candidates.push_back(best_trackid);
            trial_tracks[best_idx]->AddTPCHit(new TPCLTrackHit(hit));
          }
        }
      } //for(icl : clusters_fortest)

      if(trial_tracks.size() != candidates_trackid.size())
        std::cout << FUNC_NAME + " FATAL Error #of tracks is not matched" << std::endl;

      std::sort(reassign_candidates.begin(), reassign_candidates.end());
      reassign_candidates.erase(
        std::unique(reassign_candidates.begin(), reassign_candidates.end()),
        reassign_candidates.end());

      if(!reassign){
        del::ClearContainer(trial_tracks);
        tracks[0]->SetClustersHoughFlag(GoodForTracking);
        tracks[1]->SetClustersHoughFlag(GoodForTracking);
        continue;
      }

      // Commit only tracks listed in reassign_candidates; DoFit is the final gate
      Bool_t any_committed = false;
      for(Int_t j = 0, nt = trial_tracks.size(); j < nt; ++j){
        const Int_t trackid = candidates_trackid[j];
        T* trial = trial_tracks[j];

        if(std::find(reassign_candidates.begin(), reassign_candidates.end(), trackid)
           == reassign_candidates.end()){
          delete trial;
          continue;
        }

        T* prev_track = TrackCont[trackid];
        if(trial->DoFit(MinNumOfHits)){
          trial->SetClustersHoughFlag(GoodForTracking);
          trial->Calculate();
          if(Exclusive){
            trial->DoFitExclusive();
            trial->CalculateExclusive();
          }
          TrackCont[trackid] = trial;
          delete prev_track;
          any_committed = true;
        }
        else{
          prev_track->SetClustersHoughFlag(GoodForTracking);
          delete trial;
        }
      }

      if(!any_committed){
        tracks[0]->SetClustersHoughFlag(GoodForTracking);
        tracks[1]->SetClustersHoughFlag(GoodForTracking);
        continue;
      }

      del::ClearContainer(VertexCont);
      VertexSearch(TrackCont, VertexCont);
      candidate_pairs.push_back(std::make_pair(trackids[0], trackids[1]));
      reassign_pending = true;  // rescan VertexCont after geometry update
      break;
    } //for(vertex : VertexCont)
  } //while(reassign_pending)

  ResetHoughFlag(ClCont);

#if FragmentedTrackTest
  //Merged fragmented tracks
  RestoreFragmentedTracks(ClCont, TrackCont, TrackContFailed, VertexCont, Exclusive, MinNumOfHits);
#endif
  if(!BeamThroughTPC) MarkingBeamTracks(TrackCont);
  if(!BeamThroughTPC) MarkingAccidentalTracks(TrackCont);

}

//_____________________________________________________________________________
template <typename T> void
FindAccidentalCoincidenceTracks(std::vector<T*>& TrackCont,
				std::vector<TPCVertex*>& VertexCont,
				std::vector<TPCVertex*>& ClusteredVertexCont)
{
  TVector3 tgt(0., 0., tpc::Z_TARGET);
  std::vector<std::vector<TPCVertex*>> clustered_vertices;
  std::vector<std::vector<Int_t>> clustered_tracks;
  std::vector<TVector3> accidental_vertices;

  std::vector<Int_t> candidates_trackids;

  const Int_t n_cluster_stages = 2;
  const Int_t ntracks = TrackCont.size();

  for(Int_t stage = 0; stage < n_cluster_stages; ++stage){
    // stage 0: beam track as seed — 1st pass
    // stage 1: non-beam seeds not yet in candidates — 2nd pass (tighter GetClosestDist cut)
    if(stage >= 1){
      std::sort(candidates_trackids.begin(), candidates_trackids.end());
      candidates_trackids.erase(
        std::unique(candidates_trackids.begin(), candidates_trackids.end()),
        candidates_trackids.end());
    }

    // empirical max GetClosestDist() [mm]; stage-specific, not shared with other cuts
    const Double_t max_vertex_closest_dist = (stage == 0) ? 20. : 15.;

    for(Int_t seed_trackid = 0; seed_trackid < ntracks; ++seed_trackid){
      T* seed_track = TrackCont[seed_trackid];
      if(seed_track->GetIsK18() == 1) continue;

      // seed track filters (depend on stage)
      if(stage == 0 && seed_track->GetIsBeam() != 1) continue;
      if(stage >= 1){
        if(std::find(candidates_trackids.begin(), candidates_trackids.end(), seed_trackid)
           != candidates_trackids.end()) continue;
        if(seed_track->GetIsBeam() == 1) continue;
        TVector3 seed_residual_tgt = seed_track->GetClosestPositionTgtXZ() - tgt;
        if(TMath::Hypot(seed_residual_tgt.x(), seed_residual_tgt.z()) > tpc::TARGET_RADIUS) continue;
      }

      std::vector<TPCVertex*> candidate_clustered_vertices;
      std::vector<Int_t> candidate_clustered_tracks;
      TVector3 vertex_point(0, 0, 0);

      candidate_clustered_tracks.push_back(seed_trackid);
      for(auto& vertex : VertexCont){
        // skip unless this vertex involves the seed track
        const Int_t trackid0 = vertex->GetTrackId(0);
        const Int_t trackid1 = vertex->GetTrackId(1);
        Int_t coin_trackid = -1;
        if(trackid0 == seed_trackid) coin_trackid = trackid1;
        else if(trackid1 == seed_trackid) coin_trackid = trackid0;
        else continue;

        if(vertex->GetClosestDist() > max_vertex_closest_dist) continue;
        TVector3 vtx = vertex->GetVertex() - tgt;
        if(TMath::Hypot(vtx.x(), vtx.z()) > tpc::TARGET_RADIUS) continue;

        T* coin_track = TrackCont[coin_trackid];
        if(coin_track->GetIsK18() == 1) continue;
        if(stage >= 1 && coin_track->GetIsBeam() == 1) continue;

        TVector3 coin_residual_tgt = coin_track->GetClosestPositionTgtXZ() - tgt;
        if(TMath::Hypot(coin_residual_tgt.x(), coin_residual_tgt.z()) > tpc::TARGET_RADIUS) continue;

        candidate_clustered_vertices.push_back(vertex);
        candidate_clustered_tracks.push_back(coin_trackid);
        vertex_point += coin_residual_tgt;
        candidates_trackids.push_back(coin_trackid);
      }

      if(candidate_clustered_tracks.size() < 2) continue;
      candidates_trackids.push_back(seed_trackid);

      const Double_t ncltrk = static_cast<Double_t>(candidate_clustered_tracks.size());
      clustered_vertices.push_back(candidate_clustered_vertices);
      clustered_tracks.push_back(candidate_clustered_tracks);
      accidental_vertices.push_back(
        TVector3(vertex_point.x()/ncltrk,
                 vertex_point.y()/ncltrk,
                 vertex_point.z()/ncltrk));
    }
  }

  for(Int_t id=0; id<accidental_vertices.size(); id++){
    if(TMath::Abs(accidental_vertices[id].y()) < tpc::TARGET_HALF_Y) continue;
    TPCVertex *clustered_vtx = new TPCVertex(accidental_vertices[id], clustered_tracks[id]);
    ClusteredVertexCont.push_back(clustered_vtx);
  }

  MarkingClusteredAccidentalTracks(TrackCont, ClusteredVertexCont);

}

//_____________________________________________________________________________
template <typename T> void
TestingCharge(std::vector<T*>& TrackCont,
	      std::vector<T*>& TrackContInvertedCharge,
	      std::vector<TPCVertex*>& VertexCont,
	      Bool_t Exclusive)
{

  // Try charge inversion per track: inclusive chi2, else Lambda rescue (ppi_distcut ~ 10 mm).
  // TrackContInvertedCharge[i] is set only for charge-tested tracks; else nullptr.
  // VertexSearch runs only if at least one track adopts the inverted charge.
  TrackContInvertedCharge.assign(TrackCont.size(), nullptr);
  Bool_t any_inverted = false;
  for(Int_t trackid=0; trackid<TrackCont.size(); trackid++){
    T *track = TrackCont[trackid];
    if(track->GetIsAccidental()==1 || track->GetIsK18()==1 ||
       track->GetIsBeam()==1) continue;

    T *InvertedTrack = new T(track);
    Bool_t test = InvertedTrack->TestInvertCharge();
    if(!test || (track->GetCharge() == InvertedTrack->GetCharge())){
      delete InvertedTrack;
      continue;
    }

    InvertedTrack->Calculate();
    if(Exclusive){
      InvertedTrack->DoFitExclusive();
      InvertedTrack->CalculateExclusive();
    }

    Bool_t invert = false;
    if(track->GetChiSquare() > InvertedTrack->GetChiSquare()) invert = true;
    else{
      for(Int_t i=0; i<VertexCont.size(); i++){
        TPCVertex *vertex = VertexCont[i];
        Double_t closedist = vertex->GetClosestDist();
        if(closedist < ppi_distcut){
          Int_t trackid1 = vertex->GetTrackId(0);
          Int_t trackid2 = vertex->GetTrackId(1);
          if(trackid1 != trackid && trackid2 != trackid) continue;

          T *track1 = TrackCont[trackid1];
          T *track2 = TrackCont[trackid2];
          if(trackid1 == trackid) track1 = InvertedTrack;
          else if(trackid2 == trackid) track2 = InvertedTrack;

          TPCVertex *newvertex = new TPCVertex(trackid1, trackid2);
          newvertex->Calculate(track1, track2);
          if(TPCReconstructor::HasLambdaCandidate(newvertex)){
            invert = true;
            delete newvertex;
            break;
          }
          delete newvertex;
        }
      }
    }

    if(invert){
      TrackCont[trackid] = InvertedTrack;
      TrackContInvertedCharge[trackid] = track;
      any_inverted = true;
    }
    else{
      TrackCont[trackid] = track;
      TrackContInvertedCharge[trackid] = InvertedTrack;
    }
  }

  if(any_inverted){
    del::ClearContainer(VertexCont);
    VertexSearch(TrackCont, VertexCont);
  }
}

//_____________________________________________________________________________
template <typename T> void
RestoreFragmentedAccidentalTracks(const std::vector<TPCClusterContainer>& ClCont,
				  std::vector<T*>& TrackCont,
				  std::vector<T*>& TrackContFailed,
				  std::vector<TPCVertex*>& VertexCont,
				  Bool_t Exclusive,
				  Int_t MinNumOfHits)
{

  Bool_t reassign = false;

  const Int_t max_fragment_hits = 4;
  const Double_t min_mom_negative = 0.5;
  const Double_t min_mom_positive = 1.0;

  // case 1: merge a short upstream fragment with its downstream parent track at a vertex
  std::vector<Int_t> erase_trackid;
  std::vector<Int_t> new_trackid;

  for(auto& vertex : VertexCont){
    Int_t id0 = vertex->GetTrackId(0);
    Int_t id1 = vertex->GetTrackId(1);

    if(std::find(erase_trackid.begin(), erase_trackid.end(), id0) != erase_trackid.end()) continue;
    if(std::find(erase_trackid.begin(), erase_trackid.end(), id1) != erase_trackid.end()) continue;
    if(std::find(new_trackid.begin(), new_trackid.end(), id0) != new_trackid.end()) continue;
    if(std::find(new_trackid.begin(), new_trackid.end(), id1) != new_trackid.end()) continue;

    T *track0 = TrackCont[id0];
    T *track1 = TrackCont[id1];

    // assign parent (longer) and fragment (shorter) by hit count
    if(track0->GetNHit() == track1->GetNHit()) continue;

    Int_t parent_id = (track0->GetNHit() > track1->GetNHit()) ? id0 : id1;
    Int_t fragment_id = (track0->GetNHit() > track1->GetNHit()) ? id1 : id0;
    T *parent = TrackCont[parent_id];
    T *fragment = TrackCont[fragment_id];

    // parent: downstream body candidate (slope + charge-dependent mom cuts)
    if(parent->GetIsAccidental()==1 || parent->GetIsBeam()==1 || parent->GetIsK18()==1) continue;
    if(TMath::Abs(parent->Getdz()) > MaxCandidateAbsDz) continue;

    Double_t helix_mom = parent->GetMom0().Mag();
    const Double_t min_mom = (parent->GetCharge() > 0) ? min_mom_positive : min_mom_negative;
    if(TMath::Abs(helix_mom) < min_mom) continue;

    // fragment: short upstream stub
    if(fragment->GetIsK18()==1) continue;
    if(fragment->GetNHit() > max_fragment_hits) continue;

    Bool_t cl_added = false;
    fragment->SetClustersHoughFlag(0);
    TPCLocalTrackHelix *merged_track = new TPCLocalTrackHelix(parent);
    for(Int_t hit_idx=0; hit_idx<fragment->GetNHit(); hit_idx++){
      TPCHit *cl = fragment->GetHitInOrder(hit_idx)->GetHit();
      Double_t resi = 0.;
      Bool_t no_limitation = true;
      if(merged_track->IsGoodHitToAdd(cl, resi, no_limitation)){
        merged_track->AddTPCHit(new TPCLTrackHit(cl));
        cl_added = true;
      }
    }

    // Check whether tracks belong to the same track or not
    if(!cl_added || (fragment->GetNHit() + parent->GetNHit() - merged_track->GetNHit() > 2)){
      fragment->SetClustersHoughFlag(GoodForTracking);
      delete merged_track;
      continue;
    }

    merged_track->SetIsAccidental(); // needed to turn off VertexAtTarget()

    Int_t n_tracks_before = TrackCont.size();
    FitTrack(merged_track, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
    Int_t n_tracks_after = TrackCont.size();
    if(n_tracks_before + 1 == n_tracks_after){
      merged_track = TrackCont[n_tracks_after - 1];
      TrackCont.pop_back();
      merged_track->Calculate();
      merged_track->CheckIsAccidental();
      if(Exclusive){
        merged_track->DoFitExclusive();
        merged_track->CalculateExclusive();
      }

      if(merged_track->GetNHit() > parent->GetNHit() &&
         merged_track->GetIsAccidental()==1){
        TrackCont.push_back(merged_track);
        new_trackid.push_back(n_tracks_after - 1);
        erase_trackid.push_back(parent_id);
        erase_trackid.push_back(fragment_id);
        reassign = true;
      }
      else{
        parent->SetClustersHoughFlag(GoodForTracking);
        fragment->SetClustersHoughFlag(GoodForTracking);
        delete merged_track;
      }
    }
    else{ // Fitting is failed
      parent->SetClustersHoughFlag(GoodForTracking);
      fragment->SetClustersHoughFlag(GoodForTracking);
    }
  }

  std::sort(erase_trackid.begin(), erase_trackid.end(), std::greater<Int_t>());
  for(Int_t id : erase_trackid){
    T *track = TrackCont[id];
    TrackCont.erase(TrackCont.begin() + id);
    delete track;
  }

  // case 2: accidental track is divided into clusters on the upstream of the target and a track on the downstream of the target
  std::vector<TPCClusterContainer> candidate_cl_cont(10);
  for(Int_t layer=0; layer<10; layer++){ //inner layers 0..9
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(hit->GetHoughFlag()==GoodForTracking ||
         hit->GetHoughFlag()==K18Tracks) continue;
      TVector3 pos = cl->GetPosition();
      // upstream orphan clusters: z gate + naive x-only transverse pre-filter
      if(pos.Z() > tpc::Z_TARGET) continue;
      if(TMath::Abs(pos.x()) > tpc::TARGET_HALF_X) continue;
      candidate_cl_cont[layer].push_back(cl);
    } //ci
  } //layer

  for(Int_t trackid=0; trackid<TrackCont.size(); trackid++){
    T *parent = TrackCont[trackid];
    if(parent->GetIsAccidental()==1 || parent->GetIsBeam()==1 || parent->GetIsK18()==1) continue;
    if(TMath::Abs(parent->Getdz()) > MaxCandidateAbsDz) continue;

    Double_t helix_mom = parent->GetMom0().Mag();
    const Double_t min_mom = (parent->GetCharge() > 0) ? min_mom_positive : min_mom_negative;
    if(TMath::Abs(helix_mom) < min_mom) continue;

    T *merged_track = new T(parent);
    std::vector<TPCHit*> clusters_added;
    Bool_t cl_added = false;
    for(Int_t layer=0; layer<10; layer++){ //inner layers 0..9
      for(Int_t ci=0, n=candidate_cl_cont[layer].size(); ci<n; ci++){
        auto cl = candidate_cl_cont[layer][ci];
        TPCHit* hit = cl->GetMeanHit();

        Double_t resi = 0.;
        Bool_t no_limitation = true;
        if(parent->IsGoodHitToAdd(hit, resi, no_limitation)){
          merged_track->AddTPCHit(new TPCLTrackHit(hit));
          clusters_added.push_back(hit);
          cl_added = true;
        }
      }
    }

    if(!cl_added){
      delete merged_track;
      continue;
    }

    merged_track->SetIsAccidental(); // needed to turn off VertexAtTarget()

    Int_t n_tracks_before = TrackCont.size();
    FitTrack(merged_track, GoodForTracking, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
    Int_t n_tracks_after = TrackCont.size();

    Bool_t accept = false;
    if(n_tracks_before + 1 == n_tracks_after){
      merged_track = TrackCont[n_tracks_after - 1];
      TrackCont.pop_back();
      if(Exclusive){
        merged_track->DoFitExclusive();
        merged_track->CalculateExclusive();
      }
      merged_track->Calculate();
      merged_track->CheckIsAccidental();
      if(merged_track->GetNHit() > parent->GetNHit() &&
         merged_track->GetIsAccidental()==1){
        accept = true;
      }
    }

    if(accept){
      TrackCont[trackid] = merged_track;
      delete parent;
      reassign = true;
    }
    else{
      for(TPCHit *cl : clusters_added) cl->SetHoughFlag(0);
      parent->SetClustersHoughFlag(GoodForTracking);
      if(n_tracks_before + 1 == n_tracks_after) delete merged_track;
    }
  }

  if(reassign){
    //Vertex finding again with new tracks
    del::ClearContainer(VertexCont);
    VertexSearch(TrackCont, VertexCont);
  } //reassigning process
}

} //namespace tpc
