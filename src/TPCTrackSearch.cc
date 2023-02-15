// -*- C++ -*-

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

#define DebugDisp      0
#define RemainHitTest  0

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

  // Houghflags
  const Int_t    GoodForTracking = 100;
  const Int_t    K18Tracks = 200;
  const Int_t    KuramaTracks = 300;
  const Int_t    AccidentalBeams = 400;
  const Int_t    BadForTracking = 500;
  const Int_t    Candidate = 1000;

  // Tracks in the Hough-Space(circle)
  std::vector<Double_t> hough_x;
  std::vector<Double_t> hough_y;
  std::vector<Double_t> hough_z;

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
  Double_t GetMagneticField()
  {
    return HS_field_0*(valueHSHall/valueHSCalc);
  }
}

//_____________________________________________________________________________
namespace tpc
{

//_____________________________________________________________________________
void
FitLinearTrack(TPCLocalTrack *Track,
	       const std::vector<TPCClusterContainer>& ClCont,
	       std::vector<TPCLocalTrack*>& TrackCont,
	       std::vector<TPCLocalTrack*>& TrackContFailed,
	       Int_t MinNumOfHits)
{

  std::chrono::milliseconds sec;
  auto fit_start = std::chrono::high_resolution_clock::now();

#if DebugDisp
  Track->Print(FUNC_NAME+" Initial track", true);
#endif

  if(Track->DoFit(MinNumOfHits)){ //fitting succeeded
    auto good = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(good - fit_start);
    Track->SetFitTime(sec.count());
    Track->SetHoughFlag(GoodForTracking + TrackCont.size());
    Track->SetFitFlag(1);

#if DebugDisp
    Track->Print(FUNC_NAME+" Fitting is succeeded");
#endif
    TrackCont.push_back(Track);
  }
  else{ //fitting failed
    auto bad = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(bad - fit_start);
    Track->SetFitTime(sec.count());
    Track->SetHoughFlag(BadForTracking);
    Track->SetParamUsingHoughParam();
    Track->SetFitFlag(4);

#if DebugDisp
    Track->Print(FUNC_NAME+" Fitting is failed");
#endif
    if(Track->GetNHit() > 0.5*MinNumOfHits) TrackContFailed.push_back(Track);
    else delete Track;
  }
}

//_____________________________________________________________________________
Int_t
LocalTrackSearch(const std::vector<TPCClusterContainer>& ClCont,
		 std::vector<TPCLocalTrack*>& TrackCont,
		 std::vector<TPCLocalTrack*>& TrackContFailed,
		 Int_t MinNumOfHits /*=8*/)
{
  const Double_t MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");

  hough_x.clear();
  hough_y.clear();
  hough_z.clear();

  bool prev_add = true;
  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    if(!prev_add) continue;
    prev_add = false;

#if DebugDisp
    std::cout<<FUNC_NAME+" tracki : "<<tracki<<std::endl;
#endif

    std::chrono::milliseconds sec;
    auto before_hough = std::chrono::high_resolution_clock::now();

    //Line Hough-transform
    Double_t LinearPar[2]; Int_t MaxBin[3];
    if(!tpc::HoughTransformLineXZ(ClCont, MaxBin, LinearPar, MinNumOfHits)){
#if DebugDisp
      std::cout<<FUNC_NAME+" No more track candiate! tracki : "<<tracki<<std::endl;
#endif
      break;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<hough_x.size(); ++i){
      Int_t bindiff = TMath::Abs(MaxBin[0] - hough_x[i]) + TMath::Abs(MaxBin[1] - hough_y[i]) + TMath::Abs(MaxBin[2] - hough_z[i]);
      if(bindiff<=2){
	hough_flag = false;
#if DebugDisp
	std::cout<<"Previous hough bin "<<i<<"th x: "<<hough_x[i]
		 <<", y: "<<hough_y[i]<<", z: "<<hough_z[i]<<std::endl;
	std::cout<<"Current hough bin "<<i<<"th x: "<<MaxBin[0]
		 <<", y: "<<MaxBin[1]<<", z: "<<MaxBin[2]<<std::endl;
#endif
      }
    }
    hough_x.push_back(MaxBin[0]);
    hough_y.push_back(MaxBin[1]);
    hough_z.push_back(MaxBin[2]);
    if(!hough_flag){
#if DebugDisp
      std::cout<<FUNC_NAME+" The same track is found by Hough-Transform : tracki : "<<tracki<<" hough_x size : "<<hough_x.size()<<std::endl;
#endif
      continue;
    }

    TPCLocalTrack *track = new TPCLocalTrack;
    //(x0, y0) are position at Target position
    track->SetAx(LinearPar[0] + LinearPar[1]*tpc::ZTarget);
    track->SetAu(LinearPar[1]);

    //HoughDistCheck
    prev_add = LinearHoughDistCheck(track, ClCont, LinearPar, MaxHoughWindowY);
    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_hough - before_hough);
    track->SetSearchTime(sec.count());

    //Track fitting processes
    FitLinearTrack(track, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  }// tracki

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

  CalcTracks(TrackCont);
  CalcTracks(TrackContFailed);

  return TrackCont.size();
}

//_____________________________________________________________________________
void
ExclusiveTracking(std::vector<TPCLocalTrack*>& TrackCont)
{

  for(auto& track: TrackCont){
    track->DoLinearFitExclusive();
    track->CalculateExclusive();
  }
}

//_____________________________________________________________________________
void
ExclusiveTrackingHelix(std::vector<TPCLocalTrackHelix*>& TrackCont)
{

  for(auto& track: TrackCont){
    track->DoHelixFitExclusive();
    track->CalculateExclusive();
  }
}

//_____________________________________________________________________________
Bool_t
LinearHoughDistCheck(TPCLocalTrack *track,
		     const std::vector<TPCClusterContainer>& ClCont,
		     Double_t *LinearPar, Double_t MaxHoughWindowY){

  bool status = false;
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(hit->GetHoughFlag()>0) continue;
      TVector3 pos = cl->GetPosition();
      Double_t dist = TMath::Abs(LinearPar[1]*pos.Z()-pos.X() +
				 LinearPar[0])/TMath::Sqrt(TMath::Sq(LinearPar[1])+1.);
      if(dist < MaxHoughWindowY){
	hit->SetHoughDistY(dist);
	track->AddTPCHit(new TPCLTrackHit(hit));
	status = true;
      }
    } //ci
  } //layer

  return status;
}

//_____________________________________________________________________________
Bool_t
HelixHoughDistCheck(TPCLocalTrackHelix *track,
		    const std::vector<TPCClusterContainer>& ClCont,
		    Double_t *HelixPar, Double_t MaxHoughWindow,
		    Double_t MaxHoughWindowY)
{

  bool status = false;
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
	  track->AddTPCHit(new TPCLTrackHit(hit));
	  status = true;
	} //distY
      } //dist
    } //ci
  } //layer

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
FitHelixTrack(TPCLocalTrackHelix *Track, Int_t Houghflag,
	      const std::vector<TPCClusterContainer>& ClCont,
	      std::vector<TPCLocalTrackHelix*>& TrackCont,
	      std::vector<TPCLocalTrackHelix*>& TrackContFailed,
	      Int_t MinNumOfHits)
{

  std::chrono::milliseconds sec;
  auto fit_start = std::chrono::high_resolution_clock::now();

#if DebugDisp
  Track->Print(FUNC_NAME+" Initial track", true);
#endif

  if(Track->DoFit(MinNumOfHits)){ //first fitting
    Track->SetHoughFlag(Candidate);
#if DebugDisp
    Track->Print(FUNC_NAME+" First fitting is succeeded");
#endif
    auto first_fit = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(first_fit - fit_start);
    Track->SetFitTime(sec.count());

    TPCLocalTrackHelix *copied_track = new TPCLocalTrackHelix(Track);
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
	if(copied_track->ResidualCheck(pos,res,resi)){
	  copied_track->AddTPCHit(new TPCLTrackHit(hit));
	  Add_rescheck = true;
	}
      }
    } //for (layer)
    if(!Add_rescheck){
      Track->SetHoughFlag(Houghflag);
      Track->SetFitFlag(1);

#if DebugDisp
      Track->Print(FUNC_NAME+" No added cluster");
#endif
      TrackCont.push_back(Track);
      delete copied_track;
    }
    else{
      auto hitadd = std::chrono::high_resolution_clock::now();
#if DebugDisp
      copied_track->Print(FUNC_NAME+" After cluster adding");
#endif
      if(copied_track->DoFit(MinNumOfHits)){ //2nd fitting success
	auto second_fit = std::chrono::high_resolution_clock::now();
	sec = std::chrono::duration_cast<std::chrono::milliseconds>(second_fit - hitadd);
	copied_track->SetFitTime(sec.count());
	copied_track->SetHoughFlag(Houghflag);
	copied_track->SetFitFlag(2);

#if DebugDisp
	copied_track->Print(FUNC_NAME+" Fitting success after adding hits (residual check)", true);
#endif
	TrackCont.push_back(copied_track);
	delete Track;
      }
      else{ //2nd fitting failure
	auto hitadd_fail = std::chrono::high_resolution_clock::now();
	sec = std::chrono::duration_cast<std::chrono::milliseconds>(hitadd_fail - hitadd);
	copied_track->SetFitTime(sec.count());
	copied_track->SetParamUsingHoughParam();
	copied_track->SetHoughFlag(0);
	copied_track->SetFitFlag(3);
	copied_track->SetIsK18(0);
	copied_track->SetIsKurama(0);

#if DebugDisp
	copied_track->Print(FUNC_NAME+" 2nd fitting failure!", true);
#endif
	if(copied_track->GetNHit() > 0.5*MinNumOfHits){
	  TrackContFailed.push_back(copied_track);
	}
	else delete copied_track;

	//push 1st fitting track
	Track->SetFitTime(sec.count());
	Track->SetHoughFlag(Houghflag);
	Track->SetFitFlag(3);
	TrackCont.push_back(Track);
      }
    }
  }
  else{ //first fitting failed
    auto after_1stfit_fail = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_1stfit_fail - fit_start);
    Track->SetFitTime(sec.count());
    Track->SetHoughFlag(0);
    Track->SetParamUsingHoughParam();
    Track->SetFitFlag(4);
    Track->SetIsK18(0);
    Track->SetIsKurama(0);

#if DebugDisp
    Track->Print(FUNC_NAME+" First fitting is failed", true);
#endif

    if(Track->GetNHit() > 0.5*MinNumOfHits){
      Track->SetHoughFlag(BadForTracking);
      TrackContFailed.push_back(Track);
    }
    else delete Track;
  }

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

}

//_____________________________________________________________________________
void
HelixTrackSearch(Int_t Beamflag, Int_t Houghflag,
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

    //Circle Hough-transform
    std::chrono::milliseconds sec;
    auto before_hough = std::chrono::high_resolution_clock::now();
    Double_t HelixPar[5]; Int_t MaxBin[3];
    if(!tpc::HoughTransformCircleXZ(ClCont, MaxBin, HelixPar, MinNumOfHits)){
#if DebugDisp
      std::cout<<FUNC_NAME+" No more circle candiate! tracki : "<<tracki<<std::endl;
#endif
      break;
    }

    //Check for duplicates
    Bool_t hough_flag = true;
    for(Int_t i=0; i<hough_x.size(); ++i){
      Int_t bindiff = TMath::Abs(MaxBin[0] - hough_x[i]) + TMath::Abs(MaxBin[1] - hough_y[i]) + TMath::Abs(MaxBin[2] - hough_z[i]);
      if(bindiff<=2){
	hough_flag = false;
#if DebugDisp
	std::cout<<"Previous hough bin "<<i<<"th x: "<<hough_x[i]
		 <<", y: "<<hough_y[i]<<", z: "<<hough_z[i]<<std::endl;
	std::cout<<"Current hough bin "<<i<<"th x: "<<MaxBin[0]
		 <<", y: "<<MaxBin[1]<<", z: "<<MaxBin[2]<<std::endl;
#endif
      }
    }
    if(!hough_flag){
#if DebugDisp
      std::cout<<FUNC_NAME+" The same track is found by Hough-Transform : tracki : "<<tracki<<" hough_x size : "<<hough_x.size()<<std::endl;
#endif
      continue;
    }
    hough_x.push_back(MaxBin[0]);
    hough_y.push_back(MaxBin[1]);
    hough_z.push_back(MaxBin[2]);

    //Linear Hough-transform
    tpc::HoughTransformLineYPhi(ClCont, HelixPar, MaxHoughWindow);

    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    track->SetAcx(HelixPar[0]);
    track->SetAcy(HelixPar[1]);
    track->SetAz0(HelixPar[2]);
    track->SetAr(HelixPar[3]);
    track->SetAdz(HelixPar[4]);
    track->SetFlag(Beamflag);

    //HoughDistCheck
    prev_add = HelixHoughDistCheck(track, ClCont, HelixPar, MaxHoughWindow, MaxHoughWindowY);

    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::milliseconds>(after_hough - before_hough);
    track->SetSearchTime(sec.count());

    //Track fitting processes
    FitHelixTrack(track, Houghflag, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  }//tracki
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

    Int_t Beamflag = 1*0 + 2*0 + 4*1 + 8*0; // isBeam, isK18, isKurama, isAccidental
    HelixTrackSearch(Beamflag, KuramaTracks, ClContKurama, TrackCont, TrackContFailed, MinNumOfHits);
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
    if(BeforeTGTHits<3) break;

#if 0 //undev
    //case 1 : Beam-through
    if(AfterTGTHits>MinNumOfHits){
      Int_t Beamflag = 1*1 + 2*1 + 4*0 + 8*1; // isBeam, isK18, isKurama, isAccidental
      HelixTrackSearch(Beamflag, K18Tracks, ClContK18, TrackCont, TrackContFailed, MinNumOfHits+3);
      ResetHoughFlag(ClContK18);
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

    Int_t Beamflag = 1*1 + 2*1 + 4*0 + 8*0; // isBeam, isK18, isKurama, isAccidental
    track->SetFlag(Beamflag);
    track->SetSearchTime(0.);

    //Track fitting processes
    FitHelixTrack(track, K18Tracks, ClContK18, TrackCont, TrackContFailed, 3);
    ResetHoughFlag(ClContK18);
  } //nt

}

//_____________________________________________________________________________
void
AccidentalBeamSearch(const std::vector<TPCClusterContainer>& ClCont,
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
    Int_t Beamflag = 1*1 + 2*0 + 4*0 + 8*1; // isBeam, isK18, isKurama, isAccidental
    HelixTrackSearch(Beamflag, AccidentalBeams, ClContBeam, TrackCont, TrackContFailed, MinNumOfHits);
    ResetHoughFlag(ClContBeam);
  } //candidates

}

//_____________________________________________________________________________
Int_t
LocalTrackSearchHelix(std::vector<std::vector<TVector3>> K18VPs,
		      std::vector<std::vector<TVector3>> KuramaVPs,
		      const std::vector<TPCClusterContainer>& ClCont,
		      std::vector<TPCLocalTrackHelix*>& TrackCont,
		      std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		      Int_t MinNumOfHits /*=8*/)
{

  hough_x.clear();
  hough_y.clear();
  hough_z.clear();

  //Track searching and fitting
  //for Accidental beams and K1.8 & Kurama tracks
  AccidentalBeamSearch(ClCont, TrackCont, TrackContFailed, 12);
  K18TrackSearch(K18VPs, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  KuramaTrackSearch(KuramaVPs, ClCont, TrackCont, TrackContFailed, MinNumOfHits);
  //for scattered helix tracks
  HelixTrackSearch(0, GoodForTracking , ClCont, TrackCont, TrackContFailed, MinNumOfHits);

#if RemainHitTest
  TestRemainingHits(ClCont, TrackCont, MinNumOfHits);
#endif

#if DebugDisp
  std::cout<<FUNC_NAME+" #track : "<<TrackCont.size()<<std::endl;
  std::cout<<FUNC_NAME+" #failed track : "<<TrackContFailed.size()<<std::endl;
#endif

  CalcTracks(TrackCont);
  CalcTracks(TrackContFailed);

  return TrackCont.size();
}

} //namespace tpc
