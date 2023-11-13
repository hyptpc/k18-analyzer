//  Authors: Wooseung Jung

//GenKEK
#include "HypTPCTrack.hh"
#include "HypTPCHit.hh"

//k18-analyzer
#include "TPCLTrackHit.hh"
#include "TPCHit.hh"
#include "DeleteUtility.hh"

//GenFit
#include <AbsFitterInfo.h>
#include <AbsHMatrix.h>
#include <AbsTrackRep.h>
#include <DetPlane.h>
#include <Exception.h>
#include <FitStatus.h>
#include <KalmanFittedStateOnPlane.h>
#include <KalmanFitterInfo.h>
#include <MeasuredStateOnPlane.h>
#include <MeasurementOnPlane.h>
#include <SharedPlanePtr.h>
#include <Tools.h>
#include <TrackPoint.h>
#include <RKTrackRep.h>

//ROOT
#include <TMatrixDfwd.h>                      // for TMatrixD
#include <TMatrixDSymfwd.h>                   // for TMatrixDSym
#include <TMatrixT.h>                         // for TMatrixT
#include <TMatrixTSym.h>                      // for TMatrixTSym
#include <TVectorDfwd.h>                      // for TVectorD
#include <TVector3.h>                         // for TVector3
#include <TVectorT.h>                         // for TVectorT, operator-

//STL
#include <iostream>

ClassImp(HypTPCTrack)

TClonesArray *HypTPCTrack::_hitClusterArray = nullptr;
TClonesArray *HypTPCTrack::_genfitTrackArray = nullptr;

HypTPCTrack::HypTPCTrack(){

  _hitClusterArray = new TClonesArray("HypTPCHit");
  _genfitTrackArray = new TClonesArray("genfit::Track");
  _measurementProducer = new genfit::MeasurementProducer<HypTPCHit, genfit::HypTPCSpacepointMeasurement>(_hitClusterArray);
  _measurementFactory = new genfit::MeasurementFactory<genfit::AbsMeasurement>();
  _measurementFactory -> addProducer(TPCDetID, _measurementProducer);
  std::cout<<"GenFit : HypTPCTrack container is ready"<<std::endl;
}

void HypTPCTrack::Clear(){
  _genfitTrackArray -> Delete();
  _hitClusterArray -> Delete();
}

genfit::Track* HypTPCTrack::GetTrack(int ith) const {

  return (genfit::Track*) _genfitTrackArray -> ConstructedAt(ith);
}

void HypTPCTrack::AddReps(int ith, int pdg){

  genfit::Track* tr = GetTrack(ith);
  tr -> addTrackRep(new genfit::RKTrackRep(pdg));
}

void HypTPCTrack::AddHelixTrack(int pdg, TPCLocalTrackHelix *tp){

  _hitClusterArray -> Delete();
  genfit::TrackCand trackCand;

  int nMeasurement = tp -> GetNHit();
  for(int i=0; i<nMeasurement; i++){
    TPCLTrackHit *point = tp -> GetHitInOrder(i);
    if(!point) continue;
    new ((*_hitClusterArray)[i]) HypTPCHit(*point);
    trackCand.addHit(TPCDetID, i);
  }

  //GenFit Units : GeV/c, ns, cm, kGauss
  //K1.8Ana Units : GeV/c, ns, mm, T
  Double_t charge = tp -> GetCharge();
  TVector3 posSeed; TVector3 momSeed;
  posSeed = tp -> GetHitInOrder(0) -> GetLocalCalPosHelix();
  posSeed.SetMag(posSeed.Mag()/10.); //mm -> cm
  momSeed = tp -> GetHitInOrder(0) -> GetMomentumHelix(charge); //GeV/c

  double resX = 0; double resY = 0; double resZ = 0;
  for(int i=0; i<nMeasurement; i++){
    auto point = tp -> GetHitInOrder(i);
    const TVector3& res_vect = point -> GetResolutionVect();
    if(res_vect.x() > 0.9e+10 && res_vect.y() > 0.9e+10 && res_vect.z() > 0.9e+10) continue;
    resX = 0.1*res_vect.X(); //mm -> cm
    resY = 0.1*res_vect.Y(); //mm -> cm
    resZ = 0.1*res_vect.Z(); //mm -> cm
    //break;
  }

  TMatrixDSym covSeed(6);
  covSeed.Zero();
  covSeed(0, 0) = resX*resX;
  covSeed(1, 1) = resY*resY;
  covSeed(2, 2) = resZ*resZ;

  double resT = 0.01*TMath::Sqrt(resX*resX + resZ*resZ); //cm -> m
  double L = 0.001*tp -> GetTransversePath(); //mm -> m
  double Pt2 = momSeed.x()*momSeed.x() + momSeed.z()*momSeed.z();
  double dPt2 = 720./(nMeasurement+4.)*Pt2*Pt2/(0.09*L*L*L*L);
  covSeed(3, 3) = dPt2*resT*resT/2.;
  covSeed(4, 4) = dPt2*resT*resT/2.;
  covSeed(5, 5) = dPt2*resT*resT/2.;
  //covSeed.Print();

  // set start values and pdg to cand
  trackCand.setCovSeed(covSeed);
  trackCand.setPosMomSeedAndPdgCode(posSeed, momSeed, pdg);
  trackCand.setTimeSeed(0.); //set defualt _time=0.;
  //trackCand.setPosMomSeed(posSeed, momSeed, helixTrack -> Charge());

  new ((*_genfitTrackArray)[_genfitTrackArray -> GetEntriesFast()]) genfit::Track(trackCand, *_measurementFactory, new genfit::RKTrackRep(pdg));
}

void HypTPCTrack::AddHelixTrack(std::vector<int> pdg, TPCLocalTrackHelix *tp){

  int Ntracks = GetNTrack();
  int Nreps = pdg.size();
  for(int i=0;i<Nreps;i++){
    int pid = pdg.at(i);
    if(i==0) AddHelixTrack(pid, tp);
    else AddReps(Ntracks, pid);
  }
}

void HypTPCTrack::AddLinearTrack(int pdg, TPCLocalTrack *tp, double momentum){

  _hitClusterArray -> Delete();
  genfit::TrackCand trackCand;

  int nMeasurement = tp -> GetNHit();
  for(int i=0; i<nMeasurement; i++){
    TPCLTrackHit *point = tp -> GetHit(i);
    if(!point) continue;
    new ((*_hitClusterArray)[i]) HypTPCHit(*point);
    trackCand.addHit(TPCDetID, i);
  }

  //GenFit Units : GeV/c, ns, cm, kGauss
  //K1.8Ana Units : GeV/c, ns, mm, T
  const TVector3& res_vect = tp -> GetHit(0) -> GetResolutionVect();
  TVector3 posSeed = tp -> GetHit(0) -> GetLocalCalPos();
  posSeed.SetMag(posSeed.Mag()/10.); //mm -> cm
  double u = tp -> GetU0(); double v = tp -> GetV0();
  TVector3 momSeed(u,v,1.0);
  momSeed.SetMag(momentum);

  TMatrixDSym covSeed(6);
  covSeed.Zero();
  double resX = 0.1*res_vect.X(); //mm -> cm
  double resY = 0.1*res_vect.Y(); //mm -> cm
  double resZ = 0.1*res_vect.Z(); //mm -> cm
  covSeed(0, 0) = resX*resX;
  covSeed(1, 1) = resY*resY;
  covSeed(2, 2) = resZ*resZ;

  double dPt2 = TMath::Power(momSeed.Mag()*0.034,2); //dp/p ~ 3.4%
  covSeed(3, 3) = dPt2/2.;
  covSeed(4, 4) = dPt2/2.;
  covSeed(5, 5) = dPt2/2.;

  // set start values and pdg to cand
  trackCand.setCovSeed(covSeed);
  trackCand.setPosMomSeedAndPdgCode(posSeed, momSeed, pdg);
  trackCand.setTimeSeed(0.); //set defualt _time=0.;
  //trackCand.setPosMomSeed(posSeed, momSeed, helixTrack -> Charge());
  new ((*_genfitTrackArray)[_genfitTrackArray -> GetEntriesFast()]) genfit::Track(trackCand, *_measurementFactory, new genfit::RKTrackRep(pdg));
}

void HypTPCTrack::AddLinearTrack(std::vector<int> pdg, TPCLocalTrack *tp, double momentum){

  int Ntracks = GetNTrack();
  int Nreps = pdg.size();
  for(int i=0;i<Nreps;i++){
    int pid = pdg.at(i);
    if(i==0) AddLinearTrack(pid, tp, momentum);
    else AddReps(Ntracks, pid);
  }
}

void HypTPCTrack::AddReconstructedTrack(int pdg, TVector3 posSeed, TVector3 momSeed){

  _hitClusterArray -> Delete();
  genfit::TrackCand trackCand;
  TVector3 resVect(0.0001, 0.0001, 0.0001);
  TPCHit *hit = new TPCHit(0, 0);
  hit -> AddHit(0, 0);

  TPCLTrackHit *lhit = new TPCLTrackHit(hit);
  lhit -> SetLocalHitPos(posSeed);
  lhit -> SetCalPosition(posSeed);
  lhit -> SetResolution(resVect);
  new ((*_hitClusterArray)[0]) HypTPCHit(*lhit);
  trackCand.addHit(TPCDetID, 0);

  TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
  track->AddTPCHit(lhit);

  TMatrixDSym covSeed(6);
  covSeed.Zero();
  covSeed(0, 0) = 0.0001;
  covSeed(1, 1) = 0.0001;
  covSeed(2, 2) = 0.0001;
  covSeed(3, 3) = 0.0001;
  covSeed(4, 4) = 0.0001;
  covSeed(5, 5) = 0.0001;
  //covSeed.Print();

  // set start values and pdg to cand
  posSeed.SetMag(posSeed.Mag()/10.); //mm -> cm
  trackCand.setCovSeed(covSeed);
  trackCand.setPosMomSeedAndPdgCode(posSeed, momSeed, pdg);
  trackCand.setTimeSeed(0.); //set defualt _time=0.;

  new ((*_genfitTrackArray)[_genfitTrackArray -> GetEntriesFast()]) genfit::Track(trackCand, *_measurementFactory, new genfit::RKTrackRep(pdg));
  delete track;
  delete hit;
}
