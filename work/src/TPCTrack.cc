#include "TPCTrack.hh"

namespace {
  HypTPCTask &fHypTPCTask = HypTPCTask::GetInstance();
}

TPCTrack::TPCTrack() {
  Initialize();
}

TPCTrack::TPCTrack(TPCTrack *track) {

  Initialize();

  fMomentum = track->GetMomentum();
  fMomentumTarget = track->GetMomentumTarget();
  fMomentumHTOF = track->GetMomentumHTOF();
  
  fPosAtTarget = track->GetPosAtTarget();
  fPosAtHTOF = track->GetPosAtHTOF();
  fVertex = track->GetVertex();

  fCharge = track->GetCharge();
  fGFCharge = track->GetGFCharge();

  fTrackLength = track->GetTrackLength();

  fHelixID = track->GetHelixTrackID();
  fHelixTrack = track->GetHelixTrack();
  
  fGenFitID = track->GetGFTrackID();
  fGFTrack = track->GetGenFitTrack();

  fPID = track->GetPID();
  fPIDProb = track->GetPIDProbability();

  fNumCluster = track->GetNumCluster();
}
  
void TPCTrack::Initialize() {

  fMomentum.SetXYZ(-9999,-9999,-9999);
  fMomentumTarget.SetXYZ(-9999,-9999,-9999);
  fMomentumHTOF.SetXYZ(-9999,-9999,-9999);

  fPosAtTarget.SetXYZ(-999,-999,-999);
  fPosAtHTOF.SetXYZ(-999,-999,-999);
  fVertex.SetXYZ(-999,-999,-999);

  fCharge = 0;
  fGFCharge = 0;

  fTrackLength = -9999;
  
  fHelixID = -9999;
  fGenFitID = -9999;
  fVertexID = -9999;

  fPID = TPCPID::kNotDetermined;
  fPIDProb = 0.;

  fNumCluster = 0;
}


  
  
  
