#include "TPCTrack.hh"

namespace {
  HypTPCTask &fHypTPCTask = HypTPCTask::GetInstance();
}

TPCTrack::TPCTrack() {
  Initialize();
}

TPCTrack::TPCTrack(TPCTrack *track) {

  Initialize();

  fMomentumTarget = track->GetMomentumTarget();
  fMomentumHTOF = track->GetMomentumHTOF();
  
  fPosAtTarget = track->GetPosAtTarget();
  fPosAtHTOF = track->GetPosAtHTOF();

  fCharge = track->GetCharge();
  fGFCharge = track->GetGFCharge();

  fTrackLength = track->GetTrackLength();
  fTrackLengthHTOF = track->GetTrackLengthHTOF();

  fHelixID = track->GetHelixTrackID();
  fHelixTrack = track->GetHelixTrack();
  
  fGenFitID = track->GetGFTrackID();
  fGFTrack = track->GetGenFitTrack();

  fPID = track->GetPID();
  fPIDProb = track->GetPIDProbability();

  fNumCluster = track->GetNumCluster();

  fStatus = track->GetStatus();

  // For RAVE 
  fMomentum = track->GetMomentum();
  fVertex = track->GetVertex();
  fVertexID = track->GetVertexID();

  fChisqr = track->GetChisqr();
  fPval = track->GetPvalue();
  fTOF = track->GetTrackTOF();
  fTOFHTOF = track->GetTrackTOFHTOF();

  fVertexWeight = track->GetVertexWeight();
  
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
  fTrackLengthHTOF = -9999;
  
  fHelixID = -9999;
  fGenFitID = -9999;
  fVertexID = -9999;

  fPID = TPCPID::kNotDetermined;
  fPIDProb = 0.;

  fNumCluster = 0;

  fChisqr = -9999;
  fPval = -9999;
  fTOF = -9999;
  fTOFHTOF = -9999;

  fVertexWeight = -9999;
  
  fHelixTrack = NULL;
  fGFTrack = NULL;
}


  
  
  
