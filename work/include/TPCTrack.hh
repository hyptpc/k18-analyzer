#ifndef TPC_TRACK_HH
#define TPC_TRACK_HH

#include "TVector3.h"

#include "TPCPID.hh"
#include "TPCRiemannTrack.hh"

#include "HypTPCTask.hh"

class HypTPCTask;

class TPCTrack 
{
public:

  TPCTrack();
  TPCTrack(TPCTrack *track);
  virtual ~TPCTrack() {}

  virtual void Initialize();

  void SetPID(TPCPID::PID pdg) {fPID = pdg;}
  TPCPID::PID GetPID() const {return fPID;}
  
  void SetPIDProbability(Double_t val)  {fPIDProb = val;}
  Double_t GetPIDProbability() const {return fPIDProb;}

  void SetCharge(Int_t q) { fCharge = q;}
  Int_t GetCharge() const {return fCharge;}

  void SetGFCharge(Int_t q) { fGFCharge = q;}
  Int_t GetGFCharge() const {return fGFCharge;}

  void SetHelixTrackID(Int_t id) {fHelixID = id;}
  Int_t GetHelixTrackID() const { return fHelixID;}

  void SetGFTrackID(Int_t id) {fGenFitID = id;}
  Int_t GetGFTrackID() const { return fGenFitID;}

  void SetVertexID(Int_t id) {fVertexID = id;}
  Int_t GetVertexID() const {return fVertexID;}

  void SetMomentum(TVector3 mom) {fMomentum = mom;}
  TVector3 GetMomentum() const {return fMomentum;}

  void SetMomentumTarget(TVector3 mom) {fMomentumTarget = mom;}
  TVector3 GetMomentumTarget() const {return fMomentumTarget;}
  
  void SetPosAtTarget(TVector3 pos) {fPosAtTarget = pos;}
  TVector3 GetPosAtTarget() const {return fPosAtTarget;}

  void SetMomentumHTOF(TVector3 mom) {fMomentumHTOF = mom;}
  TVector3 GetMomentumHTOF() const {return fMomentumHTOF;}  
  
  void SetPosAtHTOF(TVector3 pos) {fPosAtHTOF = pos;}
  TVector3 GetPosAtHTOF() const {return fPosAtHTOF;}

  void SetVertex(TVector3 pos) {fVertex = pos;}
  TVector3 GetVertex() const {return fVertex;}  

  void SetVertexWeight(Double_t weight) {fVertexWeight = weight;}
  Double_t GetVertexWeight() const {return fVertexWeight;}
  
  void SetGenFitTrack(genfit::Track *track) {fGFTrack = track;}
  genfit::Track *GetGenFitTrack() {return fGFTrack;}

  void SetHelixTrack(TPCRiemannTrack *track) {fHelixTrack = track;}
  TPCRiemannTrack *GetHelixTrack(){return fHelixTrack;}

  void SetStatus(Int_t n) {fStatus = n;}
  Int_t GetStatus() const {return fStatus;}

  void SetNumCluster(Int_t n) {fNumCluster = n;}
  Int_t GetNumCluster() const {return fNumCluster;}

  void SetTrackLength(Double_t len) {fTrackLength = len;}
  Double_t GetTrackLength() const {return fTrackLength;}

  void SetTrackLengthHTOF(Double_t len) {fTrackLengthHTOF = len;}
  Double_t GetTrackLengthHTOF() const {return fTrackLengthHTOF;}

  void SetChisqr(Double_t chisqr) {fChisqr= chisqr;}
  Double_t GetChisqr() const {return fChisqr;}

  void SetPvalue(Double_t pval) {fPval = pval;}
  Double_t GetPvalue() const {return fPval;}

  void SetTrackTOF(Double_t tof) {fTOF = tof;}
  Double_t GetTrackTOF() const {return fTOF;}
  
  void SetTrackTOFHTOF(Double_t tof) {fTOFHTOF = tof;}
  Double_t GetTrackTOFHTOF() const {return fTOFHTOF;}
  
private:

  TPCRiemannTrack *fHelixTrack;
  genfit::Track *fGFTrack;
  
  TVector3 fMomentum;
  TVector3 fMomentumTarget;
  TVector3 fMomentumHTOF;

  TVector3 fPosAtTarget;
  TVector3 fPosAtHTOF;
  TVector3 fVertex;

  Int_t fCharge;
  Int_t fGFCharge;

  Double_t fTrackLength;
  Double_t fTrackLengthHTOF;
  
  Int_t fHelixID;
  Int_t fGenFitID;
  Int_t fVertexID;

  TPCPID::PID fPID;
  Double_t fPIDProb;

  Int_t fNumCluster;

  Int_t fStatus;
  
  Double_t fChisqr;
  Double_t fPval;
  Double_t fTOF;
  Double_t fTOFHTOF;

  Double_t fVertexWeight;
};


#endif
