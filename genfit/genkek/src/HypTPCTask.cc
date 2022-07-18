
//  Authors: Wooseung Jung

//GenKEK
#include "HypTPCTask.hh"

//GenFit
#include <Exception.h>
#include <TrackPoint.h>
#include <KalmanFitterInfo.h>
#include <RKTrackRep.h>

//ROOT
#include <TMath.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>

//k18-analyzer
#include <TPCPadHelper.hh>

//STL
#include <cstddef>
#include <iostream>
#include <string>

ClassImp(HypTPCTask)

#define LogWARNING(exp) std::cout<<"WARNING: "<< __FILE__<<"::"<<__func__<<" "<<exp<<std::endl

namespace{
  const double qnan = TMath::QuietNaN();
  //const TVector3 targetsize(30,20,20);
  const TVector3 targetsize(35,25,25);
}

//GenFit Units : GeV/c, ns, cm, kGauss
//K1.8Ana Units : GeV/c, ns, mm, T

HypTPCTask& HypTPCTask::GetInstance(){

  static HypTPCTask s_instance;
  return s_instance;
}

genfit::Track* HypTPCTask::GetFittedTrack(int trackid) const{

  genfit::Track* fittedTrack = nullptr;
  if(FitCheck(trackid)) fittedTrack = GetTrack(trackid);
  return fittedTrack;
}

genfit::AbsTrackRep* HypTPCTask::GetTrackRep(int trackid) const{

  genfit::AbsTrackRep* rep = nullptr;
  genfit::Track* fittedTrack = GetFittedTrack(trackid);
  if(fittedTrack) rep = fittedTrack -> getCardinalRep();
  return rep;
}

genfit::FitStatus* HypTPCTask::GetFitStatus(int trackid) const{

  genfit::FitStatus *fitStatus = nullptr;
  genfit::Track* fittedTrack = GetFittedTrack(trackid);
  genfit::AbsTrackRep* rep = GetTrackRep(trackid);
  if(rep) fitStatus = fittedTrack -> getFitStatus(rep);
  return fitStatus;
}

genfit::MeasuredStateOnPlane HypTPCTask::GetFitState(int trackid) const{

  genfit::MeasuredStateOnPlane fitState;
  genfit::Track* fittedTrack = GetFittedTrack(trackid);
  if(fittedTrack) fitState = fittedTrack -> getFittedState();
  return fitState;
}

double HypTPCTask::GetChi2(int trackid) const{

  double chi2 = qnan;
  genfit::FitStatus *fitStatus = GetFitStatus(trackid);
  if(fitStatus) chi2 = fitStatus -> getChi2();
  return chi2;
}

double HypTPCTask::GetNDF(int trackid) const{

  double ndf = qnan;
  genfit::FitStatus *fitStatus = GetFitStatus(trackid);
  if(fitStatus) ndf = fitStatus -> getNdf();
  return ndf;
}

double HypTPCTask::GetChi2NDF(int trackid) const{

  double chi2ndf = qnan;
  double chi2 = GetChi2(trackid);
  double ndf = GetNDF(trackid);
  if(!TMath::IsNaN(chi2) && !TMath::IsNaN(ndf)) chi2ndf = chi2/ndf;
  return chi2ndf;
}

double HypTPCTask::GetCharge(int trackid) const{

  double charge = qnan;
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid);
  if(FitCheck(trackid)) charge = fitState.getCharge();
  return charge;
}

TVector3 HypTPCTask::GetMom(int trackid) const{

  TVector3 mom(qnan, qnan, qnan);
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid);
  if(FitCheck(trackid)) mom = fitState.getMom();
  return mom;
}

TVector3 HypTPCTask::GetPos0(int trackid) const{

  TVector3 pos(qnan, qnan, qnan);
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid);
  if(FitCheck(trackid)) pos = 10.*fitState.getPos(); //cm -> mm
  return pos;
}

//GenFit Units : GeV/c, ns, cm, kGauss
//K1.8Ana Units : GeV/c, ns, mm, T
double HypTPCTask::GetTrackLength(int trackid, int start, int end) const{

  double length;
  genfit::Track* fittedTrack = GetFittedTrack(trackid);
  try{length = 10.*fittedTrack -> getTrackLen(nullptr,start,end);} //cm -> mm
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return length=qnan;
  }
  return length;
}

double HypTPCTask::GetTrackTOF(int trackid, int start, int end) const{

  double TOF = qnan;
  genfit::Track* fittedTrack = GetFittedTrack(trackid);
  if(end==-1) end = -2;
  if(fittedTrack) TOF = fittedTrack -> getTOF(nullptr,start,end);
  return TOF;
}

bool HypTPCTask::ExtrapolateTrack(int trackid, double distance, TVector3 &pos) const{

  Double_t d = 0.1*distance; //mm -> cm
  if(distance>=0.) d += 0.1*GetTrackLength(trackid);
  genfit::RKTrackRep *rep = (genfit::RKTrackRep *) GetTrackRep(trackid);
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid);
  if(!rep) return false;
  try{rep -> extrapolateBy(fitState, d);}
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return false;
  }

  pos = 10.*fitState.getPos(); //cm -> mm
  return true;
}

bool HypTPCTask::ExtrapolateToPoint(int trackid, TVector3 point, TVector3 &pos) const{

  genfit::RKTrackRep *rep = (genfit::RKTrackRep *) GetTrackRep(trackid);
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid);
  if(!rep) return false;
  try{rep -> extrapolateToPoint(fitState, 0.1*point);} //mm -> cm
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return false;
  }

  pos = 10.*fitState.getPos(); //cm -> mm
  return true;
}

bool HypTPCTask::GetPosOnPlane(int trackid, genfit::SharedPlanePtr plane, TVector3 &pos) const{

  genfit::RKTrackRep *rep = (genfit::RKTrackRep *) GetTrackRep(trackid);
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid);
  if(!rep) return false;
  try{rep -> extrapolateToPlane(fitState, plane);}
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return false;
  }

  pos = 10.*fitState.getPos(); //cm -> mm
  return true;
}

bool HypTPCTask::IsInsideTarget(int trackid) const{

  double ztgt = tpc::ZTarget;
  TVector3 pos; TVector3 target(0.,0.,ztgt);
  if(ExtrapolateToPoint(trackid, target, pos) &&
     (TMath::Abs(pos.x()) < targetsize.x()/2.0) &&
     (TMath::Abs(pos.y()) < targetsize.y()/2.0) &&
     (TMath::Abs(pos.z()-ztgt) < targetsize.z()/2.0)) return true;
  else return false;
}
