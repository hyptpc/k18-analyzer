
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
  const double htof_l = 34.86; //center to HTOF downstream
}

//GenFit Units : GeV/c, ns, cm, kGauss
//K1.8Ana Units : GeV/c, ns, mm, T

HypTPCTask::HypTPCTask() : HypTPCFitProcess() {

  TVector3 pointRef(0,0,htof_l);
  TVector3 normalRef(0,0,1.);

  for(int i=0; i<8; i++){
    if(i!=0){
      double angle = 0.25*TMath::Pi();
      pointRef.RotateY(angle);
      normalRef.RotateY(angle);
    }
    HTOFPlane[i] = genfit::SharedPlanePtr(new genfit::DetPlane(pointRef, normalRef));
  }
}

HypTPCTask::~HypTPCTask(){}

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

bool HypTPCTask::GetTrackPull(int trackid, int pdg, TVector3 g4mom, TVector3 g4pos, double *residual, double *pull) const{

  for(int i=0;i<5;i++){
    residual[i] = qnan;
    pull[i] = qnan;
  }
  //ideal state(generated state)
  genfit::AbsTrackRep* g4rep = new genfit::RKTrackRep(pdg);
  genfit::MeasuredStateOnPlane g4state(g4rep);
  TMatrixDSym covM(6);
  covM.Zero();
  g4rep -> setPosMomCov(g4state, g4pos, g4mom, covM);
  genfit::StateOnPlane g4stateRef(g4state);
  const TVectorD &referenceState = g4stateRef.getState();

  //fitted track's state
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid);
  if(!FitCheck(trackid)) return false;
  TVectorD &state = fitState.getState();
  TMatrixDSym &cov = fitState.getCov();

  residual[0] = GetCharge(trackid)/state[0]-g4mom.Mag();
  residual[1] = state[1]-referenceState[1];
  residual[2] = state[2]-referenceState[2];
  residual[3] = state[3]-referenceState[3];
  residual[4] = state[4]-referenceState[4];

  pull[0] = (state[0]-referenceState[0]) / sqrt(cov[0][0]);
  pull[1] = (state[1]-referenceState[1]) / sqrt(cov[1][1]);
  pull[2] = (state[2]-referenceState[2]) / sqrt(cov[2][2]);
  pull[3] = (state[3]-referenceState[3]) / sqrt(cov[3][3]);
  pull[4] = (state[4]-referenceState[4]) / sqrt(cov[4][4]);
  return true;
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

double HypTPCTask::GetPvalue(int trackid) const{

  double pval = qnan;
  genfit::FitStatus *fitStatus = GetFitStatus(trackid);
  if(fitStatus) pval = fitStatus -> getPVal();
  return pval;
}

double HypTPCTask::GetCharge(int trackid) const{

  double charge = qnan;
  genfit::AbsTrackRep* rep = GetTrackRep(trackid);
  if(rep) charge = rep -> getPDGCharge();
  return charge;
}

int HypTPCTask::GetPDGcode(int trackid) const{

  int code = 0;
  genfit::AbsTrackRep* rep = GetTrackRep(trackid);
  if(rep) code = rep -> getPDG();
  return code;
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

double HypTPCTask::GetTrackLength(int trackid, int start, int end) const{

  double length;
  genfit::Track* fittedTrack = GetFittedTrack(trackid);
  try{length = 10.*fittedTrack -> getTrackLen(nullptr,start,end);} //cm -> mm
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return qnan;
  }
  return length;
}

double HypTPCTask::GetTrackTOF(int trackid, int start, int end) const{

  double TOF;
  genfit::Track* fittedTrack = GetFittedTrack(trackid);
  if(end==-1) end = -2;
  try{TOF = fittedTrack -> getTOF(nullptr,start,end);}
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return qnan;
  }
  return TOF;
}

bool HypTPCTask::ExtrapolateTrack(int trackid, double distance, TVector3 &pos) const{

  double tracklength = GetTrackLength(trackid);
  if(TMath::IsNaN(tracklength)) return false;
  Double_t d = 0.1*distance; //mm -> cm
  if(distance>=0.) d += 0.1*tracklength;
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

bool HypTPCTask::ExtrapolateToPlane(int trackid, genfit::SharedPlanePtr plane, TVector3 &pos, double &tracklen, double &tof) const{

  double tracklength = GetTrackLength(trackid);
  if(TMath::IsNaN(tracklength)) return false;
  tracklen = 0.1*tracklength; //mm -> cm
  genfit::RKTrackRep *rep = (genfit::RKTrackRep *) GetTrackRep(trackid);
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid);
  if(!rep) return false;

  double time0 = fitState.getTime();
  try{tracklen += rep -> extrapolateToPlane(fitState, plane);}
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return false;
  }

#if 0
  std::vector<genfit::MatStep> steps = rep->getSteps();
  std::cout<<"RK number of steps : "<<steps.size()<<std::endl;
  for(int i=0;i<steps.size();i++){
    std::cout<<"Step : "<<i<<" Size of step : "<<steps.at(i).stepSize_<<std::endl;
    std::cout<<"Density [g/cm3] : "<<steps.at(i).material_.density<<std::endl;
    std::cout<<"Radiation length [cm] : "<<steps.at(i).material_.radiationLength<<std::endl;
    std::cout<<"Atomic number Z : "<<steps.at(i).material_.Z<<std::endl;
    std::cout<<"Mass number A : "<<steps.at(i).material_.A<<std::endl;
    std::cout<<"mean excitation energy [eV] : "<<steps.at(i).material_.mEE<<std::endl;
  }
#endif

  pos = 10.*fitState.getPos(); //cm -> mm
  tof = fitState.getTime() - time0;
  tracklen *= 10.; //cm -> mm
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

bool HypTPCTask::ExtrapolateToHTOF(int trackid, TVector3 &pos, double &tracklen, double &tof) const{

  bool flag = false;
  for(int i=0;i<8;i++){
    const TVector3 center = 10*HTOFPlane[i] -> getO(); //cm -> mm
    if(ExtrapolateToPlane(trackid, HTOFPlane[i], pos, tracklen, tof)){
      if(tracklen<0||tof<0) continue;
      double xzdist = TMath::Sqrt((center.x()-pos.x())*(center.x()-pos.x()) +
				  (center.z()-pos.z())*(center.z()-pos.z()));
      if(TMath::Abs(pos.y()-12)>400. || xzdist>142) continue;
      if(i==4 && TMath::Abs(pos.x())<71 && pos.y()<62 && pos.y()>-50) break; // window
      else flag=true;
      break;
    }
  }
  return flag;
}
