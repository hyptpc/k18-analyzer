
//  Authors: Wooseung Jung

//GenKEK
#include "HypTPCTask.hh"

//GenFit
#include <Exception.h>
#include <TrackPoint.h>
#include <KalmanFitterInfo.h>
#include <RKTrackRep.h>
#include <KalmanFitStatus.h>

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
  const double ztgt = tpc::ZTarget;
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
  if(TrackCheck(trackid)) fittedTrack = GetTrack(trackid);
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

genfit::MeasuredStateOnPlane HypTPCTask::GetFitState(int trackid, int pointid) const{

  genfit::MeasuredStateOnPlane fitState = 0;
  genfit::Track* fittedTrack = GetFittedTrack(trackid);
  if(fittedTrack) fitState = fittedTrack -> getFittedState(pointid);
  return fitState;
}

bool HypTPCTask::GetTrackPull(int trackid, int pdg, TVector3 res_vect, double tracklen, TVector3 g4mom, TVector3 g4pos, double *residual, double *pull, double *residual6D, double *pull6D) const{

  for(int i=0;i<6;i++){
    if(i!=5){
      residual[i] = qnan;
      pull[i] = qnan;
    }
    residual6D[i] = qnan;
    pull6D[i] = qnan;
  }
  //ideal state(generated state)
  genfit::AbsTrackRep* g4rep = new genfit::RKTrackRep(pdg);
  genfit::MeasuredStateOnPlane g4state(g4rep);
  TMatrixDSym covM(6);
  covM.Zero();
  g4rep -> setPosMomCov(g4state, 0.1*g4pos, g4mom, covM); //mm -> cm
  genfit::StateOnPlane g4stateRef(g4state);
  TVectorD &referenceState = g4stateRef.getState();
  TVectorD referenceState6D = g4stateRef.get6DState();
  //fitted track's state
  if(!TrackCheck(trackid)) return false;
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid);
  genfit::AbsTrackRep* rep = GetTrackRep(trackid);
  try{rep -> extrapolateToPlane(fitState, g4stateRef.getPlane());} //mm -> cm
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return false;
  }

  TVectorD &state = fitState.getState();
  TMatrixDSym &cov = fitState.getCov();
  TVectorD state6D = fitState.get6DState();
  TMatrixDSym cov6D = fitState.get6DCov();

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

  residual6D[0] = state6D[0]-referenceState6D[0];
  residual6D[1] = state6D[1]-referenceState6D[1];
  residual6D[2] = state6D[2]-referenceState6D[2];
  residual6D[3] = state6D[3]-referenceState6D[3];
  residual6D[4] = state6D[4]-referenceState6D[4];
  residual6D[5] = state6D[5]-referenceState6D[5];

  pull6D[0] = (state6D[0]-referenceState6D[0]) / sqrt(cov6D[0][0]);
  pull6D[1] = (state6D[1]-referenceState6D[1]) / sqrt(cov6D[1][1]);
  pull6D[2] = (state6D[2]-referenceState6D[2]) / sqrt(cov6D[2][2]);
  pull6D[3] = (state6D[3]-referenceState6D[3]) / sqrt(cov6D[3][3]);
  pull6D[4] = (state6D[4]-referenceState6D[4]) / sqrt(cov6D[4][4]);
  pull6D[5] = (state6D[5]-referenceState6D[5]) / sqrt(cov6D[5][5]);

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

int HypTPCTask::GetNumOfIterations(int trackid) const{

  int itr = 0;
  genfit::KalmanFitStatus *fitStatus = (genfit::KalmanFitStatus *) GetFitStatus(trackid);
  if(fitStatus) itr = fitStatus -> getNumIterations();
  return itr;
}

double HypTPCTask::GetCharge(int trackid) const{

  double charge = qnan;
  //genfit::AbsTrackRep* rep = GetTrackRep(trackid);
  //if(rep) charge = rep -> getPDGCharge();
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid);
  if(TrackCheck(trackid)) charge = fitState.getCharge();
  return charge;
}

int HypTPCTask::GetPDGcode(int trackid) const{

  int code = 0;
  genfit::AbsTrackRep* rep = GetTrackRep(trackid);
  if(rep) code = rep -> getPDG();
  return code;
}

TVector3 HypTPCTask::GetMom(int trackid, int pointid) const{

  TVector3 mom(qnan, qnan, qnan);
  if(GetNHits(trackid)<=pointid){
    std::cout<<"Error : # of hits <= point Id"<<std::endl;
    return mom;
  }
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid, pointid);
  if(TrackCheck(trackid)) mom = fitState.getMom();
  return mom;
}

TVector3 HypTPCTask::GetPos(int trackid, int pointid) const{

  TVector3 pos(qnan, qnan, qnan);
  if(GetNHits(trackid)<=pointid){
    std::cout<<"Error : # of hits <= point Id"<<std::endl;
    return pos;
  }
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid, pointid);
  if(TrackCheck(trackid)) pos = 10.*fitState.getPos(); //cm -> mm
  return pos;
}

int HypTPCTask::GetNHits(int trackid) const{

  int nhits;
  genfit::Track* fittedTrack = GetFittedTrack(trackid);
  try{nhits = fittedTrack -> getNumPointsWithMeasurement(); }
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return 0;
  }
  return nhits;
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

bool HypTPCTask::ExtrapolateTrack(int trackid, double distance, TVector3 &pos, TVector3 &mom) const{

  double tracklength = GetTrackLength(trackid);
  if(TMath::IsNaN(tracklength)) return false;
  Double_t d = 0.1*distance; //mm -> cm
  if(distance>=0.) d += 0.1*tracklength;
  genfit::RKTrackRep *rep = (genfit::RKTrackRep *) GetTrackRep(trackid);
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid);
  if(!rep || !TrackCheck(trackid)) return false;
  try{rep -> extrapolateBy(fitState, d);}
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return false;
  }

  pos = 10.*fitState.getPos(); //cm -> mm
  mom = fitState.getMom(); //cm -> mm
  return true;
}

bool HypTPCTask::ExtrapolateToPoint(int trackid, TVector3 point, TVector3 &pos, TVector3 &mom, double &tracklen, double &tof) const{

  double tracklength = GetTrackLength(trackid);
  if(TMath::IsNaN(tracklength)) return false;
  tracklen = 0.1*tracklength; //mm -> cm
  genfit::RKTrackRep *rep = (genfit::RKTrackRep *) GetTrackRep(trackid);
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid);
  if(!rep || !TrackCheck(trackid)) return false;

  double time0 = fitState.getTime();
  double len;
  try{len = rep -> extrapolateToPoint(fitState, 0.1*point);} //mm -> cm
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return false;
  }
  if(len >= 0.) tracklen += len;
  else tracklen = len;
  pos = 10.*fitState.getPos(); //cm -> mm
  mom = fitState.getMom(); //cm -> mm
  tof = fitState.getTime() - time0;
  tracklen *= 10.; //cm -> mm
  return true;
}

bool HypTPCTask::ExtrapolateToPlane(int trackid, genfit::SharedPlanePtr plane, TVector3 &pos, TVector3 &mom, double &tracklen, double &tof) const{

  double tracklength = GetTrackLength(trackid);
  if(TMath::IsNaN(tracklength)) return false;
  tracklen = 0.1*tracklength; //mm -> cm
  genfit::RKTrackRep *rep = (genfit::RKTrackRep *) GetTrackRep(trackid);
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid);
  if(!rep || !TrackCheck(trackid)) return false;

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
  mom = fitState.getMom(); //cm -> mm
  tof = fitState.getTime() - time0;
  tracklen *= 10.; //cm -> mm
  return true;
}

bool HypTPCTask::ExtrapolateToTarget(int trackid, TVector3 &pos, TVector3 &mom, double &tracklen, double &tof) const{

  TVector3 target(0.,0.,ztgt);
  return ExtrapolateToPoint(trackid, target, pos, mom, tracklen, tof);
}

bool HypTPCTask::IsInsideTarget(int trackid) const{

  TVector3 pos; TVector3 mom; double tracklen; double tof;
  if(ExtrapolateToTarget(trackid, pos, mom, tracklen, tof) &&
     (TMath::Abs(pos.x()) < targetsize.x()/2.0) &&
     (TMath::Abs(pos.y()) < targetsize.y()/2.0) &&
     (TMath::Abs(pos.z()-ztgt) < targetsize.z()/2.0)) return true;
  else return false;
}

bool HypTPCTask::ExtrapolateToHTOF(int trackid, TVector3 &pos, TVector3 &mom, double &tracklen, double &tof) const{

  bool flag = false;
  for(int i=0;i<8;i++){
    const TVector3 center = 10*HTOFPlane[i] -> getO(); //cm -> mm
    if(ExtrapolateToPlane(trackid, HTOFPlane[i], pos, mom, tracklen, tof)){
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
