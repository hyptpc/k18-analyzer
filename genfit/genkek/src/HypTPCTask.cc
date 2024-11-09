
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
#include <DeleteUtility.hh>

//STL
#include <cstddef>
#include <iostream>
#include <string>

ClassImp(HypTPCTask)

#define LogWARNING(exp) std::cout<<"WARNING: "<< __FILE__<<"::"<<__func__<<" "<<exp<<std::endl

namespace{
  const double qnan = TMath::QuietNaN();
  const double htof_l = 34.86; //center to HTOF downstream
  const double ztgt = tpc::ZTarget;
  const TVector3 tgtcenter(0,0,0.1*tpc::ZTarget); //cm
}

//GenFit Units : GeV/c, ns, cm, kGauss
//K1.8Ana Units : GeV/c, ns, mm, T

HypTPCTask::HypTPCTask() : HypTPCFitProcess() {

  TVector3 pointRef(0,0,-htof_l);
  TVector3 normalRef(0,0,-1.);

  for(int i=0; i<8; i++){
    if(i!=0){
      double angle = 0.25*TMath::Pi();
      pointRef.RotateY(angle);
      normalRef.RotateY(angle);
    }
    HTOFPlane[i] = genfit::SharedPlanePtr(new genfit::DetPlane(pointRef, normalRef));
  }

  TVector3 tgtnormal(0,0,1.);
  TgtPlane = genfit::SharedPlanePtr(new genfit::DetPlane(tgtcenter, tgtnormal));
}

HypTPCTask::~HypTPCTask(){
  Clear(); //Clear the track container
}

HypTPCTask& HypTPCTask::GetInstance(){

  static HypTPCTask s_instance;
  s_instance.Clear();
  return s_instance;
}

genfit::Track* HypTPCTask::GetFittedTrack(int trackid) const{

  genfit::Track* fittedTrack = nullptr;
  if(TrackCheck(trackid)) fittedTrack = GetTrack(trackid);
  return fittedTrack;
}

genfit::AbsTrackRep* HypTPCTask::GetTrackRep(int trackid, int repid) const{

  genfit::AbsTrackRep* rep = nullptr;
  genfit::Track* fittedTrack = GetFittedTrack(trackid);
  if(fittedTrack&&repid==-1) rep = fittedTrack -> getCardinalRep();
  else if(fittedTrack&&repid!=-1) rep = fittedTrack -> getTrackRep(repid);
  return rep;
}

genfit::FitStatus* HypTPCTask::GetFitStatus(int trackid, int repid) const{

  genfit::FitStatus *fitStatus = nullptr;
  genfit::Track* fittedTrack = GetFittedTrack(trackid);
  genfit::AbsTrackRep* rep = GetTrackRep(trackid, repid);
  if(rep) fitStatus = fittedTrack -> getFitStatus(rep);
  return fitStatus;
}

genfit::MeasuredStateOnPlane HypTPCTask::GetFitState(int trackid, int pointid, int repid) const{

  genfit::MeasuredStateOnPlane fitState = 0;
  genfit::Track* fittedTrack = GetFittedTrack(trackid);
  if(fittedTrack){
    genfit::AbsTrackRep* rep = GetTrackRep(trackid, repid);
    fitState = fittedTrack -> getFittedState(pointid, rep);
  }
  return fitState;
}

bool HypTPCTask::GetTrackPull(int trackid, int pdg, TVector3 res_vect, double tracklen, TVector3 g4mom, TVector3 g4pos, double *residual, double *pull, double *residual6D, double *pull6D) const{ //for Geant4

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
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid, 0);
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

double HypTPCTask::GetChi2(int trackid, int repid) const{

  double chi2 = qnan;
  genfit::FitStatus *fitStatus = GetFitStatus(trackid, repid);
  if(fitStatus) chi2 = fitStatus -> getChi2();
  return chi2;
}

double HypTPCTask::GetNDF(int trackid, int repid) const{

  double ndf = qnan;
  genfit::FitStatus *fitStatus = GetFitStatus(trackid, repid);
  if(fitStatus) ndf = fitStatus -> getNdf();
  return ndf;
}

double HypTPCTask::GetChi2NDF(int trackid, int repid) const{

  double chi2ndf = qnan;
  double chi2 = GetChi2(trackid, repid);
  double ndf = GetNDF(trackid, repid);
  if(!TMath::IsNaN(chi2) && !TMath::IsNaN(ndf)) chi2ndf = chi2/ndf;
  return chi2ndf;
}

double HypTPCTask::GetPvalue(int trackid, int repid) const{

  double pval = qnan;
  genfit::FitStatus *fitStatus = GetFitStatus(trackid, repid);
  if(fitStatus) pval = fitStatus -> getPVal();
  return pval;
}

int HypTPCTask::GetNumOfIterations(int trackid, int repid) const{

  int itr = 0;
  genfit::KalmanFitStatus *fitStatus = (genfit::KalmanFitStatus *) GetFitStatus(trackid, repid);
  if(fitStatus) itr = fitStatus -> getNumIterations();
  return itr;
}

double HypTPCTask::GetCharge(int trackid, int repid) const{

  double charge = qnan;
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid, 0, repid);
  if(TrackCheck(trackid)) charge = fitState.getCharge();
  return charge;
}

int HypTPCTask::GetPDGcode(int trackid, int repid) const{

  int code = 0;
  genfit::AbsTrackRep* rep = GetTrackRep(trackid, repid);
  if(rep) code = rep -> getPDG();
  return code;
}

TVector3 HypTPCTask::GetMom(int trackid, int pointid, int repid) const{

  TVector3 mom(qnan, qnan, qnan);
  if(GetNHits(trackid)<=pointid){
    std::cout<<"Error : # of hits <= point Id"<<std::endl;
    return mom;
  }
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid, pointid, repid);
  if(TrackCheck(trackid)) mom = fitState.getMom();
  return mom;
}

TVector3 HypTPCTask::GetPos(int trackid, int pointid, int repid) const{

  TVector3 pos(qnan, qnan, qnan);
  if(GetNHits(trackid)<=pointid){
    std::cout<<"Error : # of hits <= point Id"<<std::endl;
    return pos;
  }
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid, pointid, repid);
  if(TrackCheck(trackid)) pos = 10.*fitState.getPos(); //cm -> mm
  return pos;
}

int HypTPCTask::GetNHits(int trackid) const{

  int nhits = 0;
  if(!TrackCheck(trackid)) return nhits;
  genfit::Track* fittedTrack = GetFittedTrack(trackid);
  try{nhits = fittedTrack -> getNumPointsWithMeasurement(); }
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return 0;
  }
  return nhits;
}

double HypTPCTask::GetTrackLength(int trackid, int start, int end, int repid) const{

  double length = qnan;
  if(!TrackCheck(trackid)) return length;
  genfit::Track* fittedTrack = GetFittedTrack(trackid);
  genfit::AbsTrackRep* rep = GetTrackRep(trackid, repid);
  try{length = 10.*fittedTrack -> getTrackLen(rep, start, end);} //cm -> mm
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return qnan;
  }
  return length;
}

double HypTPCTask::GetTrackTOF(int trackid, int start, int end, int repid) const{

  double TOF = qnan;
  if(!TrackCheck(trackid, repid)) return TOF;
  genfit::Track* fittedTrack = GetFittedTrack(trackid);
  genfit::AbsTrackRep* rep = GetTrackRep(trackid, repid);
  GetFitState(trackid, 0, repid);
  try{TOF = fittedTrack -> getTOF(rep, start, end);} //ns
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return qnan;
  }
  return TOF;
}

bool HypTPCTask::ExtrapolateTrack(int trackid, double distance, TVector3 &pos, TVector3 &mom, int repid) const{

  double tracklength = GetTrackLength(trackid, 0, -1, repid);
  if(TMath::IsNaN(tracklength)) return false;
  Double_t d = 0.1*distance; //mm -> cm
  if(distance>=0.) d += 0.1*tracklength;
  genfit::RKTrackRep *rep = (genfit::RKTrackRep *) GetTrackRep(trackid, repid);
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid, 0, repid);
  if(!rep || !TrackCheck(trackid)) return false;
  try{rep -> extrapolateBy(fitState, d);}
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return false;
  }

  pos = 10.*fitState.getPos(); //cm -> mm
  mom = fitState.getMom();
  return true;
}

bool HypTPCTask::ExtrapolateToPoint(int trackid, TVector3 point, TVector3 &pos, TVector3 &mom, double &tracklen, double &tof, int repid) const{

  double tracklength = GetTrackLength(trackid, 0, -1, repid);
  if(TMath::IsNaN(tracklength)) return false;
  genfit::RKTrackRep *rep = (genfit::RKTrackRep *) GetTrackRep(trackid, repid);
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid, 0, repid);
  if(!rep || !TrackCheck(trackid)) return false;
  double time0 = fitState.getTime();

  double len;
  try{len = 10.*rep -> extrapolateToPoint(fitState, 0.1*point);} //cm -> mm
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    LogWARNING("failed!");
    std::cerr << e.what();
    return false;
  }

#if 0
  std::vector<genfit::MatStep> steps = rep->getSteps();
  std::cout<<"RK number of steps : "<<steps.size()<<std::endl;
  for(int i=0;i<steps.size();i++){
    std::cout<<"Step : "<<i<<" Size of step : "<<steps.at(i).stepSize_<<" cm"<<std::endl;
    std::cout<<"Density [g/cm3] : "<<steps.at(i).material_.density<<std::endl;
    std::cout<<"Radiation length [cm] : "<<steps.at(i).material_.radiationLength<<std::endl;
    std::cout<<"Atomic number Z : "<<steps.at(i).material_.Z<<std::endl;
    std::cout<<"Mass number A : "<<steps.at(i).material_.A<<std::endl;
    std::cout<<"mean excitation energy [eV] : "<<steps.at(i).material_.mEE<<std::endl;
  }
#endif

  pos = 10.*fitState.getPos(); //cm -> mm
  mom = fitState.getMom();
  tof = fitState.getTime() - time0;
  tracklen = len; //mm
  return true;
}

bool HypTPCTask::ExtrapolateToPlane(int trackid, genfit::SharedPlanePtr plane, TVector3 &pos, TVector3 &mom, double &tracklen, double &tof, int repid) const{

  double tracklength = GetTrackLength(trackid, 0, -1, repid);
  if(TMath::IsNaN(tracklength)) return false;
  genfit::RKTrackRep *rep = (genfit::RKTrackRep *) GetTrackRep(trackid, repid);
  genfit::MeasuredStateOnPlane fitState0 = GetFitState(trackid, 0, repid);
  genfit::MeasuredStateOnPlane fitState = GetFitState(trackid, GetNHits(trackid)-1, repid);
  if(!rep || !TrackCheck(trackid)) return false;
  double time0 = fitState0.getTime();
  double len;
  try{len = 10.*rep -> extrapolateToPlane(fitState, plane);}
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("failed!");
    if(verbosity>=1) std::cerr << e.what();
    return false;
  }

#if 0
  std::vector<genfit::MatStep> steps = rep->getSteps();
  std::cout<<"RK number of steps : "<<steps.size()<<std::endl;
  for(int i=0;i<steps.size();i++){
    std::cout<<"Step : "<<i<<" Size of step : "<<steps.at(i).stepSize_<<" cm"<<std::endl;
    std::cout<<"Density [g/cm3] : "<<steps.at(i).material_.density<<std::endl;
    std::cout<<"Radiation length [cm] : "<<steps.at(i).material_.radiationLength<<std::endl;
    std::cout<<"Atomic number Z : "<<steps.at(i).material_.Z<<std::endl;
    std::cout<<"Mass number A : "<<steps.at(i).material_.A<<std::endl;
    std::cout<<"mean excitation energy [eV] : "<<steps.at(i).material_.mEE<<std::endl;
  }
#endif

  pos = 10.*fitState.getPos(); //cm -> mm
  mom = fitState.getMom();
  tof = fitState.getTime() - time0;
  tracklen = len + tracklength; //cm -> mm
  return true;
}

bool HypTPCTask::ExtrapolateToTarget(int trackid, TVector3 &pos, TVector3 &mom, double &tracklen, double &tof, int repid) const{

  TVector3 target(0.,0.,ztgt);
  return ExtrapolateToPoint(trackid, target, pos, mom, tracklen, tof, repid);
}

bool HypTPCTask::ExtrapolateToTargetCenter(int trackid, TVector3 &pos, TVector3 &mom, double &tracklen, double &tof, int repid) const{

  return ExtrapolateToPlane(trackid, TgtPlane, pos, mom, tracklen, tof, repid);
}

bool HypTPCTask::IsInsideTarget(int trackid, int repid) const{

  TVector3 pos; TVector3 mom; double tracklen; double tof;
  if(ExtrapolateToTarget(trackid, pos, mom, tracklen, tof, repid) &&
     TMath::Abs(pos.x()) < 50. &&
     TMath::Abs(pos.y()) < 50. &&
     TMath::Abs(pos.z()-ztgt) < 50. &&
     TMath::Sqrt(pos.x()*pos.x()+pos.y()*pos.y()+(pos.z()-ztgt)*(pos.z()-ztgt)) < 50. &&
     -100. < tracklen && tracklen < 10.) return true;
  else return false;
}

bool HypTPCTask::ExtrapolateToHTOF(int trackid, int &candidates, int *ID, TVector3 *pos, TVector3 *mom, double *tracklen, double *tof, int repid) const{

  for(int i=0;i<8;i++){
    ID[i] = -1;
    pos[i] = TVector3(qnan, qnan, qnan);
    mom[i] = TVector3(qnan, qnan, qnan);
    tracklen[i] = qnan;
    tof[i] = qnan;
  }

  bool flag = false;
  candidates = 0;
  for(int i=0;i<8;i++){
    const TVector3 center = 10.*HTOFPlane[i] -> getO(); //cm -> mm
    TVector3 pos0; TVector3 mom0;
    double tracklen0; double tof0;
    if(ExtrapolateToPlane(trackid, HTOFPlane[i], pos0, mom0, tracklen0, tof0, repid)){
      TVector3 diff = pos0 - center;
      double xzdist = TMath::Sqrt(diff.x()*diff.x() + diff.z()*diff.z());
      TVector3 cross = center.Cross(diff);
      if(tof0<0) continue; //Opposite direction
      if(TMath::Abs(pos0.y()-12)>400. || xzdist>142.) continue;
      if(i==0){
	if(TMath::Abs(pos0.x())<71. && pos0.y()<62. && pos0.y()>-50.) continue; //window
	else if(cross.y()<0){
	  if(TMath::Abs(xzdist) >= 71.) ID[candidates] = 1;
	  else if(pos0.y() < 0.) ID[candidates] = 3;
	  else ID[candidates] = 2;
	}
	else{
	  if(TMath::Abs(xzdist) >= 71.) ID[candidates] = 6;
	  else if(pos0.y() < 0.) ID[candidates] = 5;
	  else ID[candidates] = 4;
	}
      }
      else{
	if(cross.y()<0){
	  if(TMath::Abs(xzdist) >= 71.) ID[candidates] = 3 + 4*i;
	  else ID[candidates] = 4 + 4*i;
	}
	else{
	  if(TMath::Abs(xzdist) < 71.) ID[candidates] = 5 + 4*i;
	  else ID[candidates] = 6 + 4*i;
	}
      }
      pos[candidates] = pos0;
      mom[candidates] = mom0;
      tracklen[candidates] = tracklen0;
      tof[candidates] = tof0;

      candidates++;
    }
  }

  if(candidates > 0) flag = true;
  return flag;
}

bool HypTPCTask::FindVertex(int trackid1, int trackid2, int repid1, int repid2, double &extrap_dist1, double &extrap_dist2, TVector3 &mom_vertex1, TVector3 &mom_vertex2, double &distance, TVector3 &vertex, double scan_range) const{

  distance = 10000.;
  int MaxStep  = 80000;
  double StepSize = 1.; // mm (1st naive scanning)
  double fineStepSize = StepSize*0.01; // (2nd fine scanning)
  int iStep1 = 0; int iStep2 = 0;
  double dist1 = 0.; double dist2 = 0.;
  if(!TrackCheck(trackid1) || !TrackCheck(trackid2)) return false;
  while(iStep2 < MaxStep){
    iStep1 = 0;
    while(iStep1 < MaxStep){
      //std::cout<<"istep "<<iStep1<<" "<<iStep2<<std::endl;
      TVector3 pos1; TVector3 mom1; TVector3 pos2; TVector3 mom2;
      bool extrapol1 = ExtrapolateTrack(trackid1, -iStep1*StepSize, pos1, mom1, repid1);
      //if(extrapol1) std::cout<<"extrapolate id "<<trackid1<<", step size "<<-iStep1*StepSize<<" mm, mom "<<mom1.Mag()<<" GeV/c, pos "<<pos1<<std::endl;
      bool extrapol2 = ExtrapolateTrack(trackid2, -iStep2*StepSize, pos2, mom2, repid2);
      //if(extrapol2) std::cout<<"extrapolate id "<<trackid2<<", step size "<<-iStep2*StepSize<<" mm, mom "<<mom2.Mag()<<" GeV/c, pos "<<pos2<<std::endl;
      if(!extrapol1 || !extrapol2){
	std::cout<<"extrapolation failed"<<std::endl;
	return false;
      }

      //Search the closest point(Vertex)
      TVector3 diff = pos1 - pos2;
      if(distance > diff.Mag()){
	distance = diff.Mag();
	dist1 = -iStep1*StepSize;
	dist2 = -iStep2*StepSize;
	vertex = pos1 + pos2;
	vertex *= 0.5;
	mom_vertex1 = mom1;
	mom_vertex2 = mom2;
      }

      //std::cout<<" dist "<<iStep1*StepSize<<" < scan range "<<scan_range<<std::endl;
      if(iStep1*StepSize >= scan_range) break;
      iStep1++;
    } //while(++iStep1 < MaxStep)
    if(iStep2*StepSize >= scan_range) break;
    iStep2++;
  } //while(++iStep2 < MaxStep){

  iStep1 = 0; iStep2 = 0;
  while(iStep2 < MaxStep){
    iStep1 = 0;
    while(iStep1 < MaxStep){
      //std::cout<<"istep "<<iStep1<<" "<<iStep2<<std::endl;
      double finedist1 = dist1 - StepSize + iStep1*fineStepSize;
      double finedist2 = dist2 - StepSize + iStep2*fineStepSize;

      TVector3 pos1; TVector3 mom1; TVector3 pos2; TVector3 mom2;
      bool extrapol1 = ExtrapolateTrack(trackid1, finedist1, pos1, mom1, repid1);
      //if(extrapol1) std::cout<<"extrapolate id "<<trackid1<<", step size "<<fineStepSize<<" mm, mom "<<mom1.Mag()<<" GeV/c, pos "<<pos1<<std::endl;
      bool extrapol2 = ExtrapolateTrack(trackid2, finedist2, pos2, mom2, repid2);
      //if(extrapol2) std::cout<<"extrapolate id "<<trackid2<<", step size "<<fineStepSize<<" mm, mom "<<mom2.Mag()<<" GeV/c, pos "<<pos2<<std::endl;

      if(!extrapol1 || !extrapol2){
	std::cout<<"extrapolation failed"<<std::endl;
	return false;
      }

      //Search the closest point(Vertex)
      TVector3 diff = pos1 - pos2;
      if(distance > diff.Mag()){
	distance = diff.Mag();
	vertex = pos1 + pos2;
	vertex *= 0.5;
	mom_vertex1 = mom1;
	mom_vertex2 = mom2;
	extrap_dist1 = finedist1;
	extrap_dist2 = finedist2;
      }
      //std::cout<<" dist "<<iStep1*fineStepSize<<" <= scan range "<<2.*StepSize<<std::endl;
      if(iStep1*fineStepSize >= 2.*StepSize) break;
      iStep1++;
    }
    if(iStep2*fineStepSize >= 2.*StepSize) break;
    iStep2++;
  }
  //std::cout<<"find vertex end "<<std::endl;
  //std::cout<<std::endl;

  return true;
}

bool HypTPCTask::FindVertexXi(int trackid, int repid, TVector3 decayvtx_lambda, TVector3 mom_lambda, double &tracklen_lambda, double &extrap_dist_pi, TVector3 &mom_pi_vertex, double &distance, TVector3 &vertex, double scan_range, double res1, double res2, double phi) const{ //trackid & repid : pi-

  bool status = false;

  distance = 10000.;
  int MaxStep  = 80000;
  double StepSize = 1.; // mm (1st naive scanning)
  double fineStepSize = StepSize*0.01; // (2nd fine scanning)
  int iStep = 0; double dist = 0.;
  if(!TrackCheck(trackid)) return false;
  iStep = 0;
  while(iStep < MaxStep){
    TVector3 pos; TVector3 mom;
    bool extrapol = ExtrapolateTrack(trackid, -iStep*StepSize, pos, mom, repid);
    if(!extrapol){
      std::cout<<"extrapolation failed"<<std::endl;
      return false;
    }

    //Search the closest point(Vertex)
    TVector3 AP = pos - decayvtx_lambda;
    TVector3 u = mom_lambda.Unit();
    Double_t dist_AX = u.Dot(AP);
    TVector3 AI = dist_AX*u;
    AI += decayvtx_lambda;
    Double_t lambdavtx_xivtx = (AI - decayvtx_lambda)*u;
    TVector3 diff = pos - AI;
    double r1 = sqrt(2)*res1 / hypot(res1,res2);
    double r2 = sqrt(2)*res2 / hypot(res1,res2);
    TVector3 e1(cos(phi), 0, sin(phi));
    TVector3 e2(-sin(phi), 0, cos(phi));

    double dy = diff.y();
    double dx = diff * e1 / r1;
    double dz = diff * e2 / r2;
    diff = TVector3(dx,dy,dz);
    if(distance > diff.Mag() && lambdavtx_xivtx < 0){
      distance = diff.Mag();
      dist = -iStep*StepSize;
      vertex = pos + AI;
      vertex *= 0.5;
      mom_pi_vertex = mom;
      status = true;
    }

    if(iStep*StepSize >= scan_range) break;
    iStep++;
  } //while(++iStep < MaxStep)
  if(!status) return false;

  iStep = 0;
  while(iStep < MaxStep){
    double finedist = dist - StepSize + iStep*fineStepSize;

    TVector3 pos; TVector3 mom;
    bool extrapol = ExtrapolateTrack(trackid, finedist, pos, mom, repid);
    if(!extrapol){
      std::cout<<"extrapolation failed"<<std::endl;
      return false;
    }

    //Search the closest point(Vertex)
    TVector3 AP = pos - decayvtx_lambda;
    TVector3 u = mom_lambda.Unit();
    Double_t dist_AX = u.Dot(AP);
    TVector3 AI = dist_AX*u;
    AI += decayvtx_lambda;
    Double_t lambdavtx_xivtx = (AI - decayvtx_lambda)*u;
    TVector3 diff = pos - AI;
    double r1 = sqrt(2)*res1 / hypot(res1,res2);
    double r2 = sqrt(2)*res2 / hypot(res1,res2);
    TVector3 e1(cos(phi), 0, sin(phi));
    TVector3 e2(-sin(phi), 0, cos(phi));

    double dy = diff.y();
    double dx = diff * e1 / r1;
    double dz = diff * e2 / r2;
    diff = TVector3(dx,dy,dz);
    if(distance > diff.Mag() && lambdavtx_xivtx < 0){
      distance = diff.Mag();
      vertex = pos + AI;
      vertex *= 0.5;
      mom_pi_vertex = mom;
      tracklen_lambda = TMath::Abs(lambdavtx_xivtx);
      extrap_dist_pi = finedist;
      status = true;
    }
    if(iStep*fineStepSize >= 2.*StepSize) break;
    iStep++;
  }

  return status;
}

bool HypTPCTask::TPCHTOFTrackMatching(int trackid, int repid, TVector3 vertex, std::vector<Double_t> HtofSeg, std::vector<Double_t> posHtof, int &htofhitid, double &tof, double &tracklen, TVector3 &pos, double &vertex_dist) const{

  bool status = false;
  double PosDiffCut = 60.; //mm

  TVector3 vtx_pos; TVector3 vtx_mom; double vtx_len; double vtx_tof;
  if(!ExtrapolateToPoint(trackid, vertex, vtx_pos, vtx_mom,
			 vtx_len, vtx_tof, repid)) return status;

  TVector3 diff = vertex - vtx_pos;
  vertex_dist = diff.Mag();

  int candidates; int extrap_id[8]; TVector3 extrap_pos[8];
  TVector3 extrap_mom[8]; double extrap_len[8]; double extrap_tof[8];
  if(!ExtrapolateToHTOF(trackid, candidates, extrap_id, extrap_pos,
			extrap_mom, extrap_len, extrap_tof, repid)) return status;

  double mintof = 10000;
  for(int i=0;i<candidates;i++){
    int nhHtof = HtofSeg.size();
    for(int j=0;j<nhHtof;j++){
      if(extrap_id[i] == (int) HtofSeg[j] &&
	 TMath::Abs(posHtof[j] - extrap_pos[i].y()) < PosDiffCut){

	double tof_htof2vtx = extrap_tof[i] - vtx_tof;
	double tracklen_htof2vtx = extrap_len[i] - vtx_len;
	if(tof_htof2vtx < 0 || tracklen_htof2vtx < 0) continue;

	if(mintof > tof_htof2vtx){
	  pos = extrap_pos[i];
	  tof = tof_htof2vtx;
	  tracklen = tracklen_htof2vtx;
	  htofhitid = j;

	  mintof = tof;
	  status = true;
	}
      }
    }
  }

  return status;
}

bool HypTPCTask::TPCHTOFTrackMatching(int trackid, int repid, std::vector<Double_t> HtofSeg, std::vector<Double_t> posHtof, int &htofhitid, double &tof, double &tracklen, TVector3 &pos) const{

  bool status = false;
  double PosDiffCut = 60.; //mm
  int candidates; int extrap_id[8]; TVector3 extrap_pos[8]; TVector3 extrap_mom[8]; double extrap_len[8]; double extrap_tof[8];
  if(!ExtrapolateToHTOF(trackid, candidates, extrap_id, extrap_pos, extrap_mom, extrap_len, extrap_tof, repid))
    return status;

  double mintof = 10000;
  for(int i=0;i<candidates;i++){
    int nhHtof = HtofSeg.size();
    for(int j=0;j<nhHtof;j++){
      if(extrap_id[i] == (int) HtofSeg[j] &&
	 TMath::Abs(posHtof[j] - extrap_pos[i].y()) < PosDiffCut){
	if(mintof > extrap_tof[i]){
	  pos = extrap_pos[i];
	  tof = extrap_tof[i];
	  tracklen = extrap_len[i];
	  mintof = tof;
	  htofhitid = j;
	  status = true;
	}
      }
    }
  }
  return status;
}

double HypTPCTask::DistLambdaTarget(TVector3 decayvtx_lambda, TVector3 mom_lambda) const{

  //Search the closest point
  TVector3 target(0.,0.,ztgt);
  TVector3 AP = target - decayvtx_lambda;
  TVector3 u = mom_lambda.Unit();
  Double_t dist_AX = u.Dot(AP);
  TVector3 AI = dist_AX*u;
  AI += decayvtx_lambda;
  TVector3 diff = target - AI;
  return diff.Mag();
}

bool HypTPCTask::XiDecayToProdVertex(int trackid, TVector3 kkvtx, TVector3 &xiprodvtx, TVector3 &mom, double &tracklen, double &tof) const{

  return ExtrapolateToPoint(trackid, kkvtx, xiprodvtx, mom, tracklen, tof, -1);
  //bool extrap = ExtrapolateToPoint(trackid, kkvtx, xiprodvtx, mom, tracklen, tof, -1);
  //bool signmomloss = (tracklen < 0); //Xi track should loss momentum in the target meterial(not gain).
  //if(extrap && signmomloss) return true;
  //else return false;
}
