
//  Authors: Wooseung Jung

//GenKEK
#include "HypTPCTask.hh"

//GenFit
#include <Track.h>
#include <AbsTrackRep.h>
#include <Exception.h>
#include <TrackPoint.h>
#include <FitStatus.h>
#include <KalmanFitterInfo.h>
#include <MeasuredStateOnPlane.h>

//ROOT
#include <TMath.h>

//STL
#include <cstddef>
#include <iostream>
#include <string>

//ClassImp(HypTPCTask)

namespace{
  const double qnan = TMath::QuietNaN();
}

HypTPCTask* HypTPCTask::GetInstance(){
  return new HypTPCTask();
}

genfit::FitStatus* HypTPCTask::GetFitStatus(int trackid) const{
  //genfit::FitStatus* HypTPCTask::GetFitStatus(int trackid){

  genfit::Track* fittedTrack = GetTrack(trackid);
  genfit::FitStatus *fitStatus;
  try{ fitStatus = fittedTrack -> getFitStatus(); }
  catch (genfit::Exception &e){
    return nullptr;
  }
  return fitStatus;

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
  if(!TMath::IsNaN(chi2)
     && !TMath::IsNaN(ndf)) chi2ndf = chi2/ndf;
  return chi2ndf;
}

double HypTPCTask::GetCharge(int trackid) const
{
  double charge = qnan;
  genfit::Track* fittedTrack = GetTrack(trackid);
  if(!fittedTrack) return charge;
  try{ genfit::TrackPoint* tp_base = nullptr;
    if(fittedTrack->getNumPointsWithMeasurement() > 0) tp_base = fittedTrack->getPointWithMeasurement(0);
    if(!tp_base) return charge;

    genfit::AbsTrackRep* rep = fittedTrack->getCardinalRep();
    if(rep){
      genfit::KalmanFitterInfo* kfi = static_cast<genfit::KalmanFitterInfo*>(tp_base->getFitterInfo(rep));
      if(!kfi) return charge;
      const genfit::MeasuredStateOnPlane* state = &(kfi->getFittedState(true));
      if(state) charge = rep->getCharge(*state);
    }
  }
  catch(...){
    std::cerr << "HypTPCTask::GetCharge - Error - obtaining charge failed. Returning NAN as charge." << std::endl;
  }
  return charge;

}

TVector3 HypTPCTask::GetMom(int trackid) const
{
  TVector3 mom(0, 0, 0);

  genfit::Track* fittedTrack = GetTrack(trackid);
  if(!fittedTrack) return mom;
  genfit::TrackPoint* tp_base = nullptr;
  if(fittedTrack->getNumPointsWithMeasurement() > 0) tp_base = fittedTrack->getPointWithMeasurement(0);
  if(!tp_base) return mom;
  genfit::AbsTrackRep* rep = fittedTrack->getCardinalRep();
  if(rep){
    genfit::KalmanFitterInfo* kfi = static_cast<genfit::KalmanFitterInfo*>(tp_base->getFitterInfo(rep));
    if(!kfi) return mom;
    const genfit::MeasuredStateOnPlane* state = &(kfi->getFittedState(true));
    if(state) return state->getMom();
  }
  return mom;

}

double HypTPCTask::GetTrackLength(int trackid, int start, int end) const{

  double length = qnan;
  genfit::Track* fittedTrack = GetTrack(trackid);
  if (!fittedTrack) return length;
  length = fittedTrack -> getTrackLen(nullptr,start,end);
  return length;

}

double HypTPCTask::GetTrackTOF(int trackid, int start, int end) const{

  double TOF = qnan;
  genfit::Track* fittedTrack = GetTrack(trackid);
  if (!fittedTrack) return TOF;
  TOF = fittedTrack -> getTOF(nullptr,start,end);
  return TOF;

}
