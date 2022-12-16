//  Authors: Wooseung Jung

//GenKEK
#include "HypTPCFitter.hh"
#include "HypTPCFitProcess.hh"

//STL
#include <cstddef>
#include <iostream>
#include <string>

#define LogWARNING(exp) std::cout<<"WARNING: "<< __FILE__<<"::"<<__func__<<" line "<<__LINE__<<" "<<exp<<std::endl

ClassImp(HypTPCFitProcess)

void HypTPCFitProcess::DebugMode(){ HypTPCFitter::_fitter -> setDebugLvl(); }

bool HypTPCFitProcess::TrackCheck(genfit::Track* fittedTrack, genfit::AbsTrackRep* rep) const {

  if(rep==nullptr) return false;
  if(!fittedTrack->getFitStatus(rep)->isFitted()){
    if(verbosity>=1) LogWARNING("Fitting is failed");
    return false;
  }
  if(!fittedTrack->getFitStatus(rep)->isFitConverged()){
    if(verbosity>=1) LogWARNING("Fit is not converged");
    return false;
  }
  try{fittedTrack->checkConsistency();}
  catch(genfit::Exception& e){
    if(verbosity>=2) std::cerr << e.what();
    if(verbosity>=1) LogWARNING("genfit::Track::checkConsistency() failed!");
    return false;
  }
  int nhits = 0;
  try{nhits = fittedTrack->getNumPointsWithMeasurement();}
  catch(genfit::Exception &e){
    if(verbosity>=2) LogWARNING("genfit::Track::getNumPointsWithMeasurement() failed!");
    if(verbosity>=1) std::cerr << e.what();
    return false;
  }
  for(int pointid=0;pointid<nhits;pointid++){
    try{fittedTrack->getFittedState(pointid);}
    catch(genfit::Exception& e){
      if(verbosity>=2) std::cerr << e.what();
      if(verbosity>=1) LogWARNING("Track has no fitted state, failed!");
      return false;
    }
  }
  return true;
}

bool HypTPCFitProcess::TrackCheck(int trackid) const {

  genfit::Track* fittedTrack = GetTrack(trackid);
  genfit::AbsTrackRep* rep = fittedTrack->getCardinalRep();
  return TrackCheck(fittedTrack, rep);
}

void HypTPCFitProcess::FitTracks(){

  int nTracks = GetNTrack();
  for(int i=0; i<nTracks; i++){
    ProcessTrack(GetTrack(i));
  }
}

bool HypTPCFitProcess::ProcessTrack(genfit::Track* Track){

  try{ HypTPCFitter::_fitter->processTrack(Track); }
  catch(genfit::Exception& e){
    if(verbosity >= 1){
      std::cerr << "HypTPCFitter::processTrack::Exception: \n";
      std::cerr << e.what();
      std::cerr << "Exception, next track" << std::endl;
    }
    return false;
  }

  if(!TrackCheck(Track)) return false;
  return true;
}

bool HypTPCFitProcess::ProcessTrack(genfit::Track* Track, genfit::AbsTrackRep* rep){

  try{ HypTPCFitter::_fitter->processTrackWithRep(Track, rep, false); }
  catch(genfit::Exception& e){
    if(verbosity >= 1){
      std::cerr << "HypTPCFitter::processTrack::Exception: \n";
      std::cerr << e.what();
      std::cerr << "Exception, next track" << std::endl;
    }
    return false;
  }

  if(!TrackCheck(Track, rep)) return false;
  return true;
}
