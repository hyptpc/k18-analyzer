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

bool HypTPCFitProcess::FitCheck(genfit::Track* fittedTrack, genfit::AbsTrackRep* rep) const {

  if(rep==nullptr) return false;
  if(!fittedTrack->getFitStatus(rep)->isFitted()){
    if(verbosity>=2) LogWARNING("Fitting is failed");
    return false;
  }
  if(!fittedTrack->getFitStatus(rep)->isFitConverged()){
    if(verbosity>=2) LogWARNING("Fit is not converged");
    return false;
  }
  try{fittedTrack->getFittedState();}
  catch(genfit::Exception& e){
    if(verbosity>=2) LogWARNING("Track has no fitted state, failed!");
    if(verbosity>=1) std::cerr << e.what();
    return false;
  }

  try{fittedTrack->checkConsistency();}
  catch(genfit::Exception& e){
    if(verbosity>=2) LogWARNING("genfit::Track::checkConsistency() failed!");
    if(verbosity>=1) std::cerr << e.what();
    return false;
  }
  return true;
}
bool HypTPCFitProcess::FitCheck(int trackid) const {

  genfit::Track* fittedTrack = GetTrack(trackid);
  genfit::AbsTrackRep* rep = fittedTrack->getCardinalRep();
  return FitCheck(fittedTrack, rep);

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

  if(!FitCheck(Track)) return false;
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

  if(!FitCheck(Track, rep)) return false;
  return true;

}
