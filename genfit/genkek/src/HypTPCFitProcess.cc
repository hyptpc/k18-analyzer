//  Authors: Wooseung Jung

//GenKEK
#include "HypTPCFitter.hh"
#include "HypTPCFitProcess.hh"

//STL
#include <cstddef>
#include <iostream>
#include <string>

#define LogDEBUG(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogERROR(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogWARNING(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl

ClassImp(HypTPCFitProcess)

namespace{
  genfit::AbsKalmanFitter* _fitter = HypTPCFitter::GetFitter();
}

bool HypTPCFitProcess::FitCheck(genfit::Track* fittedTrack, genfit::AbsTrackRep* rep){

  if(!fittedTrack->getFitStatus(rep)->isFitted()){
    if (verbosity >= 2) LogWARNING("Fitting is failed");
    return false;
  }
  if(!fittedTrack->getFitStatus(rep)->isFitConverged()){
    if (verbosity >= 2) LogWARNING("Fit is not converged");
    return false;
  }
  try{fittedTrack->checkConsistency();}
  catch(genfit::Exception& e){
    if(verbosity >= 2){
      std::cerr << "genfit::Track::checkConsistency() failed!" << std::endl;
      std::cerr << e.what();
    }
    return false;
  }
  return true;

}

void HypTPCFitProcess::FitTracks(){

  int nTracks = GetNTrack();
  for(int i=0; i<nTracks; i++){
    //ProcessTrack((genfit::Track*) _genfitTrackArray -> ConstructedAt(i));
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
