//  Authors: Wooseung Jung

#ifndef HYPTPCFITPROCESS_HH
#define HYPTPCFITPROCESS_HH

//GenKEK
#include "HypTPCTrack.hh"

//GenFit
#include <Track.h>
#include <AbsTrackRep.h>
#include <Exception.h>

class HypTPCFitProcess : public HypTPCTrack{

public:

  HypTPCFitProcess(): HypTPCTrack(), verbosity(3){}
  virtual ~HypTPCFitProcess(){}
  /*!
   * Verbose control:
   * -1: Silent, 0: Minimum
   * 1: Errors only, 2: Errors and Warnings
   * 3: Verbose mode, long term debugging (default)
   */
  void DebugMode(); //for debugging mode

  int GetVerbosity() const { return verbosity; }
  void SetVerbosity(int v){
    this->verbosity = v;
    if(verbosity >= 1) genfit::Exception::quiet(false);
    else genfit::Exception::quiet(true);
  }


  bool FitCheck(int trackid) const;
  bool FitCheck(genfit::Track* fittedTrack, genfit::AbsTrackRep* rep=nullptr) const;

  //Process track with the all AbsTrackReps of the track.
  bool ProcessTrack(genfit::Track* Track);
  //Process track with one AbsTrackRep of the Track
  bool ProcessTrack(genfit::Track* Track, genfit::AbsTrackRep* rep);
  //Process all tracks with its all AbsTrackReps.
  void FitTracks();

protected:

  int verbosity;

private:

  ClassDef(HypTPCFitProcess, 1)

};  //class HypTPCFitProcess

#endif // genfit_HypTPCFitProcess_hh
