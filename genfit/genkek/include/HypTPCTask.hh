//  Authors: Wooseung Jung

#ifndef HYPTPCTASK_HH
#define HYPTPCTASK_HH

//GenKEK
#include "HypTPCFitProcess.hh"

//GenFit
#include <FitStatus.h>

//Root
#include <TVector3.h>                         // for TVector3

class HypTPCTask : public HypTPCFitProcess{

public:

  HypTPCTask() : HypTPCFitProcess() {}
  //HypTPCTask(){}
  ~HypTPCTask(){}
  static HypTPCTask* GetInstance();

  genfit::FitStatus* GetFitStatus(int trackid) const;
  double GetChi2(int trackid) const;
  double GetNDF(int trackid) const;
  double GetChi2NDF(int trackid) const;
  double GetCharge(int trackid) const;
  TVector3 GetMom(int trackid) const;
  //Track length[cm] & TOF[ns] between start point to end point
  double GetTrackLength(int trackid, int start=0, int end=-1) const;
  double GetTrackTOF(int trackid, int start=0, int end=-1) const;

  /*
  //HypTPCFitprocess
  int GetVerbosity() const;
  void SetVerbosity(int v);
  bool FitCheck(genfit::Track* fittedTrack, genfit::AbsTrackRep* rep=nullptr);
  //Process track with the all AbsTrackReps of the track.
  bool ProcessTrack(genfit::Track* Track);
  //Process track with one AbsTrackRep of the Track
  bool ProcessTrack(genfit::Track* Track, genfit::AbsTrackRep* rep);
  //Process all tracks with its all AbsTrackReps.
  void FitTracks();

  //HypTPCTrack
  void Init();
  void AddHelixTrack(int pdg, TPCLocalTrackHelix *tp);
  genfit::Track* GetTrack(int ith) const;
  int GetNTrack() const;
  */
private:

  //ClassDef(HypTPCTask, 1)

};  //class HypTPCTask

#endif // genfit_HypTPCTask_hh
