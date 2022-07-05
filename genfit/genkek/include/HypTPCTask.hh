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
  ~HypTPCTask(){}
  HypTPCTask* getInstance(){ return new HypTPCTask(); }

  genfit::FitStatus* GetFitStatus(int trackid) const;
  double GetChi2(int trackid) const;
  double GetNDF(int trackid) const;
  double GetChi2NDF(int trackid) const;
  double GetCharge(int trackid) const;
  TVector3 GetMom(int trackid) const;
  //Track length[cm] & TOF[ns] between start point to end point
  double GetTrackLength(int trackid, int start=0, int end=-1) const;
  double GetTrackTOF(int trackid, int start=0, int end=-1) const;

  //HypTPCFitprocess.hh
  int Get_verbosity() const override;
  void Set_verbosity(int v) override;
  bool FitCheck(genfit::Track* fittedTrack, genfit::AbsTrackRep* rep=nullptr) override;
  //Process track with the all AbsTrackReps of the track.
  bool ProcessTrack(genfit::Track* Track) override;
  //Process track with one AbsTrackRep of the Track
  bool ProcessTrack(genfit::Track* Track, genfit::AbsTrackRep* rep) override;
  //Process all tracks with its all AbsTrackReps.
  void FitTracks() override;

  //HypTPCTrack.hh
  void Init() override;
  void AddHelixTrack(int pdg, TPCLocalTrackHelix *tp) override;
  genfit::Track* GetTrack(int ith) const override;
  int GetNTrack() const override;

private:

  ClassDef(HypTPCTask, 1)

};  //class HypTPCTask

#endif // genfit_HypTPCTask_hh
