//  Authors: Wooseung Jung

#ifndef HYPTPCTASK_HH
#define HYPTPCTASK_HH

//GenKEK
#include "HypTPCFitProcess.hh"

//GenFit
#include <FitStatus.h>
#include <AbsTrackRep.h>
#include <MeasuredStateOnPlane.h>
#include <SharedPlanePtr.h>

//Root
#include <TVector3.h>                         // for TVector3

class HypTPCTask : public HypTPCFitProcess{

public:

  HypTPCTask() : HypTPCFitProcess() {}
  ~HypTPCTask(){}
  static HypTPCTask& GetInstance();

  genfit::Track* GetFittedTrack(int trackid) const;
  genfit::AbsTrackRep* GetTrackRep(int trackid) const;
  genfit::FitStatus* GetFitStatus(int trackid) const;
  genfit::MeasuredStateOnPlane GetFitState(int trackid) const;

  //Parameters
  double GetChi2(int trackid) const;
  double GetNDF(int trackid) const;
  double GetChi2NDF(int trackid) const;
  double GetCharge(int trackid) const;
  TVector3 GetMom(int trackid) const;
  TVector3 GetPos0(int trackid) const; //Get Vertex position
  double GetTrackLength(int trackid, int start=0, int end=-1) const;
  double GetTrackTOF(int trackid, int start=0, int end=-1) const;

  //Extrapolation
  bool ExtrapolateTrack(int trackid, double distance, TVector3 &pos) const;
  bool ExtrapolateToPoint(int trackid, TVector3 point, TVector3 &pos) const;
  bool GetPosOnPlane(int trackid, genfit::SharedPlanePtr plane, TVector3 &pos) const;
  bool IsInsideTarget(int trackid) const;

private:

  ClassDef(HypTPCTask, 1)

};  //class HypTPCTask

#endif // genfit_HypTPCTask_hh
