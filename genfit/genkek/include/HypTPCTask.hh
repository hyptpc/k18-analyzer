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

  HypTPCTask();
  ~HypTPCTask();
  static HypTPCTask& GetInstance();

  genfit::Track* GetFittedTrack(int trackid) const;
  genfit::AbsTrackRep* GetTrackRep(int trackid) const;
  genfit::FitStatus* GetFitStatus(int trackid) const;
  genfit::MeasuredStateOnPlane GetFitState(int trackid) const;

  //Parameters
  bool GetTrackPull(int trackid, int pdg, TVector3 res_vect, double tracklen, TVector3 g4mom, TVector3 g4pos, double *residual, double *pull, double *residual6D, double *pull6D) const;

  double GetChi2(int trackid) const;
  double GetNDF(int trackid) const;
  double GetChi2NDF(int trackid) const;
  double GetPvalue(int trackid) const;
  int GetNumOfIterations(int trackid) const;
  double GetCharge(int trackid) const;
  int GetPDGcode(int trackid) const;
  TVector3 GetMom(int trackid) const;
  TVector3 GetPos0(int trackid) const; //Get Vertex position
  int GetNHits(int trackid) const;
  double GetTrackLength(int trackid, int start=0, int end=-1) const;
  double GetTrackTOF(int trackid, int start=0, int end=-1) const;

  //Extrapolation
  bool ExtrapolateTrack(int trackid, double distance, TVector3 &pos, TVector3 &mom) const;
  bool ExtrapolateToPoint(int trackid, TVector3 point, TVector3 &pos, TVector3 &mom, double &tracklen, double &tof) const;
  bool ExtrapolateToPlane(int trackid, genfit::SharedPlanePtr plane, TVector3 &pos, TVector3 &mom, double &tracklen, double &tof) const;
  bool IsInsideTarget(int trackid) const;
  bool ExtrapolateToHTOF(int trackid, TVector3 &pos, TVector3 &mom, double &tracklen, double &tof) const;

private:

  genfit::SharedPlanePtr HTOFPlane[8];

  ClassDef(HypTPCTask, 1)

};  //class HypTPCTask

#endif // genfit_HypTPCTask_hh
