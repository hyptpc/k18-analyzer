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
  static HypTPCTask& GetInstance();

  genfit::FitStatus* GetFitStatus(int trackid) const;
  double GetChi2(int trackid) const;
  double GetNDF(int trackid) const;
  double GetChi2NDF(int trackid) const;
  double GetCharge(int trackid) const;
  TVector3 GetMom(int trackid) const;
  double GetTrackLength(int trackid, int start=0, int end=-1) const;
  double GetTrackTOF(int trackid, int start=0, int end=-1) const;

private:

  ClassDef(HypTPCTask, 1)

};  //class HypTPCTask

#endif // genfit_HypTPCTask_hh
