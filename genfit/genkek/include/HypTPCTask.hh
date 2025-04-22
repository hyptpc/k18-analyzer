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
  genfit::AbsTrackRep* GetTrackRep(int trackid, int repid=-1) const;
  genfit::FitStatus* GetFitStatus(int trackid, int repid=-1) const;
  genfit::MeasuredStateOnPlane GetFitState(int trackid, int pointid, int repid=-1) const;

  //Parameters
  bool GetTrackPull(int trackid, int pdg, TVector3 res_vect, double tracklen, TVector3 g4mom, TVector3 g4pos, double *residual, double *pull, double *residual6D, double *pull6D) const;

  double GetChi2(int trackid, int repid=-1) const;
  double GetNDF(int trackid, int repid=-1) const;
  double GetChi2NDF(int trackid, int repid=-1) const;
  double GetPvalue(int trackid, int repid=-1) const;
  int GetNumOfIterations(int trackid, int repid=-1) const;
  double GetCharge(int trackid, int repid=-1) const;
  int GetPDGcode(int trackid, int repid=-1) const;
  TVector3 GetMom(int trackid, int pointid, int repid=-1) const;
  TVector3 GetPos(int trackid, int pointid, int repid=-1) const;
  int GetNHits(int trackid) const;
  double GetTrackLength(int trackid, int start, int end, int repid=-1) const;
  double GetTrackTOF(int trackid, int start, int end, int repid=-1) const;

  //Extrapolation
  bool ExtrapolateTrack(int trackid, double distance, TVector3 &pos, TVector3 &mom, int repid=-1) const;
  bool ExtrapolateToPoint(int trackid, TVector3 point, TVector3 &pos, TVector3 &mom, double &tracklen, double &tof, int repid=-1) const;
  bool ExtrapolateToPlane(int trackid, genfit::SharedPlanePtr plane, TVector3 &pos, TVector3 &mom, double &tracklen, double &tof, int repid=-1) const;
  bool IsInsideTarget(int trackid, int repid=-1, bool Beamthrough = 0) const;
  bool ExtrapolateToHTOF(int trackid, int &candidates, int *ID, TVector3 *pos, TVector3 *mom, double *tracklen, double *tof, int repid=-1) const;
  bool ExtrapolateToTarget(int trackid, TVector3 &pos, TVector3 &mom, double &tracklen, double &tof, int repid=-1) const;
  bool ExtrapolateToTargetCenter(int trackid, TVector3 &pos, TVector3 &mom, double &tracklen, double &tof, int repid=-1) const;
  bool FindVertex(int trackid1, int trackid2, int repid1, int repid2, double &extrap_dist1, double &extrap_dist2, TVector3 &mom_vertex1, TVector3 &mom_vertex2, double &distance, TVector3 &vertex, double scan_range=100.) const;
  bool FindVertex(int trackid1, int trackid2, int repid1, int repid2, double &extrap_dist1, double &extrap_dist2, TVector3 &mom_vertex1, TVector3 &mom_vertex2, double &distance, TVector3 &vertex, double scan_range, TVector3 &track1_vertex, TVector3& track2_vertex) const;
  bool FindVertexXi(int trackid, int repid, TVector3 decayvtx_lambda, TVector3 mom_lambda, double &tracklen_lambda, double &extrap_dist_pi, TVector3 &mom_pi_vertex, double &distance, TVector3 &vertex, double scan_range=100., double res1 = 1,double res2 = 1, double phi = 0) const;
  bool FindVertexXi(int trackid, int repid, TVector3 decayvtx_lambda, TVector3 mom_lambda, double &tracklen_lambda, double &extrap_dist_pi, TVector3 &mom_pi_vertex, double &distance, TVector3 &vertex, double scan_range, double res1 ,double res2 , double phi , TVector3& track1_vertex,TVector3& track2_vertex) const;
  bool TPCHTOFTrackMatching(int trackid, int repid, TVector3 vertex, std::vector<Double_t> HtofSeg, std::vector<Double_t> posHtof, int &htofhitid, double &tof, double &tracklen, TVector3 &pos, double &vertex_dist) const;
  bool TPCHTOFTrackMatching(int trackid, int repid, std::vector<Double_t> HtofSeg, std::vector<Double_t> posHtof, int &htofhitid, double &tof, double &tracklen, TVector3 &pos) const;
  double DistLambdaTarget(TVector3 decayvtx_lambda, TVector3 mom_lambda) const;
  bool XiDecayToProdVertex(int trackid, TVector3 kkvtx, TVector3 &xiprodvtx, TVector3 &mom, double &tracklen, double &tof) const;
private:

  genfit::SharedPlanePtr HTOFPlane[8];
  genfit::SharedPlanePtr TgtPlane;

  ClassDef(HypTPCTask, 1)

};  //class HypTPCTask

#endif // genfit_HypTPCTask_hh
