// -*- C++ -*-

#ifndef TPC_EVENT_ANALYZER_HH
#define TPC_EVENT_ANALYZER_HH

#include <string>
#include <vector>

#include <TVector3.h>

#include "DetectorID.hh"
#include "Event.hh"
#include "HistTools.hh"

class TPCRawData;
class TPCAnalyzer;
class TPCReconstructor;
class TPCVertex;
class TPCLTrackHit;
class TPCLocalTrack;
class TPCLocalTrackHelix;

struct TPCBaselineInfo {
  Int_t row = -1;
  Int_t layer = -1;
  Double_t rms = 0.0;
  Bool_t valid = kFALSE;
};

//_____________________________________________________________________________
class TPCEventAnalyzer
{
public:
  static void SetDstCalibFlag(Bool_t dst_calib_flag) { m_dst_calib_flag = dst_calib_flag; }
  static Bool_t GetDstCalibFlag() { return m_dst_calib_flag; }

  TPCEventAnalyzer();
  ~TPCEventAnalyzer();

private:

public:
  void TPCRawHit(const TPCRawData& TPCrawData);
  TPCBaselineInfo TPCBaselineHit(const TPCRawData& TPCrawData);
  Int_t TPCCorHit(const TPCRawData& TPCrawData);

  // NOTE: Disabled for now; this name collides with the TPCHit type name.
  // void TPCHit(const TPCAnalyzer& TPCAna);

  void SetClock(const std::vector<Double_t>& clkTpc) { m_clkTpc = clkTpc; }
  void FillCoBoClockTime(const TString& prefix, Int_t layer, Int_t row,
                         Double_t ctime, const TVector3& local_pos, Double_t ref_y);

  static Bool_t ValidateCoboClocks(const std::vector<Double_t>& clk_tpc);
  static std::string DecodePidCandidates(Int_t pid_code);

  static Double_t CalcTruncatedMean(const std::vector<Double_t>& cumulative_vec, Double_t fraction);

  void FillTrkHist(const TPCLocalTrack* track);
  void FillTrkHitHist(TPCLTrackHit* hit, const TPCLocalTrack* track);
  void FillTrkParamHist(
    Int_t nhits, Double_t chisqr,
    Double_t x0, Double_t y0, Double_t u0, Double_t v0
  );
  void FillBcOutTrackHist(
    Double_t chisqr_bcout,
    Double_t x0_bcout, Double_t y0_bcout,
    Double_t u0_bcout, Double_t v0_bcout,
    Double_t xtgt_bcout, Double_t ytgt_bcout,
    Double_t utgt_bcout, Double_t vtgt_bcout
  );
  void FillTPCBcOutTgtResidualHist(
    Double_t tpc_xtgt, Double_t tpc_ytgt, Double_t tpc_utgt, Double_t tpc_vtgt,
    Double_t bcout_xtgt, Double_t bcout_ytgt, Double_t bcout_utgt, Double_t bcout_vtgt
  );
  void FillResidualHist(
    const TString& prefix,
    Int_t layer, Int_t row, Int_t pad,
    Double_t res_x, Double_t res_y,
    Double_t ref_x, Double_t ref_y,
    Bool_t in_window,
    Double_t ctime, const TVector3& local_pos
  );
  void FillTPCBcOutTrackingResidualPullHist(
    Int_t layer, Int_t center_row, Bool_t valid_tpc_resolution,
    const TVector3& trk_res_global,
    const TVector3& trk_res_local,
    const TVector3& trk_pull_global,
    const TVector3& trk_pull_local,
    Double_t cl_ref_x, Double_t cl_ref_y_tpc, Double_t cl_ref_y_bcout,
    const TVector3& cl_res
  );
  void FillHelixHitHist(TPCLTrackHit* hit, Bool_t fill_cluster_detail, Int_t track_pid);
  void FillHelixLambdaMassHist(const TPCVertex* vertex);

private:
  inline static Bool_t m_dst_calib_flag = false;

  std::vector<Double_t> m_clkTpc;
};

#endif
