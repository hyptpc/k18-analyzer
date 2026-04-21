// -*- C++ -*-

#ifndef TPC_EVENT_ANALYZER_HH
#define TPC_EVENT_ANALYZER_HH

#include <vector>

#include <TVector3.h>

#include "DetectorID.hh"
#include "Event.hh"
#include "HistTools.hh"

class TPCRawData;
class TPCAnalyzer;
class TPCLTrackHit;
class TPCLocalTrack;

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

  static Double_t CalcTruncatedMean(const std::vector<Double_t>& cumulative_vec, Double_t fraction);

  void FillTrkHist(const TPCLocalTrack* track);
  void FillTrkHitHist(TPCLTrackHit* hit, const TPCLocalTrack* track);
  void FillResidualHist(
    const TString& prefix,
    Int_t layer, Int_t row, Int_t pad,
    Double_t res_x, Double_t res_y,
    Double_t ref_x, Double_t ref_y,
    Bool_t in_window,
    Double_t ctime, const TVector3& local_pos
  );

private:
  inline static Bool_t m_dst_calib_flag = false;

  std::vector<Double_t> m_clkTpc;
};

#endif
