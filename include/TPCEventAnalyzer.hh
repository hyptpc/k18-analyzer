// -*- C++ -*-

#ifndef TPC_EVENT_ANALYZER_HH
#define TPC_EVENT_ANALYZER_HH

#include <TVector3.h>

#include "DetectorID.hh"
#include "Event.hh"
#include "HistTools.hh"

class TPCRawData;
class TPCAnalyzer;
class TPCLTrackHit;

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
  TPCEventAnalyzer();
  ~TPCEventAnalyzer();

private:

public:
  void TPCRawHit(const TPCRawData& TPCrawData);
  TPCBaselineInfo TPCBaselineHit(const TPCRawData& TPCrawData);
  Int_t TPCCorHit(const TPCRawData& TPCrawData);
  void TPCHit(const TPCAnalyzer& TPCAna);
  void SetClock(const std::vector<Double_t>& clkTpc) { m_clkTpc = clkTpc; }
  void FillCoBoClockTime(const TString& prefix, Int_t layer, Int_t row,
                         Double_t ctime, const TVector3& localPos, Double_t referenceY);

private:
  std::vector<Double_t> m_clkTpc;
};

#endif
