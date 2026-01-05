// -*- C++ -*-

#ifndef TPC_ANALYZER_HH
#define TPC_ANALYZER_HH

#include <vector>
#include <TString.h>
#include <TVector3.h>

#include "DetectorID.hh"

class TPCRawData;
class TPCHit;

typedef std::vector<TPCHit*>        TPCHitContainer;

//_____________________________________________________________________________
class TPCAnalyzer
{
public:
  static const TString& ClassName();
  TPCAnalyzer();
  ~TPCAnalyzer();

private:
  TPCAnalyzer(const TPCAnalyzer&);
  TPCAnalyzer& operator =(const TPCAnalyzer&);

private:
  enum e_type
  { kTPC, kTPCTracking, kTPCK18, n_type };
  std::vector<Bool_t>                m_is_decoded;
  std::vector<TPCHitContainer>       m_TPCHitCont;

public:

  //TPC Hit&Cluster
  Bool_t DecodeTPCHits(TPCRawData &TPCrawData, const std::vector<Double_t> clock);

  const TPCHitContainer& GetTPCHC(Int_t l) const { return m_TPCHitCont.at(l); }

private:

protected:

  void ClearTPCHits();

};

//_____________________________________________________________________________
inline const TString&
TPCAnalyzer::ClassName()
{
  static TString s_name("TPCAnalyzer");
  return s_name;
}

#endif
