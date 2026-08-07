// -*- C++ -*-

#ifndef EVENT_ANALYZER_HH
#define EVENT_ANALYZER_HH

#include <vector>

#include "DetectorID.hh"
#include "Event.hh"
#include "HistTools.hh"

class DCAnalyzer;
class HodoAnalyzer;
class RawData;
class D5Track;

//_____________________________________________________________________________
class EventAnalyzer
{
public:
  EventAnalyzer();
  ~EventAnalyzer();

private:

public:
  beam::EBeamFlag BeamFlag(const RawData& rawData);

  void HodoRawHit(const RawData& rawData, beam::EBeamFlag beam_flag=beam::kAll);
  void HodoHit(const HodoAnalyzer& hodoAna, beam::EBeamFlag beam_flag=beam::kAll);
  void HodoCluster(const HodoAnalyzer& hodoAna,
                   beam::EBeamFlag beam_flag=beam::kAll);

  void DCRawHit(const TString& dcname,const RawData& rawData,
                beam::EBeamFlag beam_flag=beam::kAll);
  void DCHit(const TString& dcname, const DCAnalyzer& dcAna,
             beam::EBeamFlag beam_flag=beam::kAll);
  void BcInTracking(DCAnalyzer& dcAna, beam::EBeamFlag beam_flag=beam::kAll);
  void BcOutTracking(DCAnalyzer& dcAna, beam::EBeamFlag beam_flag=beam::kAll);
  void BcInPullExclusive(DCAnalyzer& dcAna, beam::EBeamFlag beam_flag=beam::kAll);
  void BcOutPullExclusive(DCAnalyzer& dcAna, beam::EBeamFlag beam_flag=beam::kAll);
  void D5Tracking(const D5Track& d5tr, beam::EBeamFlag beam_flag=beam::kAll);
  void D5WireResiduals(const D5Track& d5tr, Double_t z_out);

  void TriggerFlag(const RawData& rawData);
  void DAQ(const RawData& rawData);

private:
};

#endif
