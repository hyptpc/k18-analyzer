// -*- C++ -*-

#ifndef EVENT_ANALYZER_HH
#define EVENT_ANALYZER_HH

#include "DetectorID.hh"
#include "Event.hh"
#include "HistTools.hh"

class RawData;
class HodoAnalyzer;
class DCAnalyzer;

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
  void DCRawHit(const RawData& rawData, beam::EBeamFlag beam_flag=beam::kAll);
  void DCHit(const DCAnalyzer& dcAna, beam::EBeamFlag beam_flag=beam::kAll);
  void TriggerFlag(const RawData& rawData);
  void DAQ(const RawData& rawData);
};

#endif
