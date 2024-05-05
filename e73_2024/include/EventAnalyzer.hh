// -*- C++ -*-

#ifndef EVENT_ANALYZER_HH
#define EVENT_ANALYZER_HH

#include "Event.hh"
#include "HistTools.hh"

class RawData;

//_____________________________________________________________________________
class EventAnalyzer
{
public:
  EventAnalyzer();
  ~EventAnalyzer();

private:

public:
  void HodoRawHit(const RawData& rawData, hist::EParticle particle=hist::kAll);
  void DCRawHit(const RawData& rawData, hist::EParticle particle=hist::kAll);
  void TriggerFlag(const RawData& rawData);
  void DAQ(const RawData& rawData);
};

#endif
