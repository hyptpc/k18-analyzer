// -*- C++ -*-

#include "Event.hh"

ClassImp(HitWire);
ClassImp(CDCTrackHits);
ClassImp(CDCTrackContainer);

//_____________________________________________________________________________
EventTriggerFlag::EventTriggerFlag(const TString& name)
  : TNamed(name, name),
    trigger_flag()
{
};
