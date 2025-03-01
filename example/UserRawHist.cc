// -*- C++ -*-

#include <iostream>
#include <sstream>
#include <cmath>

#include <TString.h>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "EventAnalyzer.hh"
#include "HistTools.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "VEvent.hh"
#include "UserParamMan.hh"

//_____________________________________________________________________________
Bool_t
ProcessBegin()
{
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessNormal()
{
  RawData rawData;
  rawData.DecodeHits();

  EventAnalyzer evAna;

  root::HF1("Status", 0);
  evAna.TriggerFlag(rawData);

  root::HF1("Status", 1);
  evAna.HodoRawHit(rawData);

  root::HF1("Status", 2);
  evAna.DCRawHit("BcOut", rawData);

  root::HF1("Status", 3);
  evAna.DAQ(rawData);

  root::HF1("Status", 20);

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessEnd()
{
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  hist::BuildStatus();
  hist::BuildTriggerFlag();
  hist::BuildDAQ();
  hist::BuildHodoRaw();
  hist::BuildDCRaw("BcOut");
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    InitializeParameter<UserParamMan>("USER") &&
    InitializeParameter<DCGeomMan>("DCGEO");
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
