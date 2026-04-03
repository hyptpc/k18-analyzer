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

#define JPARC2025Nov 0 // 1: Old runs (Skip T2+), 0: 2026Apr (T2+)

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
#if JPARC2025Nov
  // Skip new detectors to avoid unpacker errors in old runs
  rawData.DecodeHits("TriggerFlag");
  for(Int_t ihodo=kBHT; ihodo<kNumHodo; ++ihodo){
    if (ihodo >= kT2) continue;
    rawData.DecodeHits(NameHodo[ihodo]);
  }
  rawData.DecodeHits("BLC1a");
  rawData.DecodeHits("BLC1b");
  rawData.DecodeHits("BLC2a");
  rawData.DecodeHits("BLC2b");
#else
  rawData.DecodeHits();
#endif

  EventAnalyzer evAna;

  root::HF1("Status", 0);
  evAna.TriggerFlag(rawData);

  root::HF1("Status", 1);
  evAna.HodoRawHit(rawData);

  root::HF1("Status", 2);
  evAna.DCRawHit("BcIn", rawData);
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
  hist::BuildDCRaw("BcIn");
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
