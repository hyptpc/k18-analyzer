// -*- C++ -*-

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#include <TString.h>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "UserAnalyzer.hh"
#include "RootHelper.hh"
#include "DCAnalyzer.hh"
#include "DCCluster.hh"
#include "DCHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoRawHit.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "UserParamMan.hh"
#include "XTMapMan.hh"
#include "DCTdcCalibMan.hh"
#include "BLDCWireMapMan.hh"
#include "RawData.hh"
#include "HistTools.hh"
#include "UnpackerManager.hh"

#define DEBUG 0
namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
// UInt_t run_number;
UInt_t event_number;
}

//_____________________________________________________________________________
Bool_t
ProcessBegin()
{
  event_number = gUnpacker.get_event_number();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessNormal()
{
  using root::HF1;

  RawData rawData;
  rawData.DecodeHits();

  HodoAnalyzer hodoAna(rawData);
  hodoAna.DecodeHits<FiberHit>("BHT");
  // hodoAna.TimeCut("BHT");
  hodoAna.DecodeHits("T0");

  EventAnalyzer evAna;

  HF1("Status", 0);
  evAna.TriggerFlag(rawData);

  HF1("Status", 1);

  // BeamFlag
  beam::EBeamFlag beam_flag = evAna.BeamFlag(rawData);

  evAna.HodoRawHit(rawData);
  evAna.HodoRawHit(rawData, beam_flag);

  evAna.HodoHit(hodoAna);
  evAna.HodoHit(hodoAna, beam_flag);

  HF1("Status", 20);

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
  Bool_t beam_flag = true;
  hist::BuildStatus();
  hist::BuildTriggerFlag();
  hist::BuildHodoRaw(beam_flag);
  hist::BuildHodoHit(beam_flag);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<HodoParamMan>("HDPRM")) &&
    (InitializeParameter<HodoPHCMan>("HDPHC")) &&
    (InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
