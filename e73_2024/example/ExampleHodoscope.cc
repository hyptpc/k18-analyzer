// -*- C++ -*-

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#include <TString.h>

#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "DetectorID.hh"
#include "EventAnalyzer.hh"
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

using seg_t = std::vector<Int_t>;
using adc_t = std::vector<Double_t>;
using tdc_t = std::vector<std::vector<Double_t>>;
TTree* tree;
UInt_t run_number;
UInt_t event_number;

seg_t t0_raw_seg;
tdc_t t0_tdc_u;
tdc_t t0_tdc_d;
adc_t t0_adc_u;
adc_t t0_adc_d;
}

//_____________________________________________________________________________
Bool_t
ProcessBegin()
{
  run_number = gUnpacker.get_run_number();
  event_number = gUnpacker.get_event_number();
  t0_raw_seg.clear();
  t0_tdc_u.clear();
  t0_tdc_d.clear();
  t0_adc_u.clear();
  t0_adc_d.clear();
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
  hodoAna.DecodeHits<BH2Hit>("T0");

  EventAnalyzer evAna;

  HF1("Status", 0);
  evAna.TriggerFlag(rawData);

  HF1("Status", 1);

  // BeamFlag
  beam::EBeamFlag beam_flag = evAna.BeamFlag(rawData);

  evAna.HodoRawHit(rawData);
  evAna.HodoRawHit(rawData, beam_flag);

  HF1("Status", 2);

  evAna.HodoHit(hodoAna);
  evAna.HodoHit(hodoAna, beam_flag);

  HF1("Status", 3);

  evAna.HodoCluster(hodoAna);
  evAna.HodoCluster(hodoAna, beam_flag);

  HF1("Status", 4);

  for(const auto& hit: rawData.GetHodoRawHC("T0")){
    t0_raw_seg.push_back(hit->SegmentId());
    t0_tdc_u.push_back(hit->GetArrayTdcUp());
    t0_tdc_d.push_back(hit->GetArrayTdcDown());
    t0_adc_u.push_back(hit->GetArrayAdcUp().front());
    t0_adc_d.push_back(hit->GetArrayAdcDown().front());
  }

  HF1("Status", 20);

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessEnd()
{
  tree->Fill();
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
  hist::BuildHodoCluster(beam_flag);

  tree = new TTree("hodo", "UserHodoscope");
  tree->Branch("run_number", &run_number);
  tree->Branch("event_number", &event_number);
  tree->Branch("t0_raw_seg", &t0_raw_seg);
  tree->Branch("t0_adc_u", &t0_adc_u);
  tree->Branch("t0_adc_d", &t0_adc_d);
  tree->Branch("t0_tdc_u", &t0_tdc_u);
  tree->Branch("t0_tdc_d", &t0_tdc_d);

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
