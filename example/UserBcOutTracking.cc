// -*- C++ -*-

#include "VEvent.hh"

#include <fstream>
#include <iostream>
#include <sstream>

#include <TString.h>

#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <UnpackerXMLReadDigit.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "DCAnalyzer.hh"
#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DCRawHit.hh"
#include "DCTdcCalibMan.hh"
#include "EventAnalyzer.hh"
#include "HistTools.hh"
#include "HodoAnalyzer.hh"
#include "HodoRawHit.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"

#define BEAMONLY 0
#define TRACKING 1
#define GLOBAL 0
#define CLUSTER 0
#define SOURCE 0
#define DEBUG 0

namespace
{
const auto& gUser = UserParamMan::GetInstance();
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();

const double mm=0.1;
const double cm=10*mm;

using seg_t = std::vector<Double_t>;
using tdc_t = std::vector<std::vector<Double_t>>;
TTree* tree;
UInt_t run_number;
UInt_t event_number;
beam::EBeamFlag beam_flag;
tdc_t trig_flag;
seg_t trig_pat;

std::map<TString, seg_t> wire;
}

//_____________________________________________________________________________
Bool_t
ProcessBegin()
{
  run_number = gUnpacker.get_run_number();
  event_number = gUnpacker.get_event_number();
  beam_flag = beam::kUnknown;
  trig_flag.clear();
  trig_pat.clear();

  for(auto& p: wire) p.second.clear();

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessNormal()
{
  using root::HF1;

  RawData rawData;
  rawData.DecodeHits();

  EventAnalyzer evAna;

  HF1("Status", 0);
  evAna.TriggerFlag(rawData);

  HF1("Status", 1);
  beam_flag = evAna.BeamFlag(rawData);

  HF1("Status", 2);
  evAna.DCRawHit("BcOut", rawData);
  evAna.DCRawHit("BcOut", rawData, beam_flag);

  DCAnalyzer dcAna(rawData);
  dcAna.DecodeRawHits();

  evAna.DCHit("BcOut", dcAna);
  evAna.DCHit("BcOut", dcAna, beam_flag);

  dcAna.TotCut("BLC2a");
  dcAna.TotCut("BLC2b");
  dcAna.DriftTimeCut("BLC2a");
  dcAna.DriftTimeCut("BLC2b");

  // dcAna.TrackSearchBcOut();

  for (Int_t plane=0; plane<NumOfLayersBcOut; ++plane) {
    for (const auto& hit : dcAna.GetBcOutHC(plane)) {
      const auto& rhit = hit->GetRawHit();
      Double_t w = hit->GetWire();
      TString name = rhit->DetectorName() + "_" + rhit->PlaneName();
      name.ToLower();
      wire[name].push_back(w);
    }
  }

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
  hist::BuildStatus();
  hist::BuildTriggerFlag();
  hist::BuildDCRaw("BcOut", true);
  hist::BuildDCHit("BcOut", true);

  tree = new TTree("bcout", "UserBcOutTracking");
  tree->Branch("run_number", &run_number);
  tree->Branch("event_number", &event_number);
  tree->Branch("beam_flag", &beam_flag, "beam_flag/I");
  tree->Branch("trig_flag", &trig_flag);
  tree->Branch("trig_pat", &trig_pat);

  const auto& gUConf = hddaq::unpacker::GConfig::get_instance();
  const auto& digit_info = gUConf.get_digit_info();
  for (const auto& name_str : DCNameList.at("BcOut")) {
    const auto name = name_str.Data();
    Int_t detector_id = digit_info.get_device_id(name);
    Int_t nplane = digit_info.get_n_plane(detector_id);
    for (Int_t plane=0; plane<nplane; ++plane) {
      TString n = (name_str + "_" +
                   digit_info.get_name_list(detector_id).at(plane));
      n.ToLower();
      tree->Branch(Form("%s_wire", n.Data()), &wire[n]);
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<HodoParamMan>("HDPRM")) &&
    (InitializeParameter<HodoPHCMan>("HDPHC")) &&
    (InitializeParameter<DCTdcCalibMan>("DCTDC")) &&
    (InitializeParameter<DCDriftParamMan>("DCDRFT")) &&
    (InitializeParameter<DCGeomMan>("DCGEO")) &&
    (InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
