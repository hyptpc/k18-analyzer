// -*- C++ -*-

#include "VEvent.hh"

#include <fstream>
#include <iostream>
#include <sstream>

#include <TString.h>

#include "ConfMan.hh"
#include "DCAnalyzer.hh"
#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DCLocalTrack.hh"
#include "DCRawHit.hh"
#include "HodoRawHit.hh"
#include "DCTdcCalibMan.hh"
#include "DetectorID.hh"
#include "EventAnalyzer.hh"
#include "HistTools.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"

#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <UnpackerXMLReadDigit.hh>

#define FillPullExclusive 1

namespace
{
const auto& gUser = UserParamMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();

using seg_t = std::vector<Double_t>;
using tdc_t = std::vector<std::vector<Double_t>>;
TTree* tree;

struct Event
{
UInt_t run_number;
UInt_t event_number;
beam::EBeamFlag beam_flag;
tdc_t trig_flag;
seg_t trig_pat;

std::map<TString, seg_t> wire;

Int_t ntrack;
std::vector<Double_t> chisqr;
std::vector<Double_t> x0;
std::vector<Double_t> y0;
std::vector<Double_t> u0;
std::vector<Double_t> v0;
};
Event event;
}

//_____________________________________________________________________________
Bool_t
ProcessBegin()
{
  event.run_number = gUnpacker.get_run_number();
  event.event_number = gUnpacker.get_event_number();
  event.beam_flag = beam::kUnknown;
  event.trig_flag.clear();
  event.trig_pat.clear();

  for(auto& p: event.wire) p.second.clear();
  event.ntrack = 0;
  event.chisqr.clear();
  event.x0.clear();
  event.y0.clear();
  event.u0.clear();
  event.v0.clear();

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessNormal()
{
  using root::HF1;

  RawData rawData;
  for (const auto& name : DCNameList.at("BcIn")) rawData.DecodeHits(name);

  EventAnalyzer evAna;

  HF1("Status", 0);
  rawData.DecodeHits("TriggerFlag");
  evAna.TriggerFlag(rawData);

  for(const auto& hit: rawData.GetHodoRawHC("TriggerFlag")){
    event.trig_flag.push_back(hit->GetArrayTdc());
    event.trig_pat.push_back(hit->SegmentId());
  }

  HF1("Status", 1);
  rawData.DecodeHits("BAC"); // for beam_flag
  rawData.DecodeHits("BHT"); // for beam_flag
  event.beam_flag = evAna.BeamFlag(rawData);

  HF1("Status", 2);
  evAna.DCRawHit("BcIn", rawData);
  evAna.DCRawHit("BcIn", rawData, event.beam_flag);

  DCAnalyzer dcAna(rawData);
  dcAna.DecodeBcInHits();

  dcAna.TotCut("BLC1a");
  dcAna.TotCut("BLC1b");
  dcAna.DriftTimeCut("BLC1a");
  dcAna.DriftTimeCut("BLC1b");
  evAna.DCHit("BcIn", dcAna);
  evAna.DCHit("BcIn", dcAna, event.beam_flag);

  for (Int_t plane=0; plane<NumOfLayersBcIn; ++plane) {
    for (const auto& hit : dcAna.GetBcInHC(plane)) {
      const auto& rhit = hit->GetRawHit();
      Double_t w = hit->GetWire();
      TString name = rhit->DetectorName() + "_" + rhit->PlaneName();
      name.ToLower();
      event.wire[name].push_back(w);
    }
  }

  dcAna.TrackSearchBcIn();
  evAna.BcInTracking(dcAna);
  evAna.BcInTracking(dcAna, event.beam_flag);
#if FillPullExclusive
  evAna.BcInPullExclusive(dcAna, event.beam_flag);
#endif

  for(const auto& track : dcAna.GetBcInTrackContainer()){
    track->Print();
    event.ntrack++;
    event.chisqr.push_back(track->GetChiSquare());
    event.x0.push_back(track->GetX0());
    event.y0.push_back(track->GetY0());
    event.u0.push_back(track->GetU0());
    event.v0.push_back(track->GetV0());
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
  hist::BuildDCRaw("BcIn", true);
  hist::BuildDCHit("BcIn", true);
  hist::BuildDCTrack("BcIn", true);

  tree = new TTree("bcin", "UserBcInTracking");
  tree->Branch("run_number", &event.run_number);
  tree->Branch("event_number", &event.event_number);
  tree->Branch("beam_flag", &event.beam_flag, "beam_flag/I");
  tree->Branch("trig_flag", &event.trig_flag);
  tree->Branch("trig_pat", &event.trig_pat);

  const auto& gUConf = hddaq::unpacker::GConfig::get_instance();
  const auto& digit_info = gUConf.get_digit_info();
  for (const auto& name_str : DCNameList.at("BcIn")) {
    const auto name = name_str.Data();
    Int_t detector_id = digit_info.get_device_id(name);
    Int_t nplane = digit_info.get_n_plane(detector_id);
    for (Int_t plane=0; plane<nplane; ++plane) {
      TString n = (name_str + "_" +
                   digit_info.get_name_list(detector_id).at(plane));
      n.ToLower();
      tree->Branch(Form("%s_wire", n.Data()), &event.wire[n]);
    }
  }
  tree->Branch("ntrack", &event.ntrack);
  tree->Branch("chisqr", &event.chisqr);
  tree->Branch("x0", &event.x0);
  tree->Branch("y0", &event.y0);
  tree->Branch("u0", &event.u0);
  tree->Branch("v0", &event.v0);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
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
