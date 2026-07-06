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

#define HodoCut 0 // with BHT/BH2
#define TimeCut 1 // in cluster analysis
#define Chi2Cut 1 // for BcOut tracking
#define BH2MatchCut 1

#define KmBeam 0 // BeamA event selection
#define JPARC2025Nov 0 // 1: Old runs (Skip T2+), 0: 2026Apr (T2+)

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

//BHT - BH2 
Double_t time0;
Double_t btof0;

// BcOut
std::map<TString, seg_t> wireBcOut;
Int_t ntBcOut;
std::vector<Double_t> chisqrBcOut;
std::vector<Double_t> x0BcOut;
std::vector<Double_t> y0BcOut;
std::vector<Double_t> u0BcOut;
std::vector<Double_t> v0BcOut;

//VP HS
std::vector<std::vector<Double_t>> xvpHS;
std::vector<std::vector<Double_t>> yvpHS;
std::vector<std::vector<Double_t>> zvpHS;
std::vector<std::vector<Double_t>> uvpHS;
std::vector<std::vector<Double_t>> vvpHS;
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

  for(auto& p: event.wireBcOut) p.second.clear();
  event.ntBcOut = 0;
  event.chisqrBcOut.clear();
  event.x0BcOut.clear();
  event.y0BcOut.clear();
  event.u0BcOut.clear();
  event.v0BcOut.clear();

  event.xvpHS.clear();
  event.yvpHS.clear();
  event.zvpHS.clear();
  event.uvpHS.clear();
  event.vvpHS.clear();

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessNormal()
{
  using root::HF1;

  RawData rawData;
  
  for(Int_t ihodo=kBHT; ihodo<kNumHodo; ++ihodo){
#if JPARC2025Nov
    if (ihodo >= kT2) continue;
#endif
    auto n = NameHodo[ihodo];
    rawData.DecodeHits(n);
  }
  for (const auto& name : DCNameList.at("BcOut")) rawData.DecodeHits(name);

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
  evAna.DCRawHit("BcOut", rawData);
  evAna.DCRawHit("BcOut", rawData, event.beam_flag);

  DCAnalyzer dcAna(rawData);
  dcAna.DecodeBcOutHits();

  dcAna.TotCut("BLC2a");
  dcAna.TotCut("BLC2b");
  dcAna.DriftTimeCut("BLC2a");
  dcAna.DriftTimeCut("BLC2b");
  evAna.DCHit("BcOut", dcAna);
  evAna.DCHit("BcOut", dcAna, event.beam_flag);

  for (Int_t plane=0; plane<NumOfLayersBcOut; ++plane) {
    for (const auto& hit : dcAna.GetBcOutHC(plane)) {
      const auto& rhit = hit->GetRawHit();
      Double_t w = hit->GetWire();
      TString name = rhit->DetectorName() + "_" + rhit->PlaneName();
      name.ToLower();
      event.wireBcOut[name].push_back(w);
    }
  }

  dcAna.TrackSearchBcOut();
  evAna.BcOutTracking(dcAna);
  evAna.BcOutTracking(dcAna, event.beam_flag);

  for(const auto& track : dcAna.GetBcOutTrackContainer()){
    track->Print();
    event.ntBcOut++;
    event.chisqrBcOut.push_back(track->GetChiSquare());
    event.x0BcOut.push_back(track->GetX0());
    event.y0BcOut.push_back(track->GetY0());
    event.u0BcOut.push_back(track->GetU0());
    event.v0BcOut.push_back(track->GetV0());
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
  hist::BuildDCTrack("BcOut", true);

  tree = new TTree("bcout", "UserBcOutTracking");
  tree->Branch("run_number", &event.run_number);
  tree->Branch("event_number", &event.event_number);
  tree->Branch("beam_flag", &event.beam_flag, "beam_flag/I");
  tree->Branch("trig_flag", &event.trig_flag);
  tree->Branch("trig_pat", &event.trig_pat);

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
      tree->Branch(Form("%s_wire", n.Data()), &event.wireBcOut[n]);
    }
  }
  tree->Branch("ntBcOut", &event.ntBcOut);
  tree->Branch("chisqrBcOut", &event.chisqrBcOut);
  tree->Branch("x0BcOut", &event.x0BcOut);
  tree->Branch("y0BcOut", &event.y0BcOut);
  tree->Branch("u0BcOut", &event.u0BcOut);
  tree->Branch("v0BcOut", &event.v0BcOut);

  tree->Branch("xvpHS", &event.xvpHS);
  tree->Branch("yvpHS", &event.yvpHS);
  tree->Branch("zvpHS", &event.zvpHS);
  tree->Branch("uvpHS", &event.uvpHS);
  tree->Branch("vvpHS", &event.vvpHS);

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
