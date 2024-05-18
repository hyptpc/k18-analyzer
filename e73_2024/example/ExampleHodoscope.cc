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

using seg_t = std::vector<Double_t>;
using adc_t = std::vector<Double_t>;
using tdc_t = std::vector<std::vector<Double_t>>;
using cl_t = std::vector<Double_t>;
TTree* tree;
UInt_t run_number;
UInt_t event_number;
beam::EBeamFlag beam_flag;

std::map<TString, seg_t> raw_seg;
std::map<TString, adc_t> adc_u;
std::map<TString, adc_t> adc_d;
std::map<TString, tdc_t> tdc_u;
std::map<TString, tdc_t> tdc_d;
std::map<TString, tdc_t> trailing_u;
std::map<TString, tdc_t> trailing_d;

std::map<TString, seg_t> hit_seg;
std::map<TString, adc_t> de_u;
std::map<TString, adc_t> de_d;
std::map<TString, adc_t> de;
std::map<TString, tdc_t> time_u;
std::map<TString, tdc_t> time_d;
std::map<TString, tdc_t> mt;
std::map<TString, tdc_t> cmt;

std::map<TString, cl_t> cl_seg;
std::map<TString, cl_t> cl_de;
std::map<TString, cl_t> cl_time;
std::map<TString, cl_t> cl_size;

Double_t time0;
Double_t btof0;
}

//_____________________________________________________________________________
Bool_t
ProcessBegin()
{
  run_number = gUnpacker.get_run_number();
  event_number = gUnpacker.get_event_number();
  beam_flag = beam::kUnknown;

  for(auto& p: raw_seg) p.second.clear();
  for(auto& p: adc_u) p.second.clear();
  for(auto& p: adc_d) p.second.clear();
  for(auto& p: tdc_u) p.second.clear();
  for(auto& p: tdc_d) p.second.clear();
  for(auto& p: trailing_u) p.second.clear();
  for(auto& p: trailing_d) p.second.clear();

  for(auto& p: hit_seg) p.second.clear();
  for(auto& p: de_u) p.second.clear();
  for(auto& p: de_d) p.second.clear();
  for(auto& p: de) p.second.clear();
  for(auto& p: time_u) p.second.clear();
  for(auto& p: time_d) p.second.clear();
  for(auto& p: mt) p.second.clear();
  for(auto& p: cmt) p.second.clear();

  for(auto& p: cl_seg) p.second.clear();
  for(auto& p: cl_de) p.second.clear();
  for(auto& p: cl_time) p.second.clear();
  for(auto& p: cl_size) p.second.clear();

  time0 = TMath::QuietNaN();
  btof0 = TMath::QuietNaN();

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
  beam_flag = evAna.BeamFlag(rawData);

  evAna.HodoRawHit(rawData);
  evAna.HodoRawHit(rawData, beam_flag);

  HF1("Status", 2);

  evAna.HodoHit(hodoAna);
  evAna.HodoHit(hodoAna, beam_flag);

  HF1("Status", 3);

  evAna.HodoCluster(hodoAna);
  evAna.HodoCluster(hodoAna, beam_flag);

  HF1("Status", 4);

  for(const auto& hit: rawData.GetHodoRawHC("BHT")){
    raw_seg["BHT"].push_back(hit->SegmentId());
    tdc_u["BHT"].push_back(hit->GetArrayTdcUp());
    tdc_d["BHT"].push_back(hit->GetArrayTdcDown());
    trailing_u["BHT"].push_back(hit->GetArrayTdcTrailing(0));
    trailing_d["BHT"].push_back(hit->GetArrayTdcTrailing(1));
  }

  for(Int_t ihodo=kT1; ihodo<kNumHodo + 1; ++ihodo){
    auto n = (ihodo < kNumHodo) ? NameHodo[ihodo] : "AC";
    for(const auto& hit: rawData.GetHodoRawHC(n)){
      raw_seg[n].push_back(hit->SegmentId());
      adc_u[n].push_back(hit->GetAdcUp());
      adc_d[n].push_back(hit->GetAdcDown());
      tdc_u[n].push_back(hit->GetArrayTdcUp());
      tdc_d[n].push_back(hit->GetArrayTdcDown());
    }
  }

  HF1("Status", 5);

  for(Int_t ihodo=kBHT; ihodo<kNumHodo; ++ihodo){
    auto n = NameHodo[ihodo];
    for(Int_t i=0, nh=hodoAna.GetNHits(n); i<nh; ++i){
      const auto& hit = hodoAna.GetHit(n, i);
      hit_seg[n].push_back(hit->SegmentId());
      de_u[n].push_back(hit->GetAUp());
      de_d[n].push_back(hit->GetADown());
      de[n].push_back(hit->DeltaE());
      time_u[n].push_back(hit->GetArrayTime(0));
      time_d[n].push_back(hit->GetArrayTime(1));
      mt[n].push_back(hit->GetArrayMeanTime());
      cmt[n].push_back(hit->GetArrayCMeanTime());
    }
  }

  HF1("Status", 6);

  for(Int_t ihodo=kBHT; ihodo<kNumHodo; ++ihodo){
    auto n = NameHodo[ihodo];
    for(Int_t i=0, nh=hodoAna.GetNClusters(n); i<nh; ++i){
      const auto& cl = hodoAna.GetCluster(n, i);
      cl_seg[n].push_back(cl->MeanSeg());
      cl_de[n].push_back(cl->DeltaE());
      cl_time[n].push_back(cl->CTime());
      cl_size[n].push_back(cl->ClusterSize());
    }
  }

  time0 = hodoAna.Time0();
  btof0 = hodoAna.Btof0();

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
  hist::BuildStatus();
  hist::BuildTriggerFlag();
  hist::BuildHodoRaw(true);
  hist::BuildHodoHit(true);
  hist::BuildHodoCluster(true);

  tree = new TTree("hodo", "UserHodoscope");
  tree->Branch("run_number", &run_number);
  tree->Branch("event_number", &event_number);
  tree->Branch("beam_flag", &beam_flag, "beam_flag/I");
  tree->Branch("bht_raw_seg", &raw_seg["BHT"]);
  tree->Branch("bht_tdc_u", &tdc_u["BHT"]);
  tree->Branch("bht_tdc_d", &tdc_d["BHT"]);
  tree->Branch("bht_trailing_u", &trailing_u["BHT"]);
  tree->Branch("bht_trailing_d", &trailing_d["BHT"]);
  for(Int_t ihodo=kT1; ihodo<kNumHodo; ++ihodo){
    auto n = NameHodo[ihodo];
    n.ToLower();
    tree->Branch(Form("%s_raw_seg", n.Data()), &raw_seg[NameHodo[ihodo]]);
    tree->Branch(Form("%s_adc_u", n.Data()), &adc_u[NameHodo[ihodo]]);
    tree->Branch(Form("%s_adc_d", n.Data()), &adc_d[NameHodo[ihodo]]);
    tree->Branch(Form("%s_tdc_u", n.Data()), &tdc_u[NameHodo[ihodo]]);
    tree->Branch(Form("%s_tdc_d", n.Data()), &tdc_d[NameHodo[ihodo]]);
  }
  tree->Branch("bac_raw_seg", &raw_seg["AC"]);
  tree->Branch("bac_adc", &adc_u["AC"]);
  tree->Branch("bac_tdc", &tdc_u["AC"]);

  for(Int_t ihodo=kBHT; ihodo<kNumHodo; ++ihodo){
    auto n = NameHodo[ihodo];
    n.ToLower();
    tree->Branch(Form("%s_hit_seg", n.Data()), &hit_seg[NameHodo[ihodo]]);
    tree->Branch(Form("%s_de_u", n.Data()), &de_u[NameHodo[ihodo]]);
    tree->Branch(Form("%s_de_d", n.Data()), &de_d[NameHodo[ihodo]]);
    tree->Branch(Form("%s_de", n.Data()), &de[NameHodo[ihodo]]);
    tree->Branch(Form("%s_time_u", n.Data()), &time_u[NameHodo[ihodo]]);
    tree->Branch(Form("%s_time_d", n.Data()), &time_d[NameHodo[ihodo]]);
    tree->Branch(Form("%s_mt", n.Data()), &mt[NameHodo[ihodo]]);
    tree->Branch(Form("%s_cmt", n.Data()), &cmt[NameHodo[ihodo]]);
  }

  for(Int_t ihodo=kBHT; ihodo<kNumHodo; ++ihodo){
    auto n = NameHodo[ihodo];
    n.ToLower();
    tree->Branch(Form("%s_cl_seg", n.Data()), &cl_seg[NameHodo[ihodo]]);
    tree->Branch(Form("%s_cl_de", n.Data()), &cl_de[NameHodo[ihodo]]);
    tree->Branch(Form("%s_cl_time", n.Data()), &cl_time[NameHodo[ihodo]]);
    tree->Branch(Form("%s_cl_size", n.Data()), &cl_size[NameHodo[ihodo]]);
  }

  tree->Branch("time0", &time0);
  tree->Branch("btof0", &btof0);

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
