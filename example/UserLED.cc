// -*- C++ -*-

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#include <TString.h>

#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "EventAnalyzer.hh"
#include "RootHelper.hh"
#include "HodoAnalyzer.hh"
#include "HodoRawHit.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "UserParamMan.hh"
#include "XTMapMan.hh"
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
tdc_t trig_flag;
seg_t trig_pat;

std::map<TString, seg_t> raw_seg;
std::map<TString, adc_t> adc_u;
std::map<TString, adc_t> adc_d;
std::map<TString, adc_t> adc_s;
std::map<TString, tdc_t> tdc_u;
std::map<TString, tdc_t> tdc_d;
std::map<TString, tdc_t> tdc_s;
std::map<TString, tdc_t> trailing_u;
std::map<TString, tdc_t> trailing_d;

std::map<TString, seg_t> hit_seg;
std::map<TString, adc_t> de_u;
std::map<TString, adc_t> de_d;
std::map<TString, adc_t> de_s;
std::map<TString, adc_t> de;
std::map<TString, tdc_t> time_u;
std::map<TString, tdc_t> time_d;
std::map<TString, tdc_t> time_s;
std::map<TString, tdc_t> mt;
std::map<TString, tdc_t> cmt;

std::map<TString, cl_t> cl_seg;
std::map<TString, cl_t> cl_de;
std::map<TString, cl_t> cl_time;
std::map<TString, cl_t> cl_tdif;
std::map<TString, cl_t> cl_size;

///// KVC
std::map<TString, adc_t> adc_a;
std::map<TString, adc_t> adc_b;
std::map<TString, adc_t> adc_c;
std::map<TString, adc_t> de_a;
std::map<TString, adc_t> de_b;
std::map<TString, adc_t> de_c;

Double_t time0;
Double_t btof0;
Double_t ftof0;
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

  for(auto& p: raw_seg) p.second.clear();
  for(auto& p: adc_u) p.second.clear();
  for(auto& p: adc_d) p.second.clear();
  for(auto& p: adc_s) p.second.clear();
  for(auto& p: tdc_u) p.second.clear();
  for(auto& p: tdc_d) p.second.clear();
  for(auto& p: tdc_s) p.second.clear();
  for(auto& p: trailing_u) p.second.clear();
  for(auto& p: trailing_d) p.second.clear();

  for(auto& p: hit_seg) p.second.clear();
  for(auto& p: de_u) p.second.clear();
  for(auto& p: de_d) p.second.clear();
  for(auto& p: de_s) p.second.clear();
  for(auto& p: de) p.second.clear();
  for(auto& p: time_u) p.second.clear();
  for(auto& p: time_d) p.second.clear();
  for(auto& p: time_s) p.second.clear();
  for(auto& p: mt) p.second.clear();
  for(auto& p: cmt) p.second.clear();

  for(auto& p: cl_seg) p.second.clear();
  for(auto& p: cl_de) p.second.clear();
  for(auto& p: cl_time) p.second.clear();
  for(auto& p: cl_tdif) p.second.clear();
  for(auto& p: cl_size) p.second.clear();

  for(auto& p: adc_a) p.second.clear();
  for(auto& p: adc_b) p.second.clear();
  for(auto& p: adc_c) p.second.clear();
  for(auto& p: de_a) p.second.clear();
  for(auto& p: de_b) p.second.clear();
  for(auto& p: de_c) p.second.clear();

  time0 = TMath::QuietNaN();
  btof0 = TMath::QuietNaN();
  ftof0 = TMath::QuietNaN();

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessNormal()
{
  using root::HF1;

  RawData rawData;
  rawData.DecodeHits("BAC");
  rawData.DecodeHits("KVC");
  
  {
    auto n = NameHodo[kBAC];
    for(const auto& hit: rawData.GetHodoRawHC(n)){      
      raw_seg[n].push_back(hit->SegmentId());
      adc_u[n].push_back(hit->GetAdcUp());
      adc_d[n].push_back(hit->GetAdcDown());
      adc_s[n].push_back(hit->GetAdcExtra());
      tdc_u[n].push_back(hit->GetArrayTdcUp());
      tdc_d[n].push_back(hit->GetArrayTdcDown());
      tdc_s[n].push_back(hit->GetArrayTdcExtra());
      // hit->Print();
    }
  }
  
  { ///// KVC
    static const TString n("KVC");
    for(const auto& hit: rawData.GetHodoRawHC(n)){
      raw_seg[n].push_back(hit->SegmentId());
      adc_a[n].push_back(hit->GetAdc(0));
      adc_b[n].push_back(hit->GetAdc(1));
      adc_c[n].push_back(hit->GetAdc(2));
      adc_d[n].push_back(hit->GetAdc(3));
      adc_s[n].push_back(hit->GetAdc(4));
      tdc_s[n].push_back(hit->GetArrayTdc(4));
      // hit->Print();
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
  hist::BuildHodoRaw(true);
  hist::BuildHodoHit(true);
  hist::BuildHodoCluster(true);

  tree = new TTree("hodo", "UserHodoscope");
  tree->Branch("run_number", &run_number);
  tree->Branch("event_number", &event_number);

  {
    Int_t ihodo = kBAC;
    auto n = NameHodo[ihodo];
    n.ToLower();
    tree->Branch(Form("%s_raw_seg", n.Data()), &raw_seg[NameHodo[ihodo]]);
    tree->Branch(Form("%s_adc_u", n.Data()), &adc_u[NameHodo[ihodo]]);
    tree->Branch(Form("%s_adc_d", n.Data()), &adc_d[NameHodo[ihodo]]);
    tree->Branch(Form("%s_adc_s", n.Data()), &adc_s[NameHodo[ihodo]]);
    tree->Branch(Form("%s_tdc_u", n.Data()), &tdc_u[NameHodo[ihodo]]);
    tree->Branch(Form("%s_tdc_d", n.Data()), &tdc_d[NameHodo[ihodo]]);
    tree->Branch(Form("%s_tdc_s", n.Data()), &tdc_s[NameHodo[ihodo]]);
  }
  { ///// KVC
    const TString n("KVC");
    const Char_t* nn = "kvc";
    tree->Branch(Form("%s_raw_seg", nn), &raw_seg[n]);
    tree->Branch(Form("%s_adc_a", nn), &adc_a[n]);
    tree->Branch(Form("%s_adc_b", nn), &adc_b[n]);
    tree->Branch(Form("%s_adc_c", nn), &adc_c[n]);
    tree->Branch(Form("%s_adc_d", nn), &adc_d[n]);
    tree->Branch(Form("%s_adc_s", nn), &adc_s[n]);
    tree->Branch(Form("%s_tdc_s", nn), &tdc_s[n]);
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
    (InitializeParameter<DCGeomMan>("DCGEO")) &&
    (InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
