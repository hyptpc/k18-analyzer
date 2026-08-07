// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iostream>
#include <map>
#include <vector>

#include <TString.h>

#include "BH2Hit.hh"
#include "CherenkovHit.hh"
#include "ConfMan.hh"
#include "DetectorID.hh"
#include "EventAnalyzer.hh"
#include "FiberHit.hh"
#include "HistTools.hh"
#include "HodoAnalyzer.hh"
#include "HodoHit.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"

#include "UnpackerManager.hh"

#define DEBUG 0
#define JPARC2025Nov 0 // 1: Old runs (Skip T2+), 0: 2026Apr (T2+)

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
std::map<TString, adc_t> npe_sum_online;
std::map<TString, adc_t> npe_sum_offline;
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

// KVC
std::map<TString, adc_t> adc_a;
std::map<TString, adc_t> adc_b;
std::map<TString, adc_t> adc_c;
std::map<TString, adc_t> de_a;
std::map<TString, adc_t> de_b;
std::map<TString, adc_t> de_c;

Double_t time0;
Double_t time0_seg;
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
  for(auto& p: npe_sum_online) p.second.clear();
  for(auto& p: npe_sum_offline) p.second.clear();
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
  time0_seg = TMath::QuietNaN();
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
  for(Int_t ihodo=kBHT; ihodo<kNumHodo; ++ihodo){
#if JPARC2025Nov
    if (ihodo >= kT2) continue;
#endif
    auto n = NameHodo[ihodo];
    rawData.DecodeHits(n);
  }
  
  // hodoAna.DecodeHits<T>(name, makeCluster = true);
  HodoAnalyzer hodoAna(rawData);
  hodoAna.DecodeHits<FiberHit>("BHT");
  hodoAna.DecodeHits<BH2Hit>("BH2");
  for(Int_t ihodo=kBAC; ihodo<kNumHodo; ++ihodo){
#if JPARC2025Nov
    if (ihodo >= kT2) continue;
#endif
    auto n = NameHodo[ihodo];
    Bool_t do_not_cluster = !HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::NoCluster);
    if (HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::Cherenkov))
      hodoAna.DecodeHits<CherenkovHit>(n, do_not_cluster);
    else
      hodoAna.DecodeHits(n, do_not_cluster);
  }

  EventAnalyzer evAna;

  HF1("Status", 0);
  rawData.DecodeHits("TriggerFlag");
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

  for(const auto& hit: rawData.GetHodoRawHC("TriggerFlag")){
    trig_flag.push_back(hit->GetArrayTdcUp());
    trig_pat.push_back(hit->SegmentId());
  }

  for(const auto& hit: rawData.GetHodoRawHC("BHT")){
    raw_seg["BHT"].push_back(hit->SegmentId());
    tdc_u["BHT"].push_back(hit->GetArrayTdcUp());
    tdc_d["BHT"].push_back(hit->GetArrayTdcDown());
    trailing_u["BHT"].push_back(hit->GetArrayTdcTrailing(0));
    trailing_d["BHT"].push_back(hit->GetArrayTdcTrailing(1));
  }
  
  for(Int_t ihodo=kBH2; ihodo<kNumHodo; ++ihodo){
#if JPARC2025Nov
    if (ihodo >= kT2) continue;
#endif
    if (ihodo == kKVC) continue;
    auto n = NameHodo[ihodo];
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
  
  { // KVC
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

  HF1("Status", 5);

  for(Int_t ihodo=kBHT; ihodo<kNumHodo; ++ihodo){
#if JPARC2025Nov
    if (ihodo >= kT2) continue;
#endif
    if (ihodo == kBAC || ihodo == kKVC) continue;

    auto n = NameHodo[ihodo];
    for(Int_t i=0, nh=hodoAna.GetNHits(n); i<nh; ++i){
      const auto& hit = hodoAna.GetHit(n, i);
      auto n_ch = hit->NumOfChannel();
      hit_seg[n].push_back(hit->SegmentId());
      if (!HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::NoADC)) {
        de_u[n].push_back(hit->GetAUp());
        de[n].push_back(hit->DeltaE());
      }
      time_u[n].push_back(hit->GetArrayTime(0));
      mt[n].push_back(hit->GetArrayMeanTime());
      cmt[n].push_back(hit->GetArrayCMeanTime());
      if (n_ch > 1) {
        de_d[n].push_back(hit->GetADown());
        time_d[n].push_back(hit->GetArrayTime(1));
      }
      if (ihodo == kHTOF) {
        de_s[n].push_back(hit->GetAExtra());
        time_s[n].push_back(hit->GetArrayTime(2));
      }
    }
  }

  { // BAC
    const TString n = "BAC";
    for (auto x : hodoAna.GetOfflineNpe(n))
      npe_sum_offline[n].push_back(x);

    for(Int_t i=0, nh=hodoAna.GetNHits(n); i<nh; ++i){
      const auto* hit = hodoAna.GetHit<CherenkovHit>(n, i);
      auto seg = hit->SegmentId();
      hit_seg[n].push_back(seg);
      de_u[n].push_back(hit->Npe());
      Double_t onsum = hit->NpeSum(0);
      if (!TMath::IsNaN(onsum)) npe_sum_online[n].push_back(onsum);
      mt[n].push_back(hit->GetArrayCTime(HodoRawHit::kUp));
      cmt[n].push_back(hit->GetArrayCTime(HodoRawHit::kUp));
    }
  }

  { // KVC
    const TString n = "KVC";
    for (auto x : hodoAna.GetOfflineNpe(n))
      npe_sum_offline[n].push_back(x);

    for(Int_t i=0, nh=hodoAna.GetNHits(n); i<nh; ++i){
      const auto* hit = hodoAna.GetHit<CherenkovHit>(n, i);
      hit_seg[n].push_back(hit->SegmentId());
      de_a[n].push_back(hit->GetNpe(HodoRawHit::EChannelKVC::kA, 0));
      de_b[n].push_back(hit->GetNpe(HodoRawHit::EChannelKVC::kB, 0));
      de_c[n].push_back(hit->GetNpe(HodoRawHit::EChannelKVC::kC, 0));
      de_d[n].push_back(hit->GetNpe(HodoRawHit::EChannelKVC::kD, 0));
      npe_sum_online[n].push_back(hit->NpeSum(0));
      mt[n].push_back(hit->GetArrayCTime(HodoRawHit::kExtra));
      cmt[n].push_back(hit->GetArrayCTime(HodoRawHit::kExtra));
    }
  }

  HF1("Status", 6);

  for(Int_t ihodo=kBHT; ihodo<kNumHodo; ++ihodo){
#if JPARC2025Nov
    if (ihodo >= kT2) continue;
#endif
    if (HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::NoCluster)) continue;
    auto n = NameHodo[ihodo];
    for(Int_t i=0, nh=hodoAna.GetNClusters(n); i<nh; ++i){
      const auto& cl = hodoAna.GetCluster(n, i);
      cl_seg[n].push_back(cl->MeanSeg());
      cl_de[n].push_back(cl->DeltaE());
      cl_time[n].push_back(cl->CTime());
      cl_tdif[n].push_back(cl->TimeDiff());
      cl_size[n].push_back(cl->ClusterSize());
    }
  }

  time0 = hodoAna.Time0();
  time0_seg = hodoAna.Time0Seg();
  btof0 = hodoAna.Btof0();
  ftof0 = hodoAna.Ftof0();

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

  struct HodoFlags {
    UInt_t mask;
    Bool_t has_adc;
    Bool_t two_side;
    Bool_t has_sum; // for HTOF and KVC
    Bool_t no_cluster;
  };

  auto GetHodoFlags = [&](Int_t ihodo) -> HodoFlags {
    const UInt_t mask = HodoGroupMask[ihodo];
    HodoFlags f;
    f.mask        = mask;
    f.has_adc     = !HasHodoGroup(mask, HodoGroup::NoADC);
    f.two_side    = !HasHodoGroup(mask, HodoGroup::OneSideReadout);
    f.has_sum     = (ihodo == kHTOF || ihodo == kKVC);
    f.no_cluster  = HasHodoGroup(mask, HodoGroup::NoCluster);
    return f;
  };

  tree = new TTree("hodo", "UserHodoscope");
  tree->Branch("run_number", &run_number);
  tree->Branch("event_number", &event_number);
  tree->Branch("beam_flag", &beam_flag, "beam_flag/I");
  tree->Branch("trig_flag", &trig_flag);
  tree->Branch("trig_pat", &trig_pat);

  tree->Branch("bht_raw_seg", &raw_seg["BHT"]);
  tree->Branch("bht_tdc_u", &tdc_u["BHT"]);
  tree->Branch("bht_tdc_d", &tdc_d["BHT"]);
  tree->Branch("bht_trailing_u", &trailing_u["BHT"]);
  tree->Branch("bht_trailing_d", &trailing_d["BHT"]);

  for(Int_t ihodo=kBH2; ihodo<kNumHodo; ++ihodo){
#if JPARC2025Nov
    if (ihodo >= kT2) continue;
#endif
    if (ihodo == kKVC) continue;
    auto n = NameHodo[ihodo];
    n.ToLower();
    const auto f = GetHodoFlags(ihodo);
    tree->Branch(Form("%s_raw_seg", n.Data()), &raw_seg[NameHodo[ihodo]]);
    if (f.has_adc) {
      tree->Branch(Form("%s_adc_u", n.Data()), &adc_u[NameHodo[ihodo]]);
      if (f.two_side) {
        tree->Branch(Form("%s_adc_d", n.Data()), &adc_d[NameHodo[ihodo]]);
        if (f.has_sum) 
          tree->Branch(Form("%s_adc_s", n.Data()), &adc_s[NameHodo[ihodo]]);
      }
    }
    tree->Branch(Form("%s_tdc_u", n.Data()), &tdc_u[NameHodo[ihodo]]);
    if (f.two_side) {
      tree->Branch(Form("%s_tdc_d", n.Data()), &tdc_d[NameHodo[ihodo]]);
      if (f.has_sum || ihodo == kBH2) 
        tree->Branch(Form("%s_tdc_s", n.Data()), &tdc_s[NameHodo[ihodo]]);
    }
  }
  { // KVC (ch 0–3: indiv a,b,c,d; ch 4: SUM)
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

  for(Int_t ihodo=kBHT; ihodo<kNumHodo; ++ihodo){
#if JPARC2025Nov
    if (ihodo >= kT2) continue;
#endif
    if (ihodo == kBAC || ihodo == kKVC) continue;
    auto n = NameHodo[ihodo];
    n.ToLower();
    const auto f = GetHodoFlags(ihodo);
    Bool_t is_cherenkov = HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::Cherenkov);
    const Char_t* dex = is_cherenkov ? "npe" : "de";
    tree->Branch(Form("%s_hit_seg", n.Data()), &hit_seg[NameHodo[ihodo]]);
    if (f.has_adc) {
      tree->Branch(Form("%s_%s_u", n.Data(), dex), &de_u[NameHodo[ihodo]]);
      if (f.two_side) {
        tree->Branch(Form("%s_%s_d", n.Data(), dex), &de_d[NameHodo[ihodo]]);
        if (f.has_sum) 
          tree->Branch(Form("%s_%s_s", n.Data(), dex), &de_s[NameHodo[ihodo]]);
      }
      tree->Branch(Form("%s_%s", n.Data(), dex), &de[NameHodo[ihodo]]);
    }
    tree->Branch(Form("%s_time_u", n.Data()), &time_u[NameHodo[ihodo]]);
    if (f.two_side) {
      tree->Branch(Form("%s_time_d", n.Data()), &time_d[NameHodo[ihodo]]);
      if (f.has_sum || ihodo == kBH2)
        tree->Branch(Form("%s_time_s", n.Data()), &time_s[NameHodo[ihodo]]);
    }
    tree->Branch(Form("%s_mt", n.Data()), &mt[NameHodo[ihodo]]);
    tree->Branch(Form("%s_cmt", n.Data()), &cmt[NameHodo[ihodo]]);
  }
  { // BAC: npe_u, npe_sum_online (seg4, per hit). npe_sum_offline: 1 (event, seg0–3 raw). hit_seg/mt/cmt: per hit, index i.
    const TString n("BAC");
    const Char_t* nn = "bac";
    tree->Branch(Form("%s_hit_seg", nn), &hit_seg[n]);
    tree->Branch(Form("%s_npe_u", nn), &de_u[n]);
    tree->Branch(Form("%s_npe_sum_online", nn), &npe_sum_online[n]);
    tree->Branch(Form("%s_npe_sum_offline", nn), &npe_sum_offline[n]);
    tree->Branch(Form("%s_mt", nn), &mt[n]);
    tree->Branch(Form("%s_cmt", nn), &cmt[n]);
  }
  { // KVC: npe_a..d, npe_sum_online (per hit). npe_sum_offline: 8, [seg]=seg id; for hit i use [hit_seg[i]]. hit_seg/mt/cmt: per hit.
    const TString n("KVC");
    const Char_t* nn = "kvc";
    tree->Branch(Form("%s_hit_seg", nn), &hit_seg[n]);
    tree->Branch(Form("%s_npe_a", nn), &de_a[n]);
    tree->Branch(Form("%s_npe_b", nn), &de_b[n]);
    tree->Branch(Form("%s_npe_c", nn), &de_c[n]);
    tree->Branch(Form("%s_npe_d", nn), &de_d[n]);
    tree->Branch(Form("%s_npe_sum_online", nn), &npe_sum_online[n]);
    tree->Branch(Form("%s_npe_sum_offline", nn), &npe_sum_offline[n]);
    tree->Branch(Form("%s_mt", nn), &mt[n]);
    tree->Branch(Form("%s_cmt", nn), &cmt[n]);
  }

  for(Int_t ihodo=kBHT; ihodo<kNumHodo; ++ihodo){
#if !JPARC2026Apr
    if (ihodo >= kT2) continue;
#endif
    const auto f = GetHodoFlags(ihodo);
    if (f.no_cluster) continue;
    auto n = NameHodo[ihodo];
    n.ToLower();
    const Char_t* cldex = HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::Cherenkov) ? "npe" : "de";
    tree->Branch(Form("%s_cl_seg", n.Data()), &cl_seg[NameHodo[ihodo]]);
    tree->Branch(Form("%s_cl_%s", n.Data(), cldex), &cl_de[NameHodo[ihodo]]);
    tree->Branch(Form("%s_cl_time", n.Data()), &cl_time[NameHodo[ihodo]]);
    tree->Branch(Form("%s_cl_tdif", n.Data()), &cl_tdif[NameHodo[ihodo]]);
    tree->Branch(Form("%s_cl_size", n.Data()), &cl_size[NameHodo[ihodo]]);
  }

  tree->Branch("time0", &time0);
  tree->Branch("time0_seg", &time0_seg);
  tree->Branch("btof0", &btof0);
  tree->Branch("ftof0", &ftof0);

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
