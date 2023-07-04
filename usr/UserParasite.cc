// -*- C++ -*-

#include <cmath>
#include <iostream>
#include <sstream>

#include <UnpackerManager.hh>

#include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "DetectorID.hh"
#include "DCGeomMan.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "HodoHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoRawHit.hh"
#include "RawData.hh"
#include "RMAnalyzer.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"

#define FHitBranch 0 // make FiberHit branches (becomes heavy)

namespace
{
using namespace root;
const auto qnan = TMath::QuietNaN();
auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
auto& gRM = RMAnalyzer::GetInstance();
auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
  Int_t spill;

  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  Int_t t1nhits;
  Int_t t1hitpat[MaxHits];
  Double_t t1t[NumOfSegT1][MaxDepth];

  Int_t t2nhits;
  Int_t t2hitpat[MaxHits];
  Double_t t2t[NumOfSegT1][MaxDepth];

  Int_t bacnhits;
  Int_t bachitpat[MaxHits];
  Double_t baca[NumOfSegE72BAC];
  Double_t bact[NumOfSegE72BAC][MaxDepth];

  Int_t sacnhits;
  Int_t sachitpat[MaxHits];
  Double_t saca[NumOfSegE90SAC];
  Double_t sact[NumOfSegE90SAC][MaxDepth];

  Int_t kvcnhits;
  Int_t kvchitpat[MaxHits];
  Double_t kvcua[NumOfSegE72KVC];
  Double_t kvcut[NumOfSegE72KVC][MaxDepth];
  Double_t kvcda[NumOfSegE72KVC];
  Double_t kvcdt[NumOfSegE72KVC][MaxDepth];

  Int_t kvcsumnhits;
  Int_t kvcsumhitpat[MaxHits];
  Double_t kvcsuma[NumOfSegE72KVC];
  Double_t kvcsumt[NumOfSegE72KVC][MaxDepth];

  Int_t bh2nhits;
  Int_t bh2hitpat[MaxHits];
  Double_t bh2ua[NumOfSegE42BH2];
  Double_t bh2ut[NumOfSegE42BH2][MaxDepth];
  Double_t bh2da[NumOfSegE42BH2];
  Double_t bh2dt[NumOfSegE42BH2][MaxDepth];

  Int_t bh2mtnhits;
  Int_t bh2mthitpat[MaxHits];
  Double_t bh2mtt[NumOfSegE42BH2][MaxDepth];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum      = 0;
  spill      = 0;
  t1nhits   = 0;
  t2nhits   = 0;
  bacnhits   = 0;
  sacnhits   = 0;
  kvcnhits   = 0;
  kvcsumnhits   = 0;
  bh2nhits   = 0;
  bh2mtnhits   = 0;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<MaxHits; ++it){
    t1hitpat[it]   = -1;
    t2hitpat[it]   = -1;
    bachitpat[it]   = -1;
    sachitpat[it]   = -1;
    kvchitpat[it]   = -1;
    kvcsumhitpat[it]   = -1;
    bh2hitpat[it]   = -1;
    bh2mthitpat[it]   = -1;
  }

  for(Int_t it=0; it<NumOfSegT1; ++it){
    for(Int_t m = 0; m<MaxDepth; ++m){
      t1t[it][m] = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegT2; ++it){
    for(Int_t m = 0; m<MaxDepth; ++m){
      t2t[it][m] = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegE72BAC; it++){
    baca[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      bact[it][m] = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegE90SAC; it++){
    saca[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      sact[it][m] = qnan;
    }
  }
  for(Int_t it=0; it<NumOfSegE72KVC; ++it){
    kvcua[it] = qnan;
    kvcda[it] = qnan;
    kvcsuma[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      kvcut[it][m] = qnan;
      kvcdt[it][m] = qnan;
      kvcsumt[it][m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegE42BH2; ++it){
    bh2ua[it] = qnan;
    bh2da[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      bh2ut[it][m] = qnan;
      bh2dt[it][m] = qnan;
      bh2mtt[it][m] = qnan;
    }
  }

}


//_____________________________________________________________________________
namespace root
{
Event  event;
TH1*   h[MaxHist];
TTree* tree;
enum eDetHid {
  T1Hid     = 10000,
  T2Hid     = 20000,
  BACHid    = 30000,
  SACHid    = 40000,
  KVCHid    = 50000,
  BH2Hid    = 50000
};
}

//_____________________________________________________________________________
Bool_t
ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingNormal()
{
  static const auto MinTdcT1 = gUser.GetParameter("TdcT1", 0);
  static const auto MaxTdcT1 = gUser.GetParameter("TdcT1", 1);
  static const auto MinTdcT2 = gUser.GetParameter("TdcT2", 0);
  static const auto MaxTdcT2 = gUser.GetParameter("TdcT2", 1);
  static const auto MinTdcBAC = gUser.GetParameter("TdcE72BAC", 0);
  static const auto MaxTdcBAC = gUser.GetParameter("TdcE72BAC", 1);
  static const auto MinTdcSAC = gUser.GetParameter("TdcE90SAC", 0);
  static const auto MaxTdcSAC = gUser.GetParameter("TdcE90SAC", 1);
  static const auto MinTdcKVC = gUser.GetParameter("TdcE72KVC", 0);
  static const auto MaxTdcKVC = gUser.GetParameter("TdcE72KVC", 1);
  static const auto MinTdcBH2 = gUser.GetParameter("TdcE42BH2", 0);
  static const auto MaxTdcBH2 = gUser.GetParameter("TdcE42BH2", 1);

  RawData rawData;
  HodoAnalyzer hodoAna(rawData);

  rawData.DecodeHits("TFlag");
  rawData.DecodeHits("E72BAC");
  rawData.DecodeHits("E90SAC");
  rawData.DecodeHits("E72KVC");
  rawData.DecodeHits("E42BH2");
  rawData.DecodeHits("T1");
  rawData.DecodeHits("T2");

  gRM.Decode();

  // event.evnum = gRM.EventNumber();
  event.evnum = gUnpacker.get_event_number();
  event.spill = gRM.SpillNumber();

  HF1(1, 0);

  //**************************************************************************
  //****************** RawData

  // Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  for(const auto& hit: rawData.GetHodoRawHC("TFlag")){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc();
    if(tdc > 0){
      event.trigpat[trigger_flag.count()] = seg;
      event.trigflag[seg] = tdc;
      trigger_flag.set(seg);
    }
  }

  if(trigger_flag[trigger::kSpillOnEnd] || trigger_flag[trigger::kSpillOffEnd])
    return true;

  HF1(1, 1);

  ///// T1
  {
    Int_t t1_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("T1");
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.t1t[seg-1][m] = T;
        if(MinTdcT1 < T && T < MaxTdcT1) is_hit = true;
      }
      // Hitpat
      if(is_hit){
        event.t1hitpat[t1_nhits++] = seg;
      }
    }
    event.t1nhits = t1_nhits;
  }
  ///// T2
  {
    Int_t t2_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("T2");
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.t2t[seg-1][m] = T;
        if(MinTdcT2 < T && T < MaxTdcT2) is_hit = true;
      }
      // Hitpat
      if(is_hit){
        event.t2hitpat[t2_nhits++] = seg;
      }
    }
    event.t2nhits = t2_nhits;
  }
  ///// BAC
  {
    Int_t bac_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("E72BAC");
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t A = hit->GetAdcUp();
      event.baca[seg-1] = A;
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.bact[seg-1][m] = T;
        if(MinTdcBAC < T && T < MaxTdcBAC) is_hit = true;
      }
      // Hitpat
      if(is_hit){
        event.bachitpat[bac_nhits++] = seg;
      }
    }
    event.bacnhits = bac_nhits;
  }
  ///// SAC
  {
    Int_t sac_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("E90SAC");
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t A = hit->GetAdcUp();
      event.saca[seg-1] = A;
      Bool_t is_hit = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.sact[seg-1][m] = T;
        if(MinTdcSAC < T && T < MaxTdcSAC) is_hit = true;
      }
      // Hitpat
      if(is_hit){
        event.sachitpat[sac_nhits++] = seg;
      }
    }
    event.sacnhits = sac_nhits;
  }

  ///// KVC
  {
    Int_t kvc_nhits = 0;
    Int_t kvcsum_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("E72KVC");
    Int_t nh = cont.size();
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      const auto &hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      // Up
      Int_t Au = hit->GetAdcUp();
      event.kvcua[seg-1] = Au;
      Bool_t is_hit_u = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.kvcut[seg-1][m] = T;
        if(MinTdcKVC < T && T < MaxTdcKVC) is_hit_u = true;
      }
      // Down
      Int_t Ad = hit->GetAdcDown();
      event.kvcda[seg-1] = Ad;
      Bool_t is_hit_d = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcDown(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcDown(m);
        event.kvcdt[seg-1][m] = T;
        if(MinTdcKVC < T && T < MaxTdcKVC) is_hit_d = true;
      }
      // HitPat
      if(is_hit_u || is_hit_d){
        ++nh1;
      }
      if(is_hit_u && is_hit_d){
        event.kvchitpat[kvc_nhits++] = seg;
        ++nh2;
      }
      // SUM
      Int_t As = hit->GetAdcExtra();
      event.kvcsuma[seg-1] = As;
      Bool_t is_hit_s = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcExtra(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcExtra(m);
        event.kvcsumt[seg-1][m] = T;
        if(MinTdcKVC < T && T < MaxTdcKVC) is_hit_s = true;
      }
      if(is_hit_s){
        event.kvcsumhitpat[kvcsum_nhits++] = seg;
      }
    }
    event.kvcnhits = kvc_nhits;
    event.kvcsumnhits = kvcsum_nhits;
  }

  ///// BH2
  {
    Int_t bh2_nhits = 0;
    Int_t bh2mt_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("E42BH2");
    Int_t nh = cont.size();
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      // Up
      Int_t Au = hit->GetAdcUp();
      event.bh2ua[seg-1] = Au;
      Bool_t is_hit_u = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        event.bh2ut[seg-1][m] = T;
        if(MinTdcBH2 < T && T < MaxTdcBH2) is_hit_u = true;
      }
      // Down
      Int_t Ad = hit->GetAdcDown();
      event.bh2da[seg-1] = Ad;
      Bool_t is_hit_d = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcDown(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcDown(m);
        event.bh2dt[seg-1][m] = T;
        if(MinTdcBH2 < T && T < MaxTdcBH2) is_hit_d = true;
      }
      // HitPat
      if(is_hit_u || is_hit_d){
        ++nh1;
      }
      if(is_hit_u && is_hit_d){
        event.bh2hitpat[bh2_nhits++] = seg;
        ++nh2;
      }
      // Meantimer
      Bool_t is_hit_mt = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcExtra(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcExtra(m);
        event.bh2mtt[seg-1][m] = T;
        if(MinTdcBH2 < T && T < MaxTdcBH2) is_hit_mt = true;
      }
      // HitPat
      if(is_hit_mt){
        event.bh2mthitpat[bh2mt_nhits++] = seg;
      }
    }
    event.bh2nhits = bh2_nhits;
    event.bh2mtnhits = bh2mt_nhits;
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
namespace
{
const Int_t    NbinAdc = 4096;
const Double_t MinAdc  =    0.;
const Double_t MaxAdc  = 4096.;

const Int_t    NbinTdc = 4096;
const Double_t MinTdc  =    0.;
const Double_t MaxTdc  = 4096.;

const Int_t    NbinTdcHr = 1e6/10;
const Double_t MinTdcHr  =  0.;
const Double_t MaxTdcHr  = 1e6;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  HB1(1, "Status", 20, 0., 20.);
  HB1(10, "Trigger HitPat", NumOfSegTrig, 0., Double_t(NumOfSegTrig));
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1(10+i+1, Form("Trigger Trig %d", i+1), 0x1000, 0, 0x1000);
  }

  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("spill",     &event.spill,     "spill/I");
  //Trig
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  //T1
  tree->Branch("t1nhits",   &event.t1nhits,    "t1nhits/I");
  tree->Branch("t1hitpat",   event.t1hitpat,   Form("t1hitpat[%d]/I",NumOfSegT1));
  tree->Branch("t1t",        event.t1t,        Form("t1t[%d][%d]/D", NumOfSegT1, MaxDepth));
  //T2
  tree->Branch("t2nhits",   &event.t2nhits,    "t2nhits/I");
  tree->Branch("t2hitpat",   event.t2hitpat,   Form("t2hitpat[%d]/I",NumOfSegT2));
  tree->Branch("t2t",        event.t2t,        Form("t2t[%d][%d]/D", NumOfSegT2, MaxDepth));
  //BAC
  tree->Branch("bacnhits",   &event.bacnhits,   "bacnhits/I");
  tree->Branch("bachitpat",   event.bachitpat,  Form("bachitpat[%d]/I", NumOfSegE72BAC));
  tree->Branch("baca",        event.baca,       Form("baca[%d]/D", NumOfSegE72BAC));
  tree->Branch("bact",        event.bact,       Form("bact[%d][%d]/D", NumOfSegE72BAC, MaxDepth));
  //SAC
  tree->Branch("sacnhits",   &event.sacnhits,   "sacnhits/I");
  tree->Branch("sachitpat",   event.sachitpat,  Form("sachitpat[%d]/I", NumOfSegE90SAC));
  tree->Branch("saca",        event.saca,       Form("saca[%d]/D", NumOfSegE90SAC));
  tree->Branch("sact",        event.sact,       Form("sact[%d][%d]/D", NumOfSegE90SAC, MaxDepth));
  //KVC
  tree->Branch("kvcnhits",   &event.kvcnhits,    "kvcnhits/I");
  tree->Branch("kvchitpat",   event.kvchitpat,   Form("kvchitpat[%d]/I", NumOfSegE72KVC));
  tree->Branch("kvcua",       event.kvcua,       Form("kvcua[%d]/D", NumOfSegE72KVC));
  tree->Branch("kvcut",       event.kvcut,       Form("kvcut[%d][%d]/D", NumOfSegE72KVC, MaxDepth));
  tree->Branch("kvcda",       event.kvcda,       Form("kvcda[%d]/D", NumOfSegE72KVC));
  tree->Branch("kvcdt",       event.kvcdt,       Form("kvcdt[%d][%d]/D", NumOfSegE72KVC, MaxDepth));
  //KVCSUM
  tree->Branch("kvcsumnhits",   &event.kvcsumnhits,    "kvcsumnhits/I");
  tree->Branch("kvcsumhitpat",   event.kvcsumhitpat,   Form("kvcsumhitpat[%d]/I", NumOfSegE72KVC));
  tree->Branch("kvcsuma",        event.kvcsuma,       Form("kvcsuma[%d]/D", NumOfSegE72KVC));
  tree->Branch("kvcsumt",        event.kvcsumt,       Form("kvcsumt[%d][%d]/D", NumOfSegE72KVC, MaxDepth));
  //BH2
  tree->Branch("bh2nhits",   &event.bh2nhits,    "bh2nhits/I");
  tree->Branch("bh2hitpat",   event.bh2hitpat,   Form("bh2hitpat[%d]/I", NumOfSegE42BH2));
  tree->Branch("bh2ua",       event.bh2ua,       Form("bh2ua[%d]/D", NumOfSegE42BH2));
  tree->Branch("bh2ut",       event.bh2ut,       Form("bh2ut[%d][%d]/D", NumOfSegE42BH2, MaxDepth));
  tree->Branch("bh2da",       event.bh2da,       Form("bh2da[%d]/D", NumOfSegE42BH2));
  tree->Branch("bh2dt",       event.bh2dt,       Form("bh2dt[%d][%d]/D", NumOfSegE42BH2, MaxDepth));
  //BH2MT
  tree->Branch("bh2mtnhits",   &event.bh2mtnhits,    "bh2mtnhits/I");
  tree->Branch("bh2mthitpat",   event.bh2mthitpat,   Form("bh2mthitpat[%d]/I", NumOfSegE42BH2));
  tree->Branch("bh2mtt",        event.bh2mtt,       Form("bh2mtt[%d][%d]/D", NumOfSegE42BH2, MaxDepth));
  // HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO") &&
     InitializeParameter<HodoParamMan>("HDPRM") &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
