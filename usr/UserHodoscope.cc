// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iostream>
#include <sstream>

#include <UnpackerManager.hh>

// #include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "HodoHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "S2sLib.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "DCGeomMan.hh"

// #define TimeCut    1 // in cluster analysis
#define FHitBranch 0 // make FiberHit branches (becomes heavy)
#define HodoHitPos 0

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

  Int_t bh1nhits;
  Int_t bh1hitpat[MaxHits];
  Double_t bh1ua[NumOfSegBH1];
  Double_t bh1ut[NumOfSegBH1][MaxDepth];
  Double_t bh1da[NumOfSegBH1];
  Double_t bh1dt[NumOfSegBH1][MaxDepth];

  Int_t bh2nhits;
  Int_t bh2hitpat[MaxHits];
  Double_t bh2ua[NumOfSegBH2];
  Double_t bh2ut[NumOfSegBH2][MaxDepth];
  Double_t bh2da[NumOfSegBH2];
  Double_t bh2dt[NumOfSegBH2][MaxDepth];

  Int_t bacnhits;
  Int_t bachitpat[MaxHits];
  Double_t baca[NumOfSegBAC];
  Double_t bact[NumOfSegBAC][MaxDepth];

  Int_t bac1nhits;
  Int_t bac2nhits;


  Int_t tofnhits;
  Int_t tofhitpat[MaxHits];
  Double_t tofua[NumOfSegTOF];
  Double_t tofut[NumOfSegTOF][MaxDepth];
  Double_t tofda[NumOfSegTOF];
  Double_t tofdt[NumOfSegTOF][MaxDepth];

  Int_t ac1nhits;
  Int_t ac1hitpat[MaxHits];
  Double_t ac1a[NumOfSegAC1];
  Double_t ac1t[NumOfSegAC1][MaxDepth];

  Int_t wcnhits;
  Int_t wchitpat[MaxHits];
  Double_t wcua[NumOfSegWC];
  Double_t wcut[NumOfSegWC][MaxDepth];
  Double_t wcda[NumOfSegWC];
  Double_t wcdt[NumOfSegWC][MaxDepth];
  Int_t wcsumnhits;
  Int_t wcsumhitpat[MaxHits];
  Double_t wcsuma[NumOfSegWC];
  Double_t wcsumt[NumOfSegWC][MaxDepth];

  ////////// Normalized
  Double_t bh1mt[NumOfSegBH1][MaxDepth];
  Double_t bh1cmt[NumOfSegBH1][MaxDepth];
  Double_t bh1utime[NumOfSegBH1][MaxDepth];
  Double_t bh1uctime[NumOfSegBH1][MaxDepth];
  Double_t bh1dtime[NumOfSegBH1][MaxDepth];
  Double_t bh1dctime[NumOfSegBH1][MaxDepth];
  Double_t bh1de[NumOfSegBH1];
  Double_t bh1ude[NumOfSegBH1];
  Double_t bh1dde[NumOfSegBH1];

  Double_t bh2mt[NumOfSegBH2][MaxDepth];
  Double_t bh2cmt[NumOfSegBH2][MaxDepth];
  Double_t bh2utime[NumOfSegBH2][MaxDepth];
  Double_t bh2uctime[NumOfSegBH2][MaxDepth];
  Double_t bh2dtime[NumOfSegBH2][MaxDepth];
  Double_t bh2dctime[NumOfSegBH2][MaxDepth];
  Double_t bh2hitpos[NumOfSegBH2][MaxDepth];
  Double_t bh2de[NumOfSegBH2];
  Double_t bh2ude[NumOfSegBH2];
  Double_t bh2dde[NumOfSegBH2];

  Double_t bacmt[NumOfSegBAC][MaxDepth];
  Double_t bacde[NumOfSegBAC];

  Double_t t0[NumOfSegBH2][MaxDepth];
  Double_t ct0[NumOfSegBH2][MaxDepth];
  Double_t btof[NumOfSegBH1][NumOfSegBH2];
  Double_t cbtof[NumOfSegBH1][NumOfSegBH2];

  Double_t tofmt[NumOfSegTOF][MaxDepth]; 
  Double_t tofde[NumOfSegTOF];
  Double_t tofude[NumOfSegTOF];
  Double_t tofdde[NumOfSegTOF];
  Double_t tofctu[NumOfSegTOF][MaxDepth];
  Double_t tofctd[NumOfSegTOF][MaxDepth];
  Double_t tofcmt[NumOfSegTOF][MaxDepth];

  // Using WCSUM
  Double_t wcde[NumOfSegWC];
  Double_t wcmt[NumOfSegWC][MaxDepth];

  // Time0
  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;

  // Btof0 BH1
  Double_t Btof0Seg;
  Double_t deBtof0;
  Double_t Btof0;
  Double_t CBtof0;

  //SAC3
  //Int_t sac3nhits;
  Double_t sac3a[NumOfSegSAC3];
  Double_t sac3t[NumOfSegSAC3][MaxDepth];

  //SFV
  //Int_t sfvnhits;
  //Int_t sfvhitpat[MaxHits];
  Double_t sfvt[NumOfSegSFV][MaxDepth];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum      = 0;
  spill      = 0;
  bh1nhits   = 0;
  bacnhits   = 0;
  bac1nhits   = 0;
  bac2nhits   = 0;
  bh2nhits   = 0;
  tofnhits   = 0;
  ac1nhits   = 0;
  wcnhits    = 0;
  wcsumnhits = 0;
  Time0Seg = qnan;
  deTime0  = qnan;
  Time0    = qnan;
  CTime0   = qnan;
  Btof0Seg = qnan;
  deBtof0  = qnan;
  Btof0    = qnan;
  CBtof0   = qnan;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<MaxHits; ++it){
    bh1hitpat[it]   = -1;
    bh2hitpat[it]   = -1;
    bachitpat[it]   = -1;
    tofhitpat[it]   = -1;
    ac1hitpat[it]   = -1;
    wchitpat[it]    = -1;
    wcsumhitpat[it] = -1;
  }

  for(Int_t it=0; it<NumOfSegBH1; ++it){
    bh1ua[it] = qnan;
    bh1da[it] = qnan;
    bh1de[it] = qnan;
    bh1ude[it] = qnan;
    bh1dde[it] = qnan;
    for(Int_t that=0; that<NumOfSegBH2; ++that){
      btof[it][that]  = qnan;
      cbtof[it][that] = qnan;
    }
    for(Int_t m = 0; m<MaxDepth; ++m){
      bh1ut[it][m] = qnan;
      bh1dt[it][m] = qnan;
      bh1mt[it][m] = qnan;
      bh1utime[it][m] = qnan;
      bh1dtime[it][m] = qnan;
      bh1uctime[it][m] = qnan;
      bh1dctime[it][m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegBH2; ++it){
    bh2ua[it] = qnan;
    bh2da[it] = qnan;
    bh2de[it] = qnan;
    bh2dde[it] = qnan;
    bh2ude[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      bh2ut[it][m] = qnan;
      bh2dt[it][m] = qnan;
      bh2mt[it][m] = qnan;
      t0[it][m]    = qnan;
      ct0[it][m]   = qnan;
      bh2utime[it][m] = qnan;
      bh2dtime[it][m] = qnan;
      bh2uctime[it][m] = qnan;
      bh2dctime[it][m] = qnan;
      bh2hitpos[it][m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegBAC; it++){
    baca[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      bact[it][m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegTOF; it++){
    tofua[it] = qnan;
    tofda[it] = qnan;
    tofde[it] = qnan;
    tofude[it] = qnan;
    tofdde[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      tofut[it][m] = qnan;
      tofdt[it][m] = qnan;
      tofmt[it][m] = qnan;
      tofctu[it][m] = qnan;
      tofctd[it][m] = qnan;
      tofcmt[it][m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegAC1; it++){
    ac1a[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      ac1t[it][m]  = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegWC; ++it){
    wcua[it]   = qnan;
    wcda[it]   = qnan;
    wcsuma[it] = qnan;
    wcde[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      wcut[it][m]   = qnan;
      wcdt[it][m]   = qnan;
      wcsumt[it][m] = qnan;
      wcmt[it][m]   = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegSAC3; it++){
    sac3a[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      sac3t[it][m]  = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegSFV; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      sfvt[it][m]  = qnan;
    }
  }

}

//_____________________________________________________________________________
struct Dst
{
  Int_t evnum;
  Int_t spill;

  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  Int_t    nhBh1;
  Int_t    csBh1[NumOfSegBH1*MaxDepth];
  Double_t Bh1Seg[NumOfSegBH1*MaxDepth];
  Double_t tBh1[NumOfSegBH1*MaxDepth];
  Double_t dtBh1[NumOfSegBH1*MaxDepth];
  Double_t deBh1[NumOfSegBH1*MaxDepth];
  Double_t btof[NumOfSegBH1*MaxDepth];
  Double_t cbtof[NumOfSegBH1*MaxDepth];

  // Btof0 BH1
  Double_t Btof0Seg;
  Double_t deBtof0;
  Double_t Btof0;
  Double_t CBtof0;

  Int_t    nhBh2;
  Int_t    csBh2[NumOfSegBH2*MaxDepth];
  Double_t Bh2Seg[NumOfSegBH2*MaxDepth];
  Double_t tBh2[NumOfSegBH2*MaxDepth];
  Double_t t0Bh2[NumOfSegBH2*MaxDepth];
  Double_t dtBh2[NumOfSegBH2*MaxDepth];
  Double_t deBh2[NumOfSegBH2*MaxDepth];
  Double_t posBh2[NumOfSegBH2*MaxDepth];

  Int_t    nhBac;
  Int_t    csBac[NumOfSegBAC*MaxDepth];
  Double_t BacSeg[NumOfSegBAC*MaxDepth];
  Double_t tBac[NumOfSegBAC*MaxDepth];
  Double_t deBac[NumOfSegBAC*MaxDepth];

  // Time0
  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;

  Int_t    nhTof;
  Int_t    csTof[NumOfSegTOF*MaxDepth];
  Double_t TofSeg[NumOfSegTOF*MaxDepth];
  Double_t tTof[NumOfSegTOF*MaxDepth];
  Double_t dtTof[NumOfSegTOF*MaxDepth];
  Double_t deTof[NumOfSegTOF*MaxDepth];

  Int_t    nhAc1;
  Int_t    csAc1[NumOfSegAC1*MaxDepth];
  Double_t Ac1Seg[NumOfSegAC1*MaxDepth];
  Double_t tAc1[NumOfSegAC1*MaxDepth];

  Int_t    nhWc;
  Int_t    csWc[NumOfSegWC*MaxDepth];
  Double_t WcSeg[NumOfSegWC*MaxDepth];
  Double_t tWc[NumOfSegWC*MaxDepth];
  Double_t deWc[NumOfSegWC*MaxDepth];

  // for HodoParam
  Double_t tofua[NumOfSegTOF];
  Double_t tofut[NumOfSegTOF][MaxDepth];
  Double_t tofda[NumOfSegTOF];
  Double_t tofdt[NumOfSegTOF][MaxDepth];
  Double_t utTofSeg[NumOfSegTOF][MaxDepth];
  Double_t dtTofSeg[NumOfSegTOF][MaxDepth];
  Double_t udeTofSeg[NumOfSegTOF];
  Double_t ddeTofSeg[NumOfSegTOF];

  void clear();
};

//_____________________________________________________________________________
void
Dst::clear()
{
  nhBh1    = 0;
  nhBac    = 0;
  nhBh2    = 0;
  nhTof    = 0;
  nhAc1    = 0;
  evnum    = 0;
  spill    = 0;
  Time0Seg = qnan;
  deTime0  = qnan;
  Time0    = qnan;
  CTime0   = qnan;
  Btof0Seg = qnan;
  deBtof0  = qnan;
  Btof0    = qnan;
  CBtof0   = qnan;

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<NumOfSegBH1; ++it){
    for(Int_t m=0; m<MaxDepth; ++m){
      csBh1[MaxDepth*it + m]  = 0;
      Bh1Seg[MaxDepth*it + m] = qnan;
      tBh1[MaxDepth*it + m]   = qnan;
      dtBh1[MaxDepth*it + m]  = qnan;
      deBh1[MaxDepth*it + m]  = qnan;
      btof[MaxDepth*it + m]   = qnan;
      cbtof[MaxDepth*it + m]  = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegBH2; ++it){
    for(Int_t m=0; m<MaxDepth; ++m){
      csBh2[MaxDepth*it + m]  = 0;
      Bh2Seg[MaxDepth*it + m] = qnan;
      tBh2[MaxDepth*it + m]   = qnan;
      t0Bh2[MaxDepth*it + m]  = qnan;
      dtBh2[MaxDepth*it + m]  = qnan;
      deBh2[MaxDepth*it + m]  = qnan;
      posBh2[MaxDepth*it + m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegBAC; it++){
    BacSeg[it] = qnan;
    deBac[it]  = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      tBac[MaxDepth*it+m] = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegTOF; it++){
    tofua[it] = qnan;
    tofda[it] = qnan;
    udeTofSeg[it] = qnan;
    ddeTofSeg[it] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      tofut[it][m] = qnan;
      tofdt[it][m] = qnan;
      utTofSeg[it][m]  = qnan;
      dtTofSeg[it][m]  = qnan;
      csTof[MaxDepth*it + m]  = 0;
      TofSeg[MaxDepth*it + m] = qnan;
      tTof[MaxDepth*it + m]   = qnan;
      dtTof[MaxDepth*it + m]  = qnan;
      deTof[MaxDepth*it + m]  = qnan;
    }
  }

  for(Int_t it=0; it<NumOfSegAC1; it++){
    for(Int_t m=0; m<MaxDepth; ++m){
      csAc1[MaxDepth*it + m]  = 0;
      Ac1Seg[MaxDepth*it + m] = qnan;
      tAc1[MaxDepth*it + m]   = qnan;
    }
  }

}


//_____________________________________________________________________________
namespace root
{
Event  event;
Dst    dst;
TH1*   h[MaxHist];
TTree* tree;
TTree* hodo;
enum eDetHid {
  BH1Hid    = 10000,
  BH2Hid    = 20000,
  BACHid    = 30000,
  TOFHid    = 60000,
  AC1Hid    = 70000,
  WCHid     = 80000,
  WCSUMHid  = 90000
};
}

//_____________________________________________________________________________
Bool_t
ProcessingBegin()
{
  event.clear();
  dst.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingNormal()
{
  static const auto MinTdcBH1 = gUser.GetParameter("TdcBH1", 0);
  static const auto MaxTdcBH1 = gUser.GetParameter("TdcBH1", 1);
  static const auto MinTdcBH2 = gUser.GetParameter("TdcBH2", 0);
  static const auto MaxTdcBH2 = gUser.GetParameter("TdcBH2", 1);
  static const auto MinTdcBAC = gUser.GetParameter("TdcBAC", 0);
  static const auto MaxTdcBAC = gUser.GetParameter("TdcBAC", 1);
  static const auto MinTdcTOF = gUser.GetParameter("TdcTOF", 0);
  static const auto MaxTdcTOF = gUser.GetParameter("TdcTOF", 1);
  static const auto MinTdcAC1 = gUser.GetParameter("TdcAC1", 0);
  static const auto MaxTdcAC1 = gUser.GetParameter("TdcAC1", 1);
  static const auto MinTdcWC = gUser.GetParameter("TdcWC", 0);
  static const auto MaxTdcWC = gUser.GetParameter("TdcWC", 1);
#if HodoHitPos
  static const auto PropVelBH2 = gUser.GetParameter("PropagationBH2");
#endif

  RawData rawData;
  HodoAnalyzer hodoAna(rawData);


  gRM.Decode();

  // event.evnum = gRM.EventNumber();
  event.evnum = gUnpacker.get_event_number();
  event.spill = gRM.SpillNumber();
  dst.evnum   = gRM.EventNumber();
  dst.spill   = gRM.SpillNumber();

  HF1(1, 0);
  //**************************************************************************
  //****************** RawData

  // Trigger Flag
  rawData.DecodeHits("TFlag");
  std::bitset<NumOfSegTrig> trigger_flag;
  for(const auto& hit: rawData.GetHodoRawHC("TFlag")){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc();
    if(tdc > 0){
      event.trigpat[trigger_flag.count()] = seg;
      event.trigflag[seg] = tdc;
      dst.trigpat[trigger_flag.count()] = seg;
      dst.trigflag[seg] = tdc;
      trigger_flag.set(seg);
      HF1(10, seg);
      HF1(10+seg, tdc);
    }
  }

  if(trigger_flag[trigger::kSpillOnEnd] || trigger_flag[trigger::kSpillOffEnd])
    return true;

  HF1(1, 1);

  ///// BH1
  rawData.DecodeHits("BH1");
  {
    Int_t bh1_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("BH1");
    Int_t nh = cont.size();
    HF1(BH1Hid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(BH1Hid +1, seg-0.5);
      // Up
      Int_t Au = hit->GetAdcUp();
      HF1(BH1Hid +100*seg +1, Au);
      event.bh1ua[seg-1] = Au;
      Bool_t is_hit_u = false;
      Int_t m_u = 0;
      for(const auto& T: hit->GetArrayTdcUp()){
        HF1(BH1Hid +100*seg +3, T);
        if(m_u < MaxDepth) event.bh1ut[seg-1][m_u++] = T;
        if(MinTdcBH1 < T && T < MaxTdcBH1) is_hit_u = true;
      }
      if(is_hit_u) HF1(BH1Hid +100*seg +5, Au);
      else         HF1(BH1Hid +100*seg +7, Au);
      // Down
      Int_t Ad = hit->GetAdcDown();
      HF1(BH1Hid +100*seg +2, Ad);
      event.bh1da[seg-1] = Ad;
      Bool_t is_hit_d = false;
      Int_t m_d = 0;
      for(const auto& T: hit->GetArrayTdcDown()){
        HF1(BH1Hid +100*seg +4, T);
        if(m_d < MaxDepth) event.bh1dt[seg-1][m_d++] = T;
        if(MinTdcBH1 < T && T < MaxTdcBH1) is_hit_d = true;
      }
      if(is_hit_d) HF1(BH1Hid +100*seg +6, Ad);
      else         HF1(BH1Hid +100*seg +8, Ad);
      // HitPat
      if(is_hit_u || is_hit_d){
        ++nh1; HF1(BH1Hid +3, seg-0.5);
      }
      if(is_hit_u && is_hit_d){
        event.bh1hitpat[bh1_nhits++] = seg;
        ++nh2; HF1(BH1Hid +5, seg-0.5);
      }
    }
    HF1(BH1Hid +2, nh1); HF1(BH1Hid +4, nh2);
    event.bh1nhits = bh1_nhits;
  }

  ///// BH2
  rawData.DecodeHits("BH2");
  {
    Int_t bh2_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("BH2");
    Int_t nh = cont.size();
    HF1(BH2Hid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(BH2Hid +1, seg-0.5);
      // Up
      Int_t Au = hit->GetAdcUp();
      HF1(BH2Hid +100*seg +1, Au);
      event.bh2ua[seg-1] = Au;
      Bool_t is_hit_u = false;
      Int_t m_u = 0;
      for(const auto& T: hit->GetArrayTdcUp()){
        HF1(BH2Hid +100*seg +3, T);
        if(m_u < MaxDepth) event.bh2ut[seg-1][m_u++] = T;
        if(MinTdcBH2 < T && T < MaxTdcBH2) is_hit_u = true;
      }
      if(is_hit_u) HF1(BH2Hid +100*seg +5, Au);
      else         HF1(BH2Hid +100*seg +7, Au);
      // Down
      Int_t Ad = hit->GetAdcDown();
      HF1(BH2Hid +100*seg +2, Double_t(Ad));
      event.bh2da[seg-1] = Ad;
      Bool_t is_hit_d = false;
      Int_t m_d = 0;
      for(const auto& T: hit->GetArrayTdcDown()){
        HF1(BH2Hid +100*seg +4, Double_t(T));
        if(m_d < MaxDepth) event.bh2dt[seg-1][m_d++] = T;
        if(MinTdcBH2 < T && T < MaxTdcBH2) is_hit_d = true;
      }
      if(is_hit_d) HF1(BH2Hid +100*seg +6, Ad);
      else         HF1(BH2Hid +100*seg +8, Ad);
      // HitPat
      if(is_hit_u || is_hit_d){
        ++nh1; HF1(BH2Hid +3, seg-0.5);
      }
      if(is_hit_u && is_hit_d){
        event.bh2hitpat[bh2_nhits++] = seg;
        ++nh2; HF1(BH2Hid +5, seg-0.5);
      }
    }
    HF1(BH2Hid +2, nh1); HF1(BH2Hid +4, nh2);
    event.bh2nhits = bh2_nhits;
  }

  ///// BAC
  rawData.DecodeHits("BAC");
  {
    Int_t bac_nhits = 0;
    Int_t bac1_nhits = 0;
    Int_t bac2_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("BAC");
    Int_t nh = cont.size();
    HF1(BACHid, nh);
    Int_t nh1 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(BACHid+1, seg-0.5);
      Int_t A = hit->GetAdcUp();
      HF1(BACHid+100*seg+1, A);
      event.baca[seg-1] = A;
      Bool_t is_hit = false;
      Int_t m = 0;
      for(const auto& T: hit->GetArrayTdcLeading()){
        HF1(BACHid+100*seg+3, T);
        if(m < MaxDepth) event.bact[seg-1][m++] = T;
        if(MinTdcBAC < T && T < MaxTdcBAC){
	  is_hit = true;
	  if(seg==1)++bac1_nhits;
	  if(seg==2)++bac2_nhits;
	}
      }
      if(is_hit) HF1(BACHid+100*seg+5, A);
      else       HF1(BACHid+100*seg+7, A);
      // Hitpat
      if(is_hit){
        event.bachitpat[bac_nhits++] = seg;
        ++nh1; HF1(BACHid+3, seg-0.5);
      }
    }
    HF1(BACHid+2, nh1);
    event.bacnhits = bac_nhits;
    event.bac1nhits = bac1_nhits;
    event.bac2nhits = bac2_nhits;
  }

  ///// TOF
  rawData.DecodeHits("TOF");
  {
    Int_t tof_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("TOF");
    Int_t nh = cont.size();
    HF1(TOFHid, Double_t(nh));
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(TOFHid+1, seg-0.5);
      // Up
      Int_t Au = hit->GetAdcUp();
      HF1(TOFHid+100*seg+1, Au);
      event.tofua[seg-1] = Au;
      dst.tofua[seg-1] = Au;
      Bool_t is_hit_u = false;
      Int_t m_u = 0;
      for(const auto& T: hit->GetArrayTdcUp()){
        HF1(TOFHid +100*seg +3, T);
        if(m_u < MaxDepth){
          event.tofut[seg-1][m_u] = T;
          dst.tofut[seg-1][m_u] = T;
          ++m_u;
        }
        if(MinTdcTOF < T && T < MaxTdcTOF) is_hit_u = true;
      }
      if(is_hit_u) HF1(TOFHid+100*seg+5, Au);
      else         HF1(TOFHid+100*seg+7, Au);
      // Down
      Int_t Ad = hit->GetAdcDown();
      HF1(TOFHid+100*seg+2, Ad);
      event.tofda[seg-1] = Ad;
      dst.tofda[seg-1] = Ad;
      Bool_t is_hit_d = false;
      Int_t m_d = 0;
      for(const auto& T: hit->GetArrayTdcDown()){
        HF1(TOFHid +100*seg +4, Double_t(T));
        if(m_d < MaxDepth){
          event.tofdt[seg-1][m_d] = T;
          dst.tofdt[seg-1][m_d] = T;
          ++m_d;
        }
        if(MinTdcTOF < T && T < MaxTdcTOF)  is_hit_d = true;
      }
      if(is_hit_d) HF1(TOFHid+100*seg+6, Ad);
      else         HF1(TOFHid+100*seg+8, Ad);
      // Hitpat
      if(is_hit_u || is_hit_d){
        ++nh1; HF1(TOFHid+3, seg-0.5);
      }
      if(is_hit_u && is_hit_d){
        event.tofhitpat[tof_nhits++] = seg;
        ++nh2; HF1(TOFHid+5, seg-0.5);
      }
    }
    HF1(TOFHid+2, nh1); HF1(TOFHid+4, nh2);
    event.tofnhits = tof_nhits;
  }

  ///// AC1
  rawData.DecodeHits("AC1");
  {
    Int_t ac1_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("AC1");
    Int_t nh = cont.size();
    HF1(AC1Hid, nh);
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(AC1Hid +1, seg-0.5);
      Int_t A = hit->GetAdcUp();
      HF1(AC1Hid +100*seg +1, A);
      event.ac1a[seg-1] = A;
      Bool_t is_hit = false;
      Int_t m = 0;
      for(const auto& T: hit->GetArrayTdcLeading()){
        HF1(AC1Hid +100*seg +3, Double_t(T));
        if(m < MaxDepth) event.ac1t[seg-1][m++] = T;
        if(MinTdcAC1 < T && T < MaxTdcAC1) is_hit = true;
      }
      //HitPat
      if(is_hit){
        HF1(AC1Hid +3, seg-0.5);
        HF1(AC1Hid +5, seg-0.5);
        event.ac1hitpat[ac1_nhits++] = seg;
      }
    }
    HF1(AC1Hid +2, ac1_nhits);
    HF1(AC1Hid +4, ac1_nhits);
    event.ac1nhits = ac1_nhits;
  }

  ///// WC
  rawData.DecodeHits("WC");
  {
    Int_t wc_nhits = 0;
    Int_t wcsum_nhits = 0;
    const auto& cont = rawData.GetHodoRawHC("WC");
    Int_t nh = cont.size();
    HF1(WCHid, Double_t(nh));
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      HF1(WCHid +1, seg-0.5);
      // Up
      Int_t Au = hit->GetAdcUp();
      HF1(WCHid +100*seg +1, Au);
      event.wcua[seg-1] = Au;
      Bool_t is_hit_u = false;
      Int_t m_u = 0;
      for(const auto& T: hit->GetArrayTdcUp()){
        HF1(WCHid +100*seg +3, T);
        if(m_u < MaxDepth) event.wcut[seg-1][m_u++] = T;
        if(MinTdcWC < T && T < MaxTdcWC) is_hit_u = true;
      }
      if(is_hit_u) HF1(WCHid +100*seg +5, Au);
      else         HF1(WCHid +100*seg +7, Au);
      // Down
      Int_t Ad = hit->GetAdcDown();
      HF1(WCHid +100*seg +2, Ad);
      event.wcda[seg-1] = Ad;
      Bool_t is_hit_d = false;
      Int_t m_d = 0;
      for(const auto& T: hit->GetArrayTdcDown()){
        HF1(WCHid +100*seg +4, Double_t(T));
        if(m_d < MaxDepth) event.wcdt[seg-1][m_d++] = T;
        if(MinTdcWC < T && T < MaxTdcWC) is_hit_d = true;
      }
      if(is_hit_d) HF1(WCHid +100*seg +6, Ad);
      else         HF1(WCHid +100*seg +8, Ad);
      // HitPat
      if(is_hit_u || is_hit_d){
        ++nh1; HF1(WCHid +3, seg-0.5);
      }
      if(is_hit_u && is_hit_d){
        event.wchitpat[wc_nhits] = seg;
        wc_nhits++;
        ++nh2; HF1(WCHid +5, seg-0.5);
      }
      // Sum
      Int_t As = hit->GetAdcExtra();
      HF1(WCSUMHid +100*seg +1, As);
      event.wcsuma[seg-1] = As;
      Bool_t is_hit_s = false;
      Int_t m_s = 0;
      for(const auto& T: hit->GetArrayTdcExtra()){
        HF1(WCSUMHid +100*seg +3, T);
        if(m_s < MaxDepth) event.wcsumt[seg-1][m_s++] = T;
        if(MinTdcWC < T && T < MaxTdcWC) is_hit_s = true;
      }
      if(is_hit_s) HF1(WCSUMHid +100*seg +5, As);
      else         HF1(WCSUMHid +100*seg +7, As);
      // HitPat
      if(is_hit_s){
        event.wcsumhitpat[wcsum_nhits++] = seg;
        HF1(WCSUMHid +1, seg-0.5);
      }
    }
    HF1(WCHid +2, Double_t(nh1)); HF1(WCHid +4, Double_t(nh2));
    event.wcnhits = wc_nhits;
    HF1(WCSUMHid, wcsum_nhits);
    event.wcsumnhits = wcsum_nhits;
  }

   // SAC3
  {
    static const int device_id = gUnpacker.get_device_id("SAC3");
    for(Int_t seg = 0; seg<NumOfSegSAC3; seg++) {
      int plane = 0;
      int ch = 0;
      int data = 0; // adc
      int nhita = gUnpacker.get_entries( device_id, plane, seg, ch, data ); // adc
      if(nhita>0){
	int adc = gUnpacker.get( device_id, plane, seg, ch, data );
	event.sac3a[seg] = adc;
      }
      data = 1; // tdc leading
      int nhitt = gUnpacker.get_entries( device_id, plane, seg, ch, data ); // tdc
      for(int m=0; m<nhitt; m++){
	int tdc = gUnpacker.get( device_id, plane, seg, ch, data, m ); // tdc multihit
	event.sac3t[seg][m] = tdc;
      }
    }
  }

   // SFV
  {
    static const int device_id = gUnpacker.get_device_id("SFV");
    for(Int_t seg = 0; seg<NumOfSegSFV; seg++) {
      int plane = 0;
      int ch = 0;
      int data = 1; // tdc leading
      int nhitt = gUnpacker.get_entries( device_id, plane, seg, ch, data ); // tdc
      for(int m=0; m<nhitt; m++){
	int tdc = gUnpacker.get( device_id, plane, seg, ch, data, m ); // tdc multihit
	event.sfvt[seg][m] = tdc;
      }
    }
  }

  //**************************************************************************
  //****************** NormalizedData

  //BH1
  hodoAna.DecodeHits("BH1");
  // hodoAna.TimeCutBH1(-10, 10);
  // hodoAna.TimeCut("BH1", -2, 2);
  {
    Int_t nh = hodoAna.GetNHits("BH1");
    HF1(BH1Hid+10, Double_t(nh));
    Int_t nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit("BH1", i);
      if(!hit) continue;
      Int_t seg = hit->SegmentId()+1;

      Int_t n_mhit = hit->GetEntries();
      for(Int_t m=0; m<n_mhit; ++m){
        HF1(BH1Hid+11, seg-0.5);

        Double_t au  = hit->GetAUp(),    ad = hit->GetADown();
        Double_t tu  = hit->GetTUp(m),   td = hit->GetTDown(m);
        Double_t ctu = hit->GetCTUp(m), ctd = hit->GetCTDown(m);
        Double_t mt  = hit->MeanTime(m),cmt = hit->CMeanTime(m);
        Double_t de  = hit->DeltaE();
        Double_t ude  = hit->UDeltaE();
        Double_t dde  = hit->DDeltaE();
        event.bh1mt[seg-1][m] = mt;
        event.bh1utime[seg-1][m] = tu;
        event.bh1dtime[seg-1][m] = td;
        event.bh1uctime[seg-1][m] = ctu;
        event.bh1dctime[seg-1][m] = ctd;
        event.bh1de[seg-1]    = de;
        event.bh1ude[seg-1]    = ude;
        event.bh1dde[seg-1]    = dde;
        HF1(BH1Hid+100*seg+11, tu);      HF1(BH1Hid+100*seg+12, td);
        HF1(BH1Hid+100*seg+13, mt);
        HF1(BH1Hid+100*seg+17, ctu);     HF1(BH1Hid+100*seg+18, ctd);
        HF1(BH1Hid+100*seg+19, cmt);     HF1(BH1Hid+100*seg+20, ctu-ctd);
        HF2(BH1Hid+100*seg+21, tu, au);  HF2(BH1Hid+100*seg+22, td, ad);
        HF2(BH1Hid+100*seg+23, ctu, au); HF2(BH1Hid+100*seg+24, ctd, ad);
        HF2(BH1Hid+100*seg+31, au, tu);  HF2(BH1Hid+100*seg+32, ad, td);
        HF2(BH1Hid+100*seg+33, au, ctu);  HF2(BH1Hid+100*seg+34, ad, ctd);
        HF1(BH1Hid+12, cmt);

        if(m == 0){
          HF1(BH1Hid+100*seg+14, au);    HF1(BH1Hid+100*seg+15, ad);
          HF1(BH1Hid+100*seg+16, de);    HF1(BH1Hid+13, de);
        }

        if(de>0.5){
          ++nh2; HF1(BH1Hid+15, seg-0.5);
          HF1(BH1Hid+16, cmt);
        }
      }
    }

    HF1(BH1Hid+14, Double_t(nh2));
    for(Int_t i1=0; i1<nh; ++i1){
      const auto& hit1 = hodoAna.GetHit("BH1", i1);
      if(!hit1 || hit1->DeltaE()<=0.5) continue;
      Int_t seg1 = hit1->SegmentId()+1;
      for(Int_t i2=0; i2<nh; ++i2){
        if(i1==i2) continue;
        const auto& hit2 = hodoAna.GetHit("BH1", i2);
        if(!hit2 || hit2->DeltaE()<=0.5) continue;
        Int_t seg2 = hit2->SegmentId()+1;

        if(1 == hit1->GetEntries() && 1 == hit2->GetEntries()){
          Double_t ct1 = hit1->CMeanTime(), ct2 = hit2->CMeanTime();
          HF2(BH1Hid+21, seg1-0.5, seg2-0.5);
          HF2(BH1Hid+22, ct1, ct2);
          HF1(BH1Hid+23, ct2-ct1);
          if(std::abs(ct2-ct1)<2.0){
            HF2(BH1Hid+24, seg1-0.5, seg2-0.5);
          }
        }
      }//for(i2)
    }//for(i1)
  }

  // BH2
  hodoAna.DecodeHits<BH2Hit>("BH2");
  // hodoAna.TimeCut("BH2", -2, 2);
  {
    Int_t nh = hodoAna.GetNHits("BH2");
    HF1(BH2Hid+10, Double_t(nh));
    Int_t nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit<BH2Hit>("BH2", i);
      if(!hit) continue;
      Int_t seg = hit->SegmentId()+1;

      Int_t n_mhit = hit->GetEntries();
      for(Int_t m=0; m<n_mhit; ++m){
        HF1(BH2Hid+11, seg-0.5);
        Double_t au  = hit->GetAUp(),   ad  = hit->GetADown();
        Double_t tu  = hit->GetTUp(m),  td  = hit->GetTDown(m);
        Double_t ctu = hit->GetCTUp(m), ctd = hit->GetCTDown(m);
        Double_t mt  = hit->MeanTime(m),cmt = hit->CMeanTime(m);
        Double_t de  = hit->DeltaE();
        Double_t ude  = hit->UDeltaE();
        Double_t dde  = hit->DDeltaE();
        Double_t ut0  = hit->UTime0(m), dt0  = hit->DTime0(m);
        Double_t uct0 = hit->UCTime0(m),dct0 = hit->DCTime0(m);
        Double_t t0   = hit->Time0(m),  ct0  = hit->CTime0(m);
#if HodoHitPos
	Double_t ctdiff = hit->TimeDiff(m);
	event.bh2hitpos[seg-1][m] = 0.5*PropVelBH2*ctdiff;
#endif
        event.bh2mt[seg-1][m] = mt;
        event.bh2utime[seg-1][m] = tu;
        event.bh2dtime[seg-1][m] = td;
        event.bh2uctime[seg-1][m] = ctu;
        event.bh2dctime[seg-1][m] = ctd;
        event.bh2de[seg-1]    = de;
        event.bh2ude[seg-1]    = ude;
        event.bh2dde[seg-1]    = dde;
        event.t0[seg-1][m]    = t0;
        event.ct0[seg-1][m]   = ct0;

	HF1(BH2Hid+100*seg+11, tu);      HF1(BH2Hid+100*seg+12, td);
        HF1(BH2Hid+100*seg+13, mt);
        HF1(BH2Hid+100*seg+17, ctu);     HF1(BH2Hid+100*seg+18, ctd);
        HF1(BH2Hid+100*seg+19, cmt);     HF1(BH2Hid+100*seg+20, ctu-ctd);
        HF1(BH2Hid+100*seg+21, ut0);     HF1(BH2Hid+100*seg+22, dt0);
        HF1(BH2Hid+100*seg+23, uct0);    HF1(BH2Hid+100*seg+24, dct0);
        HF1(BH2Hid+100*seg+25, t0);      HF1(BH2Hid+100*seg+26, ct0);
        HF2(BH2Hid+100*seg+27, tu, au);  HF2(BH2Hid+100*seg+28, td, ad);
        HF2(BH2Hid+100*seg+29, ctu, au); HF2(BH2Hid+100*seg+30, ctd, ad);

        // HF1(BH2Hid+100*seg+11, tu); HF1(BH2Hid+100*seg+13, mt);
        // HF1(BH2Hid+100*seg+14, au); HF1(BH2Hid+100*seg+16, de);
        // HF1(BH2Hid+100*seg+17, ctu); HF1(BH2Hid+100*seg+19, cmt);
        // HF2(BH2Hid+100*seg+21, tu, au); HF2(BH2Hid+100*seg+23, ctu, au);
        // HF1(BH2Hid+12, cmt);

        if(m == 0){
          HF1(BH2Hid+100*seg+14, au);	   HF1(BH2Hid+100*seg+15, ad);
          HF1(BH2Hid+100*seg+16, de);    HF1(BH2Hid+13, de);
        }

        if(de>0.5){
          ++nh2; HF1(BH2Hid+15, seg-0.5); HF1(BH2Hid+16, cmt);
        }
      }
    }//for(i)
    HF1(BH2Hid+14, Double_t(nh2));
    for(Int_t i1=0; i1<nh; ++i1){
      const auto& hit1 = hodoAna.GetHit<BH2Hit>("BH2", i1);
      if(!hit1 || hit1->DeltaE()<=0.5) continue;
      Int_t seg1 = hit1->SegmentId()+1;
      for(Int_t i2=0; i2<nh; ++i2){
        if(i1==i2) continue;
        const auto& hit2 = hodoAna.GetHit<BH2Hit>("BH2", i2);
        if(!hit2 || hit2->DeltaE()<=0.5) continue;
        Int_t seg2 = hit2->SegmentId()+1;

        if(1 == hit1->GetEntries() && 1 == hit2->GetEntries()){
          Double_t ct1 = hit1->CMeanTime(), ct2 = hit2->CMeanTime();
          HF2(BH2Hid+21, seg1-0.5, seg2-0.5);
          HF2(BH2Hid+22, ct1, ct2);
          HF1(BH2Hid+23, ct2-ct1);
          if(std::abs(ct2-ct1)<2.0){
            HF2(BH2Hid+24, seg1-0.5, seg2-0.5);
          }
        }
      }//for(i2)
    }//for(i1)

    Int_t nc=hodoAna.GetNClusters("BH2");
    HF1(BH2Hid+30, Double_t(nc));
    Int_t nc2=0;

    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster("BH2", i);
      if(!cl) continue;
      Int_t cs=cl->ClusterSize();
      Double_t ms = cl->MeanSeg()+1;
      Double_t cmt= cl->CMeanTime();
      Double_t de = cl->DeltaE();
      // Double_t mt = cl->MeanTime();
      HF1(BH2Hid+31, Double_t(cs));
      HF1(BH2Hid+32, ms-0.5);
      HF1(BH2Hid+33, cmt); HF1(BH2Hid+34, de);
      if(de>0.5){
        ++nc2; HF1(BH2Hid+36, cmt);
      }

      for(Int_t i2=0; i2<nc; ++i2){
        if(i2==i) continue;
        const auto& cl2 = hodoAna.GetCluster("BH2", i2);
        if(!cl2) continue;
        Double_t ms2=cl2->MeanSeg()+1, cmt2=cl2->CMeanTime(),
          de2=cl2->DeltaE();
        if(de<=0.5 || de2<=0.5) continue;
        HF2(BH2Hid+41, ms-0.5, ms2-0.5);
        HF2(BH2Hid+42, cmt, cmt2);
        HF1(BH2Hid+43, cmt2-cmt);
        if(std::abs(cmt2-cmt)<2.0){
          HF2(BH2Hid+44, ms-0.5, ms2-0.5);
        }
      }//for(i2)
    }//for(i)
    HF1(BH2Hid+35, Double_t(nc2));

    const auto& cl_time0 = hodoAna.GetTime0BH2Cluster();
    if(cl_time0){
      event.Time0Seg = cl_time0->MeanSeg()+1;
      event.deTime0 = cl_time0->DeltaE();
      event.Time0 = cl_time0->Time0();
      event.CTime0 = cl_time0->CTime0();
      dst.Time0Seg = cl_time0->MeanSeg()+1;
      dst.deTime0 = cl_time0->DeltaE();
      dst.Time0 = cl_time0->Time0();
      dst.CTime0 = cl_time0->CTime0();
      HF1(100, cl_time0->MeanTime());
      HF1(101, cl_time0->CTime0());
    }

    // BTOF0 segment
    const auto& cl_btof0 = hodoAna.GetBtof0BH1Cluster();
    if(cl_btof0){
      event.Btof0Seg = cl_btof0->MeanSeg()+1;
      event.deBtof0 = cl_btof0->DeltaE();
      event.Btof0 = cl_btof0->MeanTime() - dst.Time0;
      event.CBtof0 = cl_btof0->CMeanTime() - dst.CTime0;
      dst.Btof0Seg = cl_btof0->MeanSeg()+1;
      dst.deBtof0 = cl_btof0->DeltaE();
      dst.Btof0 = cl_btof0->MeanTime() - dst.Time0;
      dst.CBtof0 = cl_btof0->CMeanTime() - dst.CTime0;
      HF1(102, dst.CBtof0);
      HF2(103, dst.deTime0, dst.CBtof0);
      HF2(104, dst.deBtof0, dst.CBtof0);
    }
  }

  // BH1 with BH2 gate
  {
    Int_t nhbh2 = hodoAna.GetNHits("BH2");
    if(nhbh2 > 0){
      const auto& hit2 = hodoAna.GetHit<BH2Hit>("BH2", 0);
      Int_t    seg2 = hit2->SegmentId()+1;
      Int_t n_mhit2 = hit2->GetEntries();
      for(Int_t m2 = 0; m2<n_mhit2; ++m2){
        Double_t mt2  = hit2->CTime0(m2);
        Int_t    nh   = hodoAna.GetNHits("BH1");
        for(Int_t i=0; i<nh; ++i){
          auto hit1 = hodoAna.GetHit("BH1", i);
          if(!hit1) continue;
          Int_t n_mhit1 = hit1->GetEntries();
          for(Int_t m1 = 0; m1<n_mhit1; ++m1){
            Int_t seg1 = hit1->SegmentId()+1;
            Double_t tu1 = hit1->GetTUp(m1), td1 = hit1->GetTDown(m1);
            Double_t mt1 = hit1->MeanTime(m1);
            HF1(BH1Hid+100*seg1+1100+21+seg2*10, tu1);
            HF1(BH1Hid+100*seg1+1100+22+seg2*10, td1);
            HF1(BH1Hid+100*seg1+1100+23+seg2*10, mt1);
            //For BH1vsBH2 Correlation
            HF2(BH1Hid+100*seg1+2200+21+seg2*10, mt2, tu1);
            HF2(BH1Hid+100*seg1+2200+22+seg2*10, mt2, td1);
            HF2(BH1Hid+100*seg1+2200+23+seg2*10, mt2, mt1);
          }// for(m1)
        }// for(bh1:seg)
      }// for(m2)
    }
    for(Int_t i2=0; i2<nhbh2; ++i2){
      const auto& hit2 = hodoAna.GetHit<BH2Hit>("BH2", i2);
      Int_t    seg2 = hit2->SegmentId()+1;
      Int_t n_mhit2 = hit2->GetEntries();
      for(Int_t m2 = 0; m2<n_mhit2; ++m2){
        Double_t ct0  = hit2->CTime0();
        Double_t t0   = hit2->Time0();
        Int_t nhbh1=hodoAna.GetNHits("BH1");
        for(Int_t i1=0; i1<nhbh1; ++i1){
          const auto& hit1 = hodoAna.GetHit("BH1", i1);
          if(!hit1) continue;
          Int_t n_mhit1 = hit1->GetEntries();
          for(Int_t m1 = 0; m1<n_mhit1; ++m1){
            Int_t seg1=hit1->SegmentId()+1;
            Double_t mt1 = hit1->MeanTime(m1);
            Double_t cmt1 = hit1->CMeanTime(m1);
            Double_t btof = mt1 - t0;
            Double_t cbtof = cmt1 - ct0;
            event.btof[seg1-1][seg2-1] = btof;
            event.cbtof[seg1-1][seg2-1] = cbtof;
            HF1(BH1Hid+100*seg1+1100+24+seg2*10, mt1-ct0);
            if(seg1==5){
              HF1(30+seg2, mt1-ct0);
            }
          }// for(m1)
        }// for(bh1:seg)
      }// for(m2)
    }// for(bh2:seg)
  }

  // BH1-BH2
  {
    Int_t nhbh1 = hodoAna.GetNHits("BH1");
    Int_t nhbh2 = hodoAna.GetNHits("BH2");
    for(Int_t i2=0; i2<nhbh2; ++i2){
      const auto& hit2 = hodoAna.GetHit<BH2Hit>("BH2", i2);
      if(!hit2) continue;

      Int_t n_mhit2 = hit2->GetEntries();
      Int_t    seg2 = hit2->SegmentId()+1;
      Double_t de2  = hit2->DeltaE();
      for(Int_t m2 = 0; m2<n_mhit2; ++m2){
        Double_t t0      = hit2->Time0(m2);

        for(Int_t i1=0; i1<nhbh1; ++i1){
          const auto& hit1 = hodoAna.GetHit("BH1", i1);
          if(!hit1) continue;

          Int_t n_mhit1  = hit1->GetEntries();
          Int_t    seg1  = hit1->SegmentId()+1;
          Double_t de1   = hit1->DeltaE();

          for(Int_t m1 = 0; m1<n_mhit1; ++m1){
            Double_t mt1   = hit1->MeanTime(m1);
            HF1(201, mt1-t0);
            HF2(202, seg1-0.5, seg2-0.5);
            //For BH1vsBH2 Correlation
            HF2(203, t0, mt1);
            HF1(204, mt1);
            HF1(205, t0);
            if(de1>0.5 && de2>0.5){
              HF1(211, mt1-t0);
              HF2(212, seg1-0.5, seg2-0.5);
            }
          }// for(m1)
        }// for(bh1:seg)
      }// for(m2)
    }// for(bh2:seg)
  }

#if 1
  // BH1-BH2 PHC
  {
    Int_t nh1 = hodoAna.GetNHits("BH1");
    Int_t nh2 = hodoAna.GetNHits("BH2");
    for(Int_t i2=0; i2<nh2; ++i2){
      const auto& hit2 = hodoAna.GetHit<BH2Hit>("BH2", i2);
      Int_t     seg2 = hit2->SegmentId()+1;
      Double_t  au2  = hit2->GetAUp(),  ad2  = hit2->GetADown();
      Int_t n_mhit2  = hit2->GetEntries();
      for(Int_t m2 = 0; m2<n_mhit2; ++m2){
        Double_t  tu2  = hit2->GetTUp(m2),  td2  = hit2->GetTDown(m2);
        Double_t  ctu2 = hit2->GetCTUp(m2), ctd2 = hit2->GetCTDown(m2);
        // Double_t  t0   = hit2->Time0();
        Double_t  ct0  = hit2->CTime0(m2);
        Double_t  tofs = ct0-(ctu2+ctd2)/2.;
        for(Int_t i1=0; i1<nh1; ++i1){
          const auto& hit1 = hodoAna.GetHit("BH1", i1);
          Int_t       seg1 = hit1->SegmentId()+1;
          Double_t    au1  = hit1->GetAUp(),  ad1  = hit1->GetADown();
          Int_t    n_mhit1 = hit1->GetEntries();
          for(Int_t m1 = 0; m1<n_mhit1; ++m1){
            Double_t    tu1  = hit1->GetTUp(m1),  td1  = hit1->GetTDown(m1);
            Double_t    ctu1 = hit1->GetCTUp(m1), ctd1 = hit1->GetCTDown(m1);
            Double_t    cmt1 = hit1->CMeanTime(m1);
            // if(trigger_flag[trigger::kBeamA] == 0) continue;
            // if(event.bacnhits > 0) continue;
            HF2(100*seg1+BH1Hid+81, au1, ct0-0.5*(ctu1+ctd1));
            HF2(100*seg1+BH1Hid+82, ad1, ct0-0.5*(ctu1+ctd1));
            HF2(100*seg1+BH1Hid+83, au1, ct0-tu1);
            HF2(100*seg1+BH1Hid+84, ad1, ct0-td1);
            HF2(100*seg2+BH2Hid+81, au2, (cmt1-tofs)-0.5*(ctu2+ctd2));
            HF2(100*seg2+BH2Hid+82, ad2, (cmt1-tofs)-0.5*(ctu2+ctd2));
            HF2(100*seg2+BH2Hid+83, au2, (cmt1-tofs)-tu2);
            HF2(100*seg2+BH2Hid+84, ad2, (cmt1-tofs)-td2);
          }// for(m1)
        }// for(bh1:seg)
      }// for(m2)
    }// for(bh2:seg)
  }
#endif

  // BAC
  hodoAna.DecodeHits("BAC");
  {
    Int_t nh=hodoAna.GetNHits("BAC");
    dst.nhBac = nh;
    HF1(BACHid+10, Double_t(nh));
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit("BAC", i);
      if(!hit) continue;
      Int_t seg = hit->SegmentId()+1;
      HF1(BACHid+11, seg-0.5);
      Double_t a = hit->GetAUp();
      event.bacde[seg-1]  = a;
      Int_t n_mhit = hit->GetEntries();
      for(Int_t m=0; m<n_mhit; ++m){
        Double_t t = hit->GetTUp(m);
        Double_t ct = hit->GetCTUp(m);
        event.bacmt[i][m] = ct;
        HF1(BACHid+100*seg+11, t);
        HF1(BACHid+100*seg+12, a);
        HF1(BACHid+100*seg+13, ct);
      }
    }
  }

  // TOF
  hodoAna.DecodeHits("TOF");
  {
    Int_t nh = hodoAna.GetNHits("TOF");
    HF1(TOFHid+10, Double_t(nh));
    Int_t nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit("TOF", i);
      if(!hit) continue;
      Int_t seg = hit->SegmentId()+1;
      Int_t n_mhit = hit->GetEntries();
      for(Int_t m=0; m<n_mhit; ++m){
        HF1(TOFHid+11, seg-0.5);
        Double_t au  = hit->GetAUp(),   ad  = hit->GetADown();
        Double_t tu  = hit->GetTUp(),   td  = hit->GetTDown();
        Double_t ctu = hit->GetCTUp(),  ctd = hit->GetCTDown();
        Double_t mt  = hit->MeanTime(), cmt = hit->CMeanTime();
        Double_t de  = hit->DeltaE();
        Double_t ude  = hit->UDeltaE();
        Double_t dde  = hit->DDeltaE();
        event.tofmt[seg-1][m] = mt;
        event.tofde[seg-1]    = de;
        event.tofude[seg-1]    = ude;
        event.tofdde[seg-1]    = dde;
        event.tofctu[seg-1][m] = ctu;
        event.tofctd[seg-1][m] = ctd;
        event.tofcmt[seg-1][m] = cmt;
        HF1(TOFHid+100*seg+11, tu);      HF1(TOFHid+100*seg+12, td);
        HF1(TOFHid+100*seg+13, mt);
        HF1(TOFHid+100*seg+17, ctu);     HF1(TOFHid+100*seg+18, ctd);
        HF1(TOFHid+100*seg+19, cmt);     HF1(TOFHid+100*seg+20, ctu-ctd);
        //HF2(TOFHid+100*seg+21, tu, au);  HF2(TOFHid+100*seg+22, td, ad);
        //HF2(TOFHid+100*seg+23, ctu, au); HF2(TOFHid+100*seg+24, ctd, ad);
        HF1(TOFHid+12, cmt);

        dst.utTofSeg[seg-1][m] = tu;
        dst.dtTofSeg[seg-1][m] = td;
        dst.udeTofSeg[seg-1] = au;
        dst.ddeTofSeg[seg-1] = ad;

        if(m == 0){
          HF1(TOFHid+100*seg+14, au);    HF1(TOFHid+100*seg+15, ad);
          HF1(TOFHid+100*seg+16, de);    HF1(TOFHid+13, de);
        }

        if(de>0.5){
          HF1(TOFHid+15, seg-0.5);
          ++nh2;
        }
      }
    }

    HF1(TOFHid+14, Double_t(nh2));
    for(Int_t i1=0; i1<nh; ++i1){
      const auto& hit1 = hodoAna.GetHit("TOF", i1);
      if(!hit1 || hit1->DeltaE()<=0.5) continue;
      Int_t seg1 = hit1->SegmentId()+1;
      for(Int_t i2=0; i2<nh; ++i2){
        if(i1==i2) continue;
        const auto& hit2 = hodoAna.GetHit("TOF", i2);
        if(!hit2 || hit2->DeltaE()<=0.5) continue;
        Int_t seg2 = hit2->SegmentId()+1;

        if(1 == hit1->GetEntries() && 1 == hit2->GetEntries()){
          Double_t ct1 = hit1->CMeanTime(), ct2 = hit2->CMeanTime();
          HF2(TOFHid+21, seg1-0.5, seg2-0.5);
          HF2(TOFHid+22, ct1, ct2);
          HF1(TOFHid+23, ct2-ct1);
          if(std::abs(ct2-ct1)<3.0){
            HF2(TOFHid+24, seg1-0.5, seg2-0.5);
          }
        }//for(i2)
      }//for(i1)
    }

    Int_t nc = hodoAna.GetNClusters("TOF");
    HF1(TOFHid+30, Double_t(nc));
    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster("TOF", i);
      if(!cl) continue;
      Int_t cs = cl->ClusterSize();
      Double_t ms  = cl->MeanSeg()+1;
      Double_t cmt = cl->CMeanTime();
      Double_t de  = cl->DeltaE();
      HF1(TOFHid+31, Double_t(cs));
      HF1(TOFHid+32, ms-0.5);
      HF1(TOFHid+33, cmt); HF1(TOFHid+34, de);
    }
  }

  ////////// Dst
  {
    Int_t nc = hodoAna.GetNClusters("BH1");
    dst.nhBh1 = nc;
    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster("BH1", i);
      if(!cl) continue;
      dst.csBh1[i]  = cl->ClusterSize();
      dst.Bh1Seg[i] = cl->MeanSeg()+1;
      dst.tBh1[i]   = cl->CMeanTime();
      dst.dtBh1[i]  = cl->TimeDiff();
      dst.deBh1[i]  = cl->DeltaE();
      if(cl->ClusterSize() >1) cl->Print();

      Int_t nc2 = hodoAna.GetNClusters("BH2");
      for(Int_t i2=0; i2<nc2; ++i2){
        const auto& cl2 = hodoAna.GetCluster("BH2", i2);
        if(!cl2) continue;
        Double_t btof = cl->MeanTime() - dst.Time0;
        Double_t cbtof = cl->CMeanTime() - dst.CTime0;
        dst.btof[i] = btof;
        dst.cbtof[i] = cbtof;
      }
    }
  }

  {
    Int_t nc = hodoAna.GetNClusters("BH2");
    dst.nhBh2 = nc;
    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster<BH2Cluster>("BH2", i);
      if(!cl) continue;
      dst.csBh2[i]  = cl->ClusterSize();
      dst.Bh2Seg[i] = cl->MeanSeg()+1;
      dst.tBh2[i]   = cl->CMeanTime();
      dst.t0Bh2[i]  = cl->CTime0();
      dst.dtBh2[i]  = cl->TimeDiff();
      dst.deBh2[i]  = cl->DeltaE();
#if HodoHitPos
      dst.posBh2[i]  = 0.5*PropVelBH2*cl->TimeDiff();
#endif
    }
  }

  {
    Int_t nc = hodoAna.GetNClusters("BAC");
    dst.nhBac = nc;
    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster("BAC", i);
      if(!cl) continue;
      dst.csBac[i]  = cl->ClusterSize();
      dst.BacSeg[i] = cl->MeanSeg()+1;
      dst.tBac[i]   = cl->CMeanTime();
      dst.deBac[i]  = cl->DeltaE();
    }
  }

  {
    Int_t nc = hodoAna.GetNClusters("TOF");
    dst.nhTof = nc;
    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster("TOF", i);
      if(!cl) continue;
      dst.csTof[i]  = cl->ClusterSize();
      dst.TofSeg[i] = cl->MeanSeg()+1;
      dst.tTof[i]   = cl->CMeanTime();
      dst.dtTof[i]  = cl->TimeDiff();
      dst.deTof[i]  = cl->DeltaE();
    }
  }

  hodoAna.DecodeHits("AC1");
  {
    Int_t nc = hodoAna.GetNClusters("AC1");
    dst.nhAc1 = nc;
    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster("AC1", i);
      if(!cl) continue;
      dst.csAc1[i]  = cl->ClusterSize();
      dst.Ac1Seg[i] = cl->MeanSeg()+1;
      dst.tAc1[i]   = cl->CMeanTime();
    }
  }

  hodoAna.DecodeHits("WC");
  {
    Int_t nc = hodoAna.GetNClusters("WC");
    dst.nhWc = nc;
    for(Int_t i=0; i<nc; ++i){
      const auto& cl = hodoAna.GetCluster("WC", i);
      if(!cl) continue;
      dst.csWc[i]  = cl->ClusterSize();
      dst.WcSeg[i] = cl->MeanSeg()+1;
      dst.tWc[i]   = cl->CMeanTime();
      dst.deWc[i]  = cl->DeltaE();
    }
  }


#if 0
  // BH1 (for parameter tuning)
  if(dst.Time0Seg==4){
    const HodoRHitContainer& cont = rawData.GetBH1HodoRawHitContainer();
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      auto hit = cont[i];
      Int_t seg = hit->SegmentId()+1;

      //Up
      {
        Int_t n_mhit = hit->GetSizeTdcUp();
        for(Int_t m=0; m<n_mhit; ++m){
          Int_t T = hit->GetTdcUp(m);
          if(T > 0) HF1(BH1Hid +100*seg +9, Double_t(T));
        }// for(m)
      }

      //Down
      {
        Int_t n_mhit = hit->GetSizeTdcDown();
        for(Int_t m=0; m<n_mhit; ++m){
          Int_t T = hit->GetTdcDown(m);
          if(T > 0) HF1(BH1Hid +100*seg +10, Double_t(T));
        }// for(m)
      }
    }
  }// (Raw BH1 TDC for parameter tuning)
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingEnd()
{
  tree->Fill();
  hodo->Fill();
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
    HB1(10+i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000);
  }

  // BH1
  // Rawdata
  HB1(BH1Hid +0, "#Hits BH1",        NumOfSegBH1+1, 0., Double_t(NumOfSegBH1+1));
  HB1(BH1Hid +1, "Hitpat BH1",       NumOfSegBH1,   0., Double_t(NumOfSegBH1));
  HB1(BH1Hid +2, "#Hits BH1(Tor)",   NumOfSegBH1+1, 0., Double_t(NumOfSegBH1+1));
  HB1(BH1Hid +3, "Hitpat BH1(Tor)",  NumOfSegBH1,   0., Double_t(NumOfSegBH1));
  HB1(BH1Hid +4, "#Hits BH1(Tand)",  NumOfSegBH1+1, 0., Double_t(NumOfSegBH1+1));
  HB1(BH1Hid +5, "Hitpat BH1(Tand)", NumOfSegBH1,   0., Double_t(NumOfSegBH1));

  for(Int_t i=1; i<=NumOfSegBH1; ++i){
    TString title1 = Form("BH1-%d UpAdc", i);
    TString title2 = Form("BH1-%d DownAdc", i);
    TString title3 = Form("BH1-%d UpTdc", i);
    TString title4 = Form("BH1-%d DownTdc", i);
    TString title5 = Form("BH1-%d UpAdc(w Tdc)", i);
    TString title6 = Form("BH1-%d DownAdc(w Tdc)", i);
    TString title7 = Form("BH1-%d UpAdc(w/o Tdc)", i);
    TString title8 = Form("BH1-%d DownAdc(w/o Tdc)", i);
    TString title9 = Form("BH1-%d UpTdc (Time0Seg==4)", i);
    TString title10= Form("BH1-%d DownTdc (Time0Seg==4)", i);
    HB1(BH1Hid +100*i +1, title1, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH1Hid +100*i +2, title2, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH1Hid +100*i +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(BH1Hid +100*i +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(BH1Hid +100*i +5, title5, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH1Hid +100*i +6, title6, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH1Hid +100*i +7, title7, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH1Hid +100*i +8, title8, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH1Hid +100*i +9, title9, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(BH1Hid +100*i +10,title10,NbinTdcHr, MinTdcHr, MaxTdcHr);
  }

  //BH1 Normalized
  HB1(BH1Hid +10, "#Hits BH1[Hodo]",  NumOfSegBH1+1, 0., Double_t(NumOfSegBH1+1));
  HB1(BH1Hid +11, "Hitpat BH1[Hodo]", NumOfSegBH1,   0., Double_t(NumOfSegBH1));
  HB1(BH1Hid +12, "CMeanTime BH1", 200, -10., 10.);
  HB1(BH1Hid +13, "dE BH1", 200, -0.5, 4.5);
  HB1(BH1Hid +14, "#Hits BH1[HodoGood]",  NumOfSegBH1+1, 0., Double_t(NumOfSegBH1+1));
  HB1(BH1Hid +15, "Hitpat BH1[HodoGood]", NumOfSegBH1,   0., Double_t(NumOfSegBH1));
  HB1(BH1Hid +16, "CMeanTime BH1[HodoGood]", 200, -10., 10.);

  HB2(BH1Hid +21, "BH1HitPat%BH1HitPat[HodoGood]",
      NumOfSegBH1, 0., Double_t(NumOfSegBH1),
      NumOfSegBH1, 0., Double_t(NumOfSegBH1));
  HB2(BH1Hid +22, "CMeanTimeBH1%CMeanTimeBH1[HodoGood]",
      100, -5., 5., 100, -5., 5.);
  HB1(BH1Hid +23, "TDiff BH1[HodoGood]", 200, -10., 10.);
  HB2(BH1Hid +24, "BH1HitPat%BH1HitPat[HodoGood2]",
      NumOfSegBH1, 0., Double_t(NumOfSegBH1),
      NumOfSegBH1, 0., Double_t(NumOfSegBH1));
  HB1(BH1Hid +30, "#Clusters BH1", NumOfSegBH1+1, 0., Double_t(NumOfSegBH1+1));
  HB1(BH1Hid +31, "ClusterSize BH1", 5, 0., 5.);
  HB1(BH1Hid +32, "HitPat Cluster BH1", 2*NumOfSegBH1, 0., Double_t(NumOfSegBH1));
  HB1(BH1Hid +33, "CMeamTime Cluster BH1", 200, -10., 10.);
  HB1(BH1Hid +34, "DeltaE Cluster BH1", 100, -0.5, 4.5);
  HB1(BH1Hid +35, "#Clusters BH1(AdcGood)", NumOfSegBH1+1, 0., Double_t(NumOfSegBH1+1));
  HB1(BH1Hid +36, "CMeamTime Cluster BH1(AdcGood)", 200, -10., 10.);

  HB2(BH1Hid +41, "BH1ClP%BH1ClP",
      NumOfSegBH1, 0., Double_t(NumOfSegBH1),
      NumOfSegBH1, 0., Double_t(NumOfSegBH1));
  HB2(BH1Hid +42, "CMeanTimeBH1%CMeanTimeBH1[Cluster]",
      100, -5., 5., 100, -5., 5.);
  HB1(BH1Hid +43, "TDiff BH1[Cluster]", 200, -10., 10.);
  HB2(BH1Hid +44, "BH1ClP%BH1ClP[AdcGood]",
      NumOfSegBH1, 0., Double_t(NumOfSegBH1),
      NumOfSegBH1, 0., Double_t(NumOfSegBH1));

  for(Int_t i=1; i<=NumOfSegBH1; ++i){
    TString title11 = Form("BH1-%d Up Time", i);
    TString title12 = Form("BH1-%d Down Time", i);
    TString title13 = Form("BH1-%d MeanTime", i);
    TString title14 = Form("BH1-%d Up dE", i);
    TString title15 = Form("BH1-%d Down dE", i);
    TString title16 = Form("BH1-%d dE", i);
    TString title17 = Form("BH1-%d Up CTime", i);
    TString title18 = Form("BH1-%d Down CTime", i);
    TString title19 = Form("BH1-%d CMeanTime", i);
    TString title20 = Form("BH1-%d Tup-Tdown", i);
    TString title21 = Form("BH1-%d Up dE%%Time", i);
    TString title22 = Form("BH1-%d Down dE%%Time", i);
    TString title23 = Form("BH1-%d Up dE%%CTime", i);
    TString title24 = Form("BH1-%d Down dE%%CTime", i);
    TString title31 = Form("BH1-%d Up Time%%dE", i);
    TString title32 = Form("BH1-%d Down Time%%dE", i);
    TString title33 = Form("BH1-%d Up CTime%%dE", i);
    TString title34 = Form("BH1-%d Down CTime%%dE", i);
    HB1(BH1Hid +100*i +11, title11, 200, -10., 10.);
    HB1(BH1Hid +100*i +12, title12, 200, -10., 10.);
    HB1(BH1Hid +100*i +13, title13, 200, -10., 10.);
    HB1(BH1Hid +100*i +14, title14, 200, -0.5, 4.5);
    HB1(BH1Hid +100*i +15, title15, 200, -0.5, 4.5);
    HB1(BH1Hid +100*i +16, title16, 200, -0.5, 4.5);
    HB1(BH1Hid +100*i +17, title17, 200, -10., 10.);
    HB1(BH1Hid +100*i +18, title18, 200, -10., 10.);
    HB1(BH1Hid +100*i +19, title19, 200, -10., 10.);
    HB1(BH1Hid +100*i +20, title20, 200, -5.0, 5.0);
    HB2(BH1Hid +100*i +21, title21, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH1Hid +100*i +22, title22, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH1Hid +100*i +23, title23, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH1Hid +100*i +24, title24, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH1Hid +100*i +31, title31, 100, -0.5, 4.5, 1000,-10.,10.);
    HB2(BH1Hid +100*i +32, title32, 100, -0.5, 4.5, 1000,-10.,10.);
    HB2(BH1Hid +100*i +33, title33, 100, -0.5, 4.5, 1000,-10.,10.);
    HB2(BH1Hid +100*i +34, title34, 100, -0.5, 4.5, 1000,-10.,10.);

    for(Int_t j=1; j<=NumOfSegBH2; ++j){
      TString title1 = Form("BH1-%d Up Time [BH2-%d]", i, j);
      TString title2 = Form("BH1-%d Down Time [BH2-%d]", i, j);
      TString title3 = Form("BH1-%d MeanTime [BH2-%d]", i, j);
      TString title4 = Form("BH1-%d MeanTime-BH2MeamTime [BH2-%d]", i, j);
      HB1(BH1Hid +100*i +1100 +21 +j*10, title1, 200, -10., 10.);
      HB1(BH1Hid +100*i +1100 +22 +j*10, title2, 200, -10., 10.);
      HB1(BH1Hid +100*i +1100 +23 +j*10, title3, 200, -10., 10.);
      HB1(BH1Hid +100*i +1100 +24 +j*10, title4, 200, -10., 10.);
      //For BH1vsBH2 Correlation
      HB2(BH1Hid +100*i +2200 +21 +j*10, Form("BH1-%d Up Time %% BH2-%d MeanTime", i, j),
	  100, -10., 10., 100, -10., 10.);
      HB2(BH1Hid +100*i +2200 +22 +j*10, Form("BH1-%d Down Time %% BH2-%d MeanTime", i, j),
	  100, -10., 10., 100, -10., 10.);
      HB2(BH1Hid +100*i +2200 +23 +j*10, Form("BH1-%d MeanTime %% BH2-%d MeanTime", i, j),
	  100, -10., 10., 100, -10., 10.);
    }
  }

  // BH2
  HB1(BH2Hid +0, "#Hits BH2",        NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +1, "Hitpat BH2",       NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB1(BH2Hid +2, "#Hits BH2(Tor)",   NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +3, "Hitpat BH2(Tor)",  NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB1(BH2Hid +4, "#Hits BH2(Tand)",  NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +5, "Hitpat BH2(Tand)", NumOfSegBH2,   0., Double_t(NumOfSegBH2));

  for(Int_t i=1; i<=NumOfSegBH2; ++i){
    TString title1 = Form("BH2-%d UpAdc", i);
    TString title2 = Form("BH2-%d DownAdc", i);
    TString title3 = Form("BH2-%d UpTdc", i);
    TString title4 = Form("BH2-%d DownTdc", i);
    TString title5 = Form("BH2-%d UpAdc(w Tdc)", i);
    TString title6 = Form("BH2-%d DownAdc(w Tdc)", i);
    TString title7 = Form("BH2-%d UpAdc(w/o Tdc)", i);
    TString title8 = Form("BH2-%d DownAdc(w/o Tdc)", i);
    HB1(BH2Hid +100*i +1, title1, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH2Hid +100*i +2, title2, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH2Hid +100*i +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(BH2Hid +100*i +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(BH2Hid +100*i +5, title5, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH2Hid +100*i +6, title6, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH2Hid +100*i +7, title7, NbinAdc,   MinAdc,   MaxAdc);
    HB1(BH2Hid +100*i +8, title8, NbinAdc,   MinAdc,   MaxAdc);
  }
  HB1(BH2Hid +10, "#Hits BH2[Hodo]",  NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +11, "Hitpat BH2[Hodo]", NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB1(BH2Hid +12, "CMeanTime BH2", 200, -10., 10.);
  HB1(BH2Hid +13, "dE BH2", 200, -0.5, 4.5);
  HB1(BH2Hid +14, "#Hits BH2[HodoGood]",  NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +15, "Hitpat BH2[HodoGood]", NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB1(BH2Hid +16, "CMeanTime BH2[HodoGood]", 200, -10., 10.);

  for(Int_t i=1; i<=NumOfSegBH2; ++i){
    TString title11 = Form("BH2-%d Up Time", i);
    TString title12 = Form("BH2-%d Down Time", i);
    TString title13 = Form("BH2-%d MeanTime", i);
    TString title14 = Form("BH2-%d Up dE", i);
    TString title15 = Form("BH2-%d Down dE", i);
    TString title16 = Form("BH2-%d dE", i);
    TString title17 = Form("BH2-%d Up CTime", i);
    TString title18 = Form("BH2-%d Down CTime", i);
    TString title19 = Form("BH2-%d CMeanTime", i);
    TString title20 = Form("BH2-%d Tup-Tdown", i);
    TString title21 = Form("BH2-%d Up Time0", i);
    TString title22 = Form("BH2-%d Down Time0", i);
    TString title23 = Form("BH2-%d Up CTime", i);
    TString title24 = Form("BH2-%d Down CTime", i);
    TString title25 = Form("BH2-%d MeanTime0", i);
    TString title26 = Form("BH2-%d CMeanTime0", i);
    TString title27 = Form("BH2-%d Up dE%%Time", i);
    TString title28 = Form("BH2-%d Down dE%%Time", i);
    TString title29 = Form("BH2-%d Up dE%%CTime", i);
    TString title30 = Form("BH2-%d Down dE%%CTime", i);
    HB1(BH2Hid +100*i +11, title11, 200, -10., 10.);
    HB1(BH2Hid +100*i +12, title12, 200, -10., 10.);
    HB1(BH2Hid +100*i +13, title13, 200, -10., 10.);
    HB1(BH2Hid +100*i +14, title14, 200, -0.5, 4.5);
    HB1(BH2Hid +100*i +15, title15, 200, -0.5, 4.5);
    HB1(BH2Hid +100*i +16, title16, 200, -0.5, 4.5);
    HB1(BH2Hid +100*i +17, title17, 200, -10., 10.);
    HB1(BH2Hid +100*i +18, title18, 200, -10., 10.);
    HB1(BH2Hid +100*i +19, title19, 200, -10., 10.);
    HB1(BH2Hid +100*i +20, title20, 200, -5.0, 5.0);
    HB1(BH2Hid +100*i +21, title21, 200, -10., 10.);
    HB1(BH2Hid +100*i +22, title22, 200, -10., 10.);
    HB1(BH2Hid +100*i +23, title23, 200, -10., 10.);
    HB1(BH2Hid +100*i +24, title24, 200, -10., 10.);
    HB1(BH2Hid +100*i +25, title25, 200, -10., 10.);
    HB1(BH2Hid +100*i +26, title26, 200, -10., 10.);
    HB2(BH2Hid +100*i +27, title27, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH2Hid +100*i +28, title28, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH2Hid +100*i +29, title29, 100, -10., 10., 100, -0.5, 4.5);
    HB2(BH2Hid +100*i +30, title30, 100, -10., 10., 100, -0.5, 4.5);
  }

  HB2(BH2Hid +21, "BH2HitPat%BH2HitPat[HodoGood]", NumOfSegBH2,   0., Double_t(NumOfSegBH2),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB2(BH2Hid +22, "CMeanTimeBH2%CMeanTimeBH2[HodoGood]",
      100, -2.5, 2.5, 100, -2.5, 2.5);
  HB1(BH2Hid +23, "TDiff BH2[HodoGood]", 200, -10., 10.);
  HB2(BH2Hid +24, "BH2HitPat%BH2HitPat[HodoGood2]", NumOfSegBH2,   0., Double_t(NumOfSegBH2),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));

  HB1(BH2Hid +30, "#Clusters BH2", NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +31, "ClusterSize BH2", 5, 0., 5.);
  HB1(BH2Hid +32, "HitPat Cluster BH2", 2*NumOfSegBH2, 0., Double_t(NumOfSegBH2));
  HB1(BH2Hid +33, "CMeamTime Cluster BH2", 200, -10., 10.);
  HB1(BH2Hid +34, "DeltaE Cluster BH2", 100, -0.5, 4.5);
  HB1(BH2Hid +35, "#Clusters BH2(ADCGood)", NumOfSegBH2+1, 0., Double_t(NumOfSegBH2+1));
  HB1(BH2Hid +36, "CMeamTime Cluster BH2(ADCGood)", 200, -10., 10.);

  HB2(BH2Hid +41, "BH2ClP%BH2ClP", NumOfSegBH2,   0., Double_t(NumOfSegBH2),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  HB2(BH2Hid +42, "CMeanTimeBH2%CMeanTimeBH2[Cluster]",
      100, -2.5, 2.5, 100, -2.5, 2.5);
  HB1(BH2Hid +43, "TDiff BH2[Cluster]", 200, -10., 10.);
  HB2(BH2Hid +44, "BH2ClP%BH2ClP(ADCGood)", NumOfSegBH2,   0., Double_t(NumOfSegBH2),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));

  HB1(201, "TimeDiff BH1-BH2", 200, -10., 10.);
  HB2(202, "SegBH2%SegBH1",NumOfSegBH1,   0., Double_t(NumOfSegBH1),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));
  //For BH1vsBH2 Corr
  HB2(203, "MTBH2%MTBH1", 200, -10., 10., 200, -10., 10.);
  HB1(204, "MTBH1", 200, -10., 10.);
  HB1(205, "MTBH2", 200, -10., 10.);

  HB1(211, "TimeDiff BH1-BH2(GoodAdc)", 200, -10., 10.);
  HB2(212, "SegBH2%SegBH1(GoodAdc)",NumOfSegBH1,   0., Double_t(NumOfSegBH1),
      NumOfSegBH2,   0., Double_t(NumOfSegBH2));

  // BH1-BH2 PHC
  for(Int_t i=1; i<=NumOfSegBH1; ++i){
    TString title1 = Form("BH1-%dU CT-TOF%%dE", i);
    TString title2 = Form("BH1-%dD CT-TOF%%dE", i);
    TString title3 = Form("BH1-%dU  T-TOF%%dE", i);
    TString title4 = Form("BH1-%dD  T-TOF%%dE", i);
    HB2(BH1Hid +100*i +81, title1, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH1Hid +100*i +82, title2, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH1Hid +100*i +83, title3, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH1Hid +100*i +84, title4, 200, -0.5, 4.5, 200, -10., 10.);
  }
  for(Int_t i=1; i<=NumOfSegBH2; ++i){
    TString title1 = Form("BH2-%dU CT-TOF%%dE", i);
    TString title2 = Form("BH2-%dD CT-TOF%%dE", i);
    TString title3 = Form("BH2-%dU  T-TOF%%dE", i);
    TString title4 = Form("BH2-%dD  T-TOF%%dE", i);
    HB2(BH2Hid +100*i +81, title1, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH2Hid +100*i +82, title2, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH2Hid +100*i +83, title3, 200, -0.5, 4.5, 200, -10., 10.);
    HB2(BH2Hid +100*i +84, title4, 200, -0.5, 4.5, 200, -10., 10.);
  }

  // BTOF
  HB1(100, "BH2 MeanTime0", 400, -4, 4);
  HB1(101, "CTime0", 400, -4, 4);
  HB1(102, "CBtof0", 400, -4, 4);
  HB2(103, "CBtof0%deTime0", 200, 0, 4, 200, -4, 4);
  HB2(104, "CBtof0%deBTof0", 200, 0, 4, 200, -4, 4);

  // BAC
  HB1(BACHid +0, "#Hits BAC",        NumOfSegBAC+1, 0., Double_t(NumOfSegBAC+1));
  HB1(BACHid +1, "Hitpat BAC",       NumOfSegBAC,   0., Double_t(NumOfSegBAC));
  HB1(BACHid +2, "#Hits BAC(Tor)",   NumOfSegBAC+1, 0., Double_t(NumOfSegBAC+1));
  HB1(BACHid +3, "Hitpat BAC(Tor)",  NumOfSegBAC,   0., Double_t(NumOfSegBAC));
  HB1(BACHid +4, "#Hits BAC(Tand)",  NumOfSegBAC+1, 0., Double_t(NumOfSegBAC+1));
  HB1(BACHid +5, "Hitpat BAC(Tand)", NumOfSegBAC,   0., Double_t(NumOfSegBAC));
  for(Int_t i=1; i<=NumOfSegBAC; ++i){
    TString title1 = Form("BAC-%d UpAdc", i);
    TString title3 = Form("BAC-%d UpTdc", i);
    TString title5 = Form("BAC-%d UpAdc(w Tdc)", i);
    TString title7 = Form("BAC-%d UpAdc(w/o Tdc)", i);
    HB1(BACHid +100*i +1, title1, NbinTdc, MinTdc, MaxTdc);
    HB1(BACHid +100*i +3, title3, NbinTdc, MinTdc, MaxTdc);
    HB1(BACHid +100*i +5, title5, NbinAdc, MinAdc, MaxAdc);
    HB1(BACHid +100*i +7, title7, NbinAdc, MinAdc, MaxAdc);
  }
  HB1(BACHid +10, "#Hits BAC[Hodo]",     NumOfSegBAC+1, 0., Double_t(NumOfSegBAC+1));
  HB1(BACHid +11, "Hitpat BAC[Hodo]",    NumOfSegBAC,   0., Double_t(NumOfSegBAC));
  for(Int_t i=1; i<=NumOfSegBAC; ++i){
    TString title1 = Form("BAC-%d Time", i);
    TString title3 = Form("BAC-%d dE", i);
    TString title5 = Form("BAC-%d CTime", i);
    HB1(BACHid +100*i +11, title1, 500, -5., 45.);
    HB1(BACHid +100*i +12, title3, 200, -0.5, 4.5);
    HB1(BACHid +100*i +13, title5, 500, -5., 45.);
  }

  // TOF
  HB1(TOFHid +0, "#Hits TOF",        NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +1, "Hitpat TOF",       NumOfSegTOF,   0., Double_t(NumOfSegTOF));
  HB1(TOFHid +2, "#Hits TOF(Tor)",   NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +3, "Hitpat TOF(Tor)",  NumOfSegTOF,   0., Double_t(NumOfSegTOF));
  HB1(TOFHid +4, "#Hits TOF(Tand)",  NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +5, "Hitpat TOF(Tand)", NumOfSegTOF,   0., Double_t(NumOfSegTOF));

  for(Int_t i=1; i<=NumOfSegTOF; ++i){
    TString title1 = Form("TOF-%d UpAdc", i);
    TString title2 = Form("TOF-%d DownAdc", i);
    TString title3 = Form("TOF-%d UpTdc", i);
    TString title4 = Form("TOF-%d DownTdc", i);
    TString title5 = Form("TOF-%d UpAdc(w Tdc)", i);
    TString title6 = Form("TOF-%d DownAdc(w Tdc)", i);
    TString title7 = Form("TOF-%d UpAdc(w/o Tdc)", i);
    TString title8 = Form("TOF-%d DownAdc(w/o Tdc)", i);
    HB1(TOFHid +100*i +1, title1, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*i +2, title2, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*i +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(TOFHid +100*i +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(TOFHid +100*i +5, title5, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*i +6, title6, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*i +7, title7, NbinAdc, MinAdc, MaxAdc);
    HB1(TOFHid +100*i +8, title8, NbinAdc, MinAdc, MaxAdc);
  }

  HB1(TOFHid +10, "#Hits Tof[Hodo]",  NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +11, "Hitpat Tof[Hodo]", NumOfSegTOF,   0., Double_t(NumOfSegTOF));
  HB1(TOFHid +12, "CMeanTime Tof", 500, -5., 45.);
  HB1(TOFHid +13, "dE Tof", 200, -0.5, 4.5);
  HB1(TOFHid +14, "#Hits Tof[HodoGood]",  NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +15, "Hitpat Tof[HodoGood]", NumOfSegTOF,   0., Double_t(NumOfSegTOF));

  for(Int_t i=1; i<=NumOfSegTOF; ++i){
    TString title11 = Form("TOF-%d Up Time", i);
    TString title12 = Form("TOF-%d Down Time", i);
    TString title13 = Form("TOF-%d MeanTime", i);
    TString title14 = Form("TOF-%d Up dE", i);
    TString title15 = Form("TOF-%d Down dE", i);
    TString title16 = Form("TOF-%d dE", i);
    TString title17 = Form("TOF-%d Up CTime", i);
    TString title18 = Form("TOF-%d Down CTime", i);
    TString title19 = Form("TOF-%d CMeanTime", i);
    TString title20 = Form("TOF-%d Tup-Tdown", i);
    TString title21 = Form("TOF-%d dE (w/ TOF-HT)", i);
    HB1(TOFHid +100*i +11, title11, 500, -5., 45.);
    HB1(TOFHid +100*i +12, title12, 500, -5., 45.);
    HB1(TOFHid +100*i +13, title13, 500, -5., 45.);
    HB1(TOFHid +100*i +14, title14, 200, -0.5, 4.5);
    HB1(TOFHid +100*i +15, title15, 200, -0.5, 4.5);
    HB1(TOFHid +100*i +16, title16, 200, -0.5, 4.5);
    HB1(TOFHid +100*i +17, title17, 500, -5., 45.);
    HB1(TOFHid +100*i +18, title18, 500, -5., 45.);
    HB1(TOFHid +100*i +19, title19, 500, -5., 45.);
    HB1(TOFHid +100*i +20, title20, 200, -10.0, 10.0);
    HB1(TOFHid +100*i +21, title21, 200, -0.5, 4.5);
  }

  HB2(TOFHid +21, "TofHitPat%TofHitPat[HodoGood]", NumOfSegTOF,   0., Double_t(NumOfSegTOF),
      NumOfSegTOF,   0., Double_t(NumOfSegTOF));
  HB2(TOFHid +22, "CMeanTimeTof%CMeanTimeTof[HodoGood]",
      120, 10., 40., 120, 10., 40.);
  HB1(TOFHid +23, "TDiff Tof[HodoGood]", 200, -10., 10.);
  HB2(TOFHid +24, "TofHitPat%TofHitPat[HodoGood2]", NumOfSegTOF,   0., Double_t(NumOfSegTOF),
      NumOfSegTOF,   0., Double_t(NumOfSegTOF));

  HB1(TOFHid +30, "#Clusters Tof", NumOfSegTOF+1, 0., Double_t(NumOfSegTOF+1));
  HB1(TOFHid +31, "ClusterSize Tof", 5, 0., 5.);
  HB1(TOFHid +32, "HitPat Cluster Tof", 2*NumOfSegTOF, 0., Double_t(NumOfSegTOF));
  HB1(TOFHid +33, "CMeamTime Cluster Tof", 500, -5., 45.);
  HB1(TOFHid +34, "DeltaE Cluster Tof", 100, -0.5, 4.5);

  // AC1
  HB1(AC1Hid +0, "#Hits AC1",        NumOfSegAC1+1, 0., Double_t(NumOfSegAC1+1));
  HB1(AC1Hid +1, "Hitpat AC1",       NumOfSegAC1,   0., Double_t(NumOfSegAC1));
  HB1(AC1Hid +2, "#Hits AC1(Tor)",   NumOfSegAC1+1, 0., Double_t(NumOfSegAC1+1));
  HB1(AC1Hid +3, "Hitpat AC1(Tor)",  NumOfSegAC1,   0., Double_t(NumOfSegAC1));
  HB1(AC1Hid +4, "#Hits AC1(Tand)",  NumOfSegAC1+1, 0., Double_t(NumOfSegAC1+1));
  HB1(AC1Hid +5, "Hitpat AC1(Tand)", NumOfSegAC1,   0., Double_t(NumOfSegAC1));
  for(Int_t i=1; i<=NumOfSegAC1; ++i){
    TString title1 = Form("AC1-%d Adc", i);
    TString title3 = Form("AC1-%d Tdc", i);
    TString title5 = Form("AC1-%d Adc(w Tdc)", i);
    TString title7 = Form("AC1-%d Adc(w/o Tdc)", i);
    HB1(AC1Hid +100*i +1, title1, NbinAdc, MinAdc, MaxAdc);
    HB1(AC1Hid +100*i +3, title3, NbinTdc, MinTdc, MaxTdc);
    HB1(AC1Hid +100*i +5, title5, NbinAdc, MinAdc, MaxAdc);
    HB1(AC1Hid +100*i +7, title7, NbinAdc, MinAdc, MaxAdc);
  }
  HB1(AC1Hid +10, "#Hits Ac1[Hodo]",  NumOfSegAC1+1, 0., Double_t(NumOfSegAC1+1));
  HB1(AC1Hid +11, "Hitpat Ac1[Hodo]", NumOfSegAC1,   0., Double_t(NumOfSegAC1));
  HB1(AC1Hid +12, "CMeanTime Ac1", 500, -5., 45.);
  HB1(AC1Hid +14, "#Hits Ac1[HodoGood]",  NumOfSegAC1+1, 0., Double_t(NumOfSegAC1+1));
  HB1(AC1Hid +15, "Hitpat Ac1[HodoGood]", NumOfSegAC1,   0., Double_t(NumOfSegAC1));
  for(Int_t i=1; i<=NumOfSegAC1; ++i){
    TString title11 = Form("AC1-%d Time", i);
    HB1(AC1Hid +100*i +11, title11, 500, -5., 45.);
  }

  //WC
  HB1(WCHid +0, "#Hits WC",        NumOfSegWC+1, 0., Double_t(NumOfSegWC+1));
  HB1(WCHid +1, "Hitpat WC",       NumOfSegWC,   0., Double_t(NumOfSegWC));
  HB1(WCHid +2, "#Hits WC(Tor)",   NumOfSegWC+1, 0., Double_t(NumOfSegWC+1));
  HB1(WCHid +3, "Hitpat WC(Tor)",  NumOfSegWC,   0., Double_t(NumOfSegWC));
  HB1(WCHid +4, "#Hits WC(Tand)",  NumOfSegWC+1, 0., Double_t(NumOfSegWC+1));
  HB1(WCHid +5, "Hitpat WC(Tand)", NumOfSegWC,   0., Double_t(NumOfSegWC));

  for(Int_t i=1; i<=NumOfSegWC; ++i){
    TString title1 = Form("WC-%d UpAdc", i);
    TString title2 = Form("WC-%d DownAdc", i);
    TString title3 = Form("WC-%d UpTdc", i);
    TString title4 = Form("WC-%d DownTdc", i);
    TString title5 = Form("WC-%d UpAdc(w Tdc)", i);
    TString title6 = Form("WC-%d DownAdc(w Tdc)", i);
    TString title7 = Form("WC-%d UpAdc(w/o Tdc)", i);
    TString title8 = Form("WC-%d DownAdc(w/o Tdc)", i);
    HB1(WCHid +100*i +1, title1, NbinAdc,   MinAdc,   MaxAdc);
    HB1(WCHid +100*i +2, title2, NbinAdc,   MinAdc,   MaxAdc);
    HB1(WCHid +100*i +3, title3, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(WCHid +100*i +4, title4, NbinTdcHr, MinTdcHr, MaxTdcHr);
    HB1(WCHid +100*i +5, title5, NbinAdc,   MinAdc,   MaxAdc);
    HB1(WCHid +100*i +6, title6, NbinAdc,   MinAdc,   MaxAdc);
    HB1(WCHid +100*i +7, title7, NbinAdc,   MinAdc,   MaxAdc);
    HB1(WCHid +100*i +8, title8, NbinAdc,   MinAdc,   MaxAdc);
  }

  //WCSUM
  HB1(WCSUMHid +0, "#Hits WC SUM", NumOfSegWC+1, 0., Double_t(NumOfSegWC+1));
  HB1(WCSUMHid +1, "Hitpat WC SUM", NumOfSegWC, 0., Double_t(NumOfSegWC));
  for(Int_t i=1; i<=NumOfSegWC; ++i){
    TString title1 = Form("WCSUM-%d Adc", i);
    TString title3 = Form("WCSUM-%d Tdc", i);
    TString title5 = Form("WCSUM-%d Adc(w Tdc)", i);
    TString title7 = Form("WCSUM-%d Adc(w/o Tdc)", i);
    HB1(WCSUMHid +100*i +1, title1, NbinAdc, MinAdc, MaxAdc);
    HB1(WCSUMHid +100*i +3, title3, NbinTdc, MinTdc, MaxTdc);
    HB1(WCSUMHid +100*i +5, title5, NbinAdc, MinAdc, MaxAdc);
    HB1(WCSUMHid +100*i +7, title7, NbinAdc, MinAdc, MaxAdc);
  }
  HB1(WCSUMHid +10, "#Hits WcSum[Hodo]",  NumOfSegWC+1, 0., Double_t(NumOfSegWC+1));
  HB1(WCSUMHid +11, "Hitpat WcSum[Hodo]", NumOfSegWC,   0., Double_t(NumOfSegWC)  );
  HB1(WCSUMHid +13, "dE WcSum", 200, -0.5, 4.5);
  HB1(WCSUMHid +14, "#Hits WcSum[HodoGood]",  NumOfSegWC+1, 0., Double_t(NumOfSegWC+1));
  HB1(WCSUMHid +15, "Hitpat WcSum[HodoGood]", NumOfSegWC,   0., Double_t(NumOfSegWC)  );
  for(Int_t i=1; i<=NumOfSegWC; ++i){
    TString title11 = Form("WCSUM-%d Time", i);
    TString title14 = Form("WCSUM-%d NPE", i);
    TString title16 = Form("WCSUM-%d NPE", i);
    TString title17 = Form("WCSUM-%d Up CTime", i);
    TString title21 = Form("WCSUM-%d dE (w/ WC-HT)", i);
    HB1(WCSUMHid +100*i +11, title11, 500, -5., 45.);
    HB1(WCSUMHid +100*i +14, title14, 110, -10., 100.);
    HB1(WCSUMHid +100*i +16, title16, 110, -10., 100.);
    HB1(WCSUMHid +100*i +17, title17, 500, -5., 45.);
    HB1(WCSUMHid +100*i +21, title21, 200, -0.5, 4.5);
  }
  HB2(WCSUMHid +21, "WcSumHitPat%WcSumHitPat[HodoGood]",
      NumOfSegWC, 0., Double_t(NumOfSegWC),
      NumOfSegWC, 0., Double_t(NumOfSegWC));
  HB2(WCSUMHid +22, "CMeanTimeWcSum%CMeanTimeWcSum[HodoGood]",
      120, 10., 40., 120, 10., 40.);
  HB1(WCSUMHid +23, "TDiff WcSum[HodoGood]", 200, -10., 10.);
  HB2(WCSUMHid +24, "WcSumHitPat%WcSumHitPat[HodoGood2]",
      NumOfSegWC, 0., Double_t(NumOfSegWC),
      NumOfSegWC, 0., Double_t(NumOfSegWC));
  HB1(WCSUMHid +30, "#Clusters WcSum", NumOfSegWC+1, 0., Double_t(NumOfSegWC+1));
  HB1(WCSUMHid +31, "ClusterSize WcSum", 5, 0., 5.);
  HB1(WCSUMHid +32, "HitPat Cluster WcSum", 2*NumOfSegWC, 0., Double_t(NumOfSegWC));
  HB1(WCSUMHid +33, "CMeamTime Cluster WcSum", 500, -5., 45.);
  HB1(WCSUMHid +34, "DeltaE Cluster WcSum", 100, -0.5, 4.5);

  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("spill",     &event.spill,     "spill/I");
  //Trig
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  //BH1
  tree->Branch("bh1nhits",   &event.bh1nhits,    "bh1nhits/I");
  tree->Branch("bh1hitpat",   event.bh1hitpat,   Form("bh1hitpat[%d]/I",NumOfSegBH1));
  tree->Branch("bh1ua",       event.bh1ua,       Form("bh1ua[%d]/D", NumOfSegBH1));
  tree->Branch("bh1ut",       event.bh1ut,       Form("bh1ut[%d][%d]/D", NumOfSegBH1, MaxDepth));
  tree->Branch("bh1da",       event.bh1da,       Form("bh1da[%d]/D", NumOfSegBH1));
  tree->Branch("bh1dt",       event.bh1dt,       Form("bh1dt[%d][%d]/D", NumOfSegBH1, MaxDepth));
  //BH2
  tree->Branch("bh2nhits",   &event.bh2nhits,    "bh2nhits/I");
  tree->Branch("bh2hitpat",   event.bh2hitpat,   Form("bh2hitpat[%d]/I", NumOfSegBH2));
  tree->Branch("bh2ua",       event.bh2ua,       Form("bh2ua[%d]/D", NumOfSegBH2));
  tree->Branch("bh2ut",       event.bh2ut,       Form("bh2ut[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2da",       event.bh2da,       Form("bh2da[%d]/D", NumOfSegBH2));
  tree->Branch("bh2dt",       event.bh2dt,       Form("bh2dt[%d][%d]/D", NumOfSegBH2, MaxDepth));
  //BAC
  tree->Branch("bacnhits",   &event.bacnhits,   "bacnhits/I");
  tree->Branch("bachitpat",   event.bachitpat,  Form("bachitpat[%d]/I", NumOfSegBAC));
  tree->Branch("baca",        event.baca,       Form("baca[%d]/D", NumOfSegBAC));
  tree->Branch("bact",        event.bact,       Form("bact[%d][%d]/D", NumOfSegBAC, MaxDepth));
  tree->Branch("bac1nhits",   &event.bac1nhits,   "bac1nhits/I");
  tree->Branch("bac2nhits",   &event.bac2nhits,   "bac2nhits/I");
  //TOF
  tree->Branch("tofnhits",   &event.tofnhits,   "tofnhits/I");
  tree->Branch("tofhitpat",   event.tofhitpat,  Form("tofhitpat[%d]/I", NumOfSegTOF));
  tree->Branch("tofua",       event.tofua,      Form("tofua[%d]/D", NumOfSegTOF));
  tree->Branch("tofut",       event.tofut,      Form("tofut[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofda",       event.tofda,      Form("tofda[%d]/D", NumOfSegTOF));
  tree->Branch("tofdt",       event.tofdt,      Form("tofdt[%d][%d]/D", NumOfSegTOF, MaxDepth));
  //AC1
  tree->Branch("ac1nhits", &event.ac1nhits, "ac1nhits/I");
  tree->Branch("ac1hitpat", event.ac1hitpat, Form("ac1hitpat[%d]/I", NumOfSegAC1));
  tree->Branch("ac1a", event.ac1a, Form("ac1a[%d]/D", NumOfSegAC1));
  tree->Branch("ac1t", event.ac1t, Form("ac1t[%d][%d]/D", NumOfSegAC1, MaxDepth));
  //WC
  tree->Branch("wcnhits",   &event.wcnhits,    "wcnhits/I");
  tree->Branch("wchitpat",   event.wchitpat,   Form("wchitpat[%d]/I", NumOfSegWC));
  tree->Branch("wcua",       event.wcua,       Form("wcua[%d]/D", NumOfSegWC));
  tree->Branch("wcut",       event.wcut,       Form("wcut[%d][%d]/D", NumOfSegWC, MaxDepth));
  tree->Branch("wcda",       event.wcda,       Form("wcda[%d]/D", NumOfSegWC));
  tree->Branch("wcdt",       event.wcdt,       Form("wcdt[%d][%d]/D", NumOfSegWC, MaxDepth));
  //WCSUM
  tree->Branch("wcsumnhits", &event.wcsumnhits, "wcsumnhits/I");
  tree->Branch("wcsumhitpat", &event.wcsumhitpat, Form("wcsumhitpat[%d]/I", NumOfSegWC));
  tree->Branch("wcsuma", &event.wcsuma, Form("wcsuma[%d]/D", NumOfSegWC));
  tree->Branch("wcsumt", &event.wcsumt, Form("wcsumt[%d][%d]/D", NumOfSegWC, MaxDepth));

  //Normalized data
  tree->Branch("bh1mt",     event.bh1mt,     Form("bh1mt[%d][%d]/D", NumOfSegBH1, MaxDepth));
  tree->Branch("bh1utime",     event.bh1utime,     Form("bh1utime[%d][%d]/D", NumOfSegBH1, MaxDepth));
  tree->Branch("bh1dtime",     event.bh1dtime,     Form("bh1dtime[%d][%d]/D", NumOfSegBH1, MaxDepth));
  tree->Branch("bh1uctime",     event.bh1uctime,     Form("bh1uctime[%d][%d]/D", NumOfSegBH1, MaxDepth));
  tree->Branch("bh1dctime",     event.bh1dctime,     Form("bh1dctime[%d][%d]/D", NumOfSegBH1, MaxDepth));
  tree->Branch("bh1de",     event.bh1de,     Form("bh1de[%d]/D", NumOfSegBH1));
  tree->Branch("bh1ude",     event.bh1ude,     Form("bh1ude[%d]/D", NumOfSegBH1));
  tree->Branch("bh1dde",     event.bh1dde,     Form("bh1dde[%d]/D", NumOfSegBH1));
  tree->Branch("bh2mt",     event.bh2mt,     Form("bh2mt[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2utime",     event.bh2utime,     Form("bh2utime[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2dtime",     event.bh2dtime,     Form("bh2dtime[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2uctime",     event.bh2uctime,     Form("bh2uctime[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2dctime",     event.bh2dctime,     Form("bh2dctime[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("bh2de",     event.bh2de,     Form("bh2de[%d]/D", NumOfSegBH2));
  tree->Branch("bh2ude",     event.bh2ude,     Form("bh2ude[%d]/D", NumOfSegBH2));
  tree->Branch("bh2dde",     event.bh2dde,     Form("bh2dde[%d]/D", NumOfSegBH2));
#if HodoHitPos
  tree->Branch("bh2hitpos",     event.bh2hitpos,     Form("bh2hitpos[%d][%d]/D", NumOfSegBH2, MaxDepth));
#endif
  tree->Branch("bacmt",     event.bacmt,     Form("bacmt[%d][%d]/D", NumOfSegBAC, MaxDepth));
  tree->Branch("bacde",     event.bacde,     Form("bacde[%d]/D", NumOfSegBAC));
  tree->Branch("tofmt",     event.tofmt,     Form("tofmt[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofde",     event.tofde,     Form("tofde[%d]/D", NumOfSegTOF));
  tree->Branch("tofude",     event.tofude,     Form("tofude[%d]/D", NumOfSegTOF));
  tree->Branch("tofdde",     event.tofdde,     Form("tofdde[%d]/D", NumOfSegTOF));
  tree->Branch("tofctu",     event.tofctu,     Form("tofctu[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofctd",     event.tofctd,     Form("tofctd[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("tofcmt",     event.tofcmt,     Form("tofcmt[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("wcmt",     event.wcmt,     Form("wcmt[%d][%d]/D", NumOfSegWC, MaxDepth));
  tree->Branch("wcde",     event.wcde,     Form("wcde[%d]/D", NumOfSegWC));

  tree->Branch("t0",        event.t0,        Form("t0[%d][%d]/D",  NumOfSegBH2, MaxDepth));
  tree->Branch("ct0",       event.ct0,       Form("ct0[%d][%d]/D", NumOfSegBH2, MaxDepth));
  tree->Branch("btof",      event.btof,      Form("btof[%d][%d]/D",  NumOfSegBH1, NumOfSegBH2));
  tree->Branch("cbtof",     event.cbtof,     Form("cbtof[%d][%d]/D", NumOfSegBH1, NumOfSegBH2));

  tree->Branch("Time0Seg", &event.Time0Seg,  "Time0Seg/D");
  tree->Branch("deTime0",  &event.deTime0,   "deTime0/D");
  tree->Branch("Time0",    &event.Time0,     "Time0/D");
  tree->Branch("CTime0",   &event.CTime0,    "CTime0/D");

  tree->Branch("Btof0Seg", &event.Btof0Seg,  "Btof0Seg/D");
  tree->Branch("deBtof0",  &event.deBtof0,   "deBtof0/D");
  tree->Branch("Btof0",    &event.Btof0,     "Btof0/D");
  tree->Branch("CBtof0",   &event.CBtof0,    "CBtof0/D");

  //SAC3
  //tree->Branch("sac3nhits", &event.sac3nhits, "sac3nhits/I");
  tree->Branch("sac3a", event.sac3a, Form("sac3a[%d]/D", NumOfSegSAC3));
  tree->Branch("sac3t", event.sac3t, Form("sac3t[%d][%d]/D", NumOfSegSAC3, MaxDepth));
  //SFV
  //tree->Branch("sfvnhits", &event.sfvnhits, "sfvnhits/I");
  //tree->Branch("sfvhitpat", event.sfvhitpat, Form("sfvhitpat[%d]/I", NumOfSegSFV));
  tree->Branch("sfvt", event.sfvt, Form("sfvt[%d][%d]/D", NumOfSegSFV, MaxDepth));

  ////////////////////////////////////////////
  //Dst
  hodo = new TTree("hodo","Data Summary Table of Hodoscope");
  hodo->Branch("evnum",     &dst.evnum,     "evnum/I");
  hodo->Branch("spill",     &dst.spill,     "spill/I");
  hodo->Branch("trigpat",    dst.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  hodo->Branch("trigflag",   dst.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));
  hodo->Branch("nhBh1",     &dst.nhBh1,     "nhBh1/I");
  hodo->Branch("csBh1",      dst.csBh1,     "csBh1[nhBh1]/I");
  hodo->Branch("Bh1Seg",     dst.Bh1Seg,    "Bh1Seg[nhBh1]/D");
  hodo->Branch("tBh1",       dst.tBh1,      "tBh1[nhBh1]/D");
  hodo->Branch("dtBh1",      dst.dtBh1,     "dtBh1[nhBh1]/D");
  hodo->Branch("deBh1",      dst.deBh1,     "deBh1[nhBh1]/D");

  hodo->Branch("nhBh2",     &dst.nhBh2,     "nhBh2/I");
  hodo->Branch("csBh2",      dst.csBh2,     "csBh2[nhBh2]/I");
  hodo->Branch("Bh2Seg",     dst.Bh2Seg,    "Bh2Seg[nhBh2]/D");
  hodo->Branch("tBh2",       dst.tBh2,      "tBh2[nhBh2]/D");
  hodo->Branch("t0Bh2",      dst.t0Bh2,     "t0Bh2[nhBh2]/D");
  hodo->Branch("dtBh2",      dst.dtBh2,     "dtBh2[nhBh2]/D");
  hodo->Branch("deBh2",      dst.deBh2,     "deBh2[nhBh2]/D");
#if HodoHitPos
  hodo->Branch("posBh2",     dst.posBh2,    "posBh2[nhBh2]/D");
#endif
  hodo->Branch("btof",       dst.btof,      "btof[nhBh1]/D");
  hodo->Branch("cbtof",      dst.cbtof,     "cbtof[nhBh1]/D");

  hodo->Branch("Btof0Seg",  &dst.Btof0Seg,  "Btof0Seg/D");
  hodo->Branch("deBtof0",   &dst.deBtof0,   "deBtof0/D");
  hodo->Branch("Btof0",     &dst.Btof0,     "Btof0/D");
  hodo->Branch("CBtof0",    &dst.CBtof0,    "CBtof0/D");

  hodo->Branch("Time0Seg",  &dst.Time0Seg,  "Time0Seg/D");
  hodo->Branch("deTime0",   &dst.deTime0,   "deTime0/D");
  hodo->Branch("Time0",     &dst.Time0,     "Time0/D");
  hodo->Branch("CTime0",    &dst.CTime0,    "CTime0/D");

  hodo->Branch("nhBac",     &dst.nhBac,     "nhBac/I");
  hodo->Branch("BacSeg",     dst.BacSeg,    "BacSeg[nhBac]/D");
  hodo->Branch("tBac",       dst.tBac,      "tBac[nhBac]/D");
  hodo->Branch("deBac",      dst.deBac,     "deBac[nhBac]/D");

  hodo->Branch("nhTof",     &dst.nhTof,     "nhTof/I");
  hodo->Branch("csTof",      dst.csTof,     "csTof[nhTof]/I");
  hodo->Branch("TofSeg",     dst.TofSeg,    "TofSeg[nhTof]/D");
  hodo->Branch("tTof",       dst.tTof,      "tTof[nhTof]/D");
  hodo->Branch("dtTof",      dst.dtTof,     "dtTof[nhTof]/D");
  hodo->Branch("deTof",      dst.deTof,     "deTof[nhTof]/D");

  hodo->Branch("tofua", event.tofua, Form("tofua[%d]/D", NumOfSegTOF));
  hodo->Branch("tofut", event.tofut, Form("tofut[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("tofda", event.tofda, Form("tofda[%d]/D", NumOfSegTOF));
  hodo->Branch("tofdt", event.tofdt, Form("tofdt[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("utTofSeg", dst.utTofSeg,
               Form("utTofSeg[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("dtTofSeg", dst.dtTofSeg,
               Form("dtTofSeg[%d][%d]/D", NumOfSegTOF, MaxDepth));
  hodo->Branch("udeTofSeg",  dst.udeTofSeg,
               Form("udeTofSeg[%d]/D", NumOfSegTOF));
  hodo->Branch("ddeTofSeg",  dst.ddeTofSeg,
               Form("ddeTofSeg[%d]/D", NumOfSegTOF));

  hodo->Branch("nhAc1",     &dst.nhAc1,     "nhAc1/I");
  hodo->Branch("csAc1",      dst.csAc1,     "csAc1[nhAc1]/I");
  hodo->Branch("Ac1Seg",     dst.Ac1Seg,    "Ac1Seg[nhAc1]/D");
  hodo->Branch("tAc1",       dst.tAc1,      "tAc1[nhAc1]/D");

  hodo->Branch("nhWc",     &dst.nhWc,     "nhWc/I");
  hodo->Branch("csWc",      dst.csWc,     "csWc[nhWc]/I");
  hodo->Branch("WcSeg",     dst.WcSeg,    "WcSeg[nhWc]/D");
  hodo->Branch("tWc",       dst.tWc,      "tWc[nhWc]/D");
  hodo->Branch("deWc",      dst.deWc,     "deWc[nhWc]/D");

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
     InitializeParameter<HodoPHCMan>("HDPHC")   &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
