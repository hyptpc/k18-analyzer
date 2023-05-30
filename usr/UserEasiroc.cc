// -*- C++ -*-

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "RootHelper.hh"
#include "HodoHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
#include "UnpackerManager.hh"

#define HodoCut    0 // with BH1/BH2
#define TIME_CUT   0 // in cluster analysis
#define FHitBranch 1 // make FiberHit branches (becomes heavy)

namespace
{
enum EUorD { kU, kD, kUorD };
using namespace root;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];
  // BFT
  Int_t    bft_nhits;
  Int_t    bft_unhits;
  Int_t    bft_dnhits;
  Int_t    bft_uhitpat[NumOfSegBFT];
  Int_t    bft_dhitpat[NumOfSegBFT];
  Double_t bft_utdc[NumOfSegBFT][MaxDepth];
  Double_t bft_dtdc[NumOfSegBFT][MaxDepth];
  Double_t bft_utrailing[NumOfSegBFT][MaxDepth];
  Double_t bft_dtrailing[NumOfSegBFT][MaxDepth];
  Double_t bft_utot[NumOfSegBFT][MaxDepth];
  Double_t bft_dtot[NumOfSegBFT][MaxDepth];
  Int_t    bft_udepth[NumOfSegBFT];
  Int_t    bft_ddepth[NumOfSegBFT];
  Int_t    bft_ncl;
  Int_t    bft_clsize[NumOfSegBFT];
  Double_t bft_ctime[NumOfSegBFT];
  Double_t bft_ctot[NumOfSegBFT];
  Double_t bft_clpos[NumOfSegBFT];
  Double_t bft_clseg[NumOfSegBFT];

  // AFT raw
  Int_t    aft_nhits[NumOfPlaneAFT];
  Int_t    aft_hitpat[NumOfPlaneAFT][NumOfSegAFT];
  Double_t aft_tdc[NumOfPlaneAFT][NumOfSegAFT][kUorD][MaxDepth];
  Double_t aft_adc_high[NumOfPlaneAFT][NumOfSegAFT][kUorD];
  Double_t aft_adc_low[NumOfPlaneAFT][NumOfSegAFT][kUorD];
  // AFT normalized
  Double_t aft_mt[NumOfPlaneAFT][NumOfSegAFT][MaxDepth];
  Double_t aft_cmt[NumOfPlaneAFT][NumOfSegAFT][MaxDepth];
  Double_t aft_mtot[NumOfPlaneAFT][NumOfSegAFT][MaxDepth];
  Double_t aft_de_high[NumOfPlaneAFT][NumOfSegAFT];
  Double_t aft_de_low[NumOfPlaneAFT][NumOfSegAFT];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum      = 0;
  bft_nhits  = 0;
  bft_unhits = 0;
  bft_dnhits = 0;
  bft_ncl    = 0;

  for(Int_t it=0; it<NumOfSegTrig; it++){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<NumOfSegBFT; it++){
    bft_uhitpat[it] = -1;
    bft_dhitpat[it] = -1;
    bft_udepth[it]  = 0;
    bft_ddepth[it]  = 0;
    for(Int_t that=0; that<MaxDepth; that++){
      bft_utdc[it][that] = qnan;
      bft_dtdc[it][that] = qnan;
      bft_utrailing[it][that] = qnan;
      bft_dtrailing[it][that] = qnan;
      bft_utot[it][that] = qnan;
      bft_dtot[it][that] = qnan;
    }
    bft_clsize[it] = 0;
    bft_ctime[it]  = qnan;
    bft_ctot[it]   = qnan;
    bft_clpos[it]  = qnan;
    bft_clseg[it]  = qnan;
  }

  for(Int_t p=0; p<NumOfPlaneAFT; p++){
    aft_nhits[p] = 0;
    for(Int_t seg=0; seg<NumOfSegAFT; seg++){
      aft_hitpat[p][seg] = -1;
      for(Int_t ud=0; ud<kUorD; ud++){
        aft_adc_high[p][seg][ud] = qnan;
        aft_adc_low[p][seg][ud] = qnan;
        for(Int_t i=0; i<MaxDepth; i++){
          aft_tdc[p][seg][ud][i] = qnan;
        }
      }
      aft_de_high[p][seg] = qnan;
      aft_de_low[p][seg] = qnan;
      for(Int_t i=0; i<MaxDepth; i++){
        aft_mt[p][seg][i] = qnan;
        aft_cmt[p][seg][i] = qnan;
        aft_mtot[p][seg][i] = qnan;
      }
    }
  }
}

//_____________________________________________________________________________
namespace root
{
Event  event;
TH1   *h[MaxHist];
TTree *tree;
enum eDetHid
{
  BFTHid  =  10000,
  AFTHid  = 100000,
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
#if HodoCut
  static const Double_t MinDeBH2   = gUser.GetParameter("DeBH2", 0);
  static const Double_t MaxDeBH2   = gUser.GetParameter("DeBH2", 1);
  static const Double_t MinDeBH1   = gUser.GetParameter("DeBH1", 0);
  static const Double_t MaxDeBH1   = gUser.GetParameter("DeBH1", 1);
  static const Double_t MinBeamToF = gUser.GetParameter("BTOF",  0);
  static const Double_t MaxBeamToF = gUser.GetParameter("BTOF",  1);
#endif
  static const Double_t MinTdcBFT  = gUser.GetParameter("TdcBFT",  0);
  static const Double_t MaxTdcBFT  = gUser.GetParameter("TdcBFT",  1);
#if TIME_CUT
  static const Double_t MinTimeBFT = gUser.GetParameter("TimeBFT", 0);
  static const Double_t MaxTimeBFT = gUser.GetParameter("TimeBFT", 1);
#endif

  RawData rawData;
  rawData.DecodeHits("TFlag");
  rawData.DecodeHits("BH1");
  rawData.DecodeHits("BH2");
  rawData.DecodeHits("BFT");
  rawData.DecodeHits("AFT");
  HodoAnalyzer hodoAna(rawData);

  event.evnum = gUnpacker.get_event_number();

  HF1(1, 0);

  ///// Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  {
    for(const auto& hit: rawData.GetHodoRawHitContainer("TFlag")){
      Int_t seg = hit->SegmentId();
      Int_t tdc = hit->GetTdc();
      if(tdc > 0){
	event.trigpat[trigger_flag.count()] = seg;
	event.trigflag[seg] = tdc;
        trigger_flag.set(seg);
	HF1(10, seg);
	HF1(10+seg, tdc);
      }
    }
  }

  if(trigger_flag[trigger::kSpillEnd])
    return true;

  HF1(1, 1);

  ////////// BH2 time 0
  hodoAna.DecodeHits<BH2Hit>("BH2");
  Int_t nhBh2 = hodoAna.GetNHits("BH2");
#if HodoCut
  if(nhBh2==0) return true;
#endif
  HF1(1, 2);
  Double_t time0 = qnan;
  ////////// BH2 Analysis
  for(Int_t i=0; i<nhBh2; ++i){
    const auto& hit = hodoAna.GetHit("BH2", i);
    if(!hit) continue;
    Double_t cmt = hit->CMeanTime();
    Double_t ct0 = hit->CTime0();
    Double_t min_time = qnan;
#if HodoCut
    Double_t dE  = hit->DeltaE();
    if(dE<MinDeBH2 || MaxDeBH2<dE)
      continue;
#endif
    if(std::abs(cmt)<std::abs(min_time)){
      min_time = cmt;
      time0    = ct0;
    }
  }

  HF1(1, 3);

  ////////// BH1 Analysis
  hodoAna.DecodeHits("BH1");
  Int_t nhBh1 = hodoAna.GetNHits("BH1");
#if HodoCut
  if(nhBh1==0) return true;
#endif
  HF1(1, 4);
  Double_t btof0 = qnan;
  for(Int_t i=0; i<nhBh1; ++i){
    const auto& hit = hodoAna.GetHit("BH1", i);
    if(!hit) continue;
    Double_t cmt  = hit->CMeanTime();
    Double_t btof = cmt - time0;
#if HodoCut
    Double_t dE   = hit->DeltaE();
    if(dE<MinDeBH1 || MaxDeBH1<dE) continue;
    if(btof<MinBeamToF || MaxBeamToF<btof) continue;
#endif
    if(std::abs(btof)<std::abs(btof0)){
      btof0 = btof;
    }
  }

  HF1(1, 5);

  HF1(1, 6);

  ////////// BFT
  hodoAna.DecodeHits<FiberHit>("BFT");
  {
    const auto& U = HodoRawHit::kUp;
    Int_t nh = hodoAna.GetNHits("BFT");
    Int_t unhits = 0;
    Int_t dnhits = 0;
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit<FiberHit>("BFT", i);
      const auto& rhit = hit->GetRawHit();
      Int_t plane = hit->PlaneId();
      Int_t seg = hit->SegmentId();
      // Raw
      Int_t mh_l = rhit->GetSizeTdcLeading();
      if(plane == 0) event.bft_udepth[seg] = mh_l;
      for(Int_t j=0; j<mh_l; ++j){
        Double_t leading = rhit->GetTdcLeading(U, j);
        HF1(BFTHid +plane + 6, leading);
        HF2(BFTHid +plane + 10, seg, leading);
        HF1(BFTHid +1000*(1+plane)+seg+1, leading);
        if(plane == 0) event.bft_utdc[seg][j] = leading;
        if(plane == 1) event.bft_dtdc[seg][j] = leading;
        if(MinTdcBFT<leading && leading<MaxTdcBFT){
	  HF1(BFTHid +plane+4, seg+0.5);
	  if(plane==0) event.bft_uhitpat[unhits++] = seg;
	  if(plane==1) event.bft_dhitpat[dnhits++] = seg;
        }
      }
      Int_t mh_t = rhit->GetSizeTdcTrailing();
      if(plane == 1) event.bft_ddepth[seg] = mh_l;
      for(Int_t j=0; j<mh_t; ++j){
        Double_t trailing = rhit->GetTdcTrailing(U, j);
        if(plane == 0) event.bft_utrailing[seg][j] = trailing;
        if(plane == 1) event.bft_dtrailing[seg][j] = trailing;
      }
      // Normalized
      Int_t mh = hit->GetEntries();
      for(Int_t j=0; j<mh; ++j){
        Double_t time  = hit->MeanTime(j);
        Double_t ctime = hit->CMeanTime(j);
        Double_t tot   = hit->TOT(U, j);
        HF1(BFTHid +plane+8, tot);
        HF2(BFTHid +plane+12, seg, tot);
        HF1(BFTHid +plane+21, time);
        HF1(BFTHid +plane+31, ctime);
        HF1(BFTHid +1000*(plane+3)+seg+1, tot);
        if(plane == 0) event.bft_utot[seg][j] = tot;
        if(plane == 1) event.bft_dtot[seg][j] = tot;
        if(-10.<time && time<10.){
          HF2(BFTHid +plane+23, tot, time);
          HF2(BFTHid +plane+33, tot, ctime);
          HF2(BFTHid +1000*(plane+5)+seg+1, tot, time);
          HF2(BFTHid +1000*(plane+7)+seg+1, tot, ctime);
        }
      }
    }
    HF1(BFTHid +1, unhits);
    HF1(BFTHid +2, dnhits);
    HF1(BFTHid +3, unhits + dnhits);
    event.bft_unhits = unhits;
    event.bft_dnhits = dnhits;
    event.bft_nhits  = unhits + dnhits;

    // Fiber Cluster
#if TIME_CUT
    hodoAna.TimeCut("BFT", MinTimeBFT, MaxTimeBFT);
#endif
    Int_t ncl = hodoAna.GetNClusters("BFT");
    if(ncl > NumOfSegBFT){
      // std::cout << "#W BFT too much number of clusters" << std::endl;
      ncl = NumOfSegBFT;
    }
    event.bft_ncl = ncl;
    HF1(BFTHid +101, ncl);
    for(Int_t i=0; i<ncl; ++i){
      const auto& cl = hodoAna.GetCluster("BFT", i);
      if(!cl) continue;
      // cl->Print();
      Double_t clsize = cl->ClusterSize();
      Double_t ctime  = cl->CMeanTime();
      Double_t ctot   = cl->TOT();
      Double_t pos    = cl->MeanPosition();
      Double_t seg    = cl->MeanSeg();
      event.bft_clsize[i] = clsize;
      event.bft_ctime[i]  = ctime;
      event.bft_ctot[i]   = ctot;
      event.bft_clpos[i]  = pos;
      event.bft_clseg[i]  = seg;
      HF1(BFTHid +102, clsize);
      HF1(BFTHid +103, ctime);
      HF1(BFTHid +104, ctot);
      HF2(BFTHid +105, ctot, ctime);
      HF1(BFTHid +106, pos);
    }
  }

  ////////// AFT
  for(const auto& hit: rawData.GetHodoRawHitContainer("AFT")){
    // hit->Print();
    Int_t plane = hit->PlaneId();
    Int_t seg = hit->SegmentId();
    for(Int_t ud=0; ud<kUorD; ++ud){
      auto adc_high = hit->GetAdcHigh(ud);
      auto adc_low = hit->GetAdcLow(ud);
      event.aft_adc_high[plane][seg][ud] = adc_high;
      event.aft_adc_low[plane][seg][ud] = adc_low;
      HF1(AFTHid+plane*1000+7+ud, adc_high);
      HF1(AFTHid+plane*1000+9+ud, adc_low);
      HF2(AFTHid+plane*1000+15+ud, seg, adc_high);
      HF2(AFTHid+plane*1000+17+ud, seg, adc_low);
      for(Int_t i=0, n=hit->GetSizeTdcLeading(ud); i<n; ++i){
        auto tdc = hit->GetTdc(ud, i);
        event.aft_tdc[plane][seg][ud][i] = tdc;
        HF1(AFTHid+plane*1000+3+ud, tdc);
        HF2(AFTHid+plane*1000+11+ud, seg, tdc);
        // HF1(AFTHid+plane*1000+seg+100+ud*100, tdc);
      }
    }
  }

  hodoAna.DecodeHits<FiberHit>("AFT");
  for(Int_t i=0, n=hodoAna.GetNHits("AFT"); i<n; ++i){
    const auto& hit = hodoAna.GetHit<FiberHit>("AFT", i);
    // hit->Print();
    // const auto& rhit = hit->GetRawHit();
    // rhit->Print();
    Int_t plane = hit->PlaneId();
    Int_t seg = hit->SegmentId();
    event.aft_hitpat[plane][event.aft_nhits[plane]++] = seg;
    HF1(AFTHid+plane*1000+2, seg);
    Int_t m = hit->GetEntries();
    for(Int_t j=0; j<m; ++j){
      auto mt = hit->MeanTime(j);
      auto cmt = hit->CMeanTime(j);
      auto mtot = hit->MeanTOT(j);
      event.aft_mt[plane][seg][j] = mt;
      event.aft_cmt[plane][seg][j] = cmt;
      event.aft_mtot[plane][seg][j] = mtot;
      HF1(AFTHid+plane*1000+21, mt);
      HF1(AFTHid+plane*1000+22, cmt);
      HF1(AFTHid+plane*1000+23, mtot);
      HF2(AFTHid+plane*1000+31, seg, mt);
      HF2(AFTHid+plane*1000+32, seg, cmt);
      HF2(AFTHid+plane*1000+33, seg, mtot);
      for(Int_t ud=0; ud<kUorD; ++ud){
        auto tot = hit->TOT(ud, j);
        HF1(AFTHid+plane*1000+5+ud, tot);
        HF2(AFTHid+plane*1000+13+ud, seg, tot);
        // HF1(AFTHid+plane*1000+seg+300+ud*100, tot);
      }
    }
    auto de_high = hit->DeltaEHighGain();
    auto de_low = hit->DeltaELowGain();
    event.aft_de_high[plane][seg] = de_high;
    event.aft_de_low[plane][seg] = de_low;
    HF1(AFTHid+plane*1000+24, de_high);
    HF1(AFTHid+plane*1000+25, de_low);
    HF2(AFTHid+plane*1000+34, seg, de_high);
    HF2(AFTHid+plane*1000+35, seg, de_low);
  }
  for(Int_t plane=0; plane<NumOfPlaneAFT; ++plane){
    HF1(AFTHid+plane*1000+1, event.aft_nhits[plane]);
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
Bool_t
ConfMan::InitializeHistograms()
{
  const Int_t    NbinAdc = 4000;
  const Double_t MinAdc  =    0.;
  const Double_t MaxAdc  = 4000.;

  const Int_t    NbinTdc = 1000;
  const Double_t MinTdc  =    0.;
  const Double_t MaxTdc  = 1000.;

  const Int_t    NbinTot =  170;
  const Double_t MinTot  =  -10.;
  const Double_t MaxTot  =  160.;

  const Int_t    NbinTime = 100;
  const Double_t MinTime  = -50.;
  const Double_t MaxTime  =  50.;

  const Int_t    NbinDe = 1000;
  const Double_t MinDe  =  0.;
  const Double_t MaxDe  = 10.;

  HB1( 1, "Status",  20,   0., 20.);
  HB1(10, "Trigger HitPat", NumOfSegTrig, 0., Double_t(NumOfSegTrig));
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1(10+i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000);
  }

  //BFT
  HB1(BFTHid + 1, "BFT Nhits U",   NumOfSegBFT, 0., (Double_t)NumOfSegBFT);
  HB1(BFTHid + 2, "BFT Nhits D",   NumOfSegBFT, 0., (Double_t)NumOfSegBFT);
  HB1(BFTHid + 3, "BFT Nhits",     NumOfSegBFT, 0., (Double_t)NumOfSegBFT);
  HB1(BFTHid + 4, "BFT Hitpat U",  NumOfSegBFT, 0., (Double_t)NumOfSegBFT);
  HB1(BFTHid + 5, "BFT Hitpat D",  NumOfSegBFT, 0., (Double_t)NumOfSegBFT);
  HB1(BFTHid + 6, "BFT Tdc U",      NbinTdc, MinTdc, MaxTdc);
  HB1(BFTHid + 7, "BFT Tdc D",      NbinTdc, MinTdc, MaxTdc);
  HB1(BFTHid + 8, "BFT Tot U",      NbinTot, MinTot, MaxTot);
  HB1(BFTHid + 9, "BFT Tot D",      NbinTot, MinTot, MaxTot);
  HB2(BFTHid +10, "BFT Tdc U%Seg",
      NumOfSegBFT, 0., (Double_t)NumOfSegBFT, NbinTdc, MinTdc, MaxTdc);
  HB2(BFTHid +11, "BFT Tdc D%Seg",
      NumOfSegBFT, 0., (Double_t)NumOfSegBFT, NbinTdc, MinTdc, MaxTdc);
  HB2(BFTHid +12, "BFT Tot U%Seg",
      NumOfSegBFT, 0., (Double_t)NumOfSegBFT, NbinTot, MinTot, MaxTot);
  HB2(BFTHid +13, "BFT Tot D%Seg",
      NumOfSegBFT, 0., (Double_t)NumOfSegBFT, NbinTot, MinTot, MaxTot);
  for(Int_t i=0; i<NumOfSegBFT; i++){
    HB1(BFTHid +1000+i+1, Form("BFT Tdc U-%d", i+1), NbinTdc, MinTdc, MaxTdc);
    HB1(BFTHid +2000+i+1, Form("BFT Tdc D-%d", i+1), NbinTdc, MinTdc, MaxTdc);
    HB1(BFTHid +3000+i+1, Form("BFT Tot U-%d", i+1), NbinTot, MinTot, MaxTot);
    HB1(BFTHid +4000+i+1, Form("BFT Tot D-%d", i+1), NbinTot, MinTot, MaxTot);
    HB2(BFTHid +5000+i+1, Form("BFT Time%%Tot U-%d", i+1),
        NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
    HB2(BFTHid +6000+i+1, Form("BFT Time%%Tot D-%d", i+1),
        NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
    HB2(BFTHid +7000+i+1, Form("BFT CTime%%Tot U-%d", i+1),
        NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
    HB2(BFTHid +8000+i+1, Form("BFT CTime%%Tot D-%d", i+1),
        NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  }
  HB1(BFTHid +21, "BFT Time U",     NbinTime, MinTime, MaxTime);
  HB1(BFTHid +22, "BFT Time D",     NbinTime, MinTime, MaxTime);
  HB2(BFTHid +23, "BFT Time%Tot U", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB2(BFTHid +24, "BFT Time%Tot D", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB1(BFTHid +31, "BFT CTime U",     NbinTime, MinTime, MaxTime);
  HB1(BFTHid +32, "BFT CTime D",     NbinTime, MinTime, MaxTime);
  HB2(BFTHid +33, "BFT CTime%Tot U", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB2(BFTHid +34, "BFT CTime%Tot D", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);

  HB1(BFTHid +101, "BFT NCluster", 100, 0, 100);
  HB1(BFTHid +102, "BFT Cluster Size", 5, 0, 5);
  HB1(BFTHid +103, "BFT CTime (Cluster)", 100., -20., 30.);
  HB1(BFTHid +104, "BFT Tot (Cluster)", NbinTot, MinTot, MaxTot);
  HB2(BFTHid +105, "BFT CTime%Tot (Cluster)",
      NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB1(BFTHid +106, "BFT Cluster Position",
      NumOfSegBFT, -0.5*(Double_t)NumOfSegBFT, 0.5*(Double_t)NumOfSegBFT);

  //AFT
  for(Int_t plane=0; plane<NumOfPlaneAFT; ++plane){
    HB1(AFTHid+plane*1000+1, Form("AFT Nhits Plane#%d", plane), NumOfSegAFT, 0., NumOfSegAFT);
    HB1(AFTHid+plane*1000+2, Form("AFT Hitpat Plane#%d", plane), NumOfSegAFT, 0., NumOfSegAFT);
    HB1(AFTHid+plane*1000+21, Form("AFT MeanTime Plane#%d", plane), NbinTime, MinTime, MaxTime);
    HB1(AFTHid+plane*1000+22, Form("AFT CMeanTime Plane#%d", plane), NbinTime, MinTime, MaxTime);
    HB1(AFTHid+plane*1000+23, Form("AFT MeanTot Plane#%d", plane), NbinTot, MinTot, MaxTot);
    HB1(AFTHid+plane*1000+24, Form("AFT DeltaE HighGain Plane#%d", plane), NbinDe, MinDe, MaxDe);
    HB1(AFTHid+plane*1000+25, Form("AFT DeltaE LowGain Plane#%d", plane), NbinDe, MinDe, MaxDe);
    HB2(AFTHid+plane*1000+31, Form("AFT MeanTime Plane#%d", plane),
        NumOfSegAFT, 0., NumOfSegAFT, NbinTime, MinTime, MaxTime);
    HB2(AFTHid+plane*1000+32, Form("AFT CMeanTime Plane#%d", plane),
        NumOfSegAFT, 0., NumOfSegAFT, NbinTime, MinTime, MaxTime);
    HB2(AFTHid+plane*1000+33, Form("AFT MeanTot Plane#%d", plane),
        NumOfSegAFT, 0., NumOfSegAFT, NbinTot, MinTot, MaxTot);
    HB2(AFTHid+plane*1000+34, Form("AFT DeltaE HighGain Plane#%d", plane),
        NumOfSegAFT, 0., NumOfSegAFT, NbinDe, MinDe, MaxDe);
    HB2(AFTHid+plane*1000+35, Form("AFT DeltaE LowGain Plane#%d", plane),
        NumOfSegAFT, 0., NumOfSegAFT, NbinDe, MinDe, MaxDe);
    for(Int_t ud=0; ud<kUorD; ++ud){
      const Char_t* s = (ud == kU) ? "U" : "D";
      HB1(AFTHid+plane*1000+3+ud, Form("AFT Tdc %s Plane#%d", s, plane), NbinTdc, MinTdc, MaxTdc);
      HB1(AFTHid+plane*1000+5+ud, Form("AFT Tot %s Plane#%d", s, plane), NbinTdc, MinTdc, MaxTdc);
      HB1(AFTHid+plane*1000+7+ud, Form("AFT AdcHigh %s Plane#%d", s, plane), NbinAdc, MinAdc, MaxAdc);
      HB1(AFTHid+plane*1000+9+ud, Form("AFT AdcHigh %s Plane#%d", s, plane), NbinAdc, MinAdc, MaxAdc);
      HB2(AFTHid+plane*1000+11+ud, Form("AFT Tdc %s%%Seg Plane#%d", s, plane),
          NumOfSegAFT, 0., NumOfSegAFT, NbinTdc, MinTdc, MaxTdc);
      HB2(AFTHid+plane*1000+13+ud, Form("AFT Tot %s%%Seg Plane#%d", s, plane),
          NumOfSegAFT, 0., NumOfSegAFT, NbinTot, MinTot, MaxTot);
      HB2(AFTHid+plane*1000+15+ud, Form("AFT AdcHigh %s%%Seg Plane#%d", s, plane),
          NumOfSegAFT, 0., NumOfSegAFT, NbinAdc, MinAdc, MaxAdc);
      HB2(AFTHid+plane*1000+17+ud, Form("AFT AdcLow %s%%Seg Plane#%d", s, plane),
          NumOfSegAFT, 0., NumOfSegAFT, NbinAdc, MinAdc, MaxAdc);
    }
  }

  // HB1(AFTHid +101, "AFT NCluster", 100, 0, 100);
  // HB1(AFTHid +102, "AFT Cluster Size", 5, 0, 5);
  // HB1(AFTHid +103, "AFT CTime (Cluster)", 100., -20., 30.);
  // HB1(AFTHid +104, "AFT Tot (Cluster)", NbinTot, MinTot, MaxTot);
  // HB2(AFTHid +105, "AFT CTime%Tot (Cluster)",
  //     NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  // HB1(AFTHid +106, "AFT Cluster Position",
  //     NumOfSegAFT, -0.5*(Double_t)NumOfSegAFT, 0.5*(Double_t)NumOfSegAFT);


  //Tree
  HBTree("ea0c", "tree of Easiroc");
  //Trig
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  //BFT
#if FHitBranch
  tree->Branch("bft_nhits",     &event.bft_nhits,        "bft_nhits/I");
  tree->Branch("bft_unhits",    &event.bft_unhits,       "bft_unhits/I");
  tree->Branch("bft_dnhits",    &event.bft_dnhits,       "bft_dnhits/I");
  tree->Branch("bft_uhitpat",    event.bft_uhitpat,      "bft_uhitpat[bft_unhits]/I");
  tree->Branch("bft_dhitpat",    event.bft_dhitpat,      "bft_dhitpat[bft_dnhits]/I");
  tree->Branch("bft_utdc",       event.bft_utdc,         Form("bft_utdc[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_dtdc",       event.bft_dtdc,         Form("bft_dtdc[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_utrailing",  event.bft_utrailing,    Form("bft_utrailing[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_dtrailing",  event.bft_dtrailing,    Form("bft_dtrailing[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_utot",       event.bft_utot,         Form("bft_utot[%d][%d]/D",
                                                              NumOfSegBFT, MaxDepth));
  tree->Branch("bft_dtot",       event.bft_dtot,         Form("bft_dtot[%d][%d]/D",
                                                              NumOfSegBFT, MaxDepth));
  tree->Branch("bft_udepth",     event.bft_udepth,       Form("bft_udepth[%d]/I", NumOfSegBFT));
  tree->Branch("bft_ddepth",     event.bft_ddepth,       Form("bft_ddepth[%d]/I", NumOfSegBFT));
#endif
  tree->Branch("bft_ncl",       &event.bft_ncl,          "bft_ncl/I");
  tree->Branch("bft_clsize",     event.bft_clsize,       "bft_clsize[bft_ncl]/I");
  tree->Branch("bft_ctime",      event.bft_ctime,        "bft_ctime[bft_ncl]/D");
  tree->Branch("bft_ctot",       event.bft_ctot,         "bft_ctot[bft_ncl]/D");
  tree->Branch("bft_clpos",      event.bft_clpos,        "bft_clpos[bft_ncl]/D");
  tree->Branch("bft_clseg",      event.bft_clseg,        "bft_clseg[bft_ncl]/D");

  tree->Branch("aft_nhits", event.aft_nhits, Form("aft_nhits[%d]/I", NumOfPlaneAFT));
  tree->Branch("aft_hitpat", event.aft_hitpat,
               Form("aft_hitpat[%d][%d]/I", NumOfPlaneAFT, NumOfSegAFT));
  tree->Branch("aft_adc_high", event.aft_adc_high,
               Form("aft_adc_high[%d][%d][%d]/D", NumOfPlaneAFT, NumOfSegAFT, kUorD));
  tree->Branch("aft_adc_low", event.aft_adc_low,
               Form("aft_adc_low[%d][%d][%d]/D", NumOfPlaneAFT, NumOfSegAFT, kUorD));
  tree->Branch("aft_tdc", event.aft_tdc,
               Form("aft_tdc[%d][%d][%d][%d]/D",
                    NumOfPlaneAFT, NumOfSegAFT, kUorD, MaxDepth));
  tree->Branch("aft_mt", event.aft_mt,
               Form("aft_mt[%d][%d][%d]/D", NumOfPlaneAFT, NumOfSegAFT, MaxDepth));
  tree->Branch("aft_cmt", event.aft_cmt,
               Form("aft_cmt[%d][%d][%d]/D", NumOfPlaneAFT, NumOfSegAFT, MaxDepth));
  tree->Branch("aft_mtot", event.aft_mtot,
               Form("aft_mtot[%d][%d][%d]/D", NumOfPlaneAFT, NumOfSegAFT, MaxDepth));
  tree->Branch("aft_de_high", event.aft_de_high,
               Form("aft_de_high[%d][%d]/D", NumOfPlaneAFT, NumOfSegAFT));
  tree->Branch("aft_de_low", event.aft_de_low,
               Form("aft_de_low[%d][%d]/D", NumOfPlaneAFT, NumOfSegAFT));

  // HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")    &&
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
