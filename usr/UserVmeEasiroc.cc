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
  RawData rawData;
  // rawData.DecodeHits("TFlag");
  rawData.DecodeHits("VMEEASIROC");
  // HodoAnalyzer hodoAna(rawData);
  rawData.Print();

  event.evnum = gUnpacker.get_event_number();

  HF1(1, 0);

  ///// Trigger Flag
  // std::bitset<NumOfSegTrig> trigger_flag;
  // {
  //   for(const auto& hit: rawData.GetHodoRawHitContainer("TFlag")){
  //     Int_t seg = hit->SegmentId();
  //     Int_t tdc = hit->GetTdc();
  //     if(tdc > 0){
  //       event.trigpat[trigger_flag.count()] = seg;
  //       event.trigflag[seg] = tdc;
  //       trigger_flag.set(seg);
  //       HF1(10, seg);
  //       HF1(10+seg, tdc);
  //     }
  //   }
  // }

  // if(trigger_flag[trigger::kSpillEnd])
  //   return true;

  HF1(1, 1);

  for(const auto& hit: rawData.GetHodoRawHitContainer("VMEEASIROC")){
    hit->Print();
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
      HB1(AFTHid+plane*1000+9+ud, Form("AFT AdcLow %s Plane#%d", s, plane), NbinAdc, MinAdc, MaxAdc);
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
