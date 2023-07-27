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
#include "S2sLib.hh"
#include "RawData.hh"
#include "UnpackerManager.hh"

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

  // VMEEASIROC raw
  Int_t    vmeeasiroc_nhits[NumOfPlaneVMEEASIROC];
  Int_t    vmeeasiroc_hitpat[NumOfPlaneVMEEASIROC][NumOfSegVMEEASIROC];
  Double_t vmeeasiroc_tdc[NumOfPlaneVMEEASIROC][NumOfSegVMEEASIROC][MaxDepth];
  Double_t vmeeasiroc_adc_high[NumOfPlaneVMEEASIROC][NumOfSegVMEEASIROC];
  Double_t vmeeasiroc_adc_low[NumOfPlaneVMEEASIROC][NumOfSegVMEEASIROC];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum      = 0;
  for(Int_t it=0; it<NumOfSegTrig; it++){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t p=0; p<NumOfPlaneVMEEASIROC; p++){
    vmeeasiroc_nhits[p] = 0;
    for(Int_t seg=0; seg<NumOfSegVMEEASIROC; seg++){
      vmeeasiroc_hitpat[p][seg] = -1;
      vmeeasiroc_adc_high[p][seg] = qnan;
      vmeeasiroc_adc_low[p][seg] = qnan;
      for(Int_t i=0; i<MaxDepth; i++){
	vmeeasiroc_tdc[p][seg][i] = qnan;
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
  VMEEASIROCHid  = 200000,
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
  HodoAnalyzer hodoAna(rawData);

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

  ////////// VMEEASIROC
  for(const auto& hit: rawData.GetHodoRawHitContainer("VMEEASIROC")){
    const auto& U = HodoRawHit::kUp;
    // hit->Print();
    Int_t plane = hit->PlaneId();
    Int_t seg = hit->SegmentId();
    auto adc_high = hit->GetAdcHigh();
    auto adc_low = hit->GetAdcLow();
    event.vmeeasiroc_adc_high[plane][seg] = adc_high;
    event.vmeeasiroc_adc_low[plane][seg]  = adc_low;
    HF1(VMEEASIROCHid+plane*1000+7, adc_high);
    HF1(VMEEASIROCHid+plane*1000+9, adc_low);
    HF2(VMEEASIROCHid+plane*1000+15, seg, adc_high);
    HF2(VMEEASIROCHid+plane*1000+17, seg, adc_low);
    for(Int_t i=0, n=hit->GetSizeTdcLeading(); i<n; ++i){
      auto tdc = hit->GetTdc(U, i);
      event.vmeeasiroc_tdc[plane][seg][i] = tdc;
      HF1(VMEEASIROCHid+plane*1000+3, tdc);
      HF2(VMEEASIROCHid+plane*1000+11, seg, tdc);
      // HF1(AFTHid+plane*1000+seg+100+ud*100, tdc);
    }
  }

  // hodoAna.DecodeHits<FiberHit>("AFT");
  // for(Int_t i=0, n=hodoAna.GetNHits("AFT"); i<n; ++i){
  //   const auto& hit = hodoAna.GetHit<FiberHit>("AFT", i);
  //   hit->Print();
  //   const auto& rhit = hit->GetRawHit();
  //   rhit->Print();
  //   // Int_t plane = hit->PlaneId();
  //   // Int_t seg = hit->SegmentId();
  //   // event.aft_hitpat[plane][event.aft_nhits[plane]++] = seg;
  //   // HF1(AFTHid+plane*1000+2, seg);
  //   // Int_t m = hit->GetEntries();
  //   // for(Int_t j=0; j<m; ++j){
  //   //   auto mt = hit->MeanTime(j);
  //   //   auto cmt = hit->CMeanTime(j);
  //   //   auto mtot = hit->MeanTOT(j);
  //   //   event.aft_mt[plane][seg][j] = mt;
  //   //   event.aft_cmt[plane][seg][j] = cmt;
  //   //   event.aft_mtot[plane][seg][j] = mtot;
  //   //   HF1(AFTHid+plane*1000+21, mt);
  //   //   HF1(AFTHid+plane*1000+22, cmt);
  //   //   HF1(AFTHid+plane*1000+23, mtot);
  //   //   HF2(AFTHid+plane*1000+31, seg, mt);
  //   //   HF2(AFTHid+plane*1000+32, seg, cmt);
  //   //   HF2(AFTHid+plane*1000+33, seg, mtot);
  //   //   for(Int_t ud=0; ud<kUorD; ++ud){
  //   //     auto tot = hit->TOT(ud, j);
  //   //     HF1(AFTHid+plane*1000+5+ud, tot);
  //   //     HF2(AFTHid+plane*1000+13+ud, seg, tot);
  //   //     // HF1(AFTHid+plane*1000+seg+300+ud*100, tot);
  //   //   }
  //   // }
  //   // auto de_high = hit->DeltaEHighGain();
  //   // auto de_low = hit->DeltaELowGain();
  //   // event.aft_de_high[plane][seg] = de_high;
  //   // event.aft_de_low[plane][seg] = de_low;
  //   // HF1(AFTHid+plane*1000+24, de_high);
  //   // HF1(AFTHid+plane*1000+25, de_low);
  //   // HF2(AFTHid+plane*1000+34, seg, de_high);
  //   // HF2(AFTHid+plane*1000+35, seg, de_low);
  // }
  // // for(Int_t plane=0; plane<NumOfPlaneAFT; ++plane){
  // //   HF1(AFTHid+plane*1000+1, event.aft_nhits[plane]);
  // // }

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

  //VMEEASIROC
  for(Int_t plane=0; plane<NumOfPlaneVMEEASIROC; ++plane){
    HB1(VMEEASIROCHid+plane*1000+1, Form("VMEEASIROC Nhits Plane#%d", plane),
	NumOfSegVMEEASIROC, 0., NumOfSegVMEEASIROC);
    HB1(VMEEASIROCHid+plane*1000+2, Form("VMEEASIROC Hitpat Plane#%d", plane),
	NumOfSegVMEEASIROC, 0., NumOfSegVMEEASIROC);
    HB1(VMEEASIROCHid+plane*1000+3, Form("VMEEASIROC Tdc Plane#%d", plane),
	NbinTdc, MinTdc, MaxTdc);
    HB1(VMEEASIROCHid+plane*1000+5, Form("VMEEASIROC Tot Plane#%d", plane),
	NbinTdc, MinTdc, MaxTdc);
    HB1(VMEEASIROCHid+plane*1000+7, Form("VMEEASIROC AdcHigh Plane#%d", plane),
	NbinAdc, MinAdc, MaxAdc);
    HB1(VMEEASIROCHid+plane*1000+9, Form("VMEEASIROC AdcLow Plane#%d", plane),
	NbinAdc, MinAdc, MaxAdc);
    HB2(VMEEASIROCHid+plane*1000+11, Form("VMEEASIROC Tdc %%Seg Plane#%d", plane),
	NumOfSegVMEEASIROC, 0., NumOfSegVMEEASIROC, NbinTdc, MinTdc, MaxTdc);
    HB2(VMEEASIROCHid+plane*1000+13, Form("VMEEASIROC Tot %%Seg Plane#%d", plane),
	NumOfSegVMEEASIROC, 0., NumOfSegVMEEASIROC, NbinTot, MinTot, MaxTot);
    HB2(VMEEASIROCHid+plane*1000+15, Form("VMEEASIROC AdcHigh %%Seg Plane#%d", plane),
	NumOfSegVMEEASIROC, 0., NumOfSegVMEEASIROC, NbinAdc, MinAdc, MaxAdc);
    HB2(VMEEASIROCHid+plane*1000+17, Form("VMEEASIROC AdcLow %%Seg Plane#%d", plane),
	NumOfSegVMEEASIROC, 0., NumOfSegVMEEASIROC, NbinAdc, MinAdc, MaxAdc);
  }

  //Tree
  HBTree("ea0c", "tree of Easiroc");
  //Trig
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  tree->Branch("vmeeasiroc_nhits", event.vmeeasiroc_nhits,
	       Form("vmeeasiroc_nhits[%d]/I", NumOfPlaneVMEEASIROC));
  tree->Branch("vmeeasiroc_hitpat", event.vmeeasiroc_hitpat,
               Form("vmeeasiroc_hitpat[%d][%d]/I", NumOfPlaneVMEEASIROC, NumOfSegVMEEASIROC));
  tree->Branch("vmeeasiroc_adc_high", event.vmeeasiroc_adc_high,
               Form("vmeeasiroc_adc_high[%d][%d]/D", NumOfPlaneVMEEASIROC, NumOfSegVMEEASIROC));
  tree->Branch("vmeeasiroc_adc_low", event.vmeeasiroc_adc_low,
               Form("vmeeasiroc_adc_low[%d][%d]/D", NumOfPlaneVMEEASIROC, NumOfSegVMEEASIROC));
  tree->Branch("vmeeasiroc_tdc", event.vmeeasiroc_tdc,
               Form("vmeeasiroc_tdc[%d][%d][%d]/D",
                    NumOfPlaneVMEEASIROC, NumOfSegVMEEASIROC, MaxDepth));

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
