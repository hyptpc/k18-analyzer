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

  std::vector<Double_t> adc;
  std::vector<std::vector<Double_t>> tdc;

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum      = 0;
  spill      = 0;
  adc.clear();
  tdc.clear();
  adc.resize(NumOfSegHODO);
  tdc.resize(NumOfSegHODO);
}

//_____________________________________________________________________________
namespace root
{
Event  event;
TH1*   h[MaxHist];
TTree* tree;
enum eDetHid {
  HODOHid = 10000,
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
  // static const auto MinTdcBH1 = gUser.GetParameter("TdcBH1", 0);
  // static const auto MaxTdcBH1 = gUser.GetParameter("TdcBH1", 1);

  RawData rawData;
  // HodoAnalyzer hodoAna(rawData);

  event.evnum = gUnpacker.get_event_number();

  HF1(1, 0);

  { ///// HODO
    static const auto device_id = gUnpacker.get_device_id("HODO");
    static const auto adc_id = gUnpacker.get_data_id("HODO", "adc");
    static const auto tdc_id = gUnpacker.get_data_id("HODO", "tdc");
    for(Int_t seg=0; seg<NumOfSegHODO; ++seg){
      auto nhit = gUnpacker.get_entries(device_id, 0, seg, 0, adc_id);
      UInt_t adc = 0;
      if (nhit != 0) {
	adc = gUnpacker.get(device_id, 0, seg, 0, adc_id);
        HF1(HODOHid+1000+seg, adc);
        event.adc[seg] = adc;
      }
      Bool_t hit_flag = false;
      for(Int_t m=0, n=gUnpacker.get_entries(device_id, 0, seg, 0, tdc_id);
          m<n; ++m) {
        auto tdc = gUnpacker.get(device_id, 0, seg, 0, tdc_id, m);
        if (tdc != 0) {
          event.tdc[seg].push_back(tdc);
          HF1(HODOHid+3000+seg, tdc);
          hit_flag = true;
        }
        if (hit_flag) {
          HF1(HODOHid+2000+seg, adc);
        }
      }
    }
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

  // HODO
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1(HODOHid+1000+i, Form("HODO ADC %d", i), NbinAdc, MinAdc, MaxAdc);
    HB1(HODOHid+2000+i, Form("HODO ADCwTDC %d", i), NbinAdc, MinAdc, MaxAdc);
    HB1(HODOHid+3000+i, Form("HODO TDC %d", i), NbinTdcHr, MinTdcHr, MaxTdcHr);
  }

  //Tree
  HBTree("tree","tree of Counter");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("spill",     &event.spill,     "spill/I");
  tree->Branch("adc", &event.adc);
  tree->Branch("tdc", &event.tdc);

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
