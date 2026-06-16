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
#include "RawData.hh"
#include "TPCRawData.hh"
#include "HistTools.hh"
#include "UnpackerManager.hh"
#include "TPCAnalyzer.hh"
#include "TPCEventAnalyzer.hh"
#include "TPCHit.hh"
#include "TPCPadHelper.hh"
#include "TPCParamMan.hh"
#include "TPCRawHit.hh"
#include "TPCPositionCorrector.hh"

namespace
{
const auto qnan = TMath::QuietNaN();
using hddaq::unpacker::GUnpacker;
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser     = UserParamMan::GetInstance();

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

Int_t nhrawTpc;
std::vector<Int_t>    layerTpc;
std::vector<Int_t>    rowTpc;
std::vector<Int_t>    padTpc;
std::vector<Double_t> rawrmsTpc;
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

  nhrawTpc = 0;
  layerTpc.clear();
  rowTpc.clear();
  padTpc.clear();
  rawrmsTpc.clear();

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessNormal()
{
  using root::HF1;
  using root::HF2;

  if(event_number%10000 == 0)std::cout<<"Event : "<<event_number<<std::endl;
  
  static const Int_t MaxMultiHitTPC = gUser.GetParameter("MaxMultiHitTPC");
  static const Int_t NumOfTimeBucket = gUser.GetParameter("NumOfTimeBucket");

  TPCRawData TPCrawData;
  
  TPCrawData.DecodeTPCHits();

  HF1("Status", 0);
  
  for(Int_t layer = 0;layer<NumOfLayersTPC;++layer){
    auto hc = TPCrawData.GetTPCRawHits(layer);
    for(const auto&rhit : hc){
      Int_t row = rhit->RowId();
      Int_t asad = tpc::GetASADId(layer,row);
      Int_t pad = tpc::GetPadId(layer,row);
      Double_t rawrms = rhit->RMS(0, NumOfTimeBucket);
      layerTpc.push_back(layer);
      rowTpc.push_back(row);
      padTpc.push_back(pad);
      rawrmsTpc.push_back(rawrms);
      ++nhrawTpc;
    }
  }
  
  HF1("Status", 1);
  RawData rawData;
  
  for(const auto& hit: rawData.GetHodoRawHC("TriggerFlag")){
    trig_flag.push_back(hit->GetArrayTdc());
    trig_pat.push_back(hit->SegmentId());
  }
  HF1("Status", 2);
  
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

  tree = new TTree("tpc", "tree of TPCHit");
  tree->Branch("run_number", &run_number);
  tree->Branch("event_number", &event_number);
  tree->Branch("beam_flag", &beam_flag, "beam_flag/I");
  tree->Branch("trig_flag", &trig_flag);
  tree->Branch("trig_pat", &trig_pat);

  tree->Branch("nhrawTpc",&nhrawTpc);
  tree->Branch("layerTpc", &layerTpc);
  tree->Branch("rowTpc", &rowTpc);
  tree->Branch("padTpc", &padTpc);
  tree->Branch("rawrmsTpc", &rawrmsTpc);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO") &&
     InitializeParameter<HodoParamMan>("HDPRM") &&
     InitializeParameter<HodoPHCMan>("HDPHC") &&
     InitializeParameter<TPCParamMan>("TPCPRM") &&
     InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
