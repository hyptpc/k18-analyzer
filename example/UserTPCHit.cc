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

Int_t npadTpc;   // number of pads
Int_t nhTpc;     // number of hits

Int_t browTpc;
Int_t blayerTpc;
Double_t brmsTpc;//Baseline RMS;
// vector (size=nhTpc)
std::vector<Int_t>    layerTpc;  // layer id
std::vector<Int_t>    rowTpc;    // row id
std::vector<Int_t>    padTpc;    // pad id
std::vector<Double_t> pedTpc;    // pedestal
std::vector<Double_t> rmsTpc;    // rms
std::vector<Double_t> rawrmsTpc; // rawrms
std::vector<Double_t> deTpc;     // dE
std::vector<Double_t> sigmaTpc;  // sigma
std::vector<Double_t> tTpc;      // time
std::vector<Double_t> chisqrTpc; // chi^2 of signal fitting
std::vector<Double_t> cdeTpc;    // dE
std::vector<Double_t> ctTpc;     // time
std::vector<Double_t> dlTpc;     // time
std::vector<Double_t> clkTpc;    // clock timing
std::vector<Double_t> cobo_id;    // COBO ID
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

  npadTpc = 0;
  nhTpc   = 0;
  browTpc   = -1;
  blayerTpc = -1;
  brmsTpc   = 0;
  trig_pat.clear();
  trig_flag.clear();
  layerTpc.clear();
  rowTpc.clear();
  padTpc.clear();
  pedTpc.clear();
  rmsTpc.clear();
  rawrmsTpc.clear();
  deTpc.clear();
  sigmaTpc.clear();
  tTpc.clear();
  chisqrTpc.clear();
  cdeTpc.clear();
  ctTpc.clear();
  dlTpc.clear();
  clkTpc.clear();
  cobo_id.clear();

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessNormal()
{
  using root::HF1;
  using root::HF2;

  static const Int_t MaxMultiHitTPC = gUser.GetParameter("MaxMultiHitTPC");
  static const Int_t NumOfTimeBucket = gUser.GetParameter("NumOfTimeBucket");

  RawData rawData;
  rawData.DecodeHits("COBO");

  TPCRawData TPCrawData;
  TPCrawData.DecodeTPCHits();
  
  HodoAnalyzer hodoAna(rawData);
  hodoAna.DecodeHits("COBO", !HasHodoGroup(HodoGroupMask[kCOBO], HodoGroup::NoCluster));

  HF1("Status", 0);
  for(const auto& hit: rawData.GetHodoRawHC("TriggerFlag")){
    trig_flag.push_back(hit->GetArrayTdc());
    trig_pat.push_back(hit->SegmentId());
  }
  HF1("Status", 1);

  { ///// COBO
    static const TString n("COBO");
    for(const auto& hit: rawData.GetHodoRawHC(n)){
      for(const auto& tdc: hit->GetArrayTdc()){
        //hit->Print();
        HF1("TPC_Clock_TDC", tdc);
      }
    }

    for(Int_t i=0, nh=hodoAna.GetNHits(n); i<nh; ++i){
      const auto& hit = hodoAna.GetHit(n, i);
      cobo_id.push_back(hit->SegmentId());
      auto clock_timing = hit->GetArrayTime().at(0);
      clkTpc.push_back(clock_timing);
      HF1("TPC_Clock_Time", clock_timing);
    }
  }

  HF1("Status", 2);

  //________________________________________________________
  //___ TPCRawHit
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    auto hc = TPCrawData.GetTPCRawHits(layer);
    for(const auto& rhit : hc){
      auto mean    = rhit->Mean(0, NumOfTimeBucket);
      auto max_adc = rhit->MaxAdc(0, NumOfTimeBucket);
      auto min_adc = rhit->MinAdc(0, NumOfTimeBucket);
      auto rms     = rhit->RMS(0, NumOfTimeBucket);
      auto loc_max = rhit->LocMax(0, NumOfTimeBucket);

      HF1("TPC_FADC_Mean", mean);
      HF1("TPC_FADC_Max", max_adc);
      HF1("TPC_FADC_RMS", rms);
      HF1("TPC_FADC_LocMax", loc_max);
      HF1("TPC_FADC_Min", min_adc);

      auto fadc = rhit->Fadc();
      for(Int_t tb = 0, ntb = fadc.size(); tb < ntb; ++tb){
        HF2("TPC_FADC_Before", tb, fadc.at(tb));
      }
    }
  }

  HF1("Status", 3);

  //________________________________________________________
  //___ TPCRawHit after baseline correction
  auto baseline = TPCrawData.GetBaselineTPC();
  if(baseline){
    browTpc = baseline->RowId();
    blayerTpc = baseline->LayerId();
    brmsTpc = baseline->RMS(0, NumOfTimeBucket);

    auto fadc = baseline->Fadc();
    for(Int_t tb = 0, ntb = fadc.size(); tb < ntb; ++tb){
      HF2("TPC_FADC_Baseline", tb, fadc.at(tb));
    }
  }

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    auto hc = TPCrawData.GetTPCCorHits(layer);
    const auto nhit = hc.size();
    npadTpc += nhit;

    for(const auto& rhit : hc){
      auto mean    = rhit->Mean(0, NumOfTimeBucket);
      auto max_adc = rhit->MaxAdc(0, NumOfTimeBucket);
      auto min_adc = rhit->MinAdc(0, NumOfTimeBucket);
      auto rms     = rhit->RMS(0, NumOfTimeBucket);
      auto loc_max = rhit->LocMax(0, NumOfTimeBucket);
      auto pars    = rhit->GetParameters();

      HF1("TPC_FADC_Mean_Cor", mean);
      HF1("TPC_FADC_Max_Cor", max_adc);
      HF1("TPC_FADC_RMS_Cor", rms);
      HF1("TPC_FADC_LocMax_Cor", loc_max);
      HF1("TPC_FADC_Min_Cor", min_adc);
      HF1("TPC_FADC_Baseline_p0", pars.at(0));
      HF1("TPC_FADC_Baseline_p1", pars.at(1));
      HF1("TPC_FADC_Baseline_p2", pars.at(2));

      // 2D FADC waveform after correction
      auto fadc = rhit->Fadc();
      for(Int_t tb = 0, ntb = fadc.size(); tb < ntb; ++tb){
      	HF2("TPC_FADC_After", tb, fadc.at(tb));
      }
    }
  }

  HF1("TPC_Multiplicity_Raw", npadTpc);
  HF1("Status", 4);

  //________________________________________________________
  //___ TPCHit
  if(MaxMultiHitTPC>0 && npadTpc>MaxMultiHitTPC){
    std::cout << "#W Too many hits found, npadTpc = " << npadTpc << std::endl;
    return true;
  }

  TPCAnalyzer TPCAna;
  TPCAna.DecodeTPCHits(TPCrawData, clkTpc);
  HF1("Status", 5);

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    auto hc = TPCAna.GetTPCHC(layer);
    for(const auto& hit : hc){
      if(!hit || !hit->IsGood())
        continue;
      //Int_t layer = hit->GetLayer();
      Int_t row = hit->GetWire();
      Int_t pad = tpc::GetPadId(layer, row);
      Double_t ped = hit->GetPedestal();
      Double_t rms = hit->GetRMS();
      Double_t rawrms = hit->GetRawRMS();

      HF1("TPC_Pedestal", ped);
      HF1("TPC_RMS", rms);

      const auto& vec = tpc::getPosition(pad);
      //HF2Poly(1001, vec.Z(), vec.X());
      Int_t nhit = hit->GetNHits();
      Bool_t good_for_analysis = false;
      for(Int_t i=0; i<nhit; ++i){
        Double_t cde = hit->GetCDe(i);
        Double_t de = hit->GetDe(i);
        Double_t time = hit->GetTime(i);
        Double_t chisqr = hit->GetChisqr(i);
        Double_t ctime = hit->GetCTime(i);
        Double_t dl = hit->GetDriftLength(i);
        Double_t sigma = hit->GetSigma(i);
        layerTpc.push_back(layer);
        rowTpc.push_back(row);
        padTpc.push_back(pad);
        pedTpc.push_back(ped);
        rmsTpc.push_back(rms);
        rawrmsTpc.push_back(rawrms);
        deTpc.push_back(de);
        tTpc.push_back(time);
        chisqrTpc.push_back(chisqr);
        cdeTpc.push_back(cde);
        ctTpc.push_back(ctime);
        dlTpc.push_back(dl);
        sigmaTpc.push_back(sigma);

        HF1("TPC_DeltaE", de);
        HF1("TPC_Time", time);
        HF1("TPC_Chisqr", chisqr);
        HF1("TPC_CDeltaE", cde);
        HF1("TPC_CTime", ctime);
        HF1("TPC_DriftLength", dl);
        HF1("TPC_sigma", sigma);
        HF2("TPC_sigma%%de", de, sigma);
        HF2("TPC_time%%de", de, time);

        good_for_analysis = true;
        ++nhTpc;
      }
      if(good_for_analysis){
        auto fadc = hit->GetRawHit()->Fadc();
        for(Int_t tb = 0, ntb = fadc.size(); tb < ntb; ++tb){
          HF2("TPC_FADC_Good", tb, fadc.at(tb));
        }
      }
    }
  }
  nhTpc = nhTpc;

  HF1("Status", 6);

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
  hist::BuildTPCHit();

  tree = new TTree("tpc", "tree of TPCHit");
  tree->Branch("run_number", &run_number);
  tree->Branch("event_number", &event_number);
  tree->Branch("beam_flag", &beam_flag, "beam_flag/I");
  tree->Branch("trig_flag", &trig_flag);
  tree->Branch("trig_pat", &trig_pat);

  tree->Branch("npadTpc", &npadTpc);
  tree->Branch("nhTpc", &nhTpc);
  tree->Branch("browTpc", &browTpc);
  tree->Branch("blayerTpc", &blayerTpc);
  tree->Branch("brmsTpc", &brmsTpc);
  tree->Branch("layerTpc", &layerTpc);
  tree->Branch("rowTpc", &rowTpc);
  tree->Branch("padTpc", &padTpc);
  tree->Branch("pedTpc", &pedTpc);
  tree->Branch("rmsTpc", &rmsTpc);
  tree->Branch("rawrmsTpc", &rawrmsTpc);
  tree->Branch("deTpc", &deTpc);
  tree->Branch("tTpc", &tTpc);
  tree->Branch("chisqrTpc", &chisqrTpc);
  tree->Branch("cdeTpc", &cdeTpc);
  tree->Branch("ctTpc", &ctTpc);
  tree->Branch("dlTpc", &dlTpc);
  tree->Branch("sigmaTpc", &sigmaTpc);
  tree->Branch("clkTpc", &clkTpc);
  tree->Branch("cobo_id", &cobo_id);

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
