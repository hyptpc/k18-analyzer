// -*- C++ -*-

#if ! defined E73_2024

#include "TPCEventAnalyzer.hh"

#include <Unpacker.hh>
#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <UnpackerXMLReadDigit.hh>
#include <DAQNode.hh>

#include "DetectorID.hh"
#include "TPCHit.hh"
#include "TPCAnalyzer.hh"
#include "TPCRawHit.hh"
#include "TPCRawData.hh"
#include "TPCPadHelper.hh"
// #include "TPCParamMan.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"

namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const auto& gUConf = hddaq::unpacker::GConfig::get_instance();
const auto& gUser = UserParamMan::GetInstance();
using root::HF1;
using root::HF2;
using root::HF2Poly;

}

//_____________________________________________________________________________
TPCEventAnalyzer::TPCEventAnalyzer()
{
}

//_____________________________________________________________________________
TPCEventAnalyzer::~TPCEventAnalyzer()
{
}

//_____________________________________________________________________________
void
TPCEventAnalyzer::TPCRawHit(const TPCRawData& TPCrawData)
{
  static const Int_t NumOfTimeBucket = gUser.GetParameter("NumOfTimeBucket");
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
}

//_____________________________________________________________________________
TPCBaselineInfo
TPCEventAnalyzer::TPCBaselineHit(const TPCRawData &TPCrawData){
  static const Int_t NumOfTimeBucket = gUser.GetParameter("NumOfTimeBucket");
  TPCBaselineInfo baseout;
  
  auto baseline = TPCrawData.GetBaselineTPC();
  if(baseline){
    baseout.row = baseline->RowId();
    baseout.layer = baseline->LayerId();
    baseout.rms = baseline->RMS(0, NumOfTimeBucket);
    baseout.valid = kTRUE;

    auto fadc = baseline->Fadc();
    for(Int_t tb = 0, ntb = fadc.size(); tb < ntb; ++tb){
      HF2("TPC_FADC_Baseline", tb, fadc.at(tb));
    }
  }

  return baseout;
}


//_____________________________________________________________________________
Int_t
TPCEventAnalyzer::TPCCorHit(const TPCRawData &TPCrawData){
  static const Int_t NumOfTimeBucket = gUser.GetParameter("NumOfTimeBucket");
  Int_t npadTpc = 0;
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

  return npadTpc;
}

//_____________________________________________________________________________
void
TPCEventAnalyzer::TPCHit(const TPCAnalyzer& TPCAnalyzer){
}

#endif
