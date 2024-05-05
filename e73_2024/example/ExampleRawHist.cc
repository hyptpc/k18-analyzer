// -*- C++ -*-

#include <iostream>
#include <sstream>
#include <cmath>

#include <TString.h>

#include <Unpacker.hh>
#include <UnpackerManager.hh>
#include <DAQNode.hh>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "RootHelper.hh"
#include "MTDCRawHit.hh"
#include "HodoRawHit.hh"
#include "DCRawHit.hh"
#include "RawData.hh"
#include "UserParamMan.hh"
#include "VEvent.hh"
#include "HistTools.hh"

namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
Bool_t
ProcessBegin()
{
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessNormal()
{
  RawData rawData;
  rawData.DecodeHits();

  // Trigger flag
  {
    static const Char_t* name = "TriggerFlag";
    for(const auto& hit: rawData.GetHodoRawHC(name)){
      auto seg = hit->SegmentId();
      if(seg > NumOfSegTrigFlag)
        continue;
      Bool_t is_hit = false;
      for(const auto& tdc: hit->GetArrayTdcUp()){
	root::HF1(Form("%s_TDC_seg%d", name, seg), tdc);
        is_hit = true;
      }
      if(is_hit){
	root::HF1(Form("%s_HitPat", name), seg);
      }
    }
  }

  // BHT
  {
    static const Char_t* name = "BHT";
    Int_t multi_or = 0;
    Int_t multi_and = 0;
    for(const auto& hit: rawData.GetHodoRawHC(name)){
      auto seg = hit->SegmentId();
      Int_t ud_good = 0;
      for(Int_t ud=0; ud<2; ++ud){
        for(Int_t i=0, n=hit->GetSizeTdcLeading(ud); i<n; ++i){
          Double_t t = hit->GetTdc(ud, i);
          if(gUser.IsInRange(Form("%s_TDC", name), t))
            ++ud_good;
          root::HF1(Form("%s_TDC_seg%d%s", name, seg, UorD[ud]), t);
        }
        for(Int_t i=0, n=hit->GetSizeTdcTrailing(ud); i<n; ++i){
          Double_t t = hit->GetTdcTrailing(ud, i);
          root::HF1(Form("%s_Trailing_seg%d%s", name, seg, UorD[ud]), t);
          if(n == hit->GetSizeTdcLeading(ud)){
            Double_t l = hit->GetTdc(ud, i);
            if(gUser.IsInRange(Form("%s_TDC", name), l)){
              root::HF1(Form("%s_TOT_seg%d%s", name, seg, UorD[ud]), l - t);
            }
          }
        }
      }
      if(ud_good >= 1){
        root::HF1(Form("%s_HitPat_OR", name), seg);
        ++multi_or;
      }
      if(ud_good == 2){
        root::HF1(Form("%s_HitPat_AND", name), seg);
        ++multi_and;
      }
    }
    root::HF1(Form("%s_Multi_OR", name), multi_or);
    root::HF1(Form("%s_Multi_AND", name), multi_and);
  }

  // AC
  {
    static const Char_t* name = "AC";
    Int_t multi = 0;
    Bool_t is_good = false;
    for(const auto& hit: rawData.GetHodoRawHC(name)){
      Int_t seg = hit->SegmentId();
      for(const auto& t: hit->GetArrayTdcLeading()){
        if(gUser.IsInRange(Form("%s_TDC", name), t)){
          is_good = true;
          root::HF1(Form("%s_HitPat", name), seg);
          ++multi;
        }
	root::HF1(Form("%s_TDC_seg%d", name, seg), t);
      }
    }
    for(const auto& hit: rawData.GetHodoRawHC(name)){
      Int_t seg = hit->SegmentId();
      Double_t a = hit->GetAdc();
      root::HF1(Form("%s_ADC_seg%d", name, seg), a);
      if(is_good)
        root::HF1(Form("%s_AWT_seg%d", name, seg), a);
    }
    root::HF1(Form("%s_Multi", name), multi);
  }

  // Hodoscope
  for(Int_t ihodo=kT1; ihodo<kNumHodo;++ihodo){
    const Char_t* name = NameHodo[ihodo];
    Int_t multi_or = 0;
    Int_t multi_and = 0;
    for(const auto& hit: rawData.GetHodoRawHC(name)){
      auto seg = hit->SegmentId();
      Int_t ud_good = 0;
      for(Int_t ud=0; ud<2; ++ud){
        Bool_t is_good = false;
        for(const auto& t: hit->GetArrayTdc(ud)){
          if(gUser.IsInRange(Form("%s_TDC", name), t))
            is_good = true;
          root::HF1(Form("%s_TDC_seg%d%s", name, seg, UorD[ud]), t);
        }
        for(const auto& a: hit->GetArrayAdc(ud)){
          root::HF1(Form("%s_ADC_seg%d%s", name, seg, UorD[ud]), a);
          if(is_good)
            root::HF1(Form("%s_AWT_seg%d%s", name, seg, UorD[ud]), a);
        }
        ud_good += is_good;
      }
      if(ud_good >= 1){
        root::HF1(Form("%s_HitPat_OR", name), seg);
        ++multi_or;
      }
      if(ud_good == 2){
        root::HF1(Form("%s_HitPat_AND", name), seg);
        ++multi_and;
      }
    }
    root::HF1(Form("%s_Multi_OR", name), multi_or);
    root::HF1(Form("%s_Multi_AND", name), multi_and);
  }

  // DC
  for(Int_t idc=0; idc<=kVFT; ++idc){
    const Char_t* name = NameDC[idc].Data();
    Int_t nlayer = NumOfLayerDC[idc];
    for(Int_t layer=0; layer<nlayer; ++layer){
      const auto& cont = rawData.GetDCRawHC(DetIdDC[idc], layer);
      // hist::H1(Form("%s_Mul_layer%d",name,layer),nh,mulbins);
      Int_t multi = 0;
      Int_t cmulti = 0;
      for(Int_t i=0, n=cont.size(); i<n; ++i){
        auto hit = cont[i];
        if(!hit) continue;
        Int_t wire = hit->WireId();
        Bool_t is_good = false;
        for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
          if(m != hit->GetTrailingSize()) break;
          Double_t l = hit->GetTdc(j);
          Double_t t = hit->GetTrailing(j);
          Double_t tot = l - t;
          if(gUser.IsInRange(Form("%s_TOT", name), tot))
            is_good = true;
          for(const auto& totcut: std::vector<TString>{"", "C"}){
            if(totcut == "C" && !gUser.IsInRange(Form("%s_TOT", name), tot))
              continue;
            auto c = totcut.Data();
            root::HF1(Form("%s_%sTDC_layer%d", name, c, layer), l);
            root::HF1(Form("%s_%sTrailing_layer%d", name, c, layer), t);
            root::HF1(Form("%s_%sTOT_layer%d", name, c, layer), tot);
            if(j == 0){
              root::HF1(Form("%s_%sTDC1st_layer%d", name, c, layer), l);
              root::HF1(Form("%s_%sTrailing1st_layer%d", name, c, layer), t);
              root::HF1(Form("%s_%sTOT1st_layer%d", name, c, layer), tot);
            }
          }
        }
        if(hit->GetTdcSize() > 0){
          root::HF1(Form("%s_HitPat_layer%d", name, layer), wire);
          ++multi;
        }
        if(is_good){
          root::HF1(Form("%s_CHitPat_layer%d", name, layer), wire);
          ++cmulti;
        }
      }
      root::HF1(Form("%s_Multi_layer%d", name, layer), multi);
      root::HF1(Form("%s_CMulti_layer%d", name, layer), cmulti);
    }
  }

  // DAQ
  static const auto k_data_size = hddaq::unpacker::DAQNode::k_data_size;
  Int_t vme_index = 0;
  Int_t hul_index = 0;
  Int_t vea0c_index = 0;
  for(auto&& c : gUnpacker.get_root()->get_child_list()){
    if (!c.second) continue;
    TString name = c.second->get_name();
    auto node_id = c.second->get_id();
    auto data_size = gUnpacker.get_node_header(node_id, k_data_size);
    if(name.Contains("vme")){
      root::HF2("FE_VME_DataSize", vme_index++, data_size);
    }
    if(name.Contains("hul")){
      root::HF2("FE_HUL_DataSize", hul_index++, data_size);
    }
    if(name.Contains("veasiroc")){
      root::HF2("FE_VEASIROC_DataSize", vea0c_index++, data_size);
    }
  }

  {
    auto node_id = gUnpacker.get_fe_id("k18breb");
    auto data_size = gUnpacker.get_node_header(node_id, k_data_size);
    root::HF1("EB_DataSize", data_size);
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessEnd()
{
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  Double_t hrtdcbins1[3] = {45000, 950000, 1400000};
  Double_t hrtdcbins2[3] = {50000, 0, 1000000}; // for CVC, NC
  Double_t hrtotbins[3] = {5000, 0, 50000};
  Double_t adcbins[3] = {4096, -0.5, 4095.5};
  Double_t mhtdcbins[3] = {2000, 0, 2000};
  Double_t mhtotbins[3] = {1000, 0, 1000};

  { // TriggerFlag
    const Char_t* name = "TriggerFlag";
    Double_t patbins[3]={NumOfSegTrigFlag, -0.5, NumOfSegTrigFlag-0.5};
    for(Int_t i=0; i<NumOfSegTrigFlag; ++i){
      root::HB1(Form("%s_TDC_seg%d", name, i), mhtdcbins);
    }
    root::HB1(Form("%s_HitPat; Segment; Counts", name), patbins);
  }

  { // BHT
    const Char_t* name = "BHT";
    Int_t nseg = NumOfSegBHT;
    for(Int_t i=0; i<nseg; ++i){
      for(const auto& ud : std::vector<TString>{"U", "D"}){
        root::HB1(Form("%s_TDC_seg%d%s", name, i, ud.Data()), hrtdcbins1);
        root::HB1(Form("%s_Trailing_seg%d%s", name, i, ud.Data()), hrtdcbins1);
        root::HB1(Form("%s_TOT_seg%d%s", name, i, ud.Data()), hrtotbins);
      }
    }
    root::HB1(Form("%s_HitPat_OR", name), nseg, -0.5, nseg - 0.5);
    root::HB1(Form("%s_Multi_OR", name), nseg + 1, -0.5, nseg + 0.5);
    root::HB1(Form("%s_HitPat_AND", name), nseg, -0.5, nseg - 0.5);
    root::HB1(Form("%s_Multi_AND", name), nseg + 1, -0.5, nseg + 0.5);
  }

  { // AC
    const Char_t* name = "AC";
    Int_t nseg = NumOfSegAC;
    for(Int_t i=0; i<nseg; ++i){
      root::HB1(Form("%s_ADC_seg%d", name, i), adcbins);
      root::HB1(Form("%s_AWT_seg%d", name, i), adcbins);
      root::HB1(Form("%s_TDC_seg%d", name, i), mhtdcbins);
    }
    root::HB1(Form("%s_HitPat", name), nseg, -0.5, nseg - 0.5);
    root::HB1(Form("%s_Multi", name), nseg + 1, -0.5, nseg + 0.5);
  }

  // Hodoscope
  for(Int_t ihodo=kT1; ihodo<kNumHodo;++ihodo){
    auto name = NameHodo[ihodo].Data();
    Double_t* hrtdcbins;
    if(NameHodo[ihodo].Contains("CVC") ||
       NameHodo[ihodo].Contains("NC")){
      hrtdcbins = hrtdcbins2;
    }else{
      hrtdcbins = hrtdcbins1;
    }
    Int_t nseg = NumOfSegHodo[ihodo];
    for(const auto& uord: std::vector<TString>{"U", "D"} ){
      auto ud = uord.Data();
      for(Int_t i=0; i<nseg; ++i){
        root::HB1(Form("%s_ADC_seg%d%s", name, i, ud), adcbins);
        root::HB1(Form("%s_AWT_seg%d%s", name, i, ud), adcbins);
        root::HB1(Form("%s_TDC_seg%d%s", name, i, ud), hrtdcbins);
      }
    }
    for(const auto& uord: std::vector<TString>{"OR", "AND"}){
      auto ud = uord.Data();
      root::HB1(Form("%s_HitPat_%s", name, ud), nseg, -0.5, nseg - 0.5);
      root::HB1(Form("%s_Multi_%s", name, ud), nseg + 1, -0.5, nseg + 0.5);
    }
  }

  // DC
  for(Int_t idc=0; idc<=kVFT; ++idc){
    const Char_t* name = NameDC[idc].Data();
    Int_t nlayer = NumOfLayerDC[idc];
    Int_t nwire = NumOfWireDC[idc];
    for(Int_t layer=0; layer<nlayer; ++layer){
      for(const auto& totcut: std::vector<TString>{"", "C"}){
        auto c = totcut.Data();
        root::HB1(Form("%s_%sTDC_layer%d", name, c, layer), mhtdcbins);
        root::HB1(Form("%s_%sTDC1st_layer%d", name, c, layer), mhtdcbins);
        root::HB1(Form("%s_%sTrailing_layer%d", name, c, layer), mhtdcbins);
        root::HB1(Form("%s_%sTrailing1st_layer%d", name, c, layer), mhtdcbins);
        root::HB1(Form("%s_%sTOT_layer%d", name, c, layer), mhtotbins);
        root::HB1(Form("%s_%sTOT1st_layer%d", name, c, layer), mhtotbins);
        root::HB1(Form("%s_%sHitPat_layer%d", name, c, layer), nwire, -0.5, nwire - 0.5);
        root::HB1(Form("%s_%sMulti_layer%d", name, c, layer), nwire + 1, -0.5, nwire + 0.5);
      }
    }
  }

  // DAQ
  std::vector<Int_t> vme_fe_id;
  std::vector<Int_t> hul_fe_id;
  std::vector<Int_t> vea0c_fe_id;
  root::HB1("EB_DataSize; words", 2000, 0, 20000);
  for(auto&& c : gUnpacker.get_root()->get_child_list()){
    if (!c.second) continue;
    TString name = c.second->get_name();
    auto node_id = c.second->get_id();
    if(name.Contains("vme"))
      vme_fe_id.push_back(node_id);
    if(name.Contains("hul"))
      hul_fe_id.push_back(node_id);
    if(name.Contains("veasiroc"))
      vea0c_fe_id.push_back(node_id);
  }
  root::HB2("FE_VME_DataSize; ; words",
            vme_fe_id.size(), 0, vme_fe_id.size(), 100, 0, 1000);
  for(Int_t i=0, n=vme_fe_id.size(); i<n; ++i){
    auto h1 = gDirectory->Get<TH2>("FE_VME_DataSize");
    h1->GetXaxis()->SetBinLabel(i+1, "0x"+TString::Itoa(vme_fe_id[i], 16));
  }
  root::HB2("FE_HUL_DataSize; ; words",
            hul_fe_id.size(), 0, hul_fe_id.size(), 100, 0, 2000);
  for(Int_t i=0, n=hul_fe_id.size(); i<n; ++i){
    auto h1 = gDirectory->Get<TH2>("FE_HUL_DataSize");
    h1->GetXaxis()->SetBinLabel(i+1, "0x"+TString::Itoa(hul_fe_id[i], 16));
  }
  root::HB2("FE_VEASIROC_DataSize; ; words",
            vea0c_fe_id.size(), 0, vea0c_fe_id.size(), 100, 0, 1000);
  for(Int_t i=0, n=vea0c_fe_id.size(); i<n; ++i){
    auto h1 = gDirectory->Get<TH2>("FE_VEASIROC_DataSize");
    h1->GetXaxis()->SetBinLabel(i+1, "0x"+TString::Itoa(vea0c_fe_id[i], 16));
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    InitializeParameter<UserParamMan>("USER") &&
    InitializeParameter<DCGeomMan>("DCGEO");
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
