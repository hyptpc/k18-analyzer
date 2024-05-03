// -*- C++ -*-

#include <iostream>
#include <sstream>
#include <cmath>

#include <TString.h>

#include <Unpacker.hh>
#include <UnpackerManager.hh>
#include <DAQNode.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RootHelper.hh"
#include "MTDCRawHit.hh"
#include "HodoRawHit.hh"
#include "DCRawHit.hh"
#include "RawData.hh"
#include "CDCWireMapMan.hh"
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
    const auto& cont = rawData.GetMTDCRawHC(DetIdTrigFlag);
    for(Int_t i=0, n=cont.size(); i<n; ++i){
      auto raw = cont[i];
      auto ntu = raw->GetSizeLeading();
      auto seg = raw->SegmentId();
      for(Int_t it=0; it<ntu; ++it){
	Double_t tu = raw->GetLeading(it);
	root::HF1(Form("%s_TDC_seg%d", name, seg), tu);
      }
      if(ntu>0){
	root::HF1(Form("%s_HitPat", name), seg);
      }
    }
  }

  // BHT
  {
    static const Char_t* name = "BHT";
    const auto& cont = rawData.GetHodoRawHC(DetIdBHT);
    Int_t multi_or = 0;
    Int_t multi_and = 0;
    for(Int_t i=0, n=cont.size(); i<n; ++i){
      auto raw = cont[i];
      if(!raw) continue;
      auto seg = raw->SegmentId();
      // Up
      Bool_t u_is_good = false;
      for(Int_t j=0, m=raw->GetSizeTdcUp(); j<m; ++j){
	Double_t t = raw->GetTdcUp(j);
        if(gUser.IsInRange(Form("%s_TDC", name), t))
          u_is_good = true;
	root::HF1(Form("%s_TDC_seg%dU", name, seg), t);
      }
      for(Int_t j=0, m=raw->GetSizeAdcUp(); j<m; ++j){
	Double_t t = raw->GetAdcUp(j); // trailing
	root::HF1(Form("%s_Trailing_seg%dU", name, seg), t);
	if(m == raw->GetSizeTdcUp()){
	  Double_t l = raw->GetTdcUp(j);
          if(gUser.IsInRange(Form("%s_TDC", name), l)){
            root::HF1(Form("%s_TOT_seg%dU", name, seg), l - t);
          }
	}
      }
      // Down
      Bool_t d_is_good = false;
      for(Int_t j=0, m=raw->GetSizeTdcDown(); j<m; ++j){
	Double_t t = raw->GetTdcDown(j);
        if(gUser.IsInRange(Form("%s_TDC", name), t))
          d_is_good = true;
	root::HF1(Form("%s_TDC_seg%dD", name, seg), t);
      }
      for(Int_t j=0, m=raw->GetSizeAdcDown(); j<m; ++j){
	Double_t t = raw->GetAdcDown(j); // trailing
	root::HF1(Form("%s_Trailing_seg%dD", name, seg), t);
	if(m == raw->GetSizeTdcDown()){
	  Double_t l = raw->GetTdcDown(j);
          if(gUser.IsInRange(Form("%s_TDC", name), l)){
            root::HF1(Form("%s_TOT_seg%dD", name, seg), l - t);
          }
	}
      }
      if(u_is_good || d_is_good){
        root::HF1(Form("%s_HitPat_OR", name), seg);
        ++multi_or;
      }
      if(u_is_good && d_is_good){
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
    const auto& cont = rawData.GetHodoRawHC(DetIdAC);
    Int_t multi = 0;
    for(Int_t i=0, n=cont.size(); i<n; ++i){
      auto raw = cont[i];
      if(!raw) continue;
      Int_t seg = raw->SegmentId();
      Bool_t is_good = false;
      for(Int_t j=0, m=raw->GetSizeTdcUp(); j<m; ++j){
	Double_t t = raw->GetTdcUp(j);
        if(gUser.IsInRange(Form("%s_TDC", name), t))
          is_good = true;
	root::HF1(Form("%s_TDC_seg%d", name, seg), t);
      }
      Double_t suma = 0;
      for(Int_t j=0, m=raw->GetSizeAdcUp(); j<m; ++j){
	Double_t a = raw->GetAdcUp(j);
        suma += a;
	root::HF1(Form("%s_ADC_seg%d", name, j+1), a);
        if(is_good)
          root::HF1(Form("%s_AWT_seg%d", name, j+1), a);
      }
      // SUM
      root::HF1(Form("%s_ADC_seg%d", name, 0), suma);
      if(is_good){
        root::HF1(Form("%s_AWT_seg%d", name, 0), suma);
        root::HF1(Form("%s_HitPat", name), seg);
        multi++;
      }
    }
    root::HF1(Form("%s_Multi", name), multi);
  }

  // Hodoscope
  for(Int_t ihodo=kT1; ihodo<kNumHodo;++ihodo){
    const Char_t* name = NameHodo[ihodo];
    Int_t multi_or = 0;
    Int_t multi_and = 0;
    const auto& cont = rawData.GetHodoRawHC(DetIdHodo[ihodo]);
    for(Int_t i=0, n=cont.size(); i<n; ++i){
      auto raw = cont[i];
      if(!raw) continue;
      auto seg = raw->SegmentId();
      // Up
      Bool_t u_is_good = false;
      for(Int_t j=0, m=raw->GetSizeTdcUp(); j<m; ++j){
	Double_t t = raw->GetTdcUp(j);
        if(gUser.IsInRange(Form("%s_TDC", name), t))
          u_is_good = true;
	root::HF1(Form("%s_TDC_seg%dU", name, seg), t);
      }
      for(Int_t j=0, m=raw->GetSizeAdcUp(); j<m; ++j){
	Double_t a = raw->GetAdcUp(j);
	root::HF1(Form("%s_ADC_seg%dU", name, seg), a);
	if(u_is_good)
          root::HF1(Form("%s_AWT_seg%dU", name, seg), a);
      }
      // Down
      Bool_t d_is_good = false;
      for(Int_t j=0, m=raw->GetSizeTdcDown(); j<m; ++j){
	Double_t t = raw->GetTdcDown(j);
        if(gUser.IsInRange(Form("%s_TDC", name), t))
          d_is_good = true;
	root::HF1(Form("%s_TDC_seg%dD", name, seg), t);
      }
      for(Int_t j=0, m=raw->GetSizeAdcDown(); j<m; ++j){
	Double_t a = raw->GetAdcDown(j);
	root::HF1(Form("%s_ADC_seg%dD", name, seg), a);
	if(d_is_good)
          root::HF1(Form("%s_AWT_seg%dD", name, seg), a);
      }
      if(u_is_good || d_is_good){
        root::HF1(Form("%s_HitPat_OR", name), seg);
        ++multi_or;
      }
      if(u_is_good && d_is_good){
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
  return InitializeParameter<UserParamMan>("USER");
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
