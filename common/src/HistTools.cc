// -*- C++ -*-

#include "HistTools.hh"

#include <TString.h>

#include <Unpacker.hh>
#include <UnpackerManager.hh>
#include <DAQNode.hh>

#include "DetectorID.hh"
#include "RootHelper.hh"

namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
}

namespace hist
{
const Double_t hrtdcbins1[3] = {45000, 950000, 1400000};
const Double_t hrtdcbins2[3] = {50000, 0, 1000000}; // for CVC, NC
const Double_t hrtotbins[3] = {5000, 0, 50000};
const Double_t adcbins[3] = {4096, -0.5, 4095.5};
const Double_t mhtdcbins[3] = {2000, 0, 2000};
const Double_t mhtotbins[3] = {1000, 0, 1000};

void
BuildStatus()
{
  root::HB1("Status", 21, -0.5, 20.5);
}

void
BuildTriggerFlag()
{
  const Char_t* name = "TriggerFlag";
  Double_t patbins[3]={NumOfSegTrigFlag, -0.5, NumOfSegTrigFlag-0.5};
  for(Int_t i=0; i<NumOfSegTrigFlag; ++i){
    root::HB1(Form("%s_TDC_seg%d", name, i), mhtdcbins);
  }
  root::HB1(Form("%s_HitPat; Segment; Counts", name), patbins);
}

void
BuildHodoRaw(Bool_t flag_particle)
{
  for(const auto& particle: std::vector<TString>{"", "_Pi", "_K", "_P"}){
    const Char_t* p = particle.Data();
    { // BHT
      const Char_t* name = "BHT";
      Int_t nseg = NumOfSegBHT;
      for(Int_t i=0; i<nseg; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          root::HB1(Form("%s_TDC_seg%d%s%s", name, i, ud, p), hrtdcbins1);
          root::HB1(Form("%s_Trailing_seg%d%s%s", name, i, ud, p), hrtdcbins1);
          root::HB1(Form("%s_TOT_seg%d%s%s", name, i, ud, p), hrtotbins);
        }
      }
      for(const auto& uord: std::vector<TString>{"OR", "AND"} ){
        auto ud = uord.Data();
        root::HB1(Form("%s_HitPat_%s%s", name, ud, p), nseg, -0.5, nseg - 0.5);
        root::HB1(Form("%s_Multi_%s%s", name, ud, p), nseg + 1, -0.5, nseg + 0.5);
      }
    }
    { // AC
      const Char_t* name = "AC";
      Int_t nseg = NumOfSegAC;
      for(Int_t i=0; i<nseg; ++i){
        root::HB1(Form("%s_ADC_seg%d%s", name, i, p), adcbins);
        root::HB1(Form("%s_AWT_seg%d%s", name, i, p), adcbins);
        root::HB1(Form("%s_TDC_seg%d%s", name, i, p), mhtdcbins);
      }
      root::HB1(Form("%s_HitPat%s", name, p), nseg, -0.5, nseg - 0.5);
      root::HB1(Form("%s_Multi%s", name, p), nseg + 1, -0.5, nseg + 0.5);
    }
    // Hodoscope
    for(Int_t ihodo=kT1; ihodo<kNumHodo;++ihodo){
      auto name = NameHodo[ihodo].Data();
      const Double_t* hrtdcbins;
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
          root::HB1(Form("%s_ADC_seg%d%s%s", name, i, ud, p), adcbins);
          root::HB1(Form("%s_AWT_seg%d%s%s", name, i, ud, p), adcbins);
          root::HB1(Form("%s_TDC_seg%d%s%s", name, i, ud, p), hrtdcbins);
        }
      }
      for(const auto& uord: std::vector<TString>{"OR", "AND"} ){
        auto ud = uord.Data();
        root::HB1(Form("%s_HitPat_%s%s", name, ud, p), nseg, -0.5, nseg - 0.5);
        root::HB1(Form("%s_Multi_%s%s", name, ud, p), nseg + 1, -0.5, nseg + 0.5);
      }
    }
    if(!flag_particle) break;
  }
}

void
BuildDCRaw(Bool_t flag_particle)
{
  for(const auto& particle: std::vector<TString>{"", "_Pi", "_K", "_P"}){
    const Char_t* p = particle.Data();
    for(Int_t idc=0; idc<=kVFT; ++idc){
      const Char_t* name = NameDC[idc].Data();
      Int_t nlayer = NumOfLayerDC[idc];
      Double_t nwire = NumOfWireDC[idc];
      const Double_t patbins[3] = {nwire, -0.5, nwire - 0.5};
      const Double_t mulbins[3] = {nwire + 1, -0.5, nwire + 0.5};
      for(Int_t layer=0; layer<nlayer; ++layer){
        for(const auto& totcut: std::vector<TString>{"", "C"}){
          auto c = totcut.Data();
          root::HB1(Form("%s_%sTDC_layer%d%s", name, c, layer, p), mhtdcbins);
          root::HB1(Form("%s_%sTDC1st_layer%d%s", name, c, layer, p), mhtdcbins);
          root::HB1(Form("%s_%sTrailing_layer%d%s", name, c, layer, p), mhtdcbins);
          root::HB1(Form("%s_%sTrailing1st_layer%d%s", name, c, layer, p), mhtdcbins);
          root::HB1(Form("%s_%sTOT_layer%d%s", name, c, layer, p), mhtotbins);
          root::HB1(Form("%s_%sTOT1st_layer%d%s", name, c, layer, p), mhtotbins);
          root::HB1(Form("%s_%sHitPat_layer%d%s", name, c, layer, p), patbins);
          root::HB1(Form("%s_%sMulti_layer%d%s", name, c, layer, p), mulbins);
        }
      }
    }
    if(!flag_particle) break;
  }
}

void
BuildDAQ()
{
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
}
}
