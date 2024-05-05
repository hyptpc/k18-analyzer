// -*- C++ -*-

#include "EventAnalyzer.hh"

#include <Unpacker.hh>
#include <UnpackerManager.hh>
#include <DAQNode.hh>

#include "DCRawHit.hh"
#include "DetectorID.hh"
#include "HodoRawHit.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"

namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
using hist::EParticle;
using hist::particle_list;
}

//_____________________________________________________________________________
EventAnalyzer::EventAnalyzer()
{
}

//_____________________________________________________________________________
EventAnalyzer::~EventAnalyzer()
{
}

//_____________________________________________________________________________
void
EventAnalyzer::HodoRawHit(const RawData& rawData, EParticle particle)
{
  const Char_t* p = particle_list.at(particle).Data();
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
          root::HF1(Form("%s_TDC_seg%d%s%s", name, seg, UorD[ud], p), t);
        }
        for(Int_t i=0, n=hit->GetSizeTdcTrailing(ud); i<n; ++i){
          Double_t t = hit->GetTdcTrailing(ud, i);
          root::HF1(Form("%s_Trailing_seg%d%s%s", name, seg, UorD[ud], p), t);
          if(n == hit->GetSizeTdcLeading(ud)){
            Double_t l = hit->GetTdc(ud, i);
            if(gUser.IsInRange(Form("%s_TDC", name), l)){
              root::HF1(Form("%s_TOT_seg%d%s%s", name, seg, UorD[ud], p), l - t);
            }
          }
        }
      }
      if(ud_good >= 1){
        root::HF1(Form("%s_HitPat_OR%s", name, p), seg);
        ++multi_or;
      }
      if(ud_good == 2){
        root::HF1(Form("%s_HitPat_AND%s", name, p), seg);
        ++multi_and;
      }
    }
    root::HF1(Form("%s_Multi_OR%s", name, p), multi_or);
    root::HF1(Form("%s_Multi_AND%s", name, p), multi_and);
  }

  // AC
  {
    static const Char_t* name = "AC";
    Int_t multi = 0;
    Bool_t is_good = false;
    for(const auto& hit: rawData.GetHodoRawHC(name)){
      Int_t seg = hit->SegmentId();
      for(const auto& t: hit->GetArrayTdcLeading()){
        if(gUser.IsInRange(Form("%s_TDC", name), t))
          is_good = true;
	root::HF1(Form("%s_TDC_seg%d%s", name, seg, p), t);
      }
      if(is_good && seg == 0){
        root::HF1(Form("%s_HitPat%s", name, p), seg);
        ++multi;
      }
      Double_t a = hit->GetAdc();
      root::HF1(Form("%s_ADC_seg%d%s", name, seg, p), a);
      if(is_good)
        root::HF1(Form("%s_AWT_seg%d%s", name, seg, p), a);
    }
    root::HF1(Form("%s_Multi%s", name, p), multi);
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
          root::HF1(Form("%s_TDC_seg%d%s%s", name, seg, UorD[ud], p), t);
        }
        for(const auto& a: hit->GetArrayAdc(ud)){
          root::HF1(Form("%s_ADC_seg%d%s%s", name, seg, UorD[ud], p), a);
          if(is_good)
            root::HF1(Form("%s_AWT_seg%d%s%s", name, seg, UorD[ud], p), a);
        }
        ud_good += is_good;
      }
      if(ud_good >= 1){
        root::HF1(Form("%s_HitPat_OR%s", name, p), seg);
        ++multi_or;
      }
      if(ud_good == 2){
        root::HF1(Form("%s_HitPat_AND%s", name, p), seg);
        ++multi_and;
      }
    }
    root::HF1(Form("%s_Multi_OR%s", name, p), multi_or);
    root::HF1(Form("%s_Multi_AND%s", name, p), multi_and);
  }
}

//_____________________________________________________________________________
void
EventAnalyzer::DCRawHit(const RawData& rawData, EParticle particle)
{
  const Char_t* p = particle_list.at(particle).Data();
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
            root::HF1(Form("%s_%sTDC_layer%d%s", name, c, layer, p), l);
            root::HF1(Form("%s_%sTrailing_layer%d%s", name, c, layer, p), t);
            root::HF1(Form("%s_%sTOT_layer%d%s", name, c, layer, p), tot);
            if(j == 0){
              root::HF1(Form("%s_%sTDC1st_layer%d%s", name, c, layer, p), l);
              root::HF1(Form("%s_%sTrailing1st_layer%d%s", name, c, layer, p), t);
              root::HF1(Form("%s_%sTOT1st_layer%d%s", name, c, layer, p), tot);
            }
          }
        }
        if(hit->GetTdcSize() > 0){
          root::HF1(Form("%s_HitPat_layer%d%s", name, layer, p), wire);
          ++multi;
        }
        if(is_good){
          root::HF1(Form("%s_CHitPat_layer%d%s", name, layer, p), wire);
          ++cmulti;
        }
      }
      root::HF1(Form("%s_Multi_layer%d%s", name, layer, p), multi);
      root::HF1(Form("%s_CMulti_layer%d%s", name, layer, p), cmulti);
    }
  }
}

//_____________________________________________________________________________
void
EventAnalyzer::TriggerFlag(const RawData& rawData)
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

//_____________________________________________________________________________
void
EventAnalyzer::DAQ(const RawData& rawData)
{
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
  auto node_id = gUnpacker.get_fe_id("k18breb");
  auto data_size = gUnpacker.get_node_header(node_id, k_data_size);
  root::HF1("EB_DataSize", data_size);
}
