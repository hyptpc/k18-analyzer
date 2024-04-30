// -*- C++ -*-

#include <iostream>
#include <sstream>
#include <cmath>

#include <TString.h>

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


#include "setup.hh"

#define DEBUG 0
namespace
{
using namespace e73_2024;
using namespace root;
using namespace hddaq::unpacker;
using namespace hddaq;
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
namespace root
{
TH1   *h[MaxHist];
double tdcbins[3]={5000,0,2e6};
double totbins[3]={5000,0,5e4};
double adcbins[3]={4096,-0.5,4095.5};
double mtdcbins[3]={2000,0,2000};
double diffbins[3]={2000,-10000,10000};
}

//_____________________________________________________________________________
void
InitializeEvent()
{
}

//_____________________________________________________________________________
Bool_t
ProcessBegin()
{
  InitializeEvent();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessNormal()
{
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
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
      Bool_t u_in_range = false;
      for(Int_t j=0, m=raw->GetSizeTdcUp(); j<m; ++j){
	Double_t t = raw->GetTdcUp(j);
        if(gUser.IsInRange(Form("%s_TDC", name), t))
          u_in_range = true;
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
      Bool_t d_in_range = false;
      for(Int_t j=0, m=raw->GetSizeTdcDown(); j<m; ++j){
	Double_t t = raw->GetTdcDown(j);
        if(gUser.IsInRange(Form("%s_TDC", name), t))
          d_in_range = true;
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
      if(u_in_range || d_in_range){
        root::HF1(Form("%s_HitPat_OR", name), seg);
        ++multi_or;
      }
      if(u_in_range && d_in_range){
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
      Bool_t in_range = false;
      for(int j=0, m=raw->GetSizeTdcUp(); j<m; ++j){
	Double_t t = raw->GetTdcUp(j);
        if(gUser.IsInRange(Form("%s_TDC", name), t))
          in_range = true;
	root::HF1(Form("%s_TDC_seg%d", name, seg), t);
      }
      for(int j=0, m=raw->GetSizeAdcUp(); j<m; ++j){
	Double_t a = raw->GetAdcUp(j);
	root::HF1(Form("%s_ADC_seg%d", name, seg), a);
        if(in_range)
          root::HF1(Form("%s_AWT_seg%d", name, seg), a);
      }
      if(in_range){
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
      Bool_t u_in_range = false;
      for(Int_t j=0, m=raw->GetSizeTdcUp(); j<m; ++j){
	Double_t t = raw->GetTdcUp(j);
        if(gUser.IsInRange(Form("%s_TDC", name), t))
          u_in_range = true;
	root::HF1(Form("%s_TDC_seg%dU", name, seg), t);
      }
      for(Int_t j=0, m=raw->GetSizeAdcUp(); j<m; ++j){
	Double_t a = raw->GetAdcUp(j);
	root::HF1(Form("%s_ADC_seg%dU", name, seg), a);
	if(u_in_range)
          root::HF1(Form("%s_AWT_seg%dU", name, seg), a);
      }
      // Down
      Bool_t d_in_range = false;
      for(Int_t j=0, m=raw->GetSizeTdcDown(); j<m; ++j){
	Double_t t = raw->GetTdcDown(j);
        if(gUser.IsInRange(Form("%s_TDC", name), t))
          d_in_range = true;
	root::HF1(Form("%s_TDC_seg%dD", name, seg), t);
      }
      for(Int_t j=0, m=raw->GetSizeAdcDown(); j<m; ++j){
	Double_t a = raw->GetAdcDown(j);
	root::HF1(Form("%s_ADC_seg%dD", name, seg), a);
	if(d_in_range)
          root::HF1(Form("%s_AWT_seg%dD", name, seg), a);
      }
      if(u_in_range || d_in_range){
        root::HF1(Form("%s_HitPat_OR", name), seg);
        ++multi_or;
      }
      if(u_in_range && d_in_range){
        root::HF1(Form("%s_HitPat_AND", name), seg);
        ++multi_and;
      }
    }
    root::HF1(Form("%s_Multi_OR", name), multi_or);
    root::HF1(Form("%s_Multi_AND", name), multi_and);
  }

  // DC
  for(Int_t idc=0; idc<kCDC; ++idc){
    const Char_t* name = NameDC[idc].Data();
    Int_t nlayer = NumOfLayerDC[idc];
    for(Int_t layer=0; layer<nlayer; ++layer){
      const auto& cont = rawData.GetDCRawHC(DetIdDC[idc], layer);
      // hist::H1(Form("%s_Mul_layer%d",name,layer),nh,mulbins);
      Int_t multi = 0;
      for(Int_t i=0, n=cont.size(); i<n; ++i){
        auto hit = cont[i];
        if(!hit) continue;
        Int_t wire = hit->WireId();
        Bool_t in_range = false;
        for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
          Double_t t = hit->GetTdc(j);
          root::HF1(Form("%s_TDC_layer%d", name, layer), t);
          in_range = true;
          if(j == 0)
            root::HF1(Form("%s_TDC1st_layer%d", name, layer), t);
        }
        for(Int_t j=0, m=hit->GetTrailingSize(); j<m; ++j){
          Double_t t = hit->GetTrailing(j);
          root::HF1(Form("%s_Trailing_layer%d", name, layer), t);
          if(j == 0)
            root::HF1(Form("%s_Trailing1st_layer%d", name, layer), t);
          if(m == hit->GetTdcSize()){
            Double_t l = hit->GetTdc(j);
            Double_t tot = l - t;
            root::HF1(Form("%s_TOT_layer%d", name, layer), tot);
            if(j == 0)
              root::HF1(Form("%s_TOT1st_layer%d", name, layer), tot);
          }
        }
        if(in_range){
          root::HF1(Form("%s_HitPat_layer%d", name, layer), wire);
          ++multi;
        }
      }
      root::HF1(Form("%s_Multi_layer%d", name, layer), multi);
    }
  }

  // DAQ
  int evnum=gUnpacker.get_event_number();
  static const int k_eb      = gUnpacker.get_fe_id("k18breb");
  static const int k_vme     = gUnpacker.get_fe_id("vme_qdc1");
  static const int k_sca     = gUnpacker.get_fe_id("hulscaler-139");
  static const int k_flag    = gUnpacker.get_fe_id("hulmhtdc-138");
  static const int k_hr      = gUnpacker.get_fe_id("hulhrtdc-121");
  static const int k_cdc1    = 1650;
  static const int k_cdc2    = 1660;
  static const int k_bldc    = gUnpacker.get_fe_id("hulmhtdc-101");
  {
    int data_size = gUnpacker.get_node_header( k_eb, DAQNode::k_data_size);
    hist::H1("DataSize; Words; Counts",data_size,10000,0,20000);
    hist::H2("DataSize_vs_EvNum",evnum,data_size,300,0,3e6,200,0,20000);
  }
  { // VME node
    for( int i=0; i<2; ++i ){
      int node_id = k_vme+i;
      int data_size = gUnpacker.get_node_header( node_id, DAQNode::k_data_size);
      hist::H1(Form("DataSize_node%04d",node_id),data_size,1000,0,1000);
      hist::H2(Form("DataSize_node%04d_vs_EvNum",node_id),
	       evnum,data_size,300,0,3e6,100,0,2000);
      hist::H2(Form("DataSize_vme"),
	       i,data_size,2,-0.5,1.5,1000,0,2000);
    }
  }

  { // HUL
    for( int i=0; i<50; ++i ){
      int node_id = k_bldc+i;
      int data_size = gUnpacker.get_node_header( node_id, DAQNode::k_data_size);
      hist::H1(Form("DataSize_node%04d",node_id),data_size,1000,0,1000);
      hist::H2(Form("DataSize_node%04d_vs_EvNum",node_id),
	       evnum,data_size,300,0,3e6,100,0,2000);
      hist::H2(Form("DataSize_hul"),
	       i,data_size,50,-0.5,49.5,1000,0,2000);
    }
  }
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
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
  { // TriggerFlag
    const Char_t* name = "TriggerFlag";
    Double_t patbins[3]={NumOfSegTrigFlag, -0.5, NumOfSegTrigFlag-0.5};
    for(Int_t i=0; i<NumOfSegTrigFlag; ++i){
      root::HB1(Form("%s_TDC_seg%d", name, i), mtdcbins);
    }
    root::HB1(Form("%s_HitPat; Segment; Counts", name), patbins);
  }

  { // BHT
    const Char_t* name = "BHT";
    Int_t nseg = NumOfSegBHT;
    for(Int_t i=0; i<nseg; ++i){
      for(const auto& ud : std::vector<TString>{"U", "D"}){
        // root::HB1(Form("%s_ADC_seg%d%s", name, i, ud.Data()), adcbins);
        // root::HB1(Form("%s_AWT_seg%d%s", name, i, ud.Data()), adcbins);
        root::HB1(Form("%s_TDC_seg%d%s", name, i, ud.Data()), tdcbins);
        root::HB1(Form("%s_Trailing_seg%d%s", name, i, ud.Data()), tdcbins);
        root::HB1(Form("%s_TOT_seg%d%s", name, i, ud.Data()), totbins);
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
      root::HB1(Form("%s_TDC_seg%d", name, i), tdcbins);
    }
    root::HB1(Form("%s_HitPat", name), nseg, -0.5, nseg - 0.5);
    root::HB1(Form("%s_Multi", name), nseg + 1, -0.5, nseg + 0.5);
  }

  // Hodoscope
  for(Int_t ihodo=kT1; ihodo<kNumHodo;++ihodo){
    auto name = NameHodo[ihodo].Data();
    Int_t nseg = NumOfSegHodo[ihodo];
    for(const auto& uord: std::vector<TString>{"U", "D"} ){
      auto ud = uord.Data();
      for(Int_t i=0; i<nseg; ++i){
        root::HB1(Form("%s_ADC_seg%d%s", name, i, ud), adcbins);
        root::HB1(Form("%s_AWT_seg%d%s", name, i, ud), adcbins);
        root::HB1(Form("%s_TDC_seg%d%s", name, i, ud), tdcbins);
      }
    }
    for(const auto& uord: std::vector<TString>{"OR", "AND"} ){
      auto ud = uord.Data();
      root::HB1(Form("%s_HitPat_%s", name, ud), nseg, -0.5, nseg - 0.5);
      root::HB1(Form("%s_Multi_%s", name, ud), nseg + 1, -0.5, nseg + 0.5);
    }
  }

  // DC
  for(Int_t idc=0; idc<kCDC; ++idc){
    const Char_t* name = NameDC[idc].Data();
    Int_t nlayer = NumOfLayerDC[idc];
    Int_t nwire = NumOfWireDC[idc];
    // double mulbins[3]={nw,-0.5,nw+0.5};
    // double mulbins2[3]={10,-0.5,9.5};
    // double patbins[3]={nw,-0.5,nw-0.5};
    // double lpatbins[3]={nlayer,0.5,nlayer+0.5};
    for(Int_t layer=0; layer<nlayer; ++layer){
      root::HB1(Form("%s_TDC_layer%d", name, layer), mtdcbins);
      root::HB1(Form("%s_TDC1st_layer%d", name, layer), mtdcbins);
      root::HB1(Form("%s_Trailing_layer%d", name, layer), mtdcbins);
      root::HB1(Form("%s_Trailing1st_layer%d", name, layer), mtdcbins);
      root::HB1(Form("%s_TOT_layer%d", name, layer), mtdcbins);
      root::HB1(Form("%s_TOT1st_layer%d", name, layer), mtdcbins);
      root::HB1(Form("%s_HitPat_layer%d", name, layer), nwire, -0.5, nwire - 0.5);
      root::HB1(Form("%s_Multi_layer%d", name, layer), nwire + 1, -0.5, nwire + 0.5);
      // for(Int_t wire=0; wire<nwire; ++wire){
      // }
    }
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
