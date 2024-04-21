// -*- C++ -*-

#include <iostream>
#include <sstream>
#include <cmath>
#include "TString.h"

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RootHelper.hh"
#include "MTDCRawHit.hh"
#include "HodoRawHit.hh"
#include "DCRawHit.hh"
#include "RawData.hh"
#include "CDCWireMapMan.hh"
#include "VEvent.hh"
#include "HistTools.hh"

#include "UnpackerManager.hh"
#include "DAQNode.hh"

#include "setup.hh"

#define DEBUG 0
namespace
{
using namespace e73_2024;
using namespace root;
using namespace hddaq::unpacker;
using namespace hddaq;
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
    int cid=DetIdTrigFlag;
    double patbins[3]={16,-0.5,16-0.5};
    const MTDCRHitContainer &cont = rawData.GetMTDCRawHC(cid);
    int nh = cont.size();
    for( int i=0; i<nh; ++i ){
      MTDCRawHit *raw = cont[i];
      int ntu = raw->GetSizeLeading();
      int seg = raw->SegmentId();
      TString tmpname="TriggerFlag";
      for(int it=0;it<ntu;it++){
	double tu  = raw->GetLeading(it);
	hist::H1(Form("%s_TDC_%d",tmpname.Data(),seg),tu,mtdcbins);
      }
      if(ntu>0) {
	hist::H1(Form("%s_Pat",tmpname.Data()),seg,patbins);
      }
    }
  }

  // AC
  {
    int cid=DetIdAC;
    const HodoRHitContainer &cont = rawData.GetHodoRawHC(cid);
    int nh = cont.size();
    TString tmpname="AC";
    for( int i=0; i<nh; ++i ){
      HodoRawHit *raw = cont[i];
      if(!raw) continue;
      int seg = raw->SegmentId();
      int nt=raw->GetSizeTdcUp();
      int na=raw->GetSizeAdcUp();
      for(int it=0;it<nt;it++){
	double tu  = raw->GetTdcUp(it);
	hist::H1(Form("%s_TDC_seg%d",tmpname.Data(),seg),tu,tdcbins);
      }
      for(int it=0;it<na;it++){
	double tu  = raw->GetAdcUp(it);
	hist::H1(Form("%s_ADC_seg%d_ch%d",tmpname.Data(),seg,it),tu,adcbins);
      }
    }
  }

  // BHT
  {
    int ihodo=kBHT;
    int cid=hodoid[ihodo];
    TString tmpname=hodoname[ihodo];
    int mul=0,mulu=0,muld=0;
    double mulbins[3]={nsegs[ihodo]+1,-0.5,nsegs[ihodo]+0.5};
    double mulbins2[3]={10,-0.5,9.5};
    double patbins[3]={nsegs[ihodo],-0.5,nsegs[ihodo]-0.5};
    const HodoRHitContainer &cont = rawData.GetHodoRawHC(cid);
    int nh = cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *raw = cont[i];
      if(!raw) continue;
      int seg = raw->SegmentId();
      int ntu=raw->GetSizeTdcUp();
      int ntd=raw->GetSizeTdcDown();
      hist::H1(Form("%s_Mulu_seg%d",tmpname.Data(),seg),ntu,mulbins2);
      hist::H1(Form("%s_Muld_seg%d",tmpname.Data(),seg),ntd,mulbins2);
      hist::H2(Form("%s_Mulu_Muld_seg%d",tmpname.Data(),seg),ntu,ntd,mulbins2,mulbins2);
      hist::H2(Form("%s_Mulu_pat",tmpname.Data()),seg,ntu,patbins,mulbins2);
      hist::H2(Form("%s_Muld_pat",tmpname.Data()),seg,ntd,patbins,mulbins2);
      if(ntu>0){
	hist::H1(Form("%s_Patu",tmpname.Data()),seg,patbins);
	mulu++;
      }
      if(ntd>0){
	hist::H1(Form("%s_Patd",tmpname.Data()),seg,patbins);
	muld++;
      }
      for(int it=0;it<ntu;it++){
	double tu  = raw->GetTdcUp(it);
	hist::H1(Form("%s_TDCu_seg%d",tmpname.Data(),seg),tu,tdcbins);
      }
      for( int it=0; it<raw->GetSizeAdcUp(); ++it ){
	double tu  = raw->GetAdcUp(it);
	hist::H1(Form("%s_Trailingu_seg%d",tmpname.Data(),seg),tu,tdcbins);
	if(it<ntu){
	  double tl=raw->GetTdcUp(it);
	  double tot=tl-tu;
	  //	  std::cout<<raw->GetTdcUp(it)<<"  "<<tu<<"  "<<tot<<std::endl;
	  hist::H1(Form("%s_TOTu_seg%d",tmpname.Data(),seg),tot,totbins);
	  hist::H2(Form("%s_TOTvsTDCu_seg%d",tmpname.Data(),seg),tl,tot,200,1e6,2e6,200,0,2e4);
	}
      }
      for(int it=0;it<ntd;it++){
	double td  = raw->GetTdcDown(it);
	hist::H1(Form("%s_TDCd_seg%d",tmpname.Data(),seg),td,tdcbins);
      }
      for( int it=0; it<raw->GetSizeAdcDown(); ++it ){
	double td  = raw->GetAdcDown(it);
	hist::H1(Form("%s_Trailingd_seg%d",tmpname.Data(),seg),td,tdcbins);
	if(it<ntd){
	  double tl=raw->GetTdcDown(it);
	  double tot=tl-td;
	  hist::H1(Form("%s_TOTd_seg%d",tmpname.Data(),seg),tot,totbins);
	  hist::H2(Form("%s_TOTvsTDCd_seg%d",tmpname.Data(),seg),tl,tot,200,1e6,2e6,200,0,2e4);
	}
      }
    }//for(i)
    hist::H1(Form("%s_Mulu",tmpname.Data()),mulu,mulbins);
    hist::H1(Form("%s_Muld",tmpname.Data()),muld,mulbins);
  }

  // hodoscopes
  for(int ihodo=0;ihodo<kNumHodo;++ihodo){
    int cid=hodoid[ihodo];
    if(cid==DetIdBHT) continue;
    TString tmpname=hodoname[ihodo];
    int mul=0,mulu=0,muld=0;
    double mulbins[3]={nsegs[ihodo]+1,-0.5,nsegs[ihodo]+0.5};
    double mulbins2[3]={10,-0.5,9.5};
    double patbins[3]={nsegs[ihodo],-0.5,nsegs[ihodo]-0.5};
    const auto& cont = rawData.GetHodoRawHC(cid);
    int nh = cont.size();
    for( int i=0; i<nh; ++i ){
      auto raw = cont[i];
      if(!raw) continue;
      int seg = raw->SegmentId();
      double au  = raw->GetAdcUp();
      double ad  = raw->GetAdcDown();
      // std::cout<<tmpname<<"  "<<seg<<"  "<<au<<"  "<<ad<<std::endl;
      hist::H1(Form("%s_ADCu_seg%d",tmpname.Data(),seg),au,adcbins);
      hist::H1(Form("%s_ADCd_seg%d",tmpname.Data(),seg),ad,adcbins);
      int ntu=raw->GetSizeTdcUp();
      int ntd=raw->GetSizeTdcDown();
      hist::H1(Form("%s_Mulu_seg%d",tmpname.Data(),seg),ntu,mulbins2);
      hist::H1(Form("%s_Muld_seg%d",tmpname.Data(),seg),ntd,mulbins2);
      hist::H2(Form("%s_Mulu_Muld_seg%d",tmpname.Data(),seg),ntu,ntd,mulbins2,mulbins2);
      hist::H2(Form("%s_Mulu_pat",tmpname.Data()),seg,ntu,patbins,mulbins2);
      hist::H2(Form("%s_Muld_pat",tmpname.Data()),seg,ntd,patbins,mulbins2);
      if(ntu>0){
	hist::H1(Form("%s_Patu",tmpname.Data()),seg,patbins);
	hist::H1(Form("%s_ADCwTu_seg%d",tmpname.Data(),seg),au,adcbins);
	mulu++;
      }
      if(ntd>0){
	hist::H1(Form("%s_Patd",tmpname.Data()),seg,patbins);
	hist::H1(Form("%s_ADCwTd_seg%d",tmpname.Data(),seg),ad,adcbins);
	muld++;
      }
      for(int it=0;it<ntu;it++){
	double tu  = raw->GetTdcUp(it);
	hist::H1(Form("%s_TDCu_seg%d",tmpname.Data(),seg),tu,tdcbins);
      }
      for(int it=0;it<ntd;it++){
	double td  = raw->GetTdcDown(it);
	hist::H1(Form("%s_TDCd_seg%d",tmpname.Data(),seg),td,tdcbins);
      }
    }//for(i)
    hist::H1(Form("%s_Mulu",tmpname.Data()),mulu,mulbins);
    hist::H1(Form("%s_Muld",tmpname.Data()),muld,mulbins);
  }

  // Chamber ------------------------------------------------------------
  {
    for(int ichm=0;ichm<kCDC;ichm++){
      TString tmpname=chmname[ichm];
      int cid=chmid[ichm];
      double nl=nlayers[ichm];
      double nw=nwires[ichm];
      double mulbins[3]={nw,-0.5,nw+0.5};
      double mulbins2[3]={10,-0.5,9.5};
      double patbins[3]={nw,-0.5,nw-0.5};
      double lpatbins[3]={nl,0.5,nl+0.5};
      for( int layer=0; layer<nl; ++layer ){
	const auto& RHitCont = rawData.GetDCRawHC(cid, layer);
	int nh = RHitCont.size();
	hist::H1(Form("%s_Mul_layer%d",tmpname.Data(),layer),nh,mulbins);
	for( int i=0; i<nh; ++i ){
	  DCRawHit *hit  = RHitCont[i];
	  int hw      = hit->WireId();
	  int mul_wire= hit->GetTdcSize();
	  hist::H1(Form("%s_HitPat_layer%d",tmpname.Data(),layer),hw,patbins);
	  hist::H2(Form("%s_HitPat_2d",tmpname.Data()),layer,hw,lpatbins,patbins);
	  hist::H1(Form("%s_Mul_layer%d_wire%d",tmpname.Data(),layer,hw),mul_wire,mulbins);
	  hist::H2(Form("%s_Mul_layer%d_2d",tmpname.Data(),layer),hw,mul_wire,patbins,mulbins2);
	  for(int i=0;i<mul_wire;i++){
	    hist::H1(Form("%s_Leading_layer%d",tmpname.Data(),layer),hit->GetTdc(i),mtdcbins);
	  }
	  // for(int i=0;i<hit->GetTrailingSize();i++){
	  //   hist::H1(Form("%s_Trailing_layer%d",tmpname.Data(),layer),hit->GetTrailing(i),mtdcbins);
	  // }
	}//ihit
      } //layer
    } // ichm
  }// chamber analysis
  // DAQ
  static UnpackerManager& gUnpacker = GUnpacker::get_instance();
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
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return ( InitializeParameter<CDCWireMapMan>("CDCGeom", "CDCASD" ) );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
