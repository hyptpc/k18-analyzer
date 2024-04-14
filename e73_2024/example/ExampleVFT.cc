/**
 *  file: UserSkeleton.cc
 *  date: 2017.04.10
 *
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include "TString.h"

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "RootHelper.hh"
#include "RawData.hh"
#include "DCRawHit.hh"
#include "HodoRawHit.hh"
#include "VEvent.hh"
#include "HistTools.hh"
#include "UserParamMan.hh"
#include "UnpackerManager.hh"

#include "setup.hh"

#define DEBUG 0
namespace
{
  using namespace root;
  using namespace e73_2024;
  using namespace hddaq::unpacker;
  const std::string& classname("EventDC");
  const UnpackerManager& gUnpacker = GUnpacker::get_instance();
}

//______________________________________________________________________________
VEvent::VEvent( void )
{
}

//______________________________________________________________________________
VEvent::~VEvent( void )
{
}

//______________________________________________________________________________
class EventBeam : public VEvent
{
private:
  RawData      *rawData;

public:
  EventBeam( void );
  ~EventBeam( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventBeam::EventBeam( void )
  : VEvent(),
    rawData(0)
{
}

//______________________________________________________________________________
EventBeam::~EventBeam( void )
{
  if( rawData ) delete rawData;
}

//______________________________________________________________________________
typedef std::vector<std::vector<Double_t>> TDCVector;
namespace root
{
  TTree *tree;
  int run_number;
  int event_number;
  std::vector<Short_t> vft_layer;
  std::vector<Short_t> vft_channel;
  TDCVector vft_leading;
  TDCVector vft_trailing;
  double tdcbins[3]={5000,0, 2e6};
  double adcbins[3]={4096,-0.5,4095.5};
  double mtdcbins[3]={4000,0,4000};
  double dtbins[3]={3000,-500,2500};
  double dlbins[3]={2000,-5,15};
}

//______________________________________________________________________________
bool
EventBeam::ProcessingBegin( void )
{
  InitializeEvent();
  event_number=gUnpacker.get_event_number();
  return true;
}

//______________________________________________________________________________
bool
EventBeam::ProcessingNormal( void )
{
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  //  std::cout<<"event# "<<event_number<<std::endl;
  rawData = new RawData;
  rawData->DecodeHits();

  // CDH
  {
    int ihodo=kCDH;
    int cid=hodoid[ihodo];
    TString tmpname=hodoname[ihodo];
    int mul=0,mulu=0,muld=0;
    double mulbins[3]={nsegs[ihodo]+1,-0.5,nsegs[ihodo]+0.5};
    double mulbins2[3]={10,-0.5,9.5};
    double patbins[3]={nsegs[ihodo],-0.5,nsegs[ihodo]-0.5};
    const HodoRHitContainer &cont = rawData->GetHodoRawHC(cid);
    int nh = cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *raw = cont[i];
      if(!raw) continue;
      int seg = raw->SegmentId();
      double au  = raw->GetAdcUp();
      double ad  = raw->GetAdcDown();
      //       std::cout<<tmpname<<"  "<<seg<<"  "<<au<<"  "<<ad<<std::endl;
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
    int ichm=kVFT;
    TString tmpname=chmname[ichm];
    int cid=chmid[ichm];
    double nl=nlayers[ichm];
    double nw=nwires[ichm];
    double mulbins[3]={nw,-0.5,nw+0.5};
    double mulbins2[3]={10,-0.5,9.5};
    double patbins[3]={nw,-0.5,nw-0.5};
    double lpatbins[3]={nl,0.5,nl+0.5};
    for( int layer=0; layer<nl; ++layer ){
      const DCRHitContainer &RHitCont=rawData->GetDCRawHC(cid,layer);
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
	int tsize=hit->GetTrailingSize();
	std::vector<Double_t> tmp_leading;
	std::vector<Double_t> tmp_trailing;
	for(int i=0;i<mul_wire;i++){
	  int tdc=hit->GetTdc(i);
	  hist::H1(Form("%s_Leading_layer%d",tmpname.Data(),layer),tdc,mtdcbins);
	  tmp_leading.push_back(tdc);
	  if(i<tsize){
	    int trailing=hit->GetTrailing(i);
	    hist::H1(Form("%s_TOT_layer%d",tmpname.Data(),layer),tdc-trailing,mtdcbins);
	  }
	}
	for(int i=0;i<tsize;i++){
	  int tdc=hit->GetTrailing(i);
	  hist::H1(Form("%s_Trailing_layer%d",tmpname.Data(),layer),hit->GetTrailing(i),mtdcbins);
	  tmp_trailing.push_back(tdc);
	}
	vft_layer.push_back(layer);
	vft_channel.push_back(hw);
	vft_leading.push_back(tmp_leading);
	vft_trailing.push_back(tmp_trailing);
      }//ihit
    } //layer
  }// chamber analysis
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  return true;
}

//______________________________________________________________________________
bool
EventBeam::ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//______________________________________________________________________________
void
EventBeam::InitializeEvent( void )
{
  vft_layer.clear();
  vft_channel.clear();
  vft_leading.clear();
  vft_trailing.clear();
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventBeam;
}
//______________________________________________________________________________
bool
ConfMan::InitializeHistograms( void )
{
  tree=new TTree("tree","vft raw data");
  // tree->Branch("runnum", &run_number);
  // tree->Branch("evnum", &event_number);
  // tree->Branch("vft_layer"  ,&vft_layer);
  // tree->Branch("vft_channel",&vft_channel);
  // tree->Branch("vft_leading",&vft_leading);
  // tree->Branch("vft_trailing",&vft_trailing);
  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<UserParamMan>("USER") );
}

//______________________________________________________________________________
bool
ConfMan::BeginRunProcess()
{
  return true;
}
//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
