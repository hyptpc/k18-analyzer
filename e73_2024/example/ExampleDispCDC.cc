/**
 *  file: UserSkeleton.cc
 *  date: 2017.04.10
 *
 */
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include "TString.h"

#include <TRint.h>
#include <TSystem.h>
#include <TCanvas.h>

#include "ConfMan.hh"
#include "DetectorID.hh"
//#include "RMAnalyzer.hh"
#include "RootHelper.hh"
#include "DCAnalyzer.hh"
#include "CDCAnalyzer.hh"
#include "HodoAnalyzer.hh"
#include "HodoParamMan.hh"
#include "CDCWireMapMan.hh"
#include "DetectorList.hh"
#include "HodoPHCMan.hh"
#include "MTDCAnalyzer.hh"
#include "MTDCRawHit.hh"
#include "UserParamMan.hh"
#include "RawData.hh"
#include "UnpackerManager.hh"
#include "BeamSpectrometer.hh"
#include "VEvent.hh"
#include "HistTools.hh"
#include "Display.hh"
#include "SimTools.hh"

#define DEBUG 0
namespace
{
  using namespace root;
  using namespace hddaq::unpacker;
  const UnpackerManager& gUnpacker = GUnpacker::get_instance();
  const std::string& classname("EventDC");
  TCanvas *c1;
  TH1 *h1;
  int event_number;  
  TString pdfname;
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
class EventRaw : public VEvent
{
private:
  RawData      *rawData;
  CDCAnalyzer  *CDCAna;
  HodoAnalyzer *hodoAna;
  MTDCAnalyzer *MTDCAna;
  
public:
  EventRaw( void );
  ~EventRaw( void );
  bool  ProcessingBegin( void );
  bool  ProcessingBeginMC( DetectorData* det, MCData *mc, ReactionData* reac,const int &iev ) override;
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventRaw::EventRaw( void )
  : VEvent(),
    rawData(0),
    CDCAna( new CDCAnalyzer ),
    hodoAna(new HodoAnalyzer),  
    MTDCAna( new MTDCAnalyzer )
{
}

//______________________________________________________________________________
EventRaw::~EventRaw( void )
{
  if( hodoAna ) delete hodoAna;
  if( CDCAna )  delete CDCAna;
  if( MTDCAna ) delete MTDCAna;
  if( rawData ) delete rawData;
}

//______________________________________________________________________________
namespace root
{
  TH1   *h[MaxHist];
  const int nhodo=5;
  int Cid[nhodo]={DetIdT0new,DetIdBHD,DetIdT0,DetIdDEF,DetIdCDH};
  TString name[nhodo]={"T0new","BHD","T0","DEF","CDH"};
  double nsegs[nhodo]={1,16,5,5,36};
  enum param_hodo{kT0new,kBHD,kT0,kDEF,kCDH};
   
  // const int nchm=7;
  // const int nchm2=5;
  // int ChmCid[nchm]={DetIdBLC1a,DetIdBLC1b,DetIdBLC2a,DetIdBLC2b,DetIdBPC,DetIdBLC1,DetIdBLC2};
  // TString chmname[nchm]={"BLC1a", "BLC1b", "BLC2a", "BLC2b","BPC","BLC1","BLC2"}; 
  // int nlayers[nchm]={8,8,8,8,8,16,16};
  // int nwires[nchm]={32,32,32,32,32,32,32};
  // enum param_chm{kBLC1a,kBLC1b,kBLC2a,kBLC2b,kBPC,kBLC1,kBLC2};

}

//______________________________________________________________________________
bool
EventRaw::ProcessingBegin( void )
{
  InitializeEvent();
  rawData = new RawData;
  rawData->DecodeHits();
  hodoAna->DecodeRawHits( rawData );
  MTDCAna->DecodeRawHits( rawData );
  CDCAna ->DecodeRawHits( rawData );
  return true;
}
//______________________________________________________________________________
bool
EventRaw::ProcessingBeginMC( DetectorData *det, MCData* mc, ReactionData* reac, const int &iev )
{
  event_number=iev;
  InitializeEvent();
  sim::Convert(det,CDCAna,0,hodoAna);    
  return true;
}
//______________________________________________________________________________
bool
EventRaw::ProcessingNormal( void )
{
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  bool BEAM=false, KAON3=false;
  // Trigger flag
  {
    int cid=DetIdTrigFlag;
    int nh = MTDCAna->GetNHits(cid);
    double patbins[3]={16,-0.5,16-0.5};
    for( int i=0; i<nh; ++i ){
      MTDCHit *hit = MTDCAna->GetHit(cid,i);
      MTDCRawHit *raw = hit->GetRawHit();
      int ntu=raw->GetSizeLeading();
      int seg = hit->SegmentId();
      if(ntu>0) {
	if(seg==2) BEAM=true;
	if(seg==5) KAON3=true;
      }
    }
  }
  //  if(!BEAM) return true;
  bool GoodBeam=true;
  int hodoseg[nhodo];
  double hodode[nhodo];
  for(int ihodo=0;ihodo<nhodo;++ihodo){
    hodoseg[ihodo]=-1;
    hodode[ihodo]=-1;
    int cid=Cid[ihodo];
    int nh = hodoAna->GetNHits(cid);
    int mulgate=0;
    for( int i=0; i<nh; ++i ){
      Hodo2Hit *hit = hodoAna->GetHit(cid,i);
      if(!hit) continue;
      int seg = hit->SegmentId();
      int nind=hit->GetIndex();
      for(int it=0;it<nind;it++){
	double mt  = hit->MeanTime(it);
	if(TMath::Abs(mt)<50){
	  hodoseg[ihodo]=seg;
	  hodode[ihodo]=hit->DeltaE();
	  mulgate++;
	}
      }
    }//for(i)
    if(cid!=DetIdBHD&&mulgate!=1) GoodBeam=false;
  }

  if(hodoseg[kCDH]<0) return true; 

#if 0
  for(int layer=0;layer<15;layer++){
    const DCHitContainer &hc=CDCAna->GetDCHC(layer);
    for( int i1=0; i1<hc.size(); ++i1 ){
      DCHit *hit1=hc[i1];
      int multi1 = hit1->GetDriftLengthSize();
      std::cout<<"===layer: "<<layer
	       <<"  ,wire : "<<hit1->GetWire()
	       <<"  ,posx : "<<hit1->GetWirePosition().X()
	       <<"  ,posy : "<<hit1->GetWirePosition().Y()<<std::endl;
      for ( int m1=0; m1<multi1; ++m1 ) {
	double dt=hit1->GetDriftTime(m1);
	double tot=hit1->GetTOT(m1);
	std::cout<<"   dt : "<<dt <<"  "<<hit1->IsWithinDtRange(m1);
	std::cout<<"  ,tot: "<<tot<<"  "<<hit1->IsWithinTotRange(m1)<<std::endl;
      }
    }
  }
#endif

  CDCAna->TrackSearch();
  int ntr=CDCAna->GetNTracks();
  //  if(ntr<1) return true;
  std::cout<<"event# "<<event_number<<", ntracks = "<<ntr<<std::endl;
  //  if(!GoodBeam) return true;
  c1->Clear();
  c1->Divide(2,2);
  TVirtualPad* pad; 
  pad=c1->cd(1);
  pad->DrawFrame(-70,-70,70,70,"CDC;X;Y"); // might cause memory leak
  std::cout<<"draw cdc layers xy"<<std::endl;
  disp::DrawCDCLayersXY(pad);
  std::cout<<"draw segments xy"<<std::endl;
  disp::DrawSegmentsXY(pad,DetIdCDH);
  std::cout<<"draw cds hits xy"<<std::endl;
  disp::DrawCDSHitsXY(pad,hodoAna,DetIdCDH);
  std::cout<<"draw cdc hits xy"<<std::endl;
  disp::DrawCDCHitsXY(pad,CDCAna);
  std::cout<<"draw cdc cluster hits xy"<<std::endl;
  disp::DrawCDCClusterHitsXY(pad,CDCAna);
  std::cout<<"draw cc tracks xy"<<std::endl;
  disp::DrawCDCTracksXY(pad,CDCAna);
  //
  pad=c1->cd(2);
  pad->DrawFrame(-70,-70,70,70,"CDC;Z;Y");
  std::cout<<"draw cdc layers zy"<<std::endl;
  disp::DrawCDCLayersZY(pad);
  std::cout<<"draw segments zy"<<std::endl;
  disp::DrawSegmentsZY(pad,DetIdCDH);
  disp::DrawCDSHitsZY(pad,hodoAna,DetIdCDH);
  disp::DrawCDCHitsZY(pad,CDCAna);
  disp::DrawCDCTracksZY(pad,CDCAna);

  pad=c1->cd(3);
  pad->DrawFrame(-50,0,50,60,"CDC;Z;R");
  std::cout<<"draw cdc layers zr"<<std::endl;
  disp::DrawCDCLayersZR(pad);
  std::cout<<"draw segments zr"<<std::endl;
  disp::DrawSegmentsZR(pad,DetIdCDH);
  disp::DrawCDSHitsZR(pad,hodoAna,DetIdCDH);
  disp::DrawCDCHitsZR(pad,CDCAna);
  disp::DrawCDCTracksZR(pad,CDCAna);

  pad=c1->cd(4);
  pad->DrawFrame(-50,-3.5,50,3.5,"CDC;Z;Phi");
  std::cout<<"draw cdc layers zphi"<<std::endl;
  disp::DrawCDSHitsZPhi(pad,hodoAna,DetIdCDH);
  disp::DrawCDCHitsZPhi(pad,CDCAna);
  disp::DrawCDCTracksZPhi(pad,CDCAna);

  c1->Update();
  c1->Print(pdfname);
  if(!disp::Wait())  return false;
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif  
  return true;
}

//______________________________________________________________________________
bool
EventRaw::ProcessingEnd( void )
{
  return true;
}

//______________________________________________________________________________
void
EventRaw::InitializeEvent( void )
{
  event_number=gUnpacker.get_event_number();
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventRaw;
}
//______________________________________________________________________________
bool
ConfMan::InitializeHistograms( void )
{  
  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return 
    ( InitializeParameter<HodoParamMan>("HDPRM","CDSPRM") ) &&
    ( InitializeParameter<HodoPHCMan>("HDPHC","CDSPHC") ) &&
    ( InitializeParameter<DCTdcCalibMan>("DCTDC") )  &&
    ( InitializeParameter<XTMapMan>("XTMap","CDCXT") )   &&
    ( InitializeParameter<UserParamMan>("USER") )   &&
    ( InitializeParameter<TransferMatrixMan>("TM") )   &&
    ( InitializeParameter<BLDCWireMapMan>("BLDCWire") ) &&
    ( InitializeParameter<CDCWireMapMan>("CDCGeom", "CDCASD" ) ) &&
    ( InitializeParameter<DetectorList>("DETLIST") ) &&
    ( InitializeParameter<GeomMapMan>("GeomBL", "GeomHall", "GeomCDS" ) );
}

//______________________________________________________________________________
bool
ConfMan::BeginRunProcess()
{
  new TApplication( "theApp", 0,0 );
  c1=new TCanvas("cdc","cdc",800,800);
  pdfname="tmpdisp.pdf";
  c1->Print(pdfname+"[");
  return true;
}
//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  c1->Print(pdfname+"]");
  return true;
}
