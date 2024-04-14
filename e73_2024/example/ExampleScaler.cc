#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "RootHelper.hh"
#include "ScalerAnalyzer.hh"
#include "RawData.hh"
#include "MTDCAnalyzer.hh"
#include "MTDCRawHit.hh"
#include "HodoParamMan.hh"
#include "HistTools.hh"
#include "VEvent.hh"
#include "TFile.h"

namespace
{
  using namespace root;
  const std::string& classname("EventScaler");
  ScalerAnalyzer& gScaler = ScalerAnalyzer::GetInstance();
  int run_number;
  int event_number;
  int spill_number;
  const int ntrig=32;
  int trigger_count[ntrig];
  const int MaxColumn=3;
  //  const int MaxRow=22;//elph
  //  const int MaxRow=16; //t77,e73fisrt
  const int MaxRow=17; //t98 phase1
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
class EventScaler : public VEvent
{
private:
  RawData *rawData;
  MTDCAnalyzer *MTDCAna;

public:
        EventScaler( void );
       ~EventScaler( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventScaler::EventScaler( void )
  : VEvent(),
    rawData(0),
    MTDCAna( new MTDCAnalyzer )
{
}
//______________________________________________________________________________
EventScaler::~EventScaler( void )
{
  if( MTDCAna ) delete MTDCAna;
  if( rawData ) delete rawData;  
}
//______________________________________________________________________________
namespace root
{
  TH1   *h[MaxHist];
  TGraph *gr[MaxColumn*MaxRow];
  TString trigname[16]={"SpillStart","SpillEnd","Beam","Pion",
			"Kaon2","Kaon3","KCDH1","KCDH2",
			"KCDH3","KCDH1g","Kg","PiCDH",
			"BeamPbF2","CDHCosmic","PbF2Cosmic","p_d"};
}
//______________________________________________________________________________
bool
EventScaler::ProcessingBegin( void )
{
  //  InitializeEvent();
  return true;
}

//______________________________________________________________________________
bool
EventScaler::ProcessingNormal( void )
{
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  gScaler.Decode();
  if(gScaler.SpillIncrement()){
    std::cout<<"spill increment"<<std::endl;
    spill_number++;
    gScaler.Print();
  }
  for(int i=0;i<MaxColumn;i++){
    for(int j=0;j<MaxRow;j++){
      if( !gScaler.FlagDisp(i,j) ) continue;
      double tmpx,tmpy;
      int tmpn=gr[i*MaxRow+j]->GetN();
      int val=gScaler.GetCurrentValue(i,j);
      gr[i*MaxRow+j]->GetPoint(tmpn-1,tmpx,tmpy);
      if(tmpx==spill_number){
	gr[i*MaxRow+j]->SetPoint(tmpn-1,tmpx,val);
      }else{
	gr[i*MaxRow+j]->SetPoint(tmpn,spill_number,val);
      }
    }
  }  
  rawData = new RawData;
  rawData->DecodeMTDCHits();
  MTDCAna->DecodeRawHits( rawData );
  // Trigger flag
  {
    int cid=DetIdTrigFlag;
    int nh = MTDCAna->GetNHits(cid);
    double mtdcbins[3]={5000,0,2e4};
    double patbins[3]={16,-0.5,16-0.5};
    for( int i=0; i<nh; ++i ){
      MTDCHit *hit = MTDCAna->GetHit(cid,i);
      MTDCRawHit *raw = hit->GetRawHit();
      int ntu=raw->GetSizeLeading();
      int seg = hit->SegmentId();
      TString tmpname="TriggerFlag";
      for(int it=0;it<ntu;it++){   
	double tu  = raw->GetLeading(it);
	hist::H1(Form("%s_TDC_%d",tmpname.Data(),seg),tu,mtdcbins);
      }
      if(ntu>0) {
	hist::H1(Form("%s_Pat",tmpname.Data()),seg,patbins);
	trigger_count[seg]++;
      }
    }
  }
  return true;
}

//______________________________________________________________________________
bool
EventScaler::ProcessingEnd( void )
{
  event_number++;
  return true;
}

//______________________________________________________________________________
void
EventScaler::InitializeEvent( void )
{
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventScaler;
}

//______________________________________________________________________________
namespace
{
}

//______________________________________________________________________________
bool
ConfMan::InitializeHistograms( void )
{
  for(int i=0;i<32;i++){
    trigger_count[i]=0;
  }
  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return  ( InitializeParameter<HodoParamMan>("HDPRM") );
}

//______________________________________________________________________________
bool
ConfMan::BeginRunProcess()
{
  run_number = get_run_number();
  event_number=0;
  spill_number=0;
  //  gScaler.Setup(FilePath("param/CMAP/scaler_e73_1st.param"));
  //  TString scaconf=FilePath("param/CMAP/scaler_t77.param");
  TString scaconf=FilePath("param/CMAP/scaler_e73_202306.param");
  std::cout<<scaconf.Data()<<std::endl;;
  gScaler.Setup(scaconf.Data());
  //gScaler.Setup(run_number);
  for(int i=0;i<MaxColumn;i++){
    for(int j=0;j<MaxRow;j++){
      if( !gScaler.FlagDisp(i,j) ) continue;
      gr[i*MaxRow+j]=new TGraph();
      std::string name=gScaler.GetName(i,j);
      std::cout<<i<<"  "<<j<<"  "<<name<<std::endl;
      gr[i*MaxRow+j]->SetName(Form("Scaler_%s",name.data()));
      gr[i*MaxRow+j]->SetMarkerStyle(20);
      gr[i*MaxRow+j]->GetXaxis()->SetTitle("Spill");
      gr[i*MaxRow+j]->GetYaxis()->SetTitle("Counts");
      gr[i*MaxRow+j]->GetXaxis()->CenterTitle();
      gr[i*MaxRow+j]->GetYaxis()->CenterTitle();
    }
  }

  return true;
}
//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  gScaler.Print();
  TString outname=gFile->GetName();
  outname.Remove(outname.Last('/')+1);
  gSystem->Exec("mkdir -p "+outname+"scaler");
  outname+=Form("scaler/scaler_run%05d.dat",run_number);
  std::cout<<outname<<std::endl;
  std::ofstream ofs(outname.Data());
  gScaler.WriteToFile(ofs);
  for(int i=0;i<16;i++){
    ofs<<Form(Form("Flag%d_",i))
       <<trigname[i]<<'\t'
       <<trigger_count[i]
       <<std::endl;
  }
  ofs.close();
  for(int i=0;i<MaxColumn;i++){
    for(int j=0;j<MaxRow;j++){
      if( !gScaler.FlagDisp(i,j) ) continue;
      gr[i*MaxRow+j]->Write();
    }
  }
  return true;
}
