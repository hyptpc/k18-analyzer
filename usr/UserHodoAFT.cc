/**
 *  file: UserHodoAFT.cc
 *  date: 2017.04.10
 *
 */

#include <iostream>
#include <sstream>
#include <cmath>

#include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "RootHelper.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"
#include "CFTPedCorMan.hh"
#include "FLHit.hh"

#define HodoCut 0 // with BH1/BH2
#define TimeCut 0 // in cluster analysis

namespace
{
  using namespace root;
  const std::string& classname("EventEasiroc");
  RMAnalyzer&         gRM   = RMAnalyzer::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
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
class EventHodoAFT : public VEvent
{
private:
  RawData      *rawData;
  DCAnalyzer   *DCAna;
  HodoAnalyzer *hodoAna;

public:
  EventHodoAFT( void );
  ~EventHodoAFT( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventHodoAFT::EventHodoAFT( void )
  : VEvent(),
    rawData(0),
    //  DCAna( new DCAnalyzer ),
    hodoAna( new HodoAnalyzer )
{
}

//______________________________________________________________________________
EventHodoAFT::~EventHodoAFT( void )
{
  if ( hodoAna ){
    delete hodoAna;
    hodoAna = 0;
  }

  if ( rawData ){
    delete rawData;
    rawData = 0;
  }
}

//______________________________________________________________________________


#ifndef MaxHits2
#define MaxHits2 20
#endif


struct Event
{
  int evnum;
  int spill;
  int spillnum=0;

  /*
    int trignhits;
    int trigpat[NumOfSegTrig];
    int trigflag[NumOfSegTrig];
  */
  int AFT_AdcHi[NumOfPlaneAFT][NumOfSegAFT_MAX][2];
  int AFT_AdcLow[NumOfPlaneAFT][NumOfSegAFT_MAX][2];
  //  int    AdcHiCor[NumOfPlaneAFT][NumOfSegAFT_MAX];
  //  int    AdcLowCor[NumOfPlaneAFT][NumOfSegAFT_MAX];
  int AFT_Tdc[NumOfPlaneAFT][NumOfSegAFT_MAX][2][MaxHits2];
  int AFT_TdcT[NumOfPlaneAFT][NumOfSegAFT_MAX][2][MaxHits2];
  //  int AFT_Hit[NumOfPlaneAFT][2][NumOfSegAFT_MAX];
  //  int AFT_CHit[NumOfPlaneAFT][2][NumOfSegAFT_MAX];
  int AFT_Hit[NumOfPlaneAFT][NumOfSegAFT_MAX][2];
  //  int AFT_CHit[NumOfPlaneAFT][NumOfSegAFT_MAX][2];
  //double Time[NumOfPlaneAFT][NumOfSegAFT_MAX];
  //double CTime[NumOfPlaneAFT][NumOfSegAFT_MAX];
  int AFT_tot[NumOfPlaneAFT][NumOfSegAFT_MAX][2][MaxHits2];

  int SST_AdcHi[NumOfPlaneSST][NumOfSegSST_MAX];
  int SST_AdcLow[NumOfPlaneSST][NumOfSegSST_MAX];
  int SST_Tdc[NumOfPlaneSST][NumOfSegSST_MAX][MaxHits2];
  int SST_Hit[NumOfPlaneSST][NumOfSegSST_MAX];
  int SST_tot[NumOfPlaneSST][NumOfSegSST_MAX][MaxHits2];

  //Fiber1-4
  int    FiberHits[NumOfPlaneAFT];
  double FiberSeg [NumOfPlaneAFT][MaxHits2];
  //double FiberTime[NumOfPlaneAFT][MaxHits2];
  //double FiberEdep[NumOfPlaneAFT][MaxHits2];
  int    FiberPID [NumOfPlaneAFT][MaxHits2];

  int TS_Adc[NumOfSegTS];
  int TS_Tdc[NumOfSegTS];

  /*
    double fClEdepForScat[NumOfPlaneAFT];
    double fClEdepNormForScat[NumOfPlaneAFT];
    double fClPathLengthForScat[NumOfPlaneAFT];

    double fClEdepForPi[NumOfPlaneAFT];
    double fClEdepNormForPi[NumOfPlaneAFT];
    double fClPathLengthForPi[NumOfPlaneAFT];
  */

};

//______________________________________________________________________________
namespace root
{
  Event  event;
  TH1   *h[MaxHist];
  TTree *tree;
  /*
    enum eDetHid
    {
    BFTHid  = 10000,
    SCHHid  = 20000,
    SFTVHid = 30000, SFTUHid = 40000, SFTXHid = 50000,
    FBT1U1Hid = 60000,  FBT1D1Hid = 70000,
    FBT1U2Hid = 80000,  FBT1D2Hid = 90000,
    FBT2U1Hid = 100000, FBT2D1Hid = 110000,
    FBT2U2Hid = 120000, FBT2D2Hid = 130000
    };
  */
}

//______________________________________________________________________________
bool
EventHodoAFT::ProcessingBegin( void )
{
  InitializeEvent();
  return true;
}

//______________________________________________________________________________
bool
EventHodoAFT::ProcessingNormal( void )
{
  static const std::string funcname("["+classname+"::"+__func__+"]");

#if HodoCut
  static const double MinDeBH2   = gUser.GetParameter("DeBH2", 0);
  static const double MaxDeBH2   = gUser.GetParameter("DeBH2", 1);
  static const double MinDeBH1   = gUser.GetParameter("DeBH1", 0);
  static const double MaxDeBH1   = gUser.GetParameter("DeBH1", 1);
  static const double MinBeamToF = gUser.GetParameter("BTOF",  0);
  static const double MaxBeamToF = gUser.GetParameter("BTOF",  1);
#endif
#if 0
  static const double MinTdcBFT  = gUser.GetParameter("TdcBFT",  0);
  static const double MaxTdcBFT  = gUser.GetParameter("TdcBFT",  1);
  static const double MinTdcSCH  = gUser.GetParameter("TdcSCH",  0);
  static const double MaxTdcSCH  = gUser.GetParameter("TdcSCH",  1);
  static const double MinTdcSFT  = gUser.GetParameter("TdcSFT",  0);
  static const double MaxTdcSFT  = gUser.GetParameter("TdcSFT",  1);
  static const double MinTdcFBT1 = gUser.GetParameter("TdcFBT1", 0);
  static const double MaxTdcFBT1 = gUser.GetParameter("TdcFBT1", 1);
  static const double MinTdcFBT2 = gUser.GetParameter("TdcFBT2", 0);
  static const double MaxTdcFBT2 = gUser.GetParameter("TdcFBT2", 1);
#endif
#if TimeCut
  static const double MinTimeBFT = gUser.GetParameter("TimeBFT", 0);
  static const double MaxTimeBFT = gUser.GetParameter("TimeBFT", 1);
  static const double MinTimeSCH = gUser.GetParameter("TimeSCH", 0);
  static const double MaxTimeSCH = gUser.GetParameter("TimeSCH", 1);
  static const double MinTimeSFT = gUser.GetParameter("TimeSFT", 0);
  static const double MaxTimeSFT = gUser.GetParameter("TimeSFT", 1);
  static const double MinTimeFBT1= gUser.GetParameter("TimeFBT1", 0);
  static const double MaxTimeFBT1= gUser.GetParameter("TimeFBT1", 1);
  static const double MinTimeFBT2= gUser.GetParameter("TimeFBT2", 0);
  static const double MaxTimeFBT2= gUser.GetParameter("TimeFBT2", 1);
#endif

  rawData = new RawData;
  rawData->DecodeHits();

  gRM.Decode();
  
  event.evnum++;
  int spill=gRM.SpillNumber();
  if(spill!=event.spill){
    int sdelta;
    if(spill-event.spill<0){
      sdelta=spill-event.spill+256;
    }else{
      sdelta=spill-event.spill;
    }
    event.spillnum+=sdelta;
    event.spill=spill;
  }

  

  HF1(1, 0);
  /*
  // Trigger Flag
  {
  const HodoRHitContainer &cont=rawData->GetTrigRawHC();
  int trignhits = 0;
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
  HodoRawHit *hit = cont[i];
  int seg = hit->SegmentId()+1;
  int tdc = hit->GetTdc1();
  if( tdc ){
  event.trigpat[trignhits++] = seg;
  event.trigflag[seg-1]      = tdc;
  HF1( 10, seg-1 );
  HF1( 10+seg, tdc );
  }
  }
  event.trignhits = trignhits;
  }
  */
  // if( event.trigflag[SpillEndFlag] ) return true;

  HF1(1, 1);
  /*
////////// BH2 time 0
hodoAna->DecodeBH2Hits(rawData);
int nhBh2 = hodoAna->GetNHitsBH2();
#if HodoCut
if ( nhBh2==0 ) return true;
#endif
HF1(1, 2);
double time0 = -999;
////////// BH2 Analysis
for(int i=0; i<nhBh2; ++i){
BH2Hit *hit = hodoAna->GetHitBH2(i);
if(!hit) continue;
double cmt = hit->CMeanTime();
double ct0 = hit->CTime0();
double min_time = -999;
#if HodoCut
double dE  = hit->DeltaE();
if( dE<MinDeBH2 || MaxDeBH2<dE )
continue;
#endif
if( std::abs(cmt)<std::abs(min_time) ){
min_time = cmt;
time0    = ct0;
}
}

HF1(1, 3);

////////// BH1 Analysis
hodoAna->DecodeBH1Hits(rawData);
int nhBh1 = hodoAna->GetNHitsBH1();
#if HodoCut
if(nhBh1==0) return true;
#endif
HF1(1, 4);
double btof0 = -999;
for(int i=0; i<nhBh1; ++i){
Hodo2Hit* hit = hodoAna->GetHitBH1(i);
if(!hit) continue;
double cmt  = hit->CMeanTime();
double btof = cmt - time0;
#if HodoCut
double dE   = hit->DeltaE();
if( dE<MinDeBH1 || MaxDeBH1<dE ) continue;
if( btof<MinBeamToF || MaxBeamToF<btof ) continue;
#endif
if( std::abs(btof)<std::abs(btof0) ){
btof0 = btof;
}
}

HF1(1, 5);

HF1(1, 6);
  */
  
  //*** AFT Raw hit ***//
  // {
  //   for (int layer=0; layer<NumOfPlaneAFT; layer++) {
  //     for(int UorD=0; UorD<2; UorD++){
  // 	const HodoRHitContainer &cont=rawData->GetAFTRawHC(layer, UorD);
  // 	int nh = cont.size();
  // 	for (int i=0; i<nh; i++){
  // 	  HodoRawHit *hit=cont[i];
  // 	  int seg=hit->SegmentId();
  // 	  //	  event.AFT_Hit[layer][UorD][seg]=seg;

  // 	  int A_high=hit->GetAdcUp(), A_low=hit->GetAdcDown();
	  
  // 	  // High Gain ADC  
  // 	  int histId = ((layer+1)*1000+10*seg+UorD)*10+1;
  // 	  //HF1( histId, double(A_high) );
	  
  // 	  event.AFT_AdcHi[layer][seg][UorD] = A_high; // AFT AdcHG

  // 	  if (A_high >= 1000.)
  // 	    event.AFT_CHit[layer][UorD][seg]=seg; // AFT Cut Hit number
	  
  // 	  histId = 100 + 10 * layer + 2;
  // 	  HF2( histId, seg, double(A_high) );
	  
  // 	  // Low Gain ADC   
  // 	  histId = ((layer+1)*1000+seg)*10+2;
  // 	  //HF1( histId, double(A_low) );
  // 	  //event.AFT_AdcLow[layer][seg][UorD] = A_low;
	  
  // 	  int nhit1 = hit->SizeTdc1();
  // 	  event.AFT_Hit[layer][UorD][seg]=nhit1; // AFT Hit number
	  
  // 	  if (nhit1>3) {
  // 	    std::cout << "nhit1 : " << nhit1 << std::endl;
  // 	    continue;
  // 	  }
	  
  // 	  int tdc_1st = 0;
  // 	  bool flag_tdc = false;
  // 	  for (int j=0; j<nhit1; j++) {
  // 	    int T_lead=hit->GetTdc1(j);
  // 	    event.AFT_Tdc[layer][seg][UorD][j] = T_lead;

  // 	    if (T_lead >= 400. && T_lead <= 430.)
  // 	      //	      event.AFT_CHit[layer][seg][UorD]=seg; // AFT Cut Hit number
  // 	      event.AFT_CHit[layer][UorD][seg]=seg; // AFT Cut Hit number

  // 	    // Leading TDC   
  // 	    histId = ((layer+1)*1000+seg)*10+3;
  // 	    //HF1( histId, double(T_lead) );
	    
  // 	    histId = 100 + 10 * layer + 1;
  // 	    HF2( histId, seg, double(T_lead) );
	    
  // 	    if (T_lead > tdc_1st)
  // 	      tdc_1st = T_lead;
	    
  // 	    if (T_lead >= 400. && T_lead <= 430.)
  // 	      flag_tdc = true;
  // 	  }
	  
	  
  // 	  // Trailing TDC  
  // 	  int nhit2 = hit->SizeTdc2();
  // 	  for (int j=0; j<nhit2; j++) {
  // 	    int T_trail=hit->GetTdcDown(j);
  // 	    histId = ((layer+1)*1000+seg)*10+4;
  // 	    //HF1( histId, double(T_trail) );
  // 	  }
	  
  // 	  // ADC HG w/ TDC  
  // 	  if (flag_tdc) {
  // 	    histId = ((layer+1)*1000+seg)*10+8;
  // 	    //HF1( histId, double(A_high));
	    
  // 	    histId = ((layer+1)*1000+seg)*10+9;
  // 	    //HF2( histId, double(A_high), double(A_low));
	    
  // 	    histId = 100 + 10 * layer + 3;
  // 	    HF2( histId, seg, double(A_high) );
	    
  // 	  }
	  
  // 	  // ADC HG % Leading TDC   
  // 	  histId = ((layer+1)*1000+seg)*10+5;
  // 	  //HF2( histId, double(A_high), double(tdc_1st) );
	  
  // 	  // ToT
  // 	  if (nhit1 == nhit2) {
  // 	    for (int j=0; j<nhit2; j++) {
  // 	      int T_lead=hit->GetTdc1(j);
  // 	      int T_trail=hit->GetTdcDown(j);
  // 	      int TOT = T_lead - T_trail;
	      
  // 	      event.AFT_tot[layer][seg][UorD][j]=TOT;
	      
  // 	      histId = ((layer+1)*1000+seg)*10+6;
  // 	      //HF2( histId, double(T_lead), double(TOT) );
	      
  // 	      if (TOT > 35) {
  // 		histId = ((layer+1)*1000+seg)*10+7;
  // 		//HF1(histId, double(T_lead));
  // 	      }
  // 	    }
  // 	  }
  // 	}
  //     }
  //   }
  // }

  // // SST Raw hit 
  // {
  //   for (int layer=0; layer<NumOfPlaneSST; layer++) {
  //     const HodoRHitContainer &cont=rawData->GetSSTRawHC(layer);
  //     int nh = cont.size();
  //     for (int i=0; i<nh; i++){
  // 	HodoRawHit *hit=cont[i];
  // 	int seg=hit->SegmentId();
       
  // 	int A_high=hit->GetAdcUp(), A_low=hit->GetAdcDown();
       
  // 	// High Gain ADC  
  // 	int histId = ((layer+1)*1000+10*seg)*10+1;
  // 	//HF1( histId, double(A_high) );
       
  // 	event.SST_AdcHi[layer][seg] = A_high;
       
  // 	histId = 100 + 10 * layer + 2;
  // 	HF2( histId, seg, double(A_high) );
       
  // 	// Low Gain ADC   
  // 	histId = ((layer+1)*1000+seg)*10+2;
  // 	//HF1( histId, double(A_low) );
       
  // 	event.SST_AdcLow[layer][seg] = A_low;
       
  // 	int nhit1 = hit->SizeTdc1();
       
  // 	if (nhit1>3) {
  // 	  //std::cout << "nhit1 : " << nhit1 << std::endl;
  // 	  continue;
  // 	}
       
  // 	int tdc_1st = 0;
  // 	bool flag_tdc = false;
  // 	for (int j=0; j<nhit1; j++) {
  // 	  int T_lead=hit->GetTdc1(j);
  // 	  // Leading TDC   
  // 	  histId = ((layer+1)*1000+seg)*10+3;
  // 	  //HF1( histId, double(T_lead) );
	 
  // 	  histId = 100 + 10 * layer + 1;
  // 	  HF2( histId, seg, double(T_lead) );
	 
  // 	  if (T_lead > tdc_1st)
  // 	    tdc_1st = T_lead;
	 
  // 	  if (T_lead >= 450. && T_lead <= 650.)
  // 	    flag_tdc = true;
	 
  // 	  event.SST_Tdc[layer][seg][j] = T_lead;
	 
  // 	}
       
       
  // 	// Trailing TDC  
  // 	int nhit2 = hit->SizeTdc2();
  // 	for (int j=0; j<nhit2; j++) {
  // 	  int T_trail=hit->GetTdcDown(j);
  // 	  histId = ((layer+1)*1000+seg)*10+4;
  // 	  //HF1( histId, double(T_trail) );
  // 	}
       
  // 	// ADC HG w/ TDC  
  // 	if (flag_tdc) {
  // 	  histId = ((layer+1)*1000+seg)*10+8;
  // 	  //HF1( histId, double(A_high));
	 
  // 	  histId = ((layer+1)*1000+seg)*10+9;
  // 	  //HF2( histId, double(A_high), double(A_low));
	 
  // 	  histId = 100 + 10 * layer + 3;
  // 	  HF2( histId, seg, double(A_high) );
	 
  // 	}
       
  // 	// ADC HG % Leading TDC   
  // 	histId = ((layer+1)*1000+seg)*10+5;
  // 	//HF2( histId, double(A_high), double(tdc_1st) );
       
  // 	if (nhit1 == nhit2) {
  // 	  for (int j=0; j<nhit2; j++) {
  // 	    int T_lead=hit->GetTdc1(j);
  // 	    int T_trail=hit->GetTdcDown(j);
  // 	    int TOT = T_lead - T_trail;
	   
  // 	    event.SST_tot[layer][seg][j]=TOT;
	   
  // 	    histId = ((layer+1)*1000+seg)*10+6;
  // 	    //HF2( histId, double(T_lead), double(TOT) );
	   
  // 	    if (TOT > 35) {
  // 	      histId = ((layer+1)*1000+seg)*10+7;
  // 	      //HF1(histId, double(T_lead));
  // 	    }
  // 	  }
  // 	}
  //     }
  //   }
  // }

  //*** AFT Raw hit made by tharada 20191121***//
  {
    for (int i=0; i<NumOfPlaneAFT; i++) {
      const HodoRHitContainer &cont=rawData->GetAFTRawHC(i);
      int nh = cont.size();
      for (int j=0; j<nh; j++){
	HodoRawHit *fiber=cont[j];
	int plane = fiber->PlaneId();
	int seg = fiber->SegmentId();
	for(int UorD = 0; UorD < 2; UorD++){
	  int histId = 0;
	  if(UorD == 0){
	    //** Up **//
	    // Hit
	    int nhitTdc1 = fiber->SizeTdc1();
	    event.AFT_Hit[plane][seg][UorD] = nhitTdc1; // AFT Hit number
	      
	    // if (nhitTdc1>3) {
	    //   std::cout << "nhitTdc1 : " << nhitTdc1 << std::endl;
	    //   continue;
	    // }

	    // Leading TDC   	      
	    int tdc_1st = 0;
	    bool flag_tdc = false;
	    for (int k=0; k<nhitTdc1; k++) {
	      int T_lead=fiber->GetTdc1(k);
	      event.AFT_Tdc[plane][seg][UorD][k] = T_lead;

	      histId = 100 + 10 * plane + 1;
	      HF2( histId, seg, double(T_lead) );
	    }
	  
	    // Trailing TDC  
	    int nhitTdcT1 = fiber->SizeTdcT1();
	    for (int k=0; k<nhitTdcT1; k++) {
	      int T_trail = fiber->GetTdcT1(k);
	      event.AFT_TdcT[plane][seg][UorD][k] = T_trail;
	    }

	    // ToT
	    if (nhitTdc1 == nhitTdcT1) {
	      for (int k=0; k<nhitTdcT1; k++) {
		int T_lead=fiber->GetTdc1(k);
		int T_trail=fiber->GetTdcT1(k);
		int TOT = T_lead - T_trail;
		  
		event.AFT_tot[plane][seg][UorD][k]=TOT;
	      }
	    }

	    // AdcHi
	    // int nhitAdcHi1 = fiber->SizeAdcHi1();
	    // for(int k = 0; k < nhitAdcHi1; k++){
	    //   int A_high = fiber->GetAdcHi1(k);
	    //   event.AFT_AdcHi[plane][seg][UorD][k] = A_high;
		
	    //   histId = 100 + 10 * plane + 2;
	    //   HF2( histId, seg, double(A_high) );
	    // }
	    {
	      int A_high = fiber->GetAdcHi1(0);
	      event.AFT_AdcHi[plane][seg][UorD] = A_high;
		
	      histId = 100 + 10 * plane + 2;
	      HF2( histId, seg, double(A_high) );
	    }

	    // AdcLow
	    // int nhitAdcLow1 = fiber->SizeAdcLow1();
	    // for(int k = 0; k < nhitAdcLow1; k++){
	    //   int A_low = fiber->GetAdcLow1(k);
	    //   event.AFT_AdcLow[plane][seg][UorD][k] = A_low;
	    // }
	    {
	      int A_low = fiber->GetAdcLow1(0);
	      event.AFT_AdcLow[plane][seg][UorD] = A_low;
	    }

	  }
	  else if(UorD == 1){
	    //** Down **//
	    int nhitTdc2 = fiber->SizeTdc2();
	    event.AFT_Hit[plane][seg][UorD] = nhitTdc2; // AFT Hit number
	      
	    // if (nhitTdc2>3) {
	    //   std::cout << "nhitTdc2 : " << nhitTdc2 << std::endl;
	    //   continue;
	    // }

	    // Leading TDC   	      
	    int tdc_1st = 0;
	    bool flag_tdc = false;
	    for (int k=0; k<nhitTdc2; k++) {
	      int T_lead=fiber->GetTdc2(k);
	      event.AFT_Tdc[plane][seg][UorD][k] = T_lead;

	      histId = 100 + 10 * plane + 1;
	      HF2( histId, seg, double(T_lead) );
	    }
	  
	    // Trailing TDC  
	    int nhitTdcT2 = fiber->SizeTdcT2();
	    for (int k=0; k<nhitTdcT2; k++) {
	      int T_trail = fiber->GetTdcT2(k);
	      event.AFT_TdcT[plane][seg][UorD][k] = T_trail;
	    }

	    // ToT
	    if (nhitTdc2 == nhitTdcT2) {
	      for (int k=0; k<nhitTdcT2; k++) {
		int T_lead=fiber->GetTdc2(k);
		int T_trail=fiber->GetTdcT2(k);
		int TOT = T_lead - T_trail;
		  
		event.AFT_tot[plane][seg][UorD][k]=TOT;
	      }
	    }

	    // AdcHi
	    // int nhitAdcHi2 = fiber->SizeAdcHi2();
	    // for(int k = 0; k < nhitAdcHi2; k++){
	    //   int A_high = fiber->GetAdcHi2(k);
	    //   event.AFT_AdcHi[plane][seg][UorD][k] = A_high;
		
	    //   histId = 100 + 10 * plane + 2;
	    //   HF2( histId, seg, double(A_high) );
	    // }
	    {
	      int A_high = fiber->GetAdcHi2(0);
	      event.AFT_AdcHi[plane][seg][UorD] = A_high;
		
	      histId = 100 + 10 * plane + 2;
	      HF2( histId, seg, double(A_high) );
	    }

	    // AdcLow
	    // int nhitAdcLow2 = fiber->SizeAdcLow2();
	    // for(int k = 0; k < nhitAdcLow2; k++){
	    //   int A_low = fiber->GetAdcLow2(k);
	    //   event.AFT_AdcLow[plane][seg][UorD][k] = A_low;
	    // }
	      int A_low = fiber->GetAdcLow2(0);
	      event.AFT_AdcLow[plane][seg][UorD] = A_low;

	  }
	}
      }
    }
  }

  //*** SST Raw hit made by tharada 20191121***// 
  {
    for (int i=0; i<NumOfPlaneSST; i++) {
      const HodoRHitContainer &cont=rawData->GetSSTRawHC(i);
      int nh = cont.size();
      for (int j=0; j<nh; j++){
	HodoRawHit *fiber=cont[j];
	int plane = fiber->PlaneId();
	int seg = fiber->SegmentId();

	// Hit
	int nhitTdc1 = fiber->SizeTdc1();
	event.SST_Hit[plane][seg] = nhitTdc1; // SST Hit number
	      
	if (nhitTdc1>3) {
	  std::cout << "nhitTdc1 : " << nhitTdc1 << std::endl;
	  continue;
	}
	
	// Leading TDC   	      
	int tdc_1st = 0;
	bool flag_tdc = false;
	for (int k=0; k<nhitTdc1; k++) {
	  int T_lead=fiber->GetTdc1(k);
	  event.SST_Tdc[plane][seg][k] = T_lead;
	}
	
	// Trailing TDC  
	int nhitTdcT1 = fiber->SizeTdcT1();
	for (int k=0; k<nhitTdcT1; k++) {
	  int T_trail = fiber->GetTdcT1(k);
	}
	
	// ToT
	if (nhitTdc1 == nhitTdcT1) {
	  for (int k=0; k<nhitTdcT1; k++) {
	    int T_lead=fiber->GetTdc1(k);
	    int T_trail=fiber->GetTdcT1(k);
	    int TOT = T_lead - T_trail;
	    
	    event.SST_tot[plane][seg][k]=TOT;
	  }
	}
	
	// AdcHi
	// int nhitAdcHi1 = fiber->SizeAdcHi1();
	// for(int k = 0; k < nhitAdcHi1; k++){
	//   int A_high = fiber->GetAdcHi1(k);
	//   event.SST_AdcHi[plane][seg][k] = A_high;
	// }
	{
	  int A_high = fiber->GetAdcHi1(0);
	  event.SST_AdcHi[plane][seg] = A_high;
	}
	
	// AdcLow
	// int nhitAdcLow1 = fiber->SizeAdcLow1();
	// for(int k = 0; k < nhitAdcLow1; k++){
	//   int A_low = fiber->GetAdcLow1(k);
	//   event.SST_AdcLow[plane][seg][k] = A_low;
	// }
	{
	  int A_low = fiber->GetAdcLow1(0);
	  event.SST_AdcLow[plane][seg] = A_low;
	}
      }
    }
  }

  /*
  // FiberHitAFT    
  hodoAna->DecodeAFTHits(rawData);
  for(int p = 0; p<NumOfPlaneAFT; ++p){

  int nhit = hodoAna->GetNHitsAFT(p);          

  int nhit_t = 0;
  for(int i = 0; i<nhit; ++i){
  const FiberHit* hit = hodoAna->GetHitAFT(p, i);
  int mhit = hit->GetNumOfHit();
  int seg = hit->PairId();
  double adcHi  = hit->GetAdcHi();
  double adcLow = hit->GetAdcLow();  

  double time0 = -999;
  double ctime0 = -999;
  double tot0 = -999;
  for(int m = 0; m<mhit; ++m){
  double ctime  = hit->GetCTime(m);
  double time   = hit->GetTime(m);
  double tot    = hit->GetTot(m);

  if (fabs(ctime0) > fabs(ctime)) {
  time0 = time;
  ctime0 = ctime;
  tot0 = tot;
  }
  }

  event.AdcHiCor[p][seg] = (int)adcHi;
  event.AdcLowCor[p][seg] = (int)adcLow;
  event.Time[p][seg] = time0;
  event.CTime[p][seg] = ctime0;
  event.tot[p][seg] = tot0;

  int mh_tdc = hit->GetNLeading();
  for(int m = 0; m<mh_tdc; ++m){
  double leading  = hit->GetLeading(m);

  int histId = ((p+1)*1000+seg)*10+7;
  if (1./sqrt(event.AdcHiCor[p][seg])<0.03 && event.AdcHiCor[p][seg]>0)
  HF1( histId, leading);

  }

  }//nhit
  }  

  for ( int l = 0; l < NumOfPlaneAFT; ++l ) {
    
  int ncl = hodoAna->GetNClustersAFT( l );    
  for ( int j = 0; j < ncl; ++j ) {
            
  FiberCluster* cl = hodoAna->GetClusterAFT( l, j );
  double mean_seg = cl->MeanSeg();
  double max_seg     = cl->MaxSeg();
  double size     = cl->ClusterSize();

  double time     = cl->CMeanTime();

  double max_adcHi  = cl->MaxAdcHi();
  double max_adcLow  = cl->MaxAdcLow();
  double max_mipLow  = cl->MaxMIPLow();
  double max_dELow   = cl->MaxdELow(); 

  int seg = (int)max_seg;

  if (time>-10 && time<10) {
  int histId = ((l+1)*1000+seg)*10+8;
  HF1( histId, max_adcHi);
  }      

  int max_cluster_id = cl->GetMaxClusterId(); 
  if (max_cluster_id>0) {
  FLHit *fhit = cl->GetHit(max_cluster_id);
  double ftime = fhit->GetTime();

  int histId = ((l+1)*1000+seg)*10+9;
  HF2( histId, 1./sqrt(max_adcHi), ftime);

  }
  }
  }
  */


  //*** TS Raw hit made by tharada 20210218***// 
  {
    const HodoRHitContainer &cont = rawData->GetTSRawHC();
    int nh = cont.size();
    for( int i = 0; i < nh; i++ ){
      HodoRawHit *ts = cont[i];

      int adc = ts->GetAdcHi1();
      int tdc = ts->GetTdc1();
      event.TS_Adc[i] = adc;
      event.TS_Tdc[i] = tdc;
    }
  }


  return true;
}

//______________________________________________________________________________
bool
EventHodoAFT::ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//______________________________________________________________________________
void
EventHodoAFT::InitializeEvent( void )
{
  //event.evnum      = 0;
  /*
    for( int it=0; it<NumOfSegTrig; it++){
    event.trigpat[it]  = -1;
    event.trigflag[it] = -1;
    }
  */
  /*
    for (int i=0; i<NumOfPlaneAFT; i++) {
    for (int j=0; j<NumOfSegAFT_MAX; j++) {
    //event.AdcHiCor[i][j] = -999;
    //event.AdcLowCor[i][j] = -999;
    //event.Time[i][j] = -999.;
    //event.CTime[i][j] = -999.;
    }
    }
  */

  for (int i=0; i<NumOfPlaneAFT; i++) {
    for (int j=0; j<NumOfSegAFT_MAX; j++) {
      for(int k=0; k<2; k++){
	event.AFT_Hit[i][j][k]=-1;
	//        event.AFT_CHit[i][k][j]=-1;
	event.AFT_AdcHi[i][j][k] = -999;
	event.AFT_AdcLow[i][j][k] = -999;
	for(int l=0; l<MaxHits2; l++){
	  event.AFT_Tdc[i][j][k][l] = -999;
	  event.AFT_TdcT[i][j][k][l] = -999;
	  event.AFT_tot[i][j][k][l] = -999.;
	}
      }
    }
  }

  for (int i=0; i<NumOfPlaneSST; i++) {
    for (int j=0; j<NumOfSegSST_MAX; j++) {
      event.SST_Hit[i][j]=-1;
      event.SST_AdcHi[i][j] = -999;
      event.SST_AdcLow[i][j] = -999;
      for(int l=0; l<MaxHits2; l++){
	event.SST_Tdc[i][j][l] = -999;
	event.SST_tot[i][j][l] = -999.;
      }
    }
  }

  for( int i = 0; i < NumOfSegTS; i++ ){
    event.TS_Adc[i] = -999;
    event.TS_Tdc[i] = -999;
  }

}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventHodoAFT;
}

//______________________________________________________________________________
namespace
{
  const int    NbinTdc = 1000;
  const double MinTdc  =    0.;
  const double MaxTdc  = 1000.;

  const int    NbinTot =  170;
  const double MinTot  =  -10.;
  const double MaxTot  =  160.;

  const int    NbinTime = 1000;
  const double MinTime  = -500.;
  const double MaxTime  =  500.;
}

//______________________________________________________________________________
bool
ConfMan:: InitializeHistograms( void )
{
  HB1(  1, "Status",  20,   0., 20. );
  /*
    HB1( 10, "Trigger HitPat", NumOfSegTrig, 0., double(NumOfSegTrig) );
    for(int i=0; i<NumOfSegTrig; ++i){
    HB1( 10+i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000 );
    }
  */

  char buf[100];
  for (int l=0; l<NumOfPlaneAFT; l++) {
    int hid = 100+10*l;
    int NumOfSeg = NumOfSegAFT[l];

    sprintf(buf, "TDC layer %d", l);
    HB2( hid+1, buf, NumOfSegAFT[l], 0, NumOfSegAFT[l], 400, 0, 4000 );

    sprintf(buf, "ADC layer %d", l);
    HB2( hid+2, buf, NumOfSegAFT[l], 0, NumOfSegAFT[l], 400, 0, 4000 );

    sprintf(buf, "ADC (w/ TDC) layer %d", l);
    HB2( hid+3, buf, NumOfSegAFT[l], 0, NumOfSegAFT[l], 400, 0, 4000 );

    sprintf(buf, "dE layer %d", l);
    HB1( hid+4, buf, 1000, 0, 10 );

    sprintf(buf, "Time layer %d", l);
    HB2( hid+5, buf, NumOfSeg, 0, NumOfSeg, 500, -250, 250 );
    sprintf(buf, "dE (High) layer %d", l);
    HB2( hid+6, buf, NumOfSeg, 0, NumOfSeg, 250, 0, 50 );
    sprintf(buf, "dE (Low) layer %d", l);
    HB2( hid+7, buf, NumOfSeg, 0, NumOfSeg, 250, 0, 50 );
    sprintf(buf, "dE (High) w/ TDC layer %d", l);
    HB2( hid+8, buf, NumOfSeg, 0, NumOfSeg, 250, 0, 50 );
    sprintf(buf, "dE (Low) w/ TDC layer %d", l);
    HB2( hid+9, buf, NumOfSeg, 0, NumOfSeg, 250, 0, 50 );

  }


  // Consume too many memories 
  // use only when the each channnel study is necessary 
  // for (int l=0; l<NumOfPlaneAFT; l++) {
  //   int NumOfSeg = NumOfSegAFT[l];
  //   for (int seg=0; seg<NumOfSeg; seg++ ) {
  //     int hid = ((l+1)*1000+seg)*10;
  //     /*
  // 	sprintf(buf, "ADC High Gain %d-%d", l, seg);
  // 	HB1(hid+1, buf, 4092, 0, 4092);
  // 	sprintf(buf, "ADC Low Gain %d-%d", l, seg);
  // 	HB1(hid+2, buf, 4092, 0, 4092);

  // 	sprintf(buf, "TDC Leading %d-%d", l, seg);
  // 	HB1(hid+3, buf, 1024, 0, 1024);

  // 	sprintf(buf, "TDC Traiding %d-%d", l, seg);
  // 	HB1(hid+4, buf, 1024, 0, 1024);
  // 	sprintf(buf, "ADC HG % TDC Leading (1st) %d-%d", l, seg);
  // 	HB2(hid+5, buf, 512, 0, 4096, 200, 500, 700);
  // 	sprintf(buf, "TDC Leading % TOT %d-%d", l, seg);
  // 	HB2(hid+6, buf, 200, 500, 700, 256, 0, 256);
  // 	sprintf(buf, "TDC Leading %d-%d (w/ TOT cut)", l, seg);
  // 	HB1(hid+7, buf, 1024, 0, 1024);
  //     */
  //     sprintf(buf, "ADC High Gain (w/ TDC) %d-%d", l, seg);
  //     HB1(hid+8, buf, 4092, 0, 4092);


  //     sprintf(buf, "ADC High Gain % Low Gain (w/ TDC) %d-%d", l, seg);
  //     HB2(hid+9, buf, 300, 0, 0.1, 100, -50, 50);

  //     sprintf(buf, "TDC (ADC Cut) %d-%d", l, seg);
  //     HB1(hid+7, buf, 1024, 0, 1024);

  //   }
  // }


  //Tree
  HBTree( "ea0c", "tree of Easiroc" );
  //Trig
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("spillnum",     &event.spillnum,     "spillnum/I");
  //tree->Branch("trignhits", &event.trignhits, "trignhits/I");
  //tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  //tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  //** AFT **//
  tree->Branch("AFT_AdcHi",event.AFT_AdcHi, Form("AFT_AdcHi[%d][%d][2]/I", NumOfPlaneAFT, NumOfSegAFT_MAX));
  tree->Branch("AFT_AdcLow",event.AFT_AdcLow, Form("AFT_AdcLow[%d][%d][2]/I", NumOfPlaneAFT, NumOfSegAFT_MAX));
  tree->Branch("AFT_Tdc",event.AFT_Tdc, Form("AFT_Tdc[%d][%d][2][%d]/I", NumOfPlaneAFT, NumOfSegAFT_MAX, MaxHits2));
  tree->Branch("AFT_TdcT",event.AFT_TdcT, Form("AFT_TdcT[%d][%d][2][%d]/I", NumOfPlaneAFT, NumOfSegAFT_MAX, MaxHits2));
  tree->Branch("AFT_Hit",event.AFT_Hit, Form("AFT_Hit[%d][%d][2]/I", NumOfPlaneAFT, NumOfSegAFT_MAX));
  //  tree->Branch("AFT_CHit",event.AFT_CHit, Form("AFT_Hit[%d][2][%d]/I",NumOfPlaneAFT, NumOfSegAFT_MAX));
  tree->Branch("AFT_tot",event.AFT_tot, Form("AFT_tot[%d][%d][2][%d]/I", NumOfPlaneAFT, NumOfSegAFT_MAX, MaxHits2));

  //** SST **//
  tree->Branch("SST_AdcHi",event.SST_AdcHi, Form("SST_AdcHi[%d][%d]/I", NumOfPlaneSST, NumOfSegSST_MAX));
  tree->Branch("SST_AdcLow",event.SST_AdcLow, Form("SST_AdcLow[%d][%d]/I", NumOfPlaneSST, NumOfSegSST_MAX));
  tree->Branch("SST_Tdc",event.SST_Tdc, Form("SST_Tdc[%d][%d][%d]/I", NumOfPlaneSST, NumOfSegSST_MAX, MaxHits2));
  tree->Branch("SST_Hit",event.SST_Hit, Form("SST_Hit[%d][%d]/I", NumOfPlaneSST, NumOfSegSST_MAX));
  tree->Branch("SST_tot",event.SST_tot, Form("SST_tot[%d][%d][%d]/I", NumOfPlaneSST, NumOfSegSST_MAX, MaxHits2));

  //** TS **//
  tree->Branch("TS_Adc", event.TS_Adc, Form("TS_Adc[%d]/I", NumOfSegTS));
  tree->Branch("TS_Tdc", event.TS_Tdc, Form("TS_Tdc[%d]/I", NumOfSegTS));

  /*		 
		 for (int i=0; i<NumOfPlaneAFT; i++) {
		 char buf1[100], buf2[100];



		 /*
		 sprintf(buf1, "AdcHiCor%d", i);
		 sprintf(buf2, "AdcHiCor%d[%d]/I", i, NumOfSegAFT[i]);
		 tree->Branch(buf1,   event.AdcHiCor[i],  buf2);

		 sprintf(buf1, "AdcLowCor%d", i);
		 sprintf(buf2, "AdcLowCor%d[%d]/I", i, NumOfSegAFT[i]);
		 tree->Branch(buf1,   event.AdcLowCor[i],  buf2);
  */

  /*
    sprintf(buf1, "Time%d", i);
    sprintf(buf2, "Time%d[%d]/D", i, NumOfSegAFT[i]);
    tree->Branch(buf1,   event.Time[i],  buf2);

    sprintf(buf1, "CTime%d", i);
    sprintf(buf2, "CTime%d[%d]/D", i, NumOfSegAFT[i]);
    tree->Branch(buf1,   event.CTime[i],  buf2);
    

    }
  */

  HPrint();
  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")    &&
      InitializeParameter<HodoParamMan>("HDPRM") &&
      InitializeParameter<HodoPHCMan>("HDPHC")   &&
      //InitializeParameter<CFTPedCorMan>("CFTPED")      &&
      InitializeParameter<UserParamMan>("USER")  );
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
