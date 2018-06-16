/**
 *  file: UserHimac.cc
 *  date: 2017.04.10
 *
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "SsdCluster.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"

#define HodoCut   0
#define MaxMultiCut 0

// SSD Filter
#define SlopeFilter     0
#define DeltaEFilter    1
#define TimeFilter      1
#define ChisqrFilter    1

namespace
{
  using namespace root;
  const std::string class_name("EventHimac");
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
class EventHimac : public VEvent
{
private:
  RawData      *rawData;
  DCAnalyzer   *DCAna;
  HodoAnalyzer *hodoAna;

public:
        EventHimac( void );
       ~EventHimac( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventHimac::EventHimac( void )
  : VEvent(),
    rawData(0),
    DCAna( new DCAnalyzer ),
    hodoAna( new HodoAnalyzer )
{
}

//______________________________________________________________________________
EventHimac::~EventHimac( void )
{
  if( DCAna )   delete DCAna;
  if( hodoAna ) delete hodoAna;
  if( rawData ) delete rawData;
}

//______________________________________________________________________________
bool
EventHimac::ProcessingBegin( void )
{
  InitializeEvent();
  return true;
}

//______________________________________________________________________________
struct Event
{
  int evnum;
  
  unsigned long unixtime;

  int trignhits;
  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

  // T1
  int t1nhits;
  int t1hitpat[MaxHits];
  double t1ua[NumOfSegT1];
  double t1ut[NumOfSegT1];
  double t1da[NumOfSegT1];
  double t1dt[NumOfSegT1];

  // T2
  int t2nhits;
  int t2hitpat[MaxHits];
  double t2ua[NumOfSegT2];
  double t2ut[NumOfSegT2];
  double t2da[NumOfSegT2];
  double t2dt[NumOfSegT2];

  // T3
  int t3nhits;
  int t3hitpat[MaxHits];
  double t3ua[NumOfSegT3];
  double t3ut[NumOfSegT3];
  double t3da[NumOfSegT3];
  double t3dt[NumOfSegT3];

  // T4
  int t4nhits;
  int t4hitpat[MaxHits];
  double t4ua[NumOfSegT4];
  double t4ut[NumOfSegT4];
  double t4da[NumOfSegT4];
  double t4dt[NumOfSegT4];

  // S1
  int s1nhits;
  int s1hitpat[MaxHits];
  double s1ua[NumOfSegS1];
  double s1ut[NumOfSegS1];
  double s1da[NumOfSegS1];
  double s1dt[NumOfSegS1];

  // S2
  int s2nhits;
  int s2hitpat[MaxHits];
  double s2ua[NumOfSegS2];
  double s2ut[NumOfSegS2];
  double s2da[NumOfSegS2];
  double s2dt[NumOfSegS2];

  // SsdHit
  int    ssd1y0nhits;
  double ssd1y0hitpat[NumOfSegSSD1];
  double ssd1y0pos[NumOfSegSSD1];
  double ssd1y0ped[NumOfSegSSD1];
  double ssd1y0adc[NumOfSegSSD1];
  double ssd1y0tdc[NumOfSegSSD1];
  double ssd1y0de[NumOfSegSSD1];
  double ssd1y0time[NumOfSegSSD1];
  double ssd1y0waveform[NumOfSegSSD1][NumOfSampleSSD];
  double ssd1y0chisqr[NumOfSegSSD1];

  int    ssd1x0nhits;
  double ssd1x0hitpat[NumOfSegSSD1];
  double ssd1x0pos[NumOfSegSSD1];
  double ssd1x0ped[NumOfSegSSD1];
  double ssd1x0adc[NumOfSegSSD1];
  double ssd1x0tdc[NumOfSegSSD1];
  double ssd1x0sum[NumOfSegSSD1];
  double ssd1x0de[NumOfSegSSD1];
  double ssd1x0time[NumOfSegSSD1];
  double ssd1x0waveform[NumOfSegSSD1][NumOfSampleSSD];
  double ssd1x0chisqr[NumOfSegSSD1];

  int    ssd2x0nhits;
  double ssd2x0hitpat[NumOfSegSSD2];
  double ssd2x0pos[NumOfSegSSD2];
  double ssd2x0ped[NumOfSegSSD2];
  double ssd2x0adc[NumOfSegSSD2];
  double ssd2x0tdc[NumOfSegSSD2];
  double ssd2x0de[NumOfSegSSD2];
  double ssd2x0time[NumOfSegSSD2];
  double ssd2x0waveform[NumOfSegSSD2][NumOfSampleSSD];
  double ssd2x0chisqr[NumOfSegSSD2];

  int    ssd2y0nhits;
  double ssd2y0hitpat[NumOfSegSSD2];
  double ssd2y0pos[NumOfSegSSD2];
  double ssd2y0ped[NumOfSegSSD2];
  double ssd2y0adc[NumOfSegSSD2];
  double ssd2y0tdc[NumOfSegSSD2];
  double ssd2y0de[NumOfSegSSD2];
  double ssd2y0time[NumOfSegSSD2];
  double ssd2y0waveform[NumOfSegSSD2][NumOfSampleSSD];
  double ssd2y0chisqr[NumOfSegSSD2];

  // SsdCluster
  int    ssd1y0ncl;
  int    ssd1y0clsize[NumOfSegSSD1];
  double ssd1y0clpos[NumOfSegSSD1];
  double ssd1y0clde[NumOfSegSSD1];
  double ssd1y0cltime[NumOfSegSSD1];
  double ssd1y0cltd[NumOfSegSSD1];
  int    ssd1x0ncl;
  int    ssd1x0clsize[NumOfSegSSD1];
  double ssd1x0clpos[NumOfSegSSD1];
  double ssd1x0clde[NumOfSegSSD1];
  double ssd1x0cltime[NumOfSegSSD1];
  double ssd1x0cltd[NumOfSegSSD1];

  int    ssd2x0ncl;
  int    ssd2x0clsize[NumOfSegSSD2];
  double ssd2x0clpos[NumOfSegSSD2];
  double ssd2x0clde[NumOfSegSSD2];
  double ssd2x0cltime[NumOfSegSSD2];
  double ssd2x0cltd[NumOfSegSSD2];
  int    ssd2y0ncl;
  int    ssd2y0clsize[NumOfSegSSD2];
  double ssd2y0clpos[NumOfSegSSD2];
  double ssd2y0clde[NumOfSegSSD2];
  double ssd2y0cltime[NumOfSegSSD2];
  double ssd2y0cltd[NumOfSegSSD2];

  // SSDT
  int    ssdtnhits;
  double ssdthitpat[NumOfSegSSDT];
  double ssdtdc[NumOfSegSSDT];
  double ssdtime[NumOfSegSSDT];

};

//______________________________________________________________________________
namespace root
{
  Event  event;
  TH1   *h[MaxHist];
  TTree *tree;
  enum eDetHid  { SSD1Hid=10000, SSD2Hid=20000,
		  SSD1ClHid=30000, SSD2ClHid=40000,
		  SSDTHid=50000, 
		  T1Hid  =  60000,
		  T2Hid  =  70000,
		  T3Hid  =  80000,
		  T4Hid  =  90000,
		  S1Hid  = 100000,
		  S2Hid  = 110000,
		  nDetHid };
  enum eStatus  { All, Good, nStatus };
}

//______________________________________________________________________________
bool
EventHimac::ProcessingNormal( void )
{
  const std::string funcname("["+class_name+"::"+__func__+"]");

#if HodoCut
  static const double MinDeBH2   = gUser.GetParameter("DeBH2", 0);
  static const double MaxDeBH2   = gUser.GetParameter("DeBH2", 1);
  static const double MinDeBH1   = gUser.GetParameter("DeBH1", 0);
  static const double MaxDeBH1   = gUser.GetParameter("DeBH1", 1);
  static const double MinBeamToF = gUser.GetParameter("BTOF",  1);
  static const double MaxBeamToF = gUser.GetParameter("BTOF",  1);
#endif
  static const double MinDeSSDKaon = gUser.GetParameter("DeSSDKaon", 0);
  static const double MaxDeSSDKaon = gUser.GetParameter("DeSSDKaon", 1);
  static const double MinTimeSSD   = gUser.GetParameter("TimeSSD",   0);
  static const double MaxTimeSSD   = gUser.GetParameter("TimeSSD",   1);
  static const double MaxChisqrSSD = gUser.GetParameter("MaxChisqrSSD", 0);

  rawData = new RawData;
  rawData->DecodeHits();

  gRM.Decode();

  event.evnum = gRM.EventNumber();

  HF1(1, 0);

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
	HF1( 10,     seg-1 );
	HF1( 10+seg, tdc   );
      }
    }
    event.trignhits = trignhits;
  }

  // if( event.trigflag[SpillEndFlag] ) return true;

  HF1( 1, 1. );
  
  // T1
  {
    int t1_nhits = 0;
    const HodoRHitContainer &cont = rawData->GetT1RawHC();
    int nh = cont.size();
    HF1( T1Hid +0, double(nh) );
    int nh1 = 0, nh2 = 0;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      int seg = hit->SegmentId()+1;
      HF1( T1Hid +1, seg-0.5 );
      int Au = hit->GetAdcUp(), Ad = hit->GetAdcDown();
      int Tu = hit->GetTdcUp(), Td = hit->GetTdcDown();

      //Tree
      event.t1ua[seg-1] = Au;
      event.t1ut[seg-1] = Tu;
      event.t1da[seg-1] = Ad;
      event.t1dt[seg-1] = Td;
      if( Tu>0 && Td>0 ){
	event.t1hitpat[t1_nhits] = seg;
	t1_nhits++;
      }
      //Up
      HF1( T1Hid +100*seg +1, double(Au) );
      if( Tu>0 ){
	HF1( T1Hid +100*seg +3, double(Tu) );
	HF1( T1Hid +100*seg +5, double(Au) );
      }
      else{
	HF1( T1Hid +100*seg +7, double(Au) );
      }
      //Down
      HF1( T1Hid +100*seg +2, double(Ad) );
      if( Td>0 ){
	HF1( T1Hid +100*seg +4, double(Td) );
	HF1( T1Hid +100*seg +6, double(Ad) );
      }
      else{
	HF1( T1Hid +100*seg +8, double(Ad) );
      }
      //HitPat
      if( Tu>0 || Td>0 ){
	++nh1; HF1( T1Hid +3, seg-0.5 );
      }
      if( Tu>0 && Td>0 ){
	++nh2; HF1( T1Hid +5, seg-0.5 );
      }
    }
    HF1( T1Hid +2, double(nh1) ); HF1( T1Hid +4, double(nh2) );
    event.t1nhits = t1_nhits;
  }
  
  HF1( 1, 2. );
  
  // T2
  {
    int t2_nhits = 0;
    const HodoRHitContainer &cont = rawData->GetT2RawHC();
    int nh = cont.size();
    HF1( T2Hid +0, double(nh) );
    int nh1 = 0, nh2 = 0;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      int seg = hit->SegmentId()+1;
      HF1( T2Hid +1, seg-0.5 );
      int Au = hit->GetAdcUp(), Ad = hit->GetAdcDown();
      int Tu = hit->GetTdcUp(), Td = hit->GetTdcDown();

      //Tree
      event.t2ua[seg-1] = Au;
      event.t2ut[seg-1] = Tu;
      event.t2da[seg-1] = Ad;
      event.t2dt[seg-1] = Td;
      if( Tu>0 && Td>0 ){
	event.t2hitpat[t2_nhits] = seg;
	t2_nhits++;
      }
      //Up
      HF1( T2Hid +100*seg +1, double(Au) );
      if( Tu>0 ){
	HF1( T2Hid +100*seg +3, double(Tu) );
	HF1( T2Hid +100*seg +5, double(Au) );
      }
      else{
	HF1( T2Hid +100*seg +7, double(Au) );
      }
      //Down
      HF1( T2Hid +100*seg +2, double(Ad) );
      if( Td>0 ){
	HF1( T2Hid +100*seg +4, double(Td) );
	HF1( T2Hid +100*seg +6, double(Ad) );
      }
      else{
	HF1( T2Hid +100*seg +8, double(Ad) );
      }
      //HitPat
      if( Tu>0 || Td>0 ){
	++nh1; HF1( T2Hid +3, seg-0.5 );
      }
      if( Tu>0 && Td>0 ){
	++nh2; HF1( T2Hid +5, seg-0.5 );
      }
    }
    HF1( T2Hid +2, double(nh1) ); HF1( T2Hid +4, double(nh2) );
    event.t2nhits = t2_nhits;
  }
  
  HF1( 1, 3. );
  // T3
  {
    int t3_nhits = 0;
    const HodoRHitContainer &cont = rawData->GetT3RawHC();
    int nh = cont.size();
    HF1( T3Hid +0, double(nh) );
    int nh1 = 0, nh2 = 0;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      int seg = hit->SegmentId()+1;
      HF1( T3Hid +1, seg-0.5 );
      int Au = hit->GetAdcUp(), Ad = hit->GetAdcDown();
      int Tu = hit->GetTdcUp(), Td = hit->GetTdcDown();

      //Tree
      event.t3ua[seg-1] = Au;
      event.t3ut[seg-1] = Tu;
      event.t3da[seg-1] = Ad;
      event.t3dt[seg-1] = Td;
      if( Tu>0 && Td>0 ){
	event.t3hitpat[t3_nhits] = seg;
	t3_nhits++;
      }
      //Up
      HF1( T3Hid +100*seg +1, double(Au) );
      if( Tu>0 ){
	HF1( T3Hid +100*seg +3, double(Tu) );
	HF1( T3Hid +100*seg +5, double(Au) );
      }
      else{
	HF1( T3Hid +100*seg +7, double(Au) );
      }
      //Down
      HF1( T3Hid +100*seg +2, double(Ad) );
      if( Td>0 ){
	HF1( T3Hid +100*seg +4, double(Td) );
	HF1( T3Hid +100*seg +6, double(Ad) );
      }
      else{
	HF1( T3Hid +100*seg +8, double(Ad) );
      }
      //HitPat
      if( Tu>0 || Td>0 ){
	++nh1; HF1( T3Hid +3, seg-0.5 );
      }
      if( Tu>0 && Td>0 ){
	++nh2; HF1( T3Hid +5, seg-0.5 );
      }
    }
    HF1( T3Hid +2, double(nh1) ); HF1( T3Hid +4, double(nh2) );
    event.t3nhits = t3_nhits;
  }
  
  HF1( 1, 4. );
  // T4
  {
    int t4_nhits = 0;
    const HodoRHitContainer &cont = rawData->GetT4RawHC();
    int nh = cont.size();
    HF1( T4Hid +0, double(nh) );
    int nh1 = 0, nh2 = 0;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      int seg = hit->SegmentId()+1;
      HF1( T4Hid +1, seg-0.5 );
      int Au = hit->GetAdcUp(), Ad = hit->GetAdcDown();
      int Tu = hit->GetTdcUp(), Td = hit->GetTdcDown();

      //Tree
      event.t4ua[seg-1] = Au;
      event.t4ut[seg-1] = Tu;
      event.t4da[seg-1] = Ad;
      event.t4dt[seg-1] = Td;
      if( Tu>0 && Td>0 ){
	event.t4hitpat[t4_nhits] = seg;
	t4_nhits++;
      }
      //Up
      HF1( T4Hid +100*seg +1, double(Au) );
      if( Tu>0 ){
	HF1( T4Hid +100*seg +3, double(Tu) );
	HF1( T4Hid +100*seg +5, double(Au) );
      }
      else{
	HF1( T4Hid +100*seg +7, double(Au) );
      }
      //Down
      HF1( T4Hid +100*seg +2, double(Ad) );
      if( Td>0 ){
	HF1( T4Hid +100*seg +4, double(Td) );
	HF1( T4Hid +100*seg +6, double(Ad) );
      }
      else{
	HF1( T4Hid +100*seg +8, double(Ad) );
      }
      //HitPat
      if( Tu>0 || Td>0 ){
	++nh1; HF1( T4Hid +3, seg-0.5 );
      }
      if( Tu>0 && Td>0 ){
	++nh2; HF1( T4Hid +5, seg-0.5 );
      }
    }
    HF1( T4Hid +2, double(nh1) ); HF1( T4Hid +4, double(nh2) );
    event.t4nhits = t4_nhits;
  }
  
  HF1( 1, 5. );
  // S1
  {
    int s1_nhits = 0;
    const HodoRHitContainer &cont = rawData->GetS1RawHC();
    int nh = cont.size();
    HF1( S1Hid +0, double(nh) );
    int nh1 = 0, nh2 = 0;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      int seg = hit->SegmentId()+1;
      HF1( S1Hid +1, seg-0.5 );
      int Au = hit->GetAdcUp(), Ad = hit->GetAdcDown();
      int Tu = hit->GetTdcUp(), Td = hit->GetTdcDown();

      //Tree
      event.s1ua[seg-1] = Au;
      event.s1ut[seg-1] = Tu;
      event.s1da[seg-1] = Ad;
      event.s1dt[seg-1] = Td;
      if( Tu>0 && Td>0 ){
	event.s1hitpat[s1_nhits] = seg;
	s1_nhits++;
      }
      //Up
      HF1( S1Hid +100*seg +1, double(Au) );
      if( Tu>0 ){
	HF1( S1Hid +100*seg +3, double(Tu) );
	HF1( S1Hid +100*seg +5, double(Au) );
      }
      else{
	HF1( S1Hid +100*seg +7, double(Au) );
      }
      //Down
      HF1( S1Hid +100*seg +2, double(Ad) );
      if( Td>0 ){
	HF1( S1Hid +100*seg +4, double(Td) );
	HF1( S1Hid +100*seg +6, double(Ad) );
      }
      else{
	HF1( S1Hid +100*seg +8, double(Ad) );
      }
      //HitPat
      if( Tu>0 || Td>0 ){
	++nh1; HF1( S1Hid +3, seg-0.5 );
      }
      if( Tu>0 && Td>0 ){
	++nh2; HF1( S1Hid +5, seg-0.5 );
      }
    }
    HF1( S1Hid +2, double(nh1) ); HF1( S1Hid +4, double(nh2) );
    event.s1nhits = s1_nhits;
  }
  
  HF1( 1, 6. );
  // S2
  {
    int s2_nhits = 0;
    const HodoRHitContainer &cont = rawData->GetS2RawHC();
    int nh = cont.size();
    HF1( S2Hid +0, double(nh) );
    int nh1 = 0, nh2 = 0;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      int seg = hit->SegmentId()+1;
      HF1( S2Hid +1, seg-0.5 );
      int Au = hit->GetAdcUp(), Ad = hit->GetAdcDown();
      int Tu = hit->GetTdcUp(), Td = hit->GetTdcDown();

      //Tree
      event.s2ua[seg-1] = Au;
      event.s2ut[seg-1] = Tu;
      event.s2da[seg-1] = Ad;
      event.s2dt[seg-1] = Td;
      if( Tu>0 && Td>0 ){
	event.s2hitpat[s2_nhits] = seg;
	s2_nhits++;
      }
      //Up
      HF1( S2Hid +100*seg +1, double(Au) );
      if( Tu>0 ){
	HF1( S2Hid +100*seg +3, double(Tu) );
	HF1( S2Hid +100*seg +5, double(Au) );
      }
      else{
	HF1( S2Hid +100*seg +7, double(Au) );
      }
      //Down
      HF1( S2Hid +100*seg +2, double(Ad) );
      if( Td>0 ){
	HF1( S2Hid +100*seg +4, double(Td) );
	HF1( S2Hid +100*seg +6, double(Ad) );
      }
      else{
	HF1( S2Hid +100*seg +8, double(Ad) );
      }
      //HitPat
      if( Tu>0 || Td>0 ){
	++nh1; HF1( S2Hid +3, seg-0.5 );
      }
      if( Tu>0 && Td>0 ){
	++nh2; HF1( S2Hid +5, seg-0.5 );
      }
    }
    HF1( S2Hid +2, double(nh1) ); HF1( S2Hid +4, double(nh2) );
    event.s2nhits = s2_nhits;
  }
  
  //SSDT
  {
    const HodoRHitContainer &cont = rawData->GetSSDTRawHC();
    int ssdtnhits = 0;
    int nh = cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      int seg = hit->SegmentId()+1;
      int tdc = hit->GetTdc1();
      if( tdc>0 ){
	event.ssdthitpat[ssdtnhits++] = seg;
	event.ssdtdc[seg-1]           = tdc;
	HF1( SSDTHid +1, seg-1 );
	HF1( SSDTHid +seg*1000 +5, tdc );
	HF2( SSDTHid +6,  seg-1, tdc );
      }
    }
    event.ssdtnhits = ssdtnhits;
    HF1( SSDTHid +0, ssdtnhits );
  }

  std::vector<double> t0Ssd( NumOfSegSSDT, 0. );
  {
    hodoAna->DecodeSSDTHits( rawData );
    int nh = hodoAna->GetNHitsSSDT();
    for( int i=0; i<nh; ++i ){
      Hodo1Hit *hit = hodoAna->GetHitSSDT(i);
      if(!hit) continue;
      int    seg  = hit->SegmentId()+1;
      double time = hit->Time();
      event.ssdtime[seg-1] = time;
      HF1( SSDTHid +seg*1000 +14, time );
      HF2( SSDTHid +15,  seg-1, time );
      t0Ssd[seg-1] = time;
    }
  }

  DCAna->DecodeSsdInHits( rawData );
  DCAna->DecodeSsdOutHits( rawData );

  //////////////////// Filtering
  DCAna->DoTimeCorrectionSsd( t0Ssd );
  HF1( 1, 11. );

#if SlopeFilter
  DCAna->SlopeFilterSsd();
#endif

#if ChisqrFilter
  DCAna->ChisqrFilterSsd( MaxChisqrSSD );
#endif

  // re-clusterize
  DCAna->ClusterizeSsd();

#if DeltaEFilter
  DCAna->DeltaEFilterSsd( MinDeSSDKaon, MaxDeSSDKaon, true );
#endif

#if TimeFilter
  DCAna->TimeFilterSsd( MinTimeSSD, MaxTimeSSD, true );
#endif

  ///// SsdIn
  double multiSsdIn[nStatus]  = {};
  int nhSsdIn[NumOfLayersSSD1][nStatus];
  for( int l=0; l<NumOfLayersSSD1; ++l ){
    for( int s=0; s<nStatus; ++s ){
      nhSsdIn[l][s] = 0;
    }
  }

  for( int layer=1; layer<=NumOfLayersSSD1; ++layer ){
    const DCHitContainer &contIn = DCAna->GetSsdInHC(layer);
    int nhIn = contIn.size();
    multiSsdIn[All] += double(nhIn);
    nhSsdIn[layer-1][All] = nhIn;
    for( int i=0; i<nhIn; ++i ){
      DCHit  *hit = contIn[i];
      if(!hit) continue;
      int    seg  = hit->GetWire();
      bool   flag = hit->IsZeroSuppressed();
      bool   good = hit->IsGoodWaveForm();
      double ped  = hit->GetPedestal();
      double adc  = hit->GetAdcPeakHeight() - ped;
      double tdc  = hit->GetAdcPeakPosition();
      double pos  = hit->GetWirePosition();
      double peaktime = hit->GetPeakTime();
      double sum  = hit->GetAdcSum();
      double de   = hit->GetDe(),            chisqr   = hit->GetChisquare();
      std::vector<double> waveform = hit->GetWaveform();
      std::vector<double> time     = hit->GetTime();
      if( !flag ) continue;
      if( good ) multiSsdIn[Good]++;

      for( int status=0; status<good+1; ++status ){
	int offset = SSD1Hid+layer*1000 +status*100;
	HF1(  1 +offset, seg      ); HF1( 2    +offset, pos      );
	HF1(  3 +offset, adc      );
	HF2(  4 +offset, seg, adc );
	HF1(  5 +offset, tdc      ); HF2(  6 +offset, seg, tdc );
	HF2(  8 +offset, seg, ped );
	HF1(  9 +offset, de       ); HF2( 10 +offset, seg, de  );
	HF1( 14 +offset, peaktime ); HF2( 15 +offset, seg, peaktime );
	HF1( 17 +offset, chisqr   );
	HF2( 18 +offset, de, peaktime );
	HF2( 20 +offset, de, chisqr );
	for(size_t m=0; m<waveform.size(); ++m){
	  HF2( 16 +offset, time[m], waveform[m] );
	}
      }

      if( !good ) continue;
      int nh = nhSsdIn[layer-1][Good];
      switch(layer){
      case 1:
	event.ssd1y0ped[nh] = ped; event.ssd1y0hitpat[nh] = seg;
	event.ssd1y0pos[nh] = pos; event.ssd1y0adc[nh]    = adc;
	event.ssd1y0tdc[nh] = tdc;
	event.ssd1y0de[nh]  = de;
	event.ssd1y0time[nh]   = peaktime;
	event.ssd1y0chisqr[nh] = chisqr;
	break;
      case 2:
	event.ssd1x0ped[nh] = ped; event.ssd1x0hitpat[nh] = seg;
	event.ssd1x0pos[nh] = pos; event.ssd1x0adc[nh]    = adc;
	event.ssd1x0tdc[nh] = tdc;
	event.ssd1x0sum[nh] = sum;
	event.ssd1x0de[nh]  = de;
	event.ssd1x0time[nh]   = peaktime;
	event.ssd1x0chisqr[nh] = chisqr;
	break;
      default:
	break;
      }
      nhSsdIn[layer-1][Good]++;
    }
  }

  event.ssd1y0nhits = nhSsdIn[0][Good];
  event.ssd1x0nhits = nhSsdIn[1][Good];

  for( int l=0; l<NumOfLayersSSD1; ++l ){
    for( int s=0; s<nStatus; ++s ){
      HF1( SSD1Hid +1000*(l+1) +s*100 +0, nhSsdIn[l][s] );
    }
  }

  HF1( 1, 12. );

  ///// SsdInCluster
  for( int layer=1; layer<=NumOfLayersSSD1; ++layer ){
    const SsdClusterContainer &ClCont = DCAna->GetClusterSsdIn(layer);
    int ncl = DCAna->GetNClustersSsdIn(layer);
    switch(layer){
    case 1: event.ssd1y0ncl = ncl; break;
    case 2: event.ssd1x0ncl = ncl; break;
    default: break;
    }
    int offset = SSD1ClHid + layer*1000;
    HF1( 0 +offset, ncl );
    for( int i=0; i<ncl; ++i ){
      SsdCluster *cl = ClCont[i];
      if( !cl ) continue;
      if( !cl->GoodForAnalysis() ) continue;
      int    clsize = cl->ClusterSize();
      double seg    = cl->MeanSeg();
      double pos    = cl->Position();
      double de     = cl->DeltaE();
      double time   = cl->Time();
      HF1(  1 +offset, seg ); HF1(  2 +offset, pos );
      HF1(  9 +offset, de );
      HF2( 10 +offset, seg, de ); HF1( 14 +offset, time );
      HF2( 15 +offset, seg, time ); HF2( 18 +offset, de, time );
      HF1( 30 +offset, clsize );
      if( clsize<5  ) HF1( 30 +offset +clsize, de );
      if( clsize>=5 ) HF1( 35 +offset, de );

      switch(layer){
      case 1:
	event.ssd1y0clsize[i] = clsize;
	event.ssd1y0clpos[i]  = pos;
	event.ssd1y0clde[i]   = de;
	event.ssd1y0cltime[i] = time;
	break;
      case 2:
	event.ssd1x0clsize[i] = clsize;
	event.ssd1x0clpos[i]  = pos;
	event.ssd1x0clde[i]   = de;
	event.ssd1x0cltime[i] = time;
	break;
      default:
	break;
      }
    }
  }

  ///// SsdOut
  double multiSsdOut[nStatus] = {};
  int nhSsdOut[NumOfLayersSSD2][nStatus];
  for( int l=0; l<NumOfLayersSSD2; ++l ){
    for( int s=0; s<nStatus; ++s ){
      nhSsdOut[l][s] = 0;
    }
  }

  for( int layer=1; layer<=NumOfLayersSSD2; ++layer ){
    const DCHitContainer &contOut = DCAna->GetSsdOutHC(layer);
    int nhOut = contOut.size();
    multiSsdOut[All] += double(nhOut);
    nhSsdOut[layer-1][All] = nhOut;
    for( int i=0; i<nhOut; ++i ){
      DCHit  *hit  = contOut[i];
      if(!hit) continue;
      int    seg  = hit->GetWire();
      bool   flag = hit->IsZeroSuppressed();
      bool   good = hit->IsGoodWaveForm();
      double ped  = hit->GetPedestal();
      double adc  = hit->GetAdcPeakHeight() - ped;
      double tdc  = hit->GetAdcPeakPosition();
      double pos  = hit->GetWirePosition();
      double de  = hit->GetDe();
      double peaktime = hit->GetPeakTime(), chisqr = hit->GetChisquare();
      std::vector<double> waveform = hit->GetWaveform();
      std::vector<double> time     = hit->GetTime();
      if( !flag ) continue;
      if( good )  multiSsdOut[Good]++;

      for(int status=0; status<good+1; ++status){
	int offset = SSD2Hid + layer*1000 + status*100;
	HF1(  1 +offset, seg      ); HF1(  2 +offset, pos      );
	HF1(  3 +offset, adc      );
	HF2(  4 +offset, seg, adc );
	HF1(  5 +offset, tdc      ); HF2(  6 +offset, seg, tdc );
	HF2(  8 +offset, seg, ped );
	HF1(  9 +offset, de       ); HF2( 10 +offset, seg, de  );
	HF1( 14 +offset, peaktime ); HF2( 15 +offset, seg, peaktime );
	HF1( 17 +offset, chisqr   );
	HF2( 18 +offset, de, peaktime );
	HF2( 20 +offset, de, chisqr );
	for( size_t m=0; m<waveform.size(); ++m ){
	  HF2( 16 +offset, time[m], waveform[m] );
	}
      }

      if( !good ) continue;
      int nh = nhSsdOut[layer-1][Good];
      switch(layer){
      case 1:
	event.ssd2x0ped[nh] = ped; event.ssd2x0hitpat[nh] = seg;
	event.ssd2x0pos[nh] = pos; event.ssd2x0adc[nh]    = adc;
	event.ssd2x0tdc[nh] = tdc;
	event.ssd2x0de[nh]  = de;
	event.ssd2x0time[nh]   = peaktime;
	event.ssd2x0chisqr[nh] = chisqr;
	break;
      case 2:
	event.ssd2y0ped[nh] = ped; event.ssd2y0hitpat[nh] = seg;
	event.ssd2y0pos[nh] = pos; event.ssd2y0adc[nh]    = adc;
	event.ssd2y0tdc[nh] = tdc;
	event.ssd2y0de[nh]  = de;
	event.ssd2y0time[nh]   = peaktime;
	event.ssd2y0chisqr[nh] = chisqr;
	break;
      }
      nhSsdOut[layer-1][Good]++;
    }
  }

  event.ssd2x0nhits = nhSsdOut[0][Good];
  event.ssd2y0nhits = nhSsdOut[1][Good];

  for( int l=0; l<NumOfLayersSSD2; ++l ){
    for( int s=0; s<nStatus; ++s ){
      HF1( SSD2Hid +1000*(l+1) +s*100 +0, nhSsdOut[l][s] );
    }
  }

  HF1( 1, 13. );

  ///// SsdOutCluster
  for( int layer=1; layer<=NumOfLayersSSD2; ++layer ){
    const SsdClusterContainer &ClCont = DCAna->GetClusterSsdOut(layer);
    int ncl = DCAna->GetNClustersSsdOut(layer);
    switch(layer){
    case 1: event.ssd2x0ncl = ncl; break;
    case 2: event.ssd2y0ncl = ncl; break;
    default: break;
    }
    int offset = SSD2ClHid + layer*1000;
    HF1( 0 +offset, ncl );
    for( int i=0; i<ncl; ++i ){
      SsdCluster *cl = ClCont[i];
      if( !cl ) continue;
      if( !cl->GoodForAnalysis() ) continue;
      int    clsize = cl->ClusterSize();
      double seg    = cl->MeanSeg();
      double pos    = cl->Position();
      double de     = cl->DeltaE();
      double time   = cl->Time();
      HF1(  1 +offset, seg ); HF1(  2 +offset, pos );
      HF1(  9 +offset, de );
      HF2( 10 +offset, seg, de ); HF1( 14 +offset, time );
      HF2( 15 +offset, seg, time ); HF2( 18 +offset, de, time );
      HF1( 30 +offset, clsize );
      if( clsize<5  ) HF1( 30 +offset +clsize, de );
      if( clsize>=5 ) HF1( 35 +offset, de );

      switch(layer){
      case 1:
	event.ssd2x0clsize[i] = clsize;
	event.ssd2x0clpos[i]  = pos;
	event.ssd2x0clde[i]   = de;
	event.ssd2x0cltime[i] = time;
	break;
      case 2:
	event.ssd2y0clsize[i] = clsize;
	event.ssd2y0clpos[i]  = pos;
	event.ssd2y0clde[i]   = de;
	event.ssd2y0cltime[i] = time;
	break;
      default:
	break;
      }
    }
  }

  int ssd1nlayers[nStatus];
  int ssd2nlayers[nStatus];
  for( int s=0; s<nStatus; ++s ){
    HF1( 1000 + s*100 +0, multiSsdIn[s] );
    HF1( 2000 + s*100 +0, multiSsdOut[s] );
    ssd1nlayers[s]
      = (nhSsdIn[0][s]>0) + (nhSsdIn[1][s]>0);
    ssd2nlayers[s]
      = (nhSsdOut[0][s]>0) + (nhSsdOut[1][s]>0);
    HF1( 3000 +s*100 +0, ssd1nlayers[s] );
    HF1( 4000 +s*100 +0, ssd2nlayers[s] );
  }

  HF1( 1, 19. );

  return true;
}

//______________________________________________________________________________
void
EventHimac::InitializeEvent( void )
{
  event.evnum       = 0;
  event.trignhits   = 0;

  event.t1nhits	    = 0;
  for(int it=0; it<NumOfSegT1; it++){
    event.t1hitpat[it]    = 0;
    event.t1ua[it]    = -9999;
    event.t1ut[it]    = -9999;
    event.t1da[it]    = -9999;
    event.t1dt[it]    = -9999;
  }
  event.t2nhits	    = 0;
  for(int it=0; it<NumOfSegT2; it++){
    event.t2hitpat[it]    = 0;
    event.t2ua[it]    = -9999;
    event.t2ut[it]    = -9999;
    event.t2da[it]    = -9999;
    event.t2dt[it]    = -9999;
  }
  event.t3nhits	    = 0;
  for(int it=0; it<NumOfSegT3; it++){
    event.t3hitpat[it]    = 0;
    event.t3ua[it]    = -9999;
    event.t3ut[it]    = -9999;
    event.t3da[it]    = -9999;
    event.t3dt[it]    = -9999;
  }
  event.t4nhits	    = 0;
  for(int it=0; it<NumOfSegT4; it++){
    event.t4hitpat[it]    = 0;
    event.t4ua[it]    = -9999;
    event.t4ut[it]    = -9999;
    event.t4da[it]    = -9999;
    event.t4dt[it]    = -9999;
  }
  event.s1nhits	    = 0;
  for(int it=0; it<NumOfSegS1; it++){
    event.s1hitpat[it]    = 0;
    event.s1ua[it]    = -9999;
    event.s1ut[it]    = -9999;
    event.s1da[it]    = -9999;
    event.s1dt[it]    = -9999;
  }
  event.s2nhits	    = 0;
  for(int it=0; it<NumOfSegS2; it++){
    event.s2hitpat[it]    = 0;
    event.s2ua[it]    = -9999;
    event.s2ut[it]    = -9999;
    event.s2da[it]    = -9999;
    event.s2dt[it]    = -9999;
  }

  event.ssdtnhits   = 0;

  event.ssd1y0nhits = 0;
  event.ssd1x0nhits = 0;

  event.ssd2x0nhits = 0;
  event.ssd2y0nhits = 0;

  event.ssd1y0ncl = 0;
  event.ssd1x0ncl = 0;

  event.ssd2x0ncl = 0;
  event.ssd2y0ncl = 0;

  for(int it=0; it<NumOfSegTrig; it++){
    event.trigpat[it]  = -1;
    event.trigflag[it] = -9999;
  }

  for(int it=0; it<NumOfSegSSDT; it++){
    event.ssdthitpat[it] = -1;
    event.ssdtdc[it]  = -9999;
    event.ssdtime[it] = -9999;
  }

  for(int it=0; it<NumOfSegSSD1; it++){
    event.ssd1y0hitpat[it] = -1;
    event.ssd1y0pos[it]  = -9999;
    event.ssd1y0ped[it]  = -9999;
    event.ssd1y0adc[it]  = -9999;
    event.ssd1y0tdc[it]  = -9999;
    event.ssd1y0de[it]   = -9999;
    event.ssd1y0time[it] = -9999;
    event.ssd1y0chisqr[it] = -1.;

    event.ssd1x0hitpat[it] = -1;
    event.ssd1x0pos[it]  = -9999;
    event.ssd1x0ped[it]  = -9999;
    event.ssd1x0adc[it]  = -9999;
    event.ssd1x0tdc[it]  = -9999;
    event.ssd1x0sum[it]  = -9999;
    event.ssd1x0de[it]   = -9999;
    event.ssd1x0time[it] = -9999;
    event.ssd1x0chisqr[it] = -1.;

    event.ssd1y0clsize[it] = 0;
    event.ssd1y0clpos[it]  = -9999;
    event.ssd1y0clde[it]   = -9999;
    event.ssd1y0cltime[it] = -9999;
    event.ssd1y0cltd[it]   = -9999;
    event.ssd1x0clsize[it] = 0;
    event.ssd1x0clpos[it]  = -9999;
    event.ssd1x0clde[it]   = -9999;
    event.ssd1x0cltime[it] = -9999;
    event.ssd1x0cltd[it]   = -9999;
  }

  for(int it=0; it<NumOfSegSSD2; it++){
    event.ssd2x0hitpat[it] = -1;
    event.ssd2x0pos[it]  = -9999;
    event.ssd2x0ped[it]  = -9999;
    event.ssd2x0adc[it]  = -9999;
    event.ssd2x0tdc[it]  = -9999;
    event.ssd2x0de[it]   = -9999;
    event.ssd2x0time[it] = -9999;
    event.ssd2x0chisqr[it] = -1.;

    event.ssd2y0hitpat[it] = -1;
    event.ssd2y0pos[it]  = -9999;
    event.ssd2y0ped[it]  = -9999;
    event.ssd2y0adc[it]  = -9999;
    event.ssd2y0tdc[it]  = -9999;
    event.ssd2y0de[it]   = -9999;
    event.ssd2y0time[it] = -9999;
    event.ssd2y0chisqr[it] = -1.;

    event.ssd2x0clsize[it] = 0;
    event.ssd2x0clpos[it]  = -9999;
    event.ssd2x0clde[it]   = -9999;
    event.ssd2x0cltime[it] = -9999;
    event.ssd2x0cltd[it]   = -9999;
    event.ssd2y0clsize[it] = 0;
    event.ssd2y0clpos[it]  = -9999;
    event.ssd2y0clde[it]   = -9999;
    event.ssd2y0cltime[it] = -9999;
    event.ssd2y0cltd[it]   = -9999;
  }
}

//______________________________________________________________________________
bool
EventHimac::ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventHimac;
}

//______________________________________________________________________________
namespace
{
  const int NbinAdc   = 1000;
  const double MinAdc =    0.;
  const double MaxAdc = 1000.;

  const int NbinTdc   = 10;
  const double MinTdc =  0.;
  const double MaxTdc = 10.;

  const int NbinTime   =  300;
  const double MinTime =  -50.;
  const double MaxTime =  250.;

  const int NbinDe   =  510;
  const double MinDe = -1000.;
  const double MaxDe =  50000.;

  const int NbinX   = 300;
  const double MinX = -30.;
  const double MaxX =  30.;

  const int NbinY   = 240;
  const double MinY = -60.;
  const double MaxY =  60.;

  const int NbinChisqr   = 500;
  const double MinChisqr =   0.;
  const double MaxChisqr = 500.;

  const int NbinClsize   =  40;
  const double MinClsize =   0.;
  const double MaxClsize =  40.;

  const std::string SsdInDetName[NumOfLayersSSD1] =
    { "SSD1-X0", "SSD1-Y0" };
  const std::string SsdOutDetName[NumOfLayersSSD2] =
    { "SSD2-Y0", "SSD2-X0" };
  const std::string GoodInfo[nStatus] = { "", "[Good]" };
}

//______________________________________________________________________________
bool
ConfMan::InitializeHistograms( void )
{
  HB1( 1, "Status", 20, 0., 20. );

  HB1( 10, "Trigger HitPat", NumOfSegTrig, 0., double(NumOfSegTrig) );
  for( int i=0; i<NumOfSegTrig; ++i ){
    HB1( 10 +i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000 );
  }

  //////////////////// SSD-T
  HB1( SSDTHid +0, "#Hits SSDT",  NumOfSegSSDT, 0., (double)NumOfSegSSDT );
  HB1( SSDTHid +1, "HitPat SSDT", NumOfSegSSDT, 0., (double)NumOfSegSSDT );
  for( int i=0; i<NumOfSegSSDT; ++i ){
    HB1( SSDTHid +(i+1)*1000 + 5, Form("SSDT Tdc-%d", i+1),  0x1000, 0, 0x1000 );
    HB1( SSDTHid +(i+1)*1000 +14, Form("SSDT Time-%d", i+1), 200, -35, 35 );
  }
  HB2( SSDTHid + 6, "SSDT Tdc%Seg",
       NumOfSegSSDT, 0., (double)NumOfSegSSDT, 0x1000, 0, 0x1000 );
  HB2( SSDTHid +15,  "SSDT Time%Seg",
       NumOfSegSSDT, 0., (double)NumOfSegSSDT, 200, -35, 35 );

  for( int g=0; g<nStatus; ++g ){
    HB1( 1000 +g*100 +0, Form("SSD1 NHits %s", GoodInfo[g].c_str()), 50,  0., 50. );
    HB1( 2000 +g*100 +0, Form("SSD2 NHits %s", GoodInfo[g].c_str()), 50,  0., 50. );
    HB1( 3000 +g*100 +0, Form("SSD1 NLayers %s", GoodInfo[g].c_str()),
  	 NumOfLayersSSD1+1,  0., double(NumOfLayersSSD1+1)  );
    HB1( 4000 +g*100 +0, Form("SSD2 NLayers %s", GoodInfo[g].c_str()),
  	 NumOfLayersSSD2+1, 0., double(NumOfLayersSSD2+1) );
  }

  //////////////////// SsdIn
  for( int l=0; l<NumOfLayersSSD1; ++l ){
    const char* dname  = SsdInDetName[l].c_str();
    for( int g=0; g<nStatus; ++g ){
      int offset = SSD1Hid +(l+1)*1000 +g*100;
      const char* info = GoodInfo[g].c_str();
      HB1(  0 +offset, Form("%s #Hits %s",  dname, info), 100, 0., 100. );
      HB1(  1 +offset, Form("%s HitPat %s", dname, info),
  	    NumOfSegSSD1, 0., double(NumOfSegSSD1) );
      HB1(  2 +offset, Form("%s Position %s", dname, info), NbinX, MinX, MaxX );
      HB1(  3 +offset, Form("%s Adc %s", dname, info), NbinAdc, MinAdc, MaxAdc );
      HB2(  4 +offset, Form("%s Adc%%Seg %s", dname, info),
  	    NumOfSegSSD1, 0, double(NumOfSegSSD1), NbinAdc, MinAdc, MaxAdc );
      HB1(  5 +offset, Form("%s Tdc %s", dname, info), NbinTdc, MinTdc, MaxTdc );
      HB2(  6 +offset, Form("%s Tdc%%Seg %s", dname, info),
  	    NumOfSegSSD1, 0, double(NumOfSegSSD1), NbinTdc, MinTdc, MaxTdc );
      HB2(  8 +offset, Form("%s Pedestal%%Seg %s", dname, info),
  	    NumOfSegSSD1, 0, double(NumOfSegSSD1), NbinAdc, MinAdc, MaxAdc );
      HB1(  9 +offset, Form("%s dE %s", dname, info), NbinDe, MinDe, MaxDe );
      HB2( 10 +offset, Form("%s dE%%Seg %s", dname, info),
  	   NumOfSegSSD1, 0, double(NumOfSegSSD1), NbinDe, MinDe, MaxDe );
      HB1( 14 +offset, Form("%s Time %s", dname, info),
  	   NbinTime, MinTime, MaxTime );
      HB2( 15 +offset, Form("%s Time%%Seg %s", dname, info),
  	   NumOfSegSSD1, 0, double(NumOfSegSSD1), NbinTime, MinTime, MaxTime );
      HB2( 16 +offset, Form("%s Waveform %s", dname, info),
  	   NbinTime, MinTime, MaxTime, NbinAdc, MinAdc, MaxAdc );
      HB1( 17 +offset, Form("%s Chisqr %s", dname, info),
  	   NbinChisqr, MinChisqr, MaxChisqr);
      HB2( 18 +offset, Form("%s Time%%De %s", dname, info),
  	   NbinDe, MinDe, MaxDe, NbinTime, MinTime, MaxTime );
      HB2( 20 +offset, Form("%s Chisqr%%De %s", dname, info),
  	   NbinDe, MinDe, MaxDe, NbinChisqr, MinChisqr, MaxChisqr );
    }

    //////////////////// SsdIn (Cluster)
    int offset = SSD1ClHid +(l+1)*1000;
    HB1(  0 +offset, Form( "%s NCluster [Cluster]",  dname ), 100, 0., 100. );
    HB1(  1 +offset, Form( "%s HitPat [Cluster]", dname ),
	  NumOfSegSSD1, 0., double(NumOfSegSSD1) );
    HB1(  2 +offset, Form( "%s Position [Cluster]", dname ), NbinX, MinX, MaxX );
    HB1(  9 +offset, Form( "%s dE [Cluster]", dname ), NbinDe, MinDe, MaxDe );
    HB2( 10 +offset, Form("%s dE%%Seg [Cluster]", dname ),
	 NumOfSegSSD1, 0, double(NumOfSegSSD1), NbinDe, MinDe, MaxDe );
    HB1( 14 +offset, Form("%s Time [Cluster]", dname ),
	 NbinTime, MinTime, MaxTime );
    HB2( 15 +offset, Form("%s Time%%Seg [Cluster]", dname ),
	 NumOfSegSSD1, 0, double(NumOfSegSSD1), NbinTime, MinTime, MaxTime );
    HB2( 18 +offset, Form("%s Time%%De [Cluster]", dname ),
	 NbinDe, MinDe, MaxDe, NbinTime, MinTime, MaxTime );
    HB1( 30 +offset, Form("%s CluserSize [Cluster]", dname ),
	 NbinClsize, MinClsize, MaxClsize );
    HB1( 31 +offset, Form("%s dE [ClusterSize==1]", dname ), NbinDe, MinDe, MaxDe );
    HB1( 32 +offset, Form("%s dE [ClusterSize==2]", dname ), NbinDe, MinDe, MaxDe );
    HB1( 33 +offset, Form("%s dE [ClusterSize==3]", dname ), NbinDe, MinDe, MaxDe );
    HB1( 34 +offset, Form("%s dE [ClusterSize==4]", dname ), NbinDe, MinDe, MaxDe );
    HB1( 35 +offset, Form("%s dE [ClusterSize==5]", dname ), NbinDe, MinDe, MaxDe );
  }

  //////////////////// SsdOut
  for( int l=0; l<NumOfLayersSSD2; ++l ){
    const char* dname  = SsdOutDetName[l].c_str();
    for( int g=0; g<nStatus; ++g ){
      int offset = SSD2Hid +(l+1)*1000 +g*100;
      const char* info = GoodInfo[g].c_str();
      HB1(  0 +offset, Form("%s #Hits %s",  dname, info), 100, 0., 100. );
      HB1(  1 +offset, Form("%s HitPat %s", dname, info),
  	    NumOfSegSSD2, 0., double(NumOfSegSSD2) );
      HB1(  2 +offset, Form("%s Position %s", dname, info), NbinX, MinX, MaxX );
      HB1(  3 +offset, Form("%s Adc %s", dname, info), NbinAdc, MinAdc, MaxAdc );
      HB2(  4 +offset, Form("%s Adc%%Seg %s", dname, info),
  	    NumOfSegSSD2, 0, double(NumOfSegSSD2), NbinAdc, MinAdc, MaxAdc );
      HB1(  5 +offset, Form("%s Tdc %s", dname, info), NbinTdc, MinTdc, MaxTdc );
      HB2(  6 +offset, Form("%s Tdc%%Seg %s", dname, info),
  	    NumOfSegSSD2, 0, double(NumOfSegSSD2), NbinTdc, MinTdc, MaxTdc );
      HB2(  8 +offset, Form("%s Pedestal%%Seg %s", dname, info),
  	    NumOfSegSSD2, 0, double(NumOfSegSSD2), NbinAdc, MinAdc, MaxAdc );
      HB1(  9 +offset, Form("%s dE %s", dname, info), NbinDe, MinDe, MaxDe );
      HB2( 10 +offset, Form("%s dE%%Seg %s", dname, info),
  	   NumOfSegSSD2, 0, double(NumOfSegSSD2), NbinDe, MinDe, MaxDe );
      HB1( 14 +offset, Form("%s Time %s", dname, info),
  	   NbinTime, MinTime, MaxTime );
      HB2( 15 +offset, Form("%s Time%%Seg %s", dname, info),
  	   NumOfSegSSD2, 0, double(NumOfSegSSD2), NbinTime, MinTime, MaxTime );
      HB2( 16 +offset, Form("%s Waveform %s", dname, info),
  	   NbinTime, MinTime, MaxTime, NbinAdc, MinAdc, MaxAdc );
      HB1( 17 +offset, Form("%s Chisqr %s", dname, info),
  	   NbinChisqr, MinChisqr, MaxChisqr);
      HB2( 18 +offset, Form("%s Time%%De %s", dname, info),
  	   NbinDe, MinDe, MaxDe, NbinTime, MinTime, MaxTime );
      HB2( 20 +offset, Form("%s Chisqr%%De %s", dname, info),
  	   NbinDe, MinDe, MaxDe, NbinChisqr, MinChisqr, MaxChisqr );
    }
    //////////////////// SsdOut (Cluster)
    int offset = SSD2ClHid +(l+1)*1000;
    HB1(  0 +offset, Form( "%s NCluster [Cluster]",  dname ), 100, 0., 100. );
    HB1(  1 +offset, Form( "%s HitPat [Cluster]", dname ),
	  NumOfSegSSD2, 0., double(NumOfSegSSD2) );
    HB1(  2 +offset, Form( "%s Position [Cluster]", dname ), NbinX, MinX, MaxX );
    HB1(  9 +offset, Form( "%s dE [Cluster]", dname ), NbinDe, MinDe, MaxDe );
    HB2( 10 +offset, Form("%s dE%%Seg [Cluster]", dname ),
	 NumOfSegSSD2, 0, double(NumOfSegSSD2), NbinDe, MinDe, MaxDe );
    HB1( 14 +offset, Form("%s Time [Cluster]", dname ),
	 NbinTime, MinTime, MaxTime );
    HB2( 15 +offset, Form("%s Time%%Seg [Cluster]", dname ),
	 NumOfSegSSD2, 0, double(NumOfSegSSD2), NbinTime, MinTime, MaxTime );
    HB2( 18 +offset, Form("%s Time%%De [Cluster]", dname ),
	 NbinDe, MinDe, MaxDe, NbinTime, MinTime, MaxTime );
    HB1( 30 +offset, Form("%s CluserSize [Cluster]", dname ),
	 NbinClsize, MinClsize, MaxClsize );
    HB1( 31 +offset, Form("%s dE [ClusterSize==1]", dname ), NbinDe, MinDe, MaxDe );
    HB1( 32 +offset, Form("%s dE [ClusterSize==2]", dname ), NbinDe, MinDe, MaxDe );
    HB1( 33 +offset, Form("%s dE [ClusterSize==3]", dname ), NbinDe, MinDe, MaxDe );
    HB1( 34 +offset, Form("%s dE [ClusterSize==4]", dname ), NbinDe, MinDe, MaxDe );
    HB1( 35 +offset, Form("%s dE [ClusterSize==5]", dname ), NbinDe, MinDe, MaxDe );
  }

  ////////////////////////////////////////////
  //Tree
  HBTree( "tree", "tree of Himac" );
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trignhits", &event.trignhits, "trignhits/I");
  tree->Branch("trigpat",    event.trigpat,   "trigpat[trignhits]/I");
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));
  
  tree->Branch("unixtime",  &event.unixtime,  "unixtime/l");

  // T1
  tree->Branch("t1nhits",   &event.t1nhits,  "t1nhits/I");
  tree->Branch("t1hitpat",   event.t1hitpat,  Form("t1hitpat[%d]/I", NumOfSegT1));
  tree->Branch("t1ua", 	     event.t1ua,      Form("t1ua[%d]/D", NumOfSegT1));
  tree->Branch("t1ut",       event.t1ut,      Form("t1ut[%d]/D", NumOfSegT1));
  tree->Branch("t1da",       event.t1da,      Form("t1da[%d]/D", NumOfSegT1));
  tree->Branch("t1dt",       event.t1dt,      Form("t1dt[%d]/D", NumOfSegT1));

  // T2
  tree->Branch("t2nhits",   &event.t2nhits,  "t2nhits/I");
  tree->Branch("t2hitpat",   event.t2hitpat,  Form("t2hitpat[%d]/I", NumOfSegT2));
  tree->Branch("t2ua", 	     event.t2ua,      Form("t2ua[%d]/D", NumOfSegT2));
  tree->Branch("t2ut",       event.t2ut,      Form("t2ut[%d]/D", NumOfSegT2));
  tree->Branch("t2da",       event.t2da,      Form("t2da[%d]/D", NumOfSegT2));
  tree->Branch("t2dt",       event.t2dt,      Form("t2dt[%d]/D", NumOfSegT2));

  // T3
  tree->Branch("t3nhits",   &event.t3nhits,  "t3nhits/I");
  tree->Branch("t3hitpat",   event.t3hitpat,  Form("t3hitpat[%d]/I", NumOfSegT3));
  tree->Branch("t3ua", 	     event.t3ua,      Form("t3ua[%d]/D", NumOfSegT3));
  tree->Branch("t3ut",       event.t3ut,      Form("t3ut[%d]/D", NumOfSegT3));
  tree->Branch("t3da",       event.t3da,      Form("t3da[%d]/D", NumOfSegT3));
  tree->Branch("t3dt",       event.t3dt,      Form("t3dt[%d]/D", NumOfSegT3));

  // T4
  tree->Branch("t4nhits",   &event.t4nhits,  "t4nhits/I");
  tree->Branch("t4hitpat",   event.t4hitpat,  Form("t4hitpat[%d]/I", NumOfSegT4));
  tree->Branch("t4ua", 	     event.t4ua,      Form("t4ua[%d]/D", NumOfSegT4));
  tree->Branch("t4ut",       event.t4ut,      Form("t4ut[%d]/D", NumOfSegT4));
  tree->Branch("t4da",       event.t4da,      Form("t4da[%d]/D", NumOfSegT4));
  tree->Branch("t4dt",       event.t4dt,      Form("t4dt[%d]/D", NumOfSegT4));

  // S1
  tree->Branch("s1nhits",   &event.s1nhits,  "s1nhits/I");
  tree->Branch("s1hitpat",   event.s1hitpat,  Form("s1hitpat[%d]/I", NumOfSegS1));
  tree->Branch("s1ua", 	     event.s1ua,      Form("s1ua[%d]/D", NumOfSegS1));
  tree->Branch("s1ut",       event.s1ut,      Form("s1ut[%d]/D", NumOfSegS1));
  tree->Branch("s1da",       event.s1da,      Form("s1da[%d]/D", NumOfSegS1));
  tree->Branch("s1dt",       event.s1dt,      Form("s1dt[%d]/D", NumOfSegS1));

  // S2
  tree->Branch("s2nhits",   &event.s2nhits,  "s2nhits/I");
  tree->Branch("s2hitpat",   event.s2hitpat,  Form("s2hitpat[%d]/I", NumOfSegS2));
  tree->Branch("s2ua", 	     event.s2ua,      Form("s2ua[%d]/D", NumOfSegS2));
  tree->Branch("s2ut",       event.s2ut,      Form("s2ut[%d]/D", NumOfSegS2));
  tree->Branch("s2da",       event.s2da,      Form("s2da[%d]/D", NumOfSegS2));
  tree->Branch("s2dt",       event.s2dt,      Form("s2dt[%d]/D", NumOfSegS2));

  // SSD1-Y0
  tree->Branch("ssd1y0nhits",  &event.ssd1y0nhits,  "ssd1y0nhits/I");
  tree->Branch("ssd1y0hitpat",  event.ssd1y0hitpat, "ssd1y0hitpat[ssd1y0nhits]/D");
  tree->Branch("ssd1y0pos",     event.ssd1y0pos,    "ssd1y0pos[ssd1y0nhits]/D");
  tree->Branch("ssd1y0ped",     event.ssd1y0ped,    "ssd1y0ped[ssd1y0nhits]/D");
  tree->Branch("ssd1y0adc",     event.ssd1y0adc,    "ssd1y0adc[ssd1y0nhits]/D");
  tree->Branch("ssd1y0tdc",     event.ssd1y0tdc,    "ssd1y0tdc[ssd1y0nhits]/D");
  tree->Branch("ssd1y0de",      event.ssd1y0de,     "ssd1y0de[ssd1y0nhits]/D");
  tree->Branch("ssd1y0time",    event.ssd1y0time,   "ssd1y0time[ssd1y0nhits]/D");
  tree->Branch("ssd1y0chisqr",  event.ssd1y0chisqr, "ssd1y0chisqr[ssd1y0nhits]/D");
  // SSD1-X0
  tree->Branch("ssd1x0nhits",  &event.ssd1x0nhits,  "ssd1x0nhits/I");
  tree->Branch("ssd1x0hitpat",  event.ssd1x0hitpat, "ssd1x0hitpat[ssd1x0nhits]/D");
  tree->Branch("ssd1x0pos",     event.ssd1x0pos,    "ssd1x0pos[ssd1x0nhits]/D");
  tree->Branch("ssd1x0ped",     event.ssd1x0ped,    "ssd1x0ped[ssd1x0nhits]/D");
  tree->Branch("ssd1x0adc",     event.ssd1x0adc,    "ssd1x0adc[ssd1x0nhits]/D");
  tree->Branch("ssd1x0tdc",     event.ssd1x0tdc,    "ssd1x0tdc[ssd1x0nhits]/D");
  tree->Branch("ssd1x0sum",     event.ssd1x0sum,    "ssd1x0sum[ssd1x0nhits]/D");
  tree->Branch("ssd1x0de",      event.ssd1x0de,     "ssd1x0de[ssd1x0nhits]/D");
  tree->Branch("ssd1x0time",    event.ssd1x0time,   "ssd1x0time[ssd1x0nhits]/D");
  tree->Branch("ssd1x0chisqr",  event.ssd1x0chisqr, "ssd1x0chisqr[ssd1x0nhits]/D");

  // SSD2-X0
  tree->Branch("ssd2x0nhits",  &event.ssd2x0nhits,  "ssd2x0nhits/I");
  tree->Branch("ssd2x0hitpat",  event.ssd2x0hitpat, "ssd2x0hitpat[ssd2x0nhits]/D");
  tree->Branch("ssd2x0pos",     event.ssd2x0pos,    "ssd2x0pos[ssd2x0nhits]/D");
  tree->Branch("ssd2x0ped",     event.ssd2x0ped,    "ssd2x0ped[ssd2x0nhits]/D");
  tree->Branch("ssd2x0adc",     event.ssd2x0adc,    "ssd2x0adc[ssd2x0nhits]/D");
  tree->Branch("ssd2x0tdc",     event.ssd2x0tdc,    "ssd2x0tdc[ssd2x0nhits]/D");
  tree->Branch("ssd2x0de",      event.ssd2x0de,     "ssd2x0de[ssd2x0nhits]/D");
  tree->Branch("ssd2x0time",    event.ssd2x0time,   "ssd2x0time[ssd2x0nhits]/D");
  tree->Branch("ssd2x0chisqr",  event.ssd2x0chisqr, "ssd2x0chisqr[ssd2x0nhits]/D");
  // SSD2-Y0
  tree->Branch("ssd2y0nhits",  &event.ssd2y0nhits,  "ssd2y0nhits/I");
  tree->Branch("ssd2y0hitpat",  event.ssd2y0hitpat, "ssd2y0hitpat[ssd2y0nhits]/D");
  tree->Branch("ssd2y0pos",     event.ssd2y0pos,    "ssd2y0pos[ssd2y0nhits]/D");
  tree->Branch("ssd2y0ped",     event.ssd2y0ped,    "ssd2y0ped[ssd2y0nhits]/D");
  tree->Branch("ssd2y0adc",     event.ssd2y0adc,    "ssd2y0adc[ssd2y0nhits]/D");
  tree->Branch("ssd2y0tdc",     event.ssd2y0tdc,    "ssd2y0tdc[ssd2y0nhits]/D");
  tree->Branch("ssd2y0de",      event.ssd2y0de,     "ssd2y0de[ssd2y0nhits]/D");
  tree->Branch("ssd2y0time",    event.ssd2y0time,   "ssd2y0time[ssd2y0nhits]/D");
  tree->Branch("ssd2y0chisqr",  event.ssd2y0chisqr, "ssd2y0chisqr[ssd2y0nhits]/D");

  // SSD1-Y0 Cluster
  tree->Branch("ssd1y0ncl",    &event.ssd1y0ncl,    "ssd1y0ncl/I");
  tree->Branch("ssd1y0clsize",  event.ssd1y0clsize, "ssd1y0clsize[ssd1y0ncl]/I");
  tree->Branch("ssd1y0clpos",   event.ssd1y0clpos,  "ssd1y0clpos[ssd1y0ncl]/D");
  tree->Branch("ssd1y0clde",    event.ssd1y0clde,   "ssd1y0clde[ssd1y0ncl]/D");
  tree->Branch("ssd1y0cltime",  event.ssd1y0cltime, "ssd1y0cltime[ssd1y0ncl]/D");
  tree->Branch("ssd1y0cltd",    event.ssd1y0cltd,   "ssd1y0cltd[ssd1y0ncl]/D");
  // SSD1-X0 Cluster
  tree->Branch("ssd1x0ncl",    &event.ssd1x0ncl,    "ssd1x0ncl/I");
  tree->Branch("ssd1x0clsize",  event.ssd1x0clsize, "ssd1x0clsize[ssd1x0ncl]/I");
  tree->Branch("ssd1x0clpos",   event.ssd1x0clpos,  "ssd1x0clpos[ssd1x0ncl]/D");
  tree->Branch("ssd1x0clde",    event.ssd1x0clde,   "ssd1x0clde[ssd1x0ncl]/D");
  tree->Branch("ssd1x0cltime",  event.ssd1x0cltime, "ssd1x0cltime[ssd1x0ncl]/D");
  tree->Branch("ssd1x0cltd",    event.ssd1x0cltd,   "ssd1x0cltd[ssd1x0ncl]/D");

  // SSD2-X0 Cluster
  tree->Branch("ssd2x0ncl",    &event.ssd2x0ncl,    "ssd2x0ncl/I");
  tree->Branch("ssd2x0clsize",  event.ssd2x0clsize, "ssd2x0clsize[ssd2x0ncl]/I");
  tree->Branch("ssd2x0clpos",   event.ssd2x0clpos,  "ssd2x0clpos[ssd2x0ncl]/D");
  tree->Branch("ssd2x0clde",    event.ssd2x0clde,   "ssd2x0clde[ssd2x0ncl]/D");
  tree->Branch("ssd2x0cltime",  event.ssd2x0cltime, "ssd2x0cltime[ssd2x0ncl]/D");
  tree->Branch("ssd2x0cltd",    event.ssd2x0cltd,   "ssd2x0cltd[ssd2x0ncl]/D");
  // SSD2-Y0 Cluster
  tree->Branch("ssd2y0ncl",    &event.ssd2y0ncl,    "ssd2y0ncl/I");
  tree->Branch("ssd2y0clsize",  event.ssd2y0clsize, "ssd2y0clsize[ssd2y0ncl]/I");
  tree->Branch("ssd2y0clpos",   event.ssd2y0clpos,  "ssd2y0clpos[ssd2y0ncl]/D");
  tree->Branch("ssd2y0clde",    event.ssd2y0clde,   "ssd2y0clde[ssd2y0ncl]/D");
  tree->Branch("ssd2y0cltime",  event.ssd2y0cltime, "ssd2y0cltime[ssd2y0ncl]/D");
  tree->Branch("ssd2y0cltd",    event.ssd2y0cltd,   "ssd2y0cltd[ssd1y0ncl]/D");

  // SSDT
  tree->Branch("ssdtnhits",    &event.ssdtnhits,    "ssdtnhits/I");
  tree->Branch("ssdthitpat",    event.ssdthitpat,   "ssdthitpat[ssdtnhits]/D");
  tree->Branch("ssdtdc",        event.ssdtdc,       "ssdtdc[ssdtnhits]/D");
  tree->Branch("ssdtime",       event.ssdtime,      "ssdtime[ssdtnhits]/D");

  HPrint();

  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")        &&
      InitializeParameter<HodoParamMan>("HDPRM")     &&
      InitializeParameter<HodoPHCMan>("HDPHC")       &&
      //InitializeParameter<SsdParamMan>("SSDPRM")     &&
      InitializeParameter<UserParamMan>("USER")      );
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
