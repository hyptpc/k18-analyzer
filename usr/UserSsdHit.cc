/**
 *  file: UserSsdHit.cc
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
#define BcOutCut  0
#define SdcOutCut 0
#define MaxMultiCut 0

// SSD Filter
#define SlopeFilter     0
#define DeltaEFilter    1
#define TimeFilter      1
#define ChisqrFilter    1

namespace
{
  using namespace root;
  const std::string class_name("EventSsdHit");
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
class EventSsdHit : public VEvent
{
private:
  RawData      *rawData;
  DCAnalyzer   *DCAna;
  HodoAnalyzer *hodoAna;

public:
        EventSsdHit( void );
       ~EventSsdHit( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventSsdHit::EventSsdHit( void )
  : VEvent(),
    rawData(0),
    DCAna( new DCAnalyzer ),
    hodoAna( new HodoAnalyzer )
{
}

//______________________________________________________________________________
EventSsdHit::~EventSsdHit( void )
{
  if( DCAna )   delete DCAna;
  if( hodoAna ) delete hodoAna;
  if( rawData ) delete rawData;
}

//______________________________________________________________________________
bool
EventSsdHit::ProcessingBegin( void )
{
  InitializeEvent();
  return true;
}

//______________________________________________________________________________
struct Event
{
  int evnum;

  int trignhits;
  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

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

  int    ssd1y1nhits;
  double ssd1y1hitpat[NumOfSegSSD1];
  double ssd1y1pos[NumOfSegSSD1];
  double ssd1y1ped[NumOfSegSSD1];
  double ssd1y1adc[NumOfSegSSD1];
  double ssd1y1tdc[NumOfSegSSD1];
  double ssd1y1de[NumOfSegSSD1];
  double ssd1y1time[NumOfSegSSD1];
  double ssd1y1waveform[NumOfSegSSD1][NumOfSampleSSD];
  double ssd1y1chisqr[NumOfSegSSD1];

  int    ssd1x1nhits;
  double ssd1x1hitpat[NumOfSegSSD1];
  double ssd1x1pos[NumOfSegSSD1];
  double ssd1x1ped[NumOfSegSSD1];
  double ssd1x1adc[NumOfSegSSD1];
  double ssd1x1tdc[NumOfSegSSD1];
  double ssd1x1de[NumOfSegSSD1];
  double ssd1x1time[NumOfSegSSD1];
  double ssd1x1waveform[NumOfSegSSD1][NumOfSampleSSD];
  double ssd1x1chisqr[NumOfSegSSD1];

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

  int    ssd2x1nhits;
  double ssd2x1hitpat[NumOfSegSSD2];
  double ssd2x1pos[NumOfSegSSD2];
  double ssd2x1ped[NumOfSegSSD2];
  double ssd2x1adc[NumOfSegSSD2];
  double ssd2x1tdc[NumOfSegSSD2];
  double ssd2x1de[NumOfSegSSD2];
  double ssd2x1time[NumOfSegSSD2];
  double ssd2x1waveform[NumOfSegSSD2][NumOfSampleSSD];
  double ssd2x1chisqr[NumOfSegSSD2];

  int    ssd2y1nhits;
  double ssd2y1hitpat[NumOfSegSSD2];
  double ssd2y1pos[NumOfSegSSD2];
  double ssd2y1ped[NumOfSegSSD2];
  double ssd2y1adc[NumOfSegSSD2];
  double ssd2y1tdc[NumOfSegSSD2];
  double ssd2y1de[NumOfSegSSD2];
  double ssd2y1time[NumOfSegSSD2];
  double ssd2y1waveform[NumOfSegSSD2][NumOfSampleSSD];
  double ssd2y1chisqr[NumOfSegSSD2];

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
  int    ssd1y1ncl;
  int    ssd1y1clsize[NumOfSegSSD1];
  double ssd1y1clpos[NumOfSegSSD1];
  double ssd1y1clde[NumOfSegSSD1];
  double ssd1y1cltime[NumOfSegSSD1];
  double ssd1y1cltd[NumOfSegSSD1];
  int    ssd1x1ncl;
  int    ssd1x1clsize[NumOfSegSSD1];
  double ssd1x1clpos[NumOfSegSSD1];
  double ssd1x1clde[NumOfSegSSD1];
  double ssd1x1cltime[NumOfSegSSD1];
  double ssd1x1cltd[NumOfSegSSD1];

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
  int    ssd2x1ncl;
  int    ssd2x1clsize[NumOfSegSSD2];
  double ssd2x1clpos[NumOfSegSSD2];
  double ssd2x1clde[NumOfSegSSD2];
  double ssd2x1cltime[NumOfSegSSD2];
  double ssd2x1cltd[NumOfSegSSD2];
  int    ssd2y1ncl;
  int    ssd2y1clsize[NumOfSegSSD2];
  double ssd2y1clpos[NumOfSegSSD2];
  double ssd2y1clde[NumOfSegSSD2];
  double ssd2y1cltime[NumOfSegSSD2];
  double ssd2y1cltd[NumOfSegSSD2];

  // SSDT
  int    ssdtnhits;
  double ssdthitpat[NumOfSegSSDT];
  double ssdtdc[NumOfSegSSDT];
  double ssdtime[NumOfSegSSDT];

  // Tracking
  int ntBcOut;
  int ntSdcIn;
  int ntSdcOut;

};

//______________________________________________________________________________
namespace root
{
  Event  event;
  TH1   *h[MaxHist];
  TTree *tree;
  enum eDetHid  { SSD1Hid=10000, SSD2Hid=20000,
		  SSD1ClHid=30000, SSD2ClHid=40000,
		  SSDTHid=50000, nDetHid };
  enum eStatus  { All, Good, nStatus };
}

//______________________________________________________________________________
bool
EventSsdHit::ProcessingNormal( void )
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
#if MaxMultiCut
  static const double MaxMultiHitBcOut  = gUser.GetParameter("MaxMultiHitBcOut");
  static const double MaxMultiHitSdcOut = gUser.GetParameter("MaxMultiHitSdcOut");
#endif

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

  ////////// BH2 time 0
  hodoAna->DecodeBH2Hits(rawData);
  int nhBh2 = hodoAna->GetNHitsBH2();
#if HodoCut
  if( nhBh2==0 ) return true;
#endif
  HF1( 1, 2 );
  double time0 = -999.;
  ////////// BH2 Analysis
  for( int i=0; i<nhBh2; ++i ){
    BH2Hit *hit = hodoAna->GetHitBH2(i);
    if(!hit) continue;
    double cmt = hit->CMeanTime();
    double ct0 = hit->CTime0();
    double min_time = -999;
#if HodoCut
    double dE  = hit->DeltaE();
    if( dE<MinDeBH2 || MaxDeBH2<dE ) continue;
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
  if( nhBh1==0 ) return true;
#endif
  HF1(1, 4);

  double btof0 = -999;
  for(int i=0; i<nhBh1; ++i){
    Hodo2Hit* hit = hodoAna->GetHitBH1(i);
    if(!hit) continue;
    double btof = hit->CMeanTime() - time0;
#if HodoCut
    double dE   = hit->DeltaE();
    if( dE<MinDeBH1 || MaxDeBH1<dE ) continue;
#endif
    if( std::abs(btof)<std::abs(btof0) ){
      btof0 = btof;
    }
  }
  HF1( 100, btof0 );

  HF1( 1, 5. );

#if HodoCut
  if( MinBeamToF<btof && btof<MaxBeamToF )
    return true;
#endif

  HF1( 1, 6. );

  DCAna->DecodeBcOutHits( rawData );
  double multiBcOut = 0.;
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &cont = DCAna->GetBcOutHC(layer);
      multiBcOut += double( cont.size() );
    }
  }

#if MaxMultiCut
  if( multiBcOut/double(NumOfLayersBcOut) > MaxMultiHitBcOut )
    return true;
#endif

  HF1( 1, 7. );

  DCAna->TrackSearchBcOut();
  {
    int ntBcOut = DCAna->GetNtracksBcOut();
    event.ntBcOut = ntBcOut;
    int ntValid = 0;
    for( int it=0; it<ntBcOut; ++it ){
      DCLocalTrack *tp = DCAna->GetTrackBcOut(it);
      double chisqr = tp->GetChiSquare();
      if( chisqr<10. ) ++ntValid;
    }
    // if( ntValid==0 ) return true;
#if BcOutCut
    if( ntValid!=1 ) return true;
#endif
  }

  DCAna->DecodeSdcOutHits( rawData );
  double multiSdcOut = 0.;
  {
    for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
      const DCHitContainer &cont = DCAna->GetSdcOutHC(layer);
      multiSdcOut += double( cont.size() );
    }
  }

#if MaxMultiCut
  if( multiSdcOut/double(NumOfLayersSdcOut) > MaxMultiHitSdcOut )
    return true;
#endif

  HF1( 1, 8. );

  DCAna->TrackSearchSdcOut();
  {
    int ntSdcOut = DCAna->GetNtracksSdcOut();
    event.ntSdcOut = ntSdcOut;
    int ntValid = 0;
    for( int it=0; it<ntSdcOut; ++it ){
      DCLocalTrack *tp = DCAna->GetTrackSdcOut(it);
      double chisqr = tp->GetChiSquare();
      if( chisqr<20. ) ++ntValid;
    }
    // if( ntValid==0 ) return true;
#if SdcOutCut
    if( ntValid!=1 ) return true;
#endif
  }

  HF1( 1, 10. );

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
  int nhSsdIn[NumOfLayersSsdIn][nStatus];
  for( int l=0; l<NumOfLayersSsdIn; ++l ){
    for( int s=0; s<nStatus; ++s ){
      nhSsdIn[l][s] = 0;
    }
  }

  for( int layer=1; layer<=NumOfLayersSsdIn; ++layer ){
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
      case 3:
	event.ssd1y1ped[nh] = ped; event.ssd1y1hitpat[nh] = seg;
	event.ssd1y1pos[nh] = pos; event.ssd1y1adc[nh]    = adc;
	event.ssd1y1tdc[nh] = tdc;
	event.ssd1y1de[nh]  = de;
	event.ssd1y1time[nh]   = peaktime;
	event.ssd1y1chisqr[nh] = chisqr;
	break;
      case 4:
	event.ssd1x1ped[nh] = ped; event.ssd1x1hitpat[nh] = seg;
	event.ssd1x1pos[nh] = pos; event.ssd1x1adc[nh]    = adc;
	event.ssd1x1tdc[nh] = tdc;
	event.ssd1x1de[nh]  = de;
	event.ssd1x1time[nh]   = peaktime;
	event.ssd1x1chisqr[nh] = chisqr;
	break;
      default:
	break;
      }
      nhSsdIn[layer-1][Good]++;
    }
  }

  event.ssd1y0nhits = nhSsdIn[0][Good];
  event.ssd1x0nhits = nhSsdIn[1][Good];
  event.ssd1y1nhits = nhSsdIn[2][Good];
  event.ssd1x1nhits = nhSsdIn[3][Good];

  for( int l=0; l<NumOfLayersSsdIn; ++l ){
    for( int s=0; s<nStatus; ++s ){
      HF1( SSD1Hid +1000*(l+1) +s*100 +0, nhSsdIn[l][s] );
    }
  }

  HF1( 1, 12. );

  ///// SsdInCluster
  for( int layer=1; layer<=NumOfLayersSsdIn; ++layer ){
    const SsdClusterContainer &ClCont = DCAna->GetClusterSsdIn(layer);
    int ncl = DCAna->GetNClustersSsdIn(layer);
    switch(layer){
    case 1: event.ssd1y0ncl = ncl; break;
    case 2: event.ssd1x0ncl = ncl; break;
    case 3: event.ssd1y1ncl = ncl; break;
    case 4: event.ssd1x1ncl = ncl; break;
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
      case 3:
	event.ssd1y1clsize[i] = clsize;
	event.ssd1y1clpos[i]  = pos;
	event.ssd1y1clde[i]   = de;
	event.ssd1y1cltime[i] = time;
	break;
      case 4:
	event.ssd1x1clsize[i] = clsize;
	event.ssd1x1clpos[i]  = pos;
	event.ssd1x1clde[i]   = de;
	event.ssd1x1cltime[i] = time;
	break;
      default:
	break;
      }
    }
  }

  ///// SsdOut
  double multiSsdOut[nStatus] = {};
  int nhSsdOut[NumOfLayersSsdOut][nStatus];
  for( int l=0; l<NumOfLayersSsdOut; ++l ){
    for( int s=0; s<nStatus; ++s ){
      nhSsdOut[l][s] = 0;
    }
  }

  for( int layer=1; layer<=NumOfLayersSsdOut; ++layer ){
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
      case 3:
	event.ssd2x1ped[nh] = ped; event.ssd2x1hitpat[nh] = seg;
	event.ssd2x1pos[nh] = pos; event.ssd2x1adc[nh]    = adc;
	event.ssd2x1tdc[nh] = tdc;
	event.ssd2x1de[nh]  = de;
	event.ssd2x1time[nh]   = peaktime;
	event.ssd2x1chisqr[nh] = chisqr;
	break;
      case 4:
	event.ssd2y1ped[nh] = ped; event.ssd2y1hitpat[nh] = seg;
	event.ssd2y1pos[nh] = pos; event.ssd2y1adc[nh]    = adc;
	event.ssd2y1tdc[nh] = tdc;
	event.ssd2y1de[nh]  = de;
	event.ssd2y1time[nh]   = peaktime;
	event.ssd2y1chisqr[nh] = chisqr;
	break;
      }
      nhSsdOut[layer-1][Good]++;
    }
  }

  event.ssd2x0nhits = nhSsdOut[0][Good];
  event.ssd2y0nhits = nhSsdOut[1][Good];
  event.ssd2x1nhits = nhSsdOut[2][Good];
  event.ssd2y1nhits = nhSsdOut[3][Good];

  for( int l=0; l<NumOfLayersSsdOut; ++l ){
    for( int s=0; s<nStatus; ++s ){
      HF1( SSD2Hid +1000*(l+1) +s*100 +0, nhSsdOut[l][s] );
    }
  }

  HF1( 1, 13. );

  ///// SsdOutCluster
  for( int layer=1; layer<=NumOfLayersSsdOut; ++layer ){
    const SsdClusterContainer &ClCont = DCAna->GetClusterSsdOut(layer);
    int ncl = DCAna->GetNClustersSsdOut(layer);
    switch(layer){
    case 1: event.ssd2x0ncl = ncl; break;
    case 2: event.ssd2y0ncl = ncl; break;
    case 3: event.ssd2x1ncl = ncl; break;
    case 4: event.ssd2y1ncl = ncl; break;
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
      case 3:
	event.ssd2x1clsize[i] = clsize;
	event.ssd2x1clpos[i]  = pos;
	event.ssd2x1clde[i]   = de;
	event.ssd2x1cltime[i] = time;
	break;
      case 4:
	event.ssd2y1clsize[i] = clsize;
	event.ssd2y1clpos[i]  = pos;
	event.ssd2y1clde[i]   = de;
	event.ssd2y1cltime[i] = time;
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
      = (nhSsdIn[0][s]>0) + (nhSsdIn[1][s]>0) +
      (nhSsdIn[2][s]>0) + (nhSsdIn[3][s]>0);
    ssd2nlayers[s]
      = (nhSsdOut[0][s]>0) + (nhSsdOut[1][s]>0) +
      (nhSsdOut[2][s]>0) + (nhSsdOut[3][s]>0);
    HF1( 3000 +s*100 +0, ssd1nlayers[s] );
    HF1( 4000 +s*100 +0, ssd2nlayers[s] );
  }

  HF1( 1, 19. );

  return true;
}

//______________________________________________________________________________
void
EventSsdHit::InitializeEvent( void )
{
  event.evnum       = 0;
  event.trignhits   = 0;
  event.ssdtnhits   = 0;

  event.ssd1y0nhits = 0;
  event.ssd1x0nhits = 0;
  event.ssd1y1nhits = 0;
  event.ssd1x1nhits = 0;

  event.ssd2x0nhits = 0;
  event.ssd2y0nhits = 0;
  event.ssd2x1nhits = 0;
  event.ssd2y1nhits = 0;

  event.ssd1y0ncl = 0;
  event.ssd1x0ncl = 0;
  event.ssd1y1ncl = 0;
  event.ssd1x1ncl = 0;

  event.ssd2x0ncl = 0;
  event.ssd2y0ncl = 0;
  event.ssd2x1ncl = 0;
  event.ssd2y1ncl = 0;

  event.ntBcOut  = 0;
  event.ntSdcIn  = 0;
  event.ntSdcOut = 0;

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

    event.ssd1y1hitpat[it] = -1;
    event.ssd1y1pos[it]  = -9999;
    event.ssd1y1ped[it]  = -9999;
    event.ssd1y1adc[it]  = -9999;
    event.ssd1y1tdc[it]  = -9999;
    event.ssd1y1de[it]   = -9999;
    event.ssd1y1time[it] = -9999;
    event.ssd1y1chisqr[it] = -1.;

    event.ssd1x1hitpat[it] = -1;
    event.ssd1y1pos[it]  = -9999;
    event.ssd1x1ped[it]  = -9999;
    event.ssd1x1adc[it]  = -9999;
    event.ssd1x1tdc[it]  = -9999;
    event.ssd1y1de[it]   = -9999;
    event.ssd1y1time[it] = -9999;
    event.ssd1y1chisqr[it] = -1.;

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
    event.ssd1y1clsize[it] = 0;
    event.ssd1y1clpos[it]  = -9999;
    event.ssd1y1clde[it]   = -9999;
    event.ssd1y1cltime[it] = -9999;
    event.ssd1y1cltd[it]   = -9999;
    event.ssd1x1clsize[it] = 0;
    event.ssd1x1clpos[it]  = -9999;
    event.ssd1x1clde[it]   = -9999;
    event.ssd1x1cltime[it] = -9999;
    event.ssd1x1cltd[it]   = -9999;
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

    event.ssd2x1hitpat[it] = -1;
    event.ssd2x1pos[it]  = -9999;
    event.ssd2x1ped[it]  = -9999;
    event.ssd2x1adc[it]  = -9999;
    event.ssd2x1tdc[it]  = -9999;
    event.ssd2x1de[it]   = -9999;
    event.ssd2x1time[it] = -9999;
    event.ssd2x1chisqr[it] = -1.;

    event.ssd2y1hitpat[it] = -1;
    event.ssd2y1pos[it]  = -9999;
    event.ssd2y1ped[it]  = -9999;
    event.ssd2y1adc[it]  = -9999;
    event.ssd2y1tdc[it]  = -9999;
    event.ssd2y1de[it]   = -9999;
    event.ssd2y1time[it] = -9999;
    event.ssd2y1chisqr[it] = -1.;

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
    event.ssd2x1clsize[it] = 0;
    event.ssd2x1clpos[it]  = -9999;
    event.ssd2x1clde[it]   = -9999;
    event.ssd2x1cltime[it] = -9999;
    event.ssd2x1cltd[it]   = -9999;
    event.ssd2y1clsize[it] = 0;
    event.ssd2y1clpos[it]  = -9999;
    event.ssd2y1clde[it]   = -9999;
    event.ssd2y1cltime[it] = -9999;
    event.ssd2y1cltd[it]   = -9999;
  }
}

//______________________________________________________________________________
bool
EventSsdHit::ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventSsdHit;
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

  const std::string SsdInDetName[NumOfLayersSsdIn] =
    { "SSD1-Y0", "SSD1-X0", "SSD1-Y1", "SSD1-X1" };
  const std::string SsdOutDetName[NumOfLayersSsdOut] =
    { "SSD2-X0", "SSD2-Y0", "SSD2-X1", "SSD2-Y1" };
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

  HB1( 100, "Beam ToF", 200,  -10., 5. );

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
  	 NumOfLayersSsdIn+1,  0., double(NumOfLayersSsdIn+1)  );
    HB1( 4000 +g*100 +0, Form("SSD2 NLayers %s", GoodInfo[g].c_str()),
  	 NumOfLayersSsdOut+1, 0., double(NumOfLayersSsdOut+1) );
  }

  //////////////////// SsdIn
  for( int l=0; l<NumOfLayersSsdIn; ++l ){
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
  for( int l=0; l<NumOfLayersSsdOut; ++l ){
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
  HBTree( "ssd", "tree of SsdHit" );
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trignhits", &event.trignhits, "trignhits/I");
  tree->Branch("trigpat",    event.trigpat,   "trigpat[trignhits]/I");
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

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
  // SSD1-Y1
  tree->Branch("ssd1y1nhits",  &event.ssd1y1nhits,  "ssd1y1nhits/I");
  tree->Branch("ssd1y1hitpat",  event.ssd1y1hitpat, "ssd1y1hitpat[ssd1y1nhits]/D");
  tree->Branch("ssd1y1pos",     event.ssd1y1pos,    "ssd1y1pos[ssd1y1nhits]/D");
  tree->Branch("ssd1y1ped",     event.ssd1y1ped,    "ssd1y1ped[ssd1y1nhits]/D");
  tree->Branch("ssd1y1adc",     event.ssd1y1adc,    "ssd1y1adc[ssd1y1nhits]/D");
  tree->Branch("ssd1y1tdc",     event.ssd1y1tdc,    "ssd1y1tdc[ssd1y1nhits]/D");
  tree->Branch("ssd1y1de",      event.ssd1y1de,     "ssd1y1de[ssd1y1nhits]/D");
  tree->Branch("ssd1y1time",    event.ssd1y1time,   "ssd1y1time[ssd1y1nhits]/D");
  tree->Branch("ssd1y1chisqr",  event.ssd1y1chisqr, "ssd1y1chisqr[ssd1y1nhits]/D");
  // SSD1-X1
  tree->Branch("ssd1x1nhits",  &event.ssd1x1nhits,  "ssd1x1nhits/I");
  tree->Branch("ssd1x1hitpat",  event.ssd1x1hitpat, "ssd1x1hitpat[ssd1x1nhits]/D");
  tree->Branch("ssd1x1pos",     event.ssd1x1pos,    "ssd1x1pos[ssd1x1nhits]/D");
  tree->Branch("ssd1x1ped",     event.ssd1x1ped,    "ssd1x1ped[ssd1x1nhits]/D");
  tree->Branch("ssd1x1adc",     event.ssd1x1adc,    "ssd1x1adc[ssd1x1nhits]/D");
  tree->Branch("ssd1x1tdc",     event.ssd1x1tdc,    "ssd1x1tdc[ssd1x1nhits]/D");
  tree->Branch("ssd1x1de",      event.ssd1x1de,     "ssd1x1de[ssd1x1nhits]/D");
  tree->Branch("ssd1x1time",    event.ssd1x1time,   "ssd1x1time[ssd1x1nhits]/D");
  tree->Branch("ssd1x1chisqr",  event.ssd1x1chisqr, "ssd1x1chisqr[ssd1x1nhits]/D");

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
  // SSD2-X1
  tree->Branch("ssd2x1nhits",  &event.ssd2x1nhits,  "ssd2x1nhits/I");
  tree->Branch("ssd2x1hitpat",  event.ssd2x1hitpat, "ssd2x1hitpat[ssd2x1nhits]/D");
  tree->Branch("ssd2x1pos",     event.ssd2x1pos,    "ssd2x1pos[ssd2x1nhits]/D");
  tree->Branch("ssd2x1ped",     event.ssd2x1ped,    "ssd2x1ped[ssd2x1nhits]/D");
  tree->Branch("ssd2x1adc",     event.ssd2x1adc,    "ssd2x1adc[ssd2x1nhits]/D");
  tree->Branch("ssd2x1tdc",     event.ssd2x1tdc,    "ssd2x1tdc[ssd2x1nhits]/D");
  tree->Branch("ssd2x1de",      event.ssd2x1de,     "ssd2x1de[ssd2x1nhits]/D");
  tree->Branch("ssd2x1time",    event.ssd2x1time,   "ssd2x1time[ssd2x1nhits]/D");
  tree->Branch("ssd2x1chisqr",  event.ssd2x1chisqr, "ssd2x1chisqr[ssd2x1nhits]/D");
  // SSD2-Y1
  tree->Branch("ssd2y1nhits",  &event.ssd2y1nhits,  "ssd2y1nhits/I");
  tree->Branch("ssd2y1hitpat",  event.ssd2y1hitpat, "ssd2y1hitpat[ssd2y1nhits]/D");
  tree->Branch("ssd2y1pos",     event.ssd2y1pos,    "ssd2y1pos[ssd2y1nhits]/D");
  tree->Branch("ssd2y1ped",     event.ssd2y1ped,    "ssd2y1ped[ssd2y1nhits]/D");
  tree->Branch("ssd2y1adc",     event.ssd2y1adc,    "ssd2y1adc[ssd2y1nhits]/D");
  tree->Branch("ssd2y1tdc",     event.ssd2y1tdc,    "ssd2y1tdc[ssd2y1nhits]/D");
  tree->Branch("ssd2y1de",      event.ssd2y1de,     "ssd2y1de[ssd2y1nhits]/D");
  tree->Branch("ssd2y1time",    event.ssd2y1time,   "ssd2y1time[ssd2y1nhits]/D");
  tree->Branch("ssd2y1chisqr",  event.ssd2y1chisqr, "ssd2y1chisqr[ssd2y1nhits]/D");

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
  // SSD1-Y1 Cluster
  tree->Branch("ssd1y1ncl",    &event.ssd1y1ncl,    "ssd1y1ncl/I");
  tree->Branch("ssd1y1clsize",  event.ssd1y1clsize, "ssd1y1clsize[ssd1y1ncl]/I");
  tree->Branch("ssd1y1clpos",   event.ssd1y1clpos,  "ssd1y1clpos[ssd1y1ncl]/D");
  tree->Branch("ssd1y1clde",    event.ssd1y1clde,   "ssd1y1clde[ssd1y1ncl]/D");
  tree->Branch("ssd1y1cltime",  event.ssd1y1cltime, "ssd1y1cltime[ssd1y1ncl]/D");
  tree->Branch("ssd1y1cltd",    event.ssd1y1cltd,   "ssd1y1cltd[ssd1y1ncl]/D");
  // SSD1-X1 Cluster
  tree->Branch("ssd1x1ncl",    &event.ssd1x1ncl,    "ssd1x1ncl/I");
  tree->Branch("ssd1x1clsize",  event.ssd1x1clsize, "ssd1x1clsize[ssd1x1ncl]/I");
  tree->Branch("ssd1x1clpos",   event.ssd1x1clpos,  "ssd1x1clpos[ssd1x1ncl]/D");
  tree->Branch("ssd1x1clde",    event.ssd1x1clde,   "ssd1x1clde[ssd1x1ncl]/D");
  tree->Branch("ssd1x1cltime",  event.ssd1x1cltime, "ssd1x1cltime[ssd1x1ncl]/D");
  tree->Branch("ssd1x1cltd",    event.ssd1x1cltd,   "ssd1x1cltd[ssd1x1ncl]/D");

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
  // SSD2-X1 Cluster
  tree->Branch("ssd2x1ncl",    &event.ssd2x1ncl,    "ssd2x1ncl/I");
  tree->Branch("ssd2x1clsize",  event.ssd2x1clsize, "ssd2x1clsize[ssd2x1ncl]/I");
  tree->Branch("ssd2x1clpos",   event.ssd2x1clpos,  "ssd2x1clpos[ssd2x1ncl]/D");
  tree->Branch("ssd2x1clde",    event.ssd2x1clde,   "ssd2x1clde[ssd2x1ncl]/D");
  tree->Branch("ssd2x1cltime",  event.ssd2x1cltime, "ssd2x1cltime[ssd2x1ncl]/D");
  tree->Branch("ssd2x1cltd",    event.ssd2x1cltd,   "ssd2x1cltd[ssd2x1ncl]/D");
  // SSD2-Y1 Cluster
  tree->Branch("ssd2y1ncl",    &event.ssd2y1ncl,    "ssd2y1ncl/I");
  tree->Branch("ssd2y1clsize",  event.ssd2y1clsize, "ssd2y1clsize[ssd2y1ncl]/I");
  tree->Branch("ssd2y1clpos",   event.ssd2y1clpos,  "ssd2y1clpos[ssd2y1ncl]/D");
  tree->Branch("ssd2y1clde",    event.ssd2y1clde,   "ssd2y1clde[ssd2y1ncl]/D");
  tree->Branch("ssd2y1cltime",  event.ssd2y1cltime, "ssd2y1cltime[ssd2y1ncl]/D");
  tree->Branch("ssd2y1cltd",    event.ssd2y1cltd,   "ssd2y1cltd[ssd2y1ncl]/D");

  // SSDT
  tree->Branch("ssdtnhits",    &event.ssdtnhits,    "ssdtnhits/I");
  tree->Branch("ssdthitpat",    event.ssdthitpat,   "ssdthitpat[ssdtnhits]/D");
  tree->Branch("ssdtdc",        event.ssdtdc,       "ssdtdc[ssdtnhits]/D");
  tree->Branch("ssdtime",       event.ssdtime,      "ssdtime[ssdtnhits]/D");
  tree->Branch("ntBcOut",       event.ntBcOut,      "ntBcOut/I");
  tree->Branch("ntSdcIn",       event.ntSdcIn,      "ntSdcIn/I");
  tree->Branch("ntSdcOut",      event.ntSdcOut,     "ntSdcOut/I");

  HPrint();

  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")        &&
      InitializeParameter<DCDriftParamMan>("DCDRFT") &&
      InitializeParameter<DCTdcCalibMan>("DCTDC")    &&
      InitializeParameter<HodoParamMan>("HDPRM")     &&
      InitializeParameter<HodoPHCMan>("HDPHC")       &&
      InitializeParameter<SsdParamMan>("SSDPRM")     &&
      InitializeParameter<UserParamMan>("USER")      );
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
