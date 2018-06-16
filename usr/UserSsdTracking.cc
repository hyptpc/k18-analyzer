/**
 *  file: UserSsdTracking.cc
 *  date: 2017.04.10
 *
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "KuramaLib.hh"
#include "MathTools.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "VEvent.hh"

#define HodoCut         0
#define MaxMultiCut     0

#define SlopeFilter     0
#define DeltaEFilter    1
#define TimeFilter      1
#define ChisqrFilter    1

namespace
{
  using namespace root;
  const std::string class_name("EventSsdTracking");
  RMAnalyzer&         gRM   = RMAnalyzer::GetInstance();
  const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const double& zK18Target = gGeom.LocalZ("K18Target");
  const double& zTarget    = gGeom.LocalZ("Target");
  const double& zEmulsion  = gGeom.LocalZ("Emulsion");
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
class EventSsdTracking : public VEvent
{
private:
  RawData      *rawData;
  DCAnalyzer   *DCAna;
  HodoAnalyzer *hodoAna;

public:
        EventSsdTracking( void );
       ~EventSsdTracking( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventSsdTracking::EventSsdTracking( void )
  : VEvent(),
    rawData(0),
    DCAna( new DCAnalyzer ),
    hodoAna( new HodoAnalyzer )
{
}

//______________________________________________________________________________
EventSsdTracking::~EventSsdTracking( void )
{
  if( DCAna )   delete DCAna;
  if( hodoAna ) delete hodoAna;
  if( rawData ) delete rawData;
}

//______________________________________________________________________________
bool
EventSsdTracking::ProcessingBegin( void )
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
  double ssd1y0de[NumOfSegSSD1];
  double ssd1y0pos[NumOfSegSSD1];

  int    ssd1x0nhits;
  double ssd1x0hitpat[NumOfSegSSD1];
  double ssd1x0de[NumOfSegSSD1];
  double ssd1x0pos[NumOfSegSSD1];

  int    ssd1y1nhits;
  double ssd1y1hitpat[NumOfSegSSD1];
  double ssd1y1de[NumOfSegSSD1];
  double ssd1y1pos[NumOfSegSSD1];

  int    ssd1x1nhits;
  double ssd1x1hitpat[NumOfSegSSD1];
  double ssd1x1de[NumOfSegSSD1];
  double ssd1x1pos[NumOfSegSSD1];

  int    ssd2x0nhits;
  double ssd2x0hitpat[NumOfSegSSD2];
  double ssd2x0de[NumOfSegSSD2];
  double ssd2x0pos[NumOfSegSSD2];

  int    ssd2y0nhits;
  double ssd2y0hitpat[NumOfSegSSD2];
  double ssd2y0de[NumOfSegSSD2];
  double ssd2y0pos[NumOfSegSSD2];

  int    ssd2x1nhits;
  double ssd2x1hitpat[NumOfSegSSD2];
  double ssd2x1de[NumOfSegSSD2];
  double ssd2x1pos[NumOfSegSSD2];

  int    ssd2y1nhits;
  double ssd2y1hitpat[NumOfSegSSD2];
  double ssd2y1de[NumOfSegSSD2];
  double ssd2y1pos[NumOfSegSSD2];

  // SSDT
  int    ssdtnhits;
  double ssdthitpat[NumOfSegSSDT];
  double ssdtime[NumOfSegSSDT];

  // Tracking
  int    ntBcOut;
  double chisqrBcOut[MaxHits];
  double x0BcOut[MaxHits];
  double y0BcOut[MaxHits];
  double u0BcOut[MaxHits];
  double v0BcOut[MaxHits];
  int    ntSsdIn;
  double chisqrSsdIn[MaxHits];
  double x0SsdIn[MaxHits];
  double y0SsdIn[MaxHits];
  double u0SsdIn[MaxHits];
  double v0SsdIn[MaxHits];
  double deSsdIn[NumOfLayersSsdIn][MaxHits];
  int    ntSsdOut;
  double chisqrSsdOut[MaxHits];
  double x0SsdOut[MaxHits];
  double y0SsdOut[MaxHits];
  double u0SsdOut[MaxHits];
  double v0SsdOut[MaxHits];
  double deSsdOut[NumOfLayersSsdOut][MaxHits];
  int    ntSdcIn;
  double chisqrSdcIn[MaxHits];
  double x0SdcIn[MaxHits];
  double y0SdcIn[MaxHits];
  double u0SdcIn[MaxHits];
  double v0SdcIn[MaxHits];

};

//______________________________________________________________________________
namespace root
{
  Event  event;
  TH1   *h[MaxHist];
  TTree *tree;
}

//______________________________________________________________________________
bool
EventSsdTracking::ProcessingNormal( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"]");

#if HodoCut
  static const double MinDeBH2   = gUser.GetParameter("DeBH2", 0);
  static const double MaxDeBH2   = gUser.GetParameter("DeBH2", 1);
  static const double MinDeBH1   = gUser.GetParameter("DeBH1", 0);
  static const double MaxDeBH1   = gUser.GetParameter("DeBH1", 1);
  static const double MinBeamToF = gUser.GetParameter("BTOF",  1);
  static const double MaxBeamToF = gUser.GetParameter("BTOF",  1);
#endif
#if MaxMultiCut
  static const double MaxMultiHitBcOut  = gUser.GetParameter("MaxMultiHitBcOut");
  static const double MaxMultiHitSsdIn  = gUser.GetParameter("MaxMultiHitSsdIn");
  static const double MaxMultiHitSsdOut = gUser.GetParameter("MaxMultiHitSsdOut");
  static const double MaxMultiHitSdcIn  = gUser.GetParameter("MaxMultiHitSdcIn");
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

  double time0 = -9999.;
  ////////// BH2 Analysis
  for( int i=0; i<nhBh2; ++i ){
    BH2Hit *hit = hodoAna->GetHitBH2(i);
    if(!hit) continue;
    double cmt = hit->CMeanTime();
    double ct0 = hit->CTime0();
    double min_time = -9999;
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

  double btof0 = -9999;
  for(int i=0; i<nhBh1; ++i){
    Hodo2Hit* hit = hodoAna->GetHitBH1(i);
    if(!hit) continue;
    double cmt  = hit->CMeanTime();
    double btof = cmt - time0;
#if HodoCut
    double dE   = hit->DeltaE();
    if( dE<MinDeBH1 || MaxDeBH1<dE )
      continue;
    if( btof<MinBeamToF || MaxBeamToF<btof )
      continue;
#endif
    if( std::abs(btof)<std::abs(btof0) ){
      btof0 = btof;
    }
  }
  HF1( 100, btof0 );

  HF1(1, 5);

  HF1(1, 6);

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
      }
    }
    event.ssdtnhits = ssdtnhits;
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
      t0Ssd[seg-1]         = time;
    }
  }


  //////////////////// Filtering
  DCAna->DecodeRawHits( rawData );
  DCAna->DoTimeCorrectionSsd( t0Ssd );
  HF1( 1, 7. );

#if SlopeFilter
  DCAna->SlopeFilterSsd();
#endif

#if DeltaEFilter
  static const double MinDeSSD = gUser.GetParameter("DeSSD",0);
  static const double MaxDeSSD = gUser.GetParameter("DeSSD",1);
  DCAna->DeltaEFilterSsd( MinDeSSD, MaxDeSSD );
#endif

#if TimeFilter
  static const double MinTimeSSD = gUser.GetParameter("TimeSSD",0);
  static const double MaxTimeSSD = gUser.GetParameter("TimeSSD",1);
  DCAna->TimeFilterSsd( MinTimeSSD, MaxTimeSSD );
#endif

#if ChisqrFilter
  static const double MaxChisqrSSD = gUser.GetParameter("MaxChisqrSSD",0);
  DCAna->ChisqrFilterSsd( MaxChisqrSSD );
#endif

  double multiSsdIn  = 0.;
  double multiSsdOut = 0.;

  // SsdIn
  int ssd1y0nhits = 0;
  int ssd1x0nhits = 0;
  int ssd1y1nhits = 0;
  int ssd1x1nhits = 0;
  for( int layer=1; layer<=NumOfLayersSsdIn; ++layer ){
    const DCHitContainer &contIn = DCAna->GetSsdInHC(layer);
    int nhIn = contIn.size();
    multiSsdIn += double(nhIn);
    for( int i=0; i<nhIn; ++i ){
      DCHit  *hit  = contIn[i];
      if( !hit ) continue;
      bool    flag = hit->IsZeroSuppressed();
      bool    good = hit->IsGoodWaveForm();
      int     seg  = hit->GetWire();
      double  pos  = hit->GetWirePosition();
      double  de   = hit->GetDe();
      if( !flag ) continue;
      if( !good ) continue;

      switch(layer){
      case 1:
	event.ssd1y0de[ssd1y0nhits]     = de;
	event.ssd1y0hitpat[ssd1y0nhits] = seg;
	event.ssd1y0pos[ssd1y0nhits]    = pos;
	ssd1y0nhits++;
	break;
      case 2:
	event.ssd1x0de[ssd1x0nhits]     = de;
	event.ssd1x0hitpat[ssd1x0nhits] = seg;
	event.ssd1x0pos[ssd1x0nhits]    = pos;
	ssd1x0nhits++;
	break;
      case 3:
	event.ssd1y1de[ssd1y1nhits]     = de;
	event.ssd1y1hitpat[ssd1y1nhits] = seg;
	event.ssd1y1pos[ssd1y1nhits]    = pos;
	ssd1y1nhits++;
	break;
      case 4:
	event.ssd1x1de[ssd1x1nhits]     = de;
	event.ssd1x1hitpat[ssd1x1nhits] = seg;
	event.ssd1x1pos[ssd1x1nhits]    = pos;
	ssd1x1nhits++;
	break;
      default:
	break;
      }
    }
  }
  event.ssd1y0nhits = ssd1y0nhits;
  event.ssd1x0nhits = ssd1x0nhits;
  event.ssd1y1nhits = ssd1y1nhits;
  event.ssd1x1nhits = ssd1x1nhits;

  HF1(1, 8);

  // SsdOut
  int ssd2x0nhits = 0;
  int ssd2y0nhits = 0;
  int ssd2x1nhits = 0;
  int ssd2y1nhits = 0;
  for( int layer=1; layer<=NumOfLayersSsdOut; ++layer ){
    const DCHitContainer &contOut = DCAna->GetSsdOutHC(layer);
    int nhOut = contOut.size();
    multiSsdOut += double(nhOut);
    for( int i=0; i<nhOut; ++i ){
      DCHit  *hit  = contOut[i];
      if(!hit) continue;
      bool    flag = hit->IsZeroSuppressed();
      bool    good = hit->IsGoodWaveForm();
      int     seg  = hit->GetWire();
      double  pos  = hit->GetWirePosition();
      double  de   = hit->GetDe();
      if(!flag) continue;
      if(!good) continue;

      switch(layer){
      case 1:
	event.ssd2x0de[ssd2x0nhits]     = de;
	event.ssd2x0hitpat[ssd2x0nhits] = seg;
	event.ssd2x0pos[ssd2x0nhits]    = pos;
	ssd2x0nhits++;
	break;
      case 2:
	event.ssd2y0de[ssd2y0nhits]     = de;
	event.ssd2y0hitpat[ssd2y0nhits] = seg;
	event.ssd2y0pos[ssd2y0nhits]    = pos;
	ssd2y0nhits++;
	break;
      case 3:
	event.ssd2x1de[ssd2x1nhits]     = de;
	event.ssd2x1hitpat[ssd2x1nhits] = seg;
	event.ssd2x1pos[ssd2x1nhits]    = pos;
	ssd2x1nhits++;
	break;
      case 4:
	event.ssd2y1de[ssd2y1nhits]     = de;
	event.ssd2y1hitpat[ssd2y1nhits] = seg;
	event.ssd2y1hitpat[ssd2y1nhits] = pos;
	ssd2y1nhits++;
	break;
      }
    }
  }

  event.ssd2x0nhits = ssd2x0nhits;
  event.ssd2y0nhits = ssd2y0nhits;
  event.ssd2x1nhits = ssd2x1nhits;
  event.ssd2y1nhits = ssd2y1nhits;

  HF1(1, 9);

  int ssd1nlayers = (ssd1y0nhits>0) + (ssd1x0nhits>0)
    + (ssd1y1nhits>0) + (ssd1x1nhits>0);
  // int ssd2nlayers = (ssd2x0nhits>0) + (ssd2y0nhits>0)
  //   + (ssd2x1nhits>0) + (ssd2y1nhits>0);

  HF1( 1, 10. );

  if( ssd1nlayers < 4 ) return true;
  // if( ssd2nlayers < 4 ) return true;

  HF1( 1, 11. );

#if MaxMultiCut
  if( multiSsdIn/double(NumOfLayersSsdIn) > MaxMultiHitSsdIn )
    return true;
#endif

  HF1( 1, 12. );

#if MaxMultiCut
  if( multiSsdOut/double(NumOfLayersSsdOut) > MaxMultiHitSsdOut )
    return true;
#endif

  HF1( 1, 13. );

  {
    double multi_BcOut=0.;
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &contOut = DCAna->GetBcOutHC(layer);
      int nhOut = contOut.size();
      multi_BcOut += double(nhOut);
    }
#if MaxMultiCut
    if( multi_BcOut/double(NumOfLayersBcOut) > MaxMultiHitBcOut )
      return true;
#endif
  }

  HF1( 1, 14. );

  double xtgt_bcout = -9999.;
  double ytgt_bcout = -9999.;
  double utgt_bcout = -9999.;
  double vtgt_bcout = -9999.;
  // std::cout << "==========TrackSearch BcOut============" << std::endl;
  DCAna->TrackSearchBcOut();
  {
    int nt=DCAna->GetNtracksBcOut();
    event.ntBcOut = nt;
    if( nt > MaxHits ){
      std::cout << "#W " << func_name
		<< " too many BcOut tracks [" << nt
		<< "/" << MaxHits << "]" << std::endl;
      nt = MaxHits;
    }
    for( int it=0; it<nt; ++it ){
      DCLocalTrack *tp = DCAna->GetTrackBcOut(it);
      double chisqr=tp->GetChiSquare();
      double x0 = tp->GetX0(), y0 = tp->GetY0();
      double u0 = tp->GetU0(), v0 = tp->GetV0();
      double xtgt = tp->GetX( zK18Target );
      double ytgt = tp->GetY( zK18Target );
      double utgt = tp->GetU0(), vtgt = tp->GetV0();
      event.chisqrBcOut[it] = chisqr;
      event.x0BcOut[it]     = x0;
      event.y0BcOut[it]     = y0;
      event.u0BcOut[it]     = u0;
      event.v0BcOut[it]     = v0;
      if(nt==1){
	xtgt_bcout = xtgt;
	ytgt_bcout = ytgt;
	utgt_bcout = utgt;
	vtgt_bcout = vtgt;
      }
    }
  }
  //if(nt==0) return true;

  HF1( 1, 15. );

  int    ntSsdIn    = 0;
  double xtgt_ssdin = -9999.;
  double ytgt_ssdin = -9999.;
  double utgt_ssdin = -9999.;
  double vtgt_ssdin = -9999.;
  // std::cout << "==========TrackSearch SsdIn============" << std::endl;
  DCAna->TrackSearchSsdIn();
  {
    int nt = DCAna->GetNtracksSsdIn();
    ntSsdIn = nt;
    event.ntSsdIn = nt;
    HF1( 10000 +0, double(nt) );
    //HF1( 10000 +0, ntracks );
    if( ntSsdIn > MaxHits ){
      std::cout << "#W " << func_name
		<< " too many SsdIn tracks [" << ntSsdIn
		<< "/" << MaxHits << "]" << std::endl;
      nt = MaxHits;
    }
    for( int it=0; it<nt; ++it ){
      DCLocalTrack *tp = DCAna->GetTrackSsdIn(it);
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double x0=tp->GetX0(), y0=tp->GetY0();
      double u0=tp->GetU0(), v0=tp->GetV0();
      // double cost  = 1./std::sqrt(1.+u0*u0+v0*v0);
      // double theta = std::acos(cost)*math::Rad2Deg();
      HF1( 10000 +1, double(nh) );
      HF1( 10000 +5, chisqr );
      HF1( 10000 +10, x0 ); HF1( 10000 +11, y0 );
      HF1( 10000 +12, u0 ); HF1( 10000 +13, v0 );
      HF2( 10000 +14, x0, u0 );
      HF2( 10000 +15, y0, v0 );
      HF2( 10000 +16, x0, y0 );

      double xtgt = tp->GetX( zK18Target );
      double ytgt = tp->GetY( zK18Target );
      double utgt = u0, vtgt = v0;
      if(nt==1){
	xtgt_ssdin = xtgt;
	ytgt_ssdin = ytgt;
	utgt_ssdin = utgt;
	vtgt_ssdin = vtgt;
      }
      HF1( 10000 +20, xtgt ); HF1( 10000 +21, ytgt );
      HF1( 10000 +22, utgt ); HF1( 10000 +23, vtgt );
      HF2( 10000 +24, xtgt, utgt );
      HF2( 10000 +25, ytgt, vtgt );
      HF2( 10000 +26, xtgt, ytgt );

      double xemul=tp->GetX( zEmulsion ), yemul=tp->GetY( zEmulsion );
      double uemul=u0,                    vemul=v0;
      HF1( 10000 +30, xemul ); HF1( 10000 +31, yemul );
      HF1( 10000 +32, uemul ); HF1( 10000 +33, vemul );
      HF2( 10000 +34, xemul, uemul );
      HF2( 10000 +35, yemul, vemul );
      HF2( 10000 +36, xemul, yemul );

      event.chisqrSsdIn[it] = chisqr;
      event.x0SsdIn[it]     = x0;
      event.y0SsdIn[it]     = y0;
      event.u0SsdIn[it]     = u0;
      event.v0SsdIn[it]     = v0;

      for( int ih=0; ih<nh; ++ih ){
	DCLTrackHit *hit = tp->GetHit(ih);
	double de  = hit->GetDe();
	event.deSsdIn[ih][it]  = de;

	int layerId = hit->GetLayer();
	if( layerId>140 ) layerId = layerId - PlOffsSsd;
	HF1( 10000 +6, layerId );
	double wire = hit->GetWire();
	HF1( 10000 +layerId*1000 +2, wire-0.5 );
	double pos = hit->GetLocalHitPos();
	HF1( 10000 +layerId*1000 +3, pos );
      }
    }
  }

  // Residual SsdIn - BcOut
  if( xtgt_ssdin!=-9999. && xtgt_bcout!=-9999. &&
      ytgt_ssdin!=-9999. && ytgt_bcout!=-9999. &&
      utgt_ssdin!=-9999. && utgt_bcout!=-9999. &&
      vtgt_ssdin!=-9999. && vtgt_bcout!=-9999. ){
    double xres = xtgt_ssdin - xtgt_bcout;
    double yres = ytgt_ssdin - ytgt_bcout;
    double ures = utgt_ssdin - utgt_bcout;
    double vres = vtgt_ssdin - vtgt_bcout;
    HF1( 1000 +1, xres  );
    HF1( 1000 +2, yres  );
    HF1( 1000 +3, ures  );
    HF1( 1000 +4, vres  );
  }

  HF1( 1, 16. );

  double xtgt_ssdout = -9999.;
  double ytgt_ssdout = -9999.;
  double utgt_ssdout = -9999.;
  double vtgt_ssdout = -9999.;
  //std::cout << "==========TrackSearch SsdOut============" << std::endl;
  DCAna->TrackSearchSsdOut();
  {
    int ntSsdOut = 0;
    int nt = DCAna->GetNtracksSsdOut();
    HF1( 20000 +0, double(nt) );
    ntSsdOut = nt;
    if( ntSsdOut > MaxHits ){
      std::cout << "#W " << func_name
		<< " too many SsdOut tracks [" << ntSsdOut
		<< "/" << MaxHits << "]" << std::endl;
      nt = MaxHits;
    }
    for( int it=0; it<nt; ++it ){
      DCLocalTrack *tp = DCAna->GetTrackSsdOut(it);
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double x0=tp->GetX0(), y0=tp->GetY0();
      double u0=tp->GetU0(), v0=tp->GetV0();
      // double cost  = 1./std::sqrt(1.+u0*u0+v0*v0);
      // double theta = std::acos(cost)*math::Rad2Deg();
      HF1( 20000 +1, double(nh) );
      HF1( 20000 +5, chisqr );
      HF1( 20000 +10, x0 ); HF1( 20000 +11, y0 );
      HF1( 20000 +12, u0 ); HF1( 20000 +13, v0 );
      HF2( 20000 +14, x0, u0 );
      HF2( 20000 +15, y0, v0 );
      HF2( 20000 +16, x0, y0 );

      double xtgt = tp->GetX( zK18Target );
      double ytgt = tp->GetY( zK18Target );
      double utgt = u0, vtgt = v0;
      if(nt==1){
	xtgt_ssdout = xtgt;
	ytgt_ssdout = ytgt;
	utgt_ssdout = utgt;
	vtgt_ssdout = vtgt;
      }
      HF1( 20000 +20, xtgt ); HF1( 20000 +21, ytgt );
      HF1( 20000 +22, utgt ); HF1( 20000 +23, vtgt );
      HF2( 20000 +24, xtgt, utgt );
      HF2( 20000 +25, ytgt, vtgt );
      HF2( 20000 +26, xtgt, ytgt );

      event.chisqrSsdOut[it] = chisqr;
      event.x0SsdOut[it]     = xtgt;
      event.y0SsdOut[it]     = ytgt;
      event.u0SsdOut[it]     = utgt;
      event.v0SsdOut[it]     = vtgt;

      for( int ih=0; ih<nh; ++ih ){
	DCLTrackHit *hit = tp->GetHit(ih);
	double de  = hit->GetDe();
	event.deSsdOut[ih][it]  = de;

	int layerId = hit->GetLayer();
	if(layerId>140) layerId = layerId - PlOffsSsd;
	HF1( 20000 +6, layerId );
	double wire = hit->GetWire();
	HF1( 20000 +layerId*1000 +2, wire-0.5 );
	double pos = hit->GetLocalHitPos();
	HF1( 20000 +layerId*1000 +3, pos );
      }
    }
    event.ntSsdOut = ntSsdOut;
    HF1( 20000 +0, ntSsdOut );
  }

  // Residual SsdOut - BcOut
  if( xtgt_ssdout!=-9999. && xtgt_bcout!=-9999. &&
      ytgt_ssdout!=-9999. && ytgt_bcout!=-9999. &&
      utgt_ssdout!=-9999. && utgt_bcout!=-9999. &&
      vtgt_ssdout!=-9999. && vtgt_bcout!=-9999. ){
    double xres = xtgt_ssdout - xtgt_bcout;
    double yres = ytgt_ssdout - ytgt_bcout;
    double ures = utgt_ssdout - utgt_bcout;
    double vres = vtgt_ssdout - vtgt_bcout;
    HF1( 2000 +1, xres  );
    HF1( 2000 +2, yres  );
    HF1( 2000 +3, ures  );
    HF1( 2000 +4, vres  );
  }

  HF1( 1, 19. );

  {
    double multi_SdcIn=0.;
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCHitContainer &contIn = DCAna->GetSdcInHC(layer);
      int nhIn = contIn.size();
      multi_SdcIn += double(nhIn);
    }
#if MaxMultiCut
    if( multi_SdcIn/double(NumOfLayersSdcIn) > MaxMultiHitSdcIn )
      return true;
#endif
  }

  double xtgt_sdcin = -9999.;
  double ytgt_sdcin = -9999.;
  double utgt_sdcin = -9999.;
  double vtgt_sdcin = -9999.;
  DCAna->TrackSearchSdcIn();
  {
    int ntSdcIn=DCAna->GetNtracksSdcIn();
    event.ntSdcIn = ntSdcIn;
    if( ntSdcIn > MaxHits ){
      std::cout << "#W " << func_name
		<< " too many SdcIn tracks [" << ntSdcIn
		<< "/" << MaxHits << "]" << std::endl;
      ntSdcIn = MaxHits;
    }
    for( int it=0; it<ntSdcIn; ++it ){
      DCLocalTrack *tp = DCAna->GetTrackSdcIn(it);
      double chisqr=tp->GetChiSquare();
      double x0 = tp->GetX0(), y0 = tp->GetY0();
      double u0 = tp->GetU0(), v0 = tp->GetV0();
      double xtgt = tp->GetX( zTarget );
      double ytgt = tp->GetY( zTarget );
      double utgt = tp->GetU0(), vtgt = tp->GetV0();
      event.chisqrSdcIn[it] = chisqr;
      event.x0SdcIn[it]     = x0;
      event.y0SdcIn[it]     = y0;
      event.u0SdcIn[it]     = u0;
      event.v0SdcIn[it]     = v0;
      if( ntSdcIn==1 ){
	xtgt_sdcin = xtgt;
	ytgt_sdcin = ytgt;
	utgt_sdcin = utgt;
	vtgt_sdcin = vtgt;
      }
    }
  }
  //if(nt==0) return true;

  // Residual Sdcin - BcOut
  if( xtgt_sdcin!=-9999. && xtgt_bcout!=-9999. &&
      ytgt_sdcin!=-9999. && ytgt_bcout!=-9999. &&
      utgt_sdcin!=-9999. && utgt_bcout!=-9999. &&
      vtgt_sdcin!=-9999. && vtgt_bcout!=-9999. ){
    double xres = xtgt_sdcin - xtgt_bcout;
    double yres = ytgt_sdcin - ytgt_bcout;
    double ures = utgt_sdcin - utgt_bcout;
    double vres = vtgt_sdcin - vtgt_bcout;
    HF1( 3000 +1, xres  );
    HF1( 3000 +2, yres  );
    HF1( 3000 +3, ures  );
    HF1( 3000 +4, vres  );
  }

  return true;
}

//______________________________________________________________________________
void
EventSsdTracking::InitializeEvent( void )
{
  event.trignhits   = 0;
  event.ssdtnhits   = 0;
  event.ntBcOut     = 0;
  event.ntSsdIn     = 0;
  event.ntSsdOut    = 0;
  event.ntSdcIn     = 0;

  event.ssd1y0nhits = 0;
  event.ssd1x0nhits = 0;
  event.ssd1y1nhits = 0;
  event.ssd1x1nhits = 0;

  event.ssd2x0nhits = 0;
  event.ssd2y0nhits = 0;
  event.ssd2x1nhits = 0;
  event.ssd2y1nhits = 0;

  for(int it=0; it<NumOfSegTrig; it++){
    event.trigpat[it]  = -1;
    event.trigflag[it] = -9999;
  }

  for(int it=0; it<NumOfSegSSDT; it++){
    event.ssdthitpat[it] = -1;
    event.ssdtime[it] = -9999;
  }

  for(int it=0; it<NumOfSegSSD1; it++){
    event.ssd1y0hitpat[it] = -1;
    event.ssd1y0de[it]  = -9999;
    event.ssd1y0pos[it] = -9999;

    event.ssd1x0hitpat[it] = -1;
    event.ssd1x0de[it]  = -9999;
    event.ssd1x0pos[it] = -9999;

    event.ssd1y1hitpat[it] = -1;
    event.ssd1y1de[it]  = -9999;
    event.ssd1y1pos[it] = -9999;

    event.ssd1x1hitpat[it] = -1;
    event.ssd1y1de[it]  = -9999;
    event.ssd1y1pos[it] = -9999;
  }
  for(int it=0; it<NumOfSegSSD2; it++){
    event.ssd2x0hitpat[it] = -1;
    event.ssd2x0de[it]  = -9999;
    event.ssd2x0pos[it] = -9999;

    event.ssd2y0hitpat[it] = -1;
    event.ssd2y0de[it]  = -9999;
    event.ssd2y0pos[it] = -9999;

    event.ssd2x1hitpat[it] = -1;
    event.ssd2x1de[it]  = -9999;
    event.ssd2x1pos[it] = -9999;

    event.ssd2y1hitpat[it] = -1;
    event.ssd2y1de[it]  = -9999;
    event.ssd2y1pos[it] = -9999;
  }

  for( int it=0; it<MaxHits; it++){
    event.chisqrBcOut[it] = -1.0;
    event.x0BcOut[it] = -9999.;
    event.y0BcOut[it] = -9999.;
    event.u0BcOut[it] = -9999.;
    event.v0BcOut[it] = -9999.;

    event.chisqrSsdIn[it] = -1.0;
    event.x0SsdIn[it] = -9999.;
    event.y0SsdIn[it] = -9999.;
    event.u0SsdIn[it] = -9999.;
    event.v0SsdIn[it] = -9999.;

    for( int ih=0; ih<NumOfLayersSsdIn; ih++){
      event.deSsdIn[ih][it]  = -9999.;
    }

    event.chisqrSsdOut[it] = -1.0;
    event.x0SsdOut[it] = -9999.;
    event.y0SsdOut[it] = -9999.;
    event.u0SsdOut[it] = -9999.;
    event.v0SsdOut[it] = -9999.;
    for( int ih=0; ih<NumOfLayersSsdOut; ih++){
      event.deSsdOut[ih][it]  = -9999.;
    }

    event.chisqrSdcIn[it] = -1.0;
    event.x0SdcIn[it] = -9999.;
    event.y0SdcIn[it] = -9999.;
    event.u0SdcIn[it] = -9999.;
    event.v0SdcIn[it] = -9999.;

  }

}

//______________________________________________________________________________
bool
EventSsdTracking::ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventSsdTracking;
}

//______________________________________________________________________________
namespace
{
  const int NbinDe   = 510;
  const double MinDe = -1000.;
  const double MaxDe = 50000.;

  const int NbinX   = 300;
  const double MinX = -30.;
  const double MaxX = -30.;

  const int NbinY   = 240;
  const double MinY = -60.;
  const double MaxY = -60.;
}

//______________________________________________________________________________
bool
ConfMan::InitializeHistograms( void )
{
  HB1( 1, "Status", 20, 0., 20. );

  HB1( 10, "Trigger HitPat", NumOfSegTrig, 0., double(NumOfSegTrig) );
  for(int i=0; i<NumOfSegTrig; ++i){
    HB1( 10 +i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000 );
  }

  HB1( 100, "Beam ToF", 200,  -10., 5.);

  // SsdInTracking
  HB1( 10000 + 0, "#Tracks SsdIn", 10, 0., 10. );
  HB1( 10000 + 1, "#Hits of Track SsdIn", 20, 0., 20. );
  HB1( 10000 + 4, "Residual SsdIn", 500, 0., 50. );
  HB1( 10000 + 5, "Chisqr SsdIn", 500, 0., 50. );
  HB1( 10000 + 6, "LayerId SsdIn", 20, 0., 20. );
  HB1( 10000 +10, "X0 SsdIn", 400, -100., 100. );
  HB1( 10000 +11, "Y0 SsdIn", 400, -100., 100. );
  HB1( 10000 +12, "U0 SsdIn", 200, -0.20, 0.20 );
  HB1( 10000 +13, "V0 SsdIn", 200, -0.20, 0.20 );
  HB2( 10000 +14, "U0%X0 SsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 10000 +15, "V0%Y0 SsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 10000 +16, "Y0%X0 SsdIn", 100, -100., 100., 100, -100., 100. );
  HB1( 10000 +20, "Xtgt SsdIn", 400, -100., 100. );
  HB1( 10000 +21, "Ytgt SsdIn", 400, -100., 100. );
  HB1( 10000 +22, "Utgt SsdIn", 200, -0.20, 0.20 );
  HB1( 10000 +23, "Vtgt SsdIn", 200, -0.20, 0.20 );
  HB2( 10000 +24, "Utgt%Xtgt SsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 10000 +25, "Vtgt%Ytgt SsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 10000 +26, "Ytgt%Xtgt SsdIn", 100, -100., 100., 100, -100., 100. );
  HB1( 10000 +30, "Xemul SsdIn", 400, -100., 100. );
  HB1( 10000 +31, "Yemul SsdIn", 400, -100., 100. );
  HB1( 10000 +32, "Uemul SsdIn", 200, -0.20, 0.20 );
  HB1( 10000 +33, "Vemul SsdIn", 200, -0.20, 0.20 );
  HB2( 10000 +34, "Uemul%Xemul SsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 10000 +35, "Vemul%Yemul SsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 10000 +36, "Yemul%Xemul SsdIn", 100, -100., 100., 100, -100., 100. );

  for( int i=1; i<=NumOfLayersSsdIn; ++i ){
    HB1( 10000 +i*1000 +2, Form("HitPat Layer%02d SsdIn [Track]", i), 100, 0., 100. );
    HB1( 10000 +i*1000 +3, Form("Position SsdIn %d", i), 100, -250., 250. );
  }

  // SsdOutTracking
  HB1( 20000 + 0, "#Tracks SsdOut", 10, 0., 10. );
  HB1( 20000 + 1, "#Hits of Track SsdOut", 20, 0., 20. );
  HB1( 20000 + 4, "Residual SsdOut", 500, 0., 50. );
  HB1( 20000 + 5, "Chisqr SsdOut", 500, 0., 50. );
  HB1( 20000 + 6, "LayerId SsdOut", 20, 0., 20. );
  HB1( 20000 +10, "X0 SsdOut", 400, -100., 100. );
  HB1( 20000 +11, "Y0 SsdOut", 400, -100., 100. );
  HB1( 20000 +12, "U0 SsdOut", 200, -0.20, 0.20 );
  HB1( 20000 +13, "V0 SsdOut", 200, -0.20, 0.20 );
  HB2( 20000 +14, "U0%X0 SsdOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 20000 +15, "V0%Y0 SsdOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 20000 +16, "Y0%X0 SsdOut", 100, -100., 100., 100, -100., 100. );
  HB1( 20000 +20, "Xtgt SsdOut", 400, -100., 100. );
  HB1( 20000 +21, "Ytgt SsdOut", 400, -100., 100. );
  HB1( 20000 +22, "Utgt SsdOut", 200, -0.20, 0.20 );
  HB1( 20000 +23, "Vtgt SsdOut", 200, -0.20, 0.20 );
  HB2( 20000 +24, "Utgt%Xtgt SsdOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 20000 +25, "Vtgt%Ytgt SsdOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 20000 +26, "Ytgt%Xtgt SsdOut", 100, -100., 100., 100, -100., 100. );
  HB1( 20000 +30, "Xemul SsdOut", 400, -100., 100. );
  HB1( 20000 +31, "Yemul SsdOut", 400, -100., 100. );
  HB1( 20000 +32, "Uemul SsdOut", 200, -0.20, 0.20 );
  HB1( 20000 +33, "Vemul SsdOut", 200, -0.20, 0.20 );
  HB2( 20000 +34, "Uemul%Xemul SsdOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 20000 +35, "Vemul%Yemul SsdOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 20000 +36, "Yemul%Xemul SsdOut", 100, -100., 100., 100, -100., 100. );

  for( int i=1; i<=NumOfLayersSsdOut; ++i ){
    HB1( 20000 +i*1000 +2, Form("HitPat Layer%02d SsdOut [Track]", i), 100, 0., 100. );
    HB1( 20000 +i*1000 +3, Form("Position SsdOut %d", i), 100, -250., 250. );
  }

  // Residual
  HB1( 1000 +1, "Residual SsdIn-BcOut X [K18Target]",  500,  -50.0,  50.0 );
  HB1( 1000 +2, "Residual SsdIn-BcOut Y [K18Target]",  200,  -20.0,  20.0 );
  HB1( 1000 +3, "Residual SsdIn-BcOut U [K18Target]",  500,   -2.0,   2.0 );
  HB1( 1000 +4, "Residual SsdIn-BcOut V [K18Target]",  500,   -2.0,   2.0 );
  HB1( 2000 +1, "Residual SsdOut-BcOut X [K18Target]", 500,  -50.0,  50.0 );
  HB1( 2000 +2, "Residual SsdOut-BcOut Y [K18Target]", 200,  -20.0,  20.0 );
  HB1( 2000 +3, "Residual SsdOut-BcOut U [K18Target]", 500,   -2.0,   2.0 );
  HB1( 2000 +4, "Residual SsdOut-BcOut V [K18Target]", 500,   -2.0,   2.0 );
  HB1( 3000 +1, "Residual SdcIn-BcOut X [K18Target]",  5000,  -500.0,  500.0 );
  HB1( 3000 +2, "Residual SdcIn-BcOut Y [K18Target]",  2000,  -200.0,  200.0 );
  HB1( 3000 +3, "Residual SdcIn-BcOut U [K18Target]",  500,   -2.0,   2.0 );
  HB1( 3000 +4, "Residual SdcIn-BcOut V [K18Target]",  500,   -2.0,   2.0 );

  ////////////////////////////////////////////
  //Tree
  HBTree( "ssd","tree of SsdTracking" );
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trignhits", &event.trignhits, "trignhits/I");
  tree->Branch("trigpat",    event.trigpat,   "trigpat[trignhits]/I");
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  // SsdIn
  tree->Branch("ssd1y0nhits",  &event.ssd1y0nhits,  "ssd1y0nhits/I");
  tree->Branch("ssd1y0hitpat",  event.ssd1y0hitpat, "ssd1y0hitpat[ssd1y0nhits]/D");
  tree->Branch("ssd1y0de",      event.ssd1y0de,     "ssd1y0de[ssd1y0nhits]/D");
  tree->Branch("ssd1y0pos",     event.ssd1y0pos,    "ssd1y0pos[ssd1y0nhits]/D");
  tree->Branch("ssd1x0nhits",  &event.ssd1x0nhits,  "ssd1x0nhits/I");
  tree->Branch("ssd1x0hitpat",  event.ssd1x0hitpat, "ssd1x0hitpat[ssd1x0nhits]/D");
  tree->Branch("ssd1x0de",      event.ssd1x0de,     "ssd1x0de[ssd1x0nhits]/D");
  tree->Branch("ssd1x0pos",     event.ssd1x0pos,    "ssd1x0pos[ssd1x0nhits]/D");
  tree->Branch("ssd1y1nhits",  &event.ssd1y1nhits,  "ssd1y1nhits/I");
  tree->Branch("ssd1y1hitpat",  event.ssd1y1hitpat, "ssd1y1hitpat[ssd1y1nhits]/D");
  tree->Branch("ssd1y1de",      event.ssd1y1de,     "ssd1y1de[ssd1y1nhits]/D");
  tree->Branch("ssd1y1pos",     event.ssd1y1pos,    "ssd1y1pos[ssd1y1nhits]/D");
  tree->Branch("ssd1x1nhits",  &event.ssd1x1nhits,  "ssd1x1nhits/I");
  tree->Branch("ssd1x1hitpat",  event.ssd1x1hitpat, "ssd1x1hitpat[ssd1x1nhits]/D");
  tree->Branch("ssd1x1de",      event.ssd1x1de,     "ssd1x1de[ssd1x1nhits]/D");
  tree->Branch("ssd1x1pos",     event.ssd1x1pos,    "ssd1x1pos[ssd1x1nhits]/D");

  // SsdOut
  tree->Branch("ssd2x0nhits",  &event.ssd2x0nhits,  "ssd2x0nhits/I");
  tree->Branch("ssd2x0hitpat",  event.ssd2x0hitpat, "ssd2x0hitpat[ssd2x0nhits]/D");
  tree->Branch("ssd2x0de",      event.ssd2x0de,     "ssd2x0de[ssd2x0nhits]/D");
  tree->Branch("ssd2x0pos",     event.ssd2x0pos,    "ssd2x0pos[ssd2x0nhits]/D");
  tree->Branch("ssd2y0nhits",  &event.ssd2y0nhits,  "ssd2y0nhits/I");
  tree->Branch("ssd2y0hitpat",  event.ssd2y0hitpat, "ssd2y0hitpat[ssd2y0nhits]/D");
  tree->Branch("ssd2y0de",      event.ssd2y0de,     "ssd2y0de[ssd2y0nhits]/D");
  tree->Branch("ssd2y0pos",     event.ssd2y0pos,    "ssd2y0pos[ssd2y0nhits]/D");
  tree->Branch("ssd2x1nhits",  &event.ssd2x1nhits,  "ssd2x1nhits/I");
  tree->Branch("ssd2x1hitpat",  event.ssd2x1hitpat, "ssd2x1hitpat[ssd2x1nhits]/D");
  tree->Branch("ssd2x1de",      event.ssd2x1de,     "ssd2x1de[ssd2x1nhits]/D");
  tree->Branch("ssd2x1pos",     event.ssd2x1pos,    "ssd2x1pos[ssd2x1nhits]/D");
  tree->Branch("ssd2y1nhits",  &event.ssd2y1nhits,  "ssd2y1nhits/I");
  tree->Branch("ssd2y1hitpat",  event.ssd2y1hitpat, "ssd2y1hitpat[ssd2y1nhits]/D");
  tree->Branch("ssd2y1de",      event.ssd2y1de,     "ssd2y1de[ssd2y1nhits]/D");
  tree->Branch("ssd2y1pos",     event.ssd2y1pos,    "ssd2y1pos[ssd2y1nhits]/D");

  // SSDT
  tree->Branch("ssdtnhits",    &event.ssdtnhits,    "ssdtnhits/I");
  tree->Branch("ssdthitpat",    event.ssdthitpat,   "ssdthitpat[ssdtnhits]/D");
  tree->Branch("ssdtime",       event.ssdtime,      "ssdtime[ssdtnhits]/D");

  tree->Branch("ntBcOut",     &event.ntBcOut,     "ntBcOut/I");
  tree->Branch("chisqrBcOut",  event.chisqrBcOut, "chisqrBcOut[ntBcOut]/D");
  tree->Branch("x0BcOut",      event.x0BcOut,     "x0BcOut[ntBcOut]/D");
  tree->Branch("y0BcOut",      event.y0BcOut,     "y0BcOut[ntBcOut]/D");
  tree->Branch("u0BcOut",      event.u0BcOut,     "u0BcOut[ntBcOut]/D");
  tree->Branch("v0BcOut",      event.v0BcOut,     "v0BcOut[ntBcOut]/D");

  tree->Branch("ntSsdIn",     &event.ntSsdIn,     "ntSsdIn/I");
  tree->Branch("chisqrSsdIn",  event.chisqrSsdIn, "chisqrSsdIn[ntSsdIn]/D");
  tree->Branch("x0SsdIn",      event.x0SsdIn,     "x0SsdIn[ntSsdIn]/D");
  tree->Branch("y0SsdIn",      event.y0SsdIn,     "y0SsdIn[ntSsdIn]/D");
  tree->Branch("u0SsdIn",      event.u0SsdIn,     "u0SsdIn[ntSsdIn]/D");
  tree->Branch("v0SsdIn",      event.v0SsdIn,     "v0SsdIn[ntSsdIn]/D");
  tree->Branch("deSsdIn",      event.deSsdIn,     Form("deSsdIn[%d][%d]/D",
						       NumOfLayersSsdIn, MaxHits));
  tree->Branch("ntSsdOut",     &event.ntSsdOut,     "ntSsdOut/I");
  tree->Branch("chisqrSsdOut",  event.chisqrSsdOut, "chisqrSsdOut[ntSsdOut]/D");
  tree->Branch("x0SsdOut",      event.x0SsdOut,   "x0SsdOut[ntSsdOut]/D");
  tree->Branch("y0SsdOut",      event.y0SsdOut,   "y0SsdOut[ntSsdOut]/D");
  tree->Branch("u0SsdOut",      event.u0SsdOut,   "u0SsdOut[ntSsdOut]/D");
  tree->Branch("v0SsdOut",      event.v0SsdOut,   "v0SsdOut[ntSsdOut]/D");
  tree->Branch("deSsdOut",      event.deSsdOut,     Form("deSsdOut[%d][%d]/D",
							 NumOfLayersSsdOut, MaxHits));

  tree->Branch("ntSdcIn",     &event.ntSdcIn,     "ntSdcIn/I");
  tree->Branch("chisqrSdcIn",  event.chisqrSdcIn, "chisqrSdcIn[ntSdcIn]/D");
  tree->Branch("x0SdcIn",      event.x0SdcIn,     "x0SdcIn[ntSdcIn]/D");
  tree->Branch("y0SdcIn",      event.y0SdcIn,     "y0SdcIn[ntSdcIn]/D");
  tree->Branch("u0SdcIn",      event.u0SdcIn,     "u0SdcIn[ntSdcIn]/D");
  tree->Branch("v0SdcIn",      event.v0SdcIn,     "v0SdcIn[ntSdcIn]/D");

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
