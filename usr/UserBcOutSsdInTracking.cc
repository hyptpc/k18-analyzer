/**
 *  file: UserBcOutSsdInTracking.cc
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
#include "RootHelper.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
#include "UnpackerManager.hh"
#include "UserParamMan.hh"
#include "VEvent.hh"

#define HodoCut     0
#define BcOut       1
#define SsdIn       1
#define BcOutSsdIn  1

namespace
{
  using namespace root;
  const std::string& class_name("EventBcOutSsdInTracking");
  RMAnalyzer&         gRM   = RMAnalyzer::GetInstance();
  const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const double& zK18Target = gGeom.LocalZ("K18Target");
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
class EventBcOutSsdInTracking : public VEvent
{
private:
  RawData      *rawData;
  DCAnalyzer   *DCAna;
  HodoAnalyzer *hodoAna;

public:
        EventBcOutSsdInTracking( void );
       ~EventBcOutSsdInTracking( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventBcOutSsdInTracking::EventBcOutSsdInTracking( void )
  : VEvent(),
    rawData(0),
    DCAna( new DCAnalyzer ),
    hodoAna( new HodoAnalyzer )
{
}

//______________________________________________________________________________
EventBcOutSsdInTracking::~EventBcOutSsdInTracking( void )
{
  if ( DCAna )   delete DCAna;
  if ( hodoAna ) delete hodoAna;
  if ( rawData ) delete rawData;
}

//______________________________________________________________________________
bool
EventBcOutSsdInTracking::ProcessingBegin( void )
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
  double ssd1y0de[NumOfSegSSD1];
  double ssd1y0pos[NumOfSegSSD1];

  int    ssd1x0nhits;
  double ssd1x0de[NumOfSegSSD1];
  double ssd1x0pos[NumOfSegSSD1];

  int    ssd1y1nhits;
  double ssd1y1de[NumOfSegSSD1];
  double ssd1y1pos[NumOfSegSSD1];

  int    ssd1x1nhits;
  double ssd1x1de[NumOfSegSSD1];
  double ssd1x1pos[NumOfSegSSD1];

  // SSDT
  int    ssdtnhits;
  double ssdthitpat[NumOfSegSSDT];
  double ssdt[NumOfSegSSDT];

  // BcOutSsdIn Tracking
  int ntrack;
  double chisqr[MaxHits];
  double x0[MaxHits];
  double y0[MaxHits];
  double u0[MaxHits];
  double v0[MaxHits];

};

//______________________________________________________________________________
namespace root
{
  Event event;
  TH1   *h[MaxHist];
  TTree *tree;
  const int SSD1Hid = 20000;
  const int SSDTHid = 40000;
  const int SsdInHid      = 1000;
  const int BcOutSsdInHid = 2000;
}

//______________________________________________________________________________
bool
EventBcOutSsdInTracking::ProcessingNormal( void )
{
  static const double MaxMultiHitBcOut = gUser.GetParameter("MaxMultiHitBcOut");

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
	HF1( 10, seg-1 );
	HF1( 10+seg, tdc );
      }
    }
    event.trignhits = trignhits;
  }

  // if( event.trigflag[SpillEndFlag] ) return true;

  HF1( 1, 1. );

#if HodoCut
  ////////// BH2 time 0
  hodoAna->DecodeBH2Hits(rawData);
  int nhBh2 = hodoAna->GetNHitsBH2();
  if( nhBh2==0 ) return true;

  HF1( 1, 2 );

  double time0 = -999.;
  ////////// BH2 Analysis
  for( int i=0; i<nhBh2; ++i ){
    BH2Hit *hit = hodoAna->GetHitBH2(i);
    if(!hit) continue;
    double cmt = hit->CMeanTime();
    double ct0 = hit->CTime0();
    double dE  = hit->DeltaE();
    double min_time = -999;
    //------------------------Cut
    // if( dE<MinDeltaEBH2 || MaxDeltaEBH2<dE ) continue;
    if( std::abs(cmt)<std::abs(min_time) ){
      min_time = cmt;
      time0    = ct0;
    }
  }

  HF1(1, 3);

  ////////// BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData);
  int nhBh1 = hodoAna->GetNHitsBH1();
  if(nhBh1==0) return true;

  HF1(1, 4);

  double btof0 = -999;
  for(int i=0; i<nhBh1; ++i){
    Hodo2Hit* hit = hodoAna->GetHitBH1(i);
    if(!hit) continue;
    double cmt  = hit->CMeanTime();
    double dE   = hit->DeltaE();
    double btof = cmt - time0;
    //------------------------Cut
    //  if( dE<MinDeltaEBH1 || MaxDeltaEBH1<dE ) continue;
    //------------------------Cut
    //  if( MinBeamToF<btof && btof<MaxBeamToF ){
    if( std::abs(btof)<std::abs(btof0) ){
      btof0 = btof;
    }
  }
  HF1( 20, btof0 );

  HF1(1, 5);
#endif
  HF1(1, 6);

  //SSDT
  {
    std::string name = "SSDT";
    const HodoRHitContainer &cont = rawData->GetSSDTRawHC();
    int ssdtnhits = 0;
    int nh = cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      int seg = hit->SegmentId()+1;
      int tdc = hit->GetTdc1();
      if(tdc>0){
	event.ssdthitpat[ssdtnhits++] = seg;
	event.ssdt[seg-1]             = tdc;
	HF1( SSDTHid +1, seg-1 );
	HF2( SSDTHid +2, seg-1, tdc );
      }
    }
    event.ssdtnhits = ssdtnhits;
    HF1( SSDTHid +0, ssdtnhits );
  }

  DCAna->DecodeRawHits( rawData );
  HF1(1, 7);

  double multiSsdIn  = 0.;

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
      double  pos  = hit->GetWirePosition();
      double  de   = hit->GetDe();
      bool    good = hit->IsGoodWaveForm();
      if(!good) continue;
      switch(layer){
      case 3:
	event.ssd1y0de[ssd1y0nhits]  = de;
	event.ssd1y0pos[ssd1y0nhits] = pos;
	ssd1y0nhits++;
	HF1( SSD1Hid + 1, pos );
	HF1( SSD1Hid + 2, de  );
	HF2( SSD1Hid + 3, pos, de );
	break;
      case 4:
	event.ssd1x0de[ssd1x0nhits]  = de;
	event.ssd1x0pos[ssd1x0nhits] = pos;
	ssd1x0nhits++;
	HF1( SSD1Hid +11, pos );
	HF1( SSD1Hid +12, de  );
	HF2( SSD1Hid +13, pos, de );
	break;
      case 5:
	event.ssd1y1de[ssd1y1nhits]  = de;
	event.ssd1y1pos[ssd1y1nhits] = pos;
	ssd1y1nhits++;
	HF1( SSD1Hid +21, pos );
	HF1( SSD1Hid +22, de  );
	HF2( SSD1Hid +23, pos, de );
	break;
      case 6:
	event.ssd1x1de[ssd1x1nhits]  = de;
	event.ssd1x1pos[ssd1x1nhits] = pos;
	ssd1x1nhits++;
	HF1( SSD1Hid +31, pos );
	HF1( SSD1Hid +32, de  );
	HF2( SSD1Hid +33, pos, de  );
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
  HF1( SSD1Hid + 0, ssd1y0nhits );
  HF1( SSD1Hid +10, ssd1x0nhits );
  HF1( SSD1Hid +20, ssd1y1nhits );
  HF1( SSD1Hid +30, ssd1x1nhits );

  HF1(1, 8);

  HF1(1, 9);

  // int ssd1nhits = (ssd1y0nhits>0) + (ssd1x0nhits>0)
  //   + (ssd1y1nhits>0) + (ssd1x1nhits>0);

  HF1(1, 10);

  // if( ssd1nhits < 4 ) return true;

  HF1( 1, 11. );

  // if( multiSsdIn/double(NumOfLayersSsdIn-2) > MaxMultiHitSsdIn ){
  //   return true;
  // }

  HF1( 1, 12. );

#if BcOut
  double multi_BcOut=0.;
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &contOut = DCAna->GetBcOutHC(layer);
      int nhOut = contOut.size();
      multi_BcOut += double(nhOut);
    }
  }
  if( multi_BcOut/double(NumOfLayersBcOut) > MaxMultiHitBcOut ){
    return true;
  }

  HF1( 1, 13. );

  double xtgt_bcout = -999.;
  double ytgt_bcout = -999.;
  double utgt_bcout = -999.;
  double vtgt_bcout = -999.;
  // std::cout << "==========TrackSearch BcOut============" << std::endl;
  DCAna->TrackSearchBcOut();
  int nt=DCAna->GetNtracksBcOut();
  HF1( 10, double(nt) );
  for( int it=0; it<nt; ++it ){
    DCLocalTrack *tp = DCAna->GetTrackBcOut(it);
    double xtgt = tp->GetX(zK18Target);
    double ytgt = tp->GetY(zK18Target);
    double utgt = tp->GetU0(), vtgt = tp->GetV0();
    if(nt==1){
      xtgt_bcout = xtgt;
      ytgt_bcout = ytgt;
      utgt_bcout = utgt;
      vtgt_bcout = vtgt;
    }
  }

  //if(nt==0) return true;

#endif

  HF1( 1, 14. );

#if SsdIn
  double xtgt_ssdin = -999.;
  double ytgt_ssdin = -999.;
  double utgt_ssdin = -999.;
  double vtgt_ssdin = -999.;
  // std::cout << "==========TrackSearch SsdIn============" << std::endl;
  DCAna->TrackSearchSsdIn();
  {
    int nt = DCAna->GetNtracksSsdIn();
    HF1( SsdInHid, double(nt) );
    for( int it=0; it<nt; ++it ){
      DCLocalTrack *tp = DCAna->GetTrackSsdIn(it);
      int nh = tp->GetNHit();
      double chisqr = tp->GetChiSquare();
      double x0 = tp->GetX0(), y0 = tp->GetY0();
      double u0 = tp->GetU0(), v0 = tp->GetV0();
      HF1( SsdInHid +  1, double(nh) );
      HF1( SsdInHid +  2, chisqr );
      HF1( SsdInHid +  4, x0 ); HF1( SsdInHid +  5, y0 );
      HF1( SsdInHid +  6, u0 ); HF1( SsdInHid +  7, v0 );
      HF2( SsdInHid +  8, x0, u0 );
      HF2( SsdInHid +  9, y0, v0 );
      HF2( SsdInHid + 10, x0, y0 );
      double xtgt = tp->GetX( zK18Target ), ytgt = tp->GetY( zK18Target );
      double utgt = u0, vtgt = v0;
      if(nt==1){
	xtgt_ssdin = xtgt;
	ytgt_ssdin = ytgt;
	utgt_ssdin = utgt;
	vtgt_ssdin = vtgt;
      }
      HF1( SsdInHid + 11, xtgt ); HF1( SsdInHid + 12, ytgt );
      HF1( SsdInHid + 13, utgt ); HF1( SsdInHid + 14, vtgt );
      HF2( SsdInHid + 15, xtgt, utgt );
      HF2( SsdInHid + 16, ytgt, vtgt );
      HF2( SsdInHid + 17, xtgt, ytgt );
      double xemul = tp->GetX( zEmulsion ), yemul = tp->GetY( zEmulsion );
      double uemul = u0, vemul = v0;
      HF1( SsdInHid + 21, xemul ); HF1( SsdInHid + 22, yemul );
      HF1( SsdInHid + 23, uemul ); HF1( SsdInHid + 24, vemul );
      HF2( SsdInHid + 25, xemul, uemul );
      HF2( SsdInHid + 26, yemul, vemul );
      HF2( SsdInHid + 27, xemul, yemul );
      for( int ih=0; ih<nh; ++ih ){
	DCLTrackHit *hit = tp->GetHit(ih);
	int layerId = hit->GetLayer();
	if(layerId>140) layerId = layerId - PlOffsSsd;
	HF1( SsdInHid + 3, layerId );
	double wire=hit->GetWire();
	HF1( SsdInHid + 100*layerId+11, wire-0.5 );
	double pos=hit->GetLocalHitPos();
	HF1( SsdInHid + 100*layerId+14, pos );
      }
    }
  }
# if BcOut
  if( xtgt_ssdin!=-999. && xtgt_bcout!=-999. &&
      ytgt_ssdin!=-999. && ytgt_bcout!=-999. &&
      utgt_ssdin!=-999. && utgt_bcout!=-999. &&
      vtgt_ssdin!=-999. && vtgt_bcout!=-999. ){
    double xres = xtgt_ssdin - xtgt_bcout;
    double yres = ytgt_ssdin - ytgt_bcout;
    double ures = utgt_ssdin - utgt_bcout;
    double vres = vtgt_ssdin - vtgt_bcout;
    HF1( 100, xres  );
    HF1( 200, yres  );
    HF1( 300, ures  );
    HF1( 400, vres  );
  }
# endif

#endif

  HF1( 1, 15. );

#if BcOutSsdIn
  double xtgt_bcout_ssdin = -999.;
  double ytgt_bcout_ssdin = -999.;
  double utgt_bcout_ssdin = -999.;
  double vtgt_bcout_ssdin = -999.;
  // std::cout << "==========TrackSearch BcOutSsdIn============" << std::endl;
  DCAna->TrackSearchBcOutSsdIn();
  {
    int nhits = 0;
    int nt = DCAna->GetNtracksBcOutSsdIn();
    HF1( BcOutSsdInHid, double(nt) );
    for( int it=0; it<nt; ++it ){
      DCLocalTrack *tp = DCAna->GetTrackBcOutSsdIn(it);
      int nh = tp->GetNHit();
      double chisqr = tp->GetChiSquare();
      double x0 = tp->GetX0(), y0 = tp->GetY0();
      double u0 = tp->GetU0(), v0 = tp->GetV0();

      double cost = 1./sqrt(1.+u0*u0+v0*v0);
      double theta = acos(cost)*math::Rad2Deg();

      HF1( BcOutSsdInHid + 1, double(nh) );
      HF1( BcOutSsdInHid + 2, chisqr );
      HF1( BcOutSsdInHid + 4, x0 ); HF1( BcOutSsdInHid + 5, y0 );
      HF1( BcOutSsdInHid + 6, u0 ); HF1( BcOutSsdInHid + 7, v0 );
      HF2( BcOutSsdInHid + 8, x0, u0 );
      HF2( BcOutSsdInHid + 9, y0, v0 );
      HF2( BcOutSsdInHid + 0, x0, y0 );

      double xtgt = tp->GetX( zK18Target ), ytgt = tp->GetY( zK18Target );
      double utgt = u0, vtgt = v0;
      if(nt==1){
	xtgt_bcout_ssdin = xtgt;
	ytgt_bcout_ssdin = ytgt;
	utgt_bcout_ssdin = utgt;
	vtgt_bcout_ssdin = vtgt;
      }
      HF1( BcOutSsdInHid + 11, xtgt ); HF1( BcOutSsdInHid + 12, ytgt );
      HF1( BcOutSsdInHid + 13, utgt ); HF1( BcOutSsdInHid + 14, vtgt );
      HF2( BcOutSsdInHid + 15, xtgt, utgt );
      HF2( BcOutSsdInHid + 16, ytgt, vtgt );
      HF2( BcOutSsdInHid + 17, xtgt, ytgt );

      double xemul = tp->GetX( zEmulsion ), yemul = tp->GetY( zEmulsion );
      double uemul = u0, vemul = v0;
      HF1( BcOutSsdInHid + 21, xemul ); HF1( BcOutSsdInHid + 22, yemul );
      HF1( BcOutSsdInHid + 23, uemul ); HF1( BcOutSsdInHid + 24, vemul );
      HF2( BcOutSsdInHid + 25, xemul, uemul );
      HF2( BcOutSsdInHid + 26, yemul, vemul );
      HF2( BcOutSsdInHid + 27, xemul, yemul );

      event.chisqr[nhits] = chisqr;
      event.x0[nhits]     = x0;
      event.y0[nhits]     = y0;
      event.u0[nhits]     = u0;
      event.v0[nhits]     = v0;
      nhits++;

      for( int ih=0; ih<nh; ++ih ){
	DCLTrackHit *hit = tp->GetHit(ih);
	int layerId = hit->GetLayer();
	if( layerId>112 )
	  layerId = layerId - PlOffsBc  - 12;
	else
	  layerId = layerId - PlOffsSsd + 12;
	HF1( BcOutSsdInHid + 3, layerId );
	double wire=hit->GetWire();
	HF1( BcOutSsdInHid + 100*layerId+11, wire-0.5 );
	double xcal = hit->GetXcal(), ycal = hit->GetYcal();
	double pos = hit->GetLocalHitPos(), res = hit->GetResidual();
	HF1( BcOutSsdInHid + 100*layerId+14, pos );
	HF1( BcOutSsdInHid + 100*layerId+15, res );
	HF2( BcOutSsdInHid + 100*layerId+16, pos, res );
	HF2( BcOutSsdInHid + 100*layerId+17, xcal, ycal);

	if (theta>=0 && theta<15)
	  HF1( BcOutSsdInHid + 100*layerId+71, res );
	else if (theta>=15 && theta<30)
	  HF1( BcOutSsdInHid + 100*layerId+72, res );
	else if (theta>=30 && theta<45)
	  HF1( BcOutSsdInHid + 100*layerId+73, res );
	else if (theta>=45)
	  HF1( BcOutSsdInHid + 100*layerId+74, res );
      }
    }
    event.ntrack = nhits;

# if SsdIn
    if( xtgt_ssdin!=-999. && xtgt_bcout_ssdin!=-999. &&
	ytgt_ssdin!=-999. && ytgt_bcout_ssdin!=-999. &&
	utgt_ssdin!=-999. && utgt_bcout_ssdin!=-999. &&
	vtgt_ssdin!=-999. && vtgt_bcout_ssdin!=-999. ){
      double xres = xtgt_ssdin - xtgt_bcout_ssdin;
      double yres = ytgt_ssdin - ytgt_bcout_ssdin;
      double ures = utgt_ssdin - utgt_bcout_ssdin;
      double vres = vtgt_ssdin - vtgt_bcout_ssdin;
      HF1( 100, xres  );
      HF1( 200, yres  );
      HF1( 300, ures  );
      HF1( 400, vres  );
    }
# endif
  }
#endif

  HF1( 1, 16. );

  HF1( 1, 19. );

  return true;
}

//______________________________________________________________________________
void
EventBcOutSsdInTracking::InitializeEvent( void )
{
  event.evnum       = 0;
  event.trignhits   = 0;
  event.ssdtnhits   = 0;
  event.ntrack      = 0;
  event.ssd1y0nhits = 0;
  event.ssd1x0nhits = 0;
  event.ssd1y1nhits = 0;
  event.ssd1x1nhits = 0;

  for(int it=0; it<NumOfSegTrig; it++){
    event.trigpat[it]  = -1;
    event.trigflag[it] = -999;
  }

  for(int it=0; it<NumOfSegSSDT; it++){
    event.ssdthitpat[it] = -1;
    event.ssdt[it] = -999;
  }

  for(int it=0; it<NumOfSegSSD1; it++){
    event.ssd1y0de[it]  = -999;
    event.ssd1y0pos[it] = -999;

    event.ssd1x0de[it]  = -999;
    event.ssd1x0pos[it] = -999;

    event.ssd1y1de[it]  = -999;
    event.ssd1y1pos[it] = -999;

    event.ssd1y1de[it]  = -999;
    event.ssd1y1pos[it] = -999;
  }

  for( int it=0; it<MaxHits; it++){
    event.chisqr[it] = -1.0;
    event.x0[it] = -999.;
    event.y0[it] = -999.;
    event.u0[it] = -999.;
    event.v0[it] = -999.;
  }

}

//______________________________________________________________________________
bool
EventBcOutSsdInTracking::ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventBcOutSsdInTracking;
}

//______________________________________________________________________________
namespace
{
  const int NbinAdc   = 2000;
  const double MinAdc =    0.;
  const double MaxAdc = 2000.;

  const int NbinTdc   = 10;
  const double MinTdc =  0.;
  const double MaxTdc = 10.;

  const int NbinDe   = 300;
  const double MinDe = -50.;
  const double MaxDe = 550.;

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
    HB1( 10+i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000 );
  }

  ////////////////////
  // SsdIn
  HB1( SSD1Hid + 0, "#Hits Ssd1 Y0", 100, 0., 100. );
  HB1( SSD1Hid + 1, "Position Ssd1 Y0", NbinY, MinY, MaxY );
  HB1( SSD1Hid + 2, "dE Ssd1 Y0", NbinDe, MinDe, MaxDe );
  HB2( SSD1Hid + 3, "dE Ssd1 Y0 (2D)",
       NbinY, MinY, MaxY, NbinDe, MinDe, MaxDe );

  HB1( SSD1Hid +10, "#Hits Ssd1 X0", 100, 0., 100. );
  HB1( SSD1Hid +11, "Position Ssd1 X0", NbinX, MinX, MaxX );
  HB1( SSD1Hid +12, "dE Ssd1 X0", NbinDe, MinDe, MaxDe );
  HB2( SSD1Hid +13, "dE Ssd1 X0 (2D)",
       NbinX, MinX, MaxX, NbinDe, MinDe, MaxDe );

  HB1( SSD1Hid +20, "#Hits Ssd1 Y1 ", 100, 0., 100. );
  HB1( SSD1Hid +21, "Position Ssd1 Y1", NbinY, MinY, MaxY );
  HB1( SSD1Hid +22, "dE Ssd1 Y1", NbinDe, MinDe, MaxDe );
  HB2( SSD1Hid +23, "dE Ssd1 Y1 (2D)",
       NbinY, MinY, MaxY, NbinDe, MinDe, MaxDe );

  HB1( SSD1Hid +30, "#Hits Ssd1 X1 ", 100, 0., 100. );
  HB1( SSD1Hid +31, "Position Ssd1 X1", NbinX, MinX, MaxX );
  HB1( SSD1Hid +32, "dE Ssd1 X1", NbinDe, MinDe, MaxDe );
  HB2( SSD1Hid +33, "dE Ssd1 X1 (2D)",
       NbinX, MinX, MaxX, NbinDe, MinDe, MaxDe );

  // SSDT
  HB1( SSDTHid +0, "#Hits SSDT",  NumOfSegSSDT, 0., (double)NumOfSegSSDT );
  HB1( SSDTHid +1, "SSDT HitPat", NumOfSegSSDT, 0., (double)NumOfSegSSDT );
  HB2( SSDTHid +2, "SSDT Tdc%Seg",
       NumOfSegSSDT, 0., (double)NumOfSegSSDT, 0x1000, 0, 0x1000 );

  // SsdInTracking
  HB1( SsdInHid +  0, "#Tracks SsdIn", 10, 0., 10. );
  HB1( SsdInHid +  1, "#Hits of Track SsdIn", 20, 0., 20. );
  HB1( SsdInHid +  2, "Chisqr SsdIn", 500, 0., 50. );
  HB1( SsdInHid +  3, "LayerId SsdIn", 20, 0., 20. );
  HB1( SsdInHid +  4, "X0 SsdIn", 400, -100., 100. );
  HB1( SsdInHid +  5, "Y0 SsdIn", 400, -100., 100. );
  HB1( SsdInHid +  6, "U0 SsdIn", 200, -0.20, 0.20 );
  HB1( SsdInHid +  7, "V0 SsdIn", 200, -0.20, 0.20 );
  HB2( SsdInHid +  8, "U0%X0 SsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( SsdInHid +  9, "V0%Y0 SsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( SsdInHid + 10, "X0%Y0 SsdIn", 100, -100., 100., 100, -100, 100 );

  HB1( SsdInHid + 11, "Xtgt SsdIn", 400, -100., 100. );
  HB1( SsdInHid + 12, "Ytgt SsdIn", 400, -100., 100. );
  HB1( SsdInHid + 13, "Utgt SsdIn", 200, -0.20, 0.20 );
  HB1( SsdInHid + 14, "Vtgt SsdIn", 200, -0.20, 0.20 );
  HB2( SsdInHid + 15, "Utgt%Xtgt SsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( SsdInHid + 16, "Vtgt%Ytgt SsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( SsdInHid + 17, "Xtgt%Ytgt SsdIn", 100, -100., 100., 100, -100, 100 );

  HB1( SsdInHid + 21, "Xemul SsdIn", 400, -100., 100. );
  HB1( SsdInHid + 22, "Yemul SsdIn", 400, -100., 100. );
  HB1( SsdInHid + 23, "Uemul SsdIn", 200, -0.20, 0.20 );
  HB1( SsdInHid + 24, "Vemul SsdIn", 200, -0.20, 0.20 );
  HB2( SsdInHid + 25, "Uemul%Xemul SsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( SsdInHid + 26, "Vemul%Yemul SsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( SsdInHid + 27, "Xemul%Yemul SsdIn", 100, -100., 100., 100, -100, 100 );

  for( int i=1; i<=NumOfLayersSsdIn; ++i ){
    HB1( SsdInHid + 100*i+11, Form("wire for LayerId : %02d [Track]", i), 100, 0., 100. );
    HB1( SsdInHid + 100*i+14, Form("Position SsdIn %d", i), 100, -250., 250. );
  }

  // BcOutSsdInTracking
  HB1( BcOutSsdInHid +  0, "#Tracks BcOutSsdIn", 10, 0., 10. );
  HB1( BcOutSsdInHid +  1, "#Hits of Track BcOutSsdIn", 20, 0., 20. );
  HB1( BcOutSsdInHid +  2, "Chisqr BcOutSsdIn", 500, 0., 50. );
  HB1( BcOutSsdInHid +  3, "LayerId BcOutSsdIn", 20, 0., 20. );
  HB1( BcOutSsdInHid +  4, "X0 BcOutSsdIn", 400, -100., 100. );
  HB1( BcOutSsdInHid +  5, "Y0 BcOutSsdIn", 400, -100., 100. );
  HB1( BcOutSsdInHid +  6, "U0 BcOutSsdIn", 200, -0.20, 0.20 );
  HB1( BcOutSsdInHid +  7, "V0 BcOutSsdIn", 200, -0.20, 0.20 );
  HB2( BcOutSsdInHid +  8, "U0%X0 BcOutSsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( BcOutSsdInHid +  9, "V0%Y0 BcOutSsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( BcOutSsdInHid + 10, "X0%Y0 BcOutSsdIn", 100, -100., 100., 100, -100, 100 );

  HB1( BcOutSsdInHid + 11, "Xtgt BcOutSsdIn", 400, -100., 100. );
  HB1( BcOutSsdInHid + 12, "Ytgt BcOutSsdIn", 400, -100., 100. );
  HB1( BcOutSsdInHid + 13, "Utgt BcOutSsdIn", 200, -0.20, 0.20 );
  HB1( BcOutSsdInHid + 14, "Vtgt BcOutSsdIn", 200, -0.20, 0.20 );
  HB2( BcOutSsdInHid + 15, "Utgt%Xtgt BcOutSsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( BcOutSsdInHid + 16, "Vtgt%Ytgt BcOutSsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( BcOutSsdInHid + 17, "Xtgt%Ytgt BcOutSsdIn", 100, -100., 100., 100, -100, 100 );

  HB1( BcOutSsdInHid + 21, "Xemul BcOutSsdIn", 400, -100., 100. );
  HB1( BcOutSsdInHid + 22, "Yemul BcOutSsdIn", 400, -100., 100. );
  HB1( BcOutSsdInHid + 23, "Uemul BcOutSsdIn", 200, -0.20, 0.20 );
  HB1( BcOutSsdInHid + 24, "Vemul BcOutSsdIn", 200, -0.20, 0.20 );
  HB2( BcOutSsdInHid + 25, "Uemul%Xemul BcOutSsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( BcOutSsdInHid + 26, "Vemul%Yemul BcOutSsdIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( BcOutSsdInHid + 27, "Xemul%Yemul BcOutSsdIn", 100, -100., 100., 100, -100, 100 );

  for( int i=1; i<=NumOfLayersBcOut + NumOfLayersSsdIn; ++i ){
    HB1( BcOutSsdInHid + 100*i+11, Form("wire for LayerId : %02d [Track]", i), 100, 0., 100. );
    HB1( BcOutSsdInHid + 100*i+14, Form("Position BcOutSsdIn %d", i), 100, -250., 250. );
    double resmin, resmax;
    if(i<=NumOfLayersBcOut){ resmin = -10.0; resmax = 10.0; }
    else                   { resmin =  -0.5; resmax =  0.5; }
    HB1( BcOutSsdInHid + 100*i+15, Form("Residual BcOutSsdIn %d", i), 200, resmin, resmax );
    HB2( BcOutSsdInHid + 100*i+16, Form("Position%%Residual BcOutSsdIn %d", i),
	 100, -250, 250, 200, resmin, resmax );
    HB2( BcOutSsdInHid + 100*i+17, Form("Xcal%%Ycal BcOutSsdIn %d", i),
	 100, -250, 250, 100, -250., 250. );
    HB1( BcOutSsdInHid + 100*i+71, Form("Residual BcOutSsdIn %d (Theta:0-15deg)", i),
	 200, resmin, resmax);
    HB1( BcOutSsdInHid + 100*i+72, Form("Residual BcOutSsdIn %d (Theta:15-30deg)", i),
	 200, resmin, resmax);
    HB1( BcOutSsdInHid + 100*i+73, Form("Residual BcOutSsdIn %d (Theta:30-45deg)", i),
	 200, resmin, resmax);
    HB1( BcOutSsdInHid + 100*i+74, Form("Residual BcOutSsdIn %d (Theta:45-deg)", i),
	 200, resmin, resmax);
  }

  // Residual
  HB1( 100, "Residual SsdIn-BcOut(SsdIn) X", 200, -100.0, 100.0 );
  HB1( 200, "Residual SsdIn-BcOut(SsdIn) Y", 200,  -20.0,  20.0 );
  HB1( 300, "Residual SsdIn-BcOut(SsdIn) U", 200,   -0.5,   0.5 );
  HB1( 400, "Residual SsdIn-BcOut(SsdIn) V", 200,   -0.2,   0.2 );

  ////////////////////////////////////////////
  //Tree
  HBTree( "tree","tree of BcOutSsdInTracking" );
  tree->Branch("evnum",    &event.evnum,    "evnum/I");
  tree->Branch("trigpat",   event.trigpat,  Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",  event.trigflag, Form("trigflag[%d]/I", NumOfSegTrig));

  // SsdIn
  tree->Branch("ssd1y0nhits",  &event.ssd1y0nhits,  "ssd1y0nhits/I");
  tree->Branch("ssd1y0de",      event.ssd1y0de,     "ssd1y0de[ssd1y0nhits]/D");
  tree->Branch("ssd1y0pos",     event.ssd1y0pos,    "ssd1y0pos[ssd1y0nhits]/D");
  tree->Branch("ssd1x0nhits",  &event.ssd1x0nhits,  "ssd1x0nhits/I");
  tree->Branch("ssd1x0de",      event.ssd1x0de,     "ssd1x0de[ssd1x0nhits]/D");
  tree->Branch("ssd1x0pos",     event.ssd1x0pos,    "ssd1x0pos[ssd1x0nhits]/D");
  tree->Branch("ssd1y1nhits",  &event.ssd1y1nhits,  "ssd1y1nhits/I");
  tree->Branch("ssd1y1de",      event.ssd1y1de,     "ssd1y1de[ssd1y1nhits]/D");
  tree->Branch("ssd1y1pos",     event.ssd1y1pos,    "ssd1y1pos[ssd1y1nhits]/D");
  tree->Branch("ssd1x1nhits",  &event.ssd1x1nhits,  "ssd1x1nhits/I");
  tree->Branch("ssd1x1de",      event.ssd1x1de,     "ssd1x1de[ssd1x1nhits]/D");
  tree->Branch("ssd1x1pos",     event.ssd1x1pos,    "ssd1x1pos[ssd1x1nhits]/D");

  // SSDT
  tree->Branch("ssdtnhits",    &event.ssdtnhits,    "ssdtnhits/I");
  tree->Branch("ssdthitpat",    event.ssdthitpat,   "ssdthitpat[ssdtnhits]/D");
  tree->Branch("ssdt",          event.ssdt,         "ssdt[ssdtnhits]/D");

  tree->Branch("ntrack",   &event.ntrack,   "ntrack/I");
  tree->Branch("chisqr",    event.chisqr,   "chisqr[ntrack]/D");
  tree->Branch("x0",        event.x0,       "x0[ntrack]/D");
  tree->Branch("y0",        event.y0,       "y0[ntrack]/D");
  tree->Branch("u0",        event.u0,       "u0[ntrack]/D");
  tree->Branch("v0",        event.v0,       "v0[ntrack]/D");

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
