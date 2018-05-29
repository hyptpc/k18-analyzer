/**
 *  file: UserEasiroc.cc
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

#define HodoCut 0 // with BH1/BH2
#define TimeCut 1 // in cluster analysis

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
class EventEasiroc : public VEvent
{
private:
  RawData      *rawData;
  DCAnalyzer   *DCAna;
  HodoAnalyzer *hodoAna;

public:
        EventEasiroc( void );
       ~EventEasiroc( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventEasiroc::EventEasiroc( void )
  : VEvent(),
    rawData(0),
    DCAna( new DCAnalyzer ),
    hodoAna( new HodoAnalyzer )
{
}

//______________________________________________________________________________
EventEasiroc::~EventEasiroc( void )
{
  if ( hodoAna ){
    delete hodoAna;
    hodoAna = 0;
  }
  if ( DCAna ){
    delete DCAna;
    DCAna   = 0;
  }
  if ( rawData ){
    delete rawData;
    rawData = 0;
  }
}

//______________________________________________________________________________
struct Event
{
  int evnum;

  int trignhits;
  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

  // BFT
  int    bft_nhits;
  int    bft_unhits;
  int    bft_dnhits;
  int    bft_uhitpat[NumOfSegBFT];
  int    bft_dhitpat[NumOfSegBFT];
  double bft_utdc[NumOfSegBFT][MaxDepth];
  double bft_dtdc[NumOfSegBFT][MaxDepth];
  double bft_utrailing[NumOfSegBFT][MaxDepth];
  double bft_dtrailing[NumOfSegBFT][MaxDepth];
  double bft_utot[NumOfSegBFT][MaxDepth];
  double bft_dtot[NumOfSegBFT][MaxDepth];
  int    bft_udepth[NumOfSegBFT];
  int    bft_ddepth[NumOfSegBFT];
  int    bft_ncl;
  int    bft_clsize[NumOfSegBFT];
  double bft_ctime[NumOfSegBFT];
  double bft_ctot[NumOfSegBFT];
  double bft_clpos[NumOfSegBFT];

  // SCH
  int    sch_nhits;
  int    sch_hitpat[NumOfSegSCH];
  double sch_tdc[NumOfSegSCH][MaxDepth];
  double sch_trailing[NumOfSegSCH][MaxDepth];
  double sch_tot[NumOfSegSCH][MaxDepth];
  int    sch_depth[NumOfSegSCH];
  int    sch_ncl;
  int    sch_clsize[NumOfSegSCH];
  double sch_ctime[NumOfSegSCH];
  double sch_ctot[NumOfSegSCH];
  double sch_clpos[NumOfSegSCH];

  // SFT-V
  int    sftv_nhits;
  int    sftv_hitpat[NumOfSegSFT_UV];
  double sftv_tdc[NumOfSegSFT_UV][MaxDepth];
  double sftv_trailing[NumOfSegSFT_UV][MaxDepth];
  double sftv_tot[NumOfSegSFT_UV][MaxDepth];
  int    sftv_depth[NumOfSegSFT_UV];
  int    sftv_ncl;
  int    sftv_clsize[NumOfSegSFT_UV];
  double sftv_ctime[NumOfSegSFT_UV];
  double sftv_ctot[NumOfSegSFT_UV];
  double sftv_clpos[NumOfSegSFT_UV];

  // SFT-U
  int    sftu_nhits;
  int    sftu_hitpat[NumOfSegSFT_UV];
  double sftu_tdc[NumOfSegSFT_UV][MaxDepth];
  double sftu_trailing[NumOfSegSFT_UV][MaxDepth];
  double sftu_tot[NumOfSegSFT_UV][MaxDepth];
  int    sftu_depth[NumOfSegSFT_UV];
  int    sftu_ncl;
  int    sftu_clsize[NumOfSegSFT_UV];
  double sftu_ctime[NumOfSegSFT_UV];
  double sftu_ctot[NumOfSegSFT_UV];
  double sftu_clpos[NumOfSegSFT_UV];

  // SFT-X
  int    sftx_nhits;
  int    sftx_unhits;
  int    sftx_dnhits;
  int    sftx_uhitpat[NumOfSegSFT_X];
  int    sftx_dhitpat[NumOfSegSFT_X];
  double sftx_utdc[NumOfSegSFT_X][MaxDepth];
  double sftx_dtdc[NumOfSegSFT_X][MaxDepth];
  double sftx_utrailing[NumOfSegSFT_X][MaxDepth];
  double sftx_dtrailing[NumOfSegSFT_X][MaxDepth];
  double sftx_utot[NumOfSegSFT_X][MaxDepth];
  double sftx_dtot[NumOfSegSFT_X][MaxDepth];
  int    sftx_udepth[NumOfSegSFT_X];
  int    sftx_ddepth[NumOfSegSFT_X];
  int    sftx_ncl;
  int    sftx_clsize[NumOfSegSFT_X];
  double sftx_ctime[NumOfSegSFT_X];
  double sftx_ctot[NumOfSegSFT_X];
  double sftx_clpos[NumOfSegSFT_X];
};

//______________________________________________________________________________
namespace root
{
  Event  event;
  TH1   *h[MaxHist];
  TTree *tree;
  enum eDetHid
    {
      BFTHid  = 10000,
      SCHHid  = 20000,
      SFTVHid = 30000, SFTUHid = 40000, SFTXHid = 50000
    };
}

//______________________________________________________________________________
bool
EventEasiroc::ProcessingBegin( void )
{
  InitializeEvent();
  return true;
}

//______________________________________________________________________________
bool
EventEasiroc::ProcessingNormal( void )
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
  static const double MinTdcBFT  = gUser.GetParameter("TdcBFT",  0);
  static const double MaxTdcBFT  = gUser.GetParameter("TdcBFT",  1);
  static const double MinTdcSCH  = gUser.GetParameter("TdcSCH",  0);
  static const double MaxTdcSCH  = gUser.GetParameter("TdcSCH",  1);
  static const double MinTdcSFT  = gUser.GetParameter("TdcSFT",  0);
  static const double MaxTdcSFT  = gUser.GetParameter("TdcSFT",  1);
#if TimeCut
  static const double MinTimeBFT = gUser.GetParameter("TimeBFT", 0);
  static const double MaxTimeBFT = gUser.GetParameter("TimeBFT", 1);
  static const double MinTimeSCH = gUser.GetParameter("TimeSCH", 0);
  static const double MaxTimeSCH = gUser.GetParameter("TimeSCH", 1);
  static const double MinTimeSFT = gUser.GetParameter("TimeSFT", 0);
  static const double MaxTimeSFT = gUser.GetParameter("TimeSFT", 1);
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
	HF1( 10, seg-1 );
	HF1( 10+seg, tdc );
      }
    }
    event.trignhits = trignhits;
  }

  // if( event.trigflag[SpillEndFlag] ) return true;

  HF1(1, 1);

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

  ////////// BFT
  {
    hodoAna->DecodeBFTHits(rawData);
    // Fiber Hit
    int unhits = 0;
    int dnhits = 0;
    for(int plane = 0; plane<NumOfPlaneBFT; ++plane){
      int nh = hodoAna->GetNHitsBFT(plane);
      enum { U, D };
      for( int i=0; i<nh; ++i ){
	const FiberHit* hit = hodoAna->GetHitBFT(plane, i);
	if(!hit) continue;
	int mhit_l  = hit->GetNLeading();
	int mhit_t  = hit->GetNTrailing();
	int seg   = hit->SegmentId();
	if(plane==U) event.bft_udepth[seg] = mhit_l;
	if(plane==D) event.bft_ddepth[seg] = mhit_l;

	int  prev = 0;
	bool hit_flag = false;

	// raw leading data
	for( int m=0; m<mhit_l; ++m ){
	  if(mhit_l > MaxDepth) break;
	  double leading  = hit->GetLeading(m);

	  if(leading==prev) continue;
	  prev = leading;
	  HF1( BFTHid +plane + 6, leading );
	  HF2( BFTHid +plane + 10, seg, leading );
	  HF1( BFTHid +1000*(1+plane)+seg+1, leading );

	  if(plane==U){
	    event.bft_utdc[seg][m]      = leading;
	  }
	  if(plane==D){
	    event.bft_dtdc[seg][m]      = leading;
	  }

	  if( MinTdcBFT<leading && leading<MaxTdcBFT ){
	    hit_flag = true;
	  }
	}// for(m)

	// raw leading data
	for( int m=0; m<mhit_t; ++m ){
	  if(mhit_t > MaxDepth) break;
	  double trailing = hit->GetTrailing(m);

	  if(plane==U){
	    event.bft_utrailing[seg][m] = trailing;
	  }
	  if(plane==D){
	    event.bft_dtrailing[seg][m] = trailing;
	  }

	}// for(m)

	int mhit_pair  = hit->GetNPair();	
	// pair data
	for( int m=0; m<mhit_pair; ++m ){
	  if(mhit_pair > MaxDepth) break;
	  double time     = hit->GetTime(m);
	  double ctime    = hit->GetCTime(m);
	  double width    = hit->GetWidth(m);

	  HF1( BFTHid +plane+8, width );
	  HF2( BFTHid +plane+12, seg, width );
	  HF1( BFTHid +plane+21, time );
	  HF2( BFTHid +plane+23, width, time );
	  HF1( BFTHid +plane+31, ctime );
	  HF2( BFTHid +plane+33, width, ctime );
	  HF1( BFTHid +1000*(plane+3)+seg+1, width );
	  if( -10.<time && time<10. ){
	    HF2( BFTHid +1000*(plane+5)+seg+1, width, time );
	    HF2( BFTHid +1000*(plane+7)+seg+1, width, ctime );
	  }
	  if(plane==U){
	    event.bft_utot[seg][m]      = width;
	  }
	  if(plane==D){
	    event.bft_dtot[seg][m]      = width;
	  }

	}
	if(hit_flag){
	  HF1( BFTHid +plane+4, seg+0.5);
	  if(plane==U) event.bft_uhitpat[unhits++] = seg;
	  if(plane==D) event.bft_dhitpat[dnhits++] = seg;
	}
      }
    }
    HF1( BFTHid +1, unhits);
    HF1( BFTHid +2, dnhits);
    HF1( BFTHid +3, unhits + dnhits);
    event.bft_unhits = unhits;
    event.bft_dnhits = dnhits;
    event.bft_nhits  = unhits + dnhits;

    // Fiber Cluster
#if TimeCut
    hodoAna->TimeCutBFT(MinTimeBFT, MaxTimeBFT);
#endif
    int ncl = hodoAna->GetNClustersBFT();
    event.bft_ncl = ncl;
    HF1( BFTHid +101, ncl );
    for(int i=0; i<ncl; ++i){
      FiberCluster *cl = hodoAna->GetClusterBFT(i);
      if(!cl) continue;
      double clsize = cl->ClusterSize();
      double ctime  = cl->CMeanTime();
      double ctot   = cl->Width();
      double pos    = cl->MeanPosition();
      event.bft_clsize[i] = clsize;
      event.bft_ctime[i]  = ctime;
      event.bft_ctot[i]   = ctot;
      event.bft_clpos[i]  = pos;
      HF1( BFTHid +102, clsize );
      HF1( BFTHid +103, ctime );
      HF1( BFTHid +104, ctot );
      HF2( BFTHid +105, ctot, ctime );
      HF1( BFTHid +106, pos );
    }
  }

  ////////// SCH
  {
    hodoAna->DecodeSCHHits(rawData);
    int nh = hodoAna->GetNHitsSCH();
    int sch_nhits = 0;
    for(int i=0; i<nh; ++i){
      FiberHit* hit = hodoAna->GetHitSCH(i);
      if(!hit) continue;
      int mhit_l = hit->GetNLeading();
      int mhit_t = hit->GetNTrailing();
      int seg    = hit->SegmentId();
      event.sch_depth[seg] = mhit_l;
      int  prev     = 0;
      bool hit_flag = false;

      for( int m=0; m<mhit_l; ++m ){
	if(mhit_l > MaxDepth) break;
	double leading  = hit->GetLeading(m);
	if( leading==prev ) continue;
	prev = leading;
	HF1( SCHHid +3, leading );
	HF2( SCHHid +5, seg, leading );
	HF1( SCHHid +1000+seg+1, leading );

	event.sch_tdc[seg][m]      = leading;

	if( MinTdcSCH<leading && leading<MaxTdcSCH ){
	  hit_flag = true;
	}
      }// for(m)

      for( int m=0; m<mhit_t; ++m ){
	if(mhit_t > MaxDepth) break;
	double trailing = hit->GetTrailing(m);
	event.sch_trailing[seg][m] = trailing;
      }// for(m)

      int mhit_pair = hit->GetNPair();
      for( int m=0; m<mhit_pair; ++m ){
	if(mhit_pair > MaxDepth) break;

	double time     = hit->GetTime(m);
	double ctime    = hit->GetCTime(m);
	double width    = hit->GetWidth(m);

	HF1( SCHHid +4, width );
	HF2( SCHHid +6, seg, width );
	HF1( SCHHid +21, time );
	HF2( SCHHid +22, width, time );
	HF1( SCHHid +31, ctime );
	HF2( SCHHid +32, width, ctime );
	HF1( SCHHid +2000+seg+1, width );
	if( -10.<time && time<10. ){
	  HF2( SCHHid +3000+seg+1, width, time );
	  HF2( SCHHid +4000+seg+1, width, ctime );
	}

	event.sch_tot[seg][m]      = width;

      }//for(m)
      if(hit_flag){
	HF1( SCHHid +2, seg+0.5);
	event.sch_hitpat[sch_nhits++] = seg;
      }
    }
    HF1( SCHHid +1, sch_nhits );
    event.sch_nhits = sch_nhits;

    // Fiber Cluster
#if TimeCut
    hodoAna->TimeCutSCH( MinTimeSCH, MaxTimeSCH );
#endif
    int ncl = hodoAna->GetNClustersSCH();
    event.sch_ncl = ncl;
    HF1( SCHHid +101, ncl );
    for(int i=0; i<ncl; ++i){
      FiberCluster *cl = hodoAna->GetClusterSCH(i);
      if(!cl) continue;
      double clsize = cl->ClusterSize();
      double ctime  = cl->CMeanTime();
      double ctot   = cl->Width();
      double pos    = cl->MeanPosition();
      event.sch_clsize[i] = clsize;
      event.sch_ctime[i]  = ctime;
      event.sch_ctot[i]   = ctot;
      event.sch_clpos[i]  = pos;
      HF1( SCHHid +102, clsize );
      HF1( SCHHid +103, ctime );
      HF1( SCHHid +104, ctot );
      HF2( SCHHid +105, ctot, ctime );
      HF1( SCHHid +106, pos );
    }
  }

  hodoAna->DecodeSFTHits(rawData);
  hodoAna->WidthCutSFT( 0, 40. , 100.);
  hodoAna->WidthCutSFT( 1, 40. , 100.);
  hodoAna->WidthCutSFT( 2, 40. , 100.);

  ////////// SFT-V
  {
    // Fiber Hit
    int nhits = 0;
    int nh = hodoAna->GetNHitsSFT(SFT_V);
    for( int i=0; i<nh; ++i ){
      const FiberHit* hit = hodoAna->GetHitSFT(SFT_V, i);
      if(!hit) continue;
      int mhit_l = hit->GetNLeading();
      int mhit_t = hit->GetNTrailing();
      int seg    = hit->SegmentId();
      event.sftv_depth[seg] = mhit_l;

      int  prev = 0;
      bool hit_flag = false;

      // raw leading data
      for( int m=0; m<mhit_l; ++m ){
	if(mhit_l > MaxDepth) break;
	double leading  = hit->GetLeading(m);

	if(leading==prev) continue;
	prev = leading;
	HF1( SFTVHid +6, leading );
	HF2( SFTVHid +10, seg, leading );
	HF1( SFTVHid +1000*(1)+seg+1, leading );

	event.sftv_tdc[seg][m]      = leading;

	if( MinTdcSFT<leading && leading<MaxTdcSFT ){
	  hit_flag = true;
	}
      }// for(m)

      // raw trailing data
      for( int m=0; m<mhit_t; ++m ){
	if(mhit_t > MaxDepth) break;
	double trailing = hit->GetTrailing(m);
	event.sftv_trailing[seg][m] = trailing;
      }// for(m)

      int mhit_pair = hit->GetNPair();
      // pair data
      for( int m=0; m<mhit_pair; ++m ){
	if(mhit_pair > MaxDepth) break;
	double time     = hit->GetTime(m);
	double ctime    = hit->GetCTime(m);
	double width    = hit->GetWidth(m);

	HF1( SFTVHid +8, width );
	HF2( SFTVHid +12, seg, width );
	HF1( SFTVHid +21, time );
	HF2( SFTVHid +23, width, time );
	HF1( SFTVHid +31, ctime );
	HF2( SFTVHid +33, width, ctime );
	HF1( SFTVHid +1000*(3)+seg+1, width );

	if( -10.<time && time<10. ){
	  HF2( SFTVHid +1000*(5)+seg+1, width, time );
	  HF2( SFTVHid +1000*(7)+seg+1, width, ctime );
	}

	event.sftv_tot[seg][m]      = width;

      }// for(m)

      if(hit_flag){
	HF1( SFTVHid +4, seg+0.5);
	 event.sftv_hitpat[nhits++] = seg;
      }
    }// for(i)
    HF1( SFTVHid +1, nhits);
    event.sftv_nhits  = nhits;

    // Fiber Cluster
#if TimeCut
    hodoAna->TimeCutSFT(0, MinTimeSFT, MaxTimeSFT);
#endif

    int ncl = hodoAna->GetNClustersSFT(0);
    event.sftv_ncl = ncl;
    HF1( SFTVHid +101, ncl );
    for(int i=0; i<ncl; ++i){
      FiberCluster *cl = hodoAna->GetClusterSFT(0, i);
      if(!cl) continue;
      double clsize = cl->ClusterSize();
      double ctime  = cl->CMeanTime();
      double ctot   = cl->Width();
      double pos    = cl->MeanPosition();
      event.sftv_clsize[i] = clsize;
      event.sftv_ctime[i]  = ctime;
      event.sftv_ctot[i]   = ctot;
      event.sftv_clpos[i]  = pos;
      if(isnan(ctot)) std::cout << "nan" << std::endl;
      HF1( SFTVHid +102, clsize );
      HF1( SFTVHid +103, ctime );
      HF1( SFTVHid +104, ctot );
      HF2( SFTVHid +105, ctot, ctime );
      HF1( SFTVHid +106, pos );
    }
  }

  ////////// SFT-U
  {
    // Fiber Hit
    int nhits = 0;
    int nh = hodoAna->GetNHitsSFT(SFT_U);
    for( int i=0; i<nh; ++i ){
      const FiberHit* hit = hodoAna->GetHitSFT(SFT_U, i);
      if(!hit) continue;
      int mhit_l = hit->GetNLeading();
      int mhit_t = hit->GetNTrailing();
      int seg    = hit->SegmentId();
      event.sftu_depth[seg] = mhit_l;

      int  prev = 0;
      bool hit_flag = false;

      // raw leading data
      for( int m=0; m<mhit_l; ++m ){
	if(mhit_l > MaxDepth) break;
	double leading  = hit->GetLeading(m);

	if(leading==prev) continue;
	prev = leading;
	HF1( SFTUHid +6, leading );
	HF2( SFTUHid +10, seg, leading );
	HF1( SFTUHid +1000*(1)+seg+1, leading );

	event.sftu_tdc[seg][m]      = leading;

	if( MinTdcSFT<leading && leading<MaxTdcSFT ){
	  hit_flag = true;
	}
      }// for(m)

      // raw trailing data
      for( int m=0; m<mhit_t; ++m ){
	if(mhit_t > MaxDepth) break;
	double trailing = hit->GetTrailing(m);
	event.sftu_trailing[seg][m] = trailing;
      }// for(m)

      int mhit_pair = hit->GetNPair();
      // pair data
      for( int m=0; m<mhit_pair; ++m ){
	if(mhit_pair > MaxDepth) break;
	double time     = hit->GetTime(m);
	double ctime    = hit->GetCTime(m);
	double width    = hit->GetWidth(m);

	HF1( SFTUHid +8, width );
	HF2( SFTUHid +12, seg, width );
	HF1( SFTUHid +21, time );
	HF2( SFTUHid +23, width, time );
	HF1( SFTUHid +31, ctime );
	HF2( SFTUHid +33, width, ctime );
	HF1( SFTUHid +1000*(3)+seg+1, width );

	if( -10.<time && time<10. ){
	  HF2( SFTUHid +1000*(5)+seg+1, width, time );
	  HF2( SFTUHid +1000*(7)+seg+1, width, ctime );
	}

	event.sftu_tot[seg][m]      = width;

      }// for(m)

      if(hit_flag){
	HF1( SFTUHid +4, seg+0.5);
	event.sftu_hitpat[nhits++] = seg;
      }
    }// for(i)
    HF1( SFTUHid +1, nhits);
    event.sftu_nhits  = nhits;

    // Fiber Cluster
#if TimeCut
    hodoAna->TimeCutSFT(1, MinTimeSFT, MaxTimeSFT);
#endif
    int ncl = hodoAna->GetNClustersSFT(1);
    event.sftu_ncl = ncl;
    HF1( SFTUHid +101, ncl );
    for(int i=0; i<ncl; ++i){
      FiberCluster *cl = hodoAna->GetClusterSFT(1, i);
      if(!cl) continue;
      double clsize = cl->ClusterSize();
      double ctime  = cl->CMeanTime();
      double ctot   = cl->Width();
      double pos    = cl->MeanPosition();
      event.sftu_clsize[i] = clsize;
      event.sftu_ctime[i]  = ctime;
      event.sftu_ctot[i]   = ctot;
      event.sftu_clpos[i]  = pos;
      HF1( SFTUHid +102, clsize );
      HF1( SFTUHid +103, ctime );
      HF1( SFTUHid +104, ctot );
      HF2( SFTUHid +105, ctot, ctime );
      HF1( SFTUHid +106, pos );
    }
  }

  ////////// SFT-X
  {
    // Fiber Hit
    int unhits = 0;
    int dnhits = 0;
    for(int p = SFT_X1; p<NumOfPlaneSFT; ++p){
      int nh = hodoAna->GetNHitsSFT(p);
      enum { U, D };
      for( int i=0; i<nh; ++i ){
	const FiberHit* hit = hodoAna->GetHitSFT(p, i);
	if(!hit) continue;
	int plane  = p-SFT_X1;
	int mhit_l = hit->GetNLeading();
	int mhit_t = hit->GetNTrailing();
	int seg    = hit->SegmentId();
	if(plane==U) event.sftx_udepth[seg] = mhit_l;
	if(plane==D) event.sftx_ddepth[seg] = mhit_l;

	int  prev = 0;
	bool hit_flag = false;

	// raw leading data
	for( int m=0; m<mhit_l; ++m ){
	  if(mhit_l > MaxDepth) break;
	  double leading  = hit->GetLeading(m);

	  if(leading==prev) continue;
	  prev = leading;
	  HF1( SFTXHid +plane+6, leading );
	  HF2( SFTXHid +plane+10, seg, leading );
	  HF1( SFTXHid +1000*(plane+1)+seg+1, leading );

	  if(plane==U){
	    event.sftx_utdc[seg][m]      = leading;
	  }
	  if(plane==D){
	    event.sftx_dtdc[seg][m]      = leading;
	  }

	  if( MinTdcSFT<leading && leading<MaxTdcSFT ){
	    hit_flag = true;
	  }
	}// for(m)

	// raw trailing data
	for( int m=0; m<mhit_t; ++m ){
	  if(mhit_t > MaxDepth) break;
	  double trailing = hit->GetTrailing(m);
	  if(plane==U){
	    event.sftx_utrailing[seg][m] = trailing;
	  }
	  if(plane==D){
	    event.sftx_dtrailing[seg][m] = trailing;
	  }
	}// for(m)

	int mhit_pair = hit->GetNPair();
	// pair data
	for( int m=0; m<mhit_pair; ++m ){
	  if(mhit_pair > MaxDepth) break;
	  double time     = hit->GetTime(m);
	  double ctime    = hit->GetCTime(m);
	  double width    = hit->GetWidth(m);

	  HF1( SFTXHid +plane+8, width );

	  HF2( SFTXHid +plane+12, seg, width );
	  HF1( SFTXHid +plane+21, time );
	  HF2( SFTXHid +plane+23, width, time );
	  HF1( SFTXHid +plane+31, ctime );
	  HF2( SFTXHid +plane+33, width, ctime );

	  HF1( SFTXHid +1000*(plane+3)+seg+1, width );
	  if( -10.<time && time<10. ){
	    HF2( SFTXHid +1000*(plane+5)+seg+1, width, time );
	    HF2( SFTXHid +1000*(plane+7)+seg+1, width, ctime );
	  }
	  if(plane==U){
	    event.sftx_utot[seg][m]      = width;
	  }
	  if(plane==D){
	    event.sftx_dtot[seg][m]      = width;
	  }
	}// for(m)
	
	if(hit_flag){
	  HF1( SFTXHid +plane+4, seg+0.5);
	  if(plane==U) event.sftx_uhitpat[unhits++] = seg;
	  if(plane==D) event.sftx_dhitpat[dnhits++] = seg;
	}
      }// for(i)
    }// for(p)
    HF1( SFTXHid +1, unhits);
    HF1( SFTXHid +2, dnhits);
    HF1( SFTXHid +3, unhits + dnhits);
    event.sftx_unhits = unhits;
    event.sftx_dnhits = dnhits;
    event.sftx_nhits  = unhits + dnhits;

    // Fiber Cluster
#if TimeCut
    hodoAna->TimeCutSFT(2, MinTimeSFT, MaxTimeSFT);
#endif
    int ncl = hodoAna->GetNClustersSFT(2);
    event.sftx_ncl = ncl;
    HF1( SFTXHid +101, ncl );
    for(int i=0; i<ncl; ++i){
      FiberCluster *cl = hodoAna->GetClusterSFT(2, i);
      if(!cl) continue;
      double clsize = cl->ClusterSize();
      double ctime  = cl->CMeanTime();
      double ctot   = cl->Width();
      double pos    = cl->MeanPosition();
      event.sftx_clsize[i] = clsize;
      event.sftx_ctime[i]  = ctime;
      event.sftx_ctot[i]   = ctot;
      event.sftx_clpos[i]  = pos;
      HF1( SFTXHid +102, clsize );
      HF1( SFTXHid +103, ctime );
      HF1( SFTXHid +104, ctot );
      HF2( SFTXHid +105, ctot, ctime );
      HF1( SFTXHid +106, pos );
    }
  }

  return true;
}

//______________________________________________________________________________
bool
EventEasiroc::ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//______________________________________________________________________________
void
EventEasiroc::InitializeEvent( void )
{
  event.evnum      = 0;
  event.bft_nhits  = 0;
  event.bft_unhits = 0;
  event.bft_dnhits = 0;
  event.bft_ncl    = 0;
  event.sch_nhits  = 0;
  event.sch_ncl    = 0;

  event.sftv_nhits  = 0;
  event.sftv_nhits  = 0;
  event.sftv_ncl    = 0;

  event.sftu_nhits  = 0;
  event.sftu_nhits  = 0;
  event.sftu_ncl    = 0;

  event.sftx_nhits  = 0;
  event.sftx_unhits = 0;
  event.sftx_dnhits = 0;
  event.sftx_ncl    = 0;

  for( int it=0; it<NumOfSegTrig; it++){
    event.trigpat[it]  = -1;
    event.trigflag[it] = -1;
  }

  for( int it=0; it<NumOfSegBFT; it++){
    event.bft_uhitpat[it] = -999;
    event.bft_dhitpat[it] = -999;
    event.bft_udepth[it]  = 0;
    event.bft_ddepth[it]  = 0;
    for( int that=0; that<MaxDepth; that++){
      event.bft_utdc[it][that] = -999.;
      event.bft_dtdc[it][that] = -999.;
      event.bft_utrailing[it][that] = -999.;
      event.bft_dtrailing[it][that] = -999.;
      event.bft_utot[it][that] = -999.;
      event.bft_dtot[it][that] = -999.;
    }
    event.bft_clsize[it] = -999;
    event.bft_ctime[it]  = -999.;
    event.bft_ctot[it]   = -999.;
    event.bft_clpos[it]  = -999.;
  }

  for( int it=0; it<NumOfSegSCH; it++){
    event.sch_hitpat[it] = -999;
    event.sch_depth[it]  = 0;
    for( int that=0; that<MaxDepth; that++){
      event.sch_tdc[it][that] = -999.;
      event.sch_trailing[it][that] = -999.;
      event.sch_tot[it][that] = -999.;
    }
    event.sch_clsize[it] = -999;
    event.sch_ctime[it]  = -999.;
    event.sch_ctot[it]   = -999.;
    event.sch_clpos[it]  = -999.;
  }

  for( int it=0; it<NumOfSegSFT_UV; it++){
    event.sftv_hitpat[it] = -999;
    event.sftv_depth[it]  = 0;
    event.sftu_hitpat[it] = -999;
    event.sftu_depth[it]  = 0;
    for( int that=0; that<MaxDepth; that++){
      event.sftv_tdc[it][that] = -999.;
      event.sftv_trailing[it][that] = -999.;
      event.sftv_tot[it][that] = -999.;

      event.sftu_tdc[it][that] = -999.;
      event.sftu_trailing[it][that] = -999.;
      event.sftu_tot[it][that] = -999.;
    }
    event.sftv_clsize[it] = -999;
    event.sftv_ctime[it]  = -999.;
    event.sftv_ctot[it]   = -999.;
    event.sftv_clpos[it]  = -999.;

    event.sftu_clsize[it] = -999;
    event.sftu_ctime[it]  = -999.;
    event.sftu_ctot[it]   = -999.;
    event.sftu_clpos[it]  = -999.;
  }

  for( int it=0; it<NumOfSegSFT_X; it++){
    event.sftx_uhitpat[it] = -999;
    event.sftx_dhitpat[it] = -999;
    event.sftx_udepth[it]  = 0;
    event.sftx_ddepth[it]  = 0;
    for( int that=0; that<MaxDepth; that++){
      event.sftx_utdc[it][that] = -999.;
      event.sftx_dtdc[it][that] = -999.;
      event.sftx_utrailing[it][that] = -999.;
      event.sftx_dtrailing[it][that] = -999.;
      event.sftx_utot[it][that] = -999.;
      event.sftx_dtot[it][that] = -999.;
    }
    event.sftx_clsize[it] = -999;
    event.sftx_ctime[it]  = -999.;
    event.sftx_ctot[it]   = -999.;
    event.sftx_clpos[it]  = -999.;
  }

}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventEasiroc;
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
  HB1( 10, "Trigger HitPat", NumOfSegTrig, 0., double(NumOfSegTrig) );
  for(int i=0; i<NumOfSegTrig; ++i){
    HB1( 10+i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000 );
  }

  //BFT
  HB1( BFTHid + 1, "BFT Nhits U",   NumOfSegBFT, 0., (double)NumOfSegBFT );
  HB1( BFTHid + 2, "BFT Nhits D",   NumOfSegBFT, 0., (double)NumOfSegBFT );
  HB1( BFTHid + 3, "BFT Nhits",     NumOfSegBFT, 0., (double)NumOfSegBFT );
  HB1( BFTHid + 4, "BFT Hitpat U",  NumOfSegBFT, 0., (double)NumOfSegBFT );
  HB1( BFTHid + 5, "BFT Hitpat D",  NumOfSegBFT, 0., (double)NumOfSegBFT );
  HB1( BFTHid + 6, "BFT Tdc U",      NbinTdc, MinTdc, MaxTdc );
  HB1( BFTHid + 7, "BFT Tdc D",      NbinTdc, MinTdc, MaxTdc );
  HB1( BFTHid + 8, "BFT Tot U",      NbinTot, MinTot, MaxTot );
  HB1( BFTHid + 9, "BFT Tot D",      NbinTot, MinTot, MaxTot );
  HB2( BFTHid +10, "BFT Tdc U%Seg",
       NumOfSegBFT, 0., (double)NumOfSegBFT, NbinTdc, MinTdc, MaxTdc );
  HB2( BFTHid +11, "BFT Tdc D%Seg",
       NumOfSegBFT, 0., (double)NumOfSegBFT, NbinTdc, MinTdc, MaxTdc );
  HB2( BFTHid +12, "BFT Tot U%Seg",
       NumOfSegBFT, 0., (double)NumOfSegBFT, NbinTot, MinTot, MaxTot );
  HB2( BFTHid +13, "BFT Tot D%Seg",
       NumOfSegBFT, 0., (double)NumOfSegBFT, NbinTot, MinTot, MaxTot );
  for(int i=0; i<NumOfSegBFT; i++){
    HB1( BFTHid +1000+i+1, Form("BFT Tdc U-%d", i+1), NbinTdc, MinTdc, MaxTdc );
    HB1( BFTHid +2000+i+1, Form("BFT Tdc D-%d", i+1), NbinTdc, MinTdc, MaxTdc );
    HB1( BFTHid +3000+i+1, Form("BFT Tot U-%d", i+1), NbinTot, MinTot, MaxTot );
    HB1( BFTHid +4000+i+1, Form("BFT Tot D-%d", i+1), NbinTot, MinTot, MaxTot );
    HB2( BFTHid +5000+i+1, Form("BFT Time/Tot U-%d", i+1),
	 NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
    HB2( BFTHid +6000+i+1, Form("BFT Time/Tot D-%d", i+1),
	 NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
    HB2( BFTHid +7000+i+1, Form("BFT CTime/Tot U-%d", i+1),
	 NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
    HB2( BFTHid +8000+i+1, Form("BFT CTime/Tot D-%d", i+1),
	 NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  }
  HB1( BFTHid +21, "BFT Time U",     NbinTime, MinTime, MaxTime );
  HB1( BFTHid +22, "BFT Time D",     NbinTime, MinTime, MaxTime );
  HB2( BFTHid +23, "BFT Time/Tot U", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  HB2( BFTHid +24, "BFT Time/Tot D", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  HB1( BFTHid +31, "BFT CTime U",     NbinTime, MinTime, MaxTime );
  HB1( BFTHid +32, "BFT CTime D",     NbinTime, MinTime, MaxTime );
  HB2( BFTHid +33, "BFT CTime/Tot U", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  HB2( BFTHid +34, "BFT CTime/Tot D", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );

  HB1( BFTHid +101, "BFT NCluster", 100, 0, 100);
  HB1( BFTHid +102, "BFT Cluster Size", 5, 0, 5);
  HB1( BFTHid +103, "BFT CTime (Cluster)", NbinTime, MinTime, MaxTime );
  HB1( BFTHid +104, "BFT Tot (Cluster)", NbinTot, MinTot, MaxTot );
  HB2( BFTHid +105, "BFT CTime%Tot (Cluster)",
       NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  HB1( BFTHid +106, "BFT Cluster Position",
       NumOfSegBFT, -0.5*(double)NumOfSegBFT, 0.5*(double)NumOfSegBFT);

  //SCH
  HB1( SCHHid +1, "SCH Nhits",     NumOfSegSCH, 0., (double)NumOfSegSCH );
  HB1( SCHHid +2, "SCH Hitpat",    NumOfSegSCH, 0., (double)NumOfSegSCH );
  HB1( SCHHid +3, "SCH Tdc",      NbinTdc, MinTdc, MaxTdc );
  HB1( SCHHid +4, "SCH Tot",      NbinTot, MinTot, MaxTot );
  HB2( SCHHid +5, "SCH Tdc U%Seg",
       NumOfSegSCH, 0., (double)NumOfSegSCH, NbinTdc, MinTdc, MaxTdc );
  HB2( SCHHid +6, "SCH Tot U%Seg",
       NumOfSegSCH, 0., (double)NumOfSegSCH, NbinTot, MinTot, MaxTot );
  HB2( SCHHid +13, "SCH Tot D%Seg",
       NumOfSegSCH, 0., (double)NumOfSegSCH, NbinTot, MinTot, MaxTot );
  for(int i=0; i<NumOfSegSCH; i++){
    HB1( SCHHid +1000+i+1, Form("SCH Tdc %d", i+1), NbinTdc, MinTdc, MaxTdc );
    HB1( SCHHid +2000+i+1, Form("SCH Tot %d", i+1), NbinTot, MinTot, MaxTot );
    HB2( SCHHid +3000+i+1, Form("SCH Time/Tot %d", i+1),
	 NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
    HB2( SCHHid +4000+i+1, Form("SCH CTime/Tot %d", i+1),
	 NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  }
  HB1( SCHHid +21, "SCH Time",     NbinTime, MinTime, MaxTime );
  HB2( SCHHid +22, "SCH Time/Tot", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  HB1( SCHHid +31, "SCH CTime",     NbinTime, MinTime, MaxTime );
  HB2( SCHHid +32, "SCH CTime/Tot", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );

  HB1( SCHHid +101, "SCH NCluster", 20, 0, 20);
  HB1( SCHHid +102, "SCH Cluster Size", 5, 0, 5);
  HB1( SCHHid +103, "SCH CTime (Cluster)", NbinTime, MinTime, MaxTime );
  HB1( SCHHid +104, "SCH Tot (Cluster)", NbinTot, MinTot, MaxTot );
  HB2( SCHHid +105, "SCH CTime%Tot (Cluster)",
       NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  HB1( SCHHid +106, "SCH Cluster Position",
       NumOfSegSCH, -0.5*(double)NumOfSegSCH, 0.5*(double)NumOfSegSCH);

  //SFT-V
  HB1( SFTVHid + 1, "SFTV Nhits",     NumOfSegSFT_UV, 0., (double)NumOfSegSFT_UV );
  HB1( SFTVHid + 4, "SFTV Hitpat",    NumOfSegSFT_UV, 0., (double)NumOfSegSFT_UV );
  HB1( SFTVHid + 6, "SFTV Tdc",        NbinTdc, MinTdc, MaxTdc );
  HB1( SFTVHid + 8, "SFTV Tot",        NbinTot, MinTot, MaxTot );
  HB2( SFTVHid +10, "SFTV Tdc Seg",
       NumOfSegSFT_UV, 0., (double)NumOfSegSFT_UV, NbinTdc, MinTdc, MaxTdc );
  HB2( SFTVHid +12, "SFTV Tot Seg",
       NumOfSegSFT_UV, 0., (double)NumOfSegSFT_UV, NbinTot, MinTot, MaxTot );
  for(int i=0; i<NumOfSegSFT_UV; i++){
    HB1( SFTVHid +1000+i+1, Form("SFTV Tdc -%d", i+1), NbinTdc, MinTdc, MaxTdc );
    HB1( SFTVHid +3000+i+1, Form("SFTV Tot -%d", i+1), NbinTot, MinTot, MaxTot );
    HB2( SFTVHid +5000+i+1, Form("SFTV Time/Tot -%d", i+1),
	 NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
    HB2( SFTVHid +7000+i+1, Form("SFTV CTime/Tot -%d", i+1),
	 NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  }
  HB1( SFTVHid +21, "SFTV Time",        NbinTime, MinTime, MaxTime );
  HB2( SFTVHid +23, "SFTV Time/Tot",    NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  HB1( SFTVHid +31, "SFTV CTime",       NbinTime, MinTime, MaxTime );
  HB2( SFTVHid +33, "SFTV CTime/Tot",   NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );

  HB1( SFTVHid +101, "SFTV NCluster", 100, 0, 100);
  HB1( SFTVHid +102, "SFTV Cluster Size", 5, 0, 5);
  HB1( SFTVHid +103, "SFTV CTime (Cluster)", NbinTime, MinTime, MaxTime );
  HB1( SFTVHid +104, "SFTV Tot (Cluster)", NbinTot, MinTot, MaxTot );
  HB2( SFTVHid +105, "SFTV CTime%Tot (Cluster)",
       NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  HB1( SFTVHid +106, "SFTV Cluster Position",
       NumOfSegSFT_UV, -0.5*(double)NumOfSegSFT_UV, 0.5*(double)NumOfSegSFT_UV);

  //SFT-U
  HB1( SFTUHid + 1, "SFTU Nhits",     NumOfSegSFT_UV, 0., (double)NumOfSegSFT_UV );
  HB1( SFTUHid + 4, "SFTU Hitpat",    NumOfSegSFT_UV, 0., (double)NumOfSegSFT_UV );
  HB1( SFTUHid + 6, "SFTU Tdc",        NbinTdc, MinTdc, MaxTdc );
  HB1( SFTUHid + 8, "SFTU Tot",        NbinTot, MinTot, MaxTot );
  HB2( SFTUHid +10, "SFTU Tdc Seg",
       NumOfSegSFT_UV, 0., (double)NumOfSegSFT_UV, NbinTdc, MinTdc, MaxTdc );
  HB2( SFTUHid +12, "SFTU Tot Seg",
       NumOfSegSFT_UV, 0., (double)NumOfSegSFT_UV, NbinTot, MinTot, MaxTot );
  for(int i=0; i<NumOfSegSFT_UV; i++){
    HB1( SFTUHid +1000+i+1, Form("SFTU Tdc -%d", i+1), NbinTdc, MinTdc, MaxTdc );
    HB1( SFTUHid +3000+i+1, Form("SFTU Tot -%d", i+1), NbinTot, MinTot, MaxTot );
    HB2( SFTUHid +5000+i+1, Form("SFTU Time/Tot -%d", i+1),
	 NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
    HB2( SFTUHid +7000+i+1, Form("SFTU CTime/Tot -%d", i+1),
	 NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  }
  HB1( SFTUHid +21, "SFTU Time",        NbinTime, MinTime, MaxTime );
  HB2( SFTUHid +23, "SFTU Time/Tot",    NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  HB1( SFTUHid +31, "SFTU CTime",       NbinTime, MinTime, MaxTime );
  HB2( SFTUHid +33, "SFTU CTime/Tot",   NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );

  HB1( SFTUHid +101, "SFTU NCluster", 100, 0, 100);
  HB1( SFTUHid +102, "SFTU Cluster Size", 5, 0, 5);
  HB1( SFTUHid +103, "SFTU CTime (Cluster)", NbinTime, MinTime, MaxTime );
  HB1( SFTUHid +104, "SFTU Tot (Cluster)", NbinTot, MinTot, MaxTot );
  HB2( SFTUHid +105, "SFTU CTime%Tot (Cluster)",
       NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  HB1( SFTUHid +106, "SFTU Cluster Position",
       NumOfSegSFT_UV, -0.5*(double)NumOfSegSFT_UV, 0.5*(double)NumOfSegSFT_UV);

  //SFT-X
  HB1( SFTXHid + 1, "SFTX Nhits U",   NumOfSegSFT_X, 0., (double)NumOfSegSFT_X );
  HB1( SFTXHid + 2, "SFTX Nhits D",   NumOfSegSFT_X, 0., (double)NumOfSegSFT_X );
  HB1( SFTXHid + 3, "SFTX Nhits",     NumOfSegSFT_X, 0., (double)NumOfSegSFT_X );
  HB1( SFTXHid + 4, "SFTX Hitpat U",  NumOfSegSFT_X, 0., (double)NumOfSegSFT_X );
  HB1( SFTXHid + 5, "SFTX Hitpat D",  NumOfSegSFT_X, 0., (double)NumOfSegSFT_X );
  HB1( SFTXHid + 6, "SFTX Tdc U",      NbinTdc, MinTdc, MaxTdc );
  HB1( SFTXHid + 7, "SFTX Tdc D",      NbinTdc, MinTdc, MaxTdc );
  HB1( SFTXHid + 8, "SFTX Tot U",      NbinTot, MinTot, MaxTot );
  HB1( SFTXHid + 9, "SFTX Tot D",      NbinTot, MinTot, MaxTot );
  HB2( SFTXHid +10, "SFTX Tdc U%Seg",
       NumOfSegSFT_X, 0., (double)NumOfSegSFT_X, NbinTdc, MinTdc, MaxTdc );
  HB2( SFTXHid +11, "SFTX Tdc D%Seg",
       NumOfSegSFT_X, 0., (double)NumOfSegSFT_X, NbinTdc, MinTdc, MaxTdc );
  HB2( SFTXHid +12, "SFTX Tot U%Seg",
       NumOfSegSFT_X, 0., (double)NumOfSegSFT_X, NbinTot, MinTot, MaxTot );
  HB2( SFTXHid +13, "SFTX Tot D%Seg",
       NumOfSegSFT_X, 0., (double)NumOfSegSFT_X, NbinTot, MinTot, MaxTot );
  for(int i=0; i<NumOfSegSFT_X; i++){
    HB1( SFTXHid +1000+i+1, Form("SFTX Tdc U-%d", i+1), NbinTdc, MinTdc, MaxTdc );
    HB1( SFTXHid +2000+i+1, Form("SFTX Tdc D-%d", i+1), NbinTdc, MinTdc, MaxTdc );
    HB1( SFTXHid +3000+i+1, Form("SFTX Tot U-%d", i+1), NbinTot, MinTot, MaxTot );
    HB1( SFTXHid +4000+i+1, Form("SFTX Tot D-%d", i+1), NbinTot, MinTot, MaxTot );
    HB2( SFTXHid +5000+i+1, Form("SFTX Time/Tot U-%d", i+1),
	 NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
    HB2( SFTXHid +6000+i+1, Form("SFTX Time/Tot D-%d", i+1),
	 NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
    HB2( SFTXHid +7000+i+1, Form("SFTX CTime/Tot U-%d", i+1),
	 NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
    HB2( SFTXHid +8000+i+1, Form("SFTX CTime/Tot D-%d", i+1),
	 NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  }
  HB1( SFTXHid +21, "SFTX Time U",     NbinTime, MinTime, MaxTime );
  HB1( SFTXHid +22, "SFTX Time D",     NbinTime, MinTime, MaxTime );
  HB2( SFTXHid +23, "SFTX Time/Tot U", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  HB2( SFTXHid +24, "SFTX Time/Tot D", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  HB1( SFTXHid +31, "SFTX CTime U",     NbinTime, MinTime, MaxTime );
  HB1( SFTXHid +32, "SFTX CTime D",     NbinTime, MinTime, MaxTime );
  HB2( SFTXHid +33, "SFTX CTime/Tot U", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  HB2( SFTXHid +34, "SFTX CTime/Tot D", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );

  HB1( SFTXHid +101, "SFTX NCluster", 100, 0, 100);
  HB1( SFTXHid +102, "SFTX Cluster Size", 5, 0, 5);
  HB1( SFTXHid +103, "SFTX CTime (Cluster)", NbinTime, MinTime, MaxTime );
  HB1( SFTXHid +104, "SFTX Tot (Cluster)", NbinTot, MinTot, MaxTot );
  HB2( SFTXHid +105, "SFTX CTime%Tot (Cluster)",
       NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime );
  HB1( SFTXHid +106, "SFTX Cluster Position",
       NumOfSegSFT_X, -0.5*(double)NumOfSegSFT_X, 0.5*(double)NumOfSegSFT_X);

  //Tree
  HBTree( "ea0c", "tree of Easiroc" );
  //Trig
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trignhits", &event.trignhits, "trignhits/I");
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  //BFT
  tree->Branch("bft_nhits",     &event.bft_nhits,        "bft_nhits/I");
  tree->Branch("bft_unhits",    &event.bft_unhits,       "bft_unhits/I");
  tree->Branch("bft_dnhits",    &event.bft_dnhits,       "bft_dnhits/I");
  tree->Branch("bft_uhitpat",    event.bft_uhitpat,      "bft_uhitpat[bft_unhits]/I");
  tree->Branch("bft_dhitpat",    event.bft_dhitpat,      "bft_dhitpat[bft_dnhits]/I");
  tree->Branch("bft_utdc",       event.bft_utdc,         Form("bft_utdc[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_dtdc",       event.bft_dtdc,         Form("bft_dtdc[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_utrailing",  event.bft_utrailing,    Form("bft_utrailing[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_dtrailing",  event.bft_dtrailing,    Form("bft_dtrailing[%d][%d]/D",
							      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_utot",       event.bft_utot,         Form("bft_utot[%d][%d]/D",
						      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_dtot",       event.bft_dtot,         Form("bft_dtot[%d][%d]/D",
						      NumOfSegBFT, MaxDepth));
  tree->Branch("bft_udepth",     event.bft_udepth,       Form("bft_udepth[%d]/I", NumOfSegBFT));
  tree->Branch("bft_ddepth",     event.bft_ddepth,       Form("bft_ddepth[%d]/I", NumOfSegBFT));
  tree->Branch("bft_ncl",       &event.bft_ncl,          "bft_ncl/I");
  tree->Branch("bft_clsize",     event.bft_clsize,       "bft_clsize[bft_ncl]/I");
  tree->Branch("bft_ctime",      event.bft_ctime,        "bft_ctime[bft_ncl]/D");
  tree->Branch("bft_ctot",       event.bft_ctot,         "bft_ctot[bft_ncl]/D");
  tree->Branch("bft_clpos",      event.bft_clpos,        "bft_clpos[bft_ncl]/D");

  //SCH
  tree->Branch("sch_nhits",     &event.sch_nhits,        "sch_nhits/I");
  tree->Branch("sch_hitpat",     event.sch_hitpat,       "sch_hitpat[sch_nhits]/I");
  tree->Branch("sch_tdc",        event.sch_tdc,          Form("sch_tdc[%d][%d]/D",
							      NumOfSegSCH, MaxDepth));
  tree->Branch("sch_trailing",   event.sch_trailing,     Form("sch_trailing[%d][%d]/D",
							      NumOfSegSCH, MaxDepth));
  tree->Branch("sch_tot",        event.sch_tot,          Form("sch_tot[%d][%d]/D",
							      NumOfSegSCH, MaxDepth));
  tree->Branch("sch_depth",      event.sch_depth,        Form("sch_depth[%d]/I", NumOfSegSCH));
  tree->Branch("sch_ncl",       &event.sch_ncl,          "sch_ncl/I");
  tree->Branch("sch_clsize",     event.sch_clsize,       "sch_clsize[sch_ncl]/I");
  tree->Branch("sch_ctime",      event.sch_ctime,        "sch_ctime[sch_ncl]/D");
  tree->Branch("sch_ctot",       event.sch_ctot,         "sch_ctot[sch_ncl]/D");
  tree->Branch("sch_clpos",      event.sch_clpos,        "sch_clpos[sch_ncl]/D");

  //SFT-V
  tree->Branch("sftv_nhits",     &event.sftv_nhits,        "sftv_nhits/I");
  tree->Branch("sftv_hitpat",     event.sftv_hitpat,       "sftv_hitpat[sftv_nhits]/I");
  tree->Branch("sftv_tdc",        event.sftv_tdc,          Form("sftv_tdc[%d][%d]/D",
							      NumOfSegSFT_UV, MaxDepth));
  tree->Branch("sftv_trailing",   event.sftv_trailing,     Form("sftv_trailing[%d][%d]/D",
							      NumOfSegSFT_UV, MaxDepth));
  tree->Branch("sftv_tot",        event.sftv_tot,          Form("sftv_tot[%d][%d]/D",
						      NumOfSegSFT_UV, MaxDepth));
  tree->Branch("sftv_depth",      event.sftv_depth,        Form("sftv_depth[%d]/I", NumOfSegSFT_UV));
  tree->Branch("sftv_ncl",       &event.sftv_ncl,          "sftv_ncl/I");
  tree->Branch("sftv_clsize",     event.sftv_clsize,       "sftv_clsize[sftv_ncl]/I");
  tree->Branch("sftv_ctime",      event.sftv_ctime,        "sftv_ctime[sftv_ncl]/D");
  tree->Branch("sftv_ctot",       event.sftv_ctot,         "sftv_ctot[sftv_ncl]/D");
  tree->Branch("sftv_clpos",      event.sftv_clpos,        "sftv_clpos[sftv_ncl]/D");

  //SFT-U
  tree->Branch("sftu_nhits",     &event.sftu_nhits,        "sftu_nhits/I");
  tree->Branch("sftu_hitpat",     event.sftu_hitpat,       "sftu_hitpat[sftu_nhits]/I");
  tree->Branch("sftu_tdc",        event.sftu_tdc,          Form("sftu_tdc[%d][%d]/D",
							      NumOfSegSFT_UV, MaxDepth));
  tree->Branch("sftu_trailing",   event.sftu_trailing,     Form("sftu_trailing[%d][%d]/D",
							      NumOfSegSFT_UV, MaxDepth));
  tree->Branch("sftu_tot",        event.sftu_tot,          Form("sftu_tot[%d][%d]/D",
						      NumOfSegSFT_UV, MaxDepth));
  tree->Branch("sftu_depth",      event.sftu_depth,        Form("sftu_depth[%d]/I", NumOfSegSFT_UV));
  tree->Branch("sftu_ncl",       &event.sftu_ncl,          "sftu_ncl/I");
  tree->Branch("sftu_clsize",     event.sftu_clsize,       "sftu_clsize[sftu_ncl]/I");
  tree->Branch("sftu_ctime",      event.sftu_ctime,        "sftu_ctime[sftu_ncl]/D");
  tree->Branch("sftu_ctot",       event.sftu_ctot,         "sftu_ctot[sftu_ncl]/D");
  tree->Branch("sftu_clpos",      event.sftu_clpos,        "sftu_clpos[sftu_ncl]/D");

  //SFT-X
  tree->Branch("sftx_nhits",     &event.sftx_nhits,        "sftx_nhits/I");
  tree->Branch("sftx_unhits",    &event.sftx_unhits,       "sftx_unhits/I");
  tree->Branch("sftx_dnhits",    &event.sftx_dnhits,       "sftx_dnhits/I");
  tree->Branch("sftx_uhitpat",    event.sftx_uhitpat,      "sftx_uhitpat[sftx_unhits]/I");
  tree->Branch("sftx_dhitpat",    event.sftx_dhitpat,      "sftx_dhitpat[sftx_dnhits]/I");
  tree->Branch("sftx_utdc",       event.sftx_utdc,         Form("sftx_utdc[%d][%d]/D",
							      NumOfSegSFT_X, MaxDepth));
  tree->Branch("sftx_dtdc",       event.sftx_dtdc,         Form("sftx_dtdc[%d][%d]/D",
							      NumOfSegSFT_X, MaxDepth));
  tree->Branch("sftx_utrailing",  event.sftx_utrailing,    Form("sftx_utrailing[%d][%d]/D",
							      NumOfSegSFT_X, MaxDepth));
  tree->Branch("sftx_dtrailing",  event.sftx_dtrailing,    Form("sftx_dtrailing[%d][%d]/D",
							      NumOfSegSFT_X, MaxDepth));
  tree->Branch("sftx_utot",       event.sftx_utot,         Form("sftx_utot[%d][%d]/D",
						      NumOfSegSFT_X, MaxDepth));
  tree->Branch("sftx_dtot",       event.sftx_dtot,         Form("sftx_dtot[%d][%d]/D",
						      NumOfSegSFT_X, MaxDepth));
  tree->Branch("sftx_udepth",     event.sftx_udepth,       Form("sftx_udepth[%d]/I", NumOfSegSFT_X));
  tree->Branch("sftx_ddepth",     event.sftx_ddepth,       Form("sftx_ddepth[%d]/I", NumOfSegSFT_X));
  tree->Branch("sftx_ncl",       &event.sftx_ncl,          "sftx_ncl/I");
  tree->Branch("sftx_clsize",     event.sftx_clsize,       "sftx_clsize[sftx_ncl]/I");
  tree->Branch("sftx_ctime",      event.sftx_ctime,        "sftx_ctime[sftx_ncl]/D");
  tree->Branch("sftx_ctot",       event.sftx_ctot,         "sftx_ctot[sftx_ncl]/D");
  tree->Branch("sftx_clpos",      event.sftx_clpos,        "sftx_clpos[sftx_ncl]/D");


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
      InitializeParameter<UserParamMan>("USER")  );
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
