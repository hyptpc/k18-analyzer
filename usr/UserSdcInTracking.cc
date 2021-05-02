// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <UnpackerManager.hh>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "RMAnalyzer.hh"
#include "KuramaLib.hh"
#include "MathTools.hh"
#include "RawData.hh"
#include "RootHelper.hh"

#define HodoCut     0
#define TotCut      1
#define Chi2Cut     0
#define MaxMultiCut 0
#define BcOutCut    0

namespace
{
  using namespace root;
  auto&gRM = RMAnalyzer::GetInstance();
  const auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
VEvent::VEvent()
{
}

//_____________________________________________________________________________
VEvent::~VEvent()
{
}

//_____________________________________________________________________________
class EventSdcInTracking : public VEvent
{
private:
  RawData*      rawData;
  DCAnalyzer*   DCAna;
  HodoAnalyzer* hodoAna;

public:
  static TString ClassName();
        EventSdcInTracking();
       ~EventSdcInTracking();
  bool  ProcessingBegin();
  bool  ProcessingEnd();
  bool  ProcessingNormal();
  bool  InitializeHistograms();
  void  InitializeEvent();
};

//_____________________________________________________________________________
TString
EventSdcInTracking::ClassName()
{
  static TString s_name("EventSdcInTracking");
  return s_name;
}

//_____________________________________________________________________________
EventSdcInTracking::EventSdcInTracking()
  : VEvent(),
    rawData(new RawData),
    DCAna(new DCAnalyzer),
    hodoAna(new HodoAnalyzer)
{
}

//_____________________________________________________________________________
EventSdcInTracking::~EventSdcInTracking()
{
  if(DCAna) delete DCAna;
  if(hodoAna) delete hodoAna;
  if(rawData) delete rawData;
}

//_____________________________________________________________________________
bool
EventSdcInTracking::ProcessingBegin()
{
  InitializeEvent();
  return true;
}

//_____________________________________________________________________________
struct Event
{
  int evnum;

  int trignhits;
  int trigpat[MaxHits];
  int trigflag[NumOfSegTrig];

  int nhBh1;
  double tBh1[MaxHits];
  double deBh1[MaxHits];

  int nhBh2;
  double tBh2[MaxHits];
  double deBh2[MaxHits];
  double Bh2Seg[MaxHits];

  double Time0Seg;
  double deTime0;
  double Time0;
  double CTime0;

  double btof;

  int ntBcOut;

  int nhit[NumOfLayersSdcIn];
  int nlayer;
  double pos[NumOfLayersSdcIn][MaxHits];

  int ntrack;
  double chisqr[MaxHits];
  double x0[MaxHits];
  double y0[MaxHits];
  double u0[MaxHits];
  double v0[MaxHits];
};

//_____________________________________________________________________________
namespace root
{
  Event  event;
  TH1   *h[MaxHist];
  TTree *tree;
}

//_____________________________________________________________________________
bool
EventSdcInTracking::ProcessingNormal()
{

#if HodoCut
  static const double MinDeBH2   = gUser.GetParameter("DeBH2", 0);
  static const double MaxDeBH2   = gUser.GetParameter("DeBH2", 1);
  static const double MinDeBH1   = gUser.GetParameter("DeBH1", 0);
  static const double MaxDeBH1   = gUser.GetParameter("DeBH1", 1);
  static const double MinBeamToF = gUser.GetParameter("BTOF",  0);
  static const double MaxBeamToF = gUser.GetParameter("BTOF",  1);
#endif
#if TotCut
  static const double MinTotSDC1 = gUser.GetParameter("MinTotSDC1", 0);
  static const double MinTotSDC2 = gUser.GetParameter("MinTotSDC2", 0);
#endif
#if MaxMultiCut
  static const double MaxMultiHitBcOut = gUser.GetParameter("MaxMultiHitBcOut");
  static const double MaxMultiHitSdcIn = gUser.GetParameter("MaxMultiHitSdcIn");
#endif

  rawData->DecodeHits();

  gRM.Decode();

  event.evnum = gRM.EventNumber();

  // Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  {
    const HodoRHitContainer &cont=rawData->GetTrigRawHC();
    int trignhits = 0;
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg = hit->SegmentId()+1;
      int tdc = hit->GetTdc1();
      if( tdc>0 ){
	trigger_flag.set( seg-1 );
	event.trigpat[trignhits++] = seg;
  	event.trigflag[seg-1]      = tdc;
      }
    }
    event.trignhits = trignhits;
  }

  HF1( 1, 0. );

  if( trigger_flag[trigger::kSpillEnd] ||
      trigger_flag[trigger::kL1SpillOff] )
    return true;

  HF1( 1, 1. );


 //////////////BH2 time 0
  hodoAna->DecodeBH2Hits( rawData );
  int nhBh2 = hodoAna->GetNHitsBH2();
  event.nhBh2 = nhBh2;
#if HodoCut
  if( nhBh2==0 ) return true;
#endif
  HF1( 1, 2 );

  double time0 = -9999.;
  //////////////BH2 Analysis
  for( int i=0; i<nhBh2; ++i ){
    BH2Hit* hit = hodoAna->GetHitBH2( i );
    if( !hit ) continue;
    double seg = hit->SegmentId()+1;
    double cmt = hit->CMeanTime();
    double dE  = hit->DeltaE();

#if HodoCut
    if( dE<MinDeBH2 || MaxDeBH2<dE ) continue;
#endif
    event.tBh2[i]   = cmt;
    event.deBh2[i]  = dE;
    event.Bh2Seg[i] = seg;
  }

  BH2Cluster *cl_time0 = hodoAna->GetTime0BH2Cluster();
  if( cl_time0 ){
    event.Time0Seg = cl_time0->MeanSeg()+1;
    event.deTime0  = cl_time0->DeltaE();
    event.Time0    = cl_time0->Time0();
    event.CTime0   = cl_time0->CTime0();
    time0          = cl_time0->CTime0();
  } else {
#if HodoCut
    return true;
#endif
  }

  HF1( 1, 3. );

  //////////////BH1 Analysis
  hodoAna->DecodeBH1Hits( rawData );
  int nhBh1 = hodoAna->GetNHitsBH1();
  event.nhBh1 = nhBh1;
#if HodoCut
  if( nhBh1==0 ) return true;
#endif
  HF1( 1, 4 );

  for( int i=0; i<nhBh1; ++i ){
    Hodo2Hit *hit = hodoAna->GetHitBH1( i );
    if( !hit ) continue;
    double cmt = hit->CMeanTime();
    double dE  = hit->DeltaE();
#if HodoCut
    if( dE<MinDeBH1 || MaxDeBH1<dE ) continue;
    if( btof<MinBeamToF || MaxBeamToF<btof ) continue;
#endif
    event.tBh1[i]  = cmt;
    event.deBh1[i] = dE;
  }

  double btof0 = -999.;
  HodoCluster* cl_btof0 = event.Time0Seg > 0 ?
    hodoAna->GetBtof0BH1Cluster( event.CTime0 ) : nullptr;
  if( cl_btof0 ) btof0 = cl_btof0->CMeanTime() - time0;
  event.btof = btof0;

  HF1( 1, 5. );


  //////////////BCout
  DCAna->DecodeRawHits( rawData );
#if TotCut
  DCAna->TotCutSDC1( MinTotSDC1 );
  DCAna->TotCutSDC2( MinTotSDC2 );
#endif
  //BC3&BC4
  double multi_BcOut=0.;
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetBcOutHC(layer);
      int nhOut=contOut.size();
      multi_BcOut += double(nhOut);
    }
  }
#if MaxMultiCut
  if( multi_BcOut/double(NumOfLayersBcOut) > MaxMultiHitBcOut ){
    return true;
  }
#endif

#if 0
  {
    //  std::cout << "==========TrackSearch BcOut============" << std::endl;
    int ntOk = 0;
    DCAna->TrackSearchBcOut();
    int ntBcOut = DCAna->GetNtracksBcOut();
    if( MaxHits<ntBcOut ){
      std::cout << "#W " << FUNC_NAME << " "
		<< "too many ntBcOut " << ntBcOut << "/" << MaxHits << std::endl;
      ntBcOut = MaxHits;
    }
    event.ntBcOut = ntBcOut;
    for( int it=0; it<ntBcOut; ++it ){
      DCLocalTrack *tp=DCAna->GetTrackBcOut(it);
      double chisqr=tp->GetChiSquare();
      if( chisqr<20. ) ntOk++;
    }
#if BcOutCut
    if( ntOk==0 ) return true;
#endif
  }
#endif


  //////////////SdcIn number of hit layer
  HF1( 1, 10. );
  double multi_SdcIn=0.;
  {
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetSdcInHC(layer);
      int nhIn = contIn.size();
      event.nhit[layer-1] = nhIn;
      if( nhIn>0 ) event.nlayer++;
      multi_SdcIn += double(nhIn);
      HF1( 100*layer, nhIn );
      int plane_eff = (layer-1)*3;
      bool fl_valid_sig = false;

      for( int i=0; i<nhIn; ++i ){
	DCHit *hit=contIn[i];
	double wire=hit->GetWire();
	HF1( 100*layer+1, wire-0.5 );
	int nhtdc = hit->GetTdcSize();
	int tdc1st = -1;
	for( int k=0; k<nhtdc; k++ ){
	  int tdc = hit->GetTdcVal(k);
	  HF1( 100*layer+2, tdc );
	  HF1( 10000*layer+int(wire), tdc );
	  //	  HF2( 1000*layer, tdc, wire-0.5 );
	  if( tdc > tdc1st ){
	    tdc1st = tdc;
	    fl_valid_sig = true;
	  }
	}
	HF1( 100*layer+6, tdc1st );

	if( i<MaxHits )
	  event.pos[layer-1][i] = hit->GetWirePosition();

	int nhdt = hit->GetDriftTimeSize();
	int tot1st = -1;
	for( int k=0; k<nhdt; k++ ){
	  double dt = hit->GetDriftTime(k);
	  HF1( 100*layer+3, dt );
	  HF1( 10000*layer+1000+int(wire), dt );

	  double tot = hit->GetTot(k);
	  HF1( 100*layer+5, tot );
	  if( tot > tot1st ){
	    tot1st = tot;
	  }
	}
	HF1( 100*layer+7, tot1st );
	int nhdl = hit->GetDriftTimeSize();
	for( int k=0; k<nhdl; k++ ){
	  double dl = hit->GetDriftLength(k);
	  HF1( 100*layer+4, dl );
	}
      }
      if( fl_valid_sig ) ++plane_eff;
      HF1(38, plane_eff);
    }
  }

#if MaxMultiCut
  if( multi_SdcIn/double(NumOfLayersSdcIn) > MaxMultiHitSdcIn )
   return true;
#endif

  HF1( 1, 11. );

  // std::cout << "==========TrackSearch SdcIn============" << std::endl;
#if 1
 #if Chi2Cut
  DCAna->ChiSqrCutSdcIn(30.);
 #endif

  DCAna->TrackSearchSdcIn();
  int nt = DCAna->GetNtracksSdcIn();
  if( MaxHits<nt ){
    std::cout << "#W " << FUNC_NAME << " "
	      << "too many ntSdcIn " << nt << "/" << MaxHits << std::endl;
    nt = MaxHits;
  }
  event.ntrack = nt;
  HF1( 10, double(nt) );
  for( int it=0; it<nt; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackSdcIn(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double x0=tp->GetX0(), y0=tp->GetY0();
    double u0=tp->GetU0(), v0=tp->GetV0();
    double cost = 1./std::sqrt(1.+u0*u0+v0*v0);
    double theta = std::acos(cost)*math::Rad2Deg();
    event.chisqr[it]=chisqr;
    event.x0[it]=x0;
    event.y0[it]=y0;
    event.u0[it]=u0;
    event.v0[it]=v0;

    HF1( 11, double(nh) );
    HF1( 12, chisqr );
    HF1( 14, x0 ); HF1( 15, y0 );
    HF1( 16, u0 ); HF1( 17, v0 );
    HF2( 18, x0, u0 ); HF2( 19, y0, v0 );
    HF2( 20, x0, y0 );

    double xtgt=tp->GetX(-1289.1), ytgt=tp->GetY(-1289.1);
    double utgt=u0, vtgt=v0;
    HF1( 21, xtgt ); HF1( 22, ytgt );
    HF1( 23, utgt ); HF1( 24, vtgt );
    HF2( 25, xtgt, utgt ); HF2( 26, ytgt, vtgt );
    HF2( 27, xtgt, ytgt );

    double xbac=tp->GetX(-815.0), ybac=tp->GetY(-815.0);
    double ubac=u0, vbac=v0;
    HF1( 31, xbac ); HF1( 32, ybac );
    HF1( 33, ubac ); HF1( 34, vbac );
    HF2( 35, xbac, ubac ); HF2( 36, ybac, vbac );
    HF2( 37, xbac, ybac );

    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      if( !hit ) continue;
      int layerId = hit->GetLayer();
      HF1( 13, layerId );

      double wire=hit->GetWire();
      double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1( 100*layerId+11, wire-0.5 );
      HF1( 100*layerId+12, dt );
      HF1( 100*layerId+13, dl );
      HF1( 10000*layerId+ 5000 +int(wire), dt);
      double xcal=hit->GetXcal(), ycal=hit->GetYcal();
      double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1( 100*layerId+14, pos );
      HF1( 100*layerId+15, res );
      HF2( 100*layerId+16, pos, res );
      HF2( 100*layerId+17, xcal, ycal);
      //      HF1( 100000*layerId+50000+wire, res);
      double wp=hit->GetWirePosition();
      double sign=1.;
      if( pos-wp<0. ) sign=-1;
      HF2( 100*layerId+18, sign*dl, res );
      double xlcal=hit->GetLocalCalPos();
      HF2( 100*layerId+19, dt, xlcal-wp);

      double tot = hit->GetTot();
      HF1( 100*layerId+40, tot);

      if (theta>=0 && theta<15)
	HF1( 100*layerId+71, res );
      else if (theta>=15 && theta<30)
	HF1( 100*layerId+72, res );
      else if (theta>=30 && theta<45)
	HF1( 100*layerId+73, res );
      else if (theta>=45)
	HF1( 100*layerId+74, res );

      if (std::abs(dl-std::abs(xlcal-wp))<2.0) {
	HFProf( 100*layerId+20, dt, std::abs(xlcal-wp));
	HF2( 100*layerId+22, dt, std::abs(xlcal-wp));
	HFProf( 100000*layerId+3000+int(wire), xlcal-wp,dt);
	HF2( 100000*layerId+4000+int(wire), xlcal-wp,dt);
      }
    }
  }

  HF1( 1, 12. );
#endif

  return true;
}

//_____________________________________________________________________________
void
EventSdcInTracking::InitializeEvent()
{
  event.evnum     =  0;
  event.trignhits =  0;
  event.nlayer    =  0;
  event.ntrack    =  0;
  event.nhBh2     =  0;
  event.nhBh1     =  0;

  event.btof      = -999.;

  event.ntBcOut   =  0;

  event.Time0Seg  = -1;
  event.deTime0   = -1;
  event.Time0     = -999;
  event.CTime0    = -999;

  for( int it=0; it<MaxHits; it++){
    event.tBh1[it]   = -9999.;
    event.deBh1[it]  = -9999.;

    event.Bh2Seg[it] = -1;
    event.tBh2[it]   = -9999.;
    event.deBh2[it]  = -9999.;

    event.chisqr[it] = -1.0;
    event.x0[it]     = -9999.0;
    event.y0[it]     = -9999.0;
    event.u0[it]     = -9999.0;
    event.v0[it]     = -9999.0;

    event.trigpat[it] = -1;
  }

  for( int it=0; it<NumOfSegTrig; it++){
    event.trigflag[it] = -1;
  }

  for(int it = 0; it<NumOfLayersSdcIn; it++){
    event.nhit[it] = -1;
    for( int that=0; that<MaxHits; that++ ){
      event.pos[it][that] = -9999.;
    }
  }

}

//_____________________________________________________________________________
bool
EventSdcInTracking::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new EventSdcInTracking;
}

//_____________________________________________________________________________
namespace
{
  const int    NbinSdcInTdc = 2000;
  const double MinSdcInTdc  =    0.;
  const double MaxSdcInTdc  = 2000.;

  const int    NbinSDC1DT = 200   ;
  const double MinSDC1DT  = -30.  ;
  const double MaxSDC1DT  = 150.  ;
  const int    NbinSDC1DL =  45   ;
  const double MinSDC1DL  =  -0.5 ;
  const double MaxSDC1DL  =   4.0 ;

  const int    NbinSDC2DT = 300 ;
  const double MinSDC2DT  = -50.;
  const double MaxSDC2DT  = 200.;
  const int    NbinSDC2DL = 100 ;
  const double MinSDC2DL  =  -3.;
  const double MaxSDC2DL  =   7.;

  const int NbinRes   =  200;
  const double MinRes = -1.;
  const double MaxRes =  1.;
}

//_____________________________________________________________________________
bool
ConfMan:: InitializeHistograms()
{
  HB1( 1, "Status", 20, 0., 20. );

  // SdcInTracking
  for( int i=1; i<=NumOfLayersSdcIn; ++i ){

    std::string tag;
    int nwire = 0, nbindt = 1, nbindl = 1;
    double mindt = 0., maxdt = 1., mindl = 0., maxdl = 1.;
    if(i<=NumOfLayersSDC1){
      tag    = "SDC1";
      nwire  = MaxWireSDC1;
      nbindt = NbinSDC1DT;
      mindt  = MinSDC1DT;
      maxdt  = MaxSDC1DT;
      nbindl = NbinSDC1DL;
      mindl  = MinSDC1DL;
      maxdl  = MaxSDC1DL;
    }else if(i<=NumOfLayersSdcIn){
      tag = "SDC2";
      nwire   = ( i==7 || i==8 ) ? MaxWireSDC2X : MaxWireSDC2Y;
      nbindt = NbinSDC2DT;
      mindt  = MinSDC2DT;
      maxdt  = MaxSDC2DT;
      nbindl = NbinSDC2DL;
      mindl  = MinSDC2DL;
      maxdl  = MaxSDC2DL;
    }

    TString title0 = Form("#Hits %s#%2d", tag.c_str(), i);
    TString title1 = Form("Hitpat %s#%2d", tag.c_str(), i);
    TString title2 = Form("Tdc %s#%2d", tag.c_str(), i);
    TString title3 = Form("Drift Time %s#%2d", tag.c_str(), i);
    TString title4 = Form("Drift Length %s#%2d", tag.c_str(), i);
    TString title5 = Form("TOT %s#%2d", tag.c_str(), i);
    TString title6 = Form("Tdc 1st %s#%2d", tag.c_str(), i);
    TString title7 = Form("TOT 1st %s#%2d", tag.c_str(), i);
    HB1( 100*i+0, title0, nwire+1, 0., double(nwire+1) );
    HB1( 100*i+1, title1, nwire, 0., double(nwire) );
    HB1( 100*i+2, title2, NbinSdcInTdc, MinSdcInTdc, MaxSdcInTdc );
    HB1( 100*i+3, title3, nbindt, mindt, maxdt );
    HB1( 100*i+4, title4, nbindl, mindl, maxdl );
    HB1( 100*i+5, title5, 500,  0, 500 );
    HB1( 100*i+6, title6, NbinSdcInTdc, MinSdcInTdc, MaxSdcInTdc );
    HB1( 100*i+7, title7, 500,  0, 500 );
    for ( int wire=1; wire<=nwire; wire++ ){
      TString title11 = Form("Tdc %s#%2d  Wire#%4d", tag.c_str(), i, wire);
      TString title12 = Form("DriftTime %s#%2d Wire#%4d", tag.c_str(), i, wire);
      TString title13 = Form("DriftLength %s#%2d Wire#%d", tag.c_str(), i, wire);
      TString title14 = Form("DriftTime %s#%2d Wire#%4d [Track]", tag.c_str(), i, wire);
      HB1( 10000*i+wire, title11, NbinSdcInTdc, MinSdcInTdc, MaxSdcInTdc );
      HB1( 10000*i+1000+wire, title12, nbindt, mindt, maxdt );
      HB1( 10000*i+2000+wire, title13, nbindl, mindl, maxdl );
      HB1( 10000*i+5000+wire, title14, nbindt, mindt, maxdt );
    }


    // Tracking Histgrams
    TString title11 = Form("HitPat SdcIn%2d [Track]", i);
    TString title12 = Form("DriftTime SdcIn%2d [Track]", i);
    TString title13 = Form("DriftLength SdcIn%2d [Track]", i);
    TString title14 = Form("Position SdcIn%2d", i);
    TString title15 = Form("Residual SdcIn%2d", i);
    TString title16 = Form("Resid%%Pos SdcIn%2d", i);
    TString title17 = Form("Y%%Xcal SdcIn%2d", i);
    TString title18 = Form("Res%%DL SdcIn%2d", i);
    TString title19 = Form("HitPos%%DriftTime SdcIn%2d", i);
    TString title20 = Form("DriftLength%%DriftTime SdcIn%2d", i);
    TString title21 = title15 + " [w/o Self]";
    TString title22 = title20;
    TString title40 = Form("TOT SdcIn%2d [Track]", i);
    TString title71 = Form("Residual SdcIn%2d (0<theta<15)", i);
    TString title72 = Form("Residual SdcIn%2d (15<theta<30)", i);
    TString title73 = Form("Residual SdcIn%2d (30<theta<45)", i);
    TString title74 = Form("Residual SdcIn%2d (45<theta)", i);
    HB1( 100*i+11, title11, nwire, 0., double(nwire) );
    HB1( 100*i+12, title12, nbindt, mindt, maxdt );
    HB1( 100*i+13, title13, nbindl, mindl, maxdl );
    HB1( 100*i+14, title14, 100, -250., 250. );
    HB1( 100*i+15, title15, NbinRes, MinRes, MaxRes );
    HB2( 100*i+16, title16, 250, -250., 250., NbinRes, MinRes, MaxRes );
    HB2( 100*i+17, title17, 100, -250., 250., 100, -250., 250. );
    HB2( 100*i+18, title18, 100, -3., 3., NbinRes, MinRes, MaxRes );
    HB2( 100*i+19, title19, nbindt, mindt, maxdt, 100, -4., 4. );
    HBProf( 100*i+20, title20, nbindt, mindt, maxdt, mindl, maxdl );
    HB2( 100*i+22, title22, nbindt, mindt, maxdt, nbindl, mindl, maxdl );
    HB1( 100*i+21, title21, 200, -5.0, 5.0 );
    HB1( 100*i+40, title40, 360,    0, 300 );
    HB1( 100*i+71, title71, 200, -5.0, 5.0 );
    HB1( 100*i+72, title72, 200, -5.0, 5.0 );
    HB1( 100*i+73, title73, 200, -5.0, 5.0 );
    HB1( 100*i+74, title74, 200, -5.0, 5.0 );
    // HB2( 1000*i, Form("Wire%%Tdc for LayerId = %d", i),
    // 	 NbinSdcOutTdc/4, MinSdcOutTdc, MaxSdcOutTdc,
    // 	 MaxWire+1, 0., double(MaxWire+1) );

    for (int j=1; j<=nwire; j++) {
      TString title = Form("XT of Layer %2d Wire #%4d", i, j);
      HBProf( 100000*i+3000+j, title, 100, -4., 4., -5, 50 );
      HB2( 100000*i+4000+j, title, 100, -4., 4., 40, -5., 50. );
    }
  }

  // Tracking Histgrams
  HB1( 10, "#Tracks SdcIn", 10, 0., 10. );
  HB1( 11, "#Hits of Track SdcIn", 15, 0., 15. );
  HB1( 12, "Chisqr SdcIn", 500, 0., 50. );
  HB1( 13, "LayerId SdcIn", 15, 0., 15. );
  HB1( 14, "X0 SdcIn", 400, -100., 100. );
  HB1( 15, "Y0 SdcIn", 400, -100., 100. );
  HB1( 16, "U0 SdcIn", 200, -0.20, 0.20 );
  HB1( 17, "V0 SdcIn", 200, -0.20, 0.20 );
  HB2( 18, "U0%X0 SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 19, "V0%Y0 SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 20, "X0%Y0 SdcIn", 100, -100., 100., 100, -100, 100 );

  HB1( 21, "Xtgt SdcIn", 400, -100., 100. );
  HB1( 22, "Ytgt SdcIn", 400, -100., 100. );
  HB1( 23, "Utgt SdcIn", 200, -0.20, 0.20 );
  HB1( 24, "Vtgt SdcIn", 200, -0.20, 0.20 );
  HB2( 25, "Utgt%Xtgt SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 26, "Vtgt%Ytgt SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 27, "Xtgt%Ytgt SdcIn", 100, -100., 100., 100, -100, 100 );

  HB1( 31, "Xbac SdcIn", 400, -100., 100. );
  HB1( 32, "Ybac SdcIn", 400, -100., 100. );
  HB1( 33, "Ubac SdcIn", 200, -0.20, 0.20 );
  HB1( 34, "Vbac SdcIn", 200, -0.20, 0.20 );
  HB2( 35, "Ubac%Xbac SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 36, "Vbac%Ybac SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 37, "Xbac%Ybac SdcIn", 100, -100., 100., 100, -100, 100 );

  // Plane eff
  HB1( 38, "Plane Eff", 30, 0, 30);


  ////////////////////////////////////////////
  //Tree
  HBTree("sdcin","tree of SdcInTracking");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trignhits", &event.trignhits, "trignhits/I");
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", MaxHits));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  //Hodoscope
  tree->Branch("nhBh1",    &event.nhBh1,   "nhBh1/I");
  tree->Branch("tBh1",      event.tBh1,    Form("tBh1[%d]/D",   MaxHits));
  tree->Branch("deBh1",     event.deBh1,   Form("deBh1[%d]/D",  MaxHits));

  tree->Branch("nhBh2",    &event.nhBh2,   "nhBh2/I");
  tree->Branch("tBh2",      event.tBh2,    Form("tBh2[%d]/D",   MaxHits));
  tree->Branch("deBh2",     event.deBh2,   Form("deBh2[%d]/D",  MaxHits));
  tree->Branch("Bh2Seg",    event.Bh2Seg,  Form("Bh2Seg[%d]/D", MaxHits));

  tree->Branch("Time0Seg", &event.Time0Seg,  "Time0Seg/D");
  tree->Branch("deTime0",  &event.deTime0,   "deTime0/D");
  tree->Branch("Time0",    &event.Time0,     "Time0/D");
  tree->Branch("CTime0",   &event.CTime0,    "CTime0/D");

  tree->Branch("btof",     &event.btof,     "btof/D");

  tree->Branch("ntBcOut", &event.ntBcOut, "ntBcOut/I");

  tree->Branch("nhit",  &event.nhit,  Form( "nhit[%d]/I", NumOfLayersSdcIn ) );
  tree->Branch("nlayer",   &event.nlayer,   "nlayer/I");
  tree->Branch("ntrack",   &event.ntrack,   "ntrack/I");
  tree->Branch("pos",      &event.pos,     Form("pos[%d][%d]/D",
						 NumOfLayersSdcOut, MaxHits));
  tree->Branch("chisqr",    event.chisqr,   "chisqr[ntrack]/D");
  tree->Branch("x0",        event.x0,       "x0[ntrack]/D");
  tree->Branch("y0",        event.y0,       "y0[ntrack]/D");
  tree->Branch("u0",        event.u0,       "u0[ntrack]/D");
  tree->Branch("v0",        event.v0,       "v0[ntrack]/D");

  // TString layer_name[NumOfLayersSdcIn] =
  //   { "sdc1u1", "sdc1u0", "sdc1x1", "sdc1x0", "sdc1v1", "sdc1v0",
  //     "sdc2x0", "sdc2x1", "sdc2y0", "sdc2y1"};
  // for( int i=0; i<NumOfLayersSdcIn; ++i ){
  //   TString name = Form("%s_pos", layer_name[i].Data() );
  //   TString type = Form("%s_pos[%d]/D", layer_name[i].Data(), MaxHits );
  //   tree->Branch( name, event.pos[i], type );
  // }
  HPrint();
  return true;
}

//_____________________________________________________________________________
bool
ConfMan::InitializeParameterFiles()
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")        &&
      InitializeParameter<DCDriftParamMan>("DCDRFT") &&
      InitializeParameter<DCTdcCalibMan>("DCTDC")    &&
      InitializeParameter<HodoParamMan>("HDPRM")     &&
      InitializeParameter<HodoPHCMan>("HDPHC")       &&
      InitializeParameter<UserParamMan>("USER")      );
}

//_____________________________________________________________________________
bool
ConfMan::FinalizeProcess()
{
  return true;
}
