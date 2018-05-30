/**
 *  file: UserKKAna.cc
 *  date: 2017.04.10
 *
 */

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>

#include "ConfMan.hh"
#include "RMAnalyzer.hh"
#include "FiberCluster.hh"
#include "RootHelper.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "DCRawHit.hh"
#include "KuramaLib.hh"
#include "K18TrackD2U.hh"
#include "MathTools.hh"
#include "MatrixParamMan.hh"
#include "MsTParamMan.hh"
#include "UserParamMan.hh"

#define HodoCut 0

namespace
{
  using namespace root;
  const std::string& class_name("EventKKAna");
  const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
  RMAnalyzer&         gRM   = RMAnalyzer::GetInstance();
  const MsTParamMan&  gMsT  = MsTParamMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const double& zK18Target = gGeom.LocalZ("K18Target");

  const double pB_offset   = 1.000;
  const double pS_offset   = 0.000;
  const double pK18_offset = 0.000;
  const double pKURAMA_offset = 0.000;

  const double x_off = 0.000;
  const double y_off = 0.000;
  const double u_off = 0.000;
  const double v_off = 0.000;

  // TOF-SdcOut
  const double MinXDifTof  =  -40.;
  const double MaxXDifTof  =   70.;
  const double MinYDifTof  =  -230.;
  const double MaxYDifTof  =   200.;

  //Multi cut of DCs
  const double MaxMultiHitBcOut  = 30.;
  const double MaxMultiHitSdcIn  = 30.;
  const double MaxMultiHitSdcOut = 3.0;

  const double MaxChisqrBcOut  = 20.;
  const double MaxChisqrSdcIn  = 20.;
  const double MaxChisqrSdcOut = 20.;
  const double MaxChisqrK18    = 100.;
  const double MaxChisqrKurama    = 1000.;

  const int MaxNTrackBcOut  = 20;
  const int MaxNTrackSdcIn  = 20;
  const int MaxNTrackSdcOut = 20;
  const int MaxNTrackK18    = 4;
  const int MaxNTrackKurama    = 2;

  const double MinTgtXSdcIn  =  -100.;
  const double MaxTgtXSdcIn  =   100.;
  const double MinTgtYSdcIn  =  -50.;
  const double MaxTgtYSdcIn  =   50.;

  const double MinTgtXBcOut  =  -100.;
  const double MaxTgtXBcOut  =   100.;
  const double MinTgtYBcOut  =  -50.;
  const double MaxTgtYBcOut  =   50.;

  const double MinMassSquare = -20.0;
  const double MaxMassSquare =  40.0;

  const double AtomicMassUnit  = 0.93149432;
  const double PionMass        = 0.1395701;
  const double KaonMass        = 0.493677;
  const double ProtonMass      = 0.93827200;
  const double NeutronMass     = 0.93956563;
  const double CarbonMass      = 12.*AtomicMassUnit;
  const double DeutronMass     = 2.*AtomicMassUnit+0.01313672;
  const double DeltaC11        = 0.0106502;
  const double CoreMass        = 11.*AtomicMassUnit+DeltaC11;
  const double LambdaMass      = 1.115648;
  const double SigmaPMass      = 1.18937;
  const double XiMass          = 1.32134;
  const double ThetaMass       = 1.530;
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
class EventKKAna : public VEvent
{
private:
  RawData      *rawData;
  DCAnalyzer   *DCAna;
  HodoAnalyzer *hodoAna;

public:
        EventKKAna( void );
       ~EventKKAna( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventKKAna::EventKKAna( void )
  : VEvent(),
    rawData(0),
    DCAna( new DCAnalyzer ),
    hodoAna( new HodoAnalyzer )
{
}

//______________________________________________________________________________
EventKKAna::~EventKKAna( void )
{
  if (DCAna)   delete DCAna;
  if (hodoAna) delete hodoAna;
  if (rawData) delete rawData;
}

//______________________________________________________________________________
struct Event
{
  int evnum;
  //Trigger
  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

  //Hodoscope
  int    nhBh2;
  double Bh2Seg[MaxHits];
  double tBh2[MaxHits];
  double deBh2[MaxHits];

  int    nhBh1;
  double Bh1Seg[MaxHits];
  double tBh1[MaxHits];
  double deBh1[MaxHits];
  double btof[MaxHits];

  int    nhTof;
  double TofSeg[MaxHits];
  double tTof[MaxHits];
  double dtTof[MaxHits];
  double deTof[MaxHits];

  int    nhBac;
  double BacSeg[MaxHits];
  double tBac[MaxHits];
  double deBac[MaxHits];

  int    nhPvac;
  double PvacSeg[MaxHits];
  double tPvac[MaxHits];
  double dePvac[MaxHits];

  int    nhFac;
  double FacSeg[MaxHits];
  double tFac[MaxHits];
  double deFac[MaxHits];

  //Fiber
  int    bft_ncl;
  int    bft_clsize[NumOfSegBFT];
  double bft_ctime[NumOfSegBFT];
  double bft_ctot[NumOfSegBFT];
  double bft_clpos[NumOfSegBFT];
  int    sch_ncl;
  // int    sch_hitpat[NumOfSegSCH];
  int    sch_clsize[NumOfSegSCH];
  double sch_ctime[NumOfSegSCH];
  double sch_ctot[NumOfSegSCH];
  double sch_clpos[NumOfSegSCH];
  int    fbh_ncl;
  // int    fbh_hitpat[NumOfSegSCH];
  int    fbh_clsize[NumOfSegCFBH];
  double fbh_ctime[NumOfSegCFBH];
  double fbh_ctot[NumOfSegCFBH];
  double fbh_clpos[NumOfSegCFBH];

  //DC Beam
  int ntBcOut;
  int nlBcOut;
  int nhBcOut[MaxHits];
  double chisqrBcOut[MaxHits];
  double x0BcOut[MaxHits];
  double y0BcOut[MaxHits];
  double u0BcOut[MaxHits];
  double v0BcOut[MaxHits];
  double xtgtBcOut[MaxHits];
  double ytgtBcOut[MaxHits];
  double xbh2BcOut[MaxHits];
  double ybh2BcOut[MaxHits];

  int    ntK18;
  int    nhK18[MaxHits];
  double chisqrK18[MaxHits];
  double pK18[MaxHits];
  double xtb[MaxHits];
  double ytb[MaxHits];
  double utb[MaxHits];
  double vtb[MaxHits];
  double thetab[MaxHits];

  //DC KURAMA
  int ntSdcIn;
  int nlSdcIn;
  int nhSdcIn[MaxHits];
  double chisqrSdcIn[MaxHits];
  double x0SdcIn[MaxHits];
  double y0SdcIn[MaxHits];
  double u0SdcIn[MaxHits];
  double v0SdcIn[MaxHits];
  double xtgtSdcIn[MaxHits];
  double ytgtSdcIn[MaxHits];
  double utgtSdcIn[MaxHits];
  double vtgtSdcIn[MaxHits];

  int ntSdcOut;
  int nlSdcOut;
  int nhSdcOut[MaxHits];
  double chisqrSdcOut[MaxHits];
  double u0SdcOut[MaxHits];
  double v0SdcOut[MaxHits];
  double x0SdcOut[MaxHits];
  double y0SdcOut[MaxHits];

  int    ntKurama;
  int    nhKurama[MaxHits];
  double chisqrKurama[MaxHits];
  double path[MaxHits];
  double pKurama[MaxHits];
  double polarity[MaxHits];
  double m2[MaxHits];
  double xts[MaxHits];
  double yts[MaxHits];
  double uts[MaxHits];
  double vts[MaxHits];
  double thetas[MaxHits];

  // Mass Trigger
  int rm_event;
  int rm_spill;
  int rm_accept;
  int rm_clear;
  int mst_accept;
  int mst_tdc[NumOfSegTOF];
  int mst_tof[NumOfSegTOF];
  int mst_sch[NumOfSegSCH];

  //Reaction
  int    nKn;
  int    nKp;
  int    nKK;
  double vtx[MaxHits];
  double vty[MaxHits];
  double vtz[MaxHits];
  double closeDist[MaxHits];
  double theta[MaxHits];
  double MissMass[MaxHits];
  double MissMassCorr[MaxHits];
  double MissMassCorrDE[MaxHits];
  double BE[MaxHits];
  double theta_CM[MaxHits];
  double cost_CM[MaxHits];

  double xkn[MaxHits];
  double ykn[MaxHits];
  double ukn[MaxHits];
  double vkn[MaxHits];
  double xkp[MaxHits];
  double ykp[MaxHits];
  double ukp[MaxHits];
  double vkp[MaxHits];

  double pOrg[MaxHits];
  double pCalc[MaxHits];
  double pCorr[MaxHits];
  double pCorrDE[MaxHits];
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
EventKKAna::ProcessingBegin( void )
{
  InitializeEvent();
  return true;
}

//______________________________________________________________________________
bool
EventKKAna::ProcessingNormal( void )
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
  static const double MinTimeBFT = gUser.GetParameter("TimeBFT", 0);
  static const double MaxTimeBFT = gUser.GetParameter("TimeBFT", 1);
  static const double MinTimeSCH = gUser.GetParameter("TimeSCH", 0);
  static const double MaxTimeSCH = gUser.GetParameter("TimeSCH", 1);
  static const double MinTimeFBH = gUser.GetParameter("TimeFBH", 0);
  static const double MaxTimeFBH = gUser.GetParameter("TimeFBH", 1);
  static const double OffsetToF  = gUser.GetParameter("OffsetToF");

  rawData = new RawData;
  rawData->DecodeHits();

  gRM.Decode();

  event.evnum = gRM.EventNumber();

  // TFlag
  {
    const HodoRHitContainer &cont=rawData->GetTrigRawHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      int T=hit->GetTdc1();
      if( T>0 ){
	event.trigpat[i] = seg;
	event.trigflag[seg-1] = T;
      }
    }
  }

  HF1( 1, 0. );

  //   if( event.trigflag[SpillEndFlag]>0 ) return true;

  HF1( 1, 1. );

  //////////////BH2 Analysis
  hodoAna->DecodeBH2Hits(rawData);
  int ncBh2 = hodoAna->GetNClustersBH2();
  event.nhBh2 = ncBh2;
  HF1( 101, double(ncBh2) );
#if HodoCut
  if( ncBh2==0 ) return true;
#endif
  double time0 = -9999.;
  {
    int ncOk = 0;
    double mint = -9999.;
    for( int i=0; i<ncBh2; ++i ){
      BH2Cluster *cl = hodoAna->GetClusterBH2(i);
      double seg = cl->MeanSeg()+1;
      double cmt = cl->CMeanTime();
      double ct0 = cl->CTime0();
      double de  = cl->DeltaE();
#if HodoCut
      if( de<MinDeBH2 || MaxDeBH2<de ) continue;
#endif
      ++ncOk;
      HF1( 112, cl->ClusterSize() );
      HF1( 113, seg );
      HF1( 114, cmt );
      HF1( 115, cl->DeltaE() );
      if( i<MaxHits ){
	event.Bh2Seg[i] = seg;
	event.tBh2[i]   = cmt;
	event.deBh2[i]  = de;
      }
      if( std::abs(cmt)<std::abs(mint) ){
	mint  = cmt;
	time0 = ct0;
      }
    }
#if HodoCut
    if( ncOk==0 ) return true;
#endif
  }

  HF1( 1, 2. );

  //////////////BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData);
  int ncBh1=hodoAna->GetNClustersBH1();
  event.nhBh1=ncBh1;
#if HodoCut
  if( ncBh1==0 ) return true;
#endif
  HodoCluster *clBH1Time0=hodoAna->GetClusterBH1(0);
  {
    int ncOk=0;
    double min_tof = clBH1Time0->CMeanTime() - time0;
    HF1( 201, double(ncBh1) );
    for( int i=0; i<ncBh1; ++i ){
      HodoCluster *cl=hodoAna->GetClusterBH1(i);
      double btof= (cl->CMeanTime())-time0;
#if HodoCut
      double de = cl->DeltaE();
      if( de<MinDeBH1 || MaxDeBH1<de ) continue;
      if( btof<MinBeamToF || MaxBeamToF<btof ) continue;
#endif
      HF1( 202, cl->ClusterSize() );
      HF1( 203, cl->MeanSeg()+1-0.5 );
      HF1( 204, cl->CMeanTime() );
      HF1( 205, cl->DeltaE() );
      HF1( 206, btof );

      if( i<MaxHits ){
	event.Bh1Seg[i]=cl->MeanSeg()+1;
	event.tBh1[i]=cl->CMeanTime();
	event.deBh1[i]=cl->DeltaE();
	event.btof[i] = btof;
      }
      ++ncOk;
      HF1( 212, cl->ClusterSize() );
      HF1( 213, cl->MeanSeg()+1-0.5 );
      HF1( 214, cl->CMeanTime() );
      HF1( 215, cl->DeltaE() );
      if( std::abs(btof)<std::abs(min_tof) ){
	clBH1Time0 = cl;
	min_tof = btof;
      }
    }
    HF1( 211, double(ncOk) );
#if HodoCut
    if( ncOk==0 ) return true;
#endif
  }

  HF1( 1, 3. );

  //////////////Tof Analysis
  hodoAna->DecodeTOFHits(rawData);
  int ncTof=hodoAna->GetNClustersTOF();
  event.nhTof=ncTof;
  HF1( 301, double(ncTof) );
#if HodoCut
  if( ncTof==0 ) return true;
#endif
  {
    int ncOk=0;
    for( int i=0; i<ncTof; ++i ){
      HodoCluster *cl=hodoAna->GetClusterTOF(i);
      double de=cl->DeltaE();
      double t=cl->CMeanTime()-time0+OffsetToF, dt=cl->TimeDif();
      HF1( 302, cl->ClusterSize() );
      HF1( 303, cl->MeanSeg()+1-0.5 );
      HF1( 304, t );
      HF1( 305, de );
      if( i<MaxHits ){
	event.TofSeg[i] = cl->MeanSeg()+1;
	event.tTof[i]   = t;
	event.dtTof[i]  = dt;
	event.deTof[i]  = de;
      }
      ++ncOk;
      HF1( 312, cl->ClusterSize() );
      HF1( 313, cl->MeanSeg()+1-0.5 );
      HF1( 314, t );
      HF1( 315, de );
    }
    HF1( 311, double(ncOk) );
#if HodoCut
    if( ncOk==0 ) return true;
#endif
  }

  HF1( 1, 4. );

  ////////// SCH
  {
    hodoAna->DecodeSCHHits(rawData);
    hodoAna->TimeCutSCH( MinTimeSCH, MaxTimeSCH );
    int ncl = hodoAna->GetNClustersSCH();
    event.sch_ncl = ncl;
    for( int i=0; i<ncl; ++i ){
      FiberCluster *cl = hodoAna->GetClusterSCH(i);
      if( !cl ) continue;
      double clsize = cl->ClusterSize();
      double ctime  = cl->CMeanTime();
      double ctot   = cl->Width();
      double pos    = cl->MeanPosition();
      event.sch_clsize[i] = clsize;
      event.sch_ctime[i]  = ctime;
      event.sch_ctot[i]   = ctot;
      event.sch_clpos[i]  = pos;
    }
  }

  HF1( 1, 5. );

  HF1( 1, 6. );

  std::vector<double> xBFT;
  ////////// BFT
  {
    hodoAna->DecodeBFTHits(rawData);
    // Fiber Cluster
    hodoAna->TimeCutBFT(MinTimeBFT, MaxTimeBFT);
    int ncl = hodoAna->GetNClustersBFT();
    event.bft_ncl = ncl;
    for( int i=0; i<ncl; ++i ){
      FiberCluster *cl = hodoAna->GetClusterBFT(i);
      if( !cl ) continue;
      double clsize = cl->ClusterSize();
      double ctime  = cl->CMeanTime();
      double ctot   = cl->Width();
      double pos    = cl->MeanPosition();
      xBFT.push_back( pos );
      event.bft_clsize[i] = clsize;
      event.bft_ctime[i]  = ctime;
      event.bft_ctot[i]   = ctot;
      event.bft_clpos[i]  = pos;
    }
  }

  HF1( 1, 7.);

  if( xBFT.size()!=1 ) return true;

  HF1( 1, 10. );

  ////////////////////////////DC
  double multi_BcOut=0.;
  double multi_SdcIn=0., multi_SdcOut=0.;

  DCAna->DecodeRawHits( rawData );

  ////////////// SDC2&3 number of hit in one layer not 0
  {
    for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetSdcOutHC(layer);
      int nhOut=contOut.size();
      multi_SdcOut += double(nhOut);
      if( nhOut>0 ) event.nlSdcOut++;
    }
    if( multi_SdcOut/double(NumOfLayersSdcOut) > MaxMultiHitSdcOut )
      return true;
  }
  ////////////// BC3&4 number of hit in one layer
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetBcOutHC(layer);
      int nhOut=contOut.size();
      multi_BcOut += double(nhOut);
      if( nhOut>0 ) event.nlBcOut++;
    }
    if( multi_BcOut/double(NumOfLayersBcOut) > MaxMultiHitBcOut )
      return true;
  }
  ////////////// SDC1 number of hit in one layer not 0
  {
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetSdcInHC(layer);
      int nhIn = contIn.size();
      multi_SdcIn  += double(nhIn);
      if( nhIn>0 ) event.nlSdcIn++;
    }
    if( multi_SdcIn/double(NumOfLayersSdcIn) > MaxMultiHitSdcIn )
      return true;
  }

  HF1( 1, 11. );

  //////////////SdcIn Analysis
  DCAna->TrackSearchSdcIn();
  int ntSdcIn=DCAna->GetNtracksSdcIn();
  if( MaxHits<ntSdcIn ){
    std::cout << "#W " << func_name << " "
	      << "too many ntSdcIn " << ntSdcIn << "/" << MaxHits << std::endl;
    ntSdcIn = MaxHits;
  }
  event.ntSdcIn=ntSdcIn;
  HF1( 1001, double(ntSdcIn) );
  if( ntSdcIn<1 || ntSdcIn>MaxNTrackSdcIn )
      return true;
  HF1( 1, 12. );
  {
    int ntOk=0;
    for( int i=0; i<ntSdcIn; ++i ){
      DCLocalTrack *tp=DCAna->GetTrackSdcIn(i);
      if(!tp) continue;
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX(0.), y0=tp->GetY(0.);
      if(i<MaxHits){
	event.nhSdcIn[i]=nh;
	event.chisqrSdcIn[i]=chisqr;
	event.u0SdcIn[i]=u0;
	event.v0SdcIn[i]=v0;
	event.x0SdcIn[i]=x0;
	event.y0SdcIn[i]=y0;
      }
      HF1( 1002, double(nh) );
      HF1( 1003, chisqr );
      HF1( 1004, x0 ); HF1( 1005, y0 );
      HF1( 1006, u0 ); HF1( 1007, v0 );
      HF2( 1008, x0, u0 ); HF2( 1009, y0, v0 );
      HF2( 1010, x0, y0 );
      ++ntOk;
    }
    if( ntOk==0 ) return true;
  }

  HF1( 1, 13. );

  //////////////SdcOut Analysis
  DCAna->TrackSearchSdcOut();
  int ntSdcOut=DCAna->GetNtracksSdcOut();
  if( MaxHits<ntSdcOut ){
    std::cout << "#W " << func_name << " "
	      << "too many ntSdcOut " << ntSdcOut << "/" << MaxHits << std::endl;
    ntSdcOut = MaxHits;
  }
  event.ntSdcOut=ntSdcOut;
  HF1( 1201, double(ntSdcOut) );
  if( ntSdcOut<1 || ntSdcOut>MaxNTrackSdcOut )
    return true;

  HF1( 1, 14. );
  {
    int ntOk=0;
    for( int i=0; i<ntSdcOut; ++i ){
      DCLocalTrack *tp=DCAna->GetTrackSdcOut(i);
      if(!tp) continue;
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX0(), y0=tp->GetY0();

      if( i<MaxHits ){
	event.nhSdcOut[i]     = nh;
	event.chisqrSdcOut[i] = chisqr;
	event.u0SdcOut[i]     = u0;
	event.v0SdcOut[i]     = v0;
	event.x0SdcOut[i]     = x0;
	event.y0SdcOut[i]     = y0;
      }

      HF1( 1202, double(nh) );
      HF1( 1203, chisqr );
      HF1( 1204, x0 ); HF1( 1205, y0 );
      HF1( 1206, u0 ); HF1( 1207, v0 );
      HF2( 1208, x0, u0 ); HF2( 1209, y0, v0 );
      HF2( 1210, x0, y0 );

      ++ntOk;
    }

    if( ntOk==0 ) return true;
  }

  HF1( 1, 15. );

  HF1( 1, 20. );

  //////////////KURAMATracking Analysis
  DCAna->TrackSearchKurama();
  int ntKurama=DCAna->GetNTracksKurama();
  if( MaxHits<ntKurama ){
    std::cout << "#W " << func_name << " "
	      << "too many ntKurama " << ntKurama << "/" << MaxHits << std::endl;
    ntKurama = MaxHits;
  }
  event.ntKurama=ntKurama;
  HF1( 3001, double(ntKurama) );
  if( ntKurama<1 || ntKurama>MaxNTrackKurama )
    return true;

  HF1( 1, 21. );
  {
    int ntOk=0;
    for( int i=0; i<ntKurama; ++i ){
      KuramaTrack *track=DCAna->GetKuramaTrack(i);
      if(!track) continue;
      int nh=track->GetNHits();
      double chisqr=track->chisqr();
      ThreeVector Ppos = track->PrimaryPosition();
      ThreeVector Pmom = track->PrimaryMomentum();
      double path=track->PathLengthToTOF();
      double xt=Ppos.x(), yt=Ppos.y();
      double polarity= track->Polarity();
      double p=Pmom.Mag();
      double ut=Pmom.x()/p, vt=Pmom.y()/p;
      double cost  = 1./sqrt(1.+ut*ut+vt*vt);
      double theta = acos(cost)*math::Rad2Deg();
      if(i<MaxHits){
	event.nhKurama[i] = nh;
	event.chisqrKurama[i] = chisqr;
	event.path[i] = path;
	event.pKurama[i] = p;
	event.polarity[i] = polarity;
	event.xts[i] = xt;
	event.yts[i] = yt;
	event.uts[i] = ut;
	event.vts[i] = vt;
	event.thetas[i] = theta;
      }

      HF1( 3002, double(nh) );
      HF1( 3003, chisqr );
      HF1( 3004, xt ); HF1( 3005, yt );
      HF1( 3006, ut ); HF1( 3007, vt );
      HF2( 3008, xt, ut ); HF2( 3009, yt, vt );
      HF2( 3010, xt, yt );
      HF1( 3011, p );
      HF1( 3012, path );

      double m2;
      for( int j=0; j<ncTof; ++j ){
	HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	if( !clTof || !clTof->GoodForAnalysis() ) continue;
	//double stof = clTof->CMeanTime()-time0+OffsetToF;
	double stof = clTof->C1stTime()-time0+OffsetToF;
	m2 = Kinematics::MassSquare( p, path, stof );
	//pTof = Kinematics::MassSquare( p, path, stof );
	// 		std::cout<<"Mom="<< p <<std::endl;
	// 		std::cout<<"Path="<< path <<std::endl;
	// 		std::cout<<"Time="<< clTof->CMeanTime()-time0 <<std::endl;
	// 		std::cout<<"m2= "<< m2 <<std::endl;
	if( i<MaxHits ) event.m2[i] = m2;
	HF1( 3013, m2 );
      }
      //------------------------Cut
      //if( chisqr<MaxChisqrKurama ){
      //if( MinMassSquare<m2 && m2<MaxMassSquare && chisqr<MaxChisqrKurama ){
      ++ntOk;
      HF1( 3102, double(nh) );
      HF1( 3103, chisqr );
      HF1( 3104, xt ); HF1( 3105, yt );
      HF1( 3106, ut ); HF1( 3107, vt );
      HF2( 3108, xt, ut ); HF2( 3109, yt, vt );
      HF2( 3110, xt, yt );
      HF1( 3111, p );
      HF1( 3112, path );
	//}
      HF1( 3101, double(ntOk) );
      if( ntOk==0 ) return true;
    }
  }

  HF1( 1, 22. );

  HF1( 1, 30. );

  //////////////BcOut Analysis
  DCAna->TrackSearchBcOut();
  int ntBcOut=DCAna->GetNtracksBcOut();
  event.ntBcOut=ntBcOut;
  HF1( 1101, double(ntBcOut) );
  if( ntBcOut<1 || ntBcOut>MaxNTrackBcOut )
    return true;

  HF1( 1, 31. );
  {
    int ntOk=0;
    for( int i=0; i<ntBcOut; ++i ){
      DCLocalTrack *tp=DCAna->GetTrackBcOut(i);
      if(!tp) continue;
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX(0.), y0=tp->GetY(0.);
      double xtgt=tp->GetX( zK18Target ), ytgt=tp->GetY( zK18Target );

      if( i<MaxHits ){
	event.nhBcOut[i]=nh;
	event.chisqrBcOut[i]=chisqr;
	event.u0BcOut[i]=u0;
	event.v0BcOut[i]=v0;
	event.x0BcOut[i]=x0;
	event.y0BcOut[i]=y0;
	event.xtgtBcOut[i]=xtgt;
	event.ytgtBcOut[i]=ytgt;
      }

      HF1( 1102, double(nh) );
      HF1( 1103, chisqr );
      HF1( 1104, x0 ); HF1( 1105, y0 );
      HF1( 1106, u0 ); HF1( 1107, v0 );
      HF2( 1108, x0, u0 ); HF2( 1109, y0, v0 );
      HF2( 1110, x0, y0 );
      HF1( 1111, xtgt ); HF1( 1112, ytgt );
      HF2( 1113, xtgt, ytgt );
      if( chisqr<MaxChisqrBcOut ) ++ntOk;
    }
    if( ntOk==0 ) return true;
  }

  HF1( 1, 32. );

  // Mass Trigger
  // bool rm_accept  = false;
#if 0
  bool mst_accept = false;
  {
    const HodoRHitContainer &cont=rawData->GetMsTRMRawHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int d0=hit->GetAdcUp(); int d1=hit->GetAdcDown();
      int d2=hit->GetTdcUp(); int d3=hit->GetTdcDown();
      event.rm_event  = d0; event.rm_spill  = d1;
      event.rm_accept = d2; event.rm_clear  = d3;
      // if( d2==1 ) rm_accept = true;
    }
  }

  // CAMAC TOF TDC
  {
    const HodoRHitContainer &cont = rawData->GetMsTRawHC(0);
    int nh = cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      if( !hit ) continue;
      int seg = hit->SegmentId();
      int tdc = hit->GetAdcUp();
      event.mst_tdc[seg] = tdc;
      event.mst_tof[i]   = seg+1;
      const HodoRHitContainer &cont2 = rawData->GetMsTRawHC(1);
      int nh2 = cont2.size();
      for( int i2=0; i2<nh2; ++i2 ){
	HodoRawHit *hit2 = cont2[i2];
	if( !hit2 ) continue;
	int seg2 = hit2->SegmentId();
	int coin = hit2->GetAdcUp();
	if( coin && !mst_accept )
	  mst_accept = gMsT.IsAccept( seg, seg2, tdc );
      }
    }
  }
  event.mst_accept = mst_accept;

  // CAMAC SCH Coin
  {
    const HodoRHitContainer &cont = rawData->GetMsTRawHC(1);
    int nh = cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      int seg = hit->SegmentId();
      int coin = hit->GetAdcUp();
      if( coin ) event.mst_sch[seg] = 1;
    }
  }
#endif

  HF1( 1, 35. );

  // K18TrackingD2U
  DCAna->TrackSearchK18D2U( xBFT );
  int ntK18=DCAna->GetNTracksK18D2U();
  event.ntK18 = ntK18;
  HF1( 2201, double(ntK18) );
  if( ntK18>MaxNTrackK18 ){
    std::cout << "#W : " << func_name << " "
	      << "too many ntK18 : " << ntK18 << std::endl;
    ntK18 = MaxNTrackK18;
  }
  if( ntK18<1 || ntK18>MaxNTrackK18 )
    return true;

  HF1( 1, 36. );

  for( int i=0; i<ntK18; ++i ){
    K18TrackD2U *tp=DCAna->GetK18TrackD2U(i);
    if(!tp) continue;
    double xt=tp->Xtgt(), yt=tp->Ytgt();
    double ut=tp->Utgt(), vt=tp->Vtgt();
    double pk18=tp->P3rd();
    double cost  = 1./sqrt(1.+ut*ut+vt*vt);
    double theta = acos(cost)*math::Rad2Deg();
    event.pK18[i]   = pk18;
    event.xtb[i]    = xt;
    event.ytb[i]    = yt;
    event.utb[i]    = ut;
    event.vtb[i]    = vt;
    event.thetab[i] = theta;
    HF1( 2204, pk18 );
  }

  HF1( 1, 37. );

  HF1( 1, 40. );

  std::vector <ThreeVector> KnPCont, KnXCont;
  std::vector <ThreeVector> KpPCont, KpXCont;

  ////////// Kaon Plus
  for( int itKurama=0; itKurama<ntKurama; ++itKurama ){
    KuramaTrack *track=DCAna->GetKuramaTrack( itKurama );
    if( !track || !track->GoodForAnalysis() ) continue;
    // DCLocalTrack *trIn =track->GetLocalTrackIn();
    // DCLocalTrack *trOut=track->GetLocalTrackOut();
    HodoCluster *clTof=0;
    // double x0 = trOut->GetX0();
    // double y0 = trOut->GetY0();
    // double mindifTof = MinXDifTof;
    if( !clTof ) continue;
    int nh=track->GetNHits();
    double chisqr=track->chisqr();
    ThreeVector Ppos = track->PrimaryPosition();
    ThreeVector Pmom = track->PrimaryMomentum();
    //Calibration
    double p0 = Pmom.Mag();//+pKURAMA_offset+pS_offset;
    double px = Pmom.x();
    double py = Pmom.y();
    double pz = Pmom.z();
    double u0 = px/pz, v0 = py/pz;
    double pt = p0/sqrt(1.+u0*u0+v0*v0);

    ThreeVector PposCorr( Ppos.x()+x_off, Ppos.y()+y_off, Ppos.z() );
    ThreeVector PmomCorr( pt*(u0+u_off), pt*(v0+v_off), pt);

    double pathL = track->PathLengthToTOF();
    double xt = PposCorr.x();
    double yt = PposCorr.y();
    double zt = PposCorr.z();
    double p  = PmomCorr.Mag();
    // double q = track->Polarity();
    double ut = xt/zt;
    double vt = yt/zt;
    double stof = clTof->CMeanTime()-time0+OffsetToF;
    double m2 = Kinematics::MassSquare( p0, pathL, stof );
    HF1( 3202, double(nh) );
    HF1( 3203, chisqr );
    HF1( 3204, xt ); HF1( 3205, yt );
    HF1( 3206, ut ); HF1( 3207, vt );
    HF2( 3208, xt, ut ); HF2( 3209, yt, vt );
    HF2( 3210, xt, yt );
    HF1( 3211, p );
    HF1( 3212, pathL );
    //HF1( 3213, m2 );

    //     double xTof=(clTof->MeanSeg()-7.5)*70.;
    //     double yTof=(clTof->TimeDif())*800./12.;
    //     HF2( 4011, clTof->MeanSeg()-0.5, xtof );
    //     HF2( 4013, clTof->TimeDif(), ytof );
    //     HF1( 4015, xtof-xTof ); HF1( 4017, ytof-yTof );
    //     HF2( 4019, xtof-xTof, ytof-yTof );
    if( MinMassSquare<m2 && m2<MaxMassSquare ){
	//       double ttof=clTof->CMeanTime()-time0,
	//       HF1( 322, clTof->ClusterSize() );
	//       HF1( 323, clTof->MeanSeg()-0.5 );
	//       HF1( 324, ttof );
	//       HF1( 325, clTof->DeltaE() );
	//       double u0in=trIn->GetU0();
	//       HF2( 4001, u0in, ttof ); HF2( 4003, u0in, ttof+12.5*u0in );
	HF1( 4202, double(nh) );
	HF1( 4203, chisqr );
	HF1( 4204, xt ); HF1( 4205, yt );
	HF1( 4206, ut ); HF1( 4207, vt );
	HF2( 4208, xt, ut ); HF2( 4209, yt, vt );
	HF2( 4210, xt, yt );
	HF1( 4211, p );
	HF1( 4212, pathL );
	KpPCont.push_back(PmomCorr);
	KpXCont.push_back(PposCorr);
    }
  }

  if( KpPCont.empty() ) return true;

  HF1( 1, 41. );

  ////////// Kaon Minus
  for( int itK18=0; itK18<ntK18; ++itK18 ){
    K18TrackD2U *track = DCAna->GetK18TrackD2U( itK18 );
    if( !track || !track->GoodForAnalysis() ) continue;
    DCLocalTrack *trOut=track->TrackOut();
    //Calibration
    double p = track->P3rd()*pB_offset+pK18_offset;
    double x=track->Xtgt(), y=track->Ytgt();
    double u=track->Utgt(), v=track->Vtgt();
    // double loss_bh2 = 1.09392e-3;
    // p = p - loss_bh2;
    double pt=p/sqrt(1.+u*u+v*v);
    ThreeVector Pos( x, y, 0. );
    ThreeVector Mom( pt*u, pt*v, pt );
    double xo=trOut->GetX0(), yo=trOut->GetY0();
    HodoCluster *clBh1=0;
    for( int j=0; j<ncBh1; ++j ){
      HodoCluster *cl=hodoAna->GetClusterBH1(j);
      clBh1=cl;
    }
    if( !clBh1 ) continue;
    HF1( 4104, p );
    HF1( 4105, x ); HF1( 4106, y );
    HF1( 4107, xo ); HF1( 4108, yo ); HF1( 4109, u ); HF1( 4110, v );

    KnPCont.push_back(Mom); KnXCont.push_back(Pos);
  }

  if( KnPCont.empty() ) return true;

  HF1( 1, 42. );

#if 1
  //MissingMass
  int nKn = KnPCont.size();
  int nKp = KpPCont.size();
  event.nKn = nKn;
  event.nKp = nKp;
  event.nKK = nKn*nKp;
  HF1( 4101, double(nKn));
  HF1( 4201, double(nKp) );

  int nkk=0;
  for( int ikp=0; ikp<nKp; ++ikp ){
    ThreeVector pkp = KpPCont[ikp], xkp = KpXCont[ikp];
    for( int ikn=0; ikn<nKn; ++ikn ){
      ThreeVector pkn  = KnPCont[ikn], xkn = KnXCont[ikn];
      ThreeVector vert = Kinematics::VertexPoint( xkn, xkp, pkn, pkp );
      double closedist = Kinematics::closeDist( xkn, xkp, pkn, pkp );

      double us = pkp.x()/pkp.z(), vs = pkp.y()/pkp.z();
      double ub = pkn.x()/pkn.z(), vb = pkn.y()/pkn.z();
      double cost = pkn*pkp/(pkn.Mag()*pkp.Mag());

      double pk0   = pkp.Mag();
      double pCorr = pk0;

      ThreeVector pkpCorr( pCorr*pkp.x()/pkp.Mag(),
			   pCorr*pkp.y()/pkp.Mag(),
			   pCorr*pkp.z()/pkp.Mag() );

      ThreeVector pknCorrDE = Kinematics::CorrElossIn( pkn, xkn, vert, KaonMass );
      ThreeVector pkpCorrDE = Kinematics::CorrElossOut( pkpCorr, xkp, vert, KaonMass );

      LorentzVector LvKn( pkn, sqrt( KaonMass*KaonMass+pkn.Mag2() ) );
      LorentzVector LvKnCorrDE( pknCorrDE, sqrt( KaonMass*KaonMass+pknCorrDE.Mag2() ) );

      LorentzVector LvKp( pkp, sqrt( KaonMass*KaonMass+pkp.Mag2() ) );
      LorentzVector LvKpCorr( pkpCorr, sqrt( KaonMass*KaonMass+pkpCorr.Mag2() ) );
      LorentzVector LvKpCorrDE( pkpCorrDE, sqrt( KaonMass*KaonMass+pkpCorrDE.Mag2() ) );

      LorentzVector LvC( 0., 0., 0., ProtonMass );
      //LorentzVector LvC( 0., 0., 0., CarbonMass );
      LorentzVector LvCore( 0., 0., 0., 0. );

      LorentzVector LvRc       = LvKn+LvC-LvKp;
      LorentzVector LvRcCorr   = LvKn+LvC-LvKpCorr;
      LorentzVector LvRcCorrDE = LvKnCorrDE+LvC-LvKpCorrDE;
      double MisMass       = LvRc.Mag();//-LvC.Mag();
      double MisMassCorr   = LvRcCorr.Mag();//-LvC.Mag();
      double MisMassCorrDE = LvRcCorrDE.Mag();//-LvC.Mag();
      double BE            = LvRc.Mag()-( LvCore.Mag()+KaonMass );
      // double BECorr        = LvRcCorr.Mag()-( LvCore.Mag()+KaonMass );
      // double BECorrDE      = LvRcCorrDE.Mag()-( LvCore.Mag()+KaonMass );

      // double BE               = LvRc.Mag()-( CoreMass+LambdaMass );
      // double BECorr           = LvRcCorr.Mag()-( CoreMass+LambdaMass );
      // double BECorrDE         = LvRcCorrDE.Mag()-( CoreMass+LambdaMass );

      // std::cout<<"******* Missing Mass = "<< MisMass
      // 	       << " (" << MisMassCorr << ")"
      // 	       << " (" << MisMassCorrDE << ")" <<std::endl;
      // std::cout<<"******* Binding Energy = "<<BE
      // 	       << " (" << BECorr << ")"
      // 	       << " (" << BECorrDE <<")" <<std::endl;

      //Primary frame
      LorentzVector PrimaryLv = LvKn+LvC;
      double TotalEnergyCM = PrimaryLv.Mag();
      ThreeVector beta( 1/PrimaryLv.E()*PrimaryLv.Vect() );

      //CM
      double TotalMomCM
	= 0.5*sqrt(( TotalEnergyCM*TotalEnergyCM
		     -( KaonMass+XiMass )*( KaonMass+XiMass ))
		   *( TotalEnergyCM*TotalEnergyCM
		      -( KaonMass-XiMass )*( KaonMass-XiMass )))/TotalEnergyCM;

      double costLab = cost;
      double cottLab = costLab/sqrt(1.-costLab*costLab);
      double bt=beta.Mag(), gamma=1./sqrt(1.-bt*bt);
      double gbep=gamma*bt*sqrt(TotalMomCM*TotalMomCM+KaonMass*KaonMass)/TotalMomCM;
      double a  = gamma*gamma+cottLab*cottLab;
      double bp = gamma*gbep;
      double c  = gbep*gbep-cottLab*cottLab;
      double dd = bp*bp-a*c;

      if( dd<0. ){
	std::cerr << "dd<0." << std::endl;
	dd = 0.;
      }

      double costCM = (sqrt(dd)-bp)/a;
      if( costCM>1. || costCM<-1. ){
	std::cerr << "costCM>1. || costCM<-1." << std::endl;
	costCM=-1.;
      }
      double sintCM=sqrt(1.-costCM*costCM);
      double KaonMom = TotalMomCM*sintCM/sqrt(1.-costLab*costLab);

      if (nkk<MaxHits) {
	event.vtx[nkk]       = vert.x();
	event.vty[nkk]       = vert.y();
	event.vtz[nkk]       = vert.z();
	event.closeDist[nkk] = closedist;
	event.theta[nkk]     = acos(cost)*math::Rad2Deg();
	event.theta_CM[nkk]  = acos(costCM)*math::Rad2Deg();
	event.cost_CM[nkk]   = costCM;

	event.MissMass[nkk]=MisMass;
	event.MissMassCorr[nkk] = MisMassCorr;
	event.MissMassCorrDE[nkk] = MisMassCorrDE;

	event.BE[nkk]=BE;

	event.xkp[nkk] = xkp.x();
	event.ykp[nkk] = xkp.y();
	event.ukp[nkk] = us;
	event.vkp[nkk] = vs;

	event.xkn[nkk] = xkn.x();
	event.ykn[nkk] = xkn.y();
	event.ukn[nkk] = ub;
	event.vkn[nkk] = vb;
	event.pOrg[nkk] = pk0;
	event.pCalc[nkk] = KaonMom;
	event.pCorr[nkk] = pCorr;
	event.pCorrDE[nkk] = pkpCorrDE.Mag();
	nkk++;
      }

      HF1( 5001, vert.z() );

      HF1( 5002, MisMass );
      HF2( 5011, MisMass, us );
      HF2( 5012, MisMass, vs );
      HF2( 5013, MisMass, ub );
      HF2( 5014, MisMass, vb );
    }
  }

  HF1( 1, 50. );
#endif

  //Final Hodoscope histograms
  for( int i=1; i<ncBh2; ++i ){
    BH2Cluster *cl=hodoAna->GetClusterBH2(i);
    double t=cl->CMeanTime();
#if HodoCut
    double de = cl->DeltaE();
    if( de<MinDeBH2 || MaxDeBH2<de ) continue;
#endif
    HF1( 152, cl->ClusterSize() );
    HF1( 153, cl->MeanSeg()+1-0.5 );
    HF1( 154, t );
    HF1( 155, cl->DeltaE() );
  }

  for( int i=0; i<ncBh1; ++i ){
    HodoCluster *cl=hodoAna->GetClusterBH1(i);
    double btof= (cl->CMeanTime())-time0;
#if HodoCut
    if( !cl || !cl->GoodForAnalysis() ) continue;
#endif
    HF1( 252, cl->ClusterSize() );
    HF1( 253, cl->MeanSeg()+1-0.5 );
    HF1( 254, cl->CMeanTime() );
    HF1( 255, cl->DeltaE() );
    HF1( 256, btof );
  }

  for( int i=0; i<ncTof; ++i ){
    HodoCluster *cl=hodoAna->GetClusterTOF(i);
    double de=cl->DeltaE();
    double t=cl->CMeanTime()-time0;
#if HodoCut
    if( !cl || !cl->GoodForAnalysis() ) continue;
#endif
    HF1( 352, cl->ClusterSize() );
    HF1( 353, cl->MeanSeg()+1-0.5 );
    HF1( 354, t );
    HF1( 355, de );
  }

  return true;
}

//______________________________________________________________________________
void
EventKKAna::InitializeEvent( void )
{
  event.evnum      = 0;
  event.rm_event   = -1;
  event.rm_spill   = -1;
  event.rm_accept  = -1;
  event.rm_clear   = -1;
  event.mst_accept = -1;
  for( int i=0; i<NumOfSegTOF; ++i ){
    event.mst_tof[i] = -1;
    event.mst_tdc[i] = -1;
  }
  for( int i=0; i<NumOfSegSCH; ++i ){
    event.mst_sch[i] = -1;
  }

  //Trigger
  for( int it=0; it<NumOfSegTrig; ++it ){
    event.trigpat[it]  = -1;
    event.trigflag[it] = -1;
  }

  //Hodoscope
  event.nhBh2  = 0;
  event.nhBh1  = 0;
  event.nhTof  = 0;
  event.nhBac  = 0;
  event.nhPvac = 0;
  event.nhFac  = 0;

  for( int it=0; it<MaxHits; it++ ){
    event.Bh2Seg[it] = -1;
    event.tBh2[it]   = -9999.0;
    event.deBh2[it]  = -9999.0;

    event.Bh1Seg[it] = -1;
    event.tBh1[it]   = -9999.0;
    event.deBh1[it]  = -9999.0;
    event.btof[it]   = -1;

    event.TofSeg[it] = -1;
    event.tTof[it]   = -9999.0;
    event.dtTof[it]  = -9999.0;
    event.deTof[it]  = -9999.0;

    event.BacSeg[it] = -1;
    event.tBac[it]   = -9999.0;
    event.deBac[it]  = -9999.0;

    event.PvacSeg[it] = -1;
    event.tPvac[it]   = -9999.0;
    event.dePvac[it]  = -9999.0;

    event.FacSeg[it] = -1;
    event.tFac[it]   = -9999.0;
    event.deFac[it]  = -9999.0;
  }

  //Fiber
  event.bft_ncl = 0;
  event.sch_ncl = 0;
  event.fbh_ncl = 0;
  for( int it=0; it<NumOfSegBFT; it++ ){
    event.bft_clsize[it] = -999;
    event.bft_ctime[it]  = -999.;
    event.bft_ctot[it]   = -999.;
    event.bft_clpos[it]  = -999.;
  }
  for( int it=0; it<NumOfSegSCH; it++ ){
    event.sch_clsize[it] = -999;
    event.sch_ctime[it]  = -999.;
    event.sch_ctot[it]   = -999.;
    event.sch_clpos[it]  = -999.;
  }
  for( int it=0; it<NumOfSegFBH; it++ ){
    event.fbh_clsize[it] = -999;
    event.fbh_ctime[it]  = -999.;
    event.fbh_ctot[it]   = -999.;
    event.fbh_clpos[it]  = -999.;
  }

  //DC
  event.nlBcOut  = 0;
  event.nlSdcIn  = 0;
  event.nlSdcOut = 0;
  event.ntBcOut  = 0;
  event.ntSdcIn  = 0;
  event.ntSdcOut = 0;
  event.ntK18    = 0;
  event.ntKurama = 0;

  //Beam DC
  for( int it=0; it<MaxHits; it++){
    event.nhBcOut[it]     = 0;
    event.chisqrBcOut[it] = -1.0;
    event.x0BcOut[it] = -9999.0;
    event.y0BcOut[it] = -9999.0;
    event.u0BcOut[it] = -9999.0;
    event.v0BcOut[it] = -9999.0;

    event.xtgtBcOut[it] = -9999.0;
    event.ytgtBcOut[it] = -9999.0;
    event.xbh2BcOut[it] = -9999.0;
    event.ybh2BcOut[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhK18[it]     = 0;
    event.chisqrK18[it] = -1.0;
    event.xtb[it] = -9999.;
    event.ytb[it] = -9999.;
    event.utb[it] = -9999.;
    event.vtb[it] = -9999.;
    event.pK18[it]   = -9999.;
    event.thetab[it] = -9999.;
  }

  //KURAMA DC
  for( int it=0; it<MaxHits; it++){
    event.nhSdcIn[it]     = 0;
    event.chisqrSdcIn[it] = -1.;
    event.x0SdcIn[it] = -9999.;
    event.y0SdcIn[it] = -9999.;
    event.u0SdcIn[it] = -9999.;
    event.v0SdcIn[it] = -9999.;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhSdcOut[it]     = 0;
    event.chisqrSdcOut[it] = -1.;
    event.u0SdcOut[it] = -9999.;
    event.v0SdcOut[it] = -9999.;
    event.x0SdcOut[it] = -9999.;
    event.y0SdcOut[it] = -9999.;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhKurama[it]     = 0;
    event.chisqrKurama[it] = -1.;
    event.xts[it] = -9999.;
    event.yts[it] = -9999.;
    event.uts[it] = -9999.;
    event.vts[it] = -9999.;

    event.pKurama[it]  = -9999.;
    event.polarity[it] = 0.;
    event.path[it]     = -9999.;
    event.m2[it]       = -9999.;
    event.thetas[it]   = -9999.;
  }

  //Reaction
  event.nKn = 0;
  event.nKp = 0;
  event.nKK = 0;

  for( int it=0; it<MaxHits; it++){
    event.vtx[it]       = -9999.;
    event.vty[it]       = -9999.;
    event.vtz[it]       = -9999.;
    event.closeDist[it] = -9999.;
    event.theta[it]     = -9999.;
    event.theta_CM[it]  = -9999.;
    event.cost_CM[it]   = -9999.;
    event.MissMass[it]  = -9999.;
    event.MissMassCorr[it]  = -9999.;
    event.MissMassCorrDE[it]  = -9999.;
    event.BE[it]        = -9999.0;

    event.xkn[it] = -9999.0;
    event.ykn[it] = -9999.0;
    event.ukn[it] = -9999.0;
    event.vkn[it] = -9999.0;
    event.xkp[it] = -9999.0;
    event.ykp[it] = -9999.0;
    event.ukp[it] = -9999.0;
    event.vkp[it] = -9999.0;
    event.pOrg[it] = -9999.0;
    event.pCalc[it] = -9999.0;
    event.pCorr[it] = -9999.0;
    event.pCorrDE[it] = -9999.0;
  }

}

//______________________________________________________________________________
bool
EventKKAna::ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventKKAna;
}

//______________________________________________________________________________
bool
ConfMan::InitializeHistograms( void )
{
  HB1( 1, "Status", 60, 0., 60. );

  HB1( 101, "#Clusters BH2",   7, 0., 7. );
  HB1( 102, "ClusterSize BH2", 7, 0., 7. );
  HB1( 103, "HitPat BH2", 10, 0., 10. );
  HB1( 104, "MeanTime BH2", 200, -10., 10. );
  HB1( 105, "Delta-E BH2", 200, -0.5, 4.5 );

  HB1( 112, "ClusterSize BH2", 7, 0., 7. );
  HB1( 113, "HitPat BH2", 10, 0., 10. );
  HB1( 114, "MeanTime BH2", 200, -10., 10. );
  HB1( 115, "Delta-E BH2", 200, -0.5, 4.5 );

  HB1( 122, "ClusterSize BH2 [T0]", 7, 0., 7. );
  HB1( 123, "HitPat BH2 [T0]", 10, 0., 10. );
  HB1( 124, "MeanTime BH2 [T0]", 200, -10., 10. );
  HB1( 125, "Delta-E BH2 [T0]", 200, -0.5, 4.5 );

  HB1( 152, "ClusterSize BH2 [kk]", 7, 0., 7. );
  HB1( 153, "HitPat BH2 [kk]", 10, 0., 10. );
  HB1( 154, "MeanTime BH2 [kk]", 200, -10., 10. );
  HB1( 155, "Delta-E BH2 [kk]", 200, -0.5, 4.5 );

  HB1( 201, "#Clusters BH1",  11, 0., 11. );
  HB1( 202, "ClusterSize BH1",11, 0., 11. );
  HB1( 203, "HitPat BH1", 11, 0., 11. );
  HB1( 204, "MeanTime BH1", 200, -10., 10. );
  HB1( 205, "Delta-E BH1", 200, -0.5, 4.5 );
  HB1( 206, "Beam ToF", 200, -10., 10. );

  HB1( 211, "#Clusters BH1 [pi]",  11, 0., 11. );
  HB1( 212, "ClusterSize BH1 [pi]",11, 0., 11. );
  HB1( 213, "HitPat BH1 [pi]", 11, 0., 11. );
  HB1( 214, "MeanTime BH1 [pi]", 200, -10., 10. );
  HB1( 215, "Delta-E BH1 [pi]", 200, -0.5, 4.5 );

  HB1( 252, "ClusterSize BH1 [kk]",11, 0., 11. );
  HB1( 253, "HitPat BH1 [kk]", 11, 0., 11. );
  HB1( 254, "MeanTime BH1 [kk]", 200, -10., 10. );
  HB1( 255, "Delta-E BH1 [kk]", 200, -0.5, 4.5 );
  HB1( 256, "Beam ToF [kk]", 200, -10., 10. );

  HB1( 301, "#Clusters Tof",  32, 0., 32. );
  HB1( 302, "ClusterSize Tof",32, 0., 32. );
  HB1( 303, "HitPat Tof", 32, 0., 32. );
  HB1( 304, "TimeOfFlight Tof", 500, -50., 100. );
  HB1( 305, "Delta-E Tof", 200, -0.5, 4.5 );

  HB1( 311, "#Clusters Tof [Good]",  32, 0., 32. );
  HB1( 312, "ClusterSize Tof [Good]",32, 0., 32. );
  HB1( 313, "HitPat Tof [Good]", 32, 0., 32. );
  HB1( 314, "TimeOfFlight Tof [Good]", 500, -50., 100. );
  HB1( 315, "Delta-E Tof [Good]", 200, -0.5, 4.5 );

  HB1( 352, "ClusterSize Tof [kk]",32, 0., 32. );
  HB1( 353, "HitPat Tof [kk]", 32, 0., 32. );
  HB1( 354, "TimeOfFlight Tof [kk]", 500, -50., 100. );
  HB1( 355, "Delta-E Tof [kk]", 200, -0.5, 4.5 );

  HB1( 1001, "#Tracks SdcIn", 10, 0., 10. );
  HB1( 1002, "#Hits SdcIn", 20, 0., 20. );
  HB1( 1003, "Chisqr SdcIn", 200, 0., 100. );
  HB1( 1004, "X0 SdcIn", 500, -200., 200. );
  HB1( 1005, "Y0 SdcIn", 500, -200., 200. );
  HB1( 1006, "U0 SdcIn",  700, -0.35, 0.35 );
  HB1( 1007, "V0 SdcIn",  400, -0.20, 0.20 );
  HB2( 1008, "X0%U0 SdcIn", 120, -200., 200., 100, -0.35, 0.35 );
  HB2( 1009, "Y0%V0 SdcIn", 100, -200., 200., 100, -0.20, 0.20 );
  HB2( 1010, "X0%Y0 SdcIn", 100, -200., 200., 100, -200., 200. );

  HB1( 1101, "#Tracks BcOut", 10, 0., 10. );
  HB1( 1102, "#Hits BcOut", 20, 0., 20. );
  HB1( 1103, "Chisqr BcOut", 200, 0., 100. );
  HB1( 1104, "X0 BcOut", 500, -200., 200. );
  HB1( 1105, "Y0 BcOut", 500, -200., 200. );
  HB1( 1106, "U0 BcOut",  700, -0.35, 0.35 );
  HB1( 1107, "V0 BcOut",  400, -0.20, 0.20 );
  HB2( 1108, "X0%U0 BcOut", 120, -200., 200., 100, -0.35, 0.35 );
  HB2( 1109, "Y0%V0 BcOut", 100, -200., 200., 100, -0.20, 0.20 );
  HB2( 1110, "X0%Y0 BcOut", 100, -200., 200., 100, -200., 200. );
  HB1( 1111, "Xtgt BcOut", 500, -200., 200. );
  HB1( 1112, "Ytgt BcOut", 500, -200., 200. );
  HB2( 1113, "Xtgt%Ytgt BcOut", 100, -200., 200., 100, -200., 200. );

  HB1( 1201, "#Tracks SdcOut", 10, 0., 10. );
  HB1( 1202, "#Hits SdcOut", 20, 0., 20. );
  HB1( 1203, "Chisqr SdcOut", 200, 0., 100. );
  HB1( 1204, "X0 SdcOut", 600, -1200., 1200. );
  HB1( 1205, "Y0 SdcOut", 600, -600., 600. );
  HB1( 1206, "U0 SdcOut",  700, -0.35, 0.35 );
  HB1( 1207, "V0 SdcOut",  400, -0.20, 0.20 );
  HB2( 1208, "X0%U0 SdcOut", 120, -1200., 1200., 100, -0.35, 0.35 );
  HB2( 1209, "Y0%V0 SdcOut", 100,  -600.,  600., 100, -0.20, 0.20 );
  HB2( 1210, "X0%Y0 SdcOut", 100, -1200., 1200., 100, -600., 600. );

  HB1( 2201, "#Tracks K18", 10, 0., 10. );
  HB1( 2202, "#Hits K18", 30, 0., 30. );
  HB1( 2203, "Chisqr K18", 500, 0., 50. );
  HB1( 2204, "P K18", 1000, 0.5, 2.0 );
  HB1( 2251, "#Tracks K18 [Good]", 10, 0., 10. );

  HB1( 3001, "#Tracks Kurama", 10, 0., 10. );
  HB1( 3002, "#Hits Kurama", 30, 0., 30. );
  HB1( 3003, "Chisqr Kurama", 500, 0., 500. );
  HB1( 3004, "Xtgt Kurama", 500, -200., 200. );
  HB1( 3005, "Ytgt Kurama", 500, -100., 100. );
  HB1( 3006, "Utgt Kurama", 400, -0.35, 0.35 );
  HB1( 3007, "Vtgt Kurama", 200, -0.20, 0.20 );
  HB2( 3008, "Xtgt%U Kurama", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 3009, "Ytgt%V Kurama", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 3010, "Xtgt%Ytgt Kurama", 100, -200., 200., 100, -100., 100. );
  HB1( 3011, "P Kurama", 200, 0.0, 1.0 );
  HB1( 3012, "PathLength Kurama", 600, 3000., 6000. );
  HB1( 3013, "MassSqr", 600, -1.2, 1.2 );

  HB1( 3101, "#Tracks Kurama [Good]", 10, 0., 10. );
  HB1( 3102, "#Hits Kurama [Good]", 30, 0., 30. );
  HB1( 3103, "Chisqr Kurama [Good]", 500, 0., 500. );
  HB1( 3104, "Xtgt Kurama [Good]", 500, -200., 200. );
  HB1( 3105, "Ytgt Kurama [Good]", 500, -100., 100. );
  HB1( 3106, "Utgt Kurama [Good]", 700, -0.35, 0.35 );
  HB1( 3107, "Vtgt Kurama [Good]", 400, -0.20, 0.20 );
  HB2( 3108, "Xtgt%U Kurama [Good]", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 3109, "Ytgt%V Kurama [Good]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 3110, "Xtgt%Ytgt Kurama [Good]", 100, -200., 200., 100, -100., 100. );
  HB1( 3111, "P Kurama [Good]", 200, 0.0, 1.0 );
  HB1( 3112, "PathLength Kurama [Good]", 600, 3000., 6000. );
  HB1( 3113, "MassSqr", 600, -1.2, 1.2 );

  HB1( 3201, "#Tracks Kurama [Good2]", 10, 0., 10. );
  HB1( 3202, "#Hits Kurama [Good2]", 30, 0., 30. );
  HB1( 3203, "Chisqr Kurama [Good2]", 500, 0., 500. );
  HB1( 3204, "Xtgt Kurama [Good2]", 500, -200., 200. );
  HB1( 3205, "Ytgt Kurama [Good2]", 500, -100., 100. );
  HB1( 3206, "Utgt Kurama [Good2]", 700, -0.35, 0.35 );
  HB1( 3207, "Vtgt Kurama [Good2]", 400, -0.20, 0.20 );
  HB2( 3208, "Xtgt%U Kurama [Good2]", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 3209, "Ytgt%V Kurama [Good2]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 3210, "Xtgt%Ytgt Kurama [Good2]", 100, -200., 200., 100, -100., 100. );
  HB1( 3211, "P Kurama [Good2]", 200, 0.0, 1.0 );
  HB1( 3212, "PathLength Kurama [Good2]", 600, 3000., 6000. );
  HB1( 3213, "MassSqr", 600, -1.2, 1.2 );

  HB1( 4101, "#Tracks K18 [KuramaP]", 10, 0., 10. );
  HB1( 4102, "#Hits K18 [KuramaP]", 30, 0., 30. );
  HB1( 4103, "Chisqr K18 [KuramaP]", 500, 0., 50. );
  HB1( 4104, "P K18 [KuramaP]", 500, 0.5, 2.0 );
  HB1( 4105, "Xtgt K18 [KuramaP]", 500, -200., 200. );
  HB1( 4106, "Ytgt K18 [KuramaP]", 500, -100., 100. );
  HB1( 4107, "Xout K18 [KuramaP]", 500, -200., 200. );
  HB1( 4108, "Yout K18 [KuramaP]", 500, -100., 100. );
  HB1( 4109, "Uout K18 [KuramaP]", 400, -0.35, 0.35 );
  HB1( 4110, "Vout K18 [KuramaP]", 200, -0.20, 0.20 );
  HB1( 4111, "Xin  K18 [KuramaP]", 500, -100., 100. );
  HB1( 4112, "Yin  K18 [KuramaP]", 500, -100., 100. );
  HB1( 4113, "Uin  K18 [KuramaP]", 400, -0.35, 0.35 );
  HB1( 4114, "Vin  K18 [KuramaP]", 200, -0.20, 0.20 );

  HB1( 4201, "#Tracks Kurama [Proton]", 10, 0., 10. );
  HB1( 4202, "#Hits Kurama [Proton]", 30, 0., 30. );
  HB1( 4203, "Chisqr Kurama [Proton]", 500, 0., 500. );
  HB1( 4204, "Xtgt Kurama [Proton]", 500, -200., 200. );
  HB1( 4205, "Ytgt Kurama [Proton]", 500, -100., 100. );
  HB1( 4206, "Utgt Kurama [Proton]", 700, -0.35, 0.35 );
  HB1( 4207, "Vtgt Kurama [Proton]", 400, -0.20, 0.20 );
  HB2( 4208, "Xtgt%U Kurama [Proton]", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 4209, "Ytgt%V Kurama [Proton]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 4210, "Xtgt%Ytgt Kurama [Proton]", 100, -200., 200., 100, -100., 100. );
  HB1( 4211, "P Kurama [Proton]", 200, 0.0, 1.0 );
  HB1( 4212, "PathLength Kurama [Proton]", 600, 3000., 6000. );
  HB1( 4213, "MassSqr", 600, -1.2, 1.2 );

  HB1( 5001, "Zvert [KK]", 1000, -1000., 1000. );
  HB1( 5002, "MissingMass [KK]", 1000, 0.0, 2.0 );

  HB2( 5011, "MissingMass%Us", 200, 0.0, 2.50, 100, -0.40, 0.40 );
  HB2( 5012, "MissingMass%Vs", 200, 0.0, 2.50, 100, -0.20, 0.20 );
  HB2( 5013, "MissingMass%Ub", 200, 0.0, 2.50, 100, -0.30, 0.30 );
  HB2( 5014, "MissingMass%Vb", 200, 0.0, 2.50, 100, -0.10, 0.10 );

  ////////////////////////////////////////////
  //Tree
  HBTree("kk","tree of KKAna");
  //Trigger
  tree->Branch("evnum",   &event.evnum,    "evnum/I");
  tree->Branch("trigpat",  event.trigpat,  Form("trigpat[%d]/I",NumOfSegTrig));
  tree->Branch("trigflag", event.trigflag, Form("trigflag[%d]/I",NumOfSegTrig));

  //Hodoscope
  tree->Branch("nhBh2",   &event.nhBh2,   "nhBh2/I");
  tree->Branch("Bh2Seg",   event.Bh2Seg,  "Bh2Seg[nhBh2]/D");
  tree->Branch("tBh2",     event.tBh2,    "tBh2[nhBh2]/D");
  tree->Branch("deBh2",    event.deBh2,   "deBh2[nhBh2]/D");

  tree->Branch("nhBh1",   &event.nhBh1,   "nhBh1/I");
  tree->Branch("Bh1Seg",   event.Bh1Seg,  "Bh1Seg[nhBh1]/D");
  tree->Branch("tBh1",     event.tBh1,    "tBh1[nhBh1]/D");
  tree->Branch("deBh1",    event.deBh1,   "deBh1[nhBh1]/D");
  tree->Branch("btof",     event.btof,    "btof[nhBh1]/D");

  tree->Branch("nhTof",   &event.nhTof,   "nhTof/I");
  tree->Branch("TofSeg",   event.TofSeg,  "TofSeg[nhTof]/D");
  tree->Branch("tTof",     event.tTof,    "tTof[nhTof]/D");
  tree->Branch("dtTof",    event.dtTof,   "dtTof[nhTof]/D");
  tree->Branch("deTof",    event.deTof,   "deTof[nhTof]/D");

  tree->Branch("nhBac",   &event.nhBac,   "nhBac/I");
  tree->Branch("BacSeg",   event.BacSeg,  "BacSeg[nhBac]/D");
  tree->Branch("tBac",     event.tBac,    "tBac[nhBac]/D");
  tree->Branch("deBac",    event.deBac,   "deBac[nhBac]/D");
  tree->Branch("nhPvac",   &event.nhPvac,   "nhPvac/I");
  tree->Branch("PvacSeg",   event.PvacSeg,  "PvacSeg[nhPvac]/D");
  tree->Branch("tPvac",     event.tPvac,    "tPvac[nhPvac]/D");
  tree->Branch("dePvac",    event.dePvac,   "dePvac[nhPvac]/D");
  tree->Branch("nhFac",   &event.nhFac,   "nhFac/I");
  tree->Branch("FacSeg",   event.FacSeg,  "FacSeg[nhFac]/D");
  tree->Branch("tFac",     event.tFac,    "tFac[nhFac]/D");
  tree->Branch("deFac",    event.deFac,   "deFac[nhFac]/D");

  //Fiber
  tree->Branch("bft_ncl",    &event.bft_ncl,    "bft_ncl/I");
  tree->Branch("bft_clsize",  event.bft_clsize, "bft_clsize[bft_ncl]/I");
  tree->Branch("bft_ctime",   event.bft_ctime,  "bft_ctime[bft_ncl]/D");
  tree->Branch("bft_ctot",    event.bft_ctot,   "bft_ctot[bft_ncl]/D");
  tree->Branch("bft_clpos",   event.bft_clpos,  "bft_clpos[bft_ncl]/D");
  tree->Branch("sch_ncl",    &event.sch_ncl,    "sch_ncl/I");
  tree->Branch("sch_clsize",  event.sch_clsize, "sch_clsize[sch_ncl]/I");
  tree->Branch("sch_ctime",   event.sch_ctime,  "sch_ctime[sch_ncl]/D");
  tree->Branch("sch_ctot",    event.sch_ctot,   "sch_ctot[sch_ncl]/D");
  tree->Branch("sch_clpos",   event.sch_clpos,  "sch_clpos[sch_ncl]/D");
  tree->Branch("fbh_ncl",    &event.fbh_ncl,    "fbh_ncl/I");
  tree->Branch("fbh_clsize",  event.fbh_clsize, "fbh_clsize[fbh_ncl]/I");
  tree->Branch("fbh_ctime",   event.fbh_ctime,  "fbh_ctime[fbh_ncl]/D");
  tree->Branch("fbh_ctot",    event.fbh_ctot,   "fbh_ctot[fbh_ncl]/D");
  tree->Branch("fbh_clpos",   event.fbh_clpos,  "fbh_clpos[fbh_ncl]/D");

  //Beam DC
  tree->Branch("nlBcOut",   &event.nlBcOut,     "nlBcOut/I");
  tree->Branch("ntBcOut",   &event.ntBcOut,     "ntBcOut/I");
  tree->Branch("nhBcOut",    event.nhBcOut,     "nhBcOut[ntBcOut]/I");
  tree->Branch("chisqrBcOut",event.chisqrBcOut, "chisqrBcOut[ntBcOut]/D");
  tree->Branch("x0BcOut",    event.x0BcOut,     "x0BcOut[ntBcOut]/D");
  tree->Branch("y0BcOut",    event.y0BcOut,     "y0BcOut[ntBcOut]/D");
  tree->Branch("u0BcOut",    event.u0BcOut,     "u0BcOut[ntBcOut]/D");
  tree->Branch("v0BcOut",    event.v0BcOut,     "v0BcOut[ntBcOut]/D");
  tree->Branch("xtgtBcOut",  event.xtgtBcOut,   "xtgtBcOut[ntBcOut]/D");
  tree->Branch("ytgtBcOut",  event.ytgtBcOut,   "ytgtBcOut[ntBcOut]/D");
  tree->Branch("xbh2BcOut",  event.xbh2BcOut,   "xbh2BcOut[ntBcOut]/D");
  tree->Branch("ybh2BcOut",  event.ybh2BcOut,   "ybh2BcOut[ntBcOut]/D");

  tree->Branch("ntK18",      &event.ntK18,     "ntK18/I");
  tree->Branch("nhK18",       event.nhK18,     "nhK18[ntK18]/I");
  tree->Branch("chisqrK18",   event.chisqrK18, "chisqrK18[ntK18]/D");
  tree->Branch("pK18",        event.pK18,      "pK18[ntK18]/D");
  tree->Branch("xtb",         event.xtb,       "xtb[ntK18]/D");
  tree->Branch("ytb",         event.ytb,       "ytb[ntK18]/D");
  tree->Branch("utb",         event.utb,       "utb[ntK18]/D");
  tree->Branch("vtb",         event.vtb,       "vtb[ntK18]/D");
  tree->Branch("thetab",      event.thetab,    "thetab[ntK18]/D");

  //KURAMA
  tree->Branch("nlSdcIn",    &event.nlSdcIn,     "nlSdcIn/I");
  tree->Branch("ntSdcIn",    &event.ntSdcIn,     "ntSdcIn/I");
  tree->Branch("nhSdcIn",     event.nhSdcIn,     "nhSdcIn[ntSdcIn]/I");
  tree->Branch("chisqrSdcIn", event.chisqrSdcIn, "chisqrSdcIn[ntSdcIn]/D");
  tree->Branch("x0SdcIn",     event.x0SdcIn,     "x0SdcIn[ntSdcIn]/D");
  tree->Branch("y0SdcIn",     event.y0SdcIn,     "y0SdcIn[ntSdcIn]/D");
  tree->Branch("u0SdcIn",     event.u0SdcIn,     "u0SdcIn[ntSdcIn]/D");
  tree->Branch("v0SdcIn",     event.v0SdcIn,     "v0SdcIn[ntSdcIn]/D");

  tree->Branch("nlSdcOut",   &event.nlSdcOut,     "nlSdcOut/I");
  tree->Branch("ntSdcOut",   &event.ntSdcOut,     "ntSdcOut/I");
  tree->Branch("nhSdcOut",    event.nhSdcOut,     "nhSdcOut[ntSdcOut]/I");
  tree->Branch("chisqrSdcOut",event.chisqrSdcOut, "chisqrOut[ntSdcOut]/D");
  tree->Branch("x0SdcOut",    event.x0SdcOut,     "x0SdcOut[ntSdcOut]/D");
  tree->Branch("y0SdcOut",    event.y0SdcOut,     "y0SdcOut[ntSdcOut]/D");
  tree->Branch("u0SdcOut",    event.u0SdcOut,     "u0SdcOut[ntSdcOut]/D");
  tree->Branch("v0SdcOut",    event.v0SdcOut,     "v0SdcOut[ntSdcOut]/D");

  tree->Branch("ntKurama",      &event.ntKurama,     "ntKurama/I");
  tree->Branch("nhKurama",       event.nhKurama,     "nhKurama[ntKurama]/I");
  tree->Branch("chisqrKurama",   event.chisqrKurama, "chisqrKurama[ntKurama]/D");
  tree->Branch("path",        event.path,      "path[ntKurama]/D");
  tree->Branch("pKurama",        event.pKurama,      "pKurama[ntKurama]/D");
  tree->Branch("polarity",       event.polarity,     "polarity[ntKurama]/D");
  tree->Branch("m2",          event.m2,        "m2[ntKurama]/D");
  tree->Branch("xts",         event.xts,       "xts[ntKurama]/D");
  tree->Branch("yts",         event.yts,       "yts[ntKurama]/D");
  tree->Branch("uts",         event.uts,       "uts[ntKurama]/D");
  tree->Branch("vts",         event.vts,       "vts[ntKurama]/D");
  tree->Branch("thetas",      event.thetas,    "thetas[ntKurama]/D");

#if 0
  tree->Branch("rm_event",   &event.rm_event,   "rm_event/I");
  tree->Branch("rm_spill",   &event.rm_spill,   "rm_spill/I");
  tree->Branch("rm_accept",  &event.rm_accept,  "rm_accept/I");
  tree->Branch("rm_clear",   &event.rm_clear,   "rm_clear/I");
  tree->Branch("mst_accept", &event.mst_accept, "mst_accept/I");
  tree->Branch("mst_tdc",     event.mst_tdc,    Form( "mst_tdc[%d]/I", NumOfSegTOF ) );
  tree->Branch("mst_tof",     event.mst_tof,    Form( "mst_tof[%d]/I", NumOfSegTOF ) );
  tree->Branch("mst_sch",     event.mst_sch,    Form( "mst_sch[%d]/I", NumOfSegSCH ) );
#endif

  //Reaction
  tree->Branch("nKn",           &event.nKn,            "nKn/I");
  tree->Branch("nKp",           &event.nKp,            "nKp/I");
  tree->Branch("nKK",           &event.nKK,            "nKK/I");
  tree->Branch("vtx",            event.vtx,            "vtx[nKK]/D");
  tree->Branch("vty",            event.vty,            "vty[nKK]/D");
  tree->Branch("vtz",            event.vtz,            "vtz[nKK]/D");
  tree->Branch("closeDist",      event.closeDist,      "closeDist[nKK]/D");
  tree->Branch("theta",          event.theta,          "theta[nKK]/D");
  tree->Branch("MissMass",       event.MissMass,       "MissMass[nKK]/D");
  tree->Branch("MissMassCorr",   event.MissMassCorr,   "MissMassCorr[nKK]/D");
  tree->Branch("MissMassCorrDE", event.MissMassCorrDE, "MissMassCorrDE[nKK]/D");
  tree->Branch("BE",      event.BE,       "BE[nKK]/D");
  tree->Branch("theta_CM", event.theta_CM, "theta_CM[nKK]/D");
  tree->Branch("cost_CM",  event.cost_CM,  "cost_CM[nKK]/D");

  tree->Branch("xkn",        event.xkn,      "xkn[nKK]/D");
  tree->Branch("ykn",        event.ykn,      "ykn[nKK]/D");
  tree->Branch("ukn",        event.ukn,      "ukn[nKK]/D");
  tree->Branch("vkn",        event.vkn,      "vkn[nKK]/D");
  tree->Branch("xkp",        event.xkp,      "xkp[nKK]/D");
  tree->Branch("ykp",        event.ykp,      "ykp[nKK]/D");
  tree->Branch("ukp",        event.ukp,      "ukp[nKK]/D");
  tree->Branch("vkp",        event.vkp,      "vkp[nKK]/D");
  tree->Branch("pOrg",       event.pOrg,      "pOrg[nKK]/D");
  tree->Branch("pCalc",      event.pCalc,     "pCalc[nKK]/D");
  tree->Branch("pCorr",      event.pCorr,     "pCorr[nKK]/D");
  tree->Branch("pCorrDE",    event.pCorrDE,   "pCorrDE[nKK]/D");

  HPrint();
  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")         &&
      InitializeParameter<DCDriftParamMan>("DCDRFT")  &&
      InitializeParameter<DCTdcCalibMan>("DCTDC")     &&
      InitializeParameter<HodoParamMan>("HDPRM")      &&
      InitializeParameter<HodoPHCMan>("HDPHC")        &&
      InitializeParameter<FieldMan>("FLDMAP")         &&
      InitializeParameter<K18TransMatrix>("K18TM")    &&
      InitializeParameter<MatrixParamMan>("MATRIX2D",
					  "MATRIX3D") &&
      InitializeParameter<MsTParamMan>("MASS")        &&
      InitializeParameter<UserParamMan>("USER")       );
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
