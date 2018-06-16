/**
 *  file: UserEventDisplay.cc
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
#include "DatabasePDG.hh"
#include "DetectorID.hh"
#include "EventDisplay.hh"
#include "RMAnalyzer.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "FLHit.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
//#include "RootHelper.hh"
#include "VEvent.hh"

namespace
{
  const std::string& class_name("EventDisplay");
  const DCGeomMan&     gGeom   = DCGeomMan::GetInstance();
  EventDisplay&        gEvDisp = EventDisplay::GetInstance();
  RMAnalyzer&          gRM     = RMAnalyzer::GetInstance();
  const UserParamMan&  gUser   = UserParamMan::GetInstance();
  const double KaonMass   = pdg::KaonMass();
  const double ProtonMass = pdg::ProtonMass();
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
class UserEventDisplay : public VEvent
{
private:
  RawData      *rawData;
  DCAnalyzer   *DCAna;
  HodoAnalyzer *hodoAna;
public:
        UserEventDisplay( void );
       ~UserEventDisplay( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
};

//______________________________________________________________________________
UserEventDisplay::UserEventDisplay( void )
  : VEvent(),
    rawData(0),
    DCAna( new DCAnalyzer ),
    hodoAna( new HodoAnalyzer )
{
}

//______________________________________________________________________________
UserEventDisplay::~UserEventDisplay( void )
{
  if (DCAna)   delete DCAna;
  if (hodoAna) delete hodoAna;
  if (rawData) delete rawData;
}

//______________________________________________________________________________
bool
UserEventDisplay::ProcessingBegin( void )
{
  return true;
}

//______________________________________________________________________________
bool
UserEventDisplay::ProcessingNormal( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"]");

  static const double MaxMultiHitBcOut  = gUser.GetParameter("MaxMultiHitBcOut");
  static const double MaxMultiHitSdcIn  = gUser.GetParameter("MaxMultiHitSdcIn");
  static const double MaxMultiHitSdcOut = gUser.GetParameter("MaxMultiHitSdcOut");

  static const double MinTimeBFT = gUser.GetParameter("TimeBFT", 0);
  static const double MaxTimeBFT = gUser.GetParameter("TimeBFT", 1);
  static const double MinTdcFBH  = gUser.GetParameter("TdcFBH", 0);
  static const double MaxTdcFBH  = gUser.GetParameter("TdcFBH", 1);
  static const double MinTdcSCH  = gUser.GetParameter("TdcSCH", 0);
  static const double MaxTdcSCH  = gUser.GetParameter("TdcSCH", 1);
  static const double MinDeSSDKaon = gUser.GetParameter("DeSSDKaon", 0);
  static const double MaxDeSSDKaon = gUser.GetParameter("DeSSDKaon", 1);
  static const double MinTimeSSD = gUser.GetParameter("TimeSSD", 0);
  static const double MaxTimeSSD = gUser.GetParameter("TimeSSD", 1);

  static const double MaxChisqrSSD = gUser.GetParameter("MaxChisqrSSD", 0);

  static const double OffsetToF  = gUser.GetParameter("OffsetToF");

  // static const int IdBH2 = gGeom.GetDetectorId("BH2");
  static const int IdFBH = gGeom.GetDetectorId("FBH");
  static const int IdSCH = gGeom.GetDetectorId("SCH");
  static const int IdTOF = gGeom.GetDetectorId("TOF");

  rawData = new RawData;
  rawData->DecodeHits();

  gRM.Decode();

  gEvDisp.DrawText( 0.1, 0.3, Form("Run# %5d%4sEvent# %6d",
				    gRM.RunNumber(), "",
				    gRM.EventNumber() ) );

  // Trigger Flag
  {
    const HodoRHitContainer &cont=rawData->GetTrigRawHC();
    int trigflag[NumOfSegTrig] = {};
    int nh = cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      int seg = hit->SegmentId()+1;
      int tdc = hit->GetTdc1();
      trigflag[seg-1] = tdc;
    }
    if( trigflag[SpillEndFlag]>0 ) return true;
  }

  //SSDT
  std::vector<double> t0Ssd( NumOfSegSSDT, 0. );
  {
    hodoAna->DecodeSSDTHits( rawData );
    int nh = hodoAna->GetNHitsSSDT();
    for( int i=0; i<nh; ++i ){
      Hodo1Hit *hit = hodoAna->GetHitSSDT(i);
      if(!hit) continue;
      int    seg  = hit->SegmentId()+1;
      double time = hit->Time();
      t0Ssd[seg-1] = time;
    }
  }

  // BH2
  // {
  //   const HodoRHitContainer &cont = rawData->GetBH2RawHC();
  //   int nh=cont.size();
  //   for( int i=0; i<nh; ++i ){
  //     HodoRawHit *hit = cont[i];
  //     if( !hit ) continue;
  //     int seg=hit->SegmentId();
  //     int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();
  //     // if( Tu>0 || Td>0 )
  //     // 	gEvDisp.DrawHitHodoscope( DCGeomIdBH2, seg, Tu, Td );
  //   }
  // }
  hodoAna->DecodeBH2Hits(rawData);
  int ncBh2 = hodoAna->GetNClustersBH2();
  if( ncBh2==0 ) return true;
  BH2Cluster *clBH2 = hodoAna->GetClusterBH2(0);
  double time0 = clBH2->CTime0();

  // TOF
  {
    const HodoRHitContainer &cont = rawData->GetTOFRawHC();
    int nh = cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      if( !hit ) continue;
      int seg = hit->SegmentId();
      int Tu = hit->GetTdcUp(), Td = hit->GetTdcDown();
      if( Tu>0 || Td>0 )
	gEvDisp.DrawHitHodoscope( IdTOF, seg, Tu, Td );
    }
  }
  hodoAna->DecodeTOFHits(rawData);
  Hodo2HitContainer TOFCont;
  int nhTof = hodoAna->GetNHitsTOF();
  for( int i=0; i<nhTof; ++i ){
    Hodo2Hit *hit = hodoAna->GetHitTOF(i);
    if( !hit ) continue;
    TOFCont.push_back( hit );
  }
  if( nhTof==0 ) return true;

  // FBH
  {
    hodoAna->DecodeFBHHits(rawData);
    int nhFbh = hodoAna->GetNHitsFBHCoin();
    for( int i=0; i<nhFbh; ++i ){
      FLHit *hit = hodoAna->GetHitFBHCoin(i);
      if( !hit ) continue;
      double seg = hit->SegmentId();
      bool   hit_flag = false;
      double leading = hit->GetLeading();
      if( MinTdcFBH<leading && leading<MaxTdcFBH ){
	hit_flag = true;
      }
      if( hit_flag ){
  	gEvDisp.DrawHitHodoscope( IdFBH, seg );
      }
    }
  }

  // SCH
  {
    hodoAna->DecodeSCHHits(rawData);
    int nhSch = hodoAna->GetNHitsSCH();
    for( int i=0; i<nhSch; ++i ){
      FiberHit *hit = hodoAna->GetHitSCH(i);
      if( !hit ) continue;
      int mh  = hit->GetNumOfHit();
      int seg = hit->SegmentId();
      bool hit_flag = false;
      for( int m=0; m<mh; ++m ){
  	double leading = hit->GetLeading(m);
  	if( MinTdcSCH<leading && leading<MaxTdcSCH ){
  	  hit_flag = true;
  	}
      }
      if( hit_flag ){
  	gEvDisp.DrawHitHodoscope( IdSCH, seg );
      }
    }
  }

  DCAna->DecodeRawHits( rawData );

  // BcOut
  double multi_BcOut = 0.;
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &contIn =DCAna->GetBcOutHC(layer);
      int nhIn=contIn.size();
      for( int i=0; i<nhIn; ++i ){
	// DCHit  *hit  = contIn[i];
	// double  wire = hit->GetWire();
	// int     mhit = hit->GetTdcSize();
	// bool    goodFlag = false;
	++multi_BcOut;
	// for (int j=0; j<mhit; j++) {
	//   if (hit->IsWithinRange(j)) {
	//     goodFlag = true;
	//     break;
	//   }
	// }
	// if( goodFlag )
	//   gEvDisp.DrawHitWire( layer+112, int(wire) );
      }
    }
  }
  multi_BcOut /= (double)NumOfLayersBcOut;
  if( multi_BcOut > MaxMultiHitBcOut ) return true;

  // SdcIn
  double multi_SdcIn = 0.;
  {
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetSdcInHC(layer);
      int nhIn=contIn.size();
      if ( nhIn > MaxMultiHitSdcIn )
	continue;
      for( int i=0; i<nhIn; ++i ){
	DCHit  *hit  = contIn[i];
	double  wire = hit->GetWire();
	int     mhit = hit->GetTdcSize();
	bool    goodFlag = false;
	++multi_SdcIn;
	for (int j=0; j<mhit; j++) {
	  if (hit->IsWithinRange(j)) {
	    goodFlag = true;
	    break;
	  }
	}
	if( goodFlag )
	  gEvDisp.DrawHitWire( layer, int(wire) );
      }
    }
  }
  multi_SdcIn /= (double)NumOfLayersSdcIn;
  if( multi_SdcIn > MaxMultiHitSdcIn ) return true;

  // SdcOut
  double multi_SdcOut = 0.;
  {
    for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetSdcOutHC(layer);
      int nhOut = contOut.size();
      if ( nhOut > MaxMultiHitSdcOut )
	continue;
      for( int i=0; i<nhOut; ++i ){
	DCHit  *hit  = contOut[i];
	double  wire = hit->GetWire();
	++multi_SdcOut;
	gEvDisp.DrawHitWire( layer+30, int(wire) );
      }
    }
  }
  multi_SdcOut /= (double)NumOfLayersSdcOut;
  if( multi_SdcOut > MaxMultiHitSdcOut ) return true;

  int ntBcOut = 0;
  if( multi_BcOut<MaxMultiHitBcOut ){
    DCAna->TrackSearchBcOut();
    ntBcOut = DCAna->GetNtracksBcOut();
    for( int it=0; it<ntBcOut; ++it ){
      DCLocalTrack *tp = DCAna->GetTrackBcOut( it );
      if( tp ) gEvDisp.DrawBcOutLocalTrack( tp );
    }
  }
  if( ntBcOut==0 ) return true;

  // SsdAnalysis
  DCAna->DoTimeCorrectionSsd( t0Ssd );
  DCAna->ChisqrFilterSsd( MaxChisqrSSD );
  DCAna->ClusterizeSsd();
  // for SsdCluster
  DCAna->DeltaEFilterSsd( MinDeSSDKaon, MaxDeSSDKaon, true );
  DCAna->TimeFilterSsd( MinTimeSSD, MaxTimeSSD, true );
  // SsdIn
  for( int l=1; l<=NumOfLayersSsdIn; ++l ){
    const SsdClusterContainer &cont = DCAna->GetClusterSsdIn(l);
    const std::size_t nh = cont.size();
    for( std::size_t i=0; i<nh; ++i ){
      SsdCluster *cluster = cont[i];
      if(!cluster) continue;
      int    layer = cluster->LayerId();
      double seg   = cluster->MeanSeg();
      double de    = cluster->DeltaE();
      gEvDisp.DrawSsdHit( layer, seg, de );
    }
  }

  // SsdOut
  for( int l=1; l<=NumOfLayersSsdOut; ++l ){
    const SsdClusterContainer &cont = DCAna->GetClusterSsdOut(l);
    const std::size_t nh = cont.size();
    for( std::size_t i=0; i<nh; ++i ){
      SsdCluster *cluster = cont[i];
      if(!cluster) continue;
      int    layer = cluster->LayerId();
      double seg   = cluster->MeanSeg();
      double de    = cluster->DeltaE();
      gEvDisp.DrawSsdHit( layer, seg, de );
    }
  }

  int ntSdcIn = 0;
  if( multi_SdcIn<MaxMultiHitSdcIn ){
    DCAna->TrackSearchSdcIn();
    ntSdcIn = DCAna->GetNtracksSdcIn();
    for( int it=0; it<ntSdcIn; ++it ){
      DCLocalTrack *tp = DCAna->GetTrackSdcIn( it );
      if( tp ) gEvDisp.DrawSdcInLocalTrack( tp );
    }
  }
  if( ntSdcIn==0 ) return true;

  int ntSdcOut = 0;
  if( multi_SdcOut<MaxMultiHitSdcOut ){
    DCAna->TrackSearchSdcOut( TOFCont );
    ntSdcOut = DCAna->GetNtracksSdcOut();
    for( int it=0; it<ntSdcOut; ++it ){
      DCLocalTrack *tp = DCAna->GetTrackSdcOut( it );
      if( tp ) gEvDisp.DrawSdcOutLocalTrack( tp );
    }
  }
  if( ntSdcOut==0 ) return true;

  std::vector<ThreeVector> KnPCont, KnXCont;
  std::vector<ThreeVector> KpPCont, KpXCont;

  int ntKurama = 0;
  static int ntKurama_all = 0;
  if( ntSdcIn>0 && ntSdcOut>0 ){
    //if( ntSdcIn==1 && ntSdcOut==1 ){
    bool through_target = false;
    DCAna->TrackSearchKurama();
    ntKurama = DCAna->GetNTracksKurama();
    ntKurama_all++;
    for( int it=0; it<ntKurama; ++it ){
      KuramaTrack *tp = DCAna->GetKuramaTrack( it );
      if( !tp ) continue;
      tp->Print( "in "+func_name );
      const ThreeVector& postgt = tp->PrimaryPosition();
      const ThreeVector& momtgt = tp->PrimaryMomentum();
      double path = tp->PathLengthToTOF();
      double p    = momtgt.Mag();
      gEvDisp.DrawMomentum( p );
      if( std::abs( postgt.x() )<50. &&
	  std::abs( postgt.y() )<30. ){
	through_target = true;
      }
      // MassSquare
      double tofseg = tp->TofSeg();
      for( int j=0, n=TOFCont.size(); j<n; ++j ){
      	Hodo2Hit *hit = TOFCont[j];
      	if( !hit ) continue;
      	double seg = hit->SegmentId()+1;
	if( tofseg != seg ) continue;
      	double stof = hit->CMeanTime()-time0+OffsetToF;
      	if( stof<=0 ) continue;
      	double m2 = Kinematics::MassSquare( p, path, stof );
      	gEvDisp.DrawMassSquare( m2 );
	KpPCont.push_back( momtgt );
	KpXCont.push_back( postgt );
      }
    }
    if( through_target ) gEvDisp.DrawTarget();
    static double KuramaOk = 0.;
    KuramaOk += ( ntKurama>0 );
  }

  std::vector<double> BftXCont;
  ////////// BFT
  {
    hodoAna->DecodeBFTHits(rawData);
    // Fiber Cluster
    hodoAna->TimeCutBFT(MinTimeBFT, MaxTimeBFT);
    int ncl = hodoAna->GetNClustersBFT();
    for( int i=0; i<ncl; ++i ){
      FiberCluster *cl = hodoAna->GetClusterBFT(i);
      if( !cl ) continue;
      double pos    = cl->MeanPosition();
      BftXCont.push_back( pos );
    }
  }

  // K18TrackingD2U
  DCAna->TrackSearchK18D2U( BftXCont );
  int ntK18=DCAna->GetNTracksK18D2U();
  if( ntK18==0 ) return true;
  for( int i=0; i<ntK18; ++i ){
    K18TrackD2U *tp=DCAna->GetK18TrackD2U(i);
    if(!tp) continue;
    double x = tp->Xtgt(), y = tp->Ytgt();
    double u = tp->Utgt(), v = tp->Vtgt();
    double p = tp->P3rd();
    double pt = p/std::sqrt(1.+u*u+v*v);
    ThreeVector Pos( x, y, 0. );
    ThreeVector Mom( pt*u, pt*v, pt );
    KnPCont.push_back( Mom );
    KnXCont.push_back( Pos );
  }

#if 1
  if( KnPCont.size()==1 && KpPCont.size()==1 ){
    ThreeVector pkp = KpPCont[0];
    ThreeVector pkn = KnPCont[0];
    ThreeVector xkp = KpXCont[0];
    ThreeVector xkn = KnXCont[0];
    ThreeVector vertex = Kinematics::VertexPoint( xkn, xkp, pkn, pkp );
    LorentzVector LvKn( KnPCont[0], std::sqrt( KaonMass*KaonMass+pkn.Mag2()) );
    LorentzVector LvKp( KpPCont[0], std::sqrt( KaonMass*KaonMass+pkp.Mag2() ) );
    LorentzVector LvP( 0., 0., 0., ProtonMass );
    LorentzVector LvRp = LvKn+LvP-LvKp;
    ThreeVector MissMom = LvRp.Vect();
    gEvDisp.DrawVertex( vertex );
    gEvDisp.DrawMissingMomentum( MissMom, vertex );
  }
#endif

  return true;
}

//______________________________________________________________________________
bool
UserEventDisplay::ProcessingEnd( void )
{
  // gEvDisp.GetCommand();
  gEvDisp.EndOfEvent();
  // if( utility::UserStop() ) gEvDisp.Run();
  return true;
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new UserEventDisplay;
}

//______________________________________________________________________________
bool
ConfMan:: InitializeHistograms( void )
{
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
      InitializeParameter<FieldMan>("FLDMAP")        &&
      InitializeParameter<K18TransMatrix>("K18TM")   &&
      InitializeParameter<SsdParamMan>("SSDPRM")     &&
      InitializeParameter<UserParamMan>("USER")      &&
      InitializeParameter<EventDisplay>()            );
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
