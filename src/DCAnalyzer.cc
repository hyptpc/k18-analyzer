/**
 *  file: DCAnalyzer.cc
 *  date: 2017.04.10
 *
 */

#include "DCAnalyzer.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include "ConfMan.hh"
#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DCLocalTrack.hh"
#include "DCRawHit.hh"
#include "DCTrackSearch.hh"
#include "DebugCounter.hh"
#include "DebugTimer.hh"
#include "FiberCluster.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "K18Parameters.hh"
//#include "K18TrackU2D.hh"
#include "K18TrackD2U.hh"
#include "KuramaTrack.hh"
#include "MathTools.hh"
#include "MWPCCluster.hh"
#include "RawData.hh"
#include "SsdCluster.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"

#define DefStatic
#include "DCParameters.hh"
#undef DefStatic

// Tracking routine selection __________________________________________________
/* BcInTracking */
#define UseBcIn    0 // not supported
/* BcOutTracking */
#define BcOut_XUV  0 // XUV Tracking (slow but accerate)
#define BcOut_Pair 1 // Pair plane Tracking (fast but bad for large angle track)
/* SdcInTracking */
#define UseSsdCluster     1 // use SSD Cluster
#define SdcIn_XUV         0 // XUV Tracking (not used in KURAMA)
#define SdcIn_Pair        1 // Pair plane Tracking (fast but bad for large angle track)
#define SdcIn_SsdPreTrack 0 // SsdPreTracking for too many combinations
#define SdcIn_Deletion    1 // Deletion method for too many combinations

// SsdSlopeFilter  _____________________________________________________________
#define SlopeFilter_Tight 0 // 0->1->2->3 : Up   4->5->6->7 : Down
#define SlopeFilter_Wide  1 // 0->1->2    : Up      5->6->7 : Down

namespace
{
  using namespace K18Parameter;
  const std::string& class_name("DCAnalyzer");
  const ConfMan&      gConf = ConfMan::GetInstance();
  const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();

  //______________________________________________________________________________
  const double& pK18 = ConfMan::Get<double>("PK18");
  const int& IdTOFUX = gGeom.DetectorId("TOF-UX");
  const int& IdTOFUY = gGeom.DetectorId("TOF-UY");
  const int& IdTOFDX = gGeom.DetectorId("TOF-DX");
  const int& IdTOFDY = gGeom.DetectorId("TOF-DY");

  const double MaxChiSqrKuramaTrack = 10000.;
  const double MaxTimeDifMWPC       =   100.;

  const double kMWPCClusteringWireExtension =  1.0; // [mm]
  const double kMWPCClusteringTimeExtension = 10.0; // [nsec]

  //______________________________________________________________________________
  inline bool /* for MWPCCluster */
  isConnectable( double wire1, double leading1, double trailing1,
		 double wire2, double leading2, double trailing2,
		 double wExt,  double tExt )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");
    double w1Min = wire1 - wExt;
    double w1Max = wire1 + wExt;
    double t1Min = leading1  - tExt;
    double t1Max = trailing1 + tExt;
    double w2Min = wire2 - wExt;
    double w2Max = wire2 + wExt;
    double t2Min = leading2  - tExt;
    double t2Max = trailing2 + tExt;
    bool isWireOk = !(w1Min>w2Max || w1Max<w2Min);
    bool isTimeOk = !(t1Min>t2Max || t1Max<t2Min);
#if 0
    hddaq::cout << func_name << std::endl
		<< " w1 = " << wire1
		<< " le1 = " << leading1
		<< " tr1 = " << trailing1 << "\n"
		<< " w2 = " << wire2
		<< " le2 = " << leading2
		<< " tr2 = " << trailing2 << "\n"
		<< " w1(" << w1Min << " -- " << w1Max << "), t1("
		<< t1Min << " -- " << t1Max << ")\n"
		<< " w2(" << w2Min << " -- " << w2Max << "), t2("
		<< t2Min << " -- " << t2Max << ")\n"
		<< " wire : " << isWireOk
		<< ", time : " << isTimeOk
		<< std::endl;
#endif
    return ( isWireOk && isTimeOk );
  }

  //______________________________________________________________________________
  inline void
  printConnectionFlag( const std::vector<std::deque<bool> >& flag )
  {
    for( std::size_t i=0, n=flag.size(); i<n; ++i ){
      hddaq::cout << std::endl;
      for( std::size_t j=0, m=flag[i].size(); j<m; ++j ){
	hddaq::cout << " " << flag[i][j];
      }
    }
    hddaq::cout << std::endl;
  }

  //______________________________________________________________________________
  inline bool /* for SsdCluster */
  IsClusterable( const DCHitContainer& HitCont, DCHit* CandHit )
  {
    if( !CandHit )
      return false;
    if( !CandHit->IsGoodWaveForm() )
      return false;

    for( std::size_t i=0, n=HitCont.size(); i<n; ++i ){
      DCHit* hit      = HitCont[i];
      int    wire     = hit->GetWire();
      int    CandWire = CandHit->GetWire();
      // double time     = hit->GetPeakTime();
      // double CandTime = CandHit->GetPeakTime();
      if( std::abs( wire-CandWire )==1
	  // && std::abs( time-CandTime )<25.
	  ){
	return true;
      }
    }
    return false;
  }
}

//______________________________________________________________________________
DCAnalyzer::DCAnalyzer( void )
  : m_is_decoded(n_type),
    m_much_combi(n_type),
    m_MWPCClCont(NumOfLayersBcIn+1),
    m_TempBcInHC(NumOfLayersBcIn+1),
    m_BcInHC(NumOfLayersBcIn+1),
    m_BcOutHC(NumOfLayersBcOut+1),
    m_SdcInHC(NumOfLayersSdcIn+1),
    m_SdcOutHC(NumOfLayersSdcOut+1),
    m_SsdInHC(NumOfLayersSsdIn+1),
    m_SsdOutHC(NumOfLayersSsdOut+1),
    m_SsdInClCont(NumOfLayersSsdIn+1),
    m_SsdOutClCont(NumOfLayersSsdOut+1),
    m_SdcInExTC(NumOfLayersSdcIn+1),
    m_SdcOutExTC(NumOfLayersSdcOut+1)
{
  for( int i=0; i<n_type; ++i ){
    m_is_decoded[i] = false;
    m_much_combi[i] = 0;
  }
  debug::ObjectCounter::increase(class_name);
}

DCAnalyzer::~DCAnalyzer( void )
{
  ClearKuramaTracks();
#if UseBcIn
  ClearK18TracksU2D();
  ClearTracksBcIn();
#endif
  ClearK18TracksD2U();
  ClearTracksSsdOut();
  ClearTracksSsdIn();
  ClearTracksSsdXY();
  ClearTracksSdcOut();
  ClearTracksSdcIn();
  ClearTracksBcOut();
  ClearTracksBcOutSdcIn();
  ClearTracksBcOutSsdIn();
  ClearTracksSsdOutSdcIn();
  ClearTracksSdcInSdcOut();
  ClearDCHits();
  ClearVtxHits();
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void
DCAnalyzer::PrintKurama( const std::string& arg ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  int nn = m_KuramaTC.size();
  hddaq::cout << func_name << " " << arg << std::endl
	      << "   KuramaTC.size : " << nn << std::endl;
  for( int i=0; i<nn; ++i ){
    KuramaTrack *tp = m_KuramaTC[i];
    hddaq::cout << std::setw(3) << i
		<< " Niter=" << std::setw(3) << tp->Niteration()
		<< " ChiSqr=" << tp->chisqr()
		<< " P=" << tp->PrimaryMomentum().Mag()
		<< " PL(TOF)=" << tp->PathLengthToTOF()
		<< std::endl;
  }
}

//______________________________________________________________________________
#if UseBcIn
bool
DCAnalyzer::DecodeBcInHits( RawData *rawData )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_decoded[k_BcIn] ){
    hddaq::cout << "#D " << func_name << " "
		<< "already decoded" << std::endl;
    return true;
  }

  ClearBcInHits();

  for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
    const DCRHitContainer &RHitCont = rawData->GetBcInRawHC(layer);
    int nh = RHitCont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit  = RHitCont[i];
      DCHit    *thit  = new DCHit( rhit->PlaneId()+PlOffsBc, rhit->WireId() );
      int       nhtdc = rhit->GetTdcSize();
      if(!thit) continue;
      for( int j=0; j<nhtdc; ++j ){
	thit->SetTdcVal( rhit->GetTdc(j) );
	thit->SetTdcTrailing( rhit->GetTrailing(j) );
      }

      if(thit->CalcMWPCObservables())
	m_TempBcInHC[layer].push_back(thit);
      else
	delete thit;
    }

    // hddaq::cout<<"*************************************"<<std::endl;
    int ncl = clusterizeMWPCHit( m_TempBcInHC[layer], m_MWPCClCont[layer]);
    // hddaq::cout<<"numCl="<< ncl << std::endl;
     for( int i=0; i<ncl; ++i ){
      MWPCCluster *p = m_MWPCClCont[layer][i];
      if(!p) continue;

      const MWPCCluster::Statistics& mean  = p->GetMean();
      const MWPCCluster::Statistics& first = p->GetFirst();
      double mwire    = mean.m_wire;
      double mwirepos = mean.m_wpos;
      double mtime    = mean.m_leading;
      double mtrail   = mean.m_trailing;

      DCHit *hit = new DCHit( layer+PlOffsBc, mwire );
      if(!hit) continue;
      hit->SetClusterSize( p->GetClusterSize() );
      hit->SetMWPCFlag( true );
      hit->SetWire( mwire );
      hit->SetMeanWire( mwire );
      hit->SetMeanWirePosition( mwirepos );
      hit->SetTrailing( mtrail );
      hit->SetDummyPair();
      hit->SetTdcVal( 0 );
      hit->SetTdcTrailing( 0 );

      if( hit->CalcMWPCObservables() )
	m_BcInHC[layer].push_back(hit);
      else
	delete hit;
    }
    // hddaq::cout << "nh="<< m_BcInHC[layer].size() <<std::endl;
  }

  m_is_decoded[k_BcIn] = true;
  return true;
}
#endif

//______________________________________________________________________________
bool
DCAnalyzer::DecodeBcOutHits( RawData *rawData )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_decoded[k_BcOut] ){
    hddaq::cout << "#D " << func_name << " "
		<< "already decoded" << std::endl;
    return true;
  }

  ClearBcOutHits();

  for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
    const DCRHitContainer &RHitCont=rawData->GetBcOutRawHC(layer);
    int nh = RHitCont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit  = RHitCont[i];
      DCHit    *hit   = new DCHit( rhit->PlaneId()+PlOffsBc, rhit->WireId() );
      int       nhtdc = rhit->GetTdcSize();
      if(!hit) continue;
      for( int j=0; j<nhtdc; ++j ){
	hit->SetTdcVal( rhit->GetTdc(j) );
      }

      if( hit->CalcDCObservables() )
	m_BcOutHC[layer].push_back(hit);
      else
	delete hit;
    }
  }

  m_is_decoded[k_BcOut] = true;
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::DecodeSdcInHits( RawData *rawData )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_decoded[k_SdcIn] ){
    hddaq::cout << "#D " << func_name << " "
		<< "already decoded" << std::endl;
    return true;
  }

  ClearSdcInHits();

  // SFT
  {
    // HodoAnalyzer& hodoAna = HodoAnalyzer::GetInstance();
    HodoAnalyzer hodoAna;
    hodoAna.DecodeSFTHits( rawData );
    for ( int l = 0; l < NumOfLayersSFT; ++l ) {

      int layerId = l + PlMinSdcIn + NumOfLayersSDC1;

      //      hodoAna.TimeCutSFT( l, -15, 5 );
      int ncl = hodoAna.GetNClustersSFT( l );

      for ( int j = 0; j < ncl; ++j ) {

  	FiberCluster* cl = hodoAna.GetClusterSFT( l, j );
	double seg  = cl->MeanSeg();
  	double pos  = cl->MeanPosition();
  	double time = cl->CMeanTime();

  	// DCHit *hit = new DCHit( l + PlMinSdcIn );
  	DCHit *hit = new DCHit( layerId, seg );
  	hit->SetTdcVal( static_cast<int>( time ) );

  	if ( hit->CalcFiberObservables() ) {
  	  hit->SetWirePosition( pos );
  	  m_SdcInHC[layerId].push_back( hit );
	  // m_SdcInHC[layer - 1].push_back(hit);
  	} else {
  	  delete hit;
  	}

      }
    }
  }

  // SDC1
  {
    for( int layer=1; layer<=NumOfLayersSdcIn - NumOfLayersSFT; ++layer ){
      const DCRHitContainer &RHitCont=rawData->GetSdcInRawHC(layer);
      int nh = RHitCont.size();
      for( int i=0; i<nh; ++i ){
	DCRawHit *rhit  = RHitCont[i];
	DCHit    *hit   = new DCHit( rhit->PlaneId(), rhit->WireId() );
	// DCHit    *hit   = new DCHit( rhit->PlaneId() + NumOfLayersSFT, rhit->WireId() );
	int       nhtdc = rhit->GetTdcSize();
	if(!hit) continue;
	for( int j=0; j<nhtdc; ++j ){
	  hit->SetTdcVal( rhit->GetTdc(j) );
	}
	if( hit->CalcDCObservables() ) {
	  m_SdcInHC[layer].push_back(hit);
	} else {
	  delete hit;
	}
      }
    }
  }
  m_is_decoded[k_SdcIn] = true;
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::DecodeSdcOutHits( RawData *rawData )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_decoded[k_SdcOut] ){
    hddaq::cout << "#D " << func_name << " "
		<< "already decoded" << std::endl;
    return true;
  }

  ClearSdcOutHits();

  for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
    const DCRHitContainer &RHitCont = rawData->GetSdcOutRawHC(layer);
    int nh = RHitCont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit  = RHitCont[i];
      DCHit    *hit   = new DCHit( rhit->PlaneId(), rhit->WireId() );
      int       nhtdc      = rhit->GetTdcSize();
      int       nhtrailing = rhit->GetTrailingSize();
      if(!hit) continue;
      for( int j=0; j<nhtdc; ++j ){
	hit->SetTdcVal( rhit->GetTdc(j) );
      }
      for( int j=0; j<nhtrailing; ++j ){
	hit->SetTdcTrailing( rhit->GetTrailing(j) );
      }
      if( hit->CalcDCObservables() )
	m_SdcOutHC[layer].push_back(hit);
      else
	delete hit;
    }
  }

  m_is_decoded[k_SdcOut] = true;
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::DecodeSsdInHits( RawData *rawData )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_decoded[k_SsdIn] ){
    hddaq::cout << "#D " << func_name << " "
		<< "already decoded" << std::endl;
    return true;
  }

  ClearSsdInHits();

  for( int layer=1; layer<NumOfLayersSsdIn+1;++layer){
    const DCRHitContainer &RHitCont = rawData->GetSsdInRawHC(layer);
    int nh = RHitCont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit = RHitCont[i];
      DCHit    *hit  = new DCHit( rhit->PlaneId()+PlOffsSsd, rhit->WireId() );
      if(!hit) continue;
      // zero suppression flag
      if( rhit->GetTrailingSize()>0 )
      	hit->SetTdcTrailing( rhit->GetTrailing(0) );

      int nhadc = rhit->GetTdcSize();
      int adc[nhadc];
      int tdc[nhadc];
      for( int j=0; j<nhadc; ++j ){
	adc[j] = rhit->GetTdc(j);
	tdc[j] = j+1;
	hit->SetAdcVal( adc[j] );
	hit->SetTdcVal( tdc[j] );
      }
      if( hit->CalcSsdObservables() )
	m_SsdInHC[layer].push_back(hit);
      else
	delete hit;
    }
    ClusterizeSsd( m_SsdInHC[layer], m_SsdInClCont[layer] );
  }

  m_is_decoded[k_SsdIn] = true;
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::DecodeSsdOutHits( RawData *rawData )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_decoded[k_SsdOut] ){
    hddaq::cout << "#D " << func_name << " "
		<< "already decoded" << std::endl;
    return true;
  }

  ClearSsdOutHits();

  for( int layer=1; layer<NumOfLayersSsdOut+1;++layer){
    const DCRHitContainer &RHitCont = rawData->GetSsdOutRawHC(layer);
    int nh = RHitCont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit = RHitCont[i];
      DCHit    *hit  = new DCHit(rhit->PlaneId()+PlOffsSsd,
				 rhit->WireId());
      if(!hit) continue;
      // zero suppression flag
      if( rhit->GetTrailingSize()>0 )
	hit->SetTdcTrailing( rhit->GetTrailing(0) );

      int nhadc = rhit->GetTdcSize();
      int adc[nhadc];
      int tdc[nhadc];
      for( int j=0; j<nhadc; ++j ){
	adc[j] = rhit->GetTdc(j);
	tdc[j] = j+1;
	hit->SetAdcVal( adc[j] );
	hit->SetTdcVal( tdc[j] );
      }
      if( hit->CalcSsdObservables() )
	m_SsdOutHC[layer].push_back(hit);
      else
	delete hit;
    }
    ClusterizeSsd( m_SsdOutHC[layer], m_SsdOutClCont[layer] );
  }

  m_is_decoded[k_SsdOut] = true;
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::DecodeRawHits( RawData *rawData )
{
  ClearDCHits();
#if UseBcIn
  DecodeBcInHits( rawData );
#endif
  DecodeBcOutHits( rawData );
  DecodeSdcInHits( rawData );
  DecodeSdcOutHits( rawData );
  // DecodeSsdInHits( rawData );
  // DecodeSsdOutHits( rawData );
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::DecodeTOFHits( const Hodo2HitContainer& HitCont )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_decoded[k_TOF] ){
    hddaq::cout << "#D " << func_name << " "
		<< "already decoded" << std::endl;
    return true;
  }

  ClearTOFHits();

  static const double RA2 = gGeom.GetRotAngle2("TOF"); // for the tilting plane case
  static const ThreeVector TOFPos[2] = {
    gGeom.GetGlobalPosition("TOF-UX"),
    gGeom.GetGlobalPosition("TOF-DX") };

  const std::size_t nh = HitCont.size();
  for( std::size_t i=0; i<nh; ++i ){
    Hodo2Hit *hodo_hit = HitCont[i];
    if( !hodo_hit ) continue;
    const double seg  = hodo_hit->SegmentId()+1;
    const double dt   = hodo_hit->TimeDiff();
    int layer_x = -1;
    int layer_y = -1;
    if( (int)seg%2==0 ){
      layer_x  = IdTOFUX;
      layer_y  = IdTOFUY;
    }
    if( (int)seg%2==1 ){
      layer_x  = IdTOFDX;
      layer_y  = IdTOFDY;
    }
    double wpos = gGeom.CalcWirePosition( layer_x, seg );
    ThreeVector w( wpos, 0., 0. );
    w.RotateY( RA2*math::Deg2Rad() ); // for the tilting plane case
    const ThreeVector& hit_pos = TOFPos[(int)seg%2] + w
      + ThreeVector( 0., dt/0.01285, 0. );
    // X
    DCHit *dc_hit_x = new DCHit( layer_x, seg );
    dc_hit_x->SetWirePosition( hit_pos.x() );
    dc_hit_x->SetZ( hit_pos.z() );
    dc_hit_x->SetTiltAngle( 0. );
    dc_hit_x->SetDummyPair();
    m_TOFHC.push_back( dc_hit_x );
    // Y
    DCHit *dc_hit_y = new DCHit( layer_y, seg );
    dc_hit_y->SetWirePosition( hit_pos.y() ); // [ns] -> [mm]
    dc_hit_y->SetZ( hit_pos.z() );
    dc_hit_y->SetTiltAngle( 90. );
    dc_hit_y->SetDummyPair();
    m_TOFHC.push_back( dc_hit_y );
  }

  m_is_decoded[k_TOF] = true;
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::DecodeTOFHits( const HodoClusterContainer& ClCont )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_decoded[k_TOF] ){
    hddaq::cout << "#D " << func_name << " "
		<< "already decoded" << std::endl;
    return true;
  }

  ClearTOFHits();

  static const double RA2 = gGeom.GetRotAngle2("TOF");
  static const ThreeVector TOFPos[2] = {
    gGeom.GetGlobalPosition("TOF-UX"),
    gGeom.GetGlobalPosition("TOF-DX") };

  const std::size_t nh = ClCont.size();
  for( std::size_t i=0; i<nh; ++i ){
    HodoCluster *hodo_cluster = ClCont[i];
    if( !hodo_cluster ) continue;
    const double seg  = hodo_cluster->MeanSeg()+1;
    const double dt   = hodo_cluster->TimeDif();
    int layer_x = -1;
    int layer_y = -1;
    if( (int)seg%2==0 ){
      layer_x  = IdTOFUX;
      layer_y  = IdTOFUY;
    }
    if( (int)seg%2==1 ){
      layer_x  = IdTOFDX;
      layer_y  = IdTOFDY;
    }
    double wpos = gGeom.CalcWirePosition( layer_x, seg );
    ThreeVector w( wpos, 0., 0. );
    w.RotateY( RA2*math::Deg2Rad() );
    const ThreeVector& hit_pos = TOFPos[(int)seg%2] + w
      + ThreeVector( 0., dt/0.01285, 0. );
    // X
    DCHit *dc_hit_x = new DCHit( layer_x, seg );
    dc_hit_x->SetWirePosition( hit_pos.x() );
    dc_hit_x->SetZ( hit_pos.z() );
    dc_hit_x->SetTiltAngle( 0. );
    dc_hit_x->SetDummyPair();
    m_TOFHC.push_back( dc_hit_x );
    // Y
    DCHit *dc_hit_y = new DCHit( layer_y, seg );
    dc_hit_y->SetWirePosition( hit_pos.y() ); // [ns] -> [mm]
    dc_hit_y->SetZ( hit_pos.z() );
    dc_hit_y->SetTiltAngle( 90. );
    dc_hit_y->SetDummyPair();
    m_TOFHC.push_back( dc_hit_y );
  }

  m_is_decoded[k_TOF] = true;
  return true;
}

//______________________________________________________________________________
#if UseBcIn
bool
DCAnalyzer::TrackSearchBcIn( void )
{
  track::MWPCLocalTrackSearch( &(m_BcInHC[1]), m_BcInTC );
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::TrackSearchBcIn( const std::vector<std::vector<DCHitContainer> >& hc )
{
  track::MWPCLocalTrackSearch( hc, m_BcInTC );
  return true;
}
#endif

//______________________________________________________________________________
bool
DCAnalyzer::TrackSearchBcOut( void )
{
  static const int MinLayer = gUser.GetParameter("MinLayerBcOut");

#if BcOut_Pair //Pair Plane Tracking Routine for BcOut
  track::LocalTrackSearch( m_BcOutHC, PPInfoBcOut, NPPInfoBcOut,
			   m_BcOutTC, MinLayer );
  return true;
#endif

#if BcOut_XUV  //XUV Tracking Routine for BcOut
  track::LocalTrackSearchVUX( m_BcOutHC, PPInfoBcOut, NPPInfoBcOut,
			      m_BcOutTC, MinLayer );
  return true;
#endif

  return false;
}

//______________________________________________________________________________
// Use with BH2Filter
bool
DCAnalyzer::TrackSearchBcOut( const std::vector<std::vector<DCHitContainer> >& hc )
{
  static const int MinLayer = gUser.GetParameter("MinLayerBcOut");

#if BcOut_Pair //Pair Plane Tracking Routine for BcOut
  track::LocalTrackSearch( hc, PPInfoBcOut, NPPInfoBcOut, m_BcOutTC, MinLayer );
  return true;
#endif

#if BcOut_XUV  //XUV Tracking Routine for BcOut
  track::LocalTrackSearchVUX( hc, PPInfoBcOut, NPPInfoBcOut, m_BcOutTC, MinLayer );
  return true;
#endif

  return false;
}

//______________________________________________________________________________
bool
DCAnalyzer::TrackSearchSdcIn( void )
{
  static const int MinLayer = gUser.GetParameter("MinLayerSdcIn");

// #if SdcIn_Pair //Pair Plane Tracking Routine for SdcIn

// # if UseSsdCluster
//   int ntrack =
//     track::LocalTrackSearchSsdOutSdcIn( m_SsdInClCont, m_SsdOutClCont,
// 					m_SdcInHC, PPInfoSdcIn, NPPInfoSdcIn,
// 					m_SsdOutSdcInTC, MinLayer );
//   if( ntrack>=0 ) return true;
//   if( ntrack<0 ) m_much_combi[k_SdcIn]++;

// #  if SdcIn_SsdPreTrack
//   ClearTracksSsdOutSdcIn();

//   track::PreTrackSearchSsdXY( m_SsdInClCont, m_SsdOutClCont,
// 			      m_SsdXTC, m_SsdYTC );

//   ntrack =
//     track::LocalTrackSearchSsdOutSdcIn( m_SdcInHC, PPInfoSdcIn, NPPInfoSdcIn,
// 					m_SsdXTC, m_SsdYTC, m_SsdOutSdcInTC,
// 					MinLayer );
//   if( ntrack<0 ) m_much_combi[k_SdcIn]++;
//   else return true;

// #  endif
// #  if SdcIn_Deletion
//   ClearTracksSsdOutSdcIn();

//   track::LocalTrackSearchSsdOutSdcIn( m_SsdInClCont, m_SsdOutClCont,
// 				      m_SdcInHC, PPInfoSdcIn, NPPInfoSdcIn,
// 				      m_SsdOutSdcInTC, MinLayer, true );
// #  endif

// # else // UseSsdCluster
//   track::LocalTrackSearchSsdOutSdcIn( m_SsdInHC, m_SsdOutHC,
// 				      m_SdcInHC, PPInfoSdcIn, NPPInfoSdcIn,
// 				      m_SsdOutSdcInTC, MinLayer );
// # endif
//   // track::LocalTrackSearch( m_SdcInHC, PPInfoSdcIn, NPPInfoSdcIn, m_SdcInTC, MinLayer );
//   return true;
// #endif

// #if SdcIn_XUV  //XUV Tracking Routine for SdcIn
//   track::LocalTrackSearchVUX( m_SdcInHC, PPInfoSdcIn, NPPInfoSdcIn, m_SdcInTC, MinLayer );
//   return true;
// #endif
  // track::LocalTrackSearch( m_SdcInHC, PPInfoSdcIn, NPPInfoSdcIn, m_SdcInTC, MinLayer );
  track::LocalTrackSearchSdcInFiber( m_SdcInHC, PPInfoSdcIn, NPPInfoSdcIn, m_SdcInTC, MinLayer );
  return true;
  // return false;
}

//______________________________________________________________________________
bool
DCAnalyzer::TrackSearchSdcOut( void )
{
  static const int MinLayer = gUser.GetParameter("MinLayerSdcOut");

  track::LocalTrackSearchSdcOut( m_SdcOutHC, PPInfoSdcOut, NPPInfoSdcOut,
				 m_SdcOutTC, MinLayer );

  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::TrackSearchSdcOut( const Hodo2HitContainer& TOFCont )
{
  static const int MinLayer = gUser.GetParameter("MinLayerSdcOut");

  if( !DecodeTOFHits( TOFCont ) ) return false;

#if 0
  for( std::size_t i=0, n=m_TOFHC.size(); i<n; ++i ){
    m_TOFHC[i]->Print();
  }
#endif

  track::LocalTrackSearchSdcOut( m_TOFHC, m_SdcOutHC, PPInfoSdcOut, NPPInfoSdcOut+2,
				 m_SdcOutTC, MinLayer );

  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::TrackSearchSdcOut( const HodoClusterContainer& TOFCont )
{
  static const int MinLayer = gUser.GetParameter("MinLayerSdcOut");

  if( !DecodeTOFHits( TOFCont ) ) return false;

  track::LocalTrackSearchSdcOut( m_TOFHC, m_SdcOutHC, PPInfoSdcOut, NPPInfoSdcOut+2,
				 m_SdcOutTC, MinLayer );

  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::TrackSearchSsdIn( void )
{
  static const int MinLayer = gUser.GetParameter("MinLayerSsdIn");

#if UseSsdCluster
  track::LocalTrackSearchSsdIn( m_SsdInClCont, m_SsdInTC, MinLayer );
#else
  track::LocalTrackSearchSsdIn( m_SsdInHC, m_SsdInTC, MinLayer );
#endif

  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::TrackSearchSsdOut( void )
{
  static const int MinLayer = gUser.GetParameter("MinLayerSsdOut");

  track::LocalTrackSearchSsdOut( m_SsdOutHC, m_SsdOutTC, MinLayer );

  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::TrackSearchSsdInXY( void )
{
#if UseSsdCluster
  track::LocalTrackSearchSsdInXY( m_SsdInClCont, m_SsdInXTC, m_SsdInYTC );
#else
  track::LocalTrackSearchSsdInXY( m_SsdInHC, m_SsdInXTC, m_SsdInYTC );
#endif

  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::TrackSearchBcOutSdcIn( void )
{
  static const int MinLayer = gUser.GetParameter("MinLayerBcOutSdcIn");

  track::LocalTrackSearchBcOutSdcIn( m_BcOutHC, PPInfoBcOut,
				     m_SdcInHC, PPInfoSdcIn,
				     NPPInfoBcOut,NPPInfoSdcIn,
				     m_BcOutSdcInTC, MinLayer );

  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::TrackSearchBcOutSsdIn( void )
{
  static const int MinLayer = gUser.GetParameter("MinLayerBcOutSsdIn");

  track::LocalTrackSearchBcOutSsdIn( m_BcOutHC, PPInfoBcOut, NPPInfoBcOut,
				     m_SsdInHC,
				     m_BcOutSsdInTC, MinLayer );

  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::TrackSearchSsdOutSdcIn( void )
{
  static const int MinLayer = gUser.GetParameter("MinLayerSsdOutSdcIn");

  track::LocalTrackSearchSsdOutSdcIn( m_SsdInHC, m_SsdOutHC,
				      m_SdcInHC, PPInfoSdcIn, NPPInfoSdcIn,
				      m_SsdOutSdcInTC, MinLayer );

  // track::LocalTrackSearchSsdOutSdcIn( m_SsdOutHC,
  // 			       m_SdcInHC, PPInfoSdcIn, NPPInfoSdcIn,
  // 			       m_SsdOutSdcInTC, MinLayer );

  return true;
}

//______________________________________________________________________________
#if UseBcIn
bool
DCAnalyzer::TrackSearchK18U2D( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  ClearK18TracksU2D();

  int nIn  = m_BcInTC.size();
  int nOut = m_BcOutTC.size();

# if 0
  hddaq::cout << "**************************************" << std::endl;
  hddaq::cout << func_name << ": #TracksIn=" << std::setw(3) << nIn
	      << " #TracksOut=" << std::setw(3) << nOut << std::endl;
# endif

  if( nIn==0 || nOut==0 ) return true;

  for( int iIn=0; iIn<nIn; ++iIn ){
    DCLocalTrack *trIn = m_BcInTC[iIn];
# if 0
    hddaq::cout << "TrackIn  :" << std::setw(2) << iIn
		<< " X0=" << trIn->GetX0() << " Y0=" << trIn->GetY0()
		<< " U0=" << trIn->GetU0() << " V0=" << trIn->GetV0()
		<< std::endl;
# endif
    if( !trIn->GoodForTracking() ||
	trIn->GetX0()<MinK18InX || trIn->GetX0()>MaxK18InX ||
	trIn->GetY0()<MinK18InY || trIn->GetY0()>MaxK18InY ||
	trIn->GetU0()<MinK18InU || trIn->GetU0()>MaxK18InU ||
	trIn->GetV0()<MinK18InV || trIn->GetV0()>MaxK18InV ) continue;
    for( int iOut=0; iOut<nOut; ++iOut ){
      DCLocalTrack *trOut=m_BcOutTC[iOut];
# if 0
      hddaq::cout << "TrackOut :" << std::setw(2) << iOut
		  << " X0=" << trOut->GetX0() << " Y0=" << trOut->GetY0()
		  << " U0=" << trOut->GetU0() << " V0=" << trOut->GetV0()
		  << std::endl;
# endif
      if( !trOut->GoodForTracking() ||
	  trOut->GetX0()<MinK18OutX || trOut->GetX0()>MaxK18OutX ||
	  trOut->GetY0()<MinK18OutY || trOut->GetY0()>MaxK18OutY ||
	  trOut->GetU0()<MinK18OutU || trOut->GetU0()>MaxK18OutU ||
	  trOut->GetV0()<MinK18OutV || trOut->GetV0()>MaxK18OutV ) continue;

# if 0
      hddaq::cout << func_name << ": In -> " << trIn->GetChiSquare()
		  << " (" << std::setw(2) << trIn->GetNHit() << ") "
		  << "Out -> " << trOut->GetChiSquare()
		  << " (" << std::setw(2) << trOut->GetNHit() << ") "
		  << std::endl;
# endif

      K18TrackU2D *track=new K18TrackU2D( trIn, trOut, pK18 );
      if( track && track->DoFit() )
	m_K18U2DTC.push_back(track);
      else
	delete track;
    }
  }

# if 0
  hddaq::cout << "********************" << std::endl;
  {
    int nn = m_K18U2DTC.size();
    hddaq::cout << func_name << ": Before sorting. #Track="
		<< nn << std::endl;
    for( int i=0; i<nn; ++i ){
      K18TrackU2D *tp = m_K18U2DTC[i];
      hddaq::cout << std::setw(3) << i
		  << " ChiSqr=" << tp->chisquare()
		  << " Delta=" << tp->Delta()
		  << " P=" << tp->P() << "\n";
      //      hddaq::cout<<"********************"<<std::endl;
      //      hddaq::cout << "In :"
      // 	       << " X " << tp->Xin() << "(" << tp->TrackIn()->GetX0() << ")"
      // 	       << " Y " << tp->Yin() << "(" << tp->TrackIn()->GetY0() << ")"
      // 	       << " U " << tp->Uin() << "(" << tp->TrackIn()->GetU0() << ")"
      // 	       << " V " << tp->Vin() << "(" << tp->TrackIn()->GetV0() << ")"
      // 	       << "\n";
      //      hddaq::cout << "Out:"
      // 	       << " X " << tp->Xout() << "(" << tp->TrackOut()->GetX0() << ")"
      // 	       << " Y " << tp->Yout() << "(" << tp->TrackOut()->GetY0() << ")"
      // 	       << " U " << tp->Uout() << "(" << tp->TrackOut()->GetU0() << ")"
      // 	       << " V " << tp->Vout() << "(" << tp->TrackOut()->GetV0() << ")"
      // 	       << std::endl;
    }
  }
# endif

  std::sort( m_K18U2DTC.begin(), m_K18U2DTC.end(), K18TrackU2DComp() );

  return true;
}
#endif // UseBcIn

//______________________________________________________________________________
bool
DCAnalyzer::TrackSearchK18D2U( const std::vector<double>& XinCont )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  ClearK18TracksD2U();

  std::size_t nIn  = XinCont.size();
  std::size_t nOut = m_BcOutTC.size();

  if( nIn==0 || nOut==0 )
    return true;

  for( std::size_t iIn=0; iIn<nIn; ++iIn ){
    double LocalX = XinCont.at( iIn );
    for( std::size_t iOut=0; iOut<nOut; ++iOut ){
      DCLocalTrack *trOut = m_BcOutTC[iOut];
#if 0
      hddaq::cout << "TrackOut :" << std::setw(2) << iOut
		  << " X0=" << trOut->GetX0() << " Y0=" << trOut->GetY0()
		  << " U0=" << trOut->GetU0() << " V0=" << trOut->GetV0()
		  << std::endl;
#endif
      if( !trOut->GoodForTracking() ||
	  trOut->GetX0()<MinK18OutX || trOut->GetX0()>MaxK18OutX ||
	  trOut->GetY0()<MinK18OutY || trOut->GetY0()>MaxK18OutY ||
	  trOut->GetU0()<MinK18OutU || trOut->GetU0()>MaxK18OutU ||
	  trOut->GetV0()<MinK18OutV || trOut->GetV0()>MaxK18OutV ) continue;

#if 0
      hddaq::cout << func_name
		  << "Out -> " << trOut->GetChiSquare()
		  << " (" << std::setw(2) << trOut->GetNHit() << ") "
		  << std::endl;
#endif

      K18TrackD2U *track = new K18TrackD2U( LocalX, trOut, pK18 );
      if( !track ) continue;
      track->CalcMomentumD2U();
      m_K18D2UTC.push_back(track);
    }
  }

#if 0
  hddaq::cout<<"********************"<<std::endl;
  {
    int nn = m_K18D2UTC.size();
    hddaq::cout << func_name << ": Before sorting. #Track="
		<< nn << std::endl;
    for( int i=0; i<nn; ++i ){
      K18TrackD2U *tp = m_K18D2UTC[i];
      hddaq::cout << std::setw(3) << i
		  << " ChiSqr=" << tp->chisquare()
		  << " Delta=" << tp->Delta()
		  << " P=" << tp->P() << "\n";
//      hddaq::cout<<"********************"<<std::endl;
//      hddaq::cout << "In :"
// 	       << " X " << tp->Xin() << "(" << tp->TrackIn()->GetX0() << ")"
// 	       << " Y " << tp->Yin() << "(" << tp->TrackIn()->GetY0() << ")"
// 	       << " U " << tp->Uin() << "(" << tp->TrackIn()->GetU0() << ")"
// 	       << " V " << tp->Vin() << "(" << tp->TrackIn()->GetV0() << ")"
// 	       << "\n";
//      hddaq::cout << "Out:"
// 	       << " X " << tp->Xout() << "(" << tp->TrackOut()->GetX0() << ")"
// 	       << " Y " << tp->Yout() << "(" << tp->TrackOut()->GetY0() << ")"
// 	       << " U " << tp->Uout() << "(" << tp->TrackOut()->GetU0() << ")"
// 	       << " V " << tp->Vout() << "(" << tp->TrackOut()->GetV0() << ")"
// 	       << std::endl;
   }
 }
#endif

  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::TrackSearchKurama( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  ClearKuramaTracks();

  std::size_t nIn  = GetNtracksSdcIn();
  std::size_t nOut = GetNtracksSdcOut();

  if( nIn==0 || nOut==0 ) {
    return true;
  }
  for( std::size_t iIn=0; iIn<nIn; ++iIn ){
    DCLocalTrack *trIn = GetTrackSdcIn( iIn );
    if( !trIn->GoodForTracking() ) continue;
    for( std::size_t iOut=0; iOut<nOut; ++iOut ){
      DCLocalTrack * trOut = GetTrackSdcOut( iOut );
      if( !trOut->GoodForTracking() ) continue;

      KuramaTrack *trKurama = new KuramaTrack( trIn, trOut );
      if( !trKurama ) continue;
      double u0In    = trIn->GetU0();
      double u0Out   = trOut->GetU0();
      double bending = u0Out - u0In;
      double p[3] = { 0.08493, 0.2227, 0.01572 };
      double initial_momentum = p[0] + p[1]/( bending-p[2] );
      if( bending>0. && initial_momentum>0. ) {
	trKurama->SetInitialMomentum( initial_momentum );
      } else {
	// trKurama->SetInitialMomentum( 1.8 );
	trKurama->SetInitialMomentum( 1. );
      }
      if( trKurama->DoFit() && trKurama->chisqr()<MaxChiSqrKuramaTrack ){
	m_KuramaTC.push_back( trKurama );
      }
      else{
	// trKurama->Print( "in "+func_name );
	delete trKurama;
      }
    }// for( iOut )
  }// for( iIn )

  std::sort( m_KuramaTC.begin(), m_KuramaTC.end(), KuramaTrackComp() );

#if 0
  PrintKurama( "Before Deleting" );
#endif

  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::TrackSearchKurama( double initial_momentum )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  ClearKuramaTracks();

  int nIn  = GetNtracksSdcIn();
  int nOut = GetNtracksSdcOut();

  if( nIn==0 || nOut==0 ) return true;

  for( int iIn=0; iIn<nIn; ++iIn ){
    DCLocalTrack *trIn = GetTrackSdcIn( iIn );
    if( !trIn->GoodForTracking() ) continue;
    for( int iOut=0; iOut<nOut; ++iOut ){
      DCLocalTrack * trOut = GetTrackSdcOut( iOut );
      if( !trOut->GoodForTracking() ) continue;
      KuramaTrack *trKurama = new KuramaTrack( trIn, trOut );
      if( !trKurama ) continue;
      trKurama->SetInitialMomentum( initial_momentum );
      if( trKurama->DoFit() && trKurama->chisqr()<MaxChiSqrKuramaTrack ){
	m_KuramaTC.push_back( trKurama );
      }
      else{
	trKurama->Print( " in "+func_name );
	delete trKurama;
      }
    }// for( iOut )
  }// for( iIn )

  std::sort( m_KuramaTC.begin(), m_KuramaTC.end(), KuramaTrackComp() );

#if 0
  PrintKurama( "Before Deleting" );
#endif

  return true;
}

//______________________________________________________________________________
void
DCAnalyzer::ClearDCHits( void )
{
#if UseBcIn
  ClearBcInHits();
#endif
  ClearBcOutHits();
  ClearSdcInHits();
  ClearSdcOutHits();
  ClearSsdInHits();
  ClearSsdOutHits();
  ClearTOFHits();
}

//______________________________________________________________________________
#if UseBcIn
void
DCAnalyzer::ClearBcInHits( void )
{
  del::ClearContainerAll( m_TempBcInHC );
  del::ClearContainerAll( m_BcInHC );
  del::ClearContainerAll( m_MWPCClCont );
}
#endif

//______________________________________________________________________________
void
DCAnalyzer::ClearBcOutHits( void )
{
  del::ClearContainerAll( m_BcOutHC );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearSdcInHits( void )
{
  del::ClearContainerAll( m_SdcInHC );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearSdcOutHits( void )
{
  del::ClearContainerAll( m_SdcOutHC );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearSsdInHits( void )
{
  del::ClearContainerAll( m_SsdInHC );
  del::ClearContainerAll( m_SsdInClCont );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearSsdOutHits( void )
{
  del::ClearContainerAll( m_SsdOutHC );
  del::ClearContainerAll( m_SsdOutClCont );
}

// //______________________________________________________________________________
// void
// DCAnalyzer::ClearSftHits( void )
// {
//   del::ClearContainerAll( m_SftHC );
//   del::ClearContainerAll( m_SftClCont );
// }

//______________________________________________________________________________
void
DCAnalyzer::ClearVtxHits( void )
{
  del::ClearContainer( m_VtxPoint );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearTOFHits( void )
{
  del::ClearContainer( m_TOFHC );
}

//______________________________________________________________________________
#if UseBcIn
void
DCAnalyzer::ClearTracksBcIn( void )
{
  del::ClearContainer( m_BcInTC );
}
#endif

//______________________________________________________________________________
void
DCAnalyzer::ClearTracksBcOut( void )
{
  del::ClearContainer( m_BcOutTC );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearTracksSdcIn( void )
{
  del::ClearContainer( m_SdcInTC );
  del::ClearContainerAll( m_SdcInExTC );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearTracksSdcOut( void )
{
  del::ClearContainer( m_SdcOutTC );
  del::ClearContainerAll( m_SdcOutExTC );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearTracksSsdIn( void )
{
  del::ClearContainer( m_SsdInTC );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearTracksSsdOut( void )
{
  del::ClearContainer( m_SsdOutTC );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearTracksSsdXY( void )
{
  del::ClearContainer( m_SsdXTC );
  del::ClearContainer( m_SsdYTC );
  del::ClearContainer( m_SsdInXTC );
  del::ClearContainer( m_SsdInYTC );
}

//______________________________________________________________________________
#if UseBcIn
void
DCAnalyzer::ClearK18TracksU2D( void )
{
  del::ClearContainer( m_K18U2DTC );
}
#endif

//______________________________________________________________________________
void
DCAnalyzer::ClearK18TracksD2U( void )
{
  del::ClearContainer( m_K18D2UTC );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearKuramaTracks( void )
{
  del::ClearContainer( m_KuramaTC );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearTracksBcOutSdcIn( void )
{
  del::ClearContainer( m_BcOutSdcInTC );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearTracksBcOutSsdIn( void )
{
  del::ClearContainer( m_BcOutSsdInTC );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearTracksSsdOutSdcIn( void )
{
  del::ClearContainer( m_SsdOutSdcInTC );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearTracksSdcInSdcOut( void )
{
  del::ClearContainer( m_SdcInSdcOutTC );
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcMWPCHits( std::vector<DCHitContainer>& cont,
			    bool applyRecursively )
{
  const std::size_t n = cont.size();
  for( std::size_t l=0; l<n; ++l ){
    const std::size_t m = cont[l].size();
    for( std::size_t i=0; i<m; ++i ){
      DCHit *hit = (cont[l])[i];
      if( !hit ) continue;
      hit->ReCalcMWPC(applyRecursively);
    }
  }
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcDCHits( std::vector<DCHitContainer>& cont,
			  bool applyRecursively )
{
  const std::size_t n = cont.size();
  for( std::size_t l=0; l<n; ++l ){
    const std::size_t m = cont[l].size();
    for( std::size_t i=0; i<m; ++i ){
      DCHit *hit = (cont[l])[i];
      if( !hit ) continue;
      hit->ReCalcDC(applyRecursively);
    }
  }
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcDCHits( bool applyRecursively )
{
#if UseBcIn
  ReCalcMWPCHits( m_TempBcInHC, applyRecursively );
  ReCalcMWPCHits( m_BcInHC, applyRecursively );
#endif

  ReCalcDCHits( m_BcOutHC, applyRecursively );
  ReCalcDCHits( m_SdcInHC, applyRecursively );
  ReCalcDCHits( m_SdcOutHC, applyRecursively );
  ReCalcDCHits( m_SsdInHC, applyRecursively );
  ReCalcDCHits( m_SsdOutHC, applyRecursively );

  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcTrack( DCLocalTrackContainer& cont,
			 bool applyRecursively )
{
  const std::size_t n = cont.size();
  for( std::size_t i=0; i<n; ++i ){
    DCLocalTrack *track = cont[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcTrack( K18TrackD2UContainer& cont,
			 bool applyRecursively )
{
  const std::size_t n = cont.size();
  for( std::size_t i=0; i<n; ++i ){
    K18TrackD2U *track = cont[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcTrack( KuramaTrackContainer& cont,
			 bool applyRecursively )
{
  const std::size_t n = cont.size();
  for( std::size_t i=0; i<n; ++i ){
    KuramaTrack *track = cont[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}

//______________________________________________________________________________
#if UseBcIn
bool
DCAnalyzer::ReCalcTrackBcIn( bool applyRecursively )
{
  return ReCalcTrack( m_BcInTC );
}
#endif

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcTrackBcOut( bool applyRecursively )
{
  return ReCalcTrack( m_BcOutTC );
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcTrackSdcIn( bool applyRecursively )
{
  return ReCalcTrack( m_SdcInTC );
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcTrackSdcOut( bool applyRecursively )
{
  return ReCalcTrack( m_SdcOutTC );
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcTrackSsdIn( bool applyRecursively )
{
  return ReCalcTrack( m_SsdInTC );
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcTrackSsdOut( bool applyRecursively )
{
  return ReCalcTrack( m_SsdOutTC );
}

//______________________________________________________________________________
#if UseBcIn
bool
DCAnalyzer::ReCalcK18TrackU2D( bool applyRecursively )
{
  int n = m_K18U2DTC.size();
  for( int i=0; i<n; ++i ){
    K18TrackU2D *track = m_K18U2DTC[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}
#endif

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcK18TrackD2U( bool applyRecursively )
{
  return ReCalcTrack( m_K18D2UTC, applyRecursively );
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcKuramaTrack( bool applyRecursively )
{
  return ReCalcTrack( m_KuramaTC, applyRecursively );
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcAll( void )
{
  ReCalcDCHits();
#if UseBcIn
  ReCalcTrackBcIn();
  ReCalcK18TrackU2D();
#endif
  ReCalcTrackBcOut();
  ReCalcTrackSdcIn();
  ReCalcTrackSdcOut();
  ReCalcTrackSsdIn();
  ReCalcTrackSsdOut();

  //ReCalcK18TrackD2U();
  ReCalcKuramaTrack();

  return true;
}

//______________________________________________________________________________
// int
// clusterizeMWPCHit(const DCHitContainer& hits,
// 		  MWPCClusterContainer& clusters)
// {
//   if (!clusters.empty()){
//       std::for_each(clusters.begin(), clusters.end(), DeleteObject());
//       clusters.clear();
//   }

//   const int nhits = hits.size();
//   //   hddaq::cout << "#D " << __func__ << " " << nhits << std::endl;
//   if (nhits==0)
//     return 0;

//   int n = 0;
//   for (int i=0; i<nhits; ++i){
//     const DCHit* h = hits[i];
//     if (!h)
//       continue;
//     n += h->GetTdcSize();
//   }

//   DCHitContainer singleHits;
//   singleHits.reserve(n);
//   for (int i=0; i<nhits; ++i){
//     const DCHit* h = hits[i];
//     if (!h)
//       continue;
//     int nn = h->GetTdcSize();
//     for (int ii=0; ii<nn; ++ii){
//       DCHit* htmp = new DCHit(h->GetLayer(), h->GetWire());
//       htmp->SetTdcVal(h->GetTdcVal());
//       htmp->SetTdcTrailing(h->GetTdcTrailing());
//       htmp->SetTrailingTime(h->GetTrailingTime());
//       htmp->SetDriftTime(h->GetDriftTime());
//       htmp->SetDriftLength(h->GetDriftLength());
//       htmp->SetTiltAngle(h->GetTiltAngle());
//       htmp->SetWirePosition(h->GetWirePosition());
//       htmp->setRangeCheckStatus(h->rangecheck(), 0);
//       singleHits.push_back(htmp);
//     }
//   }

//   std::vector<std::deque<bool> > flag(n, std::deque<bool>(n, false));
//   n = singleHits.size();
//   for (int i=0;  i<n; ++i){
//     flag[i][i] = true;
//     const DCHit* h1 = singleHits[i];
//     //       h1->print("h1");
//     for (int j=i+1; j<n; ++j){
//       const DCHit* h2 = singleHits[j];
//       // 	  h2->print("h2");
//       // 	  hddaq::cout << " (i,j) = (" << i << ", " << j << ")" << std::endl;
//       bool val
// 	= isConnectable(h1->GetWirePosition(),
// 			h1->GetDriftTime(),
// 			h1->GetTrailingTime(),
// 			h2->GetWirePosition(),
// 			h2->GetDriftTime(),
// 			h2->GetTrailingTime(),
// 			kMWPCClusteringWireExtension,
// 			kMWPCClusteringTimeExtension);
//       // 	  hddaq::cout << "#D val = " << val << std::endl;
//       flag[i][j] = val;
//       flag[j][i] = val;
//     }
//   }

//   //   hddaq::cout << "#D " << __func__ << "  before " << std::endl;
//   //   printConnectionFlag(flag);

//   const int maxLoop = static_cast<int>(std::log(x)/std::log(2.))+1;
//   for (int loop=0; loop<maxLoop; ++loop){
//     std::vector<std::deque<bool> > tmp(n, std::deque<bool>(n, false));
//     for (int i=0; i<n; ++i){
//       for (int j=i; j<n; ++j){
// 	for (int k=0; k<n; ++k){
// 	  tmp[i][j] |= (flag[i][k] && flag[k][j]);
// 	  tmp[j][i] = tmp[i][j];
// 	}
//       }
//     }
//     flag = tmp;
//     //       hddaq::cout << " n iteration = " << loop << std::endl;
//     //       printConnectionFlag(flag);
//   }

//   //   hddaq::cout << "#D " << __func__ << "  after " << std::endl;
//   //   printConnectionFlag(flag);

//   std::set<int> checked;
//   for (int i=0; i<n; ++i){
//     if (checked.find(i)!=checked.end())
//       continue;
//     MWPCCluster* c = 0;
//     for (int j=i; j<n; ++j){
//       if (flag[i][j]){
// 	checked.insert(j);
// 	if (!c) {
// 	  c = new MWPCCluster;
// 	  // 		  hddaq::cout << " new cluster " << std::endl;
// 	}
// 	// 	      hddaq::cout << " " << i << "---" << j << std::endl;
// 	c->Add(singleHits[j]);
//       }
//     }

//     if (c){
//       c->Calculate();
//       clusters.push_back(c);
//     }
//   }

//   //   hddaq::cout << " end of " << __func__
//   // 	    << " : n = " << n << ", " << checked.size()
//   // 	    << std::endl;

//   //   hddaq::cout << __func__ << " n clusters = " << clusters.size() << std::endl;

//   return clusters.size();
// }

//______________________________________________________________________________
bool
DCAnalyzer::ClusterizeSsd( void )
{
  for( int l=1; l<NumOfLayersSsdIn+1; ++l )
    ClusterizeSsd( m_SsdInHC[l], m_SsdInClCont[l] );
  for( int l=1; l<NumOfLayersSsdOut+1; ++l )
    ClusterizeSsd( m_SsdOutHC[l], m_SsdOutClCont[l] );
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::ClusterizeSsd( const DCHitContainer& HitCont,
			   SsdClusterContainer& ClCont,
			   double MaxTimeDiff )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  del::ClearContainer( ClCont );

  const std::size_t nh = HitCont.size();
  if( nh==0 )
    return false;

  std::vector<int> flag( nh, 0 );

  for( std::size_t i=0; i<nh; ++i ){
    if( flag[i]>0 )
      continue;
    DCHitContainer CandCont;
    DCHit* hitA = HitCont[i];
    if( !hitA || !hitA->IsGoodWaveForm() )
      continue;
    CandCont.push_back( hitA );
    ++flag[i];

    for( std::size_t j=0; j<nh; ++j ){
      if( CandCont.size()==SsdCluster::MaxClusterSize() )
	break;
      if( i==j || flag[j]>0 )
	continue;
      DCHit* hitB = HitCont[j];
      if( !hitB || !hitB->IsGoodWaveForm() )
	continue;
      if( IsClusterable( CandCont, hitB ) ){
	CandCont.push_back( hitB );
	++flag[j];
      }
    }
    SsdCluster *cluster = new SsdCluster( CandCont );
    if( cluster ) ClCont.push_back( cluster );
  }

  return true;
}

//______________________________________________________________________________
void
DCAnalyzer::DoTimeCorrectionSsd( const std::vector<double>& t0 )
{
  DoTimeCorrectionSsdIn( t0 );
  DoTimeCorrectionSsdOut( t0 );
}

//______________________________________________________________________________
void
DCAnalyzer::DoTimeCorrectionSsdIn( const std::vector<double>& t0 )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdIn+1; ++layer ){
    const DCHitContainer& HitCont = GetSsdInHC(layer);
    const std::size_t nh = HitCont.size();
    for( std::size_t i=0; i<nh; ++i ){
      DCHit *hit = HitCont[i];
      if( !hit ) continue;
      double offset = 0.;
      switch( layer ){
      case 1: case 2: offset = t0[1]; break;
      case 3: case 4: offset = t0[1]; break;
      }
      hit->DoTimeCorrection( offset );
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::DoTimeCorrectionSsdOut( const std::vector<double>& t0 )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdOut+1; ++layer ){
    const DCHitContainer &HitCont = GetSsdOutHC(layer);
    const std::size_t nh = HitCont.size();
    for( std::size_t i=0; i<nh; ++i ){
      DCHit *hit = HitCont[i];
      if( !hit ) continue;
      double offset = 0.;
      switch(layer){
      case 1: case 2: offset = t0[1]; break;
      case 3: case 4: offset = t0[1]; break;
      }
      hit->DoTimeCorrection( offset );
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::ResetStatusSsd( void )
{
  ResetStatusSsdIn();
  ResetStatusSsdOut();
}

//______________________________________________________________________________
void
DCAnalyzer::ResetStatusSsdIn( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdIn+1;++layer){
    const DCHitContainer &HitCont = GetSsdInHC(layer);
    const std::size_t nh = HitCont.size();
    for( std::size_t i=0; i<nh; ++i ){
      DCHit* hit = HitCont[i];
      if( !hit ) continue;
      hit->SetGoodWaveForm( true );
    }
    const SsdClusterContainer &ClCont = GetClusterSsdIn(layer);
    const std::size_t ncl = ClCont.size();
    for( std::size_t i=0; i<ncl; ++i ){
      SsdCluster *cluster = ClCont[i];
      if( !cluster ) continue;
      cluster->GoodForAnalysis( true );
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::ResetStatusSsdOut( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdOut+1;++layer){
    const DCHitContainer &HitCont = GetSsdOutHC(layer);
    const std::size_t nh = HitCont.size();
    for( std::size_t i=0; i<nh; ++i ){
      DCHit *hit = HitCont[i];
      if( !hit ) continue;
      hit->SetGoodWaveForm( true );
    }
    const SsdClusterContainer &ClCont = GetClusterSsdOut(layer);
    const std::size_t ncl = ClCont.size();
    for( std::size_t i=0; i<ncl; ++i ){
      SsdCluster *cluster = ClCont[i];
      if( !cluster ) continue;
      cluster->GoodForAnalysis( true );
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::SlopeFilterSsd( void )
{
  SlopeFilterSsdIn();
  SlopeFilterSsdOut();
}

//______________________________________________________________________________
void
DCAnalyzer::SlopeFilterSsdIn( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdIn+1;++layer){
    const DCHitContainer &HitCont = GetSsdInHC(layer);
    const std::size_t nh = HitCont.size();
    for( std::size_t i=0; i<nh; ++i ){
      DCHit *hit = HitCont[i];
      if( !hit ) continue;
      if( !hit->IsGoodWaveForm() ) continue;

      std::vector<double> waveform = hit->GetWaveform();
      if( waveform.size()!=NumOfSampleSSD ) continue;

      bool slope[NumOfSampleSSD];
      for( int s=0; s<NumOfSampleSSD; ++s ){
	if( s>0 ) slope[s-1] = ( waveform[s]>waveform[s-1] );
      }

      bool status = true;
#if SlopeFilter_Tight
      status = slope[0] && slope[1] && slope[2]
      	&& !slope[4] && !slope[5] && !slope[6];
#elif SlopeFilter_Wide
      status = slope[0] && slope[1] && !slope[5] && !slope[6];
#endif
      hit->SetGoodWaveForm(status);
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::SlopeFilterSsdOut( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdOut+1;++layer){
    const DCHitContainer &HitCont = GetSsdOutHC(layer);
    const std::size_t nh = HitCont.size();
    for( std::size_t i=0; i<nh; ++i ){
      DCHit *hit = HitCont[i];
      if( !hit ) continue;
      if( !hit->IsGoodWaveForm() ) continue;

      std::vector<double> waveform = hit->GetWaveform();
      if( waveform.size()!=NumOfSampleSSD ) continue;

      bool slope[NumOfSampleSSD];
      for( int s=0; s<NumOfSampleSSD; ++s ){
	if(s>0) slope[s-1] = waveform[s] > waveform[s-1];
      }

      bool status = true;
#if SlopeFilter_Tight
      status = slope[0] && slope[1] && slope[2]
      	&& !slope[4] && !slope[5] && !slope[6];
#elif SlopeFilter_Wide
      status = slope[0] && slope[1] && !slope[5] && !slope[6];
#endif
      hit->SetGoodWaveForm(status);
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::AdcPeakHeightFilterSsd( double min, double max )
{
  AdcPeakHeightFilterSsdIn( min, max );
  AdcPeakHeightFilterSsdOut( min, max );
}

//______________________________________________________________________________
void
DCAnalyzer::AdcPeakHeightFilterSsdIn( double min, double max )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdIn+1;++layer ){
    const DCHitContainer &HitCont = GetSsdInHC(layer);
    const std::size_t nh = HitCont.size();
    for( std::size_t i=0; i<nh; ++i ){
      DCHit *hit = HitCont[i];
      if( !hit ) continue;
      if( !hit->IsGoodWaveForm() ) continue;
      double peak_height = hit->GetAdcPeakHeight();
      double pedestal    = hit->GetPedestal();
      peak_height -= pedestal;
      if( peak_height<min || max<peak_height ){
	hit->SetGoodWaveForm( false );
      }
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::AdcPeakHeightFilterSsdOut( double min, double max )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdOut+1;++layer ){
    const DCHitContainer &HitCont = GetSsdOutHC(layer);
    const std::size_t nh = HitCont.size();
    for( std::size_t i=0; i<nh; ++i ){
      DCHit *hit = HitCont[i];
      if( !hit ) continue;
      if( !hit->IsGoodWaveForm() ) continue;
      double peak_height = hit->GetAdcPeakHeight();
      double pedestal    = hit->GetPedestal();
      peak_height -= pedestal;
      if( peak_height<min || max<peak_height ){
	hit->SetGoodWaveForm( false );
      }
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::AmplitudeFilterSsd( double min, double max, bool cluster_flag /*=false*/ )
{
  AmplitudeFilterSsdIn( min, max, cluster_flag );
  AmplitudeFilterSsdOut( min, max, cluster_flag );
}

//______________________________________________________________________________
void
DCAnalyzer::AmplitudeFilterSsdIn( double min, double max, bool cluster_flag /*=false*/ )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdIn+1;++layer ){
    if( !cluster_flag ){
      const DCHitContainer &HitCont = GetSsdInHC(layer);
      const std::size_t nh = HitCont.size();
      for( std::size_t i=0; i<nh; ++i ){
	DCHit *hit = HitCont[i];
	if( !hit ) continue;
	if( !hit->IsGoodWaveForm() ) continue;
	double amplitude = hit->GetAmplitude();
	if( amplitude<min || max<amplitude ){
	  hit->SetGoodWaveForm(false);
	}
      }
    }
    if( cluster_flag ){
      const SsdClusterContainer &ClCont = GetClusterSsdIn(layer);
      const std::size_t ncl = ClCont.size();
      for( std::size_t i=0; i<ncl; ++i ){
	SsdCluster *cluster = ClCont[i];
	if( !cluster ) continue;
	if( !cluster->GoodForAnalysis() ) continue;
	double amplitude = cluster->Amplitude();
	if( amplitude<min || max<amplitude ){
	  cluster->GoodForAnalysis(false);
	}
      }
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::AmplitudeFilterSsdOut( double min, double max, bool cluster_flag /*=false*/ )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdOut+1;++layer){
    if( !cluster_flag ){
      const DCHitContainer &HitCont = GetSsdOutHC(layer);
      const std::size_t nh = HitCont.size();
      for( std::size_t i=0; i<nh; ++i ){
	DCHit *hit = HitCont[i];
	if( !hit ) continue;
	if( !hit->IsGoodWaveForm() ) continue;
	double amplitude = hit->GetAmplitude();
	if( amplitude<min || max<amplitude ){
	  hit->SetGoodWaveForm(false);
	}
      }
    }
    if( cluster_flag ){
      const SsdClusterContainer &ClCont = GetClusterSsdOut(layer);

      const std::size_t ncl = ClCont.size();
      for( std::size_t i=0; i<ncl; ++i ){
	SsdCluster *cluster = ClCont[i];
	if( !cluster ) continue;
	if( !cluster->GoodForAnalysis() ) continue;
	double amplitude = cluster->Amplitude();
	if( amplitude<min || max<amplitude ){
	  cluster->GoodForAnalysis(false);
	}
      }
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::DeltaEFilterSsd( double min, double max, bool cluster_flag /*=false*/ )
{
  DeltaEFilterSsdIn( min, max, cluster_flag );
  DeltaEFilterSsdOut( min, max, cluster_flag );
}

//______________________________________________________________________________
void
DCAnalyzer::DeltaEFilterSsdIn( double min, double max, bool cluster_flag /*=false*/ )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdIn+1;++layer){
    if( !cluster_flag ){
      const DCHitContainer &HitCont = GetSsdInHC(layer);
      const std::size_t nh = HitCont.size();
      for( std::size_t i=0; i<nh; ++i ){
	DCHit *hit = HitCont[i];
	if( !hit ) continue;
	if( !hit->IsGoodWaveForm() ) continue;
	double de = hit->GetDe();
	if( de<min || max<de ){
	  hit->SetGoodWaveForm(false);
	}
      }
    }
    if( cluster_flag ){
      const SsdClusterContainer &ClCont = GetClusterSsdIn(layer);
      const std::size_t ncl = ClCont.size();
      for( std::size_t i=0; i<ncl; ++i ){
	SsdCluster *cluster = ClCont[i];
	if( !cluster ) continue;
	if( !cluster->GoodForAnalysis() ) continue;
	double de = cluster->DeltaE();
	if( de<min || max<de ){
	  cluster->GoodForAnalysis(false);
	}
      }
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::DeltaEFilterSsdOut( double min, double max, bool cluster_flag /*=false*/ )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdOut+1;++layer){
    if( !cluster_flag ){
      const DCHitContainer &HitCont = GetSsdOutHC(layer);
      const std::size_t nh = HitCont.size();
      for( std::size_t i=0; i<nh; ++i ){
	DCHit *hit = HitCont[i];
	if( !hit ) continue;
	if( !hit->IsGoodWaveForm() ) continue;
	double de = hit->GetDe();
	if( de<min || max<de ){
	  hit->SetGoodWaveForm(false);
	}
      }
    }
    if( cluster_flag ){
      const SsdClusterContainer &ClCont = GetClusterSsdOut(layer);
      const std::size_t ncl = ClCont.size();
      for( std::size_t i=0; i<ncl; ++i ){
	SsdCluster *cluster = ClCont[i];
	if( !cluster ) continue;
	if( !cluster->GoodForAnalysis() ) continue;
	double de = cluster->DeltaE();
	if( de<min || max<de ){
	  cluster->GoodForAnalysis(false);
	}
      }
    }
  }
}
//______________________________________________________________________________
void
DCAnalyzer::RmsFilterSsd( double min, double max )
{
  RmsFilterSsdIn(min, max);
  RmsFilterSsdOut(min, max);
}

//______________________________________________________________________________
void
DCAnalyzer::RmsFilterSsdIn( double min, double max )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdIn+1;++layer){
    const DCHitContainer &HitCont = GetSsdInHC(layer);
    const std::size_t nh = HitCont.size();
    for( std::size_t i=0; i<nh; ++i ){
      DCHit *hit = HitCont[i];
      if(!hit) continue;
      if(!hit->IsGoodWaveForm()) continue;
      double de  = hit->GetDe();
      double rms = hit->GetRms();
      if( de<min*rms || max*rms<de )
	hit->SetGoodWaveForm(false);
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::RmsFilterSsdOut( double min, double max )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdOut+1;++layer){
    const DCHitContainer &HitCont = GetSsdOutHC(layer);
    const std::size_t nh = HitCont.size();
    for( std::size_t i=0; i<nh; ++i ){
      DCHit *hit = HitCont[i];
      if( !hit ) continue;
      if( !hit->IsGoodWaveForm() ) continue;
      double de  = hit->GetDe();
      double rms = hit->GetRms();
      if( de<min*rms||max*rms<de ){
	hit->SetGoodWaveForm(false);
      }
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::DeviationFilterSsd( double min, double max )
{
  DeviationFilterSsdIn( min, max );
  DeviationFilterSsdOut( min, max );
}

//______________________________________________________________________________
void
DCAnalyzer::DeviationFilterSsdIn( double min, double max )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdIn+1;++layer){
    const DCHitContainer &HitCont = GetSsdInHC(layer);
    const std::size_t nh = HitCont.size();
    for( std::size_t i=0; i<nh; ++i ){
      DCHit *hit = HitCont[i];
      if( !hit ) continue;
      if( !hit->IsGoodWaveForm() ) continue;
      double deviation = hit->GetDeviation();
      if( deviation<min || max<deviation ){
	hit->SetGoodWaveForm(false);
      }
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::DeviationFilterSsdOut( double min, double max )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdOut+1;++layer){
    const DCHitContainer &HitCont = GetSsdOutHC(layer);
    const std::size_t nh = HitCont.size();
    for( std::size_t i=0; i<nh; ++i ){
      DCHit *hit = HitCont[i];
      if( !hit ) continue;
      if( !hit->IsGoodWaveForm() ) continue;
      double deviation = hit->GetDeviation();
      if( deviation<min || max<deviation ){
	hit->SetGoodWaveForm(false);
      }
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::TimeFilterSsd( double min, double max, bool cluster_flag /*=false*/ )
{
  TimeFilterSsdIn( min, max, cluster_flag );
  TimeFilterSsdOut( min, max, cluster_flag );
}

//______________________________________________________________________________
void
DCAnalyzer::TimeFilterSsdIn( double min, double max, bool cluster_flag /*=false*/ )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdIn+1;++layer){
    if( !cluster_flag ){
      const DCHitContainer &HitCont = GetSsdInHC(layer);
      const std::size_t nh = HitCont.size();
      for( std::size_t i=0; i<nh; ++i ){
	DCHit *hit = HitCont[i];
	if( !hit ) continue;
	if( !hit->IsGoodWaveForm() ) continue;
	double peaktime = hit->GetPeakTime();
	if( peaktime<min || max<peaktime ){
	  hit->SetGoodWaveForm(false);
	}
      }
    }
    if( cluster_flag ){
      const SsdClusterContainer &ClCont = GetClusterSsdIn(layer);
      const std::size_t ncl = ClCont.size();
      for( std::size_t i=0; i<ncl; ++i ){
	SsdCluster *cluster = ClCont[i];
	if( !cluster ) continue;
	if( !cluster->GoodForAnalysis() ) continue;
	double time = cluster->Time();
	if( time<min || max<time ){
	  cluster->GoodForAnalysis(false);
	}
      }
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::TimeFilterSsdOut( double min, double max, bool cluster_flag /*=false*/ )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdOut+1;++layer){
    if( !cluster_flag ){
      const DCHitContainer &HitCont = GetSsdOutHC(layer);
      const std::size_t nh = HitCont.size();
      for( std::size_t i=0; i<nh; ++i ){
	DCHit *hit = HitCont[i];
	if( !hit ) continue;
	if( !hit->IsGoodWaveForm() ) continue;
	double peaktime = hit->GetPeakTime();
	if( peaktime<min || max<peaktime ){
	  hit->SetGoodWaveForm(false);
	}
      }
    }
    if( cluster_flag ){
      const SsdClusterContainer &ClCont = GetClusterSsdOut(layer);
      const std::size_t ncl = ClCont.size();
      for( std::size_t i=0; i<ncl; ++i ){
	SsdCluster *cluster = ClCont[i];
	if( !cluster ) continue;
	if( !cluster->GoodForAnalysis() ) continue;
	double time = cluster->Time();
	if( time<min || max<time ){
	  cluster->GoodForAnalysis(false);
	}
      }
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::ChisqrFilterSsd( double maxchi2 )
{
  ChisqrFilterSsdIn( maxchi2 );
  ChisqrFilterSsdOut( maxchi2 );
}

//______________________________________________________________________________
void
DCAnalyzer::ChisqrFilterSsdIn( double maxchi2 )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdIn+1;++layer){
    const DCHitContainer &HitCont = GetSsdInHC(layer);
    int nh = HitCont.size();
    for(int i=0;i<nh;++i){
      DCHit *hit = HitCont[i];
      if(!hit) continue;
      if(!hit->IsGoodWaveForm()) continue;
      double chisqr = hit->GetChisquare();
      if( chisqr>maxchi2 ){
	hit->SetGoodWaveForm(false);
      }
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::ChisqrFilterSsdOut( double maxchi2 )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( int layer=1; layer<NumOfLayersSsdOut+1;++layer){
    const DCHitContainer &HitCont = GetSsdOutHC(layer);
    int nh = HitCont.size();
    for(int i=0;i<nh;++i){
      DCHit *hit = HitCont[i];
      if(!hit) continue;
      if(!hit->IsGoodWaveForm()) continue;
      double chisqr = hit->GetChisquare();
      if( chisqr>maxchi2 ){
	hit->SetGoodWaveForm(false);
      }
    }
  }
}

//______________________________________________________________________________
void
DCAnalyzer::ChiSqrCutBcOut( double chisqr )
{
  ChiSqrCut( m_BcOutTC, chisqr );
}

//______________________________________________________________________________
void
DCAnalyzer::ChiSqrCutSdcIn( double chisqr )
{
  ChiSqrCut( m_SdcInTC, chisqr );
}

//______________________________________________________________________________
void
DCAnalyzer::ChiSqrCutSdcOut( double chisqr )
{
  ChiSqrCut( m_SdcOutTC, chisqr );
}

//______________________________________________________________________________
void
DCAnalyzer::ChiSqrCutSsdIn( double chisqr )
{
  ChiSqrCut( m_SsdInTC, chisqr );
}

//______________________________________________________________________________
void
DCAnalyzer::ChiSqrCutSsdOut( double chisqr )
{
  ChiSqrCut( m_SsdOutTC, chisqr );
}

//______________________________________________________________________________
void
DCAnalyzer::ChiSqrCut( DCLocalTrackContainer& TrackCont,
		       double chisqr )
{
  DCLocalTrackContainer DeleteCand;
  DCLocalTrackContainer ValidCand;
  int NofTrack = TrackCont.size();
  for(int i = NofTrack-1; i>=0; --i){
    DCLocalTrack* tempTrack = TrackCont.at(i);
    if(tempTrack->GetChiSquare() > chisqr){
      DeleteCand.push_back(tempTrack);
    }else{
      ValidCand.push_back(tempTrack);
    }
  }

  del::ClearContainer( DeleteCand );

  TrackCont.clear();
  TrackCont.resize( ValidCand.size() );
  std::copy( ValidCand.begin(), ValidCand.end(), TrackCont.begin() );
  ValidCand.clear();
}

//______________________________________________________________________________
void
DCAnalyzer::TotCutSDC2(double min_tot)
{
  for(int i = 0; i<NumOfLayersSDC2; ++i){
    TotCut(m_SdcOutHC[i + 1], min_tot, true);
  }// for(i)
}

//______________________________________________________________________________
void
DCAnalyzer::TotCutSDC3(double min_tot)
{
  for(int i = 0; i<NumOfLayersSDC3; ++i){
    TotCut(m_SdcOutHC[i + NumOfLayersSDC2 +1], min_tot, true);
  }// for(i)
}

//______________________________________________________________________________
void
DCAnalyzer::TotCut( DCHitContainer& HitCont,
		    double min_tot, bool adopt_nan )
{
  DCHitContainer ValidCand;
  DCHitContainer DeleteCand;
  for(auto *ptr : HitCont){
    ptr->TotCut(min_tot, adopt_nan);
    if(0 == ptr->GetDriftTimeSize()){
      DeleteCand.push_back(ptr);
    }else{
      ValidCand.push_back(ptr);
    }
  }

  del::ClearContainer( DeleteCand );

  HitCont.clear();
  HitCont.resize( ValidCand.size() );
  std::copy( ValidCand.begin(), ValidCand.end(), HitCont.begin() );
  ValidCand.clear();
}
