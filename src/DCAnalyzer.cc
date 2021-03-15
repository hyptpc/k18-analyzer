// -*- C++ -*-

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
#include "FuncName.hh"
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
#include "UserParamMan.hh"
#include "DeleteUtility.hh"
#include "TPCPadHelper.hh"
#include "TPCPositionCorrector.hh"
#include "TPCRawHit.hh"
#include "TPCHit.hh"
#include "TPCCluster.hh"
#include "TPCLocalTrack.hh"
#include "TPCLocalTrack_Helix.hh"

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
#define SdcIn_XUV         0 // XUV Tracking (not used in KURAMA)
#define SdcIn_Pair        1 // Pair plane Tracking (fast but bad for large angle track)
#define SdcIn_Deletion    1 // Deletion method for too many combinations
/* TPCTracking */
#define UseTpcCluster	1

namespace
{
using namespace K18Parameter;
const auto& gConf   = ConfMan::GetInstance();
const auto& gGeom   = DCGeomMan::GetInstance();
const auto& gTPCPos = TPCPositionCorrector::GetInstance();
const auto& gUser   = UserParamMan::GetInstance();

//_____________________________________________________________________________
const double& pK18 = ConfMan::Get<double>("PK18");
const int& IdTOFUX = gGeom.DetectorId("TOF-UX");
const int& IdTOFUY = gGeom.DetectorId("TOF-UY");
const int& IdTOFDX = gGeom.DetectorId("TOF-DX");
const int& IdTOFDY = gGeom.DetectorId("TOF-DY");

const double MaxChiSqrKuramaTrack = 10000.;
const double MaxTimeDifMWPC       =   100.;

const double kMWPCClusteringWireExtension =  1.0; // [mm]
const double kMWPCClusteringTimeExtension = 10.0; // [nsec]

//_____________________________________________________________________________
inline bool /* for MWPCCluster */
isConnectable( double wire1, double leading1, double trailing1,
               double wire2, double leading2, double trailing2,
               double wExt,  double tExt )
{
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
  hddaq::cout << __func__ << std::endl
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

//_____________________________________________________________________________
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
}

//_____________________________________________________________________________
DCAnalyzer::DCAnalyzer( void )
  : m_is_decoded(n_type),
    m_much_combi(n_type),
    m_MWPCClCont(NumOfLayersBcIn+1),
    m_TempBcInHC(NumOfLayersBcIn+1),
    m_BcInHC(NumOfLayersBcIn+1),
    m_BcOutHC(NumOfLayersBcOut+2),
    m_SdcInHC(NumOfLayersSdcIn+1),
    m_SdcOutHC(NumOfLayersSdcOut+1),
    m_TPCHitCont(NumOfLayersTPC+1),
    m_TempTPCHitCont(NumOfLayersTPC+1),
    m_TPCClCont(NumOfLayersTPC+1),
    m_SdcInExTC(NumOfLayersSdcIn+1),
    m_SdcOutExTC(NumOfLayersSdcOut+1)
{
  for( int i=0; i<n_type; ++i ){
    m_is_decoded[i] = false;
    m_much_combi[i] = 0;
  }
  debug::ObjectCounter::increase(ClassName());
}

DCAnalyzer::~DCAnalyzer( void )
{
  ClearKuramaTracks();
#if UseBcIn
  ClearK18TracksU2D();
  ClearTracksBcIn();
#endif
  ClearK18TracksD2U();
  ClearTracksSdcOut();
  ClearTracksSdcIn();
  ClearTracksBcOut();
  ClearTracksBcOutSdcIn();
  ClearTracksSdcInSdcOut();
  ClearTracksTPC();
  ClearDCHits();
  ClearVtxHits();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
DCAnalyzer::PrintKurama( const std::string& arg ) const
{
  int nn = m_KuramaTC.size();
  hddaq::cout << FUNC_NAME << " " << arg << std::endl
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

//_____________________________________________________________________________
#if UseBcIn
bool
DCAnalyzer::DecodeBcInHits( RawData *rawData )
{
  if( m_is_decoded[k_BcIn] ){
    hddaq::cout << "#D " << FUNC_NAME << " "
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

//_____________________________________________________________________________
bool
DCAnalyzer::DecodeBcOutHits( RawData *rawData )
{
  if( m_is_decoded[k_BcOut] ){
    hddaq::cout << "#D " << FUNC_NAME << " "
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
	m_BcOutHC[layer].push_back(hit);
      else
	delete hit;
    }
  }

  m_is_decoded[k_BcOut] = true;
  return true;
}

//_____________________________________________________________________________
bool
DCAnalyzer::ClusterizeTPC( int layerID, const TPCHitContainer& HitCont,
                           TPCClusterContainer& ClCont )
{

  static const double ClusterYCut = gUser.GetParameter("ClusterYCut");

  del::ClearContainer( ClCont );

  const std::size_t nh = HitCont.size();
  if( nh==0 ) return false;

  std::vector<int> flag( nh, 0 );

  for( std::size_t hiti=0; hiti < nh; hiti++ ) {
    if( flag[hiti] > 0 ) continue;
    TPCHitContainer CandCont;
    TPCHit* hit = HitCont[hiti];
    if( !hit || !hit->IsGood() ) continue;
    CandCont.push_back(hit);
    flag[hiti]++;

    for( std::size_t hitj=0; hitj < nh; hitj++ ) {
      if( hiti==hitj || flag[hitj]>0 ) continue;
      TPCHit* thit = HitCont[hitj];
      if( !thit || !thit->IsGood() ) continue;
      for( int ci=0; ci < CandCont.size(); ci++ ) {
	TPCHit* c_hit = CandCont[ci];
	int rowID = thit->GetRow();
	int c_rowID = c_hit->GetRow();
	// std::cout<<"clusterize TPC1 layer:"<<thit->GetLayer()<<", "
	// 	 <<"row: "<<rowID<<", "
	// 	 <<"de: "<<thit->GetCharge()<<", "
	// 	 <<"pos: "<<thit->GetPos()<<std::endl;
	// std::cout<<"clusterize TPC2 layer:"<<c_hit->GetLayer()<<", "
	// 	 <<"row: "<<c_rowID<<", "
	// 	 <<"de: "<<c_hit->GetCharge()<<", "
	// 	 <<"pos: "<<c_hit->GetPos()<<std::endl;

	if( (abs(rowID - c_rowID) <= 2 ||
	      (layerID<10 && abs(rowID - c_rowID)>=tpc::padParameter[layerID][1]-2) )
	    && fabs( thit->GetY() - c_hit->GetY() ) < ClusterYCut )
	{
	  CandCont.push_back(thit);
	  flag[hitj]++;
	  break;
	}
      }
    }
    TPCCluster* cluster = new TPCCluster( layerID, CandCont);
    // std::cout<<"After clusterize, layer:"<<layerID<<", "
    //  	     <<"pos: "<<cluster->Position()<<", "
    //  	     <<"size:"<<cluster->GetClusterSize()<<std::endl;
    if( cluster ) ClCont.push_back( cluster );
  }

  return true;
}

//_____________________________________________________________________________
bool
DCAnalyzer::DecodeTPCHits( RawData *rawData )
{
  static const Int_t TPC_Subtraction = gUser.GetParameter("TPC_Subtraction");
  static const Int_t TPC_Multi = gUser.GetParameter("TPC_Multi");
  if( m_is_decoded[k_TPC] ){
    hddaq::cout << "#D " << FUNC_NAME << " "
		<< "already decoded" << std::endl;
    return true;
  }

  ClearTPCHits();
  
  int numhit =0; 
    //multiplicity cut for partial readout mode
  for( int layer=0; layer<=NumOfLayersTPC; ++layer ){
    const auto& rhit =  rawData->GetTPCRawHC( layer ); 
    const std::size_t nh = rhit.size();
    numhit += nh;
  }
  if(TPC_Multi>0&&numhit>TPC_Multi){
    m_is_decoded[k_TPC] = true;
    return true;
  }

  for( int layer=0; layer<=NumOfLayersTPC; ++layer ){
    if(TPC_Subtraction == 1 ){
      for( const auto& rhit : rawData->GetTPCCorHC( layer ) ){
	auto hit = new TPCHit( rhit );
	if( hit->DoFit() && hit->Calculate() )
	  m_TPCHitCont[layer].push_back( hit );
	else
	  delete hit;
      }
    }
    else{
     for( const auto& rhit : rawData->GetTPCRawHC( layer ) ){
	auto hit = new TPCHit( rhit );
	if( hit->DoFit() && hit->Calculate() )
	  m_TPCHitCont[layer].push_back( hit );
	else
	  delete hit;
      }
    }

  }

  m_is_decoded[k_TPC] = true;
  return true;
}

//_____________________________________________________________________________
bool
DCAnalyzer::ReCalcTPCHits( const int nhits,
			   const std::vector<int>& padid,
			   const std::vector<double>& time,
                           const std::vector<double>& de,
                           Bool_t do_clusterize )
{
  static const Double_t Time0 = gUser.GetParameter("Time0TPC");
  static const Double_t DriftVelocity = gUser.GetParameter("DriftVelocityTPC");

  if( m_is_decoded[k_TPC] ){
    hddaq::cout << "#D " << FUNC_NAME << " "
		<< "already decoded" << std::endl;
    return true;
  }
  ClearTPCClusters();
  ClearTPCHits();
  //Temporary: these parameter should be given by different param file (ch by ch)

  for( int hiti=0; hiti<nhits; hiti++ ){
    TVector3 pos_tmp = tpc::getPosition(padid[hiti]);
    //    double y = ( time[hiti] - Time0 ) * DriftVelocity;
    //Temporary: DriftVelocity (unit mm/ch)
    double y = ( time[hiti] - Time0 ) * 80. * DriftVelocity;
    TVector3 pos(pos_tmp.x(), y, pos_tmp.z());
    TVector3 cpos = gTPCPos.Correct(pos);
    int layer = tpc::getLayerID( padid[hiti] );
    int row = tpc::getRowID( padid[hiti] );
    // std::cout<<"original padid:"<<padid[hiti]
    // 	     <<", calcpadid:"<<tpc::GetPadId(layer, row)<<std::endl;

    TPCHit* hit = new TPCHit(layer, (double)row);
    hit->SetPos(pos);
    hit->SetCharge(de[hiti]);
    // std::cout<<"Hit, layer:"<<layer<<", "
    // 	     <<"pos:"<<pos<<std::endl;
    if( hit ) m_TPCHitCont[layer].push_back( hit );
  }

  if( do_clusterize ){
    std::vector<TPCClusterContainer>  TPCClusterCont;
    TPCClusterCont.resize(NumOfLayersTPC+1);
    for( int layer=0; layer<=NumOfLayersTPC; ++layer ){
      if(m_TPCHitCont[layer].size()==0)
	continue;

      ClusterizeTPC( layer, m_TPCHitCont[layer], TPCClusterCont[layer] );

      int ncl = TPCClusterCont[layer].size();
      for(int i=0; i<ncl; ++i){
	TPCCluster *p = TPCClusterCont[layer][i];
	TVector3 pos = p->Position();
        TVector3 cpos = gTPCPos.Correct(pos);
        double charge = p->Charge();
	double mrow = p->MeanRow();
	int clusterSize = p->GetClusterSize();
	TPCHit* hit = new TPCHit(layer, mrow);
	hit->SetPos(pos);
	hit->SetCharge(charge);
	hit->SetClusterSize(clusterSize);
        // std::cout<<"Cluster, layer:"<<layer<<", "
        //  	       <<"pos:"<<pos<<", "
        //  	       <<"mrow:"<<p->MeanRow()<<", "
        //  	       <<"cluster size:"<<p->GetClusterSize()<<std::endl;
        // getchar();

        // if( hit->CalcTPCObservables() )
        //  	m_TPCHitCont[layer].push_back(hit);
        // else
        // 	delete hit;
        m_TPCClCont[layer].push_back(hit);
       }
    }
    del::ClearContainerAll( TPCClusterCont );
  }


  m_is_decoded[k_TPC] = true;
  return true;
}

//_____________________________________________________________________________
bool
DCAnalyzer::DecodeTPCHitsGeant4( const int nhits,
			         const double *x, const double *y, const double *z, const double *de )
{
  if( m_is_decoded[k_TPC] ){
    hddaq::cout << "#D " << FUNC_NAME << " "
		<< "already decoded" << std::endl;
    return true;
  }
  ClearTPCClusters();
  ClearTPCHits();

  for( int hiti=0; hiti<nhits; hiti++ ){
    int padid = tpc::findPadID( z[hiti], x[hiti] );
    int layer = tpc::getLayerID( padid );
    int row = tpc::getRowID( padid );
    TVector3 pos(x[hiti], y[hiti], z[hiti]);
    TPCHit  *hit  = new TPCHit(layer,(double)row);
    hit->SetClusterSize(1);
    hit->SetPos(pos);
    hit->SetCharge(de[hiti]);
    m_TPCClCont[layer].push_back(hit);
  }
  // for( int hiti=0; hiti<nhits; hiti++ ){
  //   TPCCluster* cluster = new TPCCluster( x[hiti], y[hiti], z[hiti], de[hiti] );
  //   int layer = tpc::getLayerID( tpc::findPadID( z[hiti], x[hiti] ) );
  //   if( cluster ) m_TPCClCont[layer].push_back( cluster );
  // }

  // for( int layer=0; layer<=NumOfLayersTPC; ++layer ){
  //   int ncl = m_TPCClCont[layer].size();
  //   for(int i=0; i<ncl; ++i){
  //     TPCCluster *p = m_TPCClCont[layer][i];
  //     int MeanPad = p->MeanPadId();
  //     TVector3 pos = p->Position();
  //     double charge = p->Charge();
  //     TPCHit  *hit  = new TPCHit( MeanPad, pos, charge);
  //     hit->SetClusterSize(1);
  //     hit->SetMRow((double)tpc::getRowID(MeanPad));//return row id

  //     // if( hit->CalcTPCObservables() )
  //     //  	m_TPCHitCont[layer].push_back(hit);
  //     // else
  //     // 	delete hit;
  //     m_TPCHitCont[layer].push_back(hit);
  //   }
  // }
  m_is_decoded[k_TPC] = true;
  return true;
}
//_____________________________________________________________________________
bool
DCAnalyzer::DecodeSdcInHits( RawData *rawData )
{
  if( m_is_decoded[k_SdcIn] ){
    hddaq::cout << "#D " << FUNC_NAME << " "
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

      hodoAna.TimeCutSFT( l, -10, 5 );
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
	int       nhtdc      = rhit->GetTdcSize();
	int       nhtrailing = rhit->GetTrailingSize();
	if(!hit) continue;
	for( int j=0; j<nhtdc; ++j ){
	  hit->SetTdcVal( rhit->GetTdc(j) );
	}
	for( int j=0; j<nhtrailing; ++j ){
	  hit->SetTdcTrailing( rhit->GetTrailing(j) );
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

//_____________________________________________________________________________
bool
DCAnalyzer::DecodeSdcOutHits( RawData *rawData , double ofs_dt)
{
  if( m_is_decoded[k_SdcOut] ){
    hddaq::cout << "#D " << FUNC_NAME << " "
		<< "already decoded" << std::endl;
    return true;
  }

  ClearSdcOutHits();

  // SdcOut
  {
    for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
      const DCRHitContainer &RHitCont = rawData->GetSdcOutRawHC(layer);
      int nh = RHitCont.size();
      for( int i=0; i<nh; ++i ){
	DCRawHit *rhit  = RHitCont[i];
	DCHit    *hit   = new DCHit( rhit->PlaneId(), rhit->WireId() );
	int       nhtdc      = rhit->GetTdcSize();
	int       nhtrailing = rhit->GetTrailingSize();
	if(!hit) continue;

	hit->SetOfsdT(ofs_dt);
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
  }
 //*
  // FBT1
  {
    HodoAnalyzer hodoAna;
    hodoAna.DecodeFBT1Hits( rawData );

    for ( int l = 0; l < 2; ++l ) {
      hodoAna.TimeCutFBT1( l, 0, -10, 5 ); //FBT1-U
      hodoAna.TimeCutFBT1( l, 1, -10, 5 ); //FBT1-D

      int nclU = hodoAna.GetNClustersFBT1( l, 0 );
      int nclD = hodoAna.GetNClustersFBT1( l, 1 );

      for ( int j = 0; j < nclU; ++j ) {
        FiberCluster* clU = hodoAna.GetClusterFBT1( l, 0, j );
        double segU  = clU->MeanSeg();
        double posU  = clU->MeanPosition();
        double timeU = clU->CMeanTime();

        DCHit *hitU = new DCHit( 1+2*l+PlOffsFht, segU );
        hitU->SetTdcVal( static_cast<int>( timeU ) );

        int layer = NumOfLayersSdcOut-7 + 1 + 2*l;
        if ( hitU->CalcFiberObservables() ) {
          hitU->SetWirePosition( posU );
          m_SdcOutHC[layer].push_back( hitU );
	} else {
          delete hitU;
        }
      }

      for ( int j = 0; j < nclD; ++j ) {
        FiberCluster* clD = hodoAna.GetClusterFBT1( l, 1, j );
	double segD  = clD->MeanSeg();
        double posD  = clD->MeanPosition();
        double timeD = clD->CMeanTime();

        DCHit *hitD = new DCHit( 2*l+PlOffsFht, segD );
        hitD->SetTdcVal( static_cast<int>( timeD ) );

        int layer = NumOfLayersSdcOut-7 + 2*l;
        if ( hitD->CalcFiberObservables() ) {
          hitD->SetWirePosition( posD );
          m_SdcOutHC[layer].push_back( hitD );
	} else {
          delete hitD;
        }
      }
    }
  }

  // FBT2
  {
    HodoAnalyzer hodoAna;
    hodoAna.DecodeFBT2Hits( rawData );

    for ( int l = 0; l < 2; ++l ) {
      hodoAna.TimeCutFBT2( l, 0, -10, 5 ); //FBT2-U
      hodoAna.TimeCutFBT2( l, 1, -10, 5 ); //FBT2-D

      int nclU = hodoAna.GetNClustersFBT2( l, 0 );
      int nclD = hodoAna.GetNClustersFBT2( l, 1 );

      for ( int j = 0; j < nclU; ++j ) {
        FiberCluster* clU = hodoAna.GetClusterFBT2( l, 0, j );
        double segU  = clU->MeanSeg();
        double posU  = clU->MeanPosition();
        double timeU = clU->CMeanTime();

        DCHit *hitU = new DCHit( 5+2*l+PlOffsFht, segU );
	hitU->SetTdcVal( timeU );

        int layer = NumOfLayersSdcOut-3 + 1 + 2*l;
        if ( hitU->CalcFiberObservables() ) {
          hitU->SetWirePosition( posU );
          m_SdcOutHC[layer].push_back( hitU );
        } else {
          delete hitU;
	}
      }

      for ( int j = 0; j < nclD; ++j ) {
        FiberCluster* clD = hodoAna.GetClusterFBT2( l, 1, j );
        double segD  = clD->MeanSeg();
	double posD  = clD->MeanPosition();
        double timeD = clD->CMeanTime();

        DCHit *hitD = new DCHit( 4+2*l+PlOffsFht, segD );
        hitD->SetTdcVal( timeD );

	int layer = NumOfLayersSdcOut-3 + 2*l;
        if ( hitD->CalcFiberObservables() ) {
          hitD->SetWirePosition( posD );
          m_SdcOutHC[layer].push_back( hitD );
        } else {
          delete hitD;
        }
      }
    }
  }

  m_is_decoded[k_SdcOut] = true;
  return true;
}

//_____________________________________________________________________________
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
  return true;
}

//_____________________________________________________________________________
bool
DCAnalyzer::DecodeTOFHits( const Hodo2HitContainer& HitCont )
{
  if( m_is_decoded[k_TOF] ){
    hddaq::cout << "#D " << FUNC_NAME << " "
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

//_____________________________________________________________________________
bool
DCAnalyzer::DecodeTOFHits( const HodoClusterContainer& ClCont )
{
  if( m_is_decoded[k_TOF] ){
    hddaq::cout << "#D " << FUNC_NAME << " "
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

//_____________________________________________________________________________
#if UseBcIn
bool
DCAnalyzer::TrackSearchBcIn( void )
{
  track::MWPCLocalTrackSearch( &(m_BcInHC[1]), m_BcInTC );
  return true;
}

//_____________________________________________________________________________
bool
DCAnalyzer::TrackSearchBcIn( const std::vector<std::vector<DCHitContainer> >& hc )
{
  track::MWPCLocalTrackSearch( hc, m_BcInTC );
  return true;
}
#endif

//_____________________________________________________________________________
bool
DCAnalyzer::TrackSearchBcOut( int T0Seg )
{
  static const int MinLayer = gUser.GetParameter("MinLayerBcOut");

#if BcOut_Pair //Pair Plane Tracking Routine for BcOut
  int ntrack = track::LocalTrackSearch( m_BcOutHC, PPInfoBcOut, NPPInfoBcOut,
					m_BcOutTC, MinLayer, T0Seg );
  return ntrack == -1 ? false : true;
#endif

#if BcOut_XUV  //XUV Tracking Routine for BcOut
  int ntrack = track::LocalTrackSearchVUX( m_BcOutHC, PPInfoBcOut, NPPInfoBcOut,
					   m_BcOutTC, MinLayer );
  return ntrack == -1 ? false : true;
#endif

  return false;
}

//_____________________________________________________________________________
// Use with BH2Filter
bool
DCAnalyzer::TrackSearchBcOut( const std::vector<std::vector<DCHitContainer> >& hc, int T0Seg )
{
  static const int MinLayer = gUser.GetParameter("MinLayerBcOut");

#if BcOut_Pair //Pair Plane Tracking Routine for BcOut
  int ntrack = track::LocalTrackSearch( hc, PPInfoBcOut, NPPInfoBcOut, m_BcOutTC, MinLayer, T0Seg );
  return ntrack == -1 ? false : true;
#endif

#if BcOut_XUV  //XUV Tracking Routine for BcOut
  int ntrack = track::LocalTrackSearchVUX( hc, PPInfoBcOut, NPPInfoBcOut, m_BcOutTC, MinLayer );
  return ntrack == -1 ? false : true;
#endif

  return false;
}

//_____________________________________________________________________________
bool
DCAnalyzer::TrackSearchSdcIn( void )
{
  static const int MinLayer = gUser.GetParameter("MinLayerSdcIn");

  // track::LocalTrackSearch( m_SdcInHC, PPInfoSdcIn, NPPInfoSdcIn, m_SdcInTC, MinLayer );
  track::LocalTrackSearchSdcInFiber( m_SdcInHC, PPInfoSdcIn, NPPInfoSdcIn, m_SdcInTC, MinLayer );
  return true;
}

//_____________________________________________________________________________
bool
DCAnalyzer::TrackSearchSdcOut( void )
{
  static const int MinLayer = gUser.GetParameter("MinLayerSdcOut");

  track::LocalTrackSearchSdcOut( m_SdcOutHC, PPInfoSdcOut, NPPInfoSdcOut,
				 m_SdcOutTC, MinLayer );

  return true;
}

//_____________________________________________________________________________
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

//_____________________________________________________________________________
bool
DCAnalyzer::TrackSearchSdcOut( const HodoClusterContainer& TOFCont )
{
  static const int MinLayer = gUser.GetParameter("MinLayerSdcOut");

  if( !DecodeTOFHits( TOFCont ) ) return false;

  track::LocalTrackSearchSdcOut( m_TOFHC, m_SdcOutHC, PPInfoSdcOut, NPPInfoSdcOut+2,
				 m_SdcOutTC, MinLayer );

  return true;
}

//_____________________________________________________________________________
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

//_____________________________________________________________________________
#if UseBcIn
bool
DCAnalyzer::TrackSearchK18U2D( void )
{
  ClearK18TracksU2D();

  int nIn  = m_BcInTC.size();
  int nOut = m_BcOutTC.size();

# if 0
  hddaq::cout << "**************************************" << std::endl;
  hddaq::cout << FUNC_NAME << ": #TracksIn=" << std::setw(3) << nIn
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
      hddaq::cout << FUNC_NAME << ": In -> " << trIn->GetChiSquare()
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
    hddaq::cout << FUNC_NAME << ": Before sorting. #Track="
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

//_____________________________________________________________________________
bool
DCAnalyzer::TrackSearchK18D2U( const std::vector<double>& XinCont )
{
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
      hddaq::cout << FUNC_NAME
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
    hddaq::cout << FUNC_NAME << ": Before sorting. #Track="
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

//_____________________________________________________________________________
bool
DCAnalyzer::TrackSearchKurama( void )
{
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
	trKurama->SetInitialMomentum( 1. );
      }
      if( trKurama->DoFit() && trKurama->chisqr()<MaxChiSqrKuramaTrack ){
	m_KuramaTC.push_back( trKurama );
      }
      else{
	//	trKurama->Print( "in "+FUNC_NAME );
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

//_____________________________________________________________________________
bool
DCAnalyzer::TrackSearchKurama( double initial_momentum )
{
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
	trKurama->Print( " in "+FUNC_NAME );
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

//_____________________________________________________________________________
bool
DCAnalyzer::TrackSearchTPC( void )
{
  static const int MinLayer = gUser.GetParameter("MinLayerTPC");

#if UseTpcCluster
  track::LocalTrackSearchTPC(m_TPCClCont, m_TPCTC, MinLayer );
#else
  track::LocalTrackSearchTPC(m_TPCHitCont, m_TPCTC, MinLayer );
#endif
  return true;
}

//_____________________________________________________________________________
bool
DCAnalyzer::TrackSearchTPC_Helix( void )
{
  static const int MinLayer = gUser.GetParameter("MinLayerTPC");

#if UseTpcCluster
  track::LocalTrackSearchTPC_Helix(m_TPCClCont, m_TPCTC_Helix, MinLayer );
#else
  track::LocalTrackSearchTPC_Helix(m_TPCHitCont, m_TPCTC_Helix, MinLayer );
#endif
  return true;
}


//_____________________________________________________________________________
void
DCAnalyzer::ClearDCHits( void )
{
#if UseBcIn
  ClearBcInHits();
#endif
  ClearBcOutHits();
  ClearSdcInHits();
  ClearSdcOutHits();
  ClearTOFHits();
  ClearTPCHits();
  ClearTPCClusters();
}

//_____________________________________________________________________________
#if UseBcIn
void
DCAnalyzer::ClearBcInHits( void )
{
  del::ClearContainerAll( m_TempBcInHC );
  del::ClearContainerAll( m_BcInHC );
  del::ClearContainerAll( m_MWPCClCont );
}
#endif

//_____________________________________________________________________________
void
DCAnalyzer::ClearBcOutHits( void )
{
  del::ClearContainerAll( m_BcOutHC );
}


//_____________________________________________________________________________
void
DCAnalyzer::ClearSdcInHits( void )
{
  del::ClearContainerAll( m_SdcInHC );
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearSdcOutHits( void )
{
  del::ClearContainerAll( m_SdcOutHC );
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearVtxHits( void )
{
  del::ClearContainer( m_VtxPoint );
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTOFHits( void )
{
  del::ClearContainer( m_TOFHC );
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTPCHits( void )
{
  del::ClearContainerAll( m_TPCHitCont );
  del::ClearContainerAll( m_TempTPCHitCont );

}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTPCClusters( void )
{
  del::ClearContainerAll( m_TPCClCont );
}

//_____________________________________________________________________________
#if UseBcIn
void
DCAnalyzer::ClearTracksBcIn( void )
{
  del::ClearContainer( m_BcInTC );
}
#endif

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksBcOut( void )
{
  del::ClearContainer( m_BcOutTC );
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksSdcIn( void )
{
  del::ClearContainer( m_SdcInTC );
  del::ClearContainerAll( m_SdcInExTC );
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksSdcOut( void )
{
  del::ClearContainer( m_SdcOutTC );
  del::ClearContainerAll( m_SdcOutExTC );
}

//_____________________________________________________________________________
#if UseBcIn
void
DCAnalyzer::ClearK18TracksU2D( void )
{
  del::ClearContainer( m_K18U2DTC );
}
#endif

//_____________________________________________________________________________
void
DCAnalyzer::ClearK18TracksD2U( void )
{
  del::ClearContainer( m_K18D2UTC );
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearKuramaTracks( void )
{
  del::ClearContainer( m_KuramaTC );
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksBcOutSdcIn( void )
{
  del::ClearContainer( m_BcOutSdcInTC );
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksSdcInSdcOut( void )
{
  del::ClearContainer( m_SdcInSdcOutTC );
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksTPC( void )
{
  del::ClearContainer( m_TPCTC );
  del::ClearContainer( m_TPCTC_Helix );
}


//_____________________________________________________________________________
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

//_____________________________________________________________________________
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

//_____________________________________________________________________________
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

  return true;
}

//_____________________________________________________________________________
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

//_____________________________________________________________________________
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

//_____________________________________________________________________________
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

//_____________________________________________________________________________
#if UseBcIn
bool
DCAnalyzer::ReCalcTrackBcIn( bool applyRecursively )
{
  return ReCalcTrack( m_BcInTC );
}
#endif

//_____________________________________________________________________________
bool
DCAnalyzer::ReCalcTrackBcOut( bool applyRecursively )
{
  return ReCalcTrack( m_BcOutTC );
}

//_____________________________________________________________________________
bool
DCAnalyzer::ReCalcTrackSdcIn( bool applyRecursively )
{
  return ReCalcTrack( m_SdcInTC );
}

//_____________________________________________________________________________
bool
DCAnalyzer::ReCalcTrackSdcOut( bool applyRecursively )
{
  return ReCalcTrack( m_SdcOutTC );
}

//_____________________________________________________________________________
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

//_____________________________________________________________________________
bool
DCAnalyzer::ReCalcK18TrackD2U( bool applyRecursively )
{
  return ReCalcTrack( m_K18D2UTC, applyRecursively );
}

//_____________________________________________________________________________
bool
DCAnalyzer::ReCalcKuramaTrack( bool applyRecursively )
{
  return ReCalcTrack( m_KuramaTC, applyRecursively );
}

//_____________________________________________________________________________
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

  //ReCalcK18TrackD2U();
  ReCalcKuramaTrack();

  return true;
}

//_____________________________________________________________________________
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

//_____________________________________________________________________________
void
DCAnalyzer::ChiSqrCutBcOut( double chisqr )
{
  ChiSqrCut( m_BcOutTC, chisqr );
}

//_____________________________________________________________________________
void
DCAnalyzer::ChiSqrCutSdcIn( double chisqr )
{
  ChiSqrCut( m_SdcInTC, chisqr );
}

//_____________________________________________________________________________
void
DCAnalyzer::ChiSqrCutSdcOut( double chisqr )
{
  ChiSqrCut( m_SdcOutTC, chisqr );
}

//_____________________________________________________________________________
void
DCAnalyzer::ChiSqrCut( DCLocalTrackContainer& TrackCont,
		       double chisqr )
{
  DCLocalTrackContainer DeleteCand;
  DCLocalTrackContainer ValidCand;
  for(auto& tempTrack : TrackCont){
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

//_____________________________________________________________________________
void
DCAnalyzer::TotCutBCOut(double min_tot)
{
  for(int i = 0; i<NumOfLayersBcOut; ++i){
    TotCut(m_BcOutHC[i + 1], min_tot, false);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCutSDC1(double min_tot)
{
  for(int i = 0; i<NumOfLayersSdcIn - NumOfLayersSFT; ++i){
    TotCut(m_SdcInHC[i + NumOfLayersSFT +1], min_tot, false);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCutSDC2(double min_tot)
{
  for(int i = 0; i<NumOfLayersSDC2; ++i){
    TotCut(m_SdcOutHC[i + 1], min_tot, false);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCutSDC3(double min_tot)
{
  for(int i = 0; i<NumOfLayersSDC3; ++i){
    TotCut(m_SdcOutHC[i + NumOfLayersSDC2 +1], min_tot, false);
  }// for(i)
}

//_____________________________________________________________________________
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

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCutBC34(double min_dt, double max_dt)
{
  for(int i = 0; i<NumOfLayersBcOut; ++i){
    DriftTimeCut(m_BcOutHC[i + 1], min_dt, max_dt, true);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCutSDC2(double min_dt, double max_dt)
{
  for(int i = 0; i<NumOfLayersSDC2; ++i){
    DriftTimeCut(m_SdcOutHC[i + 1], min_dt, max_dt, true);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCutSDC3(double min_dt, double max_dt)
{
  for(int i = 0; i<NumOfLayersSDC3; ++i){
    DriftTimeCut(m_SdcOutHC[i + NumOfLayersSDC2 + 1], min_dt, max_dt, true);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCut( DCHitContainer& HitCont,
			  double min_dt, double max_dt, bool select_1st )
{
  DCHitContainer ValidCand;
  DCHitContainer DeleteCand;
  for(auto *ptr : HitCont){
    ptr->GateDriftTime(min_dt, max_dt, select_1st);
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

//_____________________________________________________________________________
bool
DCAnalyzer::MakeBH2DCHit(int t0seg)
{
  static const double centerbh2[] = {
    -41.8, -19.3, -10.7, -3.6, 3.6, 10.7, 19.3, 41.8
  };

  bool status = true;

  double bh2pos = centerbh2[t0seg];
  DCHit *dchit = new DCHit(125, t0seg);
  dchit->SetTdcVal(0.);
  if(dchit->CalcFiberObservables()){
    dchit->SetWirePosition(bh2pos);
    m_BcOutHC[13].push_back(dchit);
  }else{
    delete dchit;
    status = false;
  }

  return status;
}
