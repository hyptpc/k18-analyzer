/**
 *  file: HodoAnalyzer.cc
 *  date: 2017.04.10
 *
 */

#include "HodoAnalyzer.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "FLHit.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoCluster.hh"
#include "RawData.hh"

namespace
{
  const std::string& class_name("HodoAnalyzer");
  const double MaxTimeDifBH1 =  2.0;
  const double MaxTimeDifBH2 =  2.0;
  const double MaxTimeDifTOF = -1.0;
  const double MaxTimeDifLC  =  3.5;
  const double MaxTimeDifBFT =  8.0;
  const double MaxTimeDifSFT =  8.0;
  const double MaxTimeDifSCH = 10.0;
  const double MaxTimeDifFBH = 20.0;
}

#define Cluster 1

//______________________________________________________________________________
HodoAnalyzer::HodoAnalyzer( void )
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
HodoAnalyzer::~HodoAnalyzer( void )
{
  ClearBH1Hits();
  ClearBH2Hits();
  ClearBACHits();
  ClearSACHits();
  ClearPVACHits();
  ClearFACHits();
  ClearSSDTHits();
  ClearTOFHits();
  ClearLCHits();
  ClearBFTHits();
  ClearSFTHits();
  ClearSCHHits();
  ClearFBHHits();
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void
HodoAnalyzer::ClearBH1Hits( void )
{
  del::ClearContainer( m_BH1Cont );
  del::ClearContainer( m_BH1ClCont );
}

//______________________________________________________________________________
void
HodoAnalyzer::ClearBH2Hits( void )
{
  del::ClearContainer( m_BH2Cont );
  del::ClearContainer( m_BH2ClCont );
}

//______________________________________________________________________________
void
HodoAnalyzer::ClearBACHits( void )
{
  del::ClearContainer( m_BACCont );
}

//______________________________________________________________________________
void
HodoAnalyzer::ClearPVACHits( void )
{
  del::ClearContainer( m_PVACCont );
}

//______________________________________________________________________________
void
HodoAnalyzer::ClearSACHits( void )
{
  del::ClearContainer( m_SACCont );
}

//______________________________________________________________________________
void
HodoAnalyzer::ClearFACHits( void )
{
  del::ClearContainer( m_FACCont );
}

//______________________________________________________________________________
void
HodoAnalyzer::ClearSSDTHits( void )
{
  del::ClearContainer( m_SSDTCont );
}

//______________________________________________________________________________
void
HodoAnalyzer::ClearTOFHits( void )
{
  del::ClearContainer( m_TOFCont );
  del::ClearContainer( m_TOFClCont );
}

//______________________________________________________________________________
void
HodoAnalyzer::ClearLCHits( void )
{
  del::ClearContainer( m_LCCont );
  del::ClearContainer( m_LCClCont );
}

//______________________________________________________________________________
void
HodoAnalyzer::ClearBFTHits( void )
{
  for ( MultiPlaneFiberHitContainer::iterator itr = m_BFTCont.begin();
  	itr != m_BFTCont.end(); ++itr ) {
    del::ClearContainer( *itr );
  }
  del::ClearContainer( m_BFTClCont );
}

//______________________________________________________________________________
void
HodoAnalyzer::ClearSFTHits( void )
{
  for ( MultiPlaneFiberHitContainer::iterator itr = m_SFTCont.begin();
  	itr != m_SFTCont.end(); ++itr ) {
    del::ClearContainer( *itr );
  }
  for ( MultiPlaneFiberClusterContainer::iterator itr = m_SFTClCont.begin();
  	itr != m_SFTClCont.end(); ++itr ) {
    del::ClearContainer( *itr );
  }
}

//______________________________________________________________________________
void
HodoAnalyzer::ClearSCHHits( void )
{
  del::ClearContainer( m_SCHCont );
  del::ClearContainer( m_SCHClCont );
}

//______________________________________________________________________________
void
HodoAnalyzer::ClearFBHHits( void )
{
  del::ClearContainer( m_FBHCont );
  // FLHit should be deleted by FiberHit
  del::ClearContainer( m_FBHCoinCont, false );
  del::ClearContainer( m_FBHClCont );
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeRawHits( RawData *rawData )
{
  DecodeBH1Hits( rawData );
  DecodeBH2Hits( rawData );
  DecodeBACHits( rawData );
  DecodePVACHits( rawData );
  DecodeSACHits( rawData );
  DecodeFACHits( rawData );
  DecodeSSDTHits( rawData );
  DecodeTOFHits( rawData );
  DecodeLCHits( rawData );
  DecodeBFTHits( rawData );
  DecodeSFTHits( rawData );
  DecodeSCHHits( rawData );
  DecodeFBHHits( rawData );
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeBH1Hits( RawData *rawData )
{
  ClearBH1Hits();
  const HodoRHitContainer &cont = rawData->GetBH1RawHC();
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    if( hit->GetTdcUp()<=0 || hit->GetTdcDown()<=0 ) continue;
    Hodo2Hit *hp = new Hodo2Hit( hit );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_BH1Cont.push_back(hp);
    else
      delete hp;
  }//for(i)

#if Cluster
  MakeUpClusters( m_BH1Cont, m_BH1ClCont, MaxTimeDifBH1 );
#endif

  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeBH2Hits( RawData *rawData )
{
  ClearBH2Hits();
  const HodoRHitContainer &cont = rawData->GetBH2RawHC();
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    if( hit->GetTdcUp()<=0 || hit->GetTdcDown()<=0 ) continue;
    BH2Hit *hp = new BH2Hit( hit );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_BH2Cont.push_back(hp);
    else
      delete hp;
  }//for(i)

#if Cluster
  MakeUpClusters( m_BH2Cont, m_BH2ClCont, MaxTimeDifBH2 );
#endif

  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeBACHits( RawData *rawData )
{
  ClearBACHits();
  const HodoRHitContainer &cont = rawData->GetBACRawHC();
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    Hodo1Hit *hp = new Hodo1Hit( hit );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_BACCont.push_back(hp);
    else
      delete hp;
  }//for(i)
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodePVACHits( RawData *rawData )
{
  ClearPVACHits();
  const HodoRHitContainer &cont = rawData->GetPVACRawHC();
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    Hodo1Hit *hp = new Hodo1Hit( hit );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_PVACCont.push_back(hp);
    else
      delete hp;
  }//for(i)
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeSACHits( RawData *rawData )
{
  ClearSACHits();
  const HodoRHitContainer &cont = rawData->GetSACRawHC();
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    Hodo1Hit *hp = new Hodo1Hit( hit );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_SACCont.push_back(hp);
    else
      delete hp;
  }//for(i)
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeFACHits( RawData *rawData )
{
  ClearFACHits();
  const HodoRHitContainer &cont = rawData->GetFACRawHC();
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    Hodo1Hit *hp = new Hodo1Hit( hit );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_FACCont.push_back(hp);
    else
      delete hp;
  }//for(i)
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeSSDTHits( RawData *rawData )
{
  ClearSSDTHits();
  const HodoRHitContainer &cont = rawData->GetSSDTRawHC();
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    Hodo1Hit *hp = new Hodo1Hit( hit );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_SSDTCont.push_back(hp);
    else
      delete hp;
  }//for(i)
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeTOFHits( RawData *rawData )
{
  ClearTOFHits();
  const HodoRHitContainer &cont = rawData->GetTOFRawHC();
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    if( hit->GetTdcUp()<=0 || hit->GetTdcDown()<=0 ) continue;
    Hodo2Hit *hp = new Hodo2Hit( hit );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_TOFCont.push_back(hp);
    else
      delete hp;
  }//for(i)

#if Cluster
  MakeUpClusters( m_TOFCont, m_TOFClCont, MaxTimeDifTOF );
#endif

  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeLCHits( RawData *rawData )
{
  ClearLCHits();
  const HodoRHitContainer &cont = rawData->GetLCRawHC();
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    Hodo1Hit *hp = new Hodo1Hit( hit );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_LCCont.push_back(hp);
    else
      delete hp;
  }//for(i)

  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeBFTHits( RawData* rawData )
{
  ClearBFTHits();

  m_BFTCont.resize( NumOfPlaneBFT );
  for(int p = 0; p<NumOfPlaneBFT; ++p){
    const HodoRHitContainer &cont = rawData->GetBFTRawHC(p);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      if( !hit ) continue;
      FiberHit *hp = new FiberHit(hit, "BFT");
      if(!hp) continue;
      if(hp->Calculate()){
	m_BFTCont.at(p).push_back(hp);
      }
      else{
	delete hp;
	hp = NULL;
      }
    }

    std::sort(m_BFTCont.at(p).begin(), m_BFTCont.at(p).end(), FiberHit::CompFiberHit);
  }


#if Cluster
  FiberHitContainer cont_merge( m_BFTCont.at( 0 ) );
  cont_merge.reserve( m_BFTCont.at( 0 ).size() + m_BFTCont.at( 1 ).size() );
  cont_merge.insert( cont_merge.end(), m_BFTCont.at( 1 ).begin(), m_BFTCont.at( 1 ).end() );
  std::sort( cont_merge.begin(), cont_merge.end(), FiberHit::CompFiberHit );
  MakeUpClusters( cont_merge, m_BFTClCont, MaxTimeDifBFT, 3 );
#endif

  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeSFTHits( RawData* rawData )
{
  ClearSFTHits();

  m_SFTCont.resize( NumOfPlaneSFT );
  m_SFTClCont.resize( NumOfLayersSFT );

  for(int p = 0; p<NumOfPlaneSFT; ++p){

    const HodoRHitContainer &cont = rawData->GetSFTRawHC(p);

    int nh=cont.size();
    for( int i=0; i<nh; ++i ){

      HodoRawHit *hit=cont[i];
      if( !hit ) continue;

      char pname[8];
      switch ( p ) {
      case 0:
	strcpy( pname, "SFT-U" ); break;
      case 1:
	strcpy( pname, "SFT-V" ); break;
      case 2: case 3:
	strcpy( pname, "SFT-X" ); break;
      default: break;
      }

      FiberHit *hp = new FiberHit( hit, pname );
      if(!hp) continue;
      if(hp->Calculate()){
	m_SFTCont.at( p ).push_back( hp );
      }
      else{
	delete hp;
	hp = NULL;
      }

    }
    std::sort( m_SFTCont.at( p ).begin(), m_SFTCont.at( p ).end(), FiberHit::CompFiberHit );
  }

#if Cluster
  FiberHitContainer cont_merge( m_SFTCont.at( 2 ) );
  cont_merge.reserve( m_SFTCont.at( 2 ).size() + m_SFTCont.at( 3 ).size() );
  cont_merge.insert( cont_merge.end(), m_SFTCont.at( 3 ).begin(), m_SFTCont.at( 3 ).end() );
  std::sort( cont_merge.begin(), cont_merge.end(), FiberHit::CompFiberHit );
  for ( int l = 0; l < NumOfLayersSFT; ++l ) {
    FiberHitContainer* p_cont;
    switch ( l ) {
    case 2:
      p_cont = &cont_merge; break;
    default:
      p_cont = &(m_SFTCont.at( l )); break;
    }
    MakeUpClusters( *p_cont, m_SFTClCont.at( l ), MaxTimeDifSFT, 3 );
  }

#endif
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeSCHHits( RawData* rawData )
{
  ClearSCHHits();
  const HodoRHitContainer &cont = rawData->GetSCHRawHC();
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    FiberHit *hp = new FiberHit( hit, "SCH" );
    if( !hp ) continue;
    if( hp->Calculate() ){
      m_SCHCont.push_back( hp );
    }else{
      delete hp;
    }
  }

  std::sort( m_SCHCont.begin(), m_SCHCont.end(), FiberHit::CompFiberHit );

#if Cluster
  MakeUpClusters( m_SCHCont, m_SCHClCont, MaxTimeDifSCH, 1 );
#endif

  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeFBHHits( RawData* rawData )
{
  ClearFBHHits();
  const HodoRHitContainer &cont = rawData->GetFBHRawHC();
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    FiberHit *hp = new FiberHit( hit, "FBH" );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_FBHCont.push_back( hp );
    else
      delete hp;
  }///for(i)

  std::sort( m_FBHCont.begin(), m_FBHCont.end(), FiberHit::CompFiberHit );

#if Cluster
    MakeUpCoincidence( m_FBHCont, m_FBHCoinCont, 15 );
    MakeUpClusters( m_FBHCoinCont, m_FBHClCont, MaxTimeDifFBH, 1);
#endif

  return true;
}

//______________________________________________________________________________
int
HodoAnalyzer::MakeUpClusters( const Hodo2HitContainer& HitCont,
			      HodoClusterContainer& ClusterCont,
			      double maxTimeDif )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  del::ClearContainer( ClusterCont );

  int nh = HitCont.size();

  //  std::vector <int> flag( nh, 0 );

  for(int i=0; i<nh; ++i){
    Hodo2Hit *hitA = HitCont[i];
    int       segA = hitA->SegmentId();
    Hodo2Hit *hitB = 0;
    int         iB = -1;
    int         mB = -1;
    int         mC = -1;
    int       segB = -1;
    if(hitA->JoinedAllMhit()) continue;
    
    int n_mhitA = hitA->GetNumOfHit();
    for(int ma = 0; ma<n_mhitA; ++ma){
      if(hitA->Joined(ma)) continue;
      double    cmtA = hitA->CMeanTime(ma);
      double    cmtB = -9999;
      
      for( int j=i+1; j<nh; ++j ){
	Hodo2Hit *hit = HitCont[j];
	int    seg  = hit->SegmentId();

	int n_mhitB = hit->GetNumOfHit();
	for(int mb = 0; mb<n_mhitB; ++mb){
	  if(hit->Joined(mb)) continue;
	  double cmt = hit->CMeanTime(mb);
	  if( std::abs(seg-segA)==1 && std::abs(cmt-cmtA)<maxTimeDif ){
	    hitB = hit; iB = j; mB = mb; segB = seg; cmtB = cmt; break;
	  }
	}// for(mb:hitB)
      }// for(j:hitB)
      if( hitB ){
	Hodo2Hit *hitC = 0;
	for( int j=i+1; j<nh; ++j ){
	  if( j==iB ) continue;
	  Hodo2Hit *hit = HitCont[j];
	  int    seg  = hit->SegmentId();
	  
	  int n_mhitC = hit->GetNumOfHit();
	  for(int mc = 0; mc<n_mhitC; ++mc){
	    if(hit->Joined(mc)) continue;
	    double cmt = hit->CMeanTime(mc);
	    if( (std::abs(seg-segA)==1 && std::abs(cmt-cmtA)<maxTimeDif) ||
		(std::abs(seg-segB)==1 && std::abs(cmt-cmtB)<maxTimeDif) ){
	      hitC=hit; mC=mc; break;
	    }
	  }
	}
	if( hitC ){
	  hitA->SetJoined(ma);
	  hitB->SetJoined(mB);
	  hitC->SetJoined(mC);
	  HodoCluster *cluster = new HodoCluster( hitA, hitB, hitC );
	  cluster->SetIndex(ma, mB, mC);
	  cluster->Calculate();
	  if( cluster ) ClusterCont.push_back( cluster );
	}
	else{
	  hitA->SetJoined(ma);
	  hitB->SetJoined(mB);
	  HodoCluster *cluster = new HodoCluster( hitA, hitB );
	  cluster->SetIndex(ma, mB);
	  cluster->Calculate();
	  if( cluster ) ClusterCont.push_back( cluster );
	}
      }
      else{
	hitA->SetJoined(ma);
	HodoCluster *cluster = new HodoCluster( hitA );
	cluster->SetIndex(ma);
	cluster->Calculate();
	if( cluster ) ClusterCont.push_back( cluster );
      }
    }// for(ma:hitA)
  }// for(i:hitA)

  return ClusterCont.size();
}

//______________________________________________________________________________
int
HodoAnalyzer::MakeUpClusters( const BH2HitContainer& HitCont,
			      BH2ClusterContainer& ClusterCont,
			      double maxTimeDif )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  del::ClearContainer( ClusterCont );

  int nh = HitCont.size();

  for(int i=0; i<nh; ++i ){
    BH2Hit *hitA = HitCont[i];
    int     segA = hitA->SegmentId();
    BH2Hit *hitB = 0;
    int       iB = -1;
    int       mB = -1;
    int       mC = -1;
    int     segB = -1;
    if(hitA->JoinedAllMhit()) continue;

    int n_mhitA = hitA->GetNumOfHit();    
    for(int ma = 0; ma<n_mhitA; ++ma){
      if(hitA->Joined(ma)) continue;
      double    cmtA = hitA->CMeanTime(ma);
      double    cmtB = -9999;    

      for( int j=i+1; j<nh; ++j ){
	BH2Hit *hit = HitCont[j];
	int     seg = hit->SegmentId();
	
	int n_mhitB = hit->GetNumOfHit();
	for(int mb = 0; mb<n_mhitB; ++mb){
	  if(hit->Joined(mb)) continue;
	  double  cmt = hit->CMeanTime(mb);
	  if( std::abs(seg-segA)==1 && std::abs(cmt-cmtA)<maxTimeDif ){
	    hitB = hit; iB = j; mB = mb; segB = seg; cmtB = cmt; break;
	  }
	}// for(mb:hitB)
      }// for(j:hitB)
      if(hitB){
	BH2Hit *hitC = 0;
	for( int j=i+1; j<nh; ++j ){
	  if( j==iB ) continue;
	  BH2Hit *hit = HitCont[j];
	  int     seg = hit->SegmentId();

	  int n_mhitC = hit->GetNumOfHit();
	  for(int mc = 0; mc<n_mhitC; ++mc){
	    if(hit->Joined(mc)) continue;
	    double  cmt = hit->CMeanTime(mc);
	    if( (std::abs(seg-segA)==1 && std::abs(cmt-cmtA)<maxTimeDif) ||
		(std::abs(seg-segB)==1 && std::abs(cmt-cmtB)<maxTimeDif) ){
	      hitC = hit; mC=mc; break;
	    }
	  }
	}
	if( hitC ){
	  hitA->SetJoined(ma);
	  hitB->SetJoined(mB);
	  hitC->SetJoined(mC);
	  BH2Cluster *cluster = new BH2Cluster( hitA, hitB, hitC );
	  cluster->SetIndex(ma, mB, mC);
	  cluster->Calculate();
	  if( cluster ) ClusterCont.push_back( cluster );
	}
	else{
	  hitA->SetJoined(ma);
	  hitB->SetJoined(mB);
	  BH2Cluster *cluster = new BH2Cluster( hitA, hitB );
	  cluster->SetIndex(ma, mB);
	  cluster->Calculate();
	  if( cluster ) ClusterCont.push_back( cluster );
	}
      }
      else{
	hitA->SetJoined(ma);
	BH2Cluster *cluster = new BH2Cluster( hitA );
	cluster->SetIndex(ma);
	cluster->Calculate();
	if( cluster ) ClusterCont.push_back( cluster );
      }
    }// for(ma:hitA)
  }// for(i:hitA)

  return ClusterCont.size();
}

//______________________________________________________________________________
int
HodoAnalyzer::MakeUpClusters( const FiberHitContainer& cont,
			      FiberClusterContainer& ClusterCont,
			      double maxTimeDif,
			      int DifPairId )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  del::ClearContainer( ClusterCont );

  int NofSeg = cont.size();
  for( int seg=0; seg<NofSeg; ++seg ){
    FiberHit* HitA = cont.at(seg);
    bool fl_ClCandA = false;
    if( seg != (NofSeg -1) ){
      if( DifPairId > (cont.at(seg+1)->PairId() - HitA->PairId())){
	fl_ClCandA = true;
      }
    }

    int NofHitA = HitA->GetNumOfHit();
    for(int mhitA = 0; mhitA<NofHitA; ++mhitA){
      if(HitA->Joined(mhitA)) continue;

      FiberCluster *cluster = new FiberCluster();
      cluster->push_back( new FLHit(HitA, mhitA) );

      if(!fl_ClCandA){
	// there is no more candidates
	if(cluster->Calculate()){
	  ClusterCont.push_back(cluster);
	}else{
	  delete cluster;
	  cluster = NULL;
	}
	continue;
      }

      // Start Search HitB
      double cmtA    = (double)HitA->GetCTime(mhitA);
      int    NofHitB = cont.at(seg+1)->GetNumOfHit();
      bool   fl_HitB = false;
      double cmtB    = -1;
      int    CurrentPair = HitA->PairId();
      for(int mhitB = 0; mhitB<NofHitB; ++mhitB){
	if(cont.at(seg+1)->Joined(mhitB)){continue;}

	FiberHit* HitB = cont.at(seg+1);

	cmtB = (double)HitB->GetCTime(mhitB);
	if(std::abs(cmtB-cmtA)<maxTimeDif){
	  cluster->push_back(new FLHit(HitB, mhitB));
	  CurrentPair = HitB->PairId();
	  fl_HitB = true;
	  break;
	}
      }

      bool fl_ClCandB  = false;
      if((seg+1) != (NofSeg -1)){
	if( DifPairId > (cont.at(seg+2)->PairId() - CurrentPair) ){
	  fl_ClCandB = true;
	}
      }

      if(!fl_ClCandB){
	// there is no more candidates
	if(cluster->Calculate()){
	  ClusterCont.push_back(cluster);
	}else{
	  delete cluster;
	  cluster = NULL;
	}
	continue;
      }

      // Start Search HitC
      int    NofHitC = cont.at(seg+2)->GetNumOfHit();
      bool   fl_HitC = false;
      double cmtC    = -1;
      for(int mhitC = 0; mhitC<NofHitC; ++mhitC){
	if(cont.at(seg+2)->Joined(mhitC)){continue;}

	FiberHit* HitC = cont.at(seg+2);

	cmtC = (double)HitC->GetCTime(mhitC);
	if(true
	   && std::abs(cmtC-cmtA)<maxTimeDif
	   && !(fl_HitB && (std::abs(cmtC-cmtB)>maxTimeDif))
	   ){
	  cluster->push_back(new FLHit(HitC, mhitC));
	  CurrentPair = HitC->PairId();
	  fl_HitC = true;
	  break;
	}
      }

      bool fl_ClCandC  = false;
      if((seg+2) != (NofSeg -1)){
	if( DifPairId > (cont.at(seg+3)->PairId() - CurrentPair) ){
	  fl_ClCandC = true;
	}
      }

      if(!fl_ClCandC){
	// there is no more candidates
	if(cluster->Calculate()){
	  ClusterCont.push_back(cluster);
	}else{
	  delete cluster;
	  cluster = NULL;
	}
	continue;
      }

      // Start Search HitD
      int    NofHitD = cont.at(seg+3)->GetNumOfHit();
      double cmtD    = -1;
      for(int mhitD = 0; mhitD<NofHitD; ++mhitD){
	if(cont.at(seg+3)->Joined(mhitD)){continue;}

	FiberHit* HitD = cont.at(seg+3);

	cmtD = (double)HitD->GetCTime(mhitD);
	if(true
	   && std::abs(cmtD-cmtA)<maxTimeDif
	   && !(fl_HitB && (std::abs(cmtD-cmtB)>maxTimeDif))
	   && !(fl_HitC && (std::abs(cmtD-cmtC)>maxTimeDif))
	   ){
	  cluster->push_back(new FLHit(HitD, mhitD));
	  break;
	}
      }

      // Finish
      if(cluster->Calculate()){
	ClusterCont.push_back(cluster);
      }else{
	delete cluster;
	cluster = NULL;
      }
    }
  }

  return ClusterCont.size();
}

//______________________________________________________________________________
int
HodoAnalyzer::MakeUpClusters( const FLHitContainer& cont,
			      FiberClusterContainer& ClusterCont,
			      double maxTimeDif,
			      int DifPairId)
{
  // Clustering function for FBH

  std::vector<FiberCluster*> DeleteCand;

  int NofSeg = cont.size();
  for(int seg = 0; seg<NofSeg; ++seg){
    FLHit* HitA = cont.at(seg);

    if(HitA->Joined()){continue;}

    HitA->SetJoined();
    FiberCluster *cluster = new FiberCluster();
    cluster->push_back(HitA);

    int CurrentPair = HitA->PairId();
    for(int smarker = seg+1; smarker<NofSeg; ++smarker){
      FLHit* HitB = cont.at(smarker);

      if(false
	 || DifPairId < (HitB->PairId() - CurrentPair)
	 || HitB->PairId() == CurrentPair
	 ) {continue;}
      if(HitB->Joined()) continue;
      HitB->SetJoined();

      double cmtB = HitB->GetCTime();
      bool   fl_all_valid = true;
      for(int cindex = 0; cindex<cluster->VectorSize(); ++cindex){
	double cmt = cluster->GetHit(cindex)->GetCTime();
	if(std::abs(cmt-cmtB) > maxTimeDif){fl_all_valid = false; break;}
      }

      if(fl_all_valid){
	// we found a cluster candidate
	cluster->push_back(HitB);
	CurrentPair = HitB->PairId();
	break;
      }
    }

    // Finish
    if(cluster->Calculate()
       ){
      ClusterCont.push_back(cluster);
    }else{
      DeleteCand.push_back(cluster);

    }
  }

  del::ClearContainer( DeleteCand );

  return ClusterCont.size();
}

// MakeUpCoincidence ------------------------------------------------------
int
HodoAnalyzer::MakeUpCoincidence( const FiberHitContainer& cont,
				 FLHitContainer& CoinCont,
				 double maxTimeDif )
{
  int NofSeg = cont.size();
  for(int seg = 0; seg<NofSeg; ++seg){
    FiberHit* HitA = cont.at(seg);

    int NofHitA = HitA->GetNumOfHit();
    for( int mhitA = 0; mhitA<NofHitA; ++mhitA ){
      if( HitA->Joined(mhitA) ) continue;
      double cmt = HitA->GetCTime(mhitA);

      int CurrentPair = HitA->PairId();
      for( int smarker = seg+1; smarker<NofSeg; ++smarker ){
	FiberHit* HitB = cont.at(smarker);

	if( 0 != (HitB->PairId() - CurrentPair) ) continue;

	int NofHitB = HitB->GetNumOfHit();
	for( int mhitB = 0; mhitB<NofHitB; ++mhitB ){
	  if( HitB->Joined(mhitB) ) continue;

	  double cmtB = HitB->GetCTime(mhitB);
	  bool   fl_all_valid = true;
	  if( std::abs(cmt-cmtB)>maxTimeDif ){
	    fl_all_valid = false;
	    break;
	  }
	  if( fl_all_valid ){
	    // we found a coin candidate
	    //	    hddaq::cout << "yes" << std::endl;
	    CoinCont.push_back(new FLHit(HitA, HitB, mhitA, mhitB));

	    break;
	  }
	}// for(mhitB)
      }// for(segB)

    }// for(mhitA)
  }// for(segA)

  return CoinCont.size();
}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcBH1Hits( bool applyRecursively )
{
  int n = m_BH1Cont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit = m_BH1Cont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcBH2Hits( bool applyRecursively )
{
  int n = m_BH2Cont.size();
  for( int i=0; i<n; ++i ){
    BH2Hit *hit = m_BH2Cont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcBACHits( bool applyRecursively )
{
  int n = m_BACCont.size();
  for( int i=0; i<n; ++i ){
    Hodo1Hit *hit = m_BACCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcPVACHits( bool applyRecursively )
{
  int n = m_PVACCont.size();
  for( int i=0; i<n; ++i ){
    Hodo1Hit *hit = m_PVACCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcFACHits( bool applyRecursively )
{
  int n = m_FACCont.size();
  for( int i=0; i<n; ++i ){
    Hodo1Hit *hit = m_FACCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcSACHits( bool applyRecursively )
{
  int n = m_SACCont.size();
  for( int i=0; i<n; ++i ){
    Hodo1Hit *hit = m_SACCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcTOFHits( bool applyRecursively )
{
  int n = m_TOFCont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit = m_TOFCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;

}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcLCHits( bool applyRecursively )
{
  int n = m_LCCont.size();
  for( int i=0; i<n; ++i ){
    Hodo1Hit *hit = m_LCCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;

}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcBH1Clusters( bool applyRecursively )
{
  int n = m_BH1ClCont.size();
  for( int i=0; i<n; ++i ){
    HodoCluster *cl = m_BH1ClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcBH2Clusters( bool applyRecursively )
{
  int n = m_BH2ClCont.size();
  for( int i=0; i<n; ++i ){
    BH2Cluster *cl = m_BH2ClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcTOFClusters( bool applyRecursively )
{
  int n = m_TOFClCont.size();
  for( int i=0; i<n; ++i ){
    HodoCluster *cl = m_TOFClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcLCClusters( bool applyRecursively )
{
  int n = m_LCClCont.size();
  for( int i=0; i<n; ++i ){
    HodoCluster *cl = m_LCClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcAll( void )
{
  ReCalcBH1Hits();
  ReCalcBH2Hits();
  ReCalcBACHits();
  ReCalcPVACHits();
  ReCalcSACHits();
  ReCalcFACHits();
  ReCalcTOFHits();
  ReCalcLCHits();
  ReCalcBH1Clusters();
  ReCalcBH2Clusters();
  ReCalcTOFClusters();
  ReCalcLCClusters();
  return true;
}

//______________________________________________________________________________
void
HodoAnalyzer::TimeCutBH1( double tmin, double tmax )
{
  TimeCut( m_BH1ClCont, tmin, tmax );
}

//______________________________________________________________________________
void
HodoAnalyzer::TimeCutBH2( double tmin, double tmax )
{
  TimeCut( m_BH2ClCont, tmin, tmax );
}

//______________________________________________________________________________
void
HodoAnalyzer::TimeCutBFT( double tmin, double tmax )
{
  TimeCut( m_BFTClCont, tmin, tmax );
}

//______________________________________________________________________________
void
HodoAnalyzer::TimeCutSFT( int layer, double tmin, double tmax )
{
  TimeCut( m_SFTClCont.at( layer ), tmin, tmax );
}

//______________________________________________________________________________
void
HodoAnalyzer::TimeCutSCH( double tmin, double tmax )
{
  TimeCut( m_SCHClCont, tmin, tmax );
}

//______________________________________________________________________________
void
HodoAnalyzer::TimeCutFBH( double tmin, double tmax )
{
  TimeCut( m_FBHClCont, tmin, tmax );
}

//______________________________________________________________________________
//Implementation of Time cut for the cluster container
template <typename TypeCluster>
void
HodoAnalyzer::TimeCut( std::vector<TypeCluster>& cont,
		       double tmin, double tmax )
{
  std::vector<TypeCluster> DeleteCand;
  std::vector<TypeCluster> ValidCand;
  std::size_t size = cont.size();
  for( std::size_t i=0; i<size; ++i ){
    double ctime = cont.at(i)->CMeanTime();
    if(tmin < ctime && ctime < tmax){
      ValidCand.push_back(cont.at(i));
    }else{
      DeleteCand.push_back(cont.at(i));
    }
  }

  del::ClearContainer( DeleteCand );

  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}

//______________________________________________________________________________
void
HodoAnalyzer::WidthCutBFT( double min_width, double max_width)
{
  WidthCut( m_BFTClCont, min_width, max_width , true);
}

//______________________________________________________________________________
void
HodoAnalyzer::WidthCutSCH( double min_width, double max_width)
{
  WidthCut( m_SCHClCont, min_width, max_width , true);
}

//______________________________________________________________________________
void
HodoAnalyzer::WidthCutSFT( int layer, double min_width, double max_width)
{
  WidthCut( m_SFTClCont.at( layer ), min_width, max_width , true);
}

//______________________________________________________________________________
//Implementation of width cut for the cluster container
template <typename TypeCluster>
void
HodoAnalyzer::WidthCut( std::vector<TypeCluster>& cont,
			double min_width, double max_width,
			bool adopt_nan)
{
  std::vector<TypeCluster> DeleteCand;
  std::vector<TypeCluster> ValidCand;
  std::size_t size = cont.size();
  for( std::size_t i=0; i<size; ++i ){
    double width = cont.at(i)->Width();
    if(isnan(width) && adopt_nan){
      ValidCand.push_back(cont.at(i));
    }else if(min_width < width && width < max_width){
      ValidCand.push_back(cont.at(i));
    }else{
      DeleteCand.push_back(cont.at(i));
    }
  }

  del::ClearContainer( DeleteCand );

  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}
