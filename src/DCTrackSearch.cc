/**
 *  file: DCTrackSearch.cc
 *  date: 2017.04.10
 *
 */

#include "DCTrackSearch.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "DCLTrackHit.hh"
#include "DCPairHitCluster.hh"
#include "DCParameters.hh"
#include "DebugTimer.hh"
#include "DetectorID.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "MathTools.hh"
#include "MWPCCluster.hh"
#include "SsdCluster.hh"
#include "TrackMaker.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"

namespace
{
  const std::string& class_name("DCTrackSearch");
  DCGeomMan& gGeom = DCGeomMan::GetInstance();
  const double& zTarget    = gGeom.LocalZ("Target");
  const double& zK18Target = gGeom.LocalZ("K18Target");
  const double MaxChisquare       = 2000.; // Set to be More than 30
  const double MaxChisquareSdcIn  = 5000.; // Set to be More than 30
  const double MaxChisquareSdcOut = 50.;   // Set to be More than 30
  const double MaxNumOfCluster = 30.;	   // Set to be Less than 30
  const double MaxCombi = 1.0e6;	   // Set to be Less than 10^6
  // SdcIn & BcOut for XUV Tracking routine
  const double MaxChisquareVXU = 50.;//
  const double ChisquareCutVXU = 50.;//
  const double MaxChisquarePreSsdXY = 1.0e5;

  //_____________________________________________________________________
  // Local Functions ____________________________________________________

  //_____________________________________________________________________
  inline void
  CalcTracks( std::vector<DCLocalTrack*>& trackCont )
  {
    for( std::size_t i=0, n=trackCont.size(); i<n; ++i )
      trackCont[i]->Calculate();
  }

  //_____________________________________________________________________
  inline void
  ClearFlags( std::vector<DCLocalTrack*>& trackCont )
  {
    for( std::size_t i=0, n=trackCont.size(); i<n; ++i){
      const DCLocalTrack* const tp = trackCont[i];
      if (!tp) continue;
      int nh = tp->GetNHit();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->QuitTrack();
    }
  }

  //_____________________________________________________________________
  inline void
  SetKaonFlags( std::vector<DCLocalTrack*>& TrackCont )
  {
    for( std::size_t i=0, n=TrackCont.size(); i<n; ++i){
      const DCLocalTrack* const tp = TrackCont[i];
      if (!tp) continue;
      for( std::size_t ii=0, nn=tp->GetNHit(); ii<nn; ++ii ){
	tp->GetHit(ii)->JoinKaonTrack();
      }
    }
  }

  //_____________________________________________________________________
  inline void
  DeleteDuplicatedTracks( std::vector<DCLocalTrack*>& trackCont, double ChisqrCut=0. )
  {
    // evaluate container size in every iteration
    for( std::size_t i=0; i<trackCont.size(); ++i ){
      const DCLocalTrack* const tp = trackCont[i];
      if (!tp) continue;
      int nh = tp->GetNHit();
      for(int j=0; j<nh; ++j) tp->GetHit(j)->JoinTrack();

      for( std::size_t i2=trackCont.size()-1; i2>i; --i2 ){
	const DCLocalTrack* tp2 = trackCont[i2];
	int nh2 = tp2->GetNHit(), flag=0;
	double chisqr = tp2->GetChiSquare();
	for( int j=0; j<nh2; ++j )
	  if( tp2->GetHit(j)->BelongToTrack() ) ++flag;
	if( flag>0 && chisqr>ChisqrCut ){
	  delete tp2;
	  tp2 = 0;
	  trackCont.erase(trackCont.begin()+i2);
	}
      }
    }
  }

  //_____________________________________________________________________
  inline void // for SSD PreTrack
  DeleteWideTracks( std::vector<DCLocalTrack*>& TrackContX,
		    std::vector<DCLocalTrack*>& TrackContY,
		    const std::size_t Nth=1 )
  {
    const double xSize = 400./2*1.; // SDC1 X-Size
    const double ySize = 250./2*1.; // SDC1 Y-Size
    const double zPos  = -658.;  // SDC1 Z-Position "V1"
    // X
    if( TrackContX.size()>Nth ){
      for( int i=TrackContX.size()-1; i>=0; --i ){
	const DCLocalTrack* tp = TrackContX[i];
	if (!tp) continue;
	double x = tp->GetVXU(zPos);
	if( std::abs(x)>xSize ){
	  delete tp;
	  tp = 0;
	  TrackContX.erase( TrackContX.begin()+i );
	}
      }
    }
    // Y
    if( TrackContY.size()>Nth ){
      for( int i=TrackContY.size()-1; i>=0; --i ){
	const DCLocalTrack* tp = TrackContY[i];
	if (!tp) continue;
	double y = tp->GetVXU(zPos);
	if( std::abs(y)>ySize ){
	  delete tp;
	  tp = 0;
	  TrackContY.erase( TrackContY.begin()+i );
	}
      }
    }
  }

  //_____________________________________________________________________
  inline void // for SSD PreTrack
  DeleteWideTracks( std::vector<DCLocalTrack*>& TrackContX,
		    const std::vector<DCHitContainer>& SdcInHC,
		    const std::size_t Nth=1 )
  {
    const double xSize =  100.;    // MaxDiff
    const double zPos  = -635.467; // SDC1 Z-Position "X1/X2"
    // X
    if( TrackContX.size()>Nth ){
      for( int i=TrackContX.size()-1; i>=0; --i ){
	const DCLocalTrack* tp = TrackContX[i];
	if (!tp) continue;
	double x = tp->GetVXU(zPos);

	bool accept = false;
	for( std::size_t j=0; j<NumOfLayersSDC1; ++j ){
	  if( j!=2 && j!=3 ) continue; // only X1/X2
	  const std::size_t nh = SdcInHC[j+1].size();
	  for( std::size_t k=0; k<nh; ++k ){
	    const DCHit *hit = SdcInHC[j+1][k];
	    if( !hit ) continue;
	    if( std::abs( x-hit->GetWirePosition() ) < xSize ){
	      accept = true;
	    }
	  }
	}
	if( !accept ){
	  delete tp;
	  tp = 0;
	  TrackContX.erase( TrackContX.begin()+i );
	}
      }
    }

    return;
  }

  //_____________________________________________________________________
  void
  DebugPrint( const IndexList& nCombi,
	      const std::string& func_name="",
	      const std::string& msg="" )
  {
    int n  =1;
    int nn =1;
    int sum=0;
    hddaq::cout << "#D " << func_name << ":" ;
    IndexList::const_iterator itr, end = nCombi.end();
    for( itr=nCombi.begin(); itr!=end; ++itr ){
      int val = *itr;
      sum += val;
      nn *= (val+1);
      hddaq::cout << " " << val;
      if ( val!=0 ) n *= val;
    }
    if (sum==0)
      n=0;
    hddaq::cout << ": total = " << n << ", " << nn << ", " << std::endl;
    return;
  }

  //_____________________________________________________________________
  void
  DebugPrint( const std::vector<DCLocalTrack*>& trackCont,
	      const std::string& arg="" )
  {
    const int nn = trackCont.size();
    hddaq::cout << arg << " " << nn << std::endl;
    for( int i=0; i<nn; ++i ){
      const DCLocalTrack * const track=trackCont[i];
      if (!track) continue;
      int    nh     = track->GetNHit();
      double chisqr = track->GetChiSquare();
      hddaq::cout << std::setw(4) << i
		  << "  #Hits : " << std::setw(2) << nh
		  << "  ChiSqr : " << chisqr << std::endl;
    }
    hddaq::cout << std::endl;
  }

  //_____________________________________________________________________
  void
  DebugPrint( const IndexList& nCombi,
	      const std::vector<ClusterList>& CandCont,
	      const std::string& arg="" )
  {
    hddaq::cout << arg << " #Hits of each group" << std::endl;
    int np = nCombi.size();
    int nn = 1;
    for( int i=0; i<np; ++i ){
      hddaq::cout << std::setw(4) << nCombi[i];
      nn *= nCombi[i] + 1;
    }
    hddaq::cout << " -> " << nn-1 << " Combinations" << std::endl;
    for( int i=0; i<np; ++i ){
      int n=CandCont[i].size();
      hddaq::cout << "[" << std::setw(3) << i << "]: "
		  << std::setw(3) << n << " ";
      for( int j=0; j<n; ++j ){
	hddaq::cout << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire() << " ";
      }
      hddaq::cout << std::endl;
    }
  }

  //______________________________________________________________________________
  template <class Functor>
  inline void
  FinalizeTrack( const std::string& arg,
		 std::vector<DCLocalTrack*>& trackCont,
		 Functor comp,
		 std::vector<ClusterList>& candCont,
		 bool delete_flag=true )
  {
    ClearFlags(trackCont);

#if 0
    DebugPrint( trackCont, arg+" Before Sorting " );
#endif

    std::stable_sort( trackCont.begin(), trackCont.end(), comp );

#if 0
    DebugPrint( trackCont, arg+" After Sorting " );
#endif

    if( delete_flag )
      DeleteDuplicatedTracks( trackCont );

#if 0
    DebugPrint( trackCont, arg+" After Deleting " );
#endif
    
    CalcTracks( trackCont );
    del::ClearContainerAll( candCont );
  }

  //______________________________________________________________________________
  // MakeCluster _________________________________________________________________

  //______________________________________________________________________________
  bool
  MakePairPlaneHitCluster( const DCHitContainer & HC1,
			   const DCHitContainer & HC2,
			   double CellSize,
			   ClusterList& Cont,
			   bool honeycomb=false )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    int nh1=HC1.size(), nh2=HC2.size();
    std::vector<int> UsedFlag(nh2,0);
    for( int i1=0; i1<nh1; ++i1 ){
      DCHit *hit1=HC1[i1];
      double wp1=hit1->GetWirePosition();
      bool flag=false;
      for( int i2=0; i2<nh2; ++i2 ){
	DCHit *hit2=HC2[i2];
	double wp2=hit2->GetWirePosition();
	if( std::abs(wp1-wp2)<CellSize ){
	  int multi1 = hit1->GetDriftLengthSize();
	  int multi2 = hit2->GetDriftLengthSize();
	  for ( int m1=0; m1<multi1; ++m1 ) {
	    if( !hit1->IsWithinRange(m1) )
	      continue;
	    for ( int m2=0; m2<multi2; ++m2 ) {
	      if( !hit2->IsWithinRange(m2) )
		continue;
	      double x1,x2;
	      if( wp1<wp2 ){
		x1=wp1+hit1->GetDriftLength(m1);
		x2=wp2-hit2->GetDriftLength(m2);
	      }
	      else {
		x1=wp1-hit1->GetDriftLength(m1);
		x2=wp2+hit2->GetDriftLength(m2);
	      }
	      DCPairHitCluster *cluster =
		new DCPairHitCluster( new DCLTrackHit(hit1,x1,m1),
				      new DCLTrackHit(hit2,x2,m2) );
	      cluster->SetHoneycomb( honeycomb );
	      Cont.push_back( cluster );
	      flag=true; ++UsedFlag[i2];
	    }
	  }
	}
      }
#if 1
      if(!flag){
	int multi1 = hit1->GetDriftLengthSize();
	for (int m1=0; m1<multi1; m1++) {
	  if( !(hit1->IsWithinRange(m1)) ) continue;
	  double dl=hit1->GetDriftLength(m1);
	  DCPairHitCluster *cluster1 = new DCPairHitCluster( new DCLTrackHit(hit1,wp1+dl,m1) );
	  DCPairHitCluster *cluster2 = new DCPairHitCluster( new DCLTrackHit(hit1,wp1-dl,m1) );
	  cluster1->SetHoneycomb( honeycomb );
	  cluster2->SetHoneycomb( honeycomb );
	  Cont.push_back( cluster1 );
	  Cont.push_back( cluster2 );
	}
      }
#endif
    }
#if 1
    for( int i2=0; i2<nh2; ++i2 ){
      if( UsedFlag[i2]==0 ) {
	DCHit *hit2=HC2[i2];
	int multi2 = hit2->GetDriftLengthSize();
	for (int m2=0; m2<multi2; m2++) {
	  if( !(hit2->IsWithinRange(m2)) ) continue;
	  double wp=hit2->GetWirePosition();
	  double dl=hit2->GetDriftLength(m2);
	  DCPairHitCluster *cluster1 = new DCPairHitCluster( new DCLTrackHit(hit2,wp+dl,m2) );
	  DCPairHitCluster *cluster2 = new DCPairHitCluster( new DCLTrackHit(hit2,wp-dl,m2) );
	  cluster1->SetHoneycomb( honeycomb );
	  cluster2->SetHoneycomb( honeycomb );
	  Cont.push_back( cluster1 );
	  Cont.push_back( cluster2 );
	}
      }
    }
#endif
    return true;
  }

  //______________________________________________________________________________
  bool
  MakeUnPairPlaneHitCluster( const DCHitContainer& HC,
			     ClusterList& Cont,
			     bool honeycomb=false )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    const std::size_t nh = HC.size();
    for( std::size_t i=0; i<nh; ++i ){
      DCHit *hit = HC[i];
      if( !hit ) continue;
      std::size_t mh = hit->GetDriftLengthSize();
      for ( std::size_t m=0; m<mh; ++m ) {
	if( !hit->IsWithinRange(m) ) continue;
	double wp = hit->GetWirePosition();
	double dl = hit->GetDriftLength(m);
	DCPairHitCluster *cluster1 =
	  new DCPairHitCluster( new DCLTrackHit(hit,wp+dl,m) );
	DCPairHitCluster *cluster2 =
	  new DCPairHitCluster( new DCLTrackHit(hit,wp-dl,m) );
	cluster1->SetHoneycomb( honeycomb );
	cluster2->SetHoneycomb( honeycomb );
	Cont.push_back( cluster1 );
	Cont.push_back( cluster2 );
      }
    }

    return true;
  }

  //______________________________________________________________________________
  bool
  MakeMWPCPairPlaneHitCluster( const DCHitContainer& HC,
			       ClusterList& Cont )
  {
    int nh=HC.size();
    for( int i=0; i<nh; ++i ){
      DCHit *hit=HC[i];
      if( hit ){
	int multi = hit->GetDriftTimeSize();
	for (int m=0; m<multi; m++) {
	  if( !(hit->IsWithinRange(m)) ) continue;
	  double wp=hit->GetWirePosition();
	  // double dl=hit->GetDriftLength(m);
	  Cont.push_back( new DCPairHitCluster( new DCLTrackHit( hit, wp, m ) ) );
	}
      }
    }
    return true;
  }

  //______________________________________________________________________________
  bool
  MakeSsdHitCluster( const DCHitContainer& HC,
		     ClusterList& Cont )
  {
    int nh = HC.size();
    for( int i=0; i<nh; ++i ){
      DCHit *hit = HC[i];
      if( !hit ) continue;
      if( !hit->IsGoodWaveForm() ) continue;
      if( hit->BelongToKaonTrack() ) continue;
      double wp = hit->GetWirePosition();
      Cont.push_back( new DCPairHitCluster( new DCLTrackHit( hit, wp, 0 ) ) );
    }
    return true;
  }

  //______________________________________________________________________________
  bool
  MakeSsdHitCluster( const SsdClusterContainer& ClCont,
		     ClusterList& Cont )
  {
    int ncl = ClCont.size();
    for( int i=0; i<ncl; ++i ){
      SsdCluster *cluster = ClCont[i];
      if( !cluster ) continue;
      if( !cluster->GoodForAnalysis() )  continue;
      if( cluster->BelongToKaonTrack() ) continue;
      DCHit  *hit = cluster->GetHit();
      double  wp  = hit->GetWirePosition();
      Cont.push_back( new DCPairHitCluster( new DCLTrackHit( hit, wp, 0 ) ) );
    }
    return true;
  }

  //______________________________________________________________________________
  bool
  MakeTOFHitCluster( const DCHitContainer& HitCont,
		     ClusterList& Cont,
		     int xy )
  {
    int nh = HitCont.size();
    for( int i=0; i<nh; ++i ){
      if( i%2!=xy ) continue;
      DCHit *hit = HitCont[i];
      if( !hit ) continue;
      double wp = hit->GetWirePosition();
      Cont.push_back( new DCPairHitCluster( new DCLTrackHit( hit, wp, 0 ) ) );
    }
    return true;
  }

  //______________________________________________________________________________
  bool
  MakePairPlaneHitClusterVUX( const DCHitContainer& HC1,
			      const DCHitContainer& HC2,
			      double CellSize,
			      ClusterList& Cont,
			      bool honeycomb=false )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    int nh1=HC1.size(), nh2=HC2.size();
    std::vector<int> UsedFlag(nh2,0);

    for( int i1=0; i1<nh1; ++i1 ){
      DCHit *hit1=HC1[i1];

      double wp1=hit1->GetWirePosition();
      bool flag=false;
      for( int i2=0; i2<nh2; ++i2 ){
	DCHit *hit2=HC2[i2];
	double wp2=hit2->GetWirePosition();
	if( std::abs(wp1-wp2)<CellSize ){

	  int multi1 = hit1->GetDriftLengthSize();
	  int multi2 = hit2->GetDriftLengthSize();
	  for (int m1=0; m1<multi1; m1++) {
	    if( !(hit1->IsWithinRange(m1)) ) continue;
	    for (int m2=0; m2<multi2; m2++) {
	      if( !(hit2->IsWithinRange(m2)) ) continue;
	      double dl1=hit1->GetDriftLength(m1);
	      double dl2=hit2->GetDriftLength(m2);

	      Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1+dl1,m1),
						    new DCLTrackHit(hit2,wp2+dl2,m2) ) );
	      Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1+dl1,m1),
						    new DCLTrackHit(hit2,wp2-dl2,m2) ) );
	      Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1-dl1,m1),
						    new DCLTrackHit(hit2,wp2+dl2,m2) ) );
	      Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1-dl1,m1),
						    new DCLTrackHit(hit2,wp2-dl2,m2) ) );

	      flag=true; ++UsedFlag[i2];
	    }
	  }
	}
      }
      if(!flag){
	int multi1 = hit1->GetDriftLengthSize();
	for (int m1=0; m1<multi1; m1++) {
	  if( !(hit1->IsWithinRange(m1)) ) continue;
	  double dl=hit1->GetDriftLength(m1);
	  Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1+dl,m1) ) );
	  Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1-dl,m1) ) );
	}
      }
    }
    for( int i2=0; i2<nh2; ++i2 ){
      if( UsedFlag[i2]==0 ) {
	DCHit *hit2=HC2[i2];
	int multi2 = hit2->GetDriftLengthSize();
	for (int m2=0; m2<multi2; m2++) {
	  if( !(hit2->IsWithinRange(m2)) ) continue;

	  double wp=hit2->GetWirePosition();
	  double dl=hit2->GetDriftLength(m2);
	  Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit2,wp+dl,m2) ) );
	  Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit2,wp-dl,m2) ) );
	}
      }
    }

    return true;
  }

  //______________________________________________________________________________
  // MakeIndex ___________________________________________________________________

  //______________________________________________________________________________
  std::vector<IndexList>
  MakeIndex( int ndim, const int *index1 )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    if( ndim==1 ){
      std::vector<IndexList> index2;
      for( int i=-1; i<index1[0]; ++i ){
	IndexList elem(1,i);
	index2.push_back(elem);
      }
      return index2;
    }

    std::vector<IndexList> index2 = MakeIndex( ndim-1, index1+1 );

    std::vector<IndexList> index;
    int n2=index2.size();
    for( int j=0; j<n2; ++j ){
      for( int i=-1; i<index1[0]; ++i ){
	IndexList elem;
	int n3=index2[j].size();
	elem.reserve(n3+1);
	elem.push_back(i);
	for( int k=0; k<n3; ++k )
	  elem.push_back(index2[j][k]);
	index.push_back(elem);
	int size1=index.size();
	if( size1>MaxCombi ){
#if 0
	  hddaq::cout << func_name << " too much combinations..." << std::endl;
#endif
	  return std::vector<IndexList>(0);
	}
      }
    }

    return index;
  }

  //______________________________________________________________________________
  std::vector<IndexList>
  MakeIndex( int ndim, const IndexList& index1 )
  {
    return MakeIndex( ndim, &(index1[0]) );
  }

  //______________________________________________________________________________
  std::vector<IndexList>
  MakeIndex_VXU( int ndim,int maximumHit, const int *index1 )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    if( ndim==1 ){
      std::vector<IndexList> index2;
      for( int i=-1; i<index1[0]; ++i ){
	IndexList elem(1,i);
	index2.push_back(elem);
      }
      return index2;
    }

    std::vector<IndexList>
      index2=MakeIndex_VXU( ndim-1, maximumHit, index1+1 );

    std::vector<IndexList> index;
    int n2=index2.size();
    for( int j=0; j<n2; ++j ){
      for( int i=-1; i<index1[0]; ++i ){
	IndexList elem;
	int validHitNum=0;
	int n3=index2[j].size();
	elem.reserve(n3+1);
	elem.push_back(i);
	if (i != -1)
	  validHitNum++;
	for( int k=0; k<n3; ++k ){
	  elem.push_back(index2[j][k]);
	  if (index2[j][k] != -1)
	    validHitNum++;
	}
	if (validHitNum <= maximumHit)
	  index.push_back(elem);
	// int size1=index.size();
      }
    }

    return index;
  }

  //______________________________________________________________________________
  std::vector<IndexList>
  MakeIndex_VXU( int ndim,int maximumHit, const IndexList& index1 )
  {
    return MakeIndex_VXU( ndim, maximumHit, &(index1[0]) );
  }

  //______________________________________________________________________________
  // MakeTrack ___________________________________________________________________

  //______________________________________________________________________________
  DCLocalTrack*
  MakeTrack( const std::vector<ClusterList>& CandCont,
	     const IndexList& combination )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    DCLocalTrack *tp = new DCLocalTrack;
    for( std::size_t i=0, n=CandCont.size(); i<n; ++i ){
      int m = combination[i];
      if( m<0 ) continue;
      DCPairHitCluster *cluster = CandCont[i][m];
      if( !cluster ) continue;
      int mm = cluster->NumberOfHits();
      for( int j=0; j<mm; ++j ){
	DCLTrackHit *hitp = cluster->GetHit(j);
	if( !hitp ) continue;
	tp->AddHit( hitp );
      }
#if 0
      hddaq::cout << func_name << ":" << std::setw(3)
		  << i << std::setw(3) << m  << " "
		  << CandCont[i][m] << " " << mm << std::endl;
#endif
    }
    return tp;
  }

}

//______________________________________________________________________________
namespace track
{
  //______________________________________________________________________________
  int /* Local Track Search without BH2Filter */
  LocalTrackSearch( const std::vector<DCHitContainer>& HC,
		    const DCPairPlaneInfo * PpInfo,
		    int npp, std::vector<DCLocalTrack*>& TrackCont,
		    int MinNumOfHits /*=6*/ )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    std::vector<ClusterList> CandCont(npp);

    for( int i=0; i<npp; ++i ){
      bool ppFlag    = PpInfo[i].pair;
      bool honeycomb = PpInfo[i].honeycomb;
      bool fiber     = PpInfo[i].fiber;
      int  layer1    = PpInfo[i].id1;
      int  layer2    = PpInfo[i].id2;
      if( ppFlag && !fiber ){
	MakePairPlaneHitCluster( HC[layer1], HC[layer2],
				 PpInfo[i].CellSize, CandCont[i], honeycomb );
      }else if( !ppFlag && fiber ){
	PpInfo[i].Print( func_name+" invalid parameter", hddaq::cerr );
      }else{
	MakeUnPairPlaneHitCluster( HC[layer1], CandCont[i], honeycomb );
      }
    }

    IndexList nCombi(npp);
    for( int i=0; i<npp; ++i ){
      int n = CandCont[i].size();
      nCombi[i] = n>MaxNumOfCluster ? 0 : n;
    }

#if 0
    DebugPrint( nCombi, CandCont, func_name );
#endif

    std::vector<IndexList> CombiIndex = MakeIndex( npp, nCombi );

    for( int i=0, n=CombiIndex.size(); i<n; ++i ){
      DCLocalTrack *track = MakeTrack( CandCont, CombiIndex[i] );
      if( !track ) continue;
      if( track->GetNHit()>=MinNumOfHits
    	  && track->DoFit()
    	  && track->GetChiSquare()<MaxChisquare ){
    	TrackCont.push_back( track );
      }
      else
    	delete track;
    }

    FinalizeTrack( func_name, TrackCont, DCLTrackComp(), CandCont );
    return TrackCont.size();
  }

  //______________________________________________________________________________
  int /* Local Track Search with BH2Filter */
  LocalTrackSearch( const std::vector< std::vector<DCHitContainer> > &hcAssemble,
		    const DCPairPlaneInfo * PpInfo,
		    int npp, std::vector<DCLocalTrack*> &trackCont,
		    int MinNumOfHits /*=6*/ )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"(BH2Filter)]");

    std::vector< std::vector<DCHitContainer> >::const_iterator
      itr, itr_end = hcAssemble.end();
    for ( itr=hcAssemble.begin(); itr!=itr_end; ++itr ){
      const std::vector<DCHitContainer>& l = *itr;
      std::vector<DCLocalTrack*> tc;
      LocalTrackSearch( l, PpInfo, npp, tc, MinNumOfHits );
      trackCont.insert( trackCont.end(), tc.begin(), tc.end() );
    }

    ClearFlags(trackCont);
    std::stable_sort( trackCont.begin(), trackCont.end(), DCLTrackComp() );
    DeleteDuplicatedTracks( trackCont );
    CalcTracks( trackCont );

    return trackCont.size();
  }

  //______________________________________________________________________________
  int
  LocalTrackSearchSsdOutSdcIn( const std::vector<DCHitContainer>& SdcInHC,
			       const DCPairPlaneInfo *PpInfo,
			       int SdcInNpp,
			       std::vector<DCLocalTrack*>& SsdXTC,
			       std::vector<DCLocalTrack*>& SsdYTC,
			       std::vector<DCLocalTrack*>& TrackCont,
			       int MinNumOfHits /*=10*/ )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"(SsdPreTrack)]");

    DeleteWideTracks( SsdXTC, SsdYTC );
    DeleteWideTracks( SsdXTC, SdcInHC );

    const int nSsdXTC = SsdXTC.size();
    const int nSsdYTC = SsdYTC.size();

    std::vector< std::vector<DCPairHitCluster*> > CandCont(SdcInNpp);
    for( int i=0; i<SdcInNpp; ++i ){
      bool ppFlag    = PpInfo[i].pair;
      bool honeycomb = PpInfo[i].honeycomb;
      bool fiber     = PpInfo[i].fiber;
      int  layer1    = PpInfo[i].id1;
      int  layer2    = PpInfo[i].id2;
      if( ppFlag && !fiber ){
	MakePairPlaneHitCluster( SdcInHC[layer1], SdcInHC[layer2],
				 PpInfo[i].CellSize, CandCont[i], honeycomb );
      }else if( !ppFlag && fiber ){
	PpInfo[i].Print( func_name+" invalid paramter", hddaq::cerr );
      }else{
	MakeUnPairPlaneHitCluster( SdcInHC[layer1], CandCont[i], honeycomb );
      }
    }

    IndexList nCombi(SdcInNpp);
    for( int i=0; i<SdcInNpp; ++i ){
      int n = CandCont[i].size();
      nCombi[i] = n>MaxNumOfCluster ? 0 : n;
    }

#if 0
    DebugPrint( nCombi, CandCont, func_name );
#endif

    std::vector<IndexList> CombiIndex = MakeIndex( SdcInNpp, nCombi );
    int nnCombi = CombiIndex.size();
    if( nnCombi==0 ){
      del::ClearContainerAll( CandCont );
      return -1;
    }

    if( nnCombi > 1.0e5 ){
      del::ClearContainerAll( CandCont );
      // hddaq::cout << func_name << " too much SdcIn combinations... "
      // 	      << nnCombi << std::endl;
      return -1;
    }

    int nnn = nnCombi * nSsdXTC * nSsdYTC;
    if( nnn > 1.0e6 || nnn < 0 ){
      del::ClearContainerAll( CandCont );
      // hddaq::cout << func_name << " too much combinations... "
      // 	      << nnn << std::endl;
      return -1;
    }

    for( int iX=0; iX<nSsdXTC; ++iX ){
      const DCLocalTrack * const tpX = SsdXTC[iX];
      if( !tpX ) continue;
      const int nhX = tpX->GetNHit();
      for( int iY=0; iY<nSsdYTC; ++iY ){
	const DCLocalTrack * const tpY = SsdYTC[iY];
	if( !tpY ) continue;
	const int nhY = tpY->GetNHit();

#if 0
	hddaq::cout << " iX : " << std::setw(3) << iX+1 << "/" << std::setw(3) << nSsdXTC << " "
		    << " iY : " << std::setw(3) << iY+1 << "/" << std::setw(3) << nSsdYTC << " "
		    << " nhX : " << nhX << " " << " nhY : " << nhY
		    << " nSdcIn : " << nnCombi << " " << " NNN : " << nnn << std::endl << "\x1b[1A";
#endif

	for( int i=0; i<nnCombi; ++i ){
	  DCLocalTrack *track = MakeTrack( CandCont, CombiIndex[i] );
	  if( !track ) continue;

	  for( int j=0; j<nhX; ++j )
	    track->AddHit( tpX->GetHit(j) );
	  for( int j=0; j<nhY; ++j )
	    track->AddHit( tpY->GetHit(j) );

	  if( track->GetNHit()>=MinNumOfHits
	      && track->DoFit()
	      && track->GetChiSquare()<MaxChisquare
	      ){
	    TrackCont.push_back( track );
	  }
	  else
	    delete track;
	}
      }// for(iY)
    }// for(iX)

    FinalizeTrack( func_name, TrackCont, DCLTrackComp(), CandCont );
    return TrackCont.size();
  }

  //______________________________________________________________________________
  int
  LocalTrackSearchSdcOut( const std::vector<DCHitContainer>& SdcOutHC,
			  const DCPairPlaneInfo *PpInfo,
			  int npp, std::vector<DCLocalTrack*>& TrackCont,
			  int MinNumOfHits /*=6*/ )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    std::vector<ClusterList> CandCont(npp);

    for( int i=0; i<npp; ++i ){
      bool ppFlag    = PpInfo[i].pair;
      bool honeycomb = PpInfo[i].honeycomb;
      int  layer1    = PpInfo[i].id1;
      int  layer2    = PpInfo[i].id2;
      if( ppFlag  )
	MakePairPlaneHitCluster( SdcOutHC[layer1], SdcOutHC[layer2],
				 PpInfo[i].CellSize, CandCont[i], honeycomb );
      else
	MakeUnPairPlaneHitCluster( SdcOutHC[layer1], CandCont[i], honeycomb );
    }

    IndexList nCombi(npp);
    for( int i=0; i<npp; ++i ){
      int n = CandCont[i].size();
      nCombi[i] = n>MaxNumOfCluster ? 0 : n;
    }

    std::vector<IndexList> CombiIndex = MakeIndex( npp, nCombi );

#if 0
    DebugPrint( nCombi, CandCont, func_name );
#endif

    for ( int i=0, n=CombiIndex.size(); i<n; ++i ) {
      DCLocalTrack *track = MakeTrack( CandCont, CombiIndex[i] );
      if ( !track ) continue;
      if ( track->GetNHit()>=MinNumOfHits     &&
	   track->DoFit()                     &&
	   track->GetChiSquare()<MaxChisquare )
	{
	  TrackCont.push_back( track );
	}
      else
	delete track;
    }

    FinalizeTrack( func_name, TrackCont, DCLTrackComp(), CandCont );
    return TrackCont.size();
  }

  //______________________________________________________________________________
  int
  LocalTrackSearchSdcOut( const DCHitContainer& TOFHC,
			  const std::vector<DCHitContainer>& SdcOutHC,
			  const DCPairPlaneInfo *PpInfo,
			  int npp,
			  std::vector<DCLocalTrack*>& TrackCont,
			  int MinNumOfHits /*=6*/ )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"(TOF)]");

    std::vector<ClusterList> CandCont(npp);

    for( int i=0; i<npp-2; ++i ){
      bool ppFlag    = PpInfo[i].pair;
      bool honeycomb = PpInfo[i].honeycomb;
      int  layer1    = PpInfo[i].id1;
      int  layer2    = PpInfo[i].id2;
      if( ppFlag  )
	MakePairPlaneHitCluster( SdcOutHC[layer1], SdcOutHC[layer2],
				 PpInfo[i].CellSize, CandCont[i], honeycomb );
      else
	MakeUnPairPlaneHitCluster( SdcOutHC[layer1], CandCont[i], honeycomb );
    }

    MakeTOFHitCluster( TOFHC, CandCont[npp-2], 0 );
    MakeTOFHitCluster( TOFHC, CandCont[npp-1], 1 );

    IndexList nCombi(npp);
    for( int i=0; i<npp; ++i ){
      int n = CandCont[i].size();
      nCombi[i] = n>MaxNumOfCluster ? 0 : n;
    }

    std::vector<IndexList> CombiIndex = MakeIndex( npp, nCombi );

#if 0
    DebugPrint( nCombi, CandCont, func_name );
#endif

    for( int i=0, n=CombiIndex.size(); i<n; ++i ){
      DCLocalTrack *track = MakeTrack( CandCont, CombiIndex[i] );
      if( !track ) continue;
      static const int IdTOF_UX = gGeom.GetDetectorId("TOF-UX");
      static const int IdTOF_UY = gGeom.GetDetectorId("TOF-UY");
      static const int IdTOF_DX = gGeom.GetDetectorId("TOF-DX");
      static const int IdTOF_DY = gGeom.GetDetectorId("TOF-DY");
      bool TOFSegXYMatching =
	( track->GetWire(IdTOF_UX)==track->GetWire(IdTOF_UY) ) ||
	( track->GetWire(IdTOF_DX)==track->GetWire(IdTOF_DY) );

      if( TOFSegXYMatching &&
	  track->GetNHit()>=MinNumOfHits+2   &&
	  track->DoFit()                     &&
	  track->GetChiSquare()<MaxChisquare ){
	TrackCont.push_back( track );
      }
      else
	delete track;
    }

    FinalizeTrack( func_name, TrackCont, DCLTrackComp(), CandCont );
    return TrackCont.size();
  }

  //______________________________________________________________________________
  int /* Local Track Search SdcIn w/Fiber */
  LocalTrackSearchSdcInFiber( const std::vector<DCHitContainer>& HC,
			      const DCPairPlaneInfo *PpInfo,
			      int npp, std::vector<DCLocalTrack*>& TrackCont,
			      int MinNumOfHits /*=6*/ )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    std::vector<ClusterList> CandCont(npp);

    for( int i=0; i<npp; ++i ){

      bool ppFlag    = PpInfo[i].pair;
      bool honeycomb = PpInfo[i].honeycomb;
      int  layer1    = PpInfo[i].id1;
      int  layer2    = PpInfo[i].id2;
      
      if(ppFlag) {
	MakePairPlaneHitCluster( HC[layer1], HC[layer2],
				 PpInfo[i].CellSize, CandCont[i], honeycomb );
      }else if(layer1==layer2){
	MakeMWPCPairPlaneHitCluster( HC[layer1], CandCont[i] );
      }else{
	MakeUnPairPlaneHitCluster( HC[layer1], CandCont[i], honeycomb );
      }
    }

    IndexList nCombi(npp);
    for ( int i=0; i<npp; ++i ) {
      int n = CandCont[i].size();
      nCombi[i] = n>MaxNumOfCluster ? 0 : n;
    }

    std::vector<IndexList> CombiIndex = MakeIndex( npp, nCombi );

    for( int i=0, n=CombiIndex.size(); i<n; ++i ){
      DCLocalTrack *track = MakeTrack( CandCont, CombiIndex[i] );
      if( !track ) continue;
      if( track->GetNHit()>=MinNumOfHits && track->DoFit() &&
	  track->GetChiSquare()<MaxChisquare ){
	TrackCont.push_back(track);
      }
      else
	delete track;
    }

    FinalizeTrack( func_name, TrackCont, DCLTrackCompSdcInFiber(), CandCont );
    return TrackCont.size();
  }

  // BC3&4, SDC1 VUX Tracking ___________________________________________
  int
  LocalTrackSearchVUX( const std::vector<DCHitContainer>& HC,
		       const DCPairPlaneInfo *PpInfo,
		       int npp, std::vector<DCLocalTrack*>& TrackCont,
		       int MinNumOfHits /*=6*/ )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    std::vector<DCLocalTrack*> TrackContV;
    std::vector<DCLocalTrack*> TrackContX;
    std::vector<DCLocalTrack*> TrackContU;
    std::vector<ClusterList>   CandCont(npp);
    std::vector<ClusterList>   CandContV(npp);
    std::vector<ClusterList>   CandContX(npp);
    std::vector<ClusterList>   CandContU(npp);

    //  int NumOfLayersDC_12 = 12 ;
    int iV=0, iX=0, iU=0;
    int nV=0, nX=0, nU=0;
    for( int i=0; i<npp; ++i ){
      bool ppFlag    = PpInfo[i].pair;
      bool honeycomb = PpInfo[i].honeycomb;
      int  layer1    = PpInfo[i].id1;
      int  layer2    = PpInfo[i].id2;
      int  nh1       = HC[layer1].size();
      int  nh2       = HC[layer2].size();
      double TiltAngle = 0.;
      if(nh1>0) TiltAngle = HC[layer1][0]->GetTiltAngle();
      if( ppFlag && nh1==0 && nh2>0 )
	TiltAngle = HC[layer2][0]->GetTiltAngle();
      if( ppFlag && nh1>0 && nh2>0 ){
	if( TiltAngle<0 ){
	  MakePairPlaneHitClusterVUX( HC[layer1], HC[layer2],
				      PpInfo[i].CellSize, CandContV[iV],
				      honeycomb );
	  ++iV; nV = nV+2;
	}
	if( TiltAngle==0 ){
	  MakePairPlaneHitClusterVUX( HC[layer1], HC[layer2],
				      PpInfo[i].CellSize, CandContX[iX],
				      honeycomb );
	  ++iX; nX = nX+2;
	}
	if( TiltAngle>0 ){
	  MakePairPlaneHitClusterVUX( HC[layer1], HC[layer2],
				      PpInfo[i].CellSize, CandContU[iU],
				      honeycomb );
	  ++iU; nU = nU+2;
	}
      }
      if( !ppFlag ){
	if( TiltAngle<0 ){
	  MakeUnPairPlaneHitCluster( HC[layer1], CandContV[iV], honeycomb );
	  ++nV; ++iV;
	}
	if( TiltAngle==0 ){
	  MakeUnPairPlaneHitCluster( HC[layer1], CandContX[iX], honeycomb );
	  ++nX; ++iX;
	}
	if( TiltAngle>0 ){
	  MakeUnPairPlaneHitCluster( HC[layer1], CandContU[iU], honeycomb );
	  ++nU; ++iU;
	}
      }
      if( ppFlag && nh1==0 && nh2>0 ){
	if( TiltAngle<0 ){
	  MakeUnPairPlaneHitCluster( HC[layer2], CandContV[iV], honeycomb );
	  ++nV; ++iV;
	}
	if( TiltAngle==0 ){
	  MakeUnPairPlaneHitCluster( HC[layer2], CandContX[iX], honeycomb );
	  ++nX; ++iX;
	}
	if( TiltAngle>0 ){
	  MakeUnPairPlaneHitCluster( HC[layer2], CandContU[iU], honeycomb );
	  ++nU; ++iU;
	}
      }
    }

    IndexList nCombi(npp);
    IndexList nCombiV(npp);
    IndexList nCombiX(npp);
    IndexList nCombiU(npp);

    for( int i=0; i<npp; ++i ){
      nCombiV[i]=(CandContV[i]).size();
      nCombiX[i]=(CandContX[i]).size();
      nCombiU[i]=(CandContU[i]).size();
    }

#if 0
    DebugPrint( nCombiV, CandContV, func_name+" V" );
    DebugPrint( nCombiX, CandContX, func_name+" X" );
    DebugPrint( nCombiU, CandContU, func_name+" U" );
#endif

    std::vector<IndexList> CombiIndexV = MakeIndex( npp, nCombiV );
    int nnCombiV=CombiIndexV.size();
    std::vector<IndexList> CombiIndexX = MakeIndex( npp, nCombiX );
    int nnCombiX=CombiIndexX.size();
    std::vector<IndexList> CombiIndexU = MakeIndex( npp, nCombiU );
    int nnCombiU=CombiIndexU.size();

    for( int i=0; i<nnCombiV; ++i ){
      DCLocalTrack *track = MakeTrack( CandContV, CombiIndexV[i] );
      if( !track ) continue;
      if( track->GetNHit()>=3 && track->DoFitVXU() &&
	  track->GetChiSquare()<MaxChisquareVXU ){
	TrackContV.push_back(track);
      }
      else{
	delete track;
      }
    }

    for( int i=0; i<nnCombiX; ++i ){
      DCLocalTrack *track = MakeTrack( CandContX, CombiIndexX[i] );
      if( !track ) continue;
      if( track->GetNHit()>=3 && track->DoFitVXU() &&
	  track->GetChiSquare()<MaxChisquareVXU ){
	TrackContX.push_back(track);
      }
      else{
	delete track;
      }
    }

    for( int i=0; i<nnCombiU; ++i ){
      DCLocalTrack *track = MakeTrack( CandContU, CombiIndexU[i] );
      if( !track ) continue;
      if( track->GetNHit()>=3 && track->DoFitVXU() &&
	  track->GetChiSquare()<MaxChisquareVXU ){
	TrackContU.push_back(track);
      }
      else{
	delete track;
      }
    }

    // Clear Flags
    if(nV>3) ClearFlags( TrackContV );
    if(nX>3) ClearFlags( TrackContX );
    if(nU>3) ClearFlags( TrackContU );

    std::stable_sort( TrackContV.begin(), TrackContV.end(), DCLTrackComp1() );
    std::stable_sort( TrackContX.begin(), TrackContX.end(), DCLTrackComp1() );
    std::stable_sort( TrackContU.begin(), TrackContU.end(), DCLTrackComp1() );

#if 0
    DebugPrint( TrackContV, func_name+" V After Sorting." );
    DebugPrint( TrackContX, func_name+" X After Sorting." );
    DebugPrint( TrackContU, func_name+" U After Sorting." );
#endif

    // Delete Duplicated Tracks (cut chisqr>100 & flag)
    double chiV = ChisquareCutVXU;
    double chiX = ChisquareCutVXU;
    double chiU = ChisquareCutVXU;

    DeleteDuplicatedTracks( TrackContV, chiV );
    DeleteDuplicatedTracks( TrackContX, chiX );
    DeleteDuplicatedTracks( TrackContU, chiU );
    CalcTracks( TrackContV );
    CalcTracks( TrackContX );
    CalcTracks( TrackContU );

#if 0
    DebugPrint( TrackContV, func_name+" V After Delete." );
    DebugPrint( TrackContX, func_name+" X After Delete." );
    DebugPrint( TrackContU, func_name+" U After Delete." );
#endif

    int nnV=1, nnX=1, nnU=1;
    int nkV=0, nkX=0, nkU=0;
    int nnVT=1, nnXT=1, nnUT=1;
    int checkV=0, checkX=0, checkU=0;

    int cV=TrackContV.size();
    if(chiV>1.5 || cV<5 ) checkV++;
    int cX=TrackContX.size();
    if(chiX>1.5 || cX<5 ) checkX++;
    int cU=TrackContU.size();
    if(chiU>1.5 || cU<5 ) checkU++;

    std::vector<IndexList> CombiIndexSV;
    std::vector<IndexList> CombiIndexSX;
    std::vector<IndexList> CombiIndexSU;

    {
      if((nV>=3) && (cV)){
	nnV = TrackContV.size();
	nnVT = TrackContV.size();
	++nkV;
      }
      if((nV>0) && checkV){
	CombiIndexSV = MakeIndex_VXU( npp, 2, nCombiV );
	nnV = nnV +  CombiIndexSV.size();
      }
    }

    {
      if((nX>=3) && (cX)){
	nnX = TrackContX.size();
	nnXT = TrackContX.size();
	++nkX;
      }
      if((nX>0) && checkX){
	CombiIndexSX = MakeIndex_VXU( npp, 2, nCombiX );
	nnX = nnX + CombiIndexSX.size();
      }
    }

    {
      if((nU>=3) && (cU)){
	nnU = TrackContU.size();
	nnUT = TrackContU.size();
	++nkU;
      }
      if((nU>0) && checkU){
	CombiIndexSU = MakeIndex_VXU( npp, 2, nCombiU );
	nnU = nnU + CombiIndexSU.size();
      }
    }

    // double DifVXU=0.0;
    // double Av=0.0;
    // double Ax=0.0;
    // double Au=0.0;

    double chiv, chix, chiu;

#if 0
    for( int i=0; i<nnV; ++i){
      for( int j=0; j<nnX; ++j){
	for( int k=0; k<nnU; ++k){

	  chiv=-1.0,chix=-1.0,chiu=-1.0;

	  DCLocalTrack *track = new DCLocalTrack();

	  /* V Plane  */
	  if(nkV){
	    DCLocalTrack *trackV=TrackContV[i];
	    //Av=trackV->GetVXU_A();
	    chiv=trackV->GetChiSquare();
	    for( int l=0; l<(trackV->GetNHit()); ++l){
	      DCLTrackHit *hitpV=trackV->GetHit(l);
	      if(hitpV){
		track->AddHit( hitpV ) ;
	      }
	    }
	    //delete trackV;
	  }
	  if((!nkV) && (nV>0)){
	    DCLocalTrack *trackV = MakeTrack( CandContV, CombiIndexV[i] );
	    for( int l=0; l<(trackV->GetNHit()); ++l){
	      DCLTrackHit *hitpV=trackV->GetHit(l);
	      if(hitpV){
		track->AddHit( hitpV ) ;
	      }
	    }
	    delete trackV;
	  }

	  /* X Plane  */
	  if(nkX){
	    DCLocalTrack *trackX=TrackContX[j];
	    //Ax=trackX->GetVXU_A();
	    chix=trackX->GetChiSquare();
	    for( int l=0; l<(trackX->GetNHit()); ++l){
	      DCLTrackHit *hitpX=trackX->GetHit(l);
	      if(hitpX){
		track->AddHit( hitpX ) ;
	      }
	    }
	    //delete trackX;
	  }
	  if((!nkX) && (nX>0)){
	    DCLocalTrack *trackX = MakeTrack( CandContX, CombiIndexX[j] );
	    for( int l=0; l<(trackX->GetNHit()); ++l){
	      DCLTrackHit *hitpX=trackX->GetHit(l);
	      if(hitpX){
		track->AddHit( hitpX ) ;
	      }
	    }
	    delete trackX;
	  }

	  /* U Plane  */
	  if(nkU){
	    DCLocalTrack *trackU=TrackContU[k];
	    //Au=trackU->GetVXU_A();
	    chiu=trackU->GetChiSquare();
	    for( int l=0; l<(trackU->GetNHit()); ++l){
	      DCLTrackHit *hitpU=trackU->GetHit(l);
	      if(hitpU){
		track->AddHit( hitpU ) ;
	      }
	    }
	    //delete trackU;
	  }
	  if((!nkU) && (nU>0)){
	    DCLocalTrack *trackU = MakeTrack( CandContU, CombiIndexU[k] );
	    for( int l=0; l<(trackU->GetNHit()); ++l){
	      DCLTrackHit *hitpU=trackU->GetHit(l);
	      if(hitpU){
		track->AddHit( hitpU ) ;
	      }
	    }
	    delete trackU;
	  }
	  //track->SetAv(Av);
	  //track->SetAx(Ax);
	  //track->SetAu(Au);
	  //DifVXU = track->GetDifVXU();

	  track->SetChiv(chiv);
	  track->SetChix(chix);
	  track->SetChiu(chiu);

	  if(!track) continue;
	  if( track->GetNHit()>=MinNumOfHits && track->DoFit() &&
	      track->GetChiSquare()<MaxChisquare ){
	    TrackCont.push_back(track);
	  }
	  else{
	    delete track;
	  }
	}
      }
    }
#endif

    for( int i=-1; i<nnV; ++i){
      for( int j=-1; j<nnX; ++j){
	for( int k=-1; k<nnU; ++k){
	  if( ((i+j)==-2) || ((j+k)==-2) || ((k+i)==-2)) continue;

	  chiv=-1.0,chix=-1.0,chiu=-1.0;

	  DCLocalTrack *track = new DCLocalTrack();

	  /* V Plane  */
	  if(i>-1){
	    if(nkV && i<nnVT){
	      DCLocalTrack *trackV=TrackContV[i];
	      // Av=trackV->GetVXU_A();
	      chiv=trackV->GetChiSquare();

	      for( int l=0; l<(trackV->GetNHit()); ++l){
		DCLTrackHit *hitpV=trackV->GetHit(l);
		if(hitpV){
		  track->AddHit( hitpV ) ;
		}
	      }
	    }
	    if((i>=nnVT) && (nV>0)){
	      DCLocalTrack *trackV = MakeTrack( CandContV, CombiIndexSV[i-nnVT] );
	      for( int l=0; l<(trackV->GetNHit()); ++l){
		DCLTrackHit *hitpV=trackV->GetHit(l);
		if(hitpV){
		  track->AddHit( hitpV ) ;
		}
	      }
	      delete trackV ;
	    }
	    // int NHitV = track->GetNHit();
	  }

	  /* X Plane  */
	  if(j>-1){
	    if(nkX && j<nnXT){
	      DCLocalTrack *trackX=TrackContX[j];
	      // Ax=trackX->GetVXU_A();
	      chix=trackX->GetChiSquare();
	      for( int l=0; l<(trackX->GetNHit()); ++l){
		DCLTrackHit *hitpX=trackX->GetHit(l);
		if(hitpX){
		  track->AddHit( hitpX ) ;
		}
	      }
	    }
	    if((j>=nnXT) && (nX>0)){
	      DCLocalTrack *trackX = MakeTrack( CandContX, CombiIndexSX[j-nnXT] );
	      for( int l=0; l<(trackX->GetNHit()); ++l){
		DCLTrackHit *hitpX=trackX->GetHit(l);
		if(hitpX){
		  track->AddHit( hitpX ) ;
		}
	      }
	      delete trackX;
	    }
	    // int NHitX = track->GetNHit();
	  }
	  /* U Plane  */
	  if(k>-1){
	    if(nkU && k<nnUT){
	      DCLocalTrack *trackU=TrackContU[k];
	      // Au=trackU->GetVXU_A();
	      chiu=trackU->GetChiSquare();
	      for( int l=0; l<(trackU->GetNHit()); ++l){
		DCLTrackHit *hitpU=trackU->GetHit(l);
		if(hitpU){
		  track->AddHit( hitpU ) ;
		}
	      }
	    }
	    if((k>=nnUT) && (nU>0)){
	      DCLocalTrack *trackU = MakeTrack( CandContU, CombiIndexSU[k-nnUT] );
	      for( int l=0; l<(trackU->GetNHit()); ++l){
		DCLTrackHit *hitpU=trackU->GetHit(l);
		if(hitpU){
		  track->AddHit( hitpU ) ;
		}
	      }
	      delete trackU;
	    }
	    // int NHitU = track->GetNHit();
	  }

	  //track->SetAv(Av);
	  //track->SetAx(Ax);
	  //track->SetAu(Au);
	  //DifVXU = track->GetDifVXU();

	  track->SetChiv(chiv);
	  track->SetChix(chix);
	  track->SetChiu(chiu);

	  if(!track) continue;
	  if(track->GetNHit()>=MinNumOfHits && track->DoFit() &&
	     track->GetChiSquare()<MaxChisquare ){//MaXChisquare
	    TrackCont.push_back(track);
	  }
	  else{
	    delete track;
	  }
	}
      }
    }

    {

      for( int i=0; i<int(TrackContV.size()); ++i ){
	DCLocalTrack *tp=TrackContV[i];
	delete tp;
	TrackContV.erase(TrackContV.begin()+i);
      }
      for( int i=0; i<int(TrackContX.size()); ++i ){
	DCLocalTrack *tp=TrackContX[i];
	delete tp;
	TrackContX.erase(TrackContX.begin()+i);
      }
      for( int i=0; i<int(TrackContU.size()); ++i ){
	DCLocalTrack *tp=TrackContU[i];
	delete tp;
	TrackContU.erase(TrackContU.begin()+i);
      }

    }

    ClearFlags( TrackCont );

    std::stable_sort( TrackCont.begin(), TrackCont.end(), DCLTrackComp1() );

#if 1
    // Delete Tracks about  (Nhit1 > Nhit2+1) (Nhit1 > Nhit2  && chi1 < chi2)
    for( int i=0; i<int(TrackCont.size()); ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      double chi=tp->GetChiSquare();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->JoinTrack();

      for( int i2=TrackCont.size()-1; i2>i; --i2 ){
	DCLocalTrack *tp2=TrackCont[i2];
	int nh2=tp2->GetNHit(), flag=0;
	double chi2=tp2->GetChiSquare();
	for( int j=0; j<nh2; ++j )
	  if( tp2->GetHit(j)->BelongToTrack() ) ++flag;
	if((flag>=2) && ((nh==nh2) || ((nh>nh2) && (chi<chi2)))){
	  //      if((flag) && ((nh>nh2+1) || ((nh==nh2) || (nh>nh2) && (chi<chi2)))){
	  //if((nh>nh2) && (chi<chi2)){
	  delete tp2;
	  TrackCont.erase(TrackCont.begin()+i2);
	}
      }
    }
#endif

    FinalizeTrack( func_name, TrackCont, DCLTrackCompSdcInFiber(), CandCont );

    std::stable_sort( TrackCont.begin(), TrackCont.end(), DCLTrackComp() );

    // Clear Flags
    {
      int nbefore=TrackCont.size();
      for( int i=0; i<nbefore; ++i ){
	DCLocalTrack *tp=TrackCont[i];
	int nh=tp->GetNHit();
	for( int j=0; j<nh; ++j ){
	  tp->GetHit(j)->QuitTrack();
	}
      }
    }

    // Delete Duplicated Tracks

    for( int i=0; i<int(TrackCont.size()); ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->JoinTrack();

      for( int i2=TrackCont.size()-1; i2>i; --i2 ){
	DCLocalTrack *tp2=TrackCont[i2];
	int nh2=tp2->GetNHit(), flag=0;
	for( int j=0; j<nh2; ++j )
	  if( tp2->GetHit(j)->BelongToTrack() ) ++flag;
	if(flag){
	  delete tp2;
	  TrackCont.erase(TrackCont.begin()+i2);
	}
      }
    }

    // Clear Flags
    {
      int nbefore=TrackCont.size();
      for( int i=0; i<nbefore; ++i ){
	DCLocalTrack *tp=TrackCont[i];
	int nh=tp->GetNHit();
	for( int j=0; j<nh; ++j ){
	  tp->GetHit(j)->QuitTrack();
	}
      }
    }

    {
      int nn=TrackCont.size();
      for(int i=0; i<nn; ++i ){
	DCLocalTrack *tp=TrackCont[i];
	int nh=tp->GetNHit();
	for( int j=0; j<nh; ++j ){
	  int lnum = tp->GetHit(j)->GetLayer();
	  double zz = gGeom.GetLocalZ( lnum );
	  tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
	}
      }
    }

#if 0
    if(TrackCont.size()>0){
      int nn=TrackCont.size();
      hddaq::cout << func_name << ": After Deleting. #Tracks = " << nn << std::endl;
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackCont[i];
	hddaq::cout << std::setw(3) << i << " #Hits="
		    << std::setw(2) << track->GetNHit()
		    << " ChiSqr=" << track->GetChiSquare()
		    << std::endl;
	hddaq::cout << std::endl;
	for( int j=0; j<(track->GetNHit()); ++j){
	  DCLTrackHit *hit = track->GetHit(j);
	  hddaq::cout << "layer = " << hit->GetLayer()+1 << " Res = " << hit->GetResidual() << std::endl ;
	  hddaq::cout << "hitp = " << hit->GetLocalHitPos() << " calp = " << hit->GetLocalCalPos() << std::endl ;
	  hddaq::cout << "X = " << hit->GetXcal() << " Y = " << hit->GetYcal() << std::endl ;
	  hddaq::cout << std::endl;
	}
	hddaq::cout << "*********************************************" << std::endl;
      }
      hddaq::cout << std::endl;
    }
#endif

    del::ClearContainerAll( CandCont );
    del::ClearContainerAll( CandContV );
    del::ClearContainerAll( CandContX );
    del::ClearContainerAll( CandContU );

    return TrackCont.size();
  }

  //______________________________________________________________________________
  int
  LocalTrackSearchSsdIn( const std::vector<DCHitContainer>& HC,
			 std::vector<DCLocalTrack*>& TrackCont,
			 int MinNumOfHits /*=4*/ )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    std::vector<ClusterList> CandCont(NumOfLayersSsdIn);

    for( int i=0; i<NumOfLayersSsdIn; ++i ){
      MakeSsdHitCluster( HC[i+1], CandCont[i] );
    }

    IndexList nCombi(NumOfLayersSsdIn);
    for( int i=0; i<NumOfLayersSsdIn; ++i ){
      nCombi[i]=(CandCont[i]).size();
    }

    std::vector<IndexList>
      CombiIndex = MakeIndex( NumOfLayersSsdIn, nCombi );
    int nnCombi=CombiIndex.size();

    for( int i=0; i<nnCombi; ++i ){
      DCLocalTrack *track = MakeTrack( CandCont, CombiIndex[i] );
      if( !track ) continue;
      if( track->GetNHit()>=MinNumOfHits && track->DoFit()
	  ){
	TrackCont.push_back(track);
      }
      else
	delete track;
    }

    if( TrackCont.empty() ) {
      del::ClearContainerAll( CandCont );
      return 0;
    }

    FinalizeTrack( func_name, TrackCont, DCLTrackCompSsd(), CandCont, false );
    return TrackCont.size();
  }

  //______________________________________________________________________________
  int
  LocalTrackSearchSsdIn( const std::vector<SsdClusterContainer>& SsdInClCont,
			 std::vector<DCLocalTrack*>& TrackCont,
			 int MinNumOfHits /*=4*/ )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    int npp = NumOfLayersSsdIn;
    std::vector< std::vector<DCPairHitCluster*> > CandCont( npp );

    for( int i=0; i<NumOfLayersSsdIn; ++i ){
      MakeSsdHitCluster( SsdInClCont[i+1], CandCont[i] );
    }

    IndexList nCombi(npp);

    for( int i=0; i<npp; ++i ){
      nCombi[i]=(CandCont[i]).size();
    }

    std::vector<IndexList> CombiIndex = MakeIndex( npp, nCombi );

    int nnCombi = CombiIndex.size();
    for( int i=0; i<nnCombi; ++i ){
      DCLocalTrack *track = MakeTrack( CandCont, CombiIndex[i] );
      if( !track ) continue;
      if( track->GetNHit()>=MinNumOfHits && track->DoFit() ){
	TrackCont.push_back( track );
      } else {
	delete track;
      }
    }

    FinalizeTrack( func_name, TrackCont, DCLTrackCompSsd(), CandCont, false );
    return TrackCont.size();
  }

  //______________________________________________________________________________
  int
  LocalTrackSearchSsdOut( const std::vector<DCHitContainer>& HC,
			  std::vector<DCLocalTrack*>& TrackCont,
			  int MinNumOfHits /*=4*/ )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    std::vector<ClusterList> CandCont(NumOfLayersSsdOut);

    for( int i=0; i<NumOfLayersSsdOut; ++i ){
      MakeSsdHitCluster( HC[i+1], CandCont[i] );
    }

    IndexList nCombi(NumOfLayersSsdOut);
    for( int i=0; i<NumOfLayersSsdOut; ++i ){
      nCombi[i]=(CandCont[i]).size();
    }

    std::vector<IndexList> CombiIndex = MakeIndex( NumOfLayersSsdOut, nCombi );
    int nnCombi=CombiIndex.size();
    for( int i=0; i<nnCombi; ++i ){
      DCLocalTrack *track = MakeTrack( CandCont, CombiIndex[i] );
      if( !track ) continue;
      if( track->GetNHit()>=MinNumOfHits && track->DoFit() ){
	TrackCont.push_back( track );
      } else {
	delete track;
      }
    }

    FinalizeTrack( func_name, TrackCont, DCLTrackCompSsd(), CandCont, false );
    return TrackCont.size();
  }

  //______________________________________________________________________________
  int
  LocalTrackSearchSsdOut( const std::vector<SsdClusterContainer>& SsdOutClCont,
			  std::vector<DCLocalTrack*>& TrackCont,
			  int MinNumOfHits /*=4*/ )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    int npp = NumOfLayersSsdOut;
    std::vector< std::vector<DCPairHitCluster*> > CandCont( npp );
    for( int i=0; i<NumOfLayersSsdOut; ++i ){
      MakeSsdHitCluster( SsdOutClCont[i+1], CandCont[i] );
    }

    IndexList nCombi(npp);
    for( int i=0; i<npp; ++i ){
      nCombi[i]=(CandCont[i]).size();
    }

    std::vector<IndexList> CombiIndex = MakeIndex( npp, nCombi );

    int nnCombi = CombiIndex.size();
    for( int i=0; i<nnCombi; ++i ){
      DCLocalTrack *track = MakeTrack( CandCont, CombiIndex[i] );
      if( !track ) continue;
      if( track->GetNHit()>=MinNumOfHits && track->DoFit() ){
	TrackCont.push_back( track );
      } else {
	delete track;
      }
    }

    FinalizeTrack( func_name, TrackCont, DCLTrackCompSsd(), CandCont, false );
    return TrackCont.size();
  }

  //______________________________________________________________________________
  int
  LocalTrackSearchSsdInXY( const std::vector<SsdClusterContainer>& SsdInClCont,
			   std::vector<DCLocalTrack*>& TrackContX,
			   std::vector<DCLocalTrack*>& TrackContY )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    int npp = 2;
    std::vector< std::vector<DCPairHitCluster*> > CandContX(npp);
    std::vector< std::vector<DCPairHitCluster*> > CandContY(npp);

    for( int i=0; i<NumOfLayersSsdIn; ++i ){
      if( SsdInClCont[i+1].size()==0 ) continue;
      SsdCluster *cluster = SsdInClCont[i+1][0];
      if( !cluster ) continue;
      int    layer = cluster->LayerId() - PlOffsSsd;
      double tilt  = cluster->TiltAngle();
      if( tilt==0. || tilt==180. ){
	MakeSsdHitCluster( SsdInClCont[layer], CandContX[layer-1] );
      }
      else{
	MakeSsdHitCluster( SsdInClCont[layer], CandContY[layer-1] );
      }
    }

    IndexList nCombiX(npp);
    IndexList nCombiY(npp);

    for( int i=0; i<npp; ++i ){
      nCombiX[i]=(CandContX[i]).size();
      nCombiY[i]=(CandContY[i]).size();
    }

    std::vector<IndexList> CombiIndexX = MakeIndex( npp, nCombiX );
    std::vector<IndexList> CombiIndexY = MakeIndex( npp, nCombiY );
    int nnCombiX = CombiIndexX.size();
    int nnCombiY = CombiIndexY.size();

    for( int i=0; i<nnCombiX; ++i ){
      DCLocalTrack *track = MakeTrack( CandContX, CombiIndexX[i] );
      if( !track ) continue;
      if( track->DoFitVXU() ){
	TrackContX.push_back( track );
      } else {
	delete track;
      }
    }
    for( int i=0; i<nnCombiY; ++i ){
      DCLocalTrack *track = MakeTrack( CandContY, CombiIndexY[i] );
      if( !track ) continue;
      if( track->DoFitVXU() ){
	TrackContY.push_back( track );
      } else {
	delete track;
      }
    }

    FinalizeTrack( func_name, TrackContX, DCLTrackCompSsd(), CandContX, false );
    FinalizeTrack( func_name, TrackContY, DCLTrackCompSsd(), CandContY, false );
    return TrackContX.size()*TrackContY.size();
  }

  //______________________________________________________________________________
  int
  PreTrackSearchSsdXY( const std::vector<SsdClusterContainer>& SsdInClCont,
		       const std::vector<SsdClusterContainer>& SsdOutClCont,
		       std::vector<DCLocalTrack*>& TrackContX,
		       std::vector<DCLocalTrack*>& TrackContY )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    static const int npp = NumOfLayersSsdIn + NumOfLayersSsdOut;

    std::vector<ClusterList> CandContX(npp);
    std::vector<ClusterList> CandContY(npp);

    int nX=0, nY=0;
    for( int i=0; i<NumOfLayersSsdIn; ++i ){
      if( SsdInClCont[i+1].size()==0 ) continue;
      DCHit *hit   = SsdInClCont[i+1][0]->GetHit(0);
      if( !hit ) continue;
      int    layer = hit->GetLayer() - PlOffsSsd;
      double tilt  = hit->GetTiltAngle();
      if( tilt==0. || tilt==180. ){
	MakeSsdHitCluster( SsdInClCont[layer], CandContX[layer-1] );
	++nX;
      }
      else{
	MakeSsdHitCluster( SsdInClCont[layer], CandContY[layer-1] );
	++nY;
      }
    }

    for( int i=0; i<NumOfLayersSsdOut; ++i ){
      if( SsdOutClCont[i+1].size()==0 ) continue;
      DCHit *hit = SsdOutClCont[i+1][0]->GetHit(0);
      if( !hit ) continue;
      int    layer = hit->GetLayer() - PlOffsSsd;
      double tilt  = hit->GetTiltAngle();
      if( tilt==0. || tilt==180. ){
	MakeSsdHitCluster( SsdOutClCont[layer-NumOfLayersSsdIn], CandContX[layer-1] );
	++nX;
      }
      else{
	MakeSsdHitCluster( SsdOutClCont[layer-NumOfLayersSsdIn], CandContY[layer-1] );
	++nY;
      }
    }

    IndexList nCombiX(npp);
    IndexList nCombiY(npp);

    for( int i=0; i<npp; ++i ){
      int n = CandContX[i].size();
      int m = CandContY[i].size();
      nCombiX[i] = n>MaxNumOfCluster ? 0 : n;
      nCombiY[i] = m>MaxNumOfCluster ? 0 : m;
    }

#if 0
    DebugPrint( nCombiX, CandContX, func_name+" X " );
    DebugPrint( nCombiY, CandContY, func_name+" Y " );
#endif

    std::vector<IndexList> CombiIndexX = MakeIndex( npp, nCombiX );
    std::vector<IndexList> CombiIndexY = MakeIndex( npp, nCombiY );
    int nnCombiX=CombiIndexX.size();
    int nnCombiY=CombiIndexY.size();

#if 0
    hddaq::cout << "X Plane [" << nX << "] " << nnCombiX
		<< " combinations will be checked..." << std::endl;
    hddaq::cout << "Y Plane [" << nY << "] " << nnCombiY
		<< " combinations will be checked..." << std::endl;
#endif

    for( int i=0; i<nnCombiX; ++i ){
      DCLocalTrack *track = MakeTrack( CandContX, CombiIndexX[i] );
      if( !track ) continue;
      if( track->GetNHit()>=2 && track->DoFitVXU()
	  && track->GetChiSquare()<MaxChisquarePreSsdXY
	  ){
	TrackContX.push_back(track);
      }else{
	delete track;
      }
    }

    for( int i=0; i<nnCombiY; ++i ){
      DCLocalTrack *track = MakeTrack( CandContY, CombiIndexY[i] );
      if( !track ) continue;
      if( track->GetNHit()>=2 && track->DoFitVXU()
	  && track->GetChiSquare()<MaxChisquarePreSsdXY
	  ){
	TrackContY.push_back(track);
      }else{
	delete track;
      }
    }

    FinalizeTrack( func_name, TrackContX, DCLTrackComp1(), CandContX, false );
    FinalizeTrack( func_name, TrackContY, DCLTrackComp1(), CandContY, false );
    return 0;
  }

  //______________________________________________________________________________
  int
  LocalTrackSearchBcOutSdcIn( const std::vector<DCHitContainer>& BcHC,
			      const DCPairPlaneInfo *BcPpInfo,
			      const std::vector<DCHitContainer>& SdcHC,
			      const DCPairPlaneInfo *SdcPpInfo,
			      int BcNpp, int SdcNpp,
			      std::vector<DCLocalTrack*>& TrackCont,
			      int MinNumOfHits )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    const int npp = BcNpp + SdcNpp;

    std::vector<ClusterList> CandCont(npp);

    for( int i=0; i<BcNpp; ++i ){
      bool ppFlag=BcPpInfo[i].pair;
      int layer1=BcPpInfo[i].id1, layer2=BcPpInfo[i].id2;
      if(ppFlag)
	MakePairPlaneHitCluster( BcHC[layer1], BcHC[layer2],
				 BcPpInfo[i].CellSize, CandCont[i] );
      else
	MakeUnPairPlaneHitCluster( BcHC[layer1], CandCont[i] );
    }

    for( int i=0; i<SdcNpp; ++i ){
      bool ppFlag=SdcPpInfo[i].pair;
      int layer1=SdcPpInfo[i].id1, layer2=SdcPpInfo[i].id2;
      if(ppFlag)
	MakePairPlaneHitCluster( SdcHC[layer1], SdcHC[layer2],
				 SdcPpInfo[i].CellSize, CandCont[i+BcNpp] );
      else
	MakeUnPairPlaneHitCluster( SdcHC[layer1], CandCont[i+BcNpp] );
    }

    IndexList nCombi( npp );
    for( int i=0; i<npp; ++i ){
      int n = CandCont[i].size();
      nCombi[i] = n>MaxNumOfCluster ? 0 : n;
    }

#if 0
    DebugPrint( nCombi, CandCont, func_name );
#endif

    std::vector<IndexList> CombiIndex = MakeIndex( npp, nCombi );
    int nnCombi = CombiIndex.size();

#if 0
    hddaq::cout << " ===> " << nnCombi << " combinations will be checked.."
		<< std::endl;
    for( int i=0; i<nnCombi; ++i ){
      for (int j=0;j<npp;j++) {
	hddaq::cout << CombiIndex[i][j] << " ";
      }
      hddaq::cout << std::endl;
    }
#endif

    for( int i=0; i<nnCombi; ++i ){
      DCLocalTrack *track = MakeTrack( CandCont, CombiIndex[i] );
      if( !track ) continue;
      if( track->GetNHit()>=MinNumOfHits && track->DoFitBcSdc() &&
	  track->GetChiSquare()<MaxChisquare )
	TrackCont.push_back(track);
      else
	delete track;
    }

    FinalizeTrack( func_name, TrackCont, DCLTrackComp(), CandCont );

    {
      int nn=TrackCont.size();
      for(int i=0; i<nn; ++i ){
	DCLocalTrack *tp = TrackCont[i];
	int nh=tp->GetNHit();
	for( int j=0; j<nh; ++j ){
	  int lnum = tp->GetHit(j)->GetLayer();
	  double zz = gGeom.GetLocalZ( lnum );
	  if ( lnum>=113 && lnum<=124 )
	    zz -= ( zK18Target - zTarget );
	  tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
	}
      }
    }

    return TrackCont.size();
  }

  //______________________________________________________________________________
  int
  LocalTrackSearchBcOutSsdIn( const std::vector<DCHitContainer>& BcHC,
			      const DCPairPlaneInfo *BcPpInfo,
			      int BcNpp,
			      const std::vector<DCHitContainer>& SsdHC,
			      std::vector<DCLocalTrack*>& TrackCont,
			      int MinNumOfHits )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    const int npp = BcNpp+NumOfLayersSsdIn;

    std::vector<ClusterList> CandCont(npp);

    for( int i=0; i<BcNpp; ++i ){
      bool ppFlag=BcPpInfo[i].pair;
      int layer1=BcPpInfo[i].id1, layer2=BcPpInfo[i].id2;
      if(ppFlag)
	MakePairPlaneHitCluster( BcHC[layer1], BcHC[layer2],
				 BcPpInfo[i].CellSize, CandCont[i] );
      else
	MakeUnPairPlaneHitCluster( BcHC[layer1], CandCont[i] );
    }

    for( int i=0; i<NumOfLayersSsdIn; ++i ){
      MakeSsdHitCluster( SsdHC[i+1], CandCont[i+BcNpp] );
    }

    IndexList nCombi(npp);

    for( int i=0; i<npp; ++i ){
      int n = CandCont[i].size();
      nCombi[i] = n>MaxNumOfCluster ? 0 : n;
    }

#if 0
    DebugPrint(nCombi, CandCont, func_name);
#endif

    std::vector<IndexList> CombiIndex = MakeIndex( npp, nCombi );
    int nnCombi=CombiIndex.size();

    for( int i=0; i<nnCombi; ++i ){
      DCLocalTrack *track = MakeTrack( CandCont, CombiIndex[i] );
      if( !track ) continue;
      if( track->GetNHit()>=MinNumOfHits && track->DoFit()
	  && track->GetChiSquare()<MaxChisquare
	  )
	TrackCont.push_back(track);
      else
	delete track;
    }

    FinalizeTrack( func_name, TrackCont, DCLTrackComp(), CandCont );
    return TrackCont.size();
  }

  //______________________________________________________________________________
  int
  LocalTrackSearchSsdInSsdOut( const std::vector<DCHitContainer>& SsdInHC,
			       const std::vector<DCHitContainer>& SsdOutHC,
			       std::vector<DCLocalTrack*>& TrackCont,
			       int MinNumOfHits )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    const int npp = NumOfLayersSsdIn + NumOfLayersSsdIn;

    std::vector<ClusterList> CandCont(npp);

    for( int i=0; i<NumOfLayersSsdIn; ++i ){
      MakeSsdHitCluster( SsdInHC[i+1], CandCont[i] );
    }
    for( int i=0; i<NumOfLayersSsdOut; ++i ){
      MakeSsdHitCluster( SsdOutHC[i+1], CandCont[i+NumOfLayersSsdIn] );
    }

    IndexList nCombi(npp);
    for( int i=0; i<npp; ++i ){
      int n = CandCont[i].size();
      nCombi[i] = n>MaxNumOfCluster ? 0 : n;
    }

#if 0
    DebugPrint( nCombi, CandCont, func_name );
#endif

    std::vector<IndexList> CombiIndex = MakeIndex( npp, nCombi );
    int nnCombi=CombiIndex.size();

    for( int i=0; i<nnCombi; ++i ){
      DCLocalTrack *track = MakeTrack( CandCont, CombiIndex[i] );
      if( !track ) continue;
      if( track->GetNHit()>=MinNumOfHits && track->DoFit()
	  //&& track->GetChiSquare()<MaxChisquare
	  ){
	TrackCont.push_back(track);
      }
      else
	delete track;
    }

    FinalizeTrack( func_name, TrackCont, DCLTrackCompSsd(), CandCont );
    return TrackCont.size();
  }

  //______________________________________________________________________________
  int
  LocalTrackSearchSsdOutSdcIn( const std::vector<DCHitContainer>& SsdInHC,
			       const std::vector<DCHitContainer>& SsdOutHC,
			       const std::vector<DCHitContainer>& SdcHC,
			       const DCPairPlaneInfo *SdcPpInfo,
			       int SdcNpp,
			       std::vector<DCLocalTrack*>& TrackCont,
			       int MinNumOfHits )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"(DCHit)]");

    const int npp = NumOfLayersSsdIn + NumOfLayersSsdOut + SdcNpp;
    std::vector<ClusterList> CandCont(npp);

    for( int i=0; i<NumOfLayersSsdIn; ++i )
      MakeSsdHitCluster( SsdInHC[i+1], CandCont[i] );
    for( int i=0; i<NumOfLayersSsdOut; ++i )
      MakeSsdHitCluster( SsdOutHC[i+1], CandCont[i+NumOfLayersSsdIn] );

    for( int i=0; i<SdcNpp; ++i ){
      bool ppFlag    = SdcPpInfo[i].pair;
      bool honeycomb = SdcPpInfo[i].honeycomb;
      int  layer1    = SdcPpInfo[i].id1;
      int  layer2    = SdcPpInfo[i].id2;
      if(ppFlag)
	MakePairPlaneHitCluster( SdcHC[layer1], SdcHC[layer2],
				 SdcPpInfo[i].CellSize, CandCont[i+NumOfLayersSsdIn+NumOfLayersSsdOut], honeycomb );
      else
	MakeUnPairPlaneHitCluster( SdcHC[layer1], CandCont[i+NumOfLayersSsdIn+NumOfLayersSsdOut], honeycomb );
    }

    IndexList nCombi(npp);
    for( int i=0; i<npp; ++i ){
      int n = CandCont[i].size();
      nCombi[i] = n>MaxNumOfCluster ? 0 : n;
    }

#if 0
    DebugPrint( nCombi, CandCont, func_name );
#endif

    std::vector<IndexList> CombiIndex = MakeIndex( npp, nCombi );

    int nnCombi = CombiIndex.size();
    if( nnCombi==0 ){
      del::ClearContainerAll( CandCont );
      return -1;
    }

    for( int i=0; i<nnCombi; ++i ){
      DCLocalTrack *track = MakeTrack( CandCont, CombiIndex[i] );
      if( !track ) continue;
      if( track->GetNHit()>=MinNumOfHits && track->DoFit()
	  && track->GetChiSquare()<MaxChisquare
	  ){
	TrackCont.push_back(track);
      } else {
	delete track;
      }
    }

    FinalizeTrack( func_name, TrackCont, DCLTrackComp(), CandCont );
    return TrackCont.size();
  }

  //______________________________________________________________________________
  inline bool
  IsDeletionTarget( const std::vector< std::pair<int,int> >& nh,
		    std::size_t NDelete, int layer )
  {
    if( NDelete==0 ) return false;

    std::vector< std::pair<int,int> >::const_iterator itr, end;
    for( itr=nh.begin(), end=nh.end(); itr!=end; ++itr ){
      if( itr->first==0 && itr->second==layer )
	return true;
    }

    if( nh[0].second == layer )
      return true;

    return false;
  }

  //______________________________________________________________________________
  int
  LocalTrackSearchSsdOutSdcIn( const std::vector<SsdClusterContainer>& SsdInClCont,
			       const std::vector<SsdClusterContainer>& SsdOutClCont,
			       const std::vector<DCHitContainer>& SdcHC,
			       const DCPairPlaneInfo *SdcPpInfo,
			       int SdcNpp,
			       std::vector<DCLocalTrack*>& TrackCont,
			       int MinNumOfHits, bool DeleteFlag /* =false */ )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"(Cluster)]");

    const int npp = NumOfLayersSsdIn + NumOfLayersSsdOut + SdcNpp;
    std::vector<ClusterList> CandCont(npp);

    std::vector<IndexList> CombiIndex;
    std::size_t nnCombi = 0;

    // first: NHits, second: LayerId
    std::vector< std::pair<int,int> >  NhitsInfo(npp);

    std::size_t NDelete = 0;
    while( nnCombi==0 ){
      del::ClearContainerAll( CandCont );

      for( int i=0; i<NumOfLayersSsdIn; ++i ){
	if( !IsDeletionTarget( NhitsInfo, NDelete, i ) ){
	  MakeSsdHitCluster( SsdInClCont[i+1], CandCont[i] );
	}
      }

      for( int i=0; i<NumOfLayersSsdOut; ++i ){
	int n = NumOfLayersSsdIn;
	if( !IsDeletionTarget( NhitsInfo, NDelete, i+n ) ){
	  MakeSsdHitCluster( SsdOutClCont[i+1], CandCont[i+n] );
	}
      }

      for( int i=0; i<SdcNpp; ++i ){
	int n = NumOfLayersSsdIn + NumOfLayersSsdOut;
	bool ppFlag    = SdcPpInfo[i].pair;
	bool honeycomb = SdcPpInfo[i].honeycomb;
	int  layer1    = SdcPpInfo[i].id1;
	int  layer2    = SdcPpInfo[i].id2;
	DCHitContainer tempHC1(0), tempHC2(0);
	if( ppFlag ){
	  if( !IsDeletionTarget( NhitsInfo, NDelete, i+n ) ){
	    MakePairPlaneHitCluster( SdcHC[layer1], SdcHC[layer2],
				     SdcPpInfo[i].CellSize, CandCont[i+n], honeycomb );
	  }
	}else{
	  if( !IsDeletionTarget( NhitsInfo, NDelete, i+n ) ){
	    MakeUnPairPlaneHitCluster( SdcHC[layer1], CandCont[i+n], honeycomb );
	  }
	}
      }

      IndexList nCombi(npp);
      for( int i=0; i<npp; ++i ){
	int n = CandCont[i].size();
	nCombi[i] = n>MaxNumOfCluster ? 0 : n;
	NhitsInfo[i] = std::make_pair( nCombi[i], i );
      }

#if 0
      if( DeleteFlag ){
	hddaq::cout << "NDelete=" << NDelete << " ===> ";
	for( std::size_t i=0, n=NhitsInfo.size(); i<n; ++i ){
	  hddaq::cout << std::setw(4) << std::right << NhitsInfo[i].first << " ";
	}
	hddaq::cout <<std::endl;
      }
#endif

      std::stable_sort( NhitsInfo.begin(), NhitsInfo.end(),
			std::greater< std::pair<int,int> >() );

#if 0
      DebugPrint( nCombi, CandCont, func_name );
#endif

      CombiIndex = MakeIndex( npp, nCombi );

      nnCombi = CombiIndex.size();

      if( nnCombi==0 && !DeleteFlag ){
	  del::ClearContainerAll( CandCont );
	  return -1;
      }

      if( !DeleteFlag ) break;
      NDelete++;
    }// while()

#if 0
    hddaq::cout << " ===> " << std::setw(7) << nnCombi
		<< " combinations will be checked..." << std::endl;
#endif

    for( std::size_t i=0; i<nnCombi; ++i ){
      DCLocalTrack *track = MakeTrack( CandCont, CombiIndex[i] );
      if( !track ) continue;
      if( track->GetNHit()>=MinNumOfHits && track->DoFit()
	  && track->GetChiSquare()<MaxChisquare
	  ){
	TrackCont.push_back(track);
      }
      else {
	delete track;
      }
    }

    FinalizeTrack( func_name, TrackCont, DCLTrackComp(), CandCont );
    return TrackCont.size();
  }

  //______________________________________________________________________________
  int
  LocalTrackSearchSdcInSdcOut( const std::vector<DCHitContainer>& SdcInHC,
			       const DCPairPlaneInfo* SdcInPpInfo,
			       const std::vector<DCHitContainer>& SdcOutHC,
			       const DCPairPlaneInfo* SdcOutPpInfo,
			       int SdcInNpp, int SdcOutNpp,
			       std::vector<DCLocalTrack*>& TrackCont,
			       int MinNumOfHits )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    const int npp = SdcInNpp+SdcOutNpp;

    std::vector<ClusterList> CandCont(npp);

    for( int i=0; i<SdcInNpp; ++i ){
      bool ppFlag = SdcInPpInfo[i].pair;
      int  layer1 = SdcInPpInfo[i].id1;
      int  layer2 = SdcInPpInfo[i].id2;
      if( ppFlag ){
	MakePairPlaneHitCluster( SdcInHC[layer1], SdcInHC[layer2],
				 SdcInPpInfo[i].CellSize, CandCont[i] );
      } else {
	MakeUnPairPlaneHitCluster( SdcInHC[layer1], CandCont[i] );
      }
    }

    for( int i=0; i<SdcOutNpp; ++i ){
      bool ppFlag = SdcOutPpInfo[i].pair;
      int  layer1 = SdcOutPpInfo[i].id1;
      int  layer2 = SdcOutPpInfo[i].id2;
      if( ppFlag ){
	MakePairPlaneHitCluster( SdcOutHC[layer1], SdcOutHC[layer2],
				 SdcOutPpInfo[i].CellSize, CandCont[i+SdcInNpp] );
      } else {
	MakeUnPairPlaneHitCluster( SdcOutHC[layer1], CandCont[i+SdcInNpp] );
      }
    }

    IndexList nCombi( npp );
    for( int i=0; i<npp; ++i ){
	int n = CandCont[i].size();
	nCombi[i] = n>MaxNumOfCluster ? 0 : n;
    }

    std::vector<IndexList> CombiIndex = MakeIndex( npp, nCombi );
    int nnCombi = CombiIndex.size();

    for( int i=0; i<nnCombi; ++i ){
      DCLocalTrack *track = MakeTrack( CandCont, CombiIndex[i] );
      if( !track ) continue;
      if( track->GetNHit()>=MinNumOfHits &&
	  track->DoFit() &&
	  track->GetChiSquare()<MaxChisquare )
	TrackCont.push_back( track );
      else
	delete track;
    }

    FinalizeTrack( func_name, TrackCont, DCLTrackComp(), CandCont );
    return TrackCont.size();
  }

  //For MWPC
  //_____________________________________________________________________________
  int MWPCLocalTrackSearch( const std::vector<DCHitContainer>& HC,
			    std::vector<DCLocalTrack*>& TrackCont )

  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    std::vector<ClusterList> CandCont(NumOfLayersBcIn);

    for( int i=0; i<NumOfLayersBcIn; ++i ){
      MakeMWPCPairPlaneHitCluster( HC[i], CandCont[i] );
    }

    //   IndexList nCombi(NumOfLayersBcIn);
    //   for( int i=0; i<NumOfLayersBcIn; ++i ){
    //     nCombi[i]=(CandCont[i]).size();

    //     // If #Cluster>MaxNumerOfCluster,  error return

    //     if(nCombi[i]>MaxNumOfCluster){
    //       hddaq::cout << "#W  too many clusters " << func_name
    // 		<< "  layer = " << i << " : " << nCombi[i] << std::endl;
    //       del::ClearContainerAll(CandCont);
    //       return 0;
    //     }
    //   }

    //   DebugPrint(nCombi, func_name);

#if 0
    DebugPrint(nCombi, CandCont, func_name);
#endif

    const int MinNumOfHitsBcIn   = 6;
    TrackMaker trackMaker(CandCont, MinNumOfHitsBcIn, MaxCombi, MaxChisquare);
    trackMaker.MakeTracks(TrackCont);

    FinalizeTrack( func_name, TrackCont, DCLTrackComp(), CandCont );
    return TrackCont.size();
  }

  //______________________________________________________________________________
  int MWPCLocalTrackSearch( const std::vector< std::vector<DCHitContainer> >& hcList,
			    std::vector<DCLocalTrack*>& trackCont )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    for ( std::vector< std::vector<DCHitContainer> >::const_iterator
	    it    = hcList.begin(),
	    itEnd = hcList.end(); it!=itEnd; ++it )
      {
	const std::vector<DCHitContainer>& l = *it;
	std::vector<DCLocalTrack*> tc;
	MWPCLocalTrackSearch( l, tc );
	trackCont.insert(trackCont.end(), tc.begin(), tc.end());
	//       hddaq::cout << " tc " << tc.size()
	// 		<< " : " << trackCont.size() << std::endl;
      }

    ClearFlags(trackCont);
    DeleteDuplicatedTracks(trackCont);
    CalcTracks( trackCont );
    return trackCont.size();
  }
}
