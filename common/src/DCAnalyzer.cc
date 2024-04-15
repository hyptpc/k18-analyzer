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
//#include "DCDriftParamMan.hh"
//#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DCRawHit.hh"
#include "DCCluster.hh"
//#include "DCTrackSearch.hh"
#include "DebugCounter.hh"
#include "DebugTimer.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "RawData.hh"
#include "UserParamMan.hh"
#include "BLDCWireMapMan.hh"
#include "DeleteUtility.hh"
#include "LocalTrack.hh"

// #define DefStatic
// #include "DCParameters.hh"
// #undef DefStatic

#define CLUSTERTIME 0

// Tracking routine selection __________________________________________________
namespace
{
// using namespace ;
const std::string& class_name("DCAnalyzer");
const auto& gConf = ConfMan::GetInstance();
const auto& gGeom  = BLDCWireMapMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const double MaxChisquare       = 100.; // Set to be More than 30
const double MaxNumOfCluster = 5.; // Set to be Less than 30
const double MaxCombi = 1.0e2; // Set to be Less than 10^6
const double tdiff_clusters = 30;
const double tdiff_xy = 20;
const double max_slope=0.3;
const int minslayers1=3;
const int minslayers2=6;

typedef std::vector<int>               IndexList;
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

}

//______________________________________________________________________________
DCAnalyzer::DCAnalyzer(const RawData& raw_data)
  : m_raw_data(&raw_data),
    m_is_decoded(n_type),
    m_much_combi(n_type),
    m_BLC1aHC(8),
    m_BLC1bHC(8),
    m_BLC2aHC(8),
    m_BLC2bHC(8),
    m_SDCHC(8),
    m_BPCHC(8),
    m_FDCHC(6)
{
  for( int i=0; i<n_type; ++i ){
    m_is_decoded[i] = false;
    m_much_combi[i] = 0;
  }
  debug::ObjectCounter::increase(class_name);
}

DCAnalyzer::~DCAnalyzer()
{
  ClearDCHits();
  ClearDCTracks();
  ClearDCClusters();
  debug::ObjectCounter::decrease(class_name);
}
//______________________________________________________________________________
Bool_t
DCAnalyzer::DecodeRawHits(e_type k_det, const int &detid, double retiming )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  if( m_is_decoded[k_det] ){
    hddaq::cout << "#D " << func_name << " "
		<< "already decoded" << std::endl;
    return true;
  }
  ClearDCHits(detid);
  for( int layer=0; layer<8; ++layer ){
    if(detid==DetIdFDC&&layer>5) break;
    const DCRHitContainer &RHitCont=m_raw_data->GetDCRawHC(detid,layer);
    int nh = RHitCont.size();
    //    if(detid==DetIdFDC)    std::cout<<"size of raw hit container "<< detid<<"  "<<layer<<"  "<<nh<<std::endl;
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit  = RHitCont[i];
      DCHit    *hit   = new DCHit( detid, rhit->PlaneId(), rhit->WireId() );
      int       nhtdc = rhit->GetTdcSize();
      int       nhtrailing = rhit->GetTrailingSize();
      if(!hit) continue;
      for( int j=0; j<nhtdc; ++j ){
	hit->SetTdcVal( rhit->GetTdc(j) );
      }
      for( int j=0; j<nhtrailing; ++j ){
	hit->SetTdcTrailing( rhit->GetTrailing(j) );
      }
      if( hit->CalcDCObservables(retiming) ){
	GetDCHC(detid,layer).push_back(hit);
      }
      else{
	delete hit;
      }
    }
  }
  //  m_is_decoded[k_det] = true;
  return true;
}


//______________________________________________________________________________
Bool_t
DCAnalyzer::DecodeRawHits(double retiming_t0, double retiming_def )
{
  ClearDCHits();
  ClearDCTracks();
#if E15
  DecodeRawHits(k_FDC, DetIdFDC );
#endif
#if E62
  DecodeRawHits(k_SDC, DetIdSDC );
  //  DecodeRawHits(k_FDC, DetIdFDC );
#elif E73_2024
  DecodeRawHits(k_BPC1, DetIdBPC1, retiming_def );
  DecodeRawHits(k_BPC2, DetIdBPC2, retiming_def );
#else
  DecodeRawHits(k_BPC, DetIdBPC, retiming_def );
#endif
  DecodeRawHits(k_BLC1a, DetIdBLC1a, retiming_t0 );
  DecodeRawHits(k_BLC1b, DetIdBLC1b, retiming_t0 );
  DecodeRawHits(k_BLC2a, DetIdBLC2a, retiming_t0 );
  DecodeRawHits(k_BLC2b, DetIdBLC2b, retiming_t0 );
  return true;
}

//______________________________________________________________________________
void
DCAnalyzer::ClearDCHits( const int &detid )
{
  switch(detid){
  case DetIdSDC:
    del::ClearContainerAll( m_SDCHC );
    break;
  case DetIdBLC1a:
    del::ClearContainerAll( m_BLC1aHC );
    break;
  case DetIdBLC1b:
    del::ClearContainerAll( m_BLC1bHC );
    break;
  case DetIdBLC2a:
    del::ClearContainerAll( m_BLC2aHC );
    break;
  case DetIdBLC2b:
    del::ClearContainerAll( m_BLC2bHC );
    break;
  case DetIdFDC:
    del::ClearContainerAll( m_FDCHC );
    break;
  case DetIdBPC:
    del::ClearContainerAll( m_BPCHC );
    break;
  }
}
//______________________________________________________________________________
void
DCAnalyzer::ClearDCTracks( const int &detid )
{
  switch(detid){
  case DetIdSDC:
    del::ClearContainer( m_SDCTC );
    break;
  case DetIdBLC1a:
    del::ClearContainer( m_BLC1aTC );
    break;
  case DetIdBLC1b:
    del::ClearContainer( m_BLC1bTC );
    break;
  case DetIdBLC1:
    del::ClearContainer( m_BLC1TC );
    break;
  case DetIdBLC2a:
    del::ClearContainer( m_BLC2aTC );
    break;
  case DetIdBLC2b:
    del::ClearContainer( m_BLC2bTC );
    break;
  case DetIdBLC2:
    del::ClearContainer( m_BLC2TC );
    break;
  case DetIdFDC:
    del::ClearContainer( m_FDCTC );
    break;
  case DetIdBPC:
    del::ClearContainer( m_BPCTC );
    break;
  }
}

//______________________________________________________________________________
void
DCAnalyzer::ClearDCHits( void)
{
  del::ClearContainerAll( m_SDCHC );
  del::ClearContainerAll( m_FDCHC );
  del::ClearContainerAll( m_BPCHC );
  del::ClearContainerAll( m_BLC1aHC );
  del::ClearContainerAll( m_BLC1bHC );
  del::ClearContainerAll( m_BLC2aHC );
  del::ClearContainerAll( m_BLC2bHC );
}

//______________________________________________________________________________
void
DCAnalyzer::ClearDCTracks( void)
{
  del::ClearContainer( m_SDCTC );
  del::ClearContainer( m_FDCTC );
  del::ClearContainer( m_BPCTC );
  del::ClearContainer( m_BLC1aTC );
  del::ClearContainer( m_BLC1bTC );
  del::ClearContainer( m_BLC2aTC );
  del::ClearContainer( m_BLC2bTC );
  del::ClearContainer( m_BLC1TC );
  del::ClearContainer( m_BLC2TC );
}
//______________________________________________________________________________
void
DCAnalyzer::ClearDCClusters(void)
{
  DCClusterListContainer::iterator itr, itr_end=m_DCCC.end();
  for( itr=m_DCCC.begin(); itr!=itr_end; ++itr )
    del::ClearContainerAll( itr->second );
  m_DCCC.clear();
}

//______________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcDCHits( std::vector<DCHitContainer>& cont,
			  Bool_t applyRecursively )
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
Bool_t
DCAnalyzer::TrackSearchAll(){
  double maxsub=999;
  //  MakePairsAll(maxsub);
  return true;
}
//______________________________________________________________________________
int
DCAnalyzer::TrackSearch( const int &detid, const Bool_t &TIMING, const int &debug ){
  int ndet=1;
  int detlist[2]={detid,detid};
  int minslayers=minslayers1;
  if(detid==DetIdBLC1){
    ndet=2;
    minslayers=minslayers2;
    detlist[0]=DetIdBLC1a;
    detlist[1]=DetIdBLC1b;
  }
  else if(detid==DetIdBLC2){
    ndet=2;
    minslayers=minslayers2;
    detlist[0]=DetIdBLC2a;
    detlist[1]=DetIdBLC2b;
  }
  else if(detid==DetIdBPC0){
    ndet=2;
    minslayers=minslayers2;
    detlist[0]=DetIdBPC1;
    detlist[1]=DetIdBPC2;
  }

  int status=0;
  LocalTrackContainer tmpcontainer[2];
  LocalTrackContainer tmpcontainer2[2];
  LocalTrackContainer tmpcontainer3;
  for(int xy=0;xy<2;xy++){
    if(debug>0) std::cout<<"============= "<<detid<<"  "<<xy<<std::endl;
    int ncomb=1;
    for(int idet=0;idet<ndet;idet++){
      int nlay=GetNClusterContainers(detlist[idet],xy);
      for(int ith=0;ith<nlay;ith++){
	int nc=GetNClusters(detlist[idet],xy,ith);
	ncomb*=nc;
	if(debug>0){
	  std::cout<<ith<<"  / "<<nlay<<"  "<<nc<<"  "<<ncomb<<std::endl;
	}
      }
    }
#if 0
    if(detid==DetIdBLC2b){
      std::cout<<"ncomb: "<<ncomb<<std::endl;
      if(ncomb==0){
	Print(detid);
	Print(DetIdBPC);
      }
    }
#endif
    for(int icomb=0;icomb<ncomb;icomb++){
      status=2;
      int tmp=icomb;
      LocalTrack *track = new LocalTrack(detlist[0]);
      for(int idet=0;idet<ndet;idet++){
	int nlay=GetNClusterContainers(detlist[idet],xy);
	for(int ith=0;ith<nlay;ith++){
	  int nc=GetNClusters(detlist[idet],xy,ith);
	  int ic=tmp%nc;
	  tmp=tmp/nc;
	  DCCluster* cl=GetCluster(detlist[idet],xy,ith,ic);
	  if(cl->ctime()>-500){
	    track->AddClusterTime(cl->ctime());
	  }
	  for(int ih=0;ih<cl->nhit();ih++){
	    track->AddHit(cl->hit(ih),cl->nth(ih),xy);
	  }
	}
      }
      if( track->nhit(xy)<minslayers){
	delete track;
	continue;
      }
      status=3;
      if( TIMING && track->GetTrackTimeRMS() > tdiff_clusters ) {
	delete track;
	continue;
      }
      status=4;
      if(!track->LeastSquareFit(xy)){
	delete track;
	continue;
      }
      if(track->chi2(xy)>1e3){
	delete track;
	continue;
      }
      double slope=track->b();
      if(xy) slope=track->e();
      if(slope>max_slope){
	delete track;
	continue;
      }
      if(debug>0){
	std::cout<<xy<<"  "
		 <<track->nhit(xy)<<"  "
		 <<track->chi2(xy)<<std::endl;
      }
      status=5;
      tmpcontainer[xy].push_back(track);
    }
    int tmpnbad=0;
    while((tmpnbad+tmpcontainer2[xy].size())<tmpcontainer[xy].size()){
      double tmpchi=9999;
      int tmpbest=-1;
      for(int itr=0;itr<tmpcontainer[xy].size();itr++){
	LocalTrack *tmptr=tmpcontainer[xy][itr];
	if(debug>0) std::cout<<xy<<"  "<<itr<<"  "<<tmptr->chi2(xy)<<std::endl;
	// if( (tmpchi < minchi &&  (tmptrrms < mintrrms || tmptrrms <10 || mintrrms > 10 || option.Contains("notiming") ) )
	//     || ( tmptrrms < mintrrms && (tmpchi < minchi + 5 ) ) ){
	if(tmptr->isgood() && tmptr->chi2(xy)<tmpchi){
	  tmpbest=itr;
	  tmpchi=tmptr->chi2(xy);
	}
      }
      if(tmpbest<0){
	break;
      }else{
	LocalTrack *tmptr=tmpcontainer[xy][tmpbest];
	tmpcontainer2[xy].push_back(tmptr);
	LocalTrackContainer::iterator it;
	for(it=tmpcontainer[xy].end()-1;it!=tmpcontainer[xy].begin()-1;--it){
	  LocalTrack* tmptr2 = *it;
	  if(!tmptr2->isgood()) continue;
	  if( tmptr->CompareTrackHit(tmptr2)){
	    tmptr2->SetBad();
	    tmpnbad++;
	    // if(debug>0) std::cout<<"erase track"<<std::endl;
	    // tmpcontainer[xy].erase(it);
	    // delete tmptr2;
	    // if(debug>0) std::cout<<"erase track done"<<std::endl;
	  }
	}
      }
    }
    if(debug>0) std::cout<<"tmpcontainer2 size: "<<tmpcontainer2[xy].size()<<std::endl;
  }
  if(tmpcontainer2[0].size()>0 && tmpcontainer2[1].size()>0)
    status=6;


  if(status==6
     &&tmpcontainer2[0].size()<10
     &&tmpcontainer2[1].size()<10){
    status=7;
    for(int itrx=0;itrx<tmpcontainer2[0].size();itrx++){
      for(int itry=0;itry<tmpcontainer2[1].size();itry++){
	LocalTrack *tmptrx=tmpcontainer2[0][itrx];
	LocalTrack *tmptry=tmpcontainer2[1][itry];
	if(TIMING&&TMath::Abs(tmptrx->GetTrackTime()-tmptry->GetTrackTime())>tdiff_xy){
	  if(debug>0) std::cout<<"tdiffxy large: "
			       <<tmptrx->GetTrackTime()<<"  "
			       <<tmptry->GetTrackTime()<<std::endl;
	  continue;
	}
	LocalTrack *track=new LocalTrack(detid);
	for(int ihit=0;ihit<tmptrx->nhit(0);ihit++){
	  track->AddHit(tmptrx->hit(0,ihit),0);
	}
	for(int ihit=0;ihit<tmptrx->nclustertimes();ihit++){
	  track->AddClusterTime(tmptrx->clustertime(ihit));
	}
	for(int ihit=0;ihit<tmptry->nhit(1);ihit++){
	  track->AddHit(tmptry->hit(1,ihit),1);
	}
	for(int ihit=0;ihit<tmptry->nclustertimes();ihit++){
	  track->AddClusterTime(tmptry->clustertime(ihit));
	}
	if( !track->DoFit() || track->chi2all()>1e3 ){
	  //	  std::cout<<"failed in fitting "<<track->chi2all()<<std::endl;
	  delete track;
	  continue;
	}
	tmpcontainer3.push_back( track );
      } //itry
    } //itrx
    if(debug>0) std::cout<<"tmpcontainer3 size: "<<tmpcontainer3.size()<<std::endl;
    int tmpnbad=0;
    while((tmpnbad+GetNTracks(detid))<tmpcontainer3.size()){
      status=8;
      double tmpchi=9999;
      int tmpbest=-1;
      for(int itr=0;itr<tmpcontainer3.size();itr++){
	LocalTrack *tmptr=tmpcontainer3[itr];
	if(debug>0){
	  std::cout<<"final "<<itr
		   <<std::setw(10)<<tmptr->chi2all()
		   <<std::setw(10)<<tmptr->GetTrackTime()
		   <<std::setw(10)<<tmptr->GetTrackTimeRMS()
		   <<std::endl;
	}
	// if( (tmpchi < minchi &&  (tmptrrms < mintrrms || tmptrrms <10 || mintrrms > 10 || option.Contains("notiming") ) )
	//     || ( tmptrrms < mintrrms && (tmpchi < minchi + 5 ) ) ){
	if(tmptr->isgood() && tmptr->chi2all()<tmpchi){
	  tmpbest=itr;
	  tmpchi=tmptr->chi2all();
	}
      }
      if(tmpbest<0){
	break; //del::ClearContainer( tmpcontainer3 );
      }else{
	LocalTrack *tmptr=tmpcontainer3[tmpbest];
	status=1;
	GetTC(detid).push_back(new LocalTrack(tmptr));
	LocalTrackContainer::iterator it;
	for(it=tmpcontainer3.end()-1;it!=tmpcontainer3.begin()-1;--it){
	  LocalTrack* tmptr2 = *it;
	  if( !tmptr2->isgood() ) continue;
	  if( tmptr->CompareTrackHit(tmptr2) ){
	    tmptr2->SetBad();
	    tmpnbad++;
	    // tmpcontainer3.erase(it);
	    // delete tmptr2;
	  }
	}
      }
    }
  }
  del::ClearContainer( tmpcontainer[0] );
  del::ClearContainer( tmpcontainer[1] );
  // tmpcontainer2[0].clear();
  // tmpcontainer2[1].clear();
  // del::ClearContainer( tmpcontainer2[0] );
  // del::ClearContainer( tmpcontainer2[1] );
  del::ClearContainer( tmpcontainer3 );
  return status;
}
//______________________________________________________________________________
int
DCAnalyzer::DeleteBadClusters( const int &detid, const int &xy,const int &icc)
{
  int good=0;
  DCClusterContainer &dccc=GetDCCL(detid,xy)[icc];
  //  int tmpn=dccc.size();
  DCClusterContainer::iterator ittr;
  for( ittr=dccc.end()-1; ittr!=dccc.begin()-1; --ittr ){
    if(!gUser.Check(Form("DCCL%d",detid),(*ittr)->time()) ) dccc.erase(ittr);
    else good++;
  }
  //  std::cout<<"before: "<<dccc.size()<<std::endl;
  // std::cout<<"after: "<<dccc.size()<<std::endl;
  // std::cout<<"detid,xy,icc,good:  "<<detid<<"  "<<xy<<"  "<<icc<<"  "<<good<<std::endl;
  return good;
}
//______________________________________________________________________________
Bool_t
DCAnalyzer::MakePairs( const int &detid, const Bool_t &isMC, const double &maxsub )
{
  int nlayers=8;
  int id1=detid;
  int id2=-1;
  if(detid==DetIdBLC1){
    id1=DetIdBLC1a;    id2=DetIdBLC1b;
  }
  if(detid==DetIdBLC2){
    id1=DetIdBLC2a;    id2=DetIdBLC2b;
  }
  for(int i=0;i<nlayers;i+=2){
    int layer1=i;
    int layer2=i+1;
    MakePairs(id1,layer1,layer2,isMC,maxsub);
    if(id2>0) MakePairs(id2,layer1,layer2,isMC,maxsub);
  }
  return true;
}
//______________________________________________________________________________
Bool_t
DCAnalyzer::MakePairsAll(const Bool_t &isMC, const double &maxsub )
{
  ClearDCClusters();
  MakePairs(DetIdBLC1a,isMC,maxsub);
  MakePairs(DetIdBLC1b,isMC,maxsub);
  MakePairs(DetIdBLC2a,isMC,maxsub);
  MakePairs(DetIdBLC2b,isMC,maxsub);
#if E62
  MakePairs(DetIdSDC,isMC,maxsub);
  //     if(layer1<6){
  //       MakePairs(DetIdFDC,layer1,layer2,maxsub);
  //     }
#elif E73_2024
  MakePairs(DetIdBPC1,isMC,maxsub);
  MakePairs(DetIdBPC2,isMC,maxsub);
#else
  MakePairs(DetIdBPC,isMC,maxsub);
#endif
  return true;
#if 0
  for(int xy=0;xy<2;xy++){
    int nlay=GetNClusterContainers(detid,xy);
    for(int ith=0;ith<nlay;ith++){
      int good=DeleteBadClusters(detid,xy,ith);
      if(good!=1) return false;
    }
  }
#endif
}
Bool_t
DCAnalyzer::MakePairs( const int &detid, const int &layer1, const int &layer2, const Bool_t &isMC, const double &maxsub)
{
  const DCHitContainer &hc1=GetDCHC(detid,layer1);
  const DCHitContainer &hc2=GetDCHC(detid,layer2);
  double cellsize=TMath::Abs(gGeom.GetWireMap(detid,layer1)->GetdXY());
  int nh1=hc1.size(), nh2=hc2.size();
  if(nh1>16) nh1=0;
  if(nh2>16) nh2=0;
  std::vector<int> UsedFlag(nh2,0);
  DCClusterContainer Cont;
  for( int i1=0; i1<nh1; ++i1 ){
    DCHit *hit1=hc1[i1];
    TVector3 wp1=hit1->GetWirePosition();
    int multi1 = hit1->GetDriftLengthSize();
    Bool_t flag=false;
    for( int i2=0; i2<nh2; ++i2 ){
      DCHit *hit2=hc2[i2];
      TVector3 wp2=hit2->GetWirePosition();
      //      std::cout<<(wp1-wp2).Perp()<<"  "<<cellsize<<"  "<<hit1->IsWithinTotRange()<<"  "<<hit2->IsWithinTotRange()<<std::endl;
      if( (wp1-wp2).Perp()<cellsize ){
	int multi2 = hit2->GetDriftLengthSize();
	for ( int m1=0; m1<multi1; ++m1 ) {
	  if( !hit1->IsWithinTotRange(m1) ) 	    continue;
	  if( !hit1->IsWithinDtRange(m1) )	    continue;
	  for ( int m2=0; m2<multi2; ++m2 ) {
	    if( !hit2->IsWithinTotRange(m2) )	      continue;
	    if( !hit2->IsWithinDtRange(m2) )	      continue;
	    if( TMath::Abs(hit1->GetDriftTime(m1)-hit2->GetDriftTime(m2))> maxsub) continue;
	    DCCluster *cluster = new DCCluster();
	    cluster->SetHit(hit1,m1);
	    cluster->SetHit(hit2,m2);
	    Cont.push_back( cluster );
	    flag=true; ++UsedFlag[i2];
	  }
	}
      }
    }
#if 1 // without hit in the 2nd layer
    if(!flag){
      for (int m1=0; m1<multi1; m1++) {
	if( !(hit1->IsWithinTotRange(m1)) ) continue;
	if( !(hit1->IsWithinDtRange(m1)) ) continue;
	DCCluster *cluster = new DCCluster();
	cluster->SetHit(hit1,m1);
	Cont.push_back( cluster );
      }
    }
#endif
  }
#if 1 // only 2nd layer hit without associated 1st layer hit
  for( int i2=0; i2<nh2; ++i2 ){
    if( UsedFlag[i2]==0 ) {
      DCHit *hit2=hc2[i2];
      int multi2 = hit2->GetDriftLengthSize();
      for (int m2=0; m2<multi2; m2++) {
	if( !(hit2->IsWithinTotRange(m2)) ) continue;
	if( !(hit2->IsWithinDtRange(m2)) ) continue;
	DCCluster *cluster = new DCCluster();
	cluster->SetHit(hit2,m2);
	Cont.push_back( cluster );
      }
    }
  }
#endif
  if(Cont.size()>30){
    std::cout<<__PRETTY_FUNCTION__<<" too many clusters "<<Cont.size()<<std::endl;
    del::ClearContainer( Cont );
  }
  for( int i=0;i<Cont.size();++i )  Cont[i]->Calc(isMC);
  int xy=gGeom.GetWireMap(detid,layer1)->GetXY();
  int key=MakeKey(detid,xy);
  m_DCCC[key].push_back(Cont);
  return true;
}
//______________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcDCHits( Bool_t applyRecursively )
{
  return true;
}

//______________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcAll()
{
  ReCalcDCHits();
  return true;
}

void DCAnalyzer::Print(const int &detid){
  std::cout<<"[DetId: "<<detid<<" ]"<<std::endl;
  for(int i=0;i<8;i++){
    const DCHitContainer &hc=GetDCHC(detid,i);
    int nhit=hc.size();
    std::cout<<"layer, nhit "<<i<<"  "<<nhit<<std::endl;
    for(int j=0; j<nhit;j++){
      std::cout<<"layer "<<i<<"  "<<j<<" / "<<nhit<<"  ";
      hc[j]->GetWirePosition().Print();
    }
  }
  for(int xy=0;xy<2;xy++){
    int nlay=GetNClusterContainers(detid,xy);
    for(int ith=0;ith<nlay;ith++){
      int nc=GetNClusters(detid,xy,ith);
      std::cout<<"=== XY: "<<xy<<"  ith: "<<ith<<"  num of clusters = "<<nc<<std::endl;
      for(int ic=0;ic<nc;ic++){
	DCCluster* cl=GetCluster(detid,xy,ith,ic);
	std::cout<<ic<<"/"<<nc<<", nhit: "<<cl->nhit()<<", time:"<<cl->GetTime()<<"  ";
	cl->pos().Print();
      }
    }
  }
}
