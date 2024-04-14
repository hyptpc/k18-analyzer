/**
 *  file: CDCAnalyzer.cc
 *  date: 2017.04.10
 *
 */

#include "CDCAnalyzer.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include "ConfMan.hh"
#include "DCHit.hh"
#include "DCRawHit.hh"
#include "DCCluster.hh"
#include "DebugCounter.hh"
#include "DebugTimer.hh"
#include "RawData.hh"
#include "UserParamMan.hh"
#include "CDCWireMapMan.hh"
#include "CDSTrack.hh"
#include "DeleteUtility.hh"
#include "MathTools.hh"

#define DEBUG 0
// Tracking routine selection __________________________________________________
namespace
{
  // using namespace ;
  const double mm=0.1;
  const double cm=10*mm;
  const std::string& class_name("CDCAnalyzer");
  const ConfMan&      gConf = ConfMan::GetInstance();
  const CDCWireMapMan&  gGeom  = CDCWireMapMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const double& cellsize=8.45*2*mm;
  const double& cellfac=1.2;
  const double& selected_chi2=1000;
  //______________________________________________________________________________
  int vector_finder(std::vector<int> vec, int number) {
    auto itr = std::find(vec.begin(), vec.end(), number);
    size_t index = std::distance( vec.begin(), itr );
    if (index != vec.size()) {
      return 1;
    }
    else {
      return 0;
    }
  }
}
//______________________________________________________________________________
CDCAnalyzer::CDCAnalyzer( void )
  : m_is_decoded(n_type),
    m_much_combi(n_type),
    m_CDCHC(15)
{
  for( int i=0; i<n_type; ++i ){
    m_is_decoded[i] = false;
    m_much_combi[i] = 0;
  }
  debug::ObjectCounter::increase(class_name);
}

CDCAnalyzer::~CDCAnalyzer( void )
{
  ClearDCHits();
  ClearDCTracks();
  ClearDCClusters();
  debug::ObjectCounter::decrease(class_name);
}
//______________________________________________________________________________
void
CDCAnalyzer::ClearDCTracks( )
{
  del::ClearContainer( m_CDCTC );
}
void
CDCAnalyzer::ClearDCClusters(void)
{
  // DCClusterList::iterator itr, itr_end=m_DCCC.end();
  // for( itr=m_DCCC.begin(); itr!=itr_end; ++itr )
  //   del::ClearContainerAll( itr );
  // m_DCCC.clear();  
  del::ClearContainerAll(m_DCCC);
}

//______________________________________________________________________________
bool
CDCAnalyzer::DecodeRawHits( RawData *rawData, e_type k_det, const int &detid, double retiming )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  if( m_is_decoded[k_det] ){
    hddaq::cout << "#D " << func_name << " "
		<< "already decoded" << std::endl;
    return true;
  }
  ClearDCHits(detid);
  for( int layer=0; layer<118; ++layer ){
    const DCRHitContainer &RHitCont=rawData->GetDCRawHC(detid,layer);
    int nh = RHitCont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit  = RHitCont[i];
      int asdnum=rhit->PlaneId();
      int asdch=rhit->WireId();
      int cdclayer=0;
      int cdcwire=0;
      if(!gGeom.GetLayerWire(asdnum,asdch,cdclayer,cdcwire)) continue;
      DCHit    *hit   = new DCHit( detid, cdclayer,cdcwire );
      if(!hit) continue;
      int       nhtdc = rhit->GetTdcSize();
      int       nhtrailing = rhit->GetTrailingSize();
      for( int j=0; j<nhtdc; ++j ){
	hit->SetTdcVal( rhit->GetTdc(j) );
      }      
      for( int j=0; j<nhtrailing; ++j ){
	hit->SetTdcTrailing( rhit->GetTrailing(j) );
      }      
      if( hit->CalcDCObservables(retiming) ){
	switch(detid){
	case DetIdCDC:
	  m_CDCHC[cdclayer].push_back(hit);
	  break;
	default:
	  std::cout<<"E# "<<func_name<<" invalid detector id "<< detid<<std::endl;
	  return false;
	}
      }
      else{
      	delete hit;
      }
    }
  }
  //  m_is_decoded[k_det] = true;
  return true;
}

  bool
CDCAnalyzer::TrackSearch()
{
  double maxsub=999;
  int nslayers=7;
  for(int i=0;i<nslayers;i++){
    MakeClusters(i,maxsub);
#if DEBUG
    std::cout<<"# of cluster in slayer "<<i<<" is : "<<GetNClusters(i)<<std::endl;
#endif
  }
  if(gUser.GetParameter("CDSFIELD")==0){
    return (FindCircleTrackCandidates(true)
	    && LineAxialFitting() 
	    && FindStereoHits()
	    && LineFitting() 
	    );    
  }
  return (FindCircleTrackCandidates()
	  && CircleFitting() 
	  && FindStereoHits()
	  && HelixFitting() 
	  );
}
bool
CDCAnalyzer::MakeClusters( const int &slayer, const double &maxsub)
{
#if DEBUG
  std::cout<<__PRETTY_FUNCTION__<<"  "<<slayer<<std::endl;
#endif
  int layer1=0;
  if(0<slayer && slayer<7) layer1=slayer*2+1;
  int layer2=layer1+1; 
  int layer3=layer1+2;
  const DCHitContainer &hc1=GetDCHC(layer1);
  const DCHitContainer &hc2=GetDCHC(layer2);
  const DCHitContainer &hc3=GetDCHC(layer3);
  int nh1=hc1.size(), nh2=hc2.size(), nh3=0;
  if(slayer==0) nh3=hc3.size();
  // if(nh1>50) nh1=0;
  // if(nh2>50) nh2=0;
  // if(nh3>50) nh3=0;
  std::vector<int> UsedFlag2(nh2+1,0);
  std::vector<int> UsedFlag3(nh3+1,0);  
  DCClusterContainer Cont;
  // with layer1 hit
  for( int i1=0; i1<nh1; ++i1 ){
    DCHit *hit1=hc1[i1];
    TVector3 wp1=hit1->GetWirePosition();
    int multi1 = hit1->GetDriftLengthSize();
    bool flag2=false;
    for( int i2=0; i2<nh2; ++i2 ){
      DCHit *hit2=hc2[i2];
      TVector3 wp2=hit2->GetWirePosition();
      if( (wp1-wp2).Perp()>cellsize*cellfac ) continue;
      int multi2 = hit2->GetDriftLengthSize();
      for ( int m1=0; m1<multi1; ++m1 ) {
	if( !hit1->IsWithinTotRange(m1) ) 	    continue;
	if( !hit1->IsWithinDtRange(m1) )	    continue;
	for ( int m2=0; m2<multi2; ++m2 ) {
	  if( !hit2->IsWithinTotRange(m2) )	      continue;
	  if( !hit2->IsWithinDtRange(m2) )	      continue;
	  if( TMath::Abs(hit1->GetDriftTime(m1)-hit2->GetDriftTime(m2)) > maxsub) continue;
	  flag2=true; ++UsedFlag2[i2];
	  bool flag3=false;
	  if(slayer==0){
	    for( int i3=0; i3<nh3; ++i3 ){
	      DCHit *hit3=hc3[i3];
	      TVector3 wp3=hit3->GetWirePosition();
	      if( (wp3-wp2).Perp()>cellsize*cellfac ) continue;
	      int multi3 = hit3->GetDriftLengthSize();
	      for ( int m3=0; m3<multi3; ++m3 ) {
		if( !hit3->IsWithinTotRange(m3) )	      continue;
		if( !hit3->IsWithinDtRange(m3) )	      continue;
		if( TMath::Abs(hit3->GetDriftTime(m3)-hit2->GetDriftTime(m2))> maxsub) continue;
		flag3=true; ++UsedFlag3[i3];
		DCCluster *cluster = new DCCluster();
		cluster->SetHit(hit1,m1);
		cluster->SetHit(hit2,m2);
		cluster->SetHit(hit3,m3);
		Cont.push_back( cluster );
	      } // m3
	    }//layer3
	  }		      
	  if(slayer!=0||!flag3)
	    {
	      DCCluster *cluster = new DCCluster();
	      cluster->SetHit(hit1,m1);
	      cluster->SetHit(hit2,m2);
	      Cont.push_back( cluster );
	    }// layer1 & layer2 not layer3
	} // m2
      } // m1
    } // layer2
    if(!flag2)
      {
	bool flag3=false;
	if(slayer==0)
	  {
	    for( int i3=0; i3<nh3; ++i3 ){
	      DCHit *hit3=hc3[i3];
	      TVector3 wp3=hit3->GetWirePosition();
	      if( (wp3-wp1).Perp()>2*cellsize*cellfac ) continue;
	      int multi3 = hit3->GetDriftLengthSize();
	      for ( int m1=0; m1<multi1; ++m1 ) {
		if( !hit1->IsWithinTotRange(m1) ) 	    continue;
		if( !hit1->IsWithinDtRange(m1) )	    continue;
		for ( int m3=0; m3<multi3; ++m3 ) {
		  if( !hit3->IsWithinTotRange(m3) )	      continue;
		  if( !hit3->IsWithinDtRange(m3) )	      continue;
		  if( TMath::Abs(hit3->GetDriftTime(m3)-hit1->GetDriftTime(m1))> maxsub) continue;
		  flag3=true; ++UsedFlag3[i3];
		  DCCluster *cluster = new DCCluster();
		  cluster->SetHit(hit1,m1);
		  cluster->SetHit(hit3,m3);
		  Cont.push_back( cluster );
		} // m3
	      }//m1
	    }
	  }
	if(slayer!=0||!flag3)
	  {
	    for ( int m1=0; m1<multi1; ++m1 ) {
	      if( !hit1->IsWithinTotRange(m1) ) 	    continue;
	      if( !hit1->IsWithinDtRange(m1) )	    continue;
	      DCCluster *cluster = new DCCluster();
	      cluster->SetHit(hit1,m1);
	      Cont.push_back( cluster );
	    }
	  }
      } //layer1 not layer2
  }
  // without layer1
  for( int i2=0; i2<nh2; ++i2 ){
    if( UsedFlag2[i2]!=0 ) continue;
    DCHit *hit2=hc2[i2];
    TVector3 wp2=hit2->GetWirePosition();
    int multi2 = hit2->GetDriftLengthSize();
    bool flag3=false;
    if(slayer==0){
      for( int i3=0; i3<nh3; ++i3 ){
	if( UsedFlag3[i3]!=0 ) continue;
	DCHit *hit3=hc3[i3];
	TVector3 wp3=hit3->GetWirePosition();
	if( (wp3-wp2).Perp()>cellsize*cellfac ) continue;
	int multi3 = hit3->GetDriftLengthSize();
	for ( int m2=0; m2<multi2; ++m2 ) {
	  if( !hit2->IsWithinTotRange(m2) ) 	    continue;
	  if( !hit2->IsWithinDtRange(m2) )	    continue;
	  for ( int m3=0;m3<multi3; ++m3 ) {
	    if( !hit3->IsWithinTotRange(m3) )	      continue;
	    if( !hit3->IsWithinDtRange(m3) )	      continue;
	    if( TMath::Abs(hit3->GetDriftTime(m3)-hit2->GetDriftTime(m2))> maxsub) continue;
	    flag3=true; ++UsedFlag3[i3];
	    DCCluster *cluster = new DCCluster();
	    cluster->SetHit(hit2,m2);
	    cluster->SetHit(hit3,m3);
	    Cont.push_back( cluster );
	  }
	}
      }
    }
    if(slayer!=0||!flag3)
      {
	for ( int m2=0; m2<multi2; ++m2 ) {
	  if( !hit2->IsWithinTotRange(m2) ) 	    continue;
	  if( !hit2->IsWithinDtRange(m2) )	    continue;
	  DCCluster *cluster = new DCCluster();
	  cluster->SetHit(hit2,m2);
	  Cont.push_back( cluster );
	}
      }// layer2 not layer3
  } // layer2
  if(slayer==0){
    for( int i3=0; i3<nh3; ++i3 ){
      if( UsedFlag3[i3]!=0 ) continue;
      DCHit *hit3=hc3[i3];
      int multi3 = hit3->GetDriftLengthSize();
      for (int m3=0; m3<multi3; m3++) {
	if( !(hit3->IsWithinTotRange(m3)) ) continue;
	if( !(hit3->IsWithinDtRange(m3)) ) continue;
	DCCluster *cluster = new DCCluster();
	cluster->SetHit(hit3,m3);
	Cont.push_back( cluster );
      }
    }  
  }
  for( int i=0;i<Cont.size();++i )  Cont[i]->Calc();
  m_DCCC.push_back(Cont); 
  //  std::cout<<slayer<<"  "<<Cont.size()<<std::endl;
  return true;
}
//______________________________________________________________________________
bool
CDCAnalyzer::FindCircleTrackCandidates(const bool &LINE)
{
#if DEBUG
  std::cout<<__PRETTY_FUNCTION__<<std::endl;
#endif
  int n1=GetNClusters(0);
  int n2=GetNClusters(3);
  int n3=GetNClusters(6);
  if(n1<1||n2<1||n3<1) return false;
  if(n1>100||n2>100||n3>100) return false;
  if(n1*n2*n3>1e4){
#if DEBUG
    std::cout<<"too many clusters in axial layers: "<<n1<<"  "<<n2<<"  "<<n3<<std::endl;
#endif
    return false;
  }
  int itr=0;
  for(int i1=0;i1<n1;i1++)
    {
      TVector3 pos1=(m_DCCC[0][i1])->pos();
      for(int i3=0;i3<n3; i3++)
	{
	  TVector3 pos3=m_DCCC[6][i3]->pos();
	  double angle13=pos1.Angle(pos3);
	  if(angle13>90*TMath::DegToRad() ) continue;
	  //#############mid axial cluster################
	  for(int i2=0;i2<n2; i2++){
	    TVector3 pos2=m_DCCC[3][i2]->pos();
	    if(LINE){
	      double tmpangle=(pos2-pos1).Angle(pos3-pos1);
	      if(tmpangle>10*TMath::DegToRad()) continue;
	    }else{
	      double angle23=pos2.Angle(pos3);
	      double angle12=pos1.Angle(pos2);
	      if(angle13<angle23) continue;
	      if(angle13<angle12) continue;
	    }
	    CDSTrack *tmptrack=new CDSTrack();
	    tmptrack->AddCluster(m_DCCC[0][i1]);
	    tmptrack->AddCluster(m_DCCC[3][i2]);
	    tmptrack->AddCluster(m_DCCC[6][i3]);
	    if(LINE){
	      double tmp[5]={pos1.X(),pos1.Y(),(pos3-pos1).Phi(),0,TMath::Pi()/2};
	      tmptrack->SetParameters(tmp);
	    }
	    tmptrack->SetTrackID(itr++);
	    m_CDCTC.push_back(tmptrack);
	  }//i2 
	}//i3
    }//i1	
#if DEBUG
  std::cout<<__PRETTY_FUNCTION__<<"  track container size: "<<m_CDCTC.size()<<std::endl;
#endif
  if(m_CDCTC.size()<1) return false; 
  return true;
}
//______________________________________________________________________________
void 
CDCAnalyzer::SelectSharedHit()
{ 
#if DEBUG
  std::cout<<__PRETTY_FUNCTION__<<std::endl;
#endif
  
  std::vector<int> goodtrack;
  int count=0;
  for(int itr=0;itr<(int)m_CDCTC.size();itr++){
    GetTrack(itr)->SetTrackID(itr);
    //    std::cout<<itr<<"  "<<GetTrack(itr)->nhits(0)<<std::endl;
  }
  while(1){
    if(count++>1000){
      std::cout<<"too many loops !! "<<__func__<<std::endl;
      return;
    }
    //###Delete track with too less Hit
    CDSTrackContainer::iterator ittr;
    for(ittr=m_CDCTC.end()-1;ittr!=m_CDCTC.begin()-1;--ittr)
      {
	CDSTrack *track=*ittr;	
	if(track->chi2()>9999 
	   ||
	   (track->naxialsuperlayers()<3)
	   ||
	   (track->FittingLevel()==2 && track->nhitlayers()<= 4)
	   ||
	   (track->FittingLevel()>2 && track->nhitlayers()<=10)
	   ||
	   (track->FittingLevel()>2 && track->nstereohits()<6)
	   ){
#if DEBUG
	  std::cout<<"Delete track: "
		   <<std::setw(10)<<track->chi2()
		   <<std::setw(5)<<track->naxialsuperlayers()
		   <<std::setw(5)<<track->FittingLevel()
		   <<std::setw(5)<<track->nhitlayers()
		   <<std::setw(5)<<track->nstereohits()
		   <<std::endl;
#endif
	  m_CDCTC.erase(ittr);
	  delete track;
	}
      }
    //    std::cout<<__func__<<"  "<<m_CDCTC.size()<<"  "<<goodtrack.size()<<std::endl;
    // select best chi2 track
    if(m_CDCTC.size()<=goodtrack.size()) break;
    double tmpchi=99999;
    int tmptr=0;
    for(int itr=0;itr<(int)m_CDCTC.size();itr++)
      {
	CDSTrack *track1=GetTrack(itr);
	if(vector_finder(goodtrack,track1->trackid())) continue;
	double chi1=track1->chi2();
	//std::cout<<itr<<"  "<<track1->trackid()<<"  "<<chi1<<"  "<<tmpchi<<"  "<<tmptr<<std::endl;
	if(chi1<tmpchi){
	  tmpchi=chi1;
	  tmptr=itr;
	}
      }
    CDSTrack *track1=GetTrack(tmptr);  
    // track1->Print();
    //    std::cout<<track1->chi2()<<std::endl;    
    goodtrack.push_back(track1->trackid());
    //      if( !(track1->FittingLevel()==2 || track1->FittingLevel()==4) ) continue;
    for(int layer=0;layer<15;layer++)
      { 
	for(int ih1=0;ih1<track1->nhits(layer);ih1++)
	  {	   
	    int wire1=track1->wire(layer,ih1);
	    for(int itr2=0;itr2<(int)m_CDCTC.size();itr2++)
	      { 
		CDSTrack *track2=GetTrack(itr2);
		if(vector_finder(goodtrack,track2->trackid())) continue;
		for(int ih2=0;ih2<track2->nhits(layer);ih2++)
		  {
		    int wire2=track2->wire(layer,ih2);
		    if(wire1==wire2)
		      track2->DeleteHit(layer,ih2);
		  }
	      }//track2
	  }//Hit1
      }//layer    
  }
  //  if(goodtrack.size()>10) exit(0);
}
//SelectSharedHit
//______________________________________________________________________________
bool 
CDCAnalyzer::FindStereoHits()
{   
#if DEBUG
  std::cout<<__PRETTY_FUNCTION__<<std::endl;
#endif
  
  //########Hit streo wire###########
  int stereo[4]={1,2,4,5};
  int ntr=m_CDCTC.size();
  for(int itr=0;itr<ntr;itr++)
    {
#if DEBUG
      std::cout<<"------------ "<<itr<<" / "<<ntr<<std::endl;
#endif
      std::vector<int> tmp[4];
      CDSTrack *track=GetTrack(itr);
      if(track->FittingLevel()!=2){
#if DEBUG
	std::cout<<"skip track because of bad fitting level"<<std::endl;
#endif
	continue;
      }
      int tmpn=0;
      int npat=1;
      for( int is=0; is<4; is++ )
	{
	  int slayer=stereo[is];
	  int ncl=m_DCCC[slayer].size();
#if DEBUG
	  std::cout<<is<<"  "<<ncl<<std::endl;
#endif
	  for( int k=0; k<ncl; k++ )
	    {	      
	      TVector3 pos=m_DCCC[slayer][k]->pos();
	      TVector3 dir=m_DCCC[slayer][k]->dir();
	      TVector3 posc=pos;//+dir*0.5;
	      double dist=999;
	      TVector3 postra=track->GetPositionatR(posc.Perp());
	      dist=(postra-posc).Mag();///dir.Mag2();
	      double wiredis = 100*mm;
	      wiredis=dir.Perp()*1.2;
	      bool tmpflag=( dist<wiredis );
#if DEBUG
	      std::cout<<"nhits  "<<m_DCCC[slayer][k]->nhit()<<std::endl;
	      std::cout<<"Pos:   "; pos.Print();
	      std::cout<<"Dir:   "; dir.Print();
	      std::cout<<"Wire:  "; posc.Print();
	      std::cout<<"Track:"; postra.Print();
	      std::cout<<"dist  "<<is<<"  "<<k<<"  "<<dist<<"  "<<wiredis<<"  "<<tmpflag<<std::endl;
#endif
	      if(tmpflag) tmp[is].push_back(k);
	    }      
	  if(tmp[is].size()>0) tmpn++;	  
	  npat*=tmp[is].size();
#if DEBUG
	  std::cout<<tmp[is].size()<<"  "<<npat<<std::endl;
#endif
	}//is

#if DEBUG
      std::cout<<"itr/ntr, tmpn, npat "<<itr<<"/"<<ntr<<"  "<<tmpn<<"  "<<npat<<std::endl;
#endif
      if(npat<1||npat>100){
#if DEBUG
	std::cout<<"skip track because of too many or less stereo combinations: "
		 <<npat<<std::endl;
#endif
	//exit(0);
	continue;
      }
      for(int is0=0;is0<(int)tmp[0].size();is0++){
	for(int is1=0;is1<(int)tmp[1].size();is1++){
	  for(int is2=0;is2<(int)tmp[2].size();is2++){
	    for(int is3=0;is3<(int)tmp[3].size();is3++){
	      CDSTrack* tmptrack=new CDSTrack(track);      
	      tmptrack->AddCluster(GetCluster(stereo[0],tmp[0][is0]) );
	      tmptrack->AddCluster(GetCluster(stereo[1],tmp[1][is1]) );
	      tmptrack->AddCluster(GetCluster(stereo[2],tmp[2][is2]) );
	      tmptrack->AddCluster(GetCluster(stereo[3],tmp[3][is3]) );
	      m_CDCTC.push_back(tmptrack);
	    }
	  }
	}
      }     
    }//track
  CDSTrackContainer::iterator ittr;
  for(ittr=m_CDCTC.end()-1;ittr!=m_CDCTC.begin()-1;--ittr)
    {
      CDSTrack *track=*ittr;	
      if(track->nstereohits()<6){
#if DEBUG
	std::cout<<__PRETTY_FUNCTION__<<"Delete track because of few stereo hits "<<track->nstereohits()<<std::endl;
#endif
	m_CDCTC.erase(ittr);
	delete track;
	continue;
      }
      //      track->Print();
    }  
#if DEBUG
  std::cout<<__PRETTY_FUNCTION__<<"  track container size: "<<m_CDCTC.size()<<std::endl;
#endif
  return true;  
}
//finding stereo hit
//____________________________________________________________________
bool
CDCAnalyzer::HelixFitting()
{
#if DEBUG
  std::cout<<__PRETTY_FUNCTION__<<std::endl;
#endif
  for(int i=0;i<(int)m_CDCTC.size();i++)
    {
      CDSTrack *track=GetTrack(i);            
      if(!track->FirstHelixFitting() ) 
	{
#if DEBUG
	  std::cout<<"failed in FirstHelixFitting()"<<std::endl;
#endif
	  track->SetGoodFlag(false);  continue;
	}
      if( track->chi2()>2*selected_chi2 ) 
	{
#if DEBUG
  std::cout<<"too large chi2"<<"  "<<track->chi2()<<std::endl;
#endif
	  track->SetGoodFlag(false);  continue;
	}
    }
  SelectSharedHit();
  for(int i=0;i<(int)m_CDCTC.size();i++)
    {
      CDSTrack *track=GetTrack(i);            
      track->SetTrackID(i);
      //      track->Print();
      if(!track->HelixFitting() ) continue;	      
      track->Calc();

    } 
  //  SelectSharedHit(cdsMan);
#if DEBUG
  std::cout<<__PRETTY_FUNCTION__<<"  track container size: "<<m_CDCTC.size()<<std::endl;
#endif
  return true;
}

bool
CDCAnalyzer::LineFitting()
{
#if DEBUG
  std::cout<<__PRETTY_FUNCTION__<<std::endl;
#endif
  for(int i=0;i<(int)m_CDCTC.size();i++)
    {
      CDSTrack *track=GetTrack(i);          
      //      track->Print();  
      if(!track->FirstLineFitting() ) 
	{
#if DEBUG
	  std::cout<<"failed in FirstLineFitting()"<<std::endl;
#endif
	  track->SetGoodFlag(false);  continue;
	}
      if( 0 && track->chi2()>3*selected_chi2 ) 
	{
#if DEBUG
	  std::cout<<"too large chi2"<<"  "<<track->chi2()<<std::endl;
#endif
	  track->SetGoodFlag(false);  continue;
	}
    }
  SelectSharedHit(); 
  for(int i=0;i<(int)m_CDCTC.size();i++)
    {
      CDSTrack *track=GetTrack(i);            
      track->SetTrackID(i);
      if(!track->LineFitting() ) continue;	      
      //      track->Print();
      //      track->Calc();
    }
  //  SelectSharedHit(cdsMan);
#if DEBUG
  std::cout<<__PRETTY_FUNCTION__<<"  track container size: "<<m_CDCTC.size()<<std::endl;
#endif
  return true;
}
//______________________________________________________________________________
bool
CDCAnalyzer::CircleFitting(){
#if DEBUG
  std::cout<<__PRETTY_FUNCTION__<<"  "<<GetNTracks()<<std::endl;
#endif
  for(int i=0;i<(int)m_CDCTC.size();i++)
    {
      CDSTrack *track=GetTrack(i);
      if(!track->CircleFitting()) 
	{
	  track->SetGoodFlag(false);
	}
    }
  SelectSharedHit();     
#if DEBUG
  std::cout<<__PRETTY_FUNCTION__<<"  track container size: "<<m_CDCTC.size()<<std::endl;
#endif
  return true;
}
//______________________________________________________________________________
bool
CDCAnalyzer::LineAxialFitting(){
#if DEBUG
  std::cout<<__PRETTY_FUNCTION__<<"  "<<GetNTracks()<<std::endl;
#endif
  for(int i=0;i<(int)m_CDCTC.size();i++)
    {
      CDSTrack *track=GetTrack(i);
      track->SetLineFlag();
      if(!track->LineAxialFitting()) 
	{
	  track->SetGoodFlag(false);
	}
    }
  SelectSharedHit();     
#if DEBUG
  std::cout<<__PRETTY_FUNCTION__<<"  track container size: "<<m_CDCTC.size()<<std::endl;
#endif
  return true;
}

//______________________________________________________________________________
void
CDCAnalyzer::SetCDCTracksFromTree( CDCTrackContainer* cont )
{
  ClearDCTracks();
  for(int it=0;it<cont->ntracks();it++){
    CDCTrackHits hits=cont->get_track(it);
    CDSTrack* tr=new CDSTrack();
    for(int ih=0;ih<hits.nhits();ih++){
      int layer=hits.get_hit(ih).layer();
      int wire=hits.get_hit(ih).wire();
      int nth=hits.get_hit(ih).nth();
      const DCHit* hit=GetDCHit(layer,wire);
      if(hit) tr->AddHit(hit,nth);
    }
    double p[5];
    hits.get_param(p);
    tr->SetParameters(p);
    if(!tr->HelixFitting() ) continue;	      
    tr->Calc();    
    tr->SetTrackID(GetNTracks());
    m_CDCTC.push_back(tr);    
  }
}
//______________________________________________________________________________
bool
CDCAnalyzer::DecodeRawHits( RawData *rawData, double retiming )
{
  ClearDCHits();
  DecodeRawHits( rawData, k_CDC, DetIdCDC, retiming );
  return true;
}

//______________________________________________________________________________
void
CDCAnalyzer::ClearDCHits( const int &detid )
{
  switch(detid){
  case DetIdCDC:
    del::ClearContainerAll( m_CDCHC );
    break;
  }
}
//______________________________________________________________________________
void
CDCAnalyzer::ClearDCHits( void)
{
  del::ClearContainerAll( m_CDCHC );
}

