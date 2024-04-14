/**
 *  file: UserSkeleton.cc
 *  date: 2017.04.10
 *
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include "TString.h"

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "RootHelper.hh"
#include "MTDCAnalyzer.hh"
#include "MTDCRawHit.hh"
#include "DCAnalyzer.hh"
#include "DCHit.hh"

#include "HodoAnalyzer.hh"
#include "HodoRawHit.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoCluster.hh"
#include "UserParamMan.hh"
#include "XTMapMan.hh"
#include "DCTdcCalibMan.hh"
#include "BLDCWireMapMan.hh"
#include "RawData.hh"
#include "VEvent.hh"
#include "HistTools.hh"
#include "UnpackerManager.hh"
#include "DeleteUtility.hh"

#include "setup.hh"

#define DEBUG 0
namespace
{
  using namespace e73_2024;
  using namespace root;
  const std::string& classname("EventHodo");
  using namespace hddaq::unpacker;
  const UnpackerManager& gUnpacker = GUnpacker::get_instance();
  const UserParamMan&   gUser = UserParamMan::GetInstance();
  int run_number;
  int event_number;
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
class EventBeam : public VEvent
{
private:
  RawData      *rawData;
  DCAnalyzer   *DCAna;
  HodoAnalyzer *hodoAna;
  MTDCAnalyzer *MTDCAna;

public:
  EventBeam( void );
  ~EventBeam( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventBeam::EventBeam( void )
  : VEvent(),
    rawData(0),
    DCAna(  new DCAnalyzer ),
    hodoAna(new HodoAnalyzer),
    MTDCAna(new MTDCAnalyzer )
{
}

//______________________________________________________________________________
EventBeam::~EventBeam( void )
{
  if( hodoAna ) delete hodoAna;
  if( DCAna )   delete DCAna;
  if( MTDCAna ) delete MTDCAna;
  if( rawData ) delete rawData;
}

//______________________________________________________________________________
namespace root
{
  TH1   *h[MaxHist];

  double tdcbins[3]={5000,0,2e6};
  double adcbins[3]={4096,-0.5,4095.5};
  double mtdcbins[3]={2000,0,2000};
  double debins[3]={5000,-2,48};
  double debins2[3]={1000,-100,900};
  double deudbins[6]={100,-2,18,100,-2,18};
}

//______________________________________________________________________________
bool
EventBeam::ProcessingBegin( void )
{
  InitializeEvent();
  event_number=gUnpacker.get_event_number();
  return true;
}

//______________________________________________________________________________
bool
EventBeam::ProcessingNormal( void )
{
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  std::cout << "# event number: "<<event_number<<std::endl;
#endif

  rawData = new RawData;
  rawData->DecodeHits();
  hodoAna->DecodeRawHits( rawData );
  MTDCAna->DecodeRawHits( rawData );
  bool BEAM  = MTDCAna->flag(kBeam);
  bool KAON2 = MTDCAna->flag(kKaon2);
  bool KAON3 = MTDCAna->flag(kKaon3);

  bool TOFK=false, TOFPi=false, TOFP=false, TOFD=false;
  bool ACHIT=false;

  // Time0
  double time0=-9999;
  double ctime0=-9999;
  {
    int cid=DetIdT0;
    int nh = hodoAna->GetNHits(cid);
    for( int i=0; i<nh; ++i ){
      Hodo2Hit *hit = hodoAna->GetHit(cid,i);
      if(!hit) continue;
      int nind=hit->GetIndex();
      for(int it=0;it<nind;it++){
	double mt  = hit->MeanTime(it);
	if(gUser.Check("Time0",mt)){
	  time0=mt;
	  ctime0=hit->CMeanTime(it);
	}
      }
    }
  }
  // K/pi by BHD MeanTime
  {
    int cid=DetIdBHT;
    int nh = hodoAna->GetNHits(cid);
    for( int i=0; i<nh; ++i ){
      Hodo2Hit *hit = hodoAna->GetHit(cid,i);
      if(!hit) continue;
      int nind=hit->GetIndex();
      for(int it=0;it<nind;it++){
	double tof  = hit->MeanTime(it)-time0;
	if(gUser.Check("TOFK", tof)) TOFK=true;
	if(gUser.Check("TOFPi",tof)) TOFPi=true;
	if(gUser.Check("TOFP", tof)) TOFP=true;
	if(gUser.Check("TOFD", tof)) TOFD=true;
      }
    }
  }
  // AC
  {
    int cid=DetIdAC;
    const HodoRHitContainer &cont = rawData->GetHodoRawHC(cid);
    int nh = cont.size();
    TString tmpname="AC";
    for( int i=0; i<nh; ++i ){
      HodoRawHit *raw = cont[i];
      if(!raw) continue;
      int ntu=raw->GetSizeTdcUp();
      for(int it=0;it<ntu;it++){
	double tu  = raw->GetTdcUp(it);
	if(gUser.Check("ACTDC",tu)) ACHIT=true;
	hist::H1(tmpname+"_TDC",tu,tdcbins);
      }
    }
  }

  std::vector<TString> trig_add;
  trig_add.push_back("");
  if(BEAM) trig_add.push_back("ifB");
  if(BEAM&&ACHIT&&TOFPi) trig_add.push_back("ifPi");
  if(TOFD&&!ACHIT)  trig_add.push_back("ifD");
  if(TOFP&&!ACHIT)  trig_add.push_back("ifP");
  //  if(KAON2) trig_add.push_back("ifK2");
  //  if(KAON3) trig_add.push_back("ifK");
  int ntrig=trig_add.size();
  // hodoscopes
  //  for(int ihodo=0;ihodo<nhodo;++ihodo){
  {
    int kHodo=kBHT;
    int cid=hodoid[kHodo];
    TString name=hodoname[kHodo];
    double mulbins[3]={nsegs[kHodo]+1,-0.5,nsegs[kHodo]+0.5};
    double patbins[3]={nsegs[kHodo],-0.5,nsegs[kHodo]-0.5};
    double mulbins2[3]={10,-0.5,9.5};
    // Decoded hit
    int mul=0,mulgate=0;
    int nh = hodoAna->GetNHits(cid);
    std::vector<double> mt_list;
    std::vector<int> seg_list;
    for( int i=0; i<nh; ++i ){
      Hodo2Hit *hit = hodoAna->GetHit(cid,i);
      if(!hit) continue;
      int seg = hit->SegmentId();
      TString segstr=Form("_seg%d",seg);
      int nind= hit->GetIndex();
      hist::H1(name+"_Mul"+segstr,nind,mulbins2);
      hist::H2(name+"_Mul_pat",seg,nind,patbins,mulbins2);
      if(nind>0){
	hist::H1(name+"_Pat",seg,patbins);
	mul++;
      }
      int idx=-1;
      for(int ii=0;ii<nind;ii++){
	double mt  = hit->MeanTime(ii);
	double cmt  =hit->CMeanTime(ii);
	double tof = mt - time0;
	double ctof = cmt - ctime0;
	hist::H1(name+"_MeanTime"        ,mt ,4000,-200,200, trig_add);
	hist::H1(name+"_MeanTime" +segstr,mt ,4000,-200,200, trig_add);
	hist::H1(name+"_CMeanTime"       ,cmt,4000,-200,200, trig_add);
	hist::H1(name+"_CMeanTime"+segstr,cmt,4000,-200,200, trig_add);
	if(time0>-9000){
	  hist::H1(name+"_TOF"        ,tof ,2000,-100,100, trig_add);
	  hist::H1(name+"_TOF" +segstr,tof ,2000,-100,100, trig_add);
	  hist::H1(name+"_cTOF"       ,ctof,2000,-100,100, trig_add);
	  hist::H1(name+"_cTOF"+segstr,ctof,2000,-100,100, trig_add);
	  double de = hit->DeltaCE(ii);
	  double demin=0,demax=4e4;
	  hist::H2(name+"_TOF_dE",tof , de, 200,-50,50, 100, demin, demax, trig_add);
	}
	if(TMath::Abs(mt)>20) continue;
	mt_list.push_back(mt);
	seg_list.push_back(seg);
	idx=ii;
	mulgate++;
      } // nindex
      if(idx<0) continue;
      double demin=0,demax=1e5,demax2=2e4;
      double au = hit->GetAUp(idx);
      double ad = hit->GetADown(idx);
      double de = hit->DeltaCE();
      double deu = hit->GetCAUp();
      double ded = hit->GetCADown();
      double tmpdebins[3]={1000,demin,demax};
      double tmpdeudbins[6]={100,demin,demax2,100,demin,demax2};
      hist::H1(name+"_dE" +segstr,de ,tmpdebins,trig_add);
      hist::H2(name+"_dEud"+segstr,deu,ded,tmpdeudbins,trig_add);
      if(BEAM&&ACHIT){
	hist::H1(name+"_ADCuwB"+segstr,au,adcbins);
	hist::H1(name+"_ADCdwB"+segstr,ad,adcbins);
	hist::H1(name+"_ADCuwB",au,adcbins);
	hist::H1(name+"_ADCdwB",ad,adcbins);
	hist::H1(name+"_dEuwB"+segstr,deu,tmpdebins);
	hist::H1(name+"_dEdwB"+segstr,ded,tmpdebins);
	hist::H1(name+"_dEuwB",deu,tmpdebins);
	hist::H1(name+"_dEdwB",ded,tmpdebins);
      }
    }//for(ihit)

    for(int i=0;i<mt_list.size();i++){
      for(int j=i+1;j<mt_list.size();j++){
	hist::H1(name+"_timediff",mt_list.at(i)-mt_list.at(j),100,-50,50);
	hist::H2(name+"_timecorr",mt_list.at(i),mt_list.at(j),100,-50,50,100,-50,50);
	hist::H1(name+"_segdiff",seg_list.at(i)-seg_list.at(j),100,-50.5,49.5);
	hist::H2(name+"_segcorr",seg_list.at(i),seg_list.at(j),63,-0.5,62.5,63,-0.5,62.5);
	hist::H2(name+"_timediff_segdiff",mt_list.at(i)-mt_list.at(j),seg_list.at(i)-seg_list.at(j),100,-50,50,100,-50.5,49.5);
      }
    }
    hist::H1(name+"_Mul",mul ,mulbins);
    hist::H1(name+"_Mulgate",mulgate,mulbins);
    // decoded hit done
#if 1
    std::vector<HodoCluster*> clusterContainer;
    std::vector<std::vector<int>> UsedFlag(nh,std::vector<int>(0));
    for( int i=0; i<nh; ++i ){
      Hodo2Hit *hit1 = hodoAna->GetHit(cid,i);
      if(!hit1) continue;
      int seg1 = hit1->SegmentId();
      int nind1= hit1->GetIndex();
      for(int i1=0;i1<nind1;i1++){
	if(std::find(UsedFlag.at(i).begin(), UsedFlag.at(i).end(), i1) != UsedFlag.at(i).end()) continue;
	double mt1  = hit1->MeanTime(i1);
	HodoCluster* clu=new HodoCluster();
	clu->AddHit(hit1,i1);
	UsedFlag.at(i).push_back(i1);
	double mt=clu->CMeanTime();
	int j=1;
	while(i+j<nh){
	  Hodo2Hit *hit2 = hodoAna->GetHit(cid,i+j);
	  if(hit2){
	    int seg2 = hit2->SegmentId();
	    //	    if( seg2!=(seg1+j) ) break;
	    if( seg2>(seg1+j+1) ) break;
	    int nind2= hit2->GetIndex();
	    for(int i2=0;i2<nind2;i2++){
	      if(std::find(UsedFlag.at(i+j).begin(), UsedFlag.at(i+j).end(), i2) != UsedFlag.at(i+j).end()) continue;	      double mt2  = hit2->MeanTime(i2);
	      if(TMath::Abs(mt-mt2)<15){
		clu->AddHit(hit2,i2);
		UsedFlag.at(i+j).push_back(i2);
		break;
	      }
	    }//i2
	  } //hit2
	  j++;
	}//while
	clusterContainer.push_back(clu);
      }//i1
    }//ihit
    int mulcluster=clusterContainer.size();
    hist::H1(name+"_Cluster_Mul",mulcluster,mulbins);
    int mulcluster_gate=0;
    bool PRINT=false;
    for(int i=0; i<mulcluster; i++){
      HodoCluster* clu=clusterContainer.at(i);
      Hodo2Hit* hit=clu->Get1stHit();
      int       nth=clu->Get1stNthHit();
      double de = hit->DeltaE(nth);
      double cde= hit->DeltaCE(nth);
      double mt = hit->MeanTime(nth);
      double cmt= hit->CMeanTime(nth);
      int seg   = hit->SegmentId();
      double tof = mt - time0;
      double ctof= cmt - ctime0;
      TString segstr=Form("_seg%d",seg);
      hist::H1(name+"_Cluster_Pat",seg,patbins);
      hist::H1(name+"_Cluster_MeanTime"        ,mt ,4000,-200,200, trig_add);
      hist::H1(name+"_Cluster_MeanTime" +segstr,mt ,4000,-200,200, trig_add);
      hist::H1(name+"_Cluster_CMeanTime"       ,cmt,4000,-200,200, trig_add);
      hist::H1(name+"_Cluster_CMeanTime"+segstr,cmt,4000,-200,200, trig_add);
      if(mt>-40&&mt<10){
	mulcluster_gate++;
	hist::H1(name+"_Cluster_Pat_gate",seg,patbins);
      }
      if(mt>-25&&mt<-15) PRINT=true;
      if(time0>-9000){
	hist::H1(name+"_Cluster_TOF"        ,tof ,2000,-100,100, trig_add);
	hist::H1(name+"_Cluster_TOF" +segstr,tof ,2000,-100,100, trig_add);
	hist::H1(name+"_Cluster_cTOF"       ,ctof,2000,-100,100, trig_add);
	hist::H1(name+"_Cluster_cTOF"+segstr,ctof,2000,-100,100, trig_add);
	hist::H2(name+"_Cluster_TOF_dE" ,tof , de , 200,-50,50, 100, 0, 5e4, trig_add);
	hist::H2(name+"_Cluster_TOF_cdE",tof , cde, 200,-50,50, 100, 0,   5, trig_add);
      }
    }
    hist::H1(name+"_Cluster_Mul_gate",mulcluster_gate,mulbins);
    if(false){
      std::cout<<std::endl<<std::endl;;
      std::cout<<"### Event number: "<<event_number<<std::endl;
      for(int i=0; i<mulcluster; i++){
	std::cout<<"================= "<<i<<" / "<<clusterContainer.size()<<std::endl;
	clusterContainer.at(i)->Print();
      }
    }
    del::ClearContainer( clusterContainer );
#endif
  }

#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  return true;
}

//______________________________________________________________________________
bool
EventBeam::ProcessingEnd( void )
{
  return true;
}

//______________________________________________________________________________
void
EventBeam::InitializeEvent( void )
{
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventBeam;
}
//______________________________________________________________________________
bool
ConfMan::InitializeHistograms( void )
{
  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<HodoParamMan>("HDPRM","CDSPRM") ) &&
    ( InitializeParameter<HodoPHCMan>  ("HDPHC","CDSPHC") ) &&
    ( InitializeParameter<UserParamMan>("USER") ) ;
}

//______________________________________________________________________________
bool
ConfMan::BeginRunProcess()
{
  run_number = get_run_number();
  return true;
}
//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
