// -*- C++ -*-

#include "VEvent.hh"

#include <fstream>
#include <iostream>
#include <sstream>

#include <TString.h>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RootHelper.hh"
#include "DCAnalyzer.hh"
#include "DCCluster.hh"
#include "DCHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoRawHit.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "MTDCAnalyzer.hh"
#include "MTDCRawHit.hh"
#include "UserParamMan.hh"
#include "GeomMapMan.hh"
#include "XTMapMan.hh"
#include "BLDCWireMapMan.hh"
#include "DCTdcCalibMan.hh"
#include "DCTimeCorrMan.hh"
#include "RawData.hh"
#include "UnpackerManager.hh"
#include "TransferMatrixMan.hh"
#include "BeamSpectrometer.hh"
#include "LocalTrack.hh"
#include "HistTools.hh"

#include "setup.hh"

#define BEAMONLY 0
#define TRACKING 1
#define GLOBAL 0
#define CLUSTER 0
#define SOURCE 0
#define DEBUG 0

namespace
{
using namespace e73_2024;
using namespace root;
using namespace hddaq::unpacker;
const auto& gUser = UserParamMan::GetInstance();
const auto& gUnpacker = GUnpacker::get_instance();

const double mm=0.1;
const double cm=10*mm;

int run_number;
int event_number;
std::ofstream ofsk;
std::ofstream ofspi;
}

//_____________________________________________________________________________
namespace root
{
TH1   *h[MaxHist];

const int nhodo=4;
int khodo[nhodo]={kT1,kBHT,kT0,kDEF};

const int nchm=9;
const int nchm2=6;
int kchm[nchm]={kBLC1a,kBLC1b,kBLC2a,kBLC2b,kBPC1,kBPC2,kBLC1,kBLC2,kBPC0};

double tdcbins[3]={5000,0,2e6};
double adcbins[3]={4096,-0.5,4095.5};
double mtdcbins[3]={2000,0,2000};
double dtbins[3]={2500,-500,2000};
double dttotbins[6]={250,-250,1000,250,0,500};

double posbins1[6]={2000,-20*cm,20*cm};
double posbins2[6]={150,-15*cm,15*cm,150,-15*cm,15*cm};
double residbins[3]={200,-5*mm,5*mm};

double evdtbins[6]={200,0,5e6,100,-100,400};
double evresidbins[6]={200,0,5e6,100,-2.5*mm,2.5*mm};
}

//_____________________________________________________________________________
void
InitializeEvent()
{
}

//_____________________________________________________________________________
Bool_t
ProcessBegin()
{
  InitializeEvent();
  event_number = gUnpacker.get_event_number();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessNormal()
{
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif

  RawData rawData;
  rawData.DecodeHits();
  MTDCAnalyzer MTDCAna(rawData);
  MTDCAna.DecodeRawHits();
  if(MTDCAna.flag(kCosmic)) return true;
  HodoAnalyzer hodoAna(rawData);
  hodoAna.DecodeRawHits();
  Bool_t BEAM =MTDCAna.flag(kBeam);
#if BEAMONLY
  if(!BEAM) return true;
#endif
  Bool_t KAON2=MTDCAna.flag(kKaon2);
  Bool_t KAON3=MTDCAna.flag(kKaon3);

  //  if(!BEAM&&!KAON3) return true;
  Bool_t TOFK=false, TOFPi=false;
  Bool_t ACHIT=false;

  // --
  {
    int cid=DetIdAC;
    const HodoRHitContainer &cont = rawData.GetHodoRawHC(cid);
    int nh = cont.size();
    TString tmpname="AC";
    for( int i=0; i<nh; ++i ){
      HodoRawHit *raw = cont[i];
      if(!raw) continue;
      //      int seg = raw->SegmentId();
      int ntu=raw->GetSizeTdcUp();
      for(int it=0;it<ntu;it++){
	double tu  = raw->GetTdcUp(it);
	if(gUser.IsInRange("ACTDC",tu)) ACHIT=true;
	hist::H1(tmpname+"_TDC",tu,tdcbins);
      }
    }
  }
  //  if(!ACHIT) return true;
  // -- hodoscopes
  int    hodoseg[kNumHodo];
  double hodode[kNumHodo];
  double hodotime[kNumHodo];
  double hodomul[kNumHodo];
  double time0=-9999;
  Bool_t GoodBeam=true;
  for(int ihodo=0;ihodo<nhodo;++ihodo){
    int kHodo=khodo[ihodo];
    hodoseg[kHodo]=-1;
    hodomul[kHodo]=0;
    int cid=DetIdHodo[kHodo];
    TString tmpname=NameHodo[kHodo];
    int mul=0,mulgate=0;
    double mulbins[3]={NumOfSegHodo[kHodo]+1,-0.5,NumOfSegHodo[kHodo]+0.5};
    double mulbins2[3]={10,-0.5,9.5};
    double patbins[3]={NumOfSegHodo[kHodo],-0.5,NumOfSegHodo[kHodo]-0.5};
    int nh = hodoAna.GetNHits(cid);
    for( int i=0; i<nh; ++i ){
      Hodo2Hit *hit = hodoAna.GetHit(cid,i);
      if(!hit) continue;
      int seg = hit->SegmentId();
      int nind=hit->GetIndex();
      TString segstr=Form("_seg%d",seg);
      hist::H1(tmpname+"_Mul"+segstr,nind,mulbins2);
      hist::H2(tmpname+"_Mul_pat",seg,nind,patbins,mulbins2);
      if(nind>0){ hist::H1(tmpname+"_Pat",seg,patbins);	mul++; }
      for(int it=0;it<nind;it++){
	double mt  = hit->MeanTime(it);
	if(cid==DetIdT1){
	  if(TMath::Abs(mt)<20){
	    time0=mt;
	  }
	}
	double tof =mt-time0;
	if(cid==DetIdBHD){
	  if(gUser.IsInRange("TOFK",tof)) TOFK=true;
	  if(gUser.IsInRange("TOFPi",tof)) TOFPi=true;
	}
	hist::H1(tmpname+"_MeanTime",mt,4000,-200,200);
	hist::H1(tmpname+"_MeanTime"+segstr,mt,4000,-200,200);
	if(time0>-9000){
	  hist::H1(tmpname+"_TOF",tof,2000,-100,100);
	  if(KAON2) hist::H1(tmpname+"_TOFifK2",tof,2000,-100,100);
	  if(KAON3) hist::H1(tmpname+"_TOFifK3",tof,2000,-100,100);
	  if(BEAM)  hist::H1(tmpname+"_TOFifB",tof,2000,-100,100);
	  if(ACHIT) hist::H1(tmpname+"_TOFifPi",tof,2000,-100,100);
	  hist::H1(tmpname+"_TOF"+segstr,tof,2000,-100,100);
	}
	if(TMath::Abs(mt)<50){
	  mulgate++;
	  hodoseg[kHodo]=seg;
	  hodotime[kHodo]=mt;
	  hodomul[kHodo]++;
	  hodode[kHodo]=hit->DeltaE();
	}
      }
    }//for(i)
    hist::H1(tmpname+"_Mul" ,mul ,mulbins);
    hist::H1(tmpname+"_Mulgate" ,mulgate,mulbins);
    if(cid!=DetIdBHD&&mulgate!=1) GoodBeam=false;
  }

  if(time0<-999) return true;

  // Chamber ------------------------------------------------------------
  DCAnalyzer DCAna(rawData);
  Bool_t Single[nchm2];
  //  DCAna.DecodeRawHits( rawData, hodotime[kT0], hodotime[kDEF] );
#if SOURCE
  DCAna.DecodeRawHits();

#else
  DCAna.DecodeRawHits(hodotime[kT0], hodotime[kT0]);
  if(GoodBeam&&(TOFK||TOFPi))
#endif
    {
      for(int ichm=0;ichm<nchm2;ichm++){
#ifndef SOURCE
	if(hodomul[kT0]!=1) continue;
	if(hodomul[kDEF]!=1) continue;
#endif
	//    DCAna.DecodeRawHits( rawData, hodotime[kT1],hodotime[kDEF] );
	int kChm=kchm[ichm];
	TString tmpname=NameDC[kChm];
	int    cid     =DetIdDC[kChm];
	double nl      =NumOfLayerDC[kChm];
	double nw      =NumOfWireDC[kChm];
	double mulbins[3] ={nw,-0.5,nw+0.5};
	double mulbins2[3]={10,-0.5,9.5};
	double patbins[3] ={nw,-0.5,nw-0.5};
	double lpatbins[3]={nl,-0.5,nl-0.5};
	Bool_t GOOD[8];
	Single[ichm]=true;
	int hw[8];
	double wpos[8];
	for( int layer=0; layer<(int)nl; ++layer ){
	  TString lstr=Form("_layer%d",layer);
	  GOOD[layer]=true;
	  const DCHitContainer &cont = DCAna.GetDCHC(cid, layer);
	  int mul=cont.size();
	  int mul_good=0;
	  hist::H1(tmpname+"_Mul"+lstr,mul,mulbins);
	  for(int ihit=0;ihit<mul;ihit++){
	    DCHit* hit=cont[ihit];
	    int mul_wire=hit->GetNHit(false,true);
	    if(mul_wire==1){
	      mul_good++;
	    }
	  }
	  if(mul_good!=1){
	    GOOD[layer]=false;
	    Single[ichm]=false;
	  }
	  hist::H1(tmpname+"_Mulgood"+lstr,mul_good,mulbins);
	}
	for( int layer=0; layer<(int)nl; ++layer ){
	  const DCHitContainer &cont = DCAna.GetDCHC(cid, layer);
	  int mul=cont.size();
	  for(int ihit=0;ihit<mul;ihit++){
	    DCHit* hit=cont[ihit];
	    int wire=hit->GetWire();
	    TString lstr=Form("_layer%d",layer);
	    TString lwstr=Form("_layer%d_wire%d",layer,wire);
	    hw[layer]=wire;
	    wpos[layer]=hit->GetWirePosition().X();
	    int mul_wire=hit->GetTdcSize();
	    int mul_wire3=hit->GetNHit(true,true);
	    hist::H1(tmpname+"_HitPat"+lstr	,hit->GetWire(),patbins);
	    hist::H1(tmpname+"_HitPatgood"+lstr	,hit->GetWire(),patbins);
	    hist::H2(tmpname+"_HitPat_2d"		,layer,hit->GetWire(),lpatbins,patbins);
	    hist::H1(tmpname+"_Mul"+lwstr		,mul_wire,mulbins);
	    hist::H1(tmpname+"_Mulgood"+lwstr	,mul_wire3,mulbins);
	    hist::H2(tmpname+"_Mul_2d"+lstr	,wire,mul_wire,patbins,mulbins2);
	    hist::H2(tmpname+"_Mulgood_2d"+lstr	,wire,mul_wire3,patbins,mulbins2);
	    if(Single[ichm]){
	      int id=hit->GetHitID(0,false,true);
	      if(id>=0){
		hist::H1(tmpname+"_Leading_single"+lstr	,hit->GetTdcVal(id)     ,mtdcbins);
		hist::H1(tmpname+"_Trailing_single"+lstr	,hit->GetTdcTrailing(id),mtdcbins);
		hist::H1(tmpname+"_TOT_single"+lstr	,hit->GetTOT(id)        ,mtdcbins);
		hist::H1(tmpname+"_dt_single"+lstr	,hit->GetDriftTime(id)  ,dtbins);
		hist::H1(tmpname+"_dt_single"+lwstr	,hit->GetDriftTime(id)  ,dtbins);
		//	    hist::H1(tmpname+"_dt_layer%d_wire%d",layer,hit->GetWire()),hit->GetDriftTime(),dtbins);
		hist::H2(tmpname+"_EvNum_dt_single"+lstr, event_number,hit->GetDriftTime(id),evdtbins);
	      }
	    }
	    if(mul_wire3==1&&GOOD[layer]){
	      int id=hit->GetHitID(0,true,true);
	      if(id>=0){
		hist::H1(tmpname+"_Leading_good"+lstr,hit->GetTdcVal(id),mtdcbins);
		hist::H1(tmpname+"_Trailing_good"+lstr,hit->GetTdcTrailing(id),mtdcbins);
		hist::H1(tmpname+"_TOT_good"+lstr,hit->GetTOT(id),mtdcbins);
		hist::H1(tmpname+"_dt_good"+lstr,hit->GetDriftTime(id),dtbins);
	      }
	    }
	    for(int i=0;i<hit->GetTdcSize();i++){
	      hist::H1(tmpname+"_Leading"+lstr,hit->GetTdcVal(i),mtdcbins);
	    }
	    for(int i=0;i<hit->GetTdcTrailingSize();i++)
	      hist::H1(tmpname+"_Trailing"+lstr,hit->GetTdcTrailing(i),mtdcbins);
	    for(int i=0;i<hit->GetTOTSize();i++){
	      hist::H1(tmpname+"_TOT"+lstr,hit->GetTOT(i),mtdcbins);
	      hist::H2(tmpname+"_dt_TOT"+lstr,hit->GetDriftTime(i),hit->GetTOT(i),dttotbins);
	      if(hit->IsWithinTotRange(i))
		hist::H1(tmpname+"_TOT_range"+lstr,hit->GetTOT(i),mtdcbins);
	    }
	    Bool_t tmpfirst=true;
	    for(int i=0;i<hit->GetDriftTimeSize();i++){
	      hist::H1(tmpname+"_dt"+lstr,hit->GetDriftTime(i),dtbins);
	      if(hit->IsWithinDtRange(i))
		hist::H1(tmpname+"_dt_range"+lstr,hit->GetDriftTime(i),dtbins);
	      if(hit->IsWithinTotRange(i)){
		hist::H1(tmpname+"_dt_tot"+lstr,hit->GetDriftTime(i),dtbins);
		if(tmpfirst){
		  hist::H1(tmpname+"_dt_totfirst"+lstr,hit->GetDriftTime(i),dtbins);
		  hist::H1(tmpname+"_dt_totfirst"+lwstr,hit->GetDriftTime(i),dtbins);
		  tmpfirst=false;
		}
	      }
	    }
	  }//ihit
	} //layer
	if(Single[ichm]){
	  for( int layer=0; layer<nl;layer+=2 ){
	    int layer2=layer+1;
	    hist::H2(tmpname+Form("_hitcorr_%d_%d",layer,layer2),hw[layer],hw[layer2],patbins,patbins);
	    hist::H2(tmpname+Form("_poscorr_%d_%d",layer,layer2),wpos[layer],wpos[layer2],posbins2);
	  }
	  for( int layer=0; layer<nl/2;layer++ ){
	    int layer2=layer+nl/2;
	    hist::H2(tmpname+Form("_hitcorr_%d_%d",layer,layer2),hw[layer],hw[layer2],patbins,patbins);
	    hist::H2(tmpname+Form("_poscorr_%d_%d",layer,layer2),wpos[layer],wpos[layer2],posbins2);
	  }
	}// single
      } // ichm
    }// chamber analysis

#if TRACKING
  DCAna.MakePairsAll();
  {
    for(int ichm=0;ichm<nchm2;ichm++){
      int kChm=kchm[ichm];
      TString tmpname=NameDC[kChm];
      int     cid    =DetIdDC[kChm];
      //      if(cid==DetIdBPC&&!GoodBeam) continue;
      //    DCAna.MakePairs(cid);
      int status=DCAna.TrackSearch(cid,true,0);
      int ntra=DCAna.GetNTracks(cid);
      hist::H1(tmpname+"_nTracks",  ntra,10,-0.5,9.5);
      hist::H1(tmpname+"_status" ,status,10,-0.5,9.5);
      if(KAON3)
	hist::H1(tmpname+"_nTracks_K3", ntra,10,-0.5,9.5);
      if(KAON3&&GoodBeam)
	hist::H1(tmpname+"_nTracks_goodK3", ntra,10,-0.5,9.5);
#if CLUSTER
      std::vector<TString> tmpadd;
      tmpadd.push_back("");
      if(Single[ichm])   tmpadd.push_back("_single");
      for( int xy=0;xy<2;xy++){
	int n1=DCAna.GetNClusterContainers(cid,xy);
	for(int i1=0;i1<n1;i1++){
	  int n2=DCAna.GetNClusters(cid,xy,i1);
	  hist::H1(tmpname+Form("_nClusters_%d",xy), n2,10,-0.5,9.5);
	  for(int i2=0;i2<n2;i2++){
	    DCCluster* cl=DCAna.GetCluster(cid,xy,i1,i2);
	    if(cl->nhit()==2){
	      TString clstr=Form("_cluster%d%d"   ,xy,i1);
	      for(int iadd=0;iadd<(int)tmpadd.size();iadd++){
		hist::H1(tmpname+"_time"         +tmpadd.at(iadd)+clstr, cl->time()   ,1000,-100,400);
		hist::H1(tmpname+"_timesub"      +tmpadd.at(iadd)+clstr, cl->timesub(),1000,-250,250);
		hist::H1(tmpname+"_ctime"        +tmpadd.at(iadd)+clstr, cl->ctime()  ,1000,-100,400);
		hist::H2(tmpname+"_timesub_time" +tmpadd.at(iadd)+clstr, cl->timesub(), cl->time() , 100,-200,200,100,-100,400);
		hist::H2(tmpname+"_timesub_ctime"+tmpadd.at(iadd)+clstr, cl->timesub(), cl->ctime(), 100,-200,200,100,-100,400);
	      }
	    }
	  }
	}
      }
#endif
      for(int itr=0;itr<ntra;itr++){
	LocalTrack* track=DCAna.GetTrack(cid,itr);
#if 0
	for(int i=0;i<track->nclustertimes();i++){
	  std::cout<<i<<"  "<<track->clustertime(i)<<std::endl;
	}
	std::cout<<"Time: "<<track->GetTrackTime()<<std::endl;
	std::cout<<"RMS:  "<<track->GetTrackTimeRMS()<<std::endl;
#endif
	hist::H1(tmpname+"_time",   track->GetTrackTime(),2000,-200,200);
	hist::H1(tmpname+"_timerms",track->GetTrackTimeRMS(),1000,0,200);
      }
      if(ntra==1){
	LocalTrack* track=DCAna.GetTrack(cid,0);
	hist::H1(tmpname+"_chi2all",track->chi2all(),1000,0,100);
	hist::H1(tmpname+"_chi2xz",track->chi2xz(),1000,0,100);
	hist::H1(tmpname+"_chi2yz",track->chi2yz(),1000,0,100);
	double x,y;
	track->XYLocalPosatZ(0,x,y);
	hist::H2(tmpname+"_XYLocal",x,y,posbins2);
	hist::H2(tmpname+"_AB",track->gdx(),track->gdy(),100,-0.1,0.1,100,-0.1,0.1);
	for(int xy=0;xy<2;xy++){
	  for(int i=0;i<track->nhit(xy);i++){
	    double resid=track->resid(xy,i);
	    int layer=track->layer(xy,i);
	    hist::H1(tmpname+"_resid"+Form("_layer%d", layer ), resid, residbins);
	    hist::H2(tmpname+"_EvNum_dE"+Form("_layer%d", layer ) ,event_number,resid,evresidbins);
	  }
	}
	switch(cid){
	case DetIdBLC1:
	case DetIdBLC1a:
	case DetIdBLC1b:
	  track->XYPosatZ(0,x,y);
	  hist::H2(tmpname+"_XY",x,y,posbins2);
	  break;
	case DetIdBLC2:
	case DetIdBLC2a:
	case DetIdBLC2b:
	  track->XYPosatZ(-130*cm,x,y);
	  hist::H2(tmpname+"_XY",x,y,posbins2);
	  if(hodoseg[kT0]>-1){
	    track->XYLocalPosatZ(20*cm,x,y);
	    hist::H2(tmpname+"_XY"+Form("_ifT0seg%d",hodoseg[kT0]),x,y,posbins2);
	  }
	  break;
	case DetIdBPC1:
	case DetIdBPC2:
	  track->XYPosatZ(0,x,y);
	  hist::H2(tmpname+"_XYatFF",x,y,posbins2);
	  if(hodoseg[kDEF]>-1){
	    track->XYLocalPosatZ(5*cm,x,y);
	    hist::H2(tmpname+"_XY"+Form("_ifDEFseg%d",hodoseg[kDEF]),x,y,posbins2);
	  }
	  break;
	}
      }
    }
  }

#if GLOBAL
  int ntra1=DCAna.GetNTracks(DetIdBLC1);
  int ntra2=DCAna.GetNTracks(DetIdBLC2);
  int nbpc=DCAna.GetNTracks(DetIdBPC1);
  if(nbpc==1&&ntra2==1){
    LocalTrack* trbpc=DCAna.GetTrack(DetIdBPC,0);
    LocalTrack* tr2=DCAna.GetTrack(DetIdBLC2,0);
    double z=-300*mm;
    double bpcx,bpcy;
    double blc2x,blc2y;
    trbpc->XYPosatZ(z,bpcx,bpcy);
    tr2->XYPosatZ(z,blc2x,blc2y);
    TVector3 dirbpc=trbpc->GetMomDir();
    TVector3 dirblc2=tr2->GetMomDir();
    TVector3 diffdir=dirblc2-dirbpc;
    double tmpbins[6]={200,-100*mm,100*mm,200,-100*mm,100*mm};
    hist::H2("CorrX_BLC2BPC",blc2x,bpcx,tmpbins);
    hist::H2("CorrY_BLC2BPC",blc2y,bpcy,tmpbins);
    hist::H2("X_diffY_BLC2BPC",bpcx,blc2y-bpcy,tmpbins);
    hist::H2("Y_diffX_BLC2BPC",bpcy,blc2x-bpcx,tmpbins);
    hist::H2("CorrAB_BLC2BPC",TMath::ATan(diffdir.X())*1000,TMath::ATan(diffdir.Y())*1000,tmpbins);
  }
  if(ntra1==1&&ntra2==1){
    LocalTrack* tr1=DCAna.GetTrack(DetIdBLC1,0);
    LocalTrack* tr2=DCAna.GetTrack(DetIdBLC2,0);
    if(tr1->chi2all()<5&&tr2->chi2all()<5){
      //     tr1->Print();
      //     tr2->Print();
      beam->TMinuitFit(tr1,tr2);
      double mom=beam->mom();
      double chi=beam->chisquare();
      hist::H1("D5Mom",mom,1000,0.5,1.0);
      hist::H1("D5chi2",chi,5000,0,1000);
      //  std::cout<<mom<<"  "<<chi<<std::endl;
      double blc1par[6];
      double parblc2[6];
      blc1par[0]=tr1->GetPosatZ(0).X(); //cm
      blc1par[1]=TMath::ATan(tr1->GetMomDir().X())*1000;//mrad
      blc1par[2]=tr1->GetPosatZ(0).Y(); //cm
      blc1par[3]=TMath::ATan(tr1->GetMomDir().Y())*1000;//mrad
      blc1par[4]=0.;
      blc1par[5]=beam->param(4);
      beam->CalcParBLC1toBLC2(blc1par,parblc2);
      double dp=beam->param(4);
      double a= tr2->gx();
      double b= tr2->gdx();
      double c= tr2->gy();
      double d= tr2->gdy();
      if(BEAM&&TOFK){
        ofsk<<a<<"\t"<<b<<"\t"<<c<<"\t"<<d<<"\t"<<dp<<std::endl;
      }
      if(BEAM&&TOFPi){
        ofspi<<a<<"\t"<<b<<"\t"<<c<<"\t"<<d<<"\t"<<dp<<std::endl;
      }
      //    std::cout<<parblc2[0]<<std::endl;
      TVector3 pos2=tr2->GetPosatZ(-130.*cm);
      TVector3 dir2=tr2->GetMomDir();
      double tmpbins[6]={200,-10*cm,10*cm,200,-10*cm,10*cm};
      hist::H2("CorrX_BLC1BLC2",pos2.X(),parblc2[0],tmpbins);
      hist::H2("CorrY_BLC1BLC2",pos2.Y(),parblc2[2],tmpbins);
      hist::H2("X_diffY_BLC1BLC2",pos2.X(),pos2.Y()-parblc2[2],tmpbins);
      hist::H2("Y_diffX_BLC1BLC2",pos2.Y(),pos2.X()-parblc2[0],tmpbins);
      hist::H2("CorrA_BLC1BLC2",TMath::ATan(dir2.X())*1000,parblc2[1],tmpbins);
      hist::H2("CorrB_BLC1BLC2",TMath::ATan(dir2.Y())*1000,parblc2[3],tmpbins);
      hist::H2("CorrXY_BLC1BLC2",pos2.X()-parblc2[0],pos2.Y()-parblc2[2],tmpbins);
      hist::H2("CorrAB_BLC1BLC2",
	       TMath::ATan(dir2.X())*1000-parblc2[1],
	       TMath::ATan(dir2.Y())*1000-parblc2[3],tmpbins);
#if 0
      std::cout<<"==========================="<<std::endl;
      std::cout<<"D5 mom, chi2: "<< mom<<"  "<<chi<<std::endl;
      std::cout<<"X: "<<blc1par[0]<<" , "<<beam->param(0)<<" ->  "<<pos2.X()<<" , "<<parblc2[0]<<std::endl;
      std::cout<<"Y: "<<blc1par[2]<<" , "<<beam->param(2)<<" ->  "<<pos2.Y()<<" , "<<parblc2[2]<<std::endl;
      std::cout<<"A: "<<blc1par[1]<<" ->  "<<TMath::ATan(dir2.X())*1000<<" , "<<parblc2[1]<<std::endl;
      std::cout<<"B: "<<blc1par[3]<<" ->  "<<TMath::ATan(dir2.Y())*1000<<" , "<<parblc2[3]<<std::endl;
      tr1->Print();
      tr2->Print();
#endif
    }
  }
#endif
#endif //Tracking
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessEnd()
{
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<HodoParamMan>("HDPRM")) &&
    (InitializeParameter<HodoPHCMan>("HDPHC")) &&
    (InitializeParameter<DCTdcCalibMan>("DCTDC")) &&
    (InitializeParameter<XTMapMan>("XTMap")) &&
    (InitializeParameter<GeomMapMan>("GeomBL","GeomHall")) &&
    (InitializeParameter<UserParamMan>("USER")) &&
    (InitializeParameter<TransferMatrixMan>("TM")) &&
    (InitializeParameter<DCTimeCorrMan>("DCTC")) &&
    (InitializeParameter<BLDCWireMapMan>("BLDCWire"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  ofsk.close();
  ofspi.close();
  return true;
}
