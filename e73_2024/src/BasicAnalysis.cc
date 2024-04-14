#include "BasicAnalysis.hh"
#include "HodoPHCMan.hh"
#include "HodoCluster.hh"
#include "HistTools.hh"
#include "GeomTools.hh"
#include "CDSAnalysis.hh"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "DatabasePDG.hh"
#include "DeleteUtility.hh"
#if TKIN // given in make file
#include "KinematicalFit.hh"
#endif

namespace{
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const HodoPHCMan& gPHC = HodoPHCMan::GetInstance();
#if TKIN
  KinematicalFit& gKin  = KinematicalFit::GetInstance();
#endif
  const double mm=0.1;
  const double cm=10*mm;

  double pi_beta=1./TMath::Hypot(1.,pdg::Mass(-211)); // at 1 GeV/c

  double tdcbins[3]={10000,0,2e6};
  double adcbins[3]={4096,-0.5,4095.5};
  double diffbins[3]={4000,-10000,10000};
  double mtdcbins[3]={2000,0,2000};
  double mulbins2[3]={10,-0.5,9.5};
  double timebins[3]={10000,-100,100};
  double tofbins[3]={10000,-100,100};
  double debins[3]={5000,-2,48};
  double evdebins[6]={200,0,5e6,200,-2,18};
  double evtimebins[6]={200,0,5e6,200,-10,10};
  double posbins2[6]={150,-15*cm,15*cm,150,-15*cm,15*cm};
  double residbins[3]={200,-5*mm,5*mm};
}

BasicAnalysis::BasicAnalysis()
{
}

bool BasicAnalysis::AnaInit(MTDCAnalyzer* mtdc, HodoAnalyzer* hodo, bool mc, bool slewing, bool yamaga){
  EventNumber=0;
  Time0=CTime0=-999;
  Time1=CTime1=-999;
  TOF_BHTT0=cTOF_BHTT0=-999;
  TOFK=TOFP=TOFPi=ACHIT=false;
  isMC=mc;
  sign_beam=isMC ? 1 : -1;
  isYamaga=yamaga;
  position[kT0]=TVector3(-999,-999,-999);
  position[kT1]=TVector3(-999,-999,-999);
  position[kDEF]=TVector3(-999,-999,-999);
  return AnaFlag(mtdc)
    && CalcTime0(hodo)
    && CheckBHTT0TOF(hodo,slewing);
}

bool BasicAnalysis::AnaFlag(MTDCAnalyzer* mtdc){
  int cid=DetIdTrigFlag;
  TString tmpname="TriggerFlag";
  int nh = mtdc->GetNHits(cid);

  double patbins[3]={kNumTrig,-0.5,kNumTrig-0.5};
  for(int i=0;i<kNumTrig;i++) TriggerFlag[i]=false;

  for( int i=0; i<nh; ++i ){
    const auto& hit = mtdc->GetHit(cid,i);
    const auto& raw = hit->GetRawHit();
    int seg = hit->SegmentId();
    TString segstr=Form("_seg%d",seg);
    for(int it=0;it<raw->GetSizeLeading();it++){
      double tu  = raw->GetLeading(it);
      hist::H1(tmpname+"_TDC"+segstr,tu,mtdcbins);
    }
    bool tmpflag=false;
    for(int it=0;it<hit->GetIndex();it++){
      double tu  = hit->Time(it);
      hist::H1(tmpname+"_Time"+segstr,tu,timebins);
      //if(gUser.Check("TriggerFlag",tu)) tmpflag=true;
      if(fabs(tu)<10) tmpflag=true;
    }
    if(tmpflag){
      hist::H1(tmpname+"_Pat",seg,patbins);
      TriggerFlag[seg]=true;
    }
  }
  for(int i=0;i<16;i++)
    for(int j=i;j<16;j++)
      if(TriggerFlag[i]&&TriggerFlag[j])
	hist::H2(tmpname+"_2DPat",i,j,patbins,patbins);
  {
    TriggerMode=0;
    if(TriggerFlag[kBeam])		TriggerMode=1;
    else if(TriggerFlag[kKGamma])	TriggerMode=2;
    else if(TriggerFlag[kKCDH1])	TriggerMode=3;
    else if(TriggerFlag[kKCDH1Gamma])	TriggerMode=4;
    else if(TriggerFlag[kKCDH2])	TriggerMode=5;
    else if(TriggerFlag[kPiCDH])	TriggerMode=6;
    else if(TriggerFlag[kBeamPbF2])	TriggerMode=7;
    hist::H1("TriggerMode",TriggerMode,mulbins2);
  }
  return true;
}

bool BasicAnalysis::CalcTime0(HodoAnalyzer* hodo){
  int cid= DetIdT1;
  int nh = hodo->GetNHits(cid);
  for( int i=0; i<nh; ++i ){
    auto hit = hodo->GetHit(cid,i);
    if(!hit) continue;
    int nind=hit->GetIndex();
    for(int it=0;it<nind;it++){
      double mt  = hit->MeanTime(it);
      double cmt = hit->CMeanTime(it);
      //std::cout<<i<<"  "<<it<<"  "<<mt<<"  "<<cmt<<std::endl;
      hist::H1("T0MeanTime" , mt,tofbins);
      hist::H1("T0CMeanTime",cmt,tofbins);
      if(TMath::Abs(mt)<10){
	Time0=mt;
	CTime0=cmt;
	hist::H1("T0MeanTime_selected" ,mt ,tofbins);
	hist::H1("T0CMeanTime_selected",cmt,tofbins);
      }
    }
  }
  if(Time0<-100) return false;
  cid= DetIdT0;
  nh = hodo->GetNHits(cid);
  for( int i=0; i<nh; ++i ){
    auto hit = hodo->GetHit(cid,i);
    if(!hit) continue;
    int nind=hit->GetIndex();
    for(int it=0;it<nind;it++){
      double mt  = hit->MeanTime(it);
      double cmt = hit->CMeanTime(it);
      //std::cout<<i<<"  "<<it<<"  "<<mt<<"  "<<cmt<<std::endl;
      hist::H1("T1MeanTime" , mt,tofbins);
      hist::H1("T1CMeanTime",cmt,tofbins);
      if(TMath::Abs(mt)<10){
	Time1=mt;
	CTime1=cmt;
	hist::H1("T1MeanTime_selected" ,mt ,tofbins);
	hist::H1("T1CMeanTime_selected",cmt,tofbins);
      }
    }
  }
  if(Time0<-100||Time1<-100) return false;
  //std::cout<<"failed in calculating time0"<<std::endl;
  return true;
}

bool BasicAnalysis::CheckBHTT0TOF(HodoAnalyzer* hodo, bool slewing){
  int cid=DetIdBHT;
  int nh = hodo->GetNHits(cid);
  for( int i=0; i<nh; ++i ){
    auto hit = hodo->GetHit(cid,i);
    int seg = hit->SegmentId();
    if(seg<3||12<seg) continue;
    if(!hit) continue;
    int nind=hit->GetIndex();
    for(int it=0;it<nind;it++){
      double mt    = hit->MeanTime(it);
      double cmt   = hit->CMeanTime(it);
      double tof   = mt-time0();
      double ctof  = cmt-ctime0();
      hist::H1("BHTMeanTime" ,mt ,tofbins);
      hist::H1("BHTCMeanTime",cmt,tofbins);
      hist::H1("BHTT0TOF" ,tof ,tofbins);
      hist::H1("BHTT0cTOF",ctof,tofbins);
      if(TMath::Abs(mt)<10){
	TOF_BHTT0=tof;
	cTOF_BHTT0=ctof;
	hist::H1("BHTMeanTime_selected" ,mt ,tofbins);
	hist::H1("BHTCMeanTime_selected",cmt,tofbins);
	hist::H1("BHTT0TOF_selected" ,tof,tofbins);
	hist::H1("BHTT0cTOF_selected",ctof,tofbins);
	if(gUser.Check("TOFK",ctof)){
	  TOFK=true;
	  hist::H1("BHTT0cTOF_TOFK",ctof,tofbins);
	  //	  return true;
	}
	if(gUser.Check("TOFP",ctof)){
	  TOFP=true;
	  hist::H1("BHTT0cTOF_TOFP",ctof,tofbins);
	}
	if(gUser.Check("TOFPi",ctof)){
	  TOFPi=true;
	  hist::H1("BHTT0cTOF_TOFPi",ctof,tofbins);
	}
      }
    }
  }
  return true;
}

bool BasicAnalysis::AnaAC(RawData* raw){
  if(!raw) return false;
  int cid=DetIdAC;
  const HodoRHitContainer &cont = raw->GetHodoRawHC(cid);
  int nh = cont.size();
  TString tmpname="AC";
  ACHIT=false;
  for( int i=0; i<nh; ++i ){
    HodoRawHit *raw = cont[i];
    if(!raw) continue;
    int ntu=raw->GetSizeTdcUp();
    for(int it=0;it<ntu;it++){
      double tu  = raw->GetTdcUp(it);
      if(gUser.Check("ACTDC",tu)) ACHIT=true;
      hist::H1(tmpname+"_TDC",tu,tdcbins);
    }
    double adcsum=0;
    for(int ia=0;ia<raw->GetSizeAdcUp();ia++)
      adcsum+=raw->GetAdcUp(ia);
    hist::H1(tmpname+"_ADC",adcsum,adcbins);
    if(trig(kBeam)){
      hist::H1(tmpname+"_ADC_ifB",adcsum,adcbins);
      if(tofk()){
	hist::H1(tmpname+"_ADC_ifB_TOFK",adcsum,adcbins);
	if(achit()){
	  hist::H1(tmpname+"_ADC_ifB_TOFK_AChit",adcsum,adcbins);
	}
      }else if(tofpi()){
	hist::H1(tmpname+"_ADC_ifB_TOFPi",adcsum,adcbins);
	if(achit()){
	  hist::H1(tmpname+"_ADC_ifB_TOFPi_AChit",adcsum,adcbins);
	}
      }
    }
  }
  return ACHIT;
}

bool BasicAnalysis::FillHodoRaw(RawData* raw,int kHodo){
  int mulu=0,muld=0;
  int cid=hodoid[kHodo];
  TString name=hodoname[kHodo];
  double mulbins[3]={nsegs[kHodo]+1,-0.5,nsegs[kHodo]+0.5};
  double patbins[3]={nsegs[kHodo],-0.5,nsegs[kHodo]-0.5};
  const HodoRHitContainer &cont = raw->GetHodoRawHC(cid);
  int nh = cont.size();
  int hitsegu=-1,hitsegd=-1;
  for( int i=0; i<nh; ++i ){
    HodoRawHit *raw = cont[i];
    if(!raw) continue;
    bool TIMING=false;
    int seg = raw->SegmentId();
    TString segstr=Form("_seg%d",seg);
    double au  = raw->GetAdcUp();
    double ad  = raw->GetAdcDown();
    int ntu=raw->GetSizeTdcUp();
    int ntd=raw->GetSizeTdcDown();
    hist::H1(name+"_ADCu"+segstr,au,adcbins);
    hist::H1(name+"_ADCd"+segstr,ad,adcbins);
    hist::H1(name+"_Mulu"+segstr,ntu,mulbins2);
    hist::H1(name+"_Muld"+segstr,ntd,mulbins2);
    hist::H2(name+"_Mulu_Muld"+segstr,ntu,ntd,mulbins2,mulbins2);
    hist::H2(name+"_Mulu_pat",seg,ntu,patbins,mulbins2);
    hist::H2(name+"_Muld_pat",seg,ntd,patbins,mulbins2);
    if(ntu>0){
      hist::H1(name+"_Patu",seg,patbins);
      hist::H1(name+"_ADCwTu"+segstr,au,adcbins);
      mulu++;
      hitsegu=seg;
    }else{
      hist::H1(name+"_ADCwoTu"+segstr,au,adcbins);
    }
    if(ntd>0){
      hist::H1(name+"_Patd",seg,patbins);
      hist::H1(name+"_ADCwTd"+segstr,ad,adcbins);
      muld++;
      hitsegd=seg;
    }else{
      hist::H1(name+"_ADCwoTd"+segstr,ad,adcbins);
    }
    for(int it=0;it<ntu;it++){
      double tu  = raw->GetTdcUp(it);
      hist::H1(name+"_TDCu"+segstr,tu,tdcbins);
      if(it<ntd){
	double td  = raw->GetTdcDown(it);
	hist::H1(name+"_TDCdiffud",tu-td,diffbins);
      }
      if(gUser.Check("HODOTDC",tu)) TIMING=true;
    }
    for(int it=0;it<ntd;it++){
      double td  = raw->GetTdcDown(it);
      hist::H1(name+"_TDCd"+segstr,td,tdcbins);
    }
    if(TIMING&&trig(kBeam)&&achit()){
      hist::H1(name+"_ADCwBu"+segstr,au,adcbins);
      hist::H1(name+"_ADCwBd"+segstr,ad,adcbins);
    }
  }//for(i)
  if(mulu==1&&muld==1){
    hist::H2(name+"_UDcorr" ,hitsegu,hitsegd ,patbins,patbins);
  }
  hist::H1(name+"_Mulu",mulu,mulbins);
  hist::H1(name+"_Muld",muld,mulbins);
  return true;
}

bool BasicAnalysis::FillHodoDecoded(HodoAnalyzer* hodo,int kHodo,std::vector<TString> trig_add)
{
  if(kHodo==kPbF2||kHodo==kPbF2)
    return FillHodo1Decoded(hodo,kHodo,trig_add);
  int mul=0;
  int cid=hodoid[kHodo];
  TString name=hodoname[kHodo];
  double mulbins[3]={nsegs[kHodo]+1,-0.5,nsegs[kHodo]+0.5};
  double patbins[3]={nsegs[kHodo],-0.5,nsegs[kHodo]-0.5};
  int nh = hodo->GetNHits(cid);
  for( int i=0; i<nh; ++i ){
    Hodo2Hit *hit = hodo->GetHit(cid,i);
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
    double de = hit->DeltaE();
    double deu = hit->GetAUp();
    double ded = hit->GetADown();
    hist::H1(name+"_dE" +segstr,de ,debins,trig_add);
    hist::H1(name+"_dEu"+segstr,deu,debins,trig_add);
    hist::H1(name+"_dEd"+segstr,ded,debins,trig_add);
    double slewbins[6]={100,0,10,400,-2.,2.};
    for(int ii=0;ii<nind;ii++){
      double mt   = hit->MeanTime(ii);
      double cmt  = hit->CMeanTime(ii);
      double tof0 = mt - time0();
      double ctof0= cmt - ctime0();
      double tof1 = mt - time1();
      double ctof1= cmt - ctime1();
      double fl0  = (pos(kHodo) - pos(kT1)).Mag();
      double ft0  = fl0/pi_beta/(TMath::C()*1e-6*mm);
      double fl1  = (pos(kHodo) - pos(kT0)).Mag();
      double ft1  = fl1/pi_beta/(TMath::C()*1e-6*mm);
      //
      double tof = cid==DetIdT1 ? tof1  : tof0;
      double ctof= cid==DetIdT1 ? ctof1 : ctof0;
      double fl  = cid==DetIdT1 ? fl1 : fl0;
      double ft  = cid==DetIdT1 ? ft1 : ft0;
      if(fabs(mt)<20){
	hist::H2(name+"_EvNum_dE" ,EventNumber,de,evdebins);
	hist::H2(name+"_EvNum_MeanTime"  ,EventNumber,mt ,evtimebins);
	hist::H2(name+"_EvNum_CMeanTime" ,EventNumber,cmt,evtimebins);
	hist::H2(name+"_EvNum_TOF"  ,EventNumber,tof ,evtimebins);
	hist::H2(name+"_EvNum_cTOF" ,EventNumber,ctof,evtimebins);
	hist::H2(name+"_EvNum_dT"  ,EventNumber,tof-ft ,evtimebins);
	hist::H2(name+"_EvNum_cdT" ,EventNumber,ctof-ft,evtimebins);
      }
      hist::H1(name+"_MeanTime" ,mt,timebins,trig_add);
      hist::H1(name+"_MeanTime" +segstr,mt,timebins,trig_add);
      hist::H1(name+"_CMeanTime",cmt,timebins,trig_add);
      hist::H1(name+"_CMeanTime"+segstr,cmt,timebins,trig_add);
      if(achit()&&tofpi()){
	hist::H2(name+"_dECMeanTime",de,cmt,slewbins);
      }
      if(time0()>-9000){
	hist::H1(name+"_TOF" ,tof,tofbins,trig_add);
	hist::H1(name+"_TOF" +segstr,tof,tofbins,trig_add);
	hist::H1(name+"_cTOF0",ctof0,tofbins,trig_add);
	hist::H1(name+"_cTOF0"+segstr,ctof0,tofbins,trig_add);
	hist::H1(name+"_cTOF1",ctof1,tofbins,trig_add);
	hist::H1(name+"_cTOF1"+segstr,ctof1,tofbins,trig_add);
	hist::H1(name+"_dT" ,tof-ft,tofbins,trig_add);
	hist::H1(name+"_dT" +segstr,tof-ft,tofbins,trig_add);
	hist::H1(name+"_cdT0",ctof0-ft0,tofbins,trig_add);
	hist::H1(name+"_cdT0"+segstr,ctof0-ft0,tofbins,trig_add);
	hist::H1(name+"_cdT1",ctof1-ft1,tofbins,trig_add);
	hist::H1(name+"_cdT1"+segstr,ctof1-ft1,tofbins,trig_add);
	if(achit()&&tofpi()){
	  hist::H2(name+"_dEcTOF"       ,de,ctof,slewbins);
	  hist::H2(name+"_dEcTOF"+segstr,de,ctof,slewbins);
	  hist::H2(name+"_dEcdT"        ,de,ctof-ft,slewbins);
	  hist::H2(name+"_dEcdT" +segstr,de ,ctof-ft,slewbins);
	  hist::H2(name+"_dEucdT"       ,deu,ctof-ft,slewbins);
	  hist::H2(name+"_dEdcdT"+segstr,ded,ctof-ft,slewbins);
	}
      }
    }
    // if(trig(kBeam)&&achit()){
    //   hist::H2("Run_"+name+"_dEuwB"+segstr,run_number,deu,400,100,900,400,-1,19);
    //   hist::H2("Run_"+name+"_dEdwB"+segstr,run_number,ded,400,100,900,400,-1,19);
    // }
  }//for(ihit)
  hist::H1(name+"_Mul",mul ,mulbins);
  //  hist::H1(name+"_Mulgate",hodomul,mulbins);
  return true;
}

bool BasicAnalysis::FillHodo1Decoded(HodoAnalyzer* hodo,int kHodo,std::vector<TString> trig_add)
{
  int mul=0;
  int cid=hodoid[kHodo];
  TString name=hodoname[kHodo];
  double mulbins[3]={nsegs[kHodo]+1,-0.5,nsegs[kHodo]+0.5};
  double patbins[3]={nsegs[kHodo],-0.5,nsegs[kHodo]-0.5};
  int nh = hodo->GetNHits(cid);
  for( int i=0; i<nh; ++i ){
    Hodo1Hit *hit = hodo->Get1Hit(cid,i);
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
    double de = hit->GetAUp();
    hist::H1(name+"_dE" +segstr,de ,debins,trig_add);
    for(int ii=0;ii<nind;ii++){
      double mt  = hit->Time(ii);
      double cmt  =hit->CTime(ii);
      double tof = mt - time0();
      double ctof = cmt - ctime0();
      hist::H1(name+"_Time" ,mt,timebins,trig_add);
      hist::H1(name+"_Time" +segstr,mt,timebins,trig_add);
      hist::H1(name+"_CTime",cmt,timebins,trig_add);
      hist::H1(name+"_CTime"+segstr,cmt,timebins,trig_add);
      if(time0()>-9000){
	hist::H1(name+"_TOF" ,tof,tofbins,trig_add);
	hist::H1(name+"_TOF" +segstr,tof,tofbins,trig_add);
	hist::H1(name+"_cTOF",ctof,tofbins,trig_add);
	hist::H1(name+"_cTOF"+segstr,ctof,tofbins,trig_add);
      }
    }
  }//for(ihit)
  hist::H1(name+"_Mul",mul ,mulbins);
  //  hist::H1(name+"_Mulgate",hodomul,mulbins);
  return true;
}

bool BasicAnalysis::AnaDC(DCAnalyzer* dc,int kDC,double tmin,double tmax,double chi2)
{
  int cid=chmid[kDC];
  TString name=chmname[kDC];
  int status=dc->TrackSearch(cid,!isMC);
  int ntra=dc->GetNTracks(cid);
  for( int i=0; i<ntra; i++ ){
    LocalTrack *tr = dc->GetTrack(cid,i);
    hist::H1( name+"_tracktime", tr->GetTrackTime(),2000,-100,100 );
    hist::H1( name+"_chi2"     , tr->chi2all(),1000,0,200 );
    if( tr->CheckRange(tmin,tmax) && tr->chi2all()<chi2 ){
      good_trackid[kDC].push_back(i);
    }
  }
  hist::H1(name+"_nTracks_all",   ntra,10,-0.5,9.5);
  hist::H1(name+"_nTracks", ngood(kDC),10,-0.5,9.5);
  hist::H1(name+"_status",      status,10,-0.5,9.5);
  if(trig(kKaon3)){
    hist::H1(name+"_nTracks_all_K3",   ntra,10,-0.5,9.5);
    hist::H1(name+"_nTracks_K3", ngood(kDC),10,-0.5,9.5);
  }
  if(ngood(kDC)==1){
    LocalTrack* track=dc->GetTrack(cid,goodid(kDC));
    double x,y;
    track->XYLocalPosatZ(0,x,y);
    hist::H2(name+"_XYLocal",x,y,posbins2);
    hist::H2(name+"_AB",track->gdx(),track->gdy(),100,-0.1,0.1,100,-0.1,0.1);
    for(int xy=0;xy<2;xy++){
      for(int i=0;i<track->nhit(xy);i++){
	double resid=track->resid(xy,i);
	int layer=track->layer(xy,i);
	hist::H1(name+"_resid"+Form("_layer%d", layer ), resid, residbins);
      }
    }
    double z=0;
    if(cid==DetIdBLC2||cid==DetIdBLC2a||cid==DetIdBLC2b){
      z=-130*cm;
      if(hodoseg(kT0)>-1){
	track->XYLocalPosatZ(20*cm,x,y);
	hist::H2(name+"_XY"+Form("_ifT0seg%d",hodoseg(kT0)),x,y,posbins2);
      }
    }else if(cid==DetIdBPC||cid==DetIdSDC){
      if(hodoseg(kDEF)>-1){
	track->XYLocalPosatZ(5*cm,x,y);
	hist::H2(name+"_XY"+Form("_ifDEFseg%d",hodoseg(kDEF)),x,y,posbins2);
      }
    }
    track->XYPosatZ(z,x,y);
    hist::H2(name+"_XY",x,y,posbins2);
  }
  return true;
}

bool BasicAnalysis::AnaHodo(HodoAnalyzer* hodo,int kHodo)
{
  if(kHodo==kPbG||kHodo==kPbF2)
    return AnaHodo1(hodo,kHodo);
  if(kHodo==kBHT)
    return AnaBHT(hodo);
  int cid=hodoid[kHodo];
  TString name=hodoname[kHodo];
  int nh = hodo->GetNHits(cid);
  //  std::cout<<kHodo<<"  "<<cid<<"  "<<name<<std::endl;
  for( int i=0; i<nh; ++i ){
    Hodo2Hit *hit = hodo->GetHit(cid,i);
    if(!hit) continue;
    int nind= hit->GetIndex();
    for(int ii=0;ii<nind;ii++){
      double mt  = hit->MeanTime(ii)-time0();
      double cmt  = hit->CMeanTime(ii)-ctime0();
      hist::H1(name+"_Time_all" ,mt,timebins);
      hist::H1(name+"_CTime_all" ,cmt,timebins);
      if(cid==DetIdCDH){
	if(!gUser.Check("CDHTOF",mt)) continue;
      }else if(cid==DetIdVeto){
	if(!gUser.Check("VetoTOF",mt)) continue;
      }else{
	if(!gUser.Check("HodoGATE",mt)) continue;
      }
      hist::H1(name+"_Time_gate" ,mt,timebins);
      hist::H1(name+"_CTime_gate",cmt,timebins);
      HodoEvent ev;
      ev.cid=cid;
      ev.seg= hit->SegmentId();
      ev.time=hit->MeanTime(ii);
      ev.ctime=hit->CMeanTime(ii);
      ev.deu= hit->GetAUp();
      ev.ded= hit->GetADown();
      ev.tu = hit->GetTUp();
      ev.td = hit->GetTDown();
      ev.ctu = hit->GetCTUp();
      ev.ctd = hit->GetCTDown();
      hodo_container[kHodo].push_back(ev);
      break;
    }
  }//for(ihit)
  return true;
}

bool BasicAnalysis::AnaBHT(HodoAnalyzer* hodo,int kHodo)
{
  int cid=hodoid[kHodo];
  TString name=hodoname[kHodo];
  int nh = hodo->GetNHits(cid);
  //Clustering
  //  std::cout<<"clustering"<<std::endl;
  std::vector<HodoCluster*> clusterContainer;
  std::vector<std::vector<int>> UsedFlag(nh,std::vector<int>(0));
  for( int i=0; i<nh; ++i ){
    Hodo2Hit *hit1 = hodo->GetHit(cid,i);
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
	Hodo2Hit *hit2 = hodo->GetHit(cid,i+j);
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

  //  std::cout<<"fill"<<std::endl;
  // Fill hits
  int mulcluster=clusterContainer.size();
  for(int i=0; i<mulcluster; i++){
    HodoCluster* clu=clusterContainer.at(i);
    Hodo2Hit* hit=clu->Get1stHit();
    int       nth=clu->Get1stNthHit();
    double mt = hit->MeanTime(nth);
    double cmt= hit->CMeanTime(nth);
    hist::H1(name+"_Time_all" ,mt,timebins);
    hist::H1(name+"_CTime_all" ,cmt,timebins);
    if(!gUser.Check("HodoGATE",mt)) continue;
    hist::H1(name+"_Time_gate" ,mt,timebins);
    hist::H1(name+"_CTime_gate",cmt,timebins);
    HodoEvent ev;
    ev.cid=cid;
    ev.seg= hit->SegmentId();
    ev.time=hit->MeanTime(nth);
    ev.ctime=hit->CMeanTime(nth);
    ev.deu= hit->GetCAUp(nth);
    ev.ded= hit->GetCADown(nth);
    ev.tu = hit->GetTUp(nth);
    ev.td = hit->GetTDown(nth);
    ev.ctu = hit->GetCTUp(nth);
    ev.ctd = hit->GetCTDown(nth);
    hodo_container[kHodo].push_back(ev);
  }
  del::ClearContainer( clusterContainer );
  return true;
}

bool BasicAnalysis::AnaHodo1(HodoAnalyzer* hodo,int kHodo)
{
  int cid=hodoid[kHodo];
  TString name=hodoname[kHodo];
  int nh = hodo->GetNHits(cid);
  for( int i=0; i<nh; ++i ){
    Hodo1Hit *hit = hodo->Get1Hit(cid,i);
    if(!hit) continue;
    int nind= hit->GetIndex();
    // if(kHodo==kVeto0){
    //   HodoEvent ev;
    //   ev.cid=cid;
    //   ev.seg= hit->SegmentId();
    //   ev.deu= hit->GetAUp();
    //   hodo_container[kHodo].push_back(ev);
    //   break;
    // }
    for(int ii=0;ii<nind;ii++){
      double mt  = hit->Time(ii);
      if(!gUser.Check("HodoGATE",mt)) continue;
      HodoEvent ev;
      ev.cid=cid;
      ev.seg= hit->SegmentId();
      ev.time=mt;
      ev.ctime=hit->CTime(ii);
      ev.deu= hit->GetAUp();
      ev.ded= hit->GetAUp();
      hodo_container[kHodo].push_back(ev);
      break;
    }
  }//for(ihit)
  return true;
}

bool BasicAnalysis::AnaPbF2(HodoAnalyzer* hodo,std::vector<TString> add)
{
  int cid=hodoid[kPbF2];
  TString name=hodoname[kPbF2];
  int nh = hodo->GetNHits(cid);
  std::vector<int> tmphit;
  PbF2_AdcSum=0;
  PbF2_EnergySum=0;
  PbF2_EnergySumwT=0;
  PbF2_EnergySum2=0;
  PbF2_EnergySum2wT=0;
  PbF2_EnergySum3=0;
  PbF2_Seg=-1;
  double tmpde[40];
  for( int i=0; i<40; ++i ) tmpde[i]=0.;
  double tmpemax=0;
  double tmpbins[3]={1000,-2,8};
  for( int i=0; i<nh; ++i ){
    Hodo1Hit *hit = hodo->Get1Hit(cid,i);
    if(!hit) continue;
    int seg = hit->SegmentId();
    TString segstr=Form("_seg%d",seg);
    double de = hit->GetAUp();
    double de2= de;
    if(!isMC){
      gPHC.DoCorrection(DetIdPbF2,0,seg,0,0,de,de2);
    }
    tmpde[seg]=de2;
    if(!isMC){
      double adc= hit->GetRawHit()->GetAdcUp();
      PbF2_AdcSum+=adc;
    }
    PbF2_EnergySum+=de;
    PbF2_EnergySum2+=de2;
    if(de2>0.03){
      PbF2_EnergySum3+=de2;
      hist::H1(name+"_Energy3"+segstr, de2 ,tmpbins,add);
    }
    hist::H1(name+"_Energy"+segstr, de ,tmpbins,add);
    hist::H1(name+"_Energy2"+segstr, de2 ,tmpbins,add);
    int nind= hit->GetIndex();
    for(int ii=0;ii<nind;ii++){
      double tof  = hit->Time(ii)-time0();
      hist::H1(name+"_TOF_all", tof ,tofbins,add);
      if(!gUser.Check("PbF2TOF",tof)) continue;
      hist::H1(name+"_TOF_gate", tof ,tofbins,add);
      //      std::cout<<"seg,adc,de,de2: "<<seg<<"  "<<adc<<"  "<<de<<"  "<<de2<<std::endl;
      PbF2_EnergySumwT+=de;
      PbF2_EnergySum2wT+=de2;
      if(de>tmpemax){
	PbF2_Seg=seg;
	PbF2_Time=tof;
	PbF2_EnergyMax=de2;
	tmpemax=de;
      }
      tmphit.push_back(i);
      break;
    }
  }//for(ihit)
  if(pbf2_seg()>-1){
    int seg=pbf2_seg();
    double de_cluster=tmpde[seg];
    if(seg%8!=0)  de_cluster+=tmpde[seg-1];
    if(seg%8!=7)  de_cluster+=tmpde[seg+1];
    if(seg>7){
      de_cluster+=tmpde[seg-8];
      if(seg%8!=0)  de_cluster+=tmpde[seg-9];
      if(seg%8!=7)  de_cluster+=tmpde[seg-7];
    }
    if(seg<32){
      de_cluster+=tmpde[seg+8];
      if(seg%8!=7)  de_cluster+=tmpde[seg+9];
      if(seg%8!=0)  de_cluster+=tmpde[seg+7];
    }
    PbF2_EnergySum4=de_cluster;
    TString segstr=Form("_seg%d",pbf2_seg());
    hist::H1(name+"_EnergySum", pbf2_enesum() ,tmpbins,add);
    hist::H1(name+"_EnergySumwT", pbf2_enesumwt() ,tmpbins,add);
    hist::H1(name+"_EnergySum2", pbf2_enesum2() ,tmpbins,add);
    hist::H1(name+"_EnergySum2wT", pbf2_enesum2wt() ,tmpbins,add);
    hist::H1(name+"_EnergySum3", pbf2_enesum3() ,tmpbins,add);
    hist::H1(name+"_EnergySum4", pbf2_enesum4() ,tmpbins,add);
    hist::H1(name+"_ADCSum", pbf2_adcsum() ,adcbins,add);

    hist::H1(name+"_EnergySum"+segstr, pbf2_enesum() ,tmpbins,add);
    hist::H1(name+"_EnergySumwT"+segstr, pbf2_enesumwt() ,tmpbins,add);
    hist::H1(name+"_EnergySum2"+segstr, pbf2_enesum2() ,tmpbins,add);
    hist::H1(name+"_EnergySum2wT"+segstr, pbf2_enesum2wt() ,tmpbins,add);
    hist::H1(name+"_EnergySum3"+segstr, pbf2_enesum3() ,tmpbins,add);
    hist::H1(name+"_EnergySum4"+segstr, pbf2_enesum4() ,tmpbins,add);
    hist::H1(name+"_ADCSum"+segstr, pbf2_adcsum() ,adcbins,add);
  }
  if(tmphit.size()==1){
    Hodo1Hit *hit = hodo->Get1Hit(cid,tmphit[0]);
    if(hit){
      int seg = hit->SegmentId();
      TString segstr=Form("_seg%d",seg);
      double de = hit->GetAUp();
      double de2= 0;
      gPHC.DoCorrection(DetIdPbF2,0,seg,0,0,de,de2);
      hist::H1(name+"_Energy_single"+segstr, de ,tmpbins,add);
      hist::H1(name+"_Energy2_single"+segstr, de2 ,tmpbins,add);
      if(!isMC){
	double adc= hit->GetRawHit()->GetAdcUp();
	hist::H1(name+"_ADC_single"+segstr, adc ,adcbins,add);
      }
    }
  }
  return true;
}
SlewEvent BasicAnalysis::get_slewevent()
{
  SlewEvent ev;
  ev.mt_t0new =hodotime( kT1,0);
  ev.cmt_t0new=hodoctime(kT1,0);
  ev.deu_t0new=hododeu(  kT1,0);
  ev.ded_t0new=hododed(  kT1,0);
  ev.seg_t0new=hodoseg(  kT1,0);

  ev.mt_t0 =hodotime( kT0,0);
  ev.cmt_t0=hodoctime(kT0,0);
  ev.deu_t0=hododeu(  kT0,0);
  ev.ded_t0=hododed(  kT0,0);
  ev.seg_t0=hodoseg(  kT0,0);

  ev.mt_bht =hodotime( kBHT,0);
  ev.cmt_bht=hodoctime(kBHT,0);
  ev.deu_bht=hododeu(  kBHT,0);
  ev.ded_bht=hododed(  kBHT,0);
  ev.seg_bht=hodoseg(  kBHT,0);

  ev.mt_def =hodotime( kDEF,0);
  ev.cmt_def=hodoctime(kDEF,0);
  ev.deu_def=hododeu(  kDEF,0);
  ev.ded_def=hododed(  kDEF,0);
  ev.seg_def=hodoseg(  kDEF,0);

  ev.mt_veto =hodotime( kVeto,0);
  ev.cmt_veto=hodoctime(kVeto,0);
  ev.deu_veto=hododeu(  kVeto,0);
  ev.ded_veto=hododed(  kVeto,0);
  ev.seg_veto=hodoseg(  kVeto,0);

  ev.mt_btc =hodotime( kBTC,0);
  ev.cmt_btc=hodoctime(kBTC,0);
  ev.deu_btc=hododeu(  kBTC,0);
  ev.ded_btc=hododed(  kBTC,0);
  ev.seg_btc=hodoseg(  kBTC,0);

  return ev;
}

DeuteronEvent BasicAnalysis::get_devent(HodoAnalyzer *hodo)
{
  DeuteronEvent ev;
  ev.time_t0new = hodotime(kT1,0);
  ev.de_t0new	= hodode(  kT1,0);
  ev.seg_t0new	= hodoseg( kT1,0);

  ev.time_t0	= hodotime(kT0,0);
  ev.de_t0	= hodode(  kT0,0);
  ev.seg_t0	= hodoseg( kT0,0);

  ev.time_bht	= hodotime(kBHT,0);
  ev.de_bht	= hodode(  kBHT,0);
  ev.seg_bht	= hodoseg( kBHT,0);

  ev.time_def	= hodotime(kDEF,0);
  ev.de_def	= hodode(  kDEF,0);
  ev.seg_def	= hodoseg( kDEF,0);

  ev.time_veto	= hodotime(kVeto,0);
  ev.de_veto	= hodode(  kVeto,0);
  ev.seg_veto	= hodoseg( kVeto,0);

  ev.time_btc	= hodotime(kBTC,0);
  ev.de_btc	= hodode(  kBTC,0);
  ev.seg_btc	= hodoseg( kBTC,0);

  return ev;
}

CDHEvent BasicAnalysis::get_cdhevent()
{
  CDHEvent ev;
  ev.tof_bhtt0=tof_bhtt0();
  ev.ctof_bhtt0=ctof_bhtt0();

  ev.mt_t0 =hodotime( kT0,0);
  ev.cmt_t0=hodoctime(kT0,0);
  ev.deu_t0=hododeu(  kT0,0);
  ev.ded_t0=hododed(  kT0,0);
  ev.seg_t0=hodoseg(  kT0,0);

  ev.mt_cdh =hodotime( kCDH,0);
  ev.cmt_cdh=hodoctime(kCDH,0);
  ev.deu_cdh=hododeu(  kCDH,0);
  ev.ded_cdh=hododed(  kCDH,0);
  ev.seg_cdh=hodoseg(  kCDH,0);
  ev.tu_cdh =hodotu(   kCDH,0);
  ev.td_cdh =hodotd(   kCDH,0);
  ev.ctu_cdh=hodoctu(  kCDH,0);
  ev.ctd_cdh=hodoctd(  kCDH,0);

  ev.mt_def =hodotime( kDEF,0);
  ev.cmt_def=hodoctime(kDEF,0);
  ev.deu_def=hododeu(  kDEF,0);
  ev.ded_def=hododed(  kDEF,0);
  ev.seg_def=hodoseg(  kDEF,0);

  return ev;
}
