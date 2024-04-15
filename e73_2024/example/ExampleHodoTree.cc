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
#include "HodoAnalyzer.hh"
#include "HodoRawHit.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "UserParamMan.hh"
#include "RawData.hh"
#include "VEvent.hh"
#include "HistTools.hh"
#include "UnpackerManager.hh"

#include "setup.hh"

#define TREE 1
#define DEBUG 0

namespace
{
using namespace e73_2024;
using namespace root;
using namespace hddaq::unpacker;
const auto& gUser = UserParamMan::GetInstance();
const auto& gUnpacker = GUnpacker::get_instance();
}

typedef std::vector<std::vector<Double_t>> TDCVector;
typedef std::vector<std::vector<Double_t>> TOTVector;
//_____________________________________________________________________________
namespace root
{
TTree *tree;
Int_t run_number;
Int_t event_number;

std::vector<Int_t>   ac_leading;

std::vector<Short_t> bht_segment;
TDCVector   bht_leading_top;
TDCVector   bht_leading_bottom;
TDCVector   bht_trailing_top;
TDCVector   bht_trailing_bottom;

TDCVector   t0new_leading_top;
TDCVector   t0new_leading_bottom;
std::vector<Short_t> t0new_qdc_top;
std::vector<Short_t> t0new_qdc_bottom;

std::vector<Short_t> t0_segment;
TDCVector   t0_leading_top;
TDCVector   t0_leading_bottom;
std::vector<Short_t> t0_qdc_top;
std::vector<Short_t> t0_qdc_bottom;

std::vector<Short_t> def_segment;
TDCVector   def_leading_top;
TDCVector   def_leading_bottom;
std::vector<Short_t> def_qdc_top;
std::vector<Short_t> def_qdc_bottom;

std::vector<Short_t> veto_segment;
TDCVector   veto_leading_left;
TDCVector   veto_leading_right;
std::vector<Short_t> veto_qdc_left;
std::vector<Short_t> veto_qdc_right;

std::vector<Short_t> btc_segment;
TDCVector   btc_leading_left;
TDCVector   btc_leading_right;
std::vector<Short_t> btc_qdc_left;
std::vector<Short_t> btc_qdc_right;

const Int_t nhodo=7;
Int_t khodo[nhodo]={kT0new,kBHT,kT0,kDEF,kVeto,kBTC,kRC};
// const Int_t nhodo=5;
// Int_t khodo[nhodo]={kT0new,kBHT,kT0,kDEF,kRC};

Double_t tdcbins[3]={5000,0,2e6};
Double_t adcbins[3]={4096,-0.5,4095.5};
Double_t mtdcbins[3]={2000,0,2000};
Double_t diffbins[3]={4000,-10000,10000};
Double_t debins[3]={5000,-2,48};
Double_t sumbins[3]={4096*4,-0.5,4096*4-0.5};
}

//_____________________________________________________________________________
void
InitializeEvent()
{
  ac_leading.clear();

  bht_segment.clear();
  bht_leading_top.clear();
  bht_leading_bottom.clear();
  bht_trailing_top.clear();
  bht_trailing_bottom.clear();

  t0new_leading_top.clear();
  t0new_leading_bottom.clear();
  t0new_qdc_top.clear();
  t0new_qdc_bottom.clear();

  t0_segment.clear();
  t0_leading_top.clear();
  t0_leading_bottom.clear();
  t0_qdc_top.clear();
  t0_qdc_bottom.clear();

  def_segment.clear();
  def_leading_top.clear();
  def_leading_bottom.clear();
  def_qdc_top.clear();
  def_qdc_bottom.clear();

  veto_segment.clear();
  veto_leading_left.clear();
  veto_leading_right.clear();
  veto_qdc_left.clear();
  veto_qdc_right.clear();

  btc_segment.clear();
  btc_leading_left.clear();
  btc_leading_right.clear();
  btc_qdc_left.clear();
  btc_qdc_right.clear();
}

//_____________________________________________________________________________
Bool_t
ProcessBegin()
{
  InitializeEvent();
  event_number=gUnpacker.get_event_number();
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
  HodoAnalyzer hodoAna(rawData);
  hodoAna.DecodeRawHits();
  MTDCAnalyzer MTDCAna(rawData);
  MTDCAna.DecodeRawHits();
  Bool_t BEAM  = MTDCAna.flag(kBeam);
  Bool_t KAON2 = MTDCAna.flag(kKaon2);
  Bool_t KAON3 = MTDCAna.flag(kKaon3);
  Bool_t KCDH1 = MTDCAna.flag(kKCDH1);

  Bool_t TOFK=false, TOFPi=false;
  Bool_t ACHIT=false;
  Bool_t VETO=false;

  // Time0
  Double_t time0=-9999;
  Double_t ctime0=-9999;
  {
    Int_t cid=DetIdT0new;
    Int_t nh = hodoAna.GetNHits(cid);
    for( Int_t i=0; i<nh; ++i ){
      Hodo2Hit *hit = hodoAna.GetHit(cid,i);
      if(!hit) continue;
      Int_t nind=hit->GetIndex();
      for(Int_t it=0;it<nind;it++){
	Double_t mt  = hit->MeanTime(it);
	if(gUser.Check("Time0",mt)){
	  time0=mt;
	  ctime0=hit->CMeanTime(it);
	}
      }
    }
  }
  // K/pi by BHD MeanTime
  {
    Int_t cid=DetIdBHT;
    Int_t nh = hodoAna.GetNHits(cid);
    for( Int_t i=0; i<nh; ++i ){
      Hodo2Hit *hit = hodoAna.GetHit(cid,i);
      if(!hit) continue;
      Int_t nind=hit->GetIndex();
      for(Int_t it=0;it<nind;it++){
	Double_t tof  = hit->MeanTime(it)-time0;
	if(gUser.Check("TOFK",tof)) TOFK=true;
	if(gUser.Check("TOFPi",tof)) TOFPi=true;
      }
    }
  }
  // AC
  {
    Int_t cid=DetIdAC;
    const HodoRHitContainer &cont = rawData.GetHodoRawHC(cid);
    Int_t nh = cont.size();
    TString tmpname="AC";
    for( Int_t i=0; i<nh; ++i ){
      HodoRawHit *raw = cont[i];
      if(!raw) continue;
      Int_t ntu=raw->GetSizeTdcUp();
      for(Int_t it=0;it<ntu;it++){
	Double_t tu  = raw->GetTdcUp(it);
	ac_leading.push_back(tu);
	if(gUser.Check("ACTDC",tu)) ACHIT=true;
	hist::H1(tmpname+"_TDC",tu,tdcbins);
      }
    }
  }

  std::vector<TString> trig_add;
  trig_add.push_back("");
  if(BEAM) trig_add.push_back("ifB");
  if(BEAM&&ACHIT) trig_add.push_back("ifPi");
  if(KAON2) trig_add.push_back("ifK2");
  if(KAON3) trig_add.push_back("ifK");
  if(KCDH1) trig_add.push_back("ifKCDH");
  Int_t ntrig=trig_add.size();

  // hodoscopes
  Int_t hodoseg[nhodo];
  Double_t hodotime[nhodo];
  for(Int_t ihodo=0;ihodo<nhodo;++ihodo){
    Int_t kHodo=khodo[ihodo];
    Int_t cid=hodoid[kHodo];
    TString name=hodoname[kHodo];
    Double_t mulbins[3]={nsegs[kHodo]+1,-0.5,nsegs[kHodo]+0.5};
    Double_t patbins[3]={nsegs[kHodo],-0.5,nsegs[kHodo]-0.5};
    Double_t mulbins2[3]={10,-0.5,9.5};

    // rawhit
    Int_t mulu=0,muld=0;
    const HodoRHitContainer &cont = rawData.GetHodoRawHC(cid);
    Int_t nh = cont.size();
    Int_t hitsegu=-1,hitsegd=-1;
    for( Int_t i=0; i<nh; ++i ){
      HodoRawHit *raw = cont[i];
      if(!raw) continue;
      Int_t seg = raw->SegmentId();
      TString segstr=Form("_seg%d",seg);
      Double_t au  = -1;
      Double_t ad  = -1;
      if(raw->GetSizeAdcUp()>0)      au=raw->GetAdcUp();
      if(raw->GetSizeAdcDown()>0)    ad=raw->GetAdcDown();
      Int_t ntu    =raw->GetSizeTdcUp();
      Int_t ntd    =raw->GetSizeTdcDown();
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
	hist::H1(name+"_ADCuwoT"+segstr,au,adcbins);
      }
      if(ntd>0){
	hist::H1(name+"_Patd",seg,patbins);
	hist::H1(name+"_ADCwTd"+segstr,ad,adcbins);
	muld++;
	hitsegd=seg;
      }else{
	hist::H1(name+"_ADCdwoT"+segstr,ad,adcbins);
      }
      std::vector<Double_t> tmp_tdcu;
      std::vector<Double_t> tmp_tdcd;
      for(Int_t it=0;it<ntu;it++){
	Double_t tu  = raw->GetTdcUp(it);
	hist::H1(name+"_TDCu"+segstr,tu,tdcbins);
	if(gUser.Check("HODOTDC",tu)) tmp_tdcu.push_back(tu);
	if(it<ntd){
	  Double_t td  = raw->GetTdcDown(it);
	  hist::H1(name+"_TDCdiffud",tu-td,diffbins);
	}
      }
      for(Int_t it=0;it<ntd;it++){
	Double_t td  = raw->GetTdcDown(it);
	if(gUser.Check("HODOTDC",td)) tmp_tdcd.push_back(td);
	hist::H1(name+"_TDCd"+segstr,td,tdcbins);
      }
      if(cid==DetIdBHT){
	std::vector<Double_t> tmp_trailingu;
	std::vector<Double_t> tmp_trailingd;
	for( Int_t it=0; it<raw->GetSizeAdcUp(); ++it ){
	  Double_t tt  = raw->GetAdcUp(it);
	  if(gUser.Check("HODOTDC",tt)) tmp_trailingu.push_back(tt);
	}
	for( Int_t it=0; it<raw->GetSizeAdcDown(); ++it ){
	  Double_t tt  = raw->GetAdcDown(it);
	  if(gUser.Check("HODOTDC",tt)) tmp_trailingd.push_back(tt);
	}
	bht_segment.push_back(seg);
	bht_leading_top.push_back(tmp_tdcu);
	bht_leading_bottom.push_back(tmp_tdcd);
	bht_trailing_top.push_back(tmp_trailingu);
	bht_trailing_bottom.push_back(tmp_trailingd);
      }else if(cid==DetIdT0new){
	t0new_leading_top.push_back(tmp_tdcu);
	t0new_leading_bottom.push_back(tmp_tdcd);
	t0new_qdc_top.push_back(au);
	t0new_qdc_bottom.push_back(ad);
      }else if(cid==DetIdT0){
	t0_segment.push_back(seg);
	t0_leading_top.push_back(tmp_tdcu);
	t0_leading_bottom.push_back(tmp_tdcd);
	t0_qdc_top.push_back(au);
	t0_qdc_bottom.push_back(ad);
      }else if(cid==DetIdDEF){
	def_segment.push_back(seg);
	def_leading_top.push_back(tmp_tdcu);
	def_leading_bottom.push_back(tmp_tdcd);
	def_qdc_top.push_back(au);
	def_qdc_bottom.push_back(ad);
      }else if(cid==DetIdVeto){
	veto_segment.push_back(seg);
	veto_leading_left.push_back(tmp_tdcu);
	veto_leading_right.push_back(tmp_tdcd);
	veto_qdc_left.push_back(au);
	veto_qdc_right.push_back(ad);
      }else if(cid==DetIdBTC){
	btc_segment.push_back(seg);
	btc_leading_left.push_back(tmp_tdcu);
	btc_leading_right.push_back(tmp_tdcd);
	btc_qdc_left.push_back(au);
	btc_qdc_right.push_back(ad);
      }
    }//for(i)
    if(mulu==1&&muld==1){
      hist::H2(name+"_UDcorr" ,hitsegu,hitsegd ,patbins,patbins);
    }
    hist::H1(name+"_Mulu",mulu,mulbins);
    hist::H1(name+"_Muld",muld,mulbins);
    // raw hit done

    // Decoded hit
    Int_t mul=0,mulgate=0;
    nh = hodoAna.GetNHits(cid);
    for( Int_t i=0; i<nh; ++i ){
      Hodo2Hit *hit = hodoAna.GetHit(cid,i);
      if(!hit) continue;
      Int_t seg = hit->SegmentId();
      TString segstr=Form("_seg%d",seg);
      Int_t nind= hit->GetIndex();
      hist::H1(name+"_Mul"+segstr,nind,mulbins2);
      hist::H2(name+"_Mul_pat",seg,nind,patbins,mulbins2);
      if(nind>0){
	hist::H1(name+"_Pat",seg,patbins);
	mul++;
      }
      Double_t de=-1,deu=-1,ded=-1,au=-1,ad=-1;
      if(cid!=DetIdBHT){
	de = hit->DeltaE();
	deu = hit->GetAUp();
	ded = hit->GetADown();
	au = hit->GetRawHit()->GetAdcUp();
	ad = hit->GetRawHit()->GetAdcDown();
	hist::H1(name+"_dE" +segstr,de ,debins,trig_add);
	hist::H1(name+"_dEu"+segstr,deu,debins,trig_add);
	hist::H1(name+"_dEd"+segstr,ded,debins,trig_add);
      }
      Bool_t TDCHIT=false;
      for(Int_t ii=0;ii<nind;ii++){
	Double_t mt  = hit->MeanTime(ii);
	Double_t cmt  =hit->CMeanTime(ii);
	//	std::cout<<ii<<"  "<<mt<<std::endl;
	Double_t tof = mt - time0;
	Double_t ctof = cmt - ctime0;
	hist::H1(name+"_MeanTime"        ,mt ,4000,-200,200, trig_add);
	hist::H1(name+"_MeanTime" +segstr,mt ,4000,-200,200, trig_add);
	hist::H1(name+"_CMeanTime"       ,cmt,4000,-200,200, trig_add);
	hist::H1(name+"_CMeanTime"+segstr,cmt,4000,-200,200, trig_add);
	if(time0>-9000){
	  hist::H1(name+"_TOF"        ,tof ,2000,-100,100, trig_add);
	  hist::H1(name+"_TOF" +segstr,tof ,2000,-100,100, trig_add);
	  hist::H1(name+"_cTOF"       ,ctof,2000,-100,100, trig_add);
	  hist::H1(name+"_cTOF"+segstr,ctof,2000,-100,100, trig_add);
	}
	if(TMath::Abs(mt)>10) continue;
	TDCHIT=true;
	mulgate++;
	hodoseg[ihodo]=seg;
	hodotime[ihodo]=mt;
      } // nindex
      if(BEAM&&ACHIT){
	if(cid!=DetIdDEF&&!TDCHIT) continue;
	if(cid==DetIdBHT) continue;
	hist::H1(name+"_ADCuwB"+segstr,au,adcbins);
	hist::H1(name+"_ADCdwB"+segstr,ad,adcbins);
	hist::H2("Run_"+name+"_dEuwB"+segstr,run_number,deu,400,100,900,400,-1,19);
	hist::H2("Run_"+name+"_dEdwB"+segstr,run_number,ded,400,100,900,400,-1,19);
      }
    }//for(ihit)
    hist::H1(name+"_Mul",mul ,mulbins);
    hist::H1(name+"_Mulgate",mulgate,mulbins);
    // decoded hit done
  }
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  tree=new TTree("tree","vft raw data");
  tree->Branch("runnum", &run_number);
  tree->Branch("evnum", &event_number);
  tree->Branch("ac_leading",    &ac_leading);

  tree->Branch("bht_segment",        &bht_segment);
  tree->Branch("bht_leading_top",    &bht_leading_top);
  tree->Branch("bht_leading_bottom", &bht_leading_bottom);
  tree->Branch("bht_trailing_top",   &bht_trailing_top);
  tree->Branch("bht_trailing_bottom",&bht_trailing_bottom);

  tree->Branch("t0_segment",        &t0_segment);
  tree->Branch("t0_leading_top",    &t0_leading_top);
  tree->Branch("t0_leading_bottom", &t0_leading_bottom);
  tree->Branch("t0_qdc_top",        &t0_qdc_top);
  tree->Branch("t0_qdc_bottom",     &t0_qdc_bottom);

  tree->Branch("t0new_leading_top",    &t0new_leading_top);
  tree->Branch("t0new_leading_bottom", &t0new_leading_bottom);
  tree->Branch("t0new_qdc_top",        &t0new_qdc_top);
  tree->Branch("t0new_qdc_bottom",     &t0new_qdc_bottom);

  tree->Branch("def_segment",        &def_segment);
  tree->Branch("def_leading_top",    &def_leading_top);
  tree->Branch("def_leading_bottom", &def_leading_bottom);
  tree->Branch("def_qdc_top",        &def_qdc_top);
  tree->Branch("def_qdc_bottom",     &def_qdc_bottom);

  tree->Branch("veto_segment",       &veto_segment);
  tree->Branch("veto_leading_left",  &veto_leading_left);
  tree->Branch("veto_leading_right", &veto_leading_right);
  tree->Branch("veto_qdc_left",      &veto_qdc_left);
  tree->Branch("veto_qdc_right",     &veto_qdc_right);

  tree->Branch("btc_segment",       &btc_segment);
  tree->Branch("btc_leading_left",  &btc_leading_left);
  tree->Branch("btc_leading_right", &btc_leading_right);
  tree->Branch("btc_qdc_left",      &btc_qdc_left);
  tree->Branch("btc_qdc_right",     &btc_qdc_right);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    ( InitializeParameter<HodoParamMan>("HDPRM","CDSPRM") ) &&
    ( InitializeParameter<HodoPHCMan>  ("HDPHC","CDSPHC") ) &&
    ( InitializeParameter<UserParamMan>("USER") ) ;
}

//_____________________________________________________________________________
Bool_t
ConfMan::BeginRunProcess()
{
  run_number = get_run_number();
  return true;
}
//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
