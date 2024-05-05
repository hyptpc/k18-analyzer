// -*- C++ -*-

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#include <TString.h>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "EventAnalyzer.hh"
#include "RootHelper.hh"
#include "DCAnalyzer.hh"
#include "DCCluster.hh"
#include "DCHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoRawHit.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "UserParamMan.hh"
#include "XTMapMan.hh"
#include "DCTdcCalibMan.hh"
#include "BLDCWireMapMan.hh"
#include "RawData.hh"
#include "HistTools.hh"
#include "UnpackerManager.hh"

#define DEBUG 0
namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
// UInt_t run_number;
UInt_t event_number;
}

//_____________________________________________________________________________
Bool_t
ProcessBegin()
{
  event_number = gUnpacker.get_event_number();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessNormal()
{
  using root::HF1;

  RawData rawData;
  rawData.DecodeHits();

  HodoAnalyzer hodoAna(rawData);
  hodoAna.DecodeHits<FiberHit>("BHT");
  hodoAna.TimeCut("BHT");

  EventAnalyzer evAna;

  HF1("Status", 0);
  evAna.TriggerFlag(rawData);

  HF1("Status", 1);

  // BeamFlag
  beam::EBeamFlag beam_flag = evAna.BeamFlag(rawData);

  evAna.HodoRawHit(rawData);
  evAna.HodoRawHit(rawData, beam_flag);

  evAna.HodoHit(hodoAna);
  evAna.HodoHit(hodoAna, beam_flag);

#if 0
  Bool_t BEAM  = MTDCAna.flag(kBeam);
  Bool_t KAON2 = MTDCAna.flag(kKaon2);
  Bool_t KAON3 = MTDCAna.flag(kKaon3);
  // Bool_t KCDH1 = MTDCAna.flag(kKCDH1);
  //  Bool_t DEUTERON = MTDCAna.flag(kDeteuteron);
  Bool_t DEUTERON = MTDCAna.flag(15);

  Bool_t TOFK=false, TOFPi=false, TOFP=false, TOFD=false;
  Bool_t ACHIT=false;
  Bool_t VETO=false;

  // Time0
  Double_t time0 = TMath::QuietNaN();
  Double_t ctime0= TMath::QuietNaN();
  {
    Int_t cid = DetIdT0new;
    for(Int_t i=0, n=hodoAna.GetNHits(cid); i<n; ++i){
      auto hit = hodoAna.GetHit(cid, i);
      if(!hit) continue;
      Int_t nind = hit->GetIndex();
      for(Int_t it=0; it<nind; ++it){
	Double_t mt = hit->MeanTime(it);
	if(gUser.IsInRange("Time0", mt)){
	  time0 = mt;
          std::cout << "time0 : " << time0 << " -> " << mt << std::endl;
	  ctime0 = hit->CMeanTime(it);
	}
      }
    }
  }

  // K/pi by BHD MeanTime
  {
    Int_t cid = DetIdBHT;
    Int_t nh = hodoAna.GetNHits(cid);
    for(Int_t i=0; i<nh; ++i){
      auto hit = hodoAna.GetHit(cid, i);
      if(!hit) continue;
      Int_t nind = hit->GetIndex();
      for(Int_t it=0; it<nind; ++it){
	Double_t tof = hit->MeanTime(it)-time0;
	if(gUser.IsInRange("TOFK",tof)) TOFK=true;
	if(gUser.IsInRange("TOFPi",tof)) TOFPi=true;
	if(gUser.IsInRange("TOFP", tof)) TOFP=true;
	if(gUser.IsInRange("TOFD", tof)) TOFD=true;
      }
    }
  }

  // AC
  {
    Int_t cid = DetIdAC;
    const auto& cont = rawData.GetHodoRawHC(cid);
    Int_t nh = cont.size();
    TString tmpname="AC";
    for( Int_t i=0; i<nh; ++i ){
      HodoRawHit *raw = cont[i];
      if(!raw) continue;
      Int_t ntu=raw->GetSizeTdcUp();
      for(Int_t it=0;it<ntu;it++){
	Double_t tu  = raw->GetTdcUp(it);
	if(gUser.IsInRange("AC_TDC",tu)) ACHIT=true;
	hist::H1(tmpname+"_TDC",tu,tdcbins);
      }
    }
  }

  std::vector<TString> trig_add;
  trig_add.push_back("");
  if(BEAM) trig_add.push_back("ifB");
  if(BEAM&&ACHIT) trig_add.push_back("ifPi");
  if(TOFD&&!ACHIT)  trig_add.push_back("ifD");
  if(TOFP&&!ACHIT)  trig_add.push_back("ifP");
  if(KAON2) trig_add.push_back("ifK2");
  //  if(KAON3) trig_add.push_back("ifK");
  //  if(KCDH1) trig_add.push_back("ifKCDH");
  Int_t ntrig=trig_add.size();

  HF1("Status", 4);

  // hodoscopes
  Int_t hodoseg[kNumHodo];
  Double_t hodotime[kNumHodo];
  for(Int_t ihodo=0;ihodo<kNumHodo;++ihodo){
    Int_t cid=DetIdHodo[ihodo];
    TString name=NameHodo[ihodo];
    Double_t mulbins[3]={NumOfSegHodo[ihodo]+1,-0.5,NumOfSegHodo[ihodo]+0.5};
    Double_t patbins[3]={NumOfSegHodo[ihodo],-0.5,NumOfSegHodo[ihodo]-0.5};
    Double_t mulbins2[3]={10,-0.5,9.5};

    // rawhit
    Int_t mulu=0,muld=0;
    const HodoRHitContainer &cont = rawData.GetHodoRawHC(cid);
    Int_t nh = cont.size();
    std::vector<Int_t> hitsegu;
    std::vector<Int_t> hitsegd;
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
      Int_t ngateu=0;
      Int_t ngated=0;
      for(Int_t it=0;it<ntu;it++){
	Double_t tu  = raw->GetTdcUp(it);
	hist::H1(name+"_TDCu"+segstr,tu,tdcbins);
	if(it+1<ntu)
	  hist::H1(name+"_TDCdiffu"+segstr,raw->GetTdcUp(it+1)-tu,tdcbins);
	if(it<ntd){
	  Double_t td  = raw->GetTdcDown(it);
	  hist::H1(name+"_TDCdiffud",tu-td,diffbins);
	}
	if(gUser.IsInRange("HODOTDC",tu)){
	  ngateu++;
	}
      }
      for(Int_t it=0;it<ntd;it++){
	Double_t td  = raw->GetTdcDown(it);
	hist::H1(name+"_TDCd"+segstr,td,tdcbins);
	if(it+1<ntd)
	  hist::H1(name+"_TDCdiffd"+segstr,raw->GetTdcDown(it+1)-td,tdcbins);
	if(gUser.IsInRange("HODOTDC",td)){
	  ngated++;
	}
      }
      if(cid==DetIdBHT){
	Double_t totbins[3]={2000,0,1e5};
	Double_t segtotbins[6]={64,-0.5,63.5,200,0,1e5};
	Double_t dttotbins[6]={200,1e6,1.4e6,200,0,1e5};
	Int_t nau = raw->GetSizeAdcUp();
	Int_t nad = raw->GetSizeAdcDown();
	for(Int_t it=0;it<ntu;it++){
	  Double_t tu  = raw->GetTdcUp(it);
	  if(it<nau){
	    Double_t au  = raw->GetAdcUp(it);
	    hist::H1(name+"_TOTu",tu-au,totbins);
	    hist::H1(name+"_TOTu"+segstr,tu-au,totbins);
	    hist::H2(name+"_TDC_TOTu",tu,tu-au,dttotbins);
	    hist::H2(name+"_Seg_TOTu",seg,tu-au,segtotbins);
	  }
	}
	for(Int_t it=0;it<ntd;it++){
	  Double_t td  = raw->GetTdcDown(it);
	  if(it<nad){
	    Double_t ad  = raw->GetAdcDown(it);
	    hist::H1(name+"_TOTd",td-ad,totbins);
	    hist::H1(name+"_TOTd"+segstr,td-ad,totbins);
	    hist::H2(name+"_TDC_TOTd",td,td-ad,dttotbins);
	    hist::H2(name+"_Seg_TOTd",seg,td-ad,segtotbins);
	  }
	}
      }
      hist::H1(name+"_ADCu"+segstr,au,adcbins);
      hist::H1(name+"_ADCd"+segstr,ad,adcbins);
      hist::H1(name+"_Mulu"+segstr,ntu,mulbins2);
      hist::H1(name+"_Muld"+segstr,ntd,mulbins2);
      hist::H2(name+"_Mulu_Muld"+segstr,ntu,ntd,mulbins2,mulbins2);
      hist::H2(name+"_Mulu_pat",seg,ntu,patbins,mulbins2);
      hist::H2(name+"_Muld_pat",seg,ntd,patbins,mulbins2);
      hist::H1(name+"_Mulgateu"+segstr,ngateu,mulbins2);
      hist::H1(name+"_Mulgated"+segstr,ngated,mulbins2);
      hist::H2(name+"_Mulgateu_Mulgated"+segstr,ngateu,ngated,mulbins2,mulbins2);
      hist::H2(name+"_Mulgateu_pat",seg,ngateu,patbins,mulbins2);
      hist::H2(name+"_Mulgated_pat",seg,ngated,patbins,mulbins2);
      if(ngateu>0){
	hist::H1(name+"_Patu",seg,patbins);
	hist::H1(name+"_ADCwTu"+segstr,au,adcbins);
	mulu++;
	hitsegu.push_back(seg);
      }else{
	hist::H1(name+"_ADCwoTu"+segstr,au,adcbins);
      }
      if(ngated>0){
	hist::H1(name+"_Patd",seg,patbins);
	hist::H1(name+"_ADCwTd"+segstr,ad,adcbins);
	muld++;
	hitsegd.push_back(seg);
      }else{
	hist::H1(name+"_ADCwoTd"+segstr,ad,adcbins);
      }
    }//for(i)
    for(Int_t iu=0;iu<hitsegu.size();iu++){
      for(Int_t id=0;id<hitsegd.size();id++){
	hist::H2(name+"_UDcorr" ,hitsegu.at(iu),hitsegd.at(id) ,patbins,patbins);
	hist::H2(name+"_EvNum_UDdiff",event_number,hitsegu.at(iu)-hitsegd.at(id),evudbins);
	if(name=="BHT"){
	  if(hitsegu.at(iu)<32&&hitsegd.at(id)<32){
	    hist::H2(name+"_EvNum_UDdiff_left",event_number,hitsegu.at(iu)-hitsegd.at(id),evudbins);
	  }else if(hitsegu.at(iu)<64&&hitsegd.at(id)<64){
	    hist::H2(name+"_EvNum_UDdiff_right",event_number,hitsegu.at(iu)-hitsegd.at(id),evudbins);
	  }
	}
      }
    }
    hist::H1(name+"_Mulu",mulu,mulbins);
    hist::H1(name+"_Muld",muld,mulbins);
    // raw hit done
    // Decoded hit
    Int_t mul=0,mulgate=0;
    nh = hodoAna.GetNHits(cid);
    Bool_t TDCHIT=false;
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
      Int_t idx=-1;
      for(Int_t ii=0;ii<nind;ii++){
	Double_t mt  = hit->MeanTime(ii);
	Double_t cmt  =hit->CMeanTime(ii);
	// if(!gUser.IsInRange("BHTTOT",hit->GetAUp(ii)))   continue;
	// if(!gUser.IsInRange("BHTTOT",hit->GetADown(ii))) continue;
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
	  Double_t de = hit->DeltaE();
	  Double_t demin=-1,demax=19;
	  if(cid==DetIdBHT){
	    de=hit->DeltaCE(ii);
	    demin=0;
	    demax=4e4;
	  }
	  hist::H2(name+"_TOF_dE",tof , de, 200,-50,50, 100, demin, demax, trig_add);
	}
	if(TMath::Abs(mt)>20) continue;
	hist::H2(name+"_EvNum_MeanTime" ,event_number,mt ,evtimebins);
	hist::H2(name+"_EvNum_CMeanTime",event_number,cmt,evtimebins);
	hist::H2(name+"_EvNum_TOF"      ,event_number,tof ,evtimebins);
	hist::H2(name+"_EvNum_cTOF"     ,event_number,ctof,evtimebins);
	TDCHIT=true;
	idx=ii;
	mulgate++;
	hodoseg[ihodo]=seg;
	hodotime[ihodo]=mt;
      } // nindex
      if(cid==DetIdBHT && idx<0) continue;
      Double_t demin=-1,demax=49,demax2=19;
      Double_t au = hit->GetRawHit()->GetAdcUp();
      Double_t ad = hit->GetRawHit()->GetAdcDown();
      Double_t de = hit->DeltaE();
      Double_t deu = hit->GetAUp();
      Double_t ded = hit->GetADown();
      if(cid==DetIdBHT){
	demin=0;demax=1e5;demax2=2e4;
	// au = hit->GetAUp(idx);
	// ad = hit->GetADown(idx);
	// de = hit->DeltaCE(idx);
	// deu = hit->GetCAUp(idx);
	// ded = hit->GetCADown(idx);
      }
      Double_t tmpdebins[3]={1000,demin,demax};
      Double_t tmpdeudbins[6]={100,demin,demax2,100,demin,demax2};

      hist::H1(name+"_dE" +segstr,de ,tmpdebins,trig_add);
      hist::H2(name+"_dEud"+segstr,deu,ded,tmpdeudbins,trig_add);
      if(TDCHIT){
	hist::H2(name+"_EvNum_dE",event_number,de,evdebins);
      }
      if(BEAM&&ACHIT){
	if(cid!=DetIdDEF&&!TDCHIT) continue;
	hist::H1(name+"_ADCuwB"+segstr,au,adcbins);
	hist::H1(name+"_ADCdwB"+segstr,ad,adcbins);
	hist::H1(name+"_ADCuwB",au,adcbins);
	hist::H1(name+"_ADCdwB",ad,adcbins);
	hist::H1(name+"_dEuwB"+segstr,deu,tmpdebins);
	hist::H1(name+"_dEdwB"+segstr,ded,tmpdebins);
	hist::H1(name+"_dEuwB",deu,tmpdebins);
	hist::H1(name+"_dEdwB",ded,tmpdebins);
	hist::H2("Run_"+name+"_dEuwB",run_number,deu,400,100,900,100,demin,demax2);
	hist::H2("Run_"+name+"_dEdwB",run_number,ded,400,100,900,100,demin,demax2);
	if(cid==DetIdBHT || cid==DetIdCDH) continue;
	hist::H2("Run_"+name+"_dEuwB"+segstr,run_number,deu,400,100,900,100,demin,demax2);
	hist::H2("Run_"+name+"_dEdwB"+segstr,run_number,ded,400,100,900,100,demin,demax2);
      }
    }//for(ihit)
    hist::H1(name+"_Mul",mul ,mulbins);
    hist::H1(name+"_Mulgate",mulgate,mulbins);
    hist::H2(name+"_EvNum_Mul"    ,event_number,mul,evmulbins);
    hist::H2(name+"_EvNum_Mulgate",event_number,mulgate,evmulbins);
    // decoded hit done
  }
  if(hodoseg[kVeto]) VETO=true;

#endif

  HF1("Status", 20);

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
  Bool_t beam_flag = true;
  hist::BuildStatus();
  hist::BuildTriggerFlag();
  hist::BuildHodoRaw(beam_flag);
  hist::BuildHodoHit(beam_flag);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<HodoParamMan>("HDPRM")) &&
    (InitializeParameter<HodoPHCMan>("HDPHC")) &&
    (InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
