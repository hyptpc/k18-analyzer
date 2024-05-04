// -*- C++ -*-

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#include <TString.h>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RootHelper.hh"
#include "MTDCAnalyzer.hh"
#include "MTDCRawHit.hh"
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
using namespace root;
using namespace hddaq::unpacker;
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
UInt_t run_number;
UInt_t event_number;
}

//_____________________________________________________________________________
namespace root
{
TH1   *h[MaxHist];
Double_t tdcbins[3]={5000,0,2e6};
Double_t adcbins[3]={4096,-0.5,4095.5};
Double_t mtdcbins[3]={2000,0,2000};
Double_t diffbins[3]={2000,-1e5,1e5};
Double_t debins[3]={5000,-2,48};
Double_t debins2[3]={1000,-100,900};
Double_t deudbins[6]={100,-2,18,100,-2,18};
Double_t sumbins[3]={4096*4,-0.5,4096*4-0.5};
Double_t evmulbins[6]={200,0,5e6,11,-0.5,10.5};
Double_t evudbins[6]={200,0,5e6,21,-10.5,10.5};
Double_t evdebins[6]={200,0,5e6,200,-2,18};
Double_t evdebins2[6]={200,0,5e6,100,-100,900};
Double_t evadcbins[6]={200,0,5e6,100,0,500};
Double_t evtimebins[6]={200,0,5e6,200,-10,10};
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
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif

  RawData rawData;
  rawData.DecodeHits();

  HodoAnalyzer hodoAna(rawData);
  hodoAna.DecodeRawHits();

  MTDCAnalyzer MTDCAna(rawData);
  MTDCAna.DecodeRawHits();

  root::HF1("Status", 0);

  // Trigger flag
  {
    static const Char_t* name = "TriggerFlag";
    const auto& cont = rawData.GetMTDCRawHC(DetIdTrigFlag);
    for(Int_t i=0, n=cont.size(); i<n; ++i){
      auto raw = cont[i];
      auto ntu = raw->GetSizeLeading();
      auto seg = raw->SegmentId();
      for(Int_t it=0; it<ntu; ++it){
	Double_t tu = raw->GetLeading(it);
	root::HF1(Form("%s_TDC_seg%d", name, seg), tu);
      }
      if(ntu>0){
	root::HF1(Form("%s_HitPat", name), seg);
      }
    }
  }

  root::HF1("Status", 1);

  // AC
  Bool_t ac_hit = false;
  {
    static const Char_t* name = "AC";
    const auto& cont = rawData.GetHodoRawHC(DetIdAC);
    for(Int_t i=0, n=cont.size(); i<n; ++i){
      auto raw = cont[i];
      if(!raw) continue;
      for(Int_t j=0, m=raw->GetSizeTdcUp(); j<m; ++j){
	Double_t t = raw->GetTdcUp(j);
        if(gUser.IsInRange(Form("%s_TDC", name), t)){
          ac_hit = true; break;
        }
      }
    }
  }

  // ParticleFlag
  Bool_t is_pion = ac_hit;
  Bool_t is_k = !ac_hit;
  Bool_t is_p = !ac_hit;

  // BHT
  for(const auto& particle: std::vector<TString>{"", "_Pi", "_K", "_P"}){
    if(particle == "_Pi" && !is_pion) continue;
    if(particle == "_K" && !is_k) continue;
    if(particle == "_P" && !is_p) continue;
    const Char_t* p = particle.Data();
    static const Char_t* name = "BHT";
    const auto& cont = rawData.GetHodoRawHC(DetIdBHT);
    Int_t multi_or = 0;
    Int_t multi_and = 0;
    for(Int_t i=0, n=cont.size(); i<n; ++i){
      auto raw = cont[i];
      if(!raw) continue;
      auto seg = raw->SegmentId();
      // Up
      Bool_t u_in_range = false;
      for(Int_t j=0, m=raw->GetSizeTdcUp(); j<m; ++j){
	Double_t t = raw->GetTdcUp(j);
        if(gUser.IsInRange(Form("%s_TDC", name), t))
          u_in_range = true;
	root::HF1(Form("%s_TDC_seg%dU%s", name, seg, p), t);
      }
      for(Int_t j=0, m=raw->GetSizeAdcUp(); j<m; ++j){
	Double_t t = raw->GetAdcUp(j); // trailing
	root::HF1(Form("%s_Trailing_seg%dU%s", name, seg, p), t);
	if(m == raw->GetSizeTdcUp()){
	  Double_t l = raw->GetTdcUp(j);
          if(gUser.IsInRange(Form("%s_TDC", name), l)){
            root::HF1(Form("%s_TOT_seg%dU%s", name, seg, p), l - t);
          }
	}
      }
      // Down
      Bool_t d_in_range = false;
      for(Int_t j=0, m=raw->GetSizeTdcDown(); j<m; ++j){
	Double_t t = raw->GetTdcDown(j);
        if(gUser.IsInRange(Form("%s_TDC", name), t))
          d_in_range = true;
	root::HF1(Form("%s_TDC_seg%dD%s", name, seg, p), t);
      }
      for(Int_t j=0, m=raw->GetSizeAdcDown(); j<m; ++j){
	Double_t t = raw->GetAdcDown(j); // trailing
	root::HF1(Form("%s_Trailing_seg%dD%s", name, seg, p), t);
	if(m == raw->GetSizeTdcDown()){
	  Double_t l = raw->GetTdcDown(j);
          if(gUser.IsInRange(Form("%s_TDC", name), l)){
            root::HF1(Form("%s_TOT_seg%dD%s", name, seg, p), l - t);
          }
	}
      }
      if(u_in_range || d_in_range){
        root::HF1(Form("%s_HitPat_OR%s", name, p), seg);
        ++multi_or;
      }
      if(u_in_range && d_in_range){
        root::HF1(Form("%s_HitPat_AND%s", name, p), seg);
        ++multi_and;
      }
    }
    root::HF1(Form("%s_Multi_OR%s", name, p), multi_or);
    root::HF1(Form("%s_Multi_AND%s", name, p), multi_and);
  }

  root::HF1("Status", 2);

  // AC
  for(const auto& particle: std::vector<TString>{"", "_Pi", "_K", "_P"}){
    if(particle == "_Pi" && !is_pion) continue;
    if(particle == "_K" && !is_k) continue;
    if(particle == "_P" && !is_p) continue;
    const Char_t* p = particle.Data();
    static const Char_t* name = "AC";
    const auto& cont = rawData.GetHodoRawHC(DetIdAC);
    Int_t multi = 0;
    for(Int_t i=0, n=cont.size(); i<n; ++i){
      auto raw = cont[i];
      if(!raw) continue;
      Int_t seg = raw->SegmentId();
      Bool_t is_good = false;
      for(Int_t j=0, m=raw->GetSizeTdcUp(); j<m; ++j){
	Double_t t = raw->GetTdcUp(j);
        if(gUser.IsInRange(Form("%s_TDC", name), t))
          is_good = true;
	root::HF1(Form("%s_TDC_seg%d%s", name, seg, p), t);
      }
      Double_t suma = 0;
      for(Int_t j=0, m=raw->GetSizeAdcUp(); j<m; ++j){
	Double_t a = raw->GetAdcUp(j);
        suma += a;
	root::HF1(Form("%s_ADC_seg%d%s", name, j+1, p), a);
        if(is_good)
          root::HF1(Form("%s_AWT_seg%d%s", name, j+1, p), a);
      }
      // SUM
      root::HF1(Form("%s_ADC_seg%d%s", name, 0, p), suma);
      if(is_good){
        root::HF1(Form("%s_AWT_seg%d%s", name, 0, p), suma);
        root::HF1(Form("%s_HitPat%s", name, p), seg);
        multi++;
      }
    }
    root::HF1(Form("%s_Multi%s", name, p), multi);
  }

  root::HF1("Status", 3);

  // Hodoscope
  for(const auto& particle: std::vector<TString>{"", "_Pi", "_K", "_P"}){
    if(particle == "_Pi" && !is_pion) continue;
    if(particle == "_K" && !is_k) continue;
    if(particle == "_P" && !is_p) continue;
    const Char_t* p = particle.Data();
    for(Int_t ihodo=kT1; ihodo<kNumHodo;++ihodo){
      const Char_t* name = NameHodo[ihodo];
      Int_t multi_or = 0;
      Int_t multi_and = 0;
      const auto& cont = rawData.GetHodoRawHC(DetIdHodo[ihodo]);
      for(Int_t i=0, n=cont.size(); i<n; ++i){
        auto raw = cont[i];
        if(!raw) continue;
        auto seg = raw->SegmentId();
        // Up
        Bool_t u_in_range = false;
        for(Int_t j=0, m=raw->GetSizeTdcUp(); j<m; ++j){
          Double_t t = raw->GetTdcUp(j);
          if(gUser.IsInRange(Form("%s_TDC", name), t))
            u_in_range = true;
          root::HF1(Form("%s_TDC_seg%dU%s", name, seg, p), t);
        }
        for(Int_t j=0, m=raw->GetSizeAdcUp(); j<m; ++j){
          Double_t a = raw->GetAdcUp(j);
          root::HF1(Form("%s_ADC_seg%dU%s", name, seg, p), a);
          if(u_in_range)
            root::HF1(Form("%s_AWT_seg%dU%s", name, seg, p), a);
        }
        // Down
        Bool_t d_in_range = false;
        for(Int_t j=0, m=raw->GetSizeTdcDown(); j<m; ++j){
          Double_t t = raw->GetTdcDown(j);
          if(gUser.IsInRange(Form("%s_TDC", name), t))
            d_in_range = true;
          root::HF1(Form("%s_TDC_seg%dD%s", name, seg, p), t);
        }
        for(Int_t j=0, m=raw->GetSizeAdcDown(); j<m; ++j){
          Double_t a = raw->GetAdcDown(j);
          root::HF1(Form("%s_ADC_seg%dD%s", name, seg, p), a);
          if(d_in_range)
            root::HF1(Form("%s_AWT_seg%dD%s", name, seg, p), a);
        }
        if(u_in_range || d_in_range){
          root::HF1(Form("%s_HitPat_OR%s", name, p), seg);
          ++multi_or;
        }
        if(u_in_range && d_in_range){
          root::HF1(Form("%s_HitPat_AND%s", name, p), seg);
          ++multi_and;
        }
      }
      root::HF1(Form("%s_Multi_OR%s", name, p), multi_or);
      root::HF1(Form("%s_Multi_AND%s", name, p), multi_and);
    }
  }

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

  root::HF1("Status", 4);

#if 0
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

  root::HF1("Status", 20);

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
  Double_t hrtdcbins1[3] = {45000, 950000, 1400000};
  Double_t hrtdcbins2[3] = {50000, 0, 1000000}; // for CVC, NC
  Double_t hrtotbins[3] = {5000, 0, 50000};
  Double_t adcbins[3] = {4096, -0.5, 4095.5};
  Double_t mhtdcbins[3] = {2000, 0, 2000};
  Double_t mhtotbins[3] = {1000, 0, 1000};

  root::HB1("Status", 21, -0.5, 20.5);

  { // TriggerFlag
    const Char_t* name = "TriggerFlag";
    Double_t patbins[3]={NumOfSegTrigFlag, -0.5, NumOfSegTrigFlag-0.5};
    for(Int_t i=0; i<NumOfSegTrigFlag; ++i){
      root::HB1(Form("%s_TDC_seg%d", name, i), mhtdcbins);
    }
    root::HB1(Form("%s_HitPat; Segment; Counts", name), patbins);
  }

  for(const auto& particle: std::vector<TString>{"", "_Pi", "_K", "_P"}){
    const Char_t* p = particle.Data();
    { // AC
      const Char_t* name = "AC";
      Int_t nseg = NumOfSegAC;
      for(Int_t i=0; i<nseg; ++i){
        root::HB1(Form("%s_ADC_seg%d%s", name, i, p), adcbins);
        root::HB1(Form("%s_AWT_seg%d%s", name, i, p), adcbins);
        root::HB1(Form("%s_TDC_seg%d%s", name, i, p), mhtdcbins);
      }
      root::HB1(Form("%s_HitPat%s", name, p), nseg, -0.5, nseg - 0.5);
      root::HB1(Form("%s_Multi%s", name, p), nseg + 1, -0.5, nseg + 0.5);
    }
    { // BHT
      const Char_t* name = "BHT";
      Int_t nseg = NumOfSegBHT;
      for(Int_t i=0; i<nseg; ++i){
        for(const auto& uord: std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          root::HB1(Form("%s_TDC_seg%d%s%s", name, i, ud, p), hrtdcbins1);
          root::HB1(Form("%s_Trailing_seg%d%s%s", name, i, ud, p), hrtdcbins1);
          root::HB1(Form("%s_TOT_seg%d%s%s", name, i, ud, p), hrtotbins);
        }
      }
      for(const auto& uord: std::vector<TString>{"OR", "AND"} ){
        auto ud = uord.Data();
        root::HB1(Form("%s_HitPat_%s%s", name, ud, p), nseg, -0.5, nseg - 0.5);
        root::HB1(Form("%s_Multi_%s%s", name, ud, p), nseg + 1, -0.5, nseg + 0.5);
      }
    }
    // Hodoscope
    for(Int_t ihodo=kT1; ihodo<kNumHodo;++ihodo){
      auto name = NameHodo[ihodo].Data();
      Double_t* hrtdcbins;
      if(NameHodo[ihodo].Contains("CVC") ||
         NameHodo[ihodo].Contains("NC")){
        hrtdcbins = hrtdcbins2;
      }else{
        hrtdcbins = hrtdcbins1;
      }
      Int_t nseg = NumOfSegHodo[ihodo];
      for(const auto& uord: std::vector<TString>{"U", "D"} ){
        auto ud = uord.Data();
        for(Int_t i=0; i<nseg; ++i){
          root::HB1(Form("%s_ADC_seg%d%s%s", name, i, ud, p), adcbins);
          root::HB1(Form("%s_AWT_seg%d%s%s", name, i, ud, p), adcbins);
          root::HB1(Form("%s_TDC_seg%d%s%s", name, i, ud, p), hrtdcbins);
        }
      }
      for(const auto& uord: std::vector<TString>{"OR", "AND"} ){
        auto ud = uord.Data();
        root::HB1(Form("%s_HitPat_%s%s", name, ud, p), nseg, -0.5, nseg - 0.5);
        root::HB1(Form("%s_Multi_%s%s", name, ud, p), nseg + 1, -0.5, nseg + 0.5);
      }
    }
  }

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
