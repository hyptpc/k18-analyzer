// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DetectorID.hh"
#include "HodoAnalyzer.hh"
#include "HodoHit.hh"
#include "RMAnalyzer.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "S2sLib.hh"
#include "RawData.hh"
#include "UnpackerManager.hh"

#define HodoCut 0
#define TdcCut  1
#define TotCut  0
#define Chi2Cut 0

namespace
{
using namespace root;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& zK18Target = gGeom.LocalZ("K18Target");
// const Double_t SdcInXoffset = 49.95; // mm
const Double_t SdcInXoffset = 0.; // mm
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  Int_t nhBh1;
  Double_t tBh1[MaxHits];
  Double_t deBh1[MaxHits];

  Int_t nhBh2;
  Double_t tBh2[MaxHits];
  Double_t deBh2[MaxHits];
  Double_t Bh2Seg[MaxHits];

  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;

  Int_t ntrack;
  Int_t nlayer[MaxHits];
  Double_t chisqr[MaxHits];
  Double_t x0[MaxHits];
  Double_t y0[MaxHits];
  Double_t u0[MaxHits];
  Double_t v0[MaxHits];
  Double_t pos[NumOfLayersBcOut+NumOfLayersSdcIn];
  Double_t res[NumOfLayersBcOut+NumOfLayersSdcIn];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum     =  0;
  ntrack    =  0;
  nhBh2     =  0;
  nhBh1     =  0;
  Time0Seg  = -1;
  deTime0   = -1;
  Time0     = qnan;
  CTime0    = qnan;

  for(Int_t it=0; it<MaxHits; it++){
    tBh1[it]   = qnan;
    deBh1[it]  = qnan;
    tBh2[it]   = qnan;
    deBh2[it]  = qnan;
    Bh2Seg[it] = -1;
    nlayer[it] =  0;
    chisqr[it] = qnan;
    x0[it]     = qnan;
    y0[it]     = qnan;
    u0[it]     = qnan;
    v0[it]     = qnan;
  }

  for(Int_t it=0; it<NumOfSegTrig; it++){
    trigpat[it] = -1;
    trigflag[it] = -1;
  }

  for (Int_t it=0; it<NumOfLayersBcOut+NumOfLayersSdcIn; it++){
      pos[it] = qnan;
      res[it] = qnan;
  }
}

//_____________________________________________________________________________
namespace root
{
Event event;
TH1   *h[MaxHist];
TTree *tree;
enum eParticle { kKaon, kPion, nParticle };
}

//_____________________________________________________________________________
Bool_t
ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingNormal()
{
#if HodoCut
  static const auto MinDeBH2 = gUser.GetParameter("DeBH2", 0);
  static const auto MaxDeBH2 = gUser.GetParameter("DeBH2", 1);
  static const auto MinDeBH1 = gUser.GetParameter("DeBH1", 0);
  static const auto MaxDeBH1 = gUser.GetParameter("DeBH1", 1);
  static const auto MinBeamToF = gUser.GetParameter("BTOF", 1);
  static const auto MaxBeamToF = gUser.GetParameter("BTOF", 1);
#endif
  static const auto MaxMultiHitBcOut = gUser.GetParameter("MaxMultiHitBcOut");
  static const auto MaxMultiHitSdcIn = gUser.GetParameter("MaxMultiHitSdcIn");
#if TotCut
  static const auto MinTotBcOut = gUser.GetParameter("MinTotBcOut");
  static const auto MinTotSDC1  = gUser.GetParameter("MinTotSDC1");
  static const auto MinTotSDC2  = gUser.GetParameter("MinTotSDC2");
#endif

  RawData rawData;
  rawData.DecodeHits("TFlag");
  rawData.DecodeHits("BH1");
  rawData.DecodeHits("BH2");
  for(const auto& name: DCNameList.at("BcOut")) rawData.DecodeHits(name);
  for(const auto& name: DCNameList.at("SdcIn")) rawData.DecodeHits(name);
  HodoAnalyzer hodoAna(rawData);
  DCAnalyzer   DCAna(rawData);

  event.evnum = gUnpacker.get_event_number();

  // Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  for(const auto& hit: rawData.GetHodoRawHitContainer("TFlag")){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc();
    if(tdc > 0){
      event.trigpat[trigger_flag.count()] = seg;
      event.trigflag[seg] = tdc;
      trigger_flag.set(seg);
    }
  }

  HF1(1, 0.);

  if(trigger_flag[trigger::kSpillOnEnd] || trigger_flag[trigger::kSpillOffEnd])
    return true;

  HF1(1, 1.);

  //////////////BH2 time 0
  hodoAna.DecodeHits("BH2");
  Int_t nhBh2 = hodoAna.GetNHits("BH2");
  event.nhBh2 = nhBh2;
#if HodoCut
  if(nhBh2==0) return true;
#endif
  HF1(1, 2);

  //////////////BH2 Analysis
  for(Int_t i=0; i<nhBh2; ++i){
    const auto& hit = hodoAna.GetHit("BH2", i);
    if(!hit) continue;
    Double_t seg = hit->SegmentId()+1;
    Double_t cmt = hit->CMeanTime();
    Double_t dE  = hit->DeltaE();

#if HodoCut
    if(dE<MinDeBH2 || MaxDeBH2<dE) continue;
#endif
    event.tBh2[i]   = cmt;
    event.deBh2[i]  = dE;
    event.Bh2Seg[i] = seg;
  }

  const auto& cl_time0 = hodoAna.GetTime0BH2Cluster();
  if(cl_time0){
    event.Time0Seg = cl_time0->MeanSeg()+1;
    event.deTime0  = cl_time0->DeltaE();
    event.Time0    = cl_time0->Time0();
    event.CTime0   = cl_time0->CTime0();
  } else {
#if HodoCut
    return true;
#endif
  }

  HF1(1, 3.);

  //////////////BH1 Analysis
  hodoAna.DecodeHits("BH1");
  Int_t nhBh1 = hodoAna.GetNHits("BH1");
  event.nhBh1 = nhBh1;
#if HodoCut
  if(nhBh1==0) return true;
#endif
  HF1(1, 4);

  for(Int_t i=0; i<nhBh1; ++i){
    const auto& hit = hodoAna.GetHit("BH1", i);
    if(!hit) continue;
    Double_t cmt = hit->CMeanTime();
    Double_t dE  = hit->DeltaE();
#if HodoCut
    if(dE<MinDeBH1 || MaxDeBH1<dE) continue;
    if(btof<MinBeamToF || MaxBeamToF<btof) continue;
#endif
    event.tBh1[i]  = cmt;
    event.deBh1[i] = dE;
  }

  // Double_t btof0 = hodoAna.Btof0();
  // event.btof = btof0;

  HF1(1, 5.);

  //////////////BC3&4 number of hit layer
#if TdcCut
  rawData.TdcCutBCOut();
#endif
  DCAna.DecodeBcOutHits();
#if TotCut
  DCAna.TotCutBCOut(MinTotBcOut);
#endif
  // DCAna.DriftTimeCutBC34(-10, 50);
  // DCAna.MakeBH2DCHit(event.Time0Seg-1);
  Double_t multi_BcOut = 0.;
  for(Int_t plane=0; plane<NumOfLayersBcOut; ++plane){
    const auto& cont =DCAna.GetBcOutHC(plane);
    Int_t nh = cont.size();
    multi_BcOut += Double_t(nh);
  }
  if(multi_BcOut/Double_t(NumOfLayersBcOut) > MaxMultiHitBcOut){
    return true;
  }

  HF1(1, 11.);

  ///// SdcIn number of hit layer
#if TdcCut
  rawData.TdcCutSDCIn();
#endif
  DCAna.DecodeSdcInHits();
#if TotCut
  DCAna.TotCutSDC1(MinTotSDC1);
  DCAna.TotCutSDC2(MinTotSDC2);
#endif
  Double_t multi_SdcIn = 0.;
  for(Int_t plane=0; plane<NumOfLayersSdcIn; ++plane){
    const auto& cont =DCAna.GetSdcInHC(plane);
    Int_t nh = cont.size();
    multi_SdcIn += Double_t(nh);
  }
  if(multi_SdcIn/Double_t(NumOfLayersSdcIn) > MaxMultiHitSdcIn){
    return true;
  }

  HF1(1, 12.);

  //  std::cout << "==========TrackSearch BcOutSdcIn============" << std::endl;
  Bool_t status_tracking = DCAna.TrackSearchBcOutSdcIn();
  Int_t nt=DCAna.GetNtracksBcOutSdcIn();
  event.ntrack=nt;
  HF1(10, Double_t(nt));
  HF1(40, status_tracking? Double_t(nt) : -1);
  for(Int_t it=0; it<nt; ++it){
    const auto& tp=DCAna.GetTrackBcOutSdcIn(it);
    Int_t nh=tp->GetNHit();
    Double_t chisqr=tp->GetChiSquare();
    Double_t x0=tp->GetX0(), y0=tp->GetY0();
    Double_t u0=tp->GetU0(), v0=tp->GetV0();
    Double_t theta = tp->GetTheta();
    event.nlayer[it] = nh;
    event.chisqr[it]=chisqr;
    event.x0[it]=x0;
    event.y0[it]=y0;
    event.u0[it]=u0;
    event.v0[it]=v0;

    HF1(11, Double_t(nh));
    HF1(12, chisqr);
    HF1(14, x0); HF1(15, y0);
    HF1(16, u0); HF1(17, v0);
    HF2(18, x0, u0); HF2(19, y0, v0);
    HF2(20, x0, y0);

    Double_t xtgt=tp->GetX(zK18Target), ytgt=tp->GetY(zK18Target);
    Double_t utgt=u0, vtgt=v0;
    HF1(21, xtgt); HF1(22, ytgt);
    HF1(23, utgt); HF1(24, vtgt);
    HF2(25, xtgt, utgt); HF2(26, ytgt, vtgt);
    HF2(27, xtgt, ytgt);

    Double_t Xangle = -1000*std::atan(u0);

    HF2(51, -tp->GetX(245.), Xangle);
    HF2(52, -tp->GetX(600.), Xangle);
    HF2(53, -tp->GetX(1200.), Xangle);
    HF2(54, -tp->GetX(1600.), Xangle);

    HF2(61, tp->GetX(600.), tp->GetX(245.));
    HF2(62, tp->GetX(1200.), tp->GetX(245.));
    HF2(63, tp->GetX(1600.), tp->GetX(245.));

    for(Int_t ih=0; ih<nh; ++ih){
      DCLTrackHit *hit=tp->GetHit(ih);
      Int_t layerId = hit->GetLayer();
      if(layerId >= PlMinBcOut) layerId -= PlMinBcOut-1;
      else layerId += NumOfLayersBcOut;
      HF1(13, layerId);
      Double_t wire=hit->GetWire();
      Double_t dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1(100*layerId+11, wire+0.5);
      HF1(100*layerId+12, dt);
      HF1(100*layerId+13, dl);
      // HF1(10000*layerId+ 5000 +(Int_t)wire, dt);
      Double_t xcal=hit->GetXcal(), ycal=hit->GetYcal();
      Double_t pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      // if(layerId > 12){
      //   res += 49.95*TMath::Cos(hit->GetTiltAngle()*TMath::DegToRad());
      // }
      event.pos[layerId-1]=pos;
      event.res[layerId-1]=res;
      HF1(100*layerId+14, pos);
      HF1(100*layerId+15, res);
      HF2(100*layerId+16, pos, res);
      HF2(100*layerId+17, xcal, ycal);
      //      HF1(100000*layerId+50000+wire, res);
      Double_t wp=hit->GetWirePosition();
      Double_t sign=1.;
      if(pos-wp<0.) sign=-1;
      HF2(100*layerId+18, sign*dl, res);
      Double_t xlcal=hit->GetLocalCalPos();
      HF2(100*layerId+19, dt, xlcal-wp);

      Double_t tot = hit->GetTot();
      HF1(100*layerId+40, tot);

      if (theta>=0 && theta<15)
	HF1(100*layerId+71, res);
      else if (theta>=15 && theta<30)
	HF1(100*layerId+72, res);
      else if (theta>=30 && theta<45)
	HF1(100*layerId+73, res);
      else if (theta>=45)
	HF1(100*layerId+74, res);

      // if (std::abs(dl-std::abs(xlcal-wp))<2.0) {
      // 	HFProf(100*layerId+20, dt, std::abs(xlcal-wp));
      // 	HF2(100*layerId+22, dt, std::abs(xlcal-wp));
      // 	HFProf(100000*layerId+3000+(Int_t)wire, xlcal-wp,dt);
      // 	HF2(100000*layerId+4000+(Int_t)wire, xlcal-wp,dt);
      // }
    }
  }

  HF1(1, 15.);

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan:: InitializeHistograms()
{
  const Int_t nl = NumOfLayersBcOut+NumOfLayersSdcIn;

  const Int_t NbinBcOutDT   =  96;
  const Double_t MinBcOutDT =  -10.;
  const Double_t MaxBcOutDT =   70.;

  const Int_t NbinBcOutDL   = 80;
  const Double_t MinBcOutDL =  -0.5;
  const Double_t MaxBcOutDL =   3.5;

  const Int_t    NbinSDC1DT = 240;
  const Double_t MinSDC1DT  = -30.;
  const Double_t MaxSDC1DT  = 170.;
  const Int_t    NbinSDC1DL = 55;
  const Double_t MinSDC1DL  = -0.5;
  const Double_t MaxSDC1DL  = 5.0;

  const Int_t    NbinSDC2DT = 360;
  const Double_t MinSDC2DT  = -50.;
  const Double_t MaxSDC2DT  = 250.;
  const Int_t    NbinSDC2DL = 85;
  const Double_t MinSDC2DL  =  -0.5;
  const Double_t MaxSDC2DL  =   8.;

  HB1(1, "Status", 20, 0., 20.);

  //***********************Chamber
  // BC34 & SDC12
  for(Int_t i=1; i<=nl; ++i){
    TString tag;
    Int_t nwire;
    Int_t nbindt;
    Double_t mindt;
    Double_t maxdt;
    Int_t nbindl;
    Double_t mindl;
    Double_t maxdl;
    if(i <= NumOfLayersBc){
      tag = "BC3";
      nwire = MaxWireBC3;
      nbindt = NbinBcOutDT;
      mindt = MinBcOutDT;
      maxdt = MaxBcOutDT;
      nbindl = NbinBcOutDL;
      mindl = MinBcOutDL;
      maxdl = MaxBcOutDL;
    }else if(i <= NumOfLayersBcOut){
      tag = "BC4";
      nwire = MaxWireBC4;
      nbindt = NbinBcOutDT;
      mindt = MinBcOutDT;
      maxdt = MaxBcOutDT;
      nbindl = NbinBcOutDL;
      mindl = MinBcOutDL;
      maxdl = MaxBcOutDL;
    }else if(i <= NumOfLayersBcOut+NumOfLayersSDC1){
      tag = "SDC1";
      nwire = MaxWireSDC1;
      nbindt = NbinSDC1DT;
      mindt = MinSDC1DT;
      maxdt = MaxSDC1DT;
      nbindl = NbinSDC1DL;
      mindl = MinSDC1DL;
      maxdl = MaxSDC1DL;
    }else if(i <= NumOfLayersBcOut+NumOfLayersSdcIn){
      tag = "SDC2";
      nwire = MaxWireSDC2;
      nbindt = NbinSDC2DT;
      mindt = MinSDC2DT;
      maxdt = MaxSDC2DT;
      nbindl = NbinSDC2DL;
      mindl = MinSDC2DL;
      maxdl = MaxSDC2DL;
    }
    // Tracking Histgrams
    TString title11 = Form("HitPat BcOutSdcIn L%2d [Track]", i);
    TString title12 = Form("DriftTime BcOutSdcIn L%2d [Track]", i);
    TString title13 = Form("DriftLength BcOutSdcIn L%2d [Track]", i);
    TString title14 = Form("Position BcOutSdcIn L%2d", i);
    TString title15 = Form("Residual BcOutSdcIn L%2d", i);
    TString title16 = Form("Resid%%Pos BcOutSdcIn L%2d", i);
    TString title17 = Form("Y%%Xcal BcOutSdcIn L%2d", i);
    TString title18 = Form("Res%%DL BcOutSdcIn L%2d", i);
    TString title19 = Form("HitPos%%DriftTime BcOutSdcIn L%2d", i);
    TString title20 = Form("DriftLength%%DriftTime BcOutSdcIn L%2d", i);
    TString title21 = title15 + " [w/o Self]";
    TString title22 = title20;
    TString title40 = Form("TOT BcOutSdcIn L%2d [Track]", i);
    TString title71 = Form("Residual BcOutSdcIn L%2d (0<theta<15)", i);
    TString title72 = Form("Residual BcOutSdcIn L%2d (15<theta<30)", i);
    TString title73 = Form("Residual BcOutSdcIn L%2d (30<theta<45)", i);
    TString title74 = Form("Residual BcOutSdcIn L%2d (45<theta)", i);
    HB1(100*i+11, title11, nwire, 0., nwire);
    HB1(100*i+12, title12, nbindt, mindt, maxdt);
    HB1(100*i+13, title13, nbindl, mindl, maxdl);
    HB1(100*i+14, title14, 100, -250., 250.);
    HB1(100*i+15, title15, 400, -2.0, 2.0);
    HB2(100*i+16, title16, 250, -250., 250., 100, -1.0, 1.0);
    HB2(100*i+17, title17, 100, -250., 250., 100, -250., 250.);
    HB2(100*i+18, title18, 50, 3., 3., 100, -1.0, 1.0);
    HB2(100*i+19, title19, nbindt, mindt, maxdt,
        Int_t(maxdl*20), -maxdl, maxdl);
    HBProf(100*i+20, title20, 100, -5., 50., 0., 12.);
    HB2(100*i+22, title22, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
    HB1(100*i+21, title21, 200, -5.0, 5.0);
    HB1(100*i+40, title40, 300,    0, 300);
    HB1(100*i+71, title71, 200, -5.0, 5.0);
    HB1(100*i+72, title72, 200, -5.0, 5.0);
    HB1(100*i+73, title73, 200, -5.0, 5.0);
    HB1(100*i+74, title74, 200, -5.0, 5.0);
    for (Int_t wire=0; wire<nwire; wire++) {
      TString title = Form("XT of Layer %2d Wire #%4d", i, wire);
      HBProf(100000*i+3000+wire, title, 100, -4., 4., -5., 40.);
      HB2(100000*i+4000+wire, title, 100, -4., 4., 100, -5., 40.);
    }
  }

  // Tracking Histgrams
  HB1(10, "#Tracks BcOutSdcIn", 10, 0., 10.);
  HB1(11, "#Hits of Track BcOutSdcIn", 24, 0., 24.);
  HB1(12, "Chisqr BcOutSdcIn", 500, 0., 50.);
  HB1(13, "LayerId BcOutSdcIn", 24, 0., 24);
  HB1(14, "X0 BcOutSdcIn", 400, -100., 100.);
  HB1(15, "Y0 BcOutSdcIn", 400, -100., 100.);
  HB1(16, "U0 BcOutSdcIn", 200, -0.20, 0.20);
  HB1(17, "V0 BcOutSdcIn", 200, -0.20, 0.20);
  HB2(18, "U0%X0 BcOutSdcIn", 100, -100., 100., 100, -0.20, 0.20);
  HB2(19, "V0%Y0 BcOutSdcIn", 100, -100., 100., 100, -0.20, 0.20);
  HB2(20, "X0%Y0 BcOutSdcIn", 100, -100., 100., 100, -100, 100);

  HB1(21, "Xtgt BcOutSdcIn", 400, -100., 100.);
  HB1(22, "Ytgt BcOutSdcIn", 400, -100., 100.);
  HB1(23, "Utgt BcOutSdcIn", 200, -0.20, 0.20);
  HB1(24, "Vtgt BcOutSdcIn", 200, -0.20, 0.20);
  HB2(25, "Utgt%Xtgt BcOutSdcIn", 100, -100., 100., 100, -0.20, 0.20);
  HB2(26, "Vtgt%Ytgt BcOutSdcIn", 100, -100., 100., 100, -0.20, 0.20);
  HB2(27, "Xtgt%Ytgt BcOutSdcIn", 100, -100., 100., 100, -100, 100);

  HB1(38, "Plane Eff", 36, 0, 36);

  HB2(51, "X-X' 245 BcOutSdcIn", 400, -100., 100., 120, -60, 60);
  HB2(52, "X-X' 600 BcOutSdcIn", 400, -100., 100., 120, -60, 60);
  HB2(53, "X-X' 1200 BcOutSdcIn", 400, -100., 100., 120, -60, 60);
  HB2(54, "X-X' 1600 BcOutSdcIn", 400, -100., 100., 120, -60, 60);

  HB2(61, "X-X 600 BcOutSdcIn", 400, -100., 100., 400, -100, 100);
  HB2(62, "X-X 1200 BcOutSdcIn", 400, -100., 100., 400, -100, 100);
  HB2(63, "X-X 1600 BcOutSdcIn", 400, -100., 100., 400, -100, 100);

  // Analysis status
  HB1(40, "Tacking status", 11, -1., 10.);
  HB2(41, "BC3X0/BC4X1", 20, 0, 20, 20, 0, 20);

  ////////////////////////////////////////////
  //Tree
  HBTree("bcsdc", "tree of BcOutSdcInSdcInTracking");
  tree->Branch("evnum", &event.evnum, "evnum/I");
  tree->Branch("trigpat", event.trigpat, Form("trigpat[%d]/I", MaxHits));
  tree->Branch("trigflag", event.trigflag, Form("trigflag[%d]/I", NumOfSegTrig));

  //Hodoscope
  tree->Branch("nhBh1",    &event.nhBh1,   "nhBh1/I");
  tree->Branch("tBh1",      event.tBh1,    Form("tBh1[%d]/D",   MaxHits));
  tree->Branch("deBh1",     event.deBh1,   Form("deBh1[%d]/D",  MaxHits));

  tree->Branch("nhBh2",    &event.nhBh2,   "nhBh2/I");
  tree->Branch("tBh2",      event.tBh2,    Form("tBh2[%d]/D",   MaxHits));
  tree->Branch("deBh2",     event.deBh2,   Form("deBh2[%d]/D",  MaxHits));
  tree->Branch("Bh2Seg",    event.Bh2Seg,  Form("Bh2Seg[%d]/D", MaxHits));

  tree->Branch("Time0Seg", &event.Time0Seg,  "Time0Seg/D");
  tree->Branch("deTime0",  &event.deTime0,   "deTime0/D");
  tree->Branch("Time0",    &event.Time0,     "Time0/D");
  tree->Branch("CTime0",   &event.CTime0,    "CTime0/D");

  tree->Branch("ntrack",   &event.ntrack,   "ntrack/I");
  tree->Branch("nlayer",   &event.nlayer,   "nlayer[ntrack]/I");
  tree->Branch("chisqr",    event.chisqr,   "chisqr[ntrack]/D");
  tree->Branch("x0",        event.x0,       "x0[ntrack]/D");
  tree->Branch("y0",        event.y0,       "y0[ntrack]/D");
  tree->Branch("u0",        event.u0,       "u0[ntrack]/D");
  tree->Branch("v0",        event.v0,       "v0[ntrack]/D");
  tree->Branch("pos",      &event.pos,      Form("pos[%d]/D", nl));
  tree->Branch("res",      &event.res,      Form("res[%d]/D", nl));

  // HPrint();

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")        &&
     InitializeParameter<DCDriftParamMan>("DCDRFT") &&
     InitializeParameter<DCTdcCalibMan>("DCTDC")    &&
     InitializeParameter<HodoParamMan>("HDPRM")     &&
     InitializeParameter<HodoPHCMan>("HDPHC")       &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
