// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iostream>
#include <sstream>

// #include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCRawHit.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "HodoHit.hh"
#include "HodoCluster.hh"
#include "HodoAnalyzer.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "K18TrackD2U.hh"
#include "BH1Match.hh"
#include "BH2Filter.hh"
#include "S2sLib.hh"
#include "MathTools.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UnpackerManager.hh"

#define HodoCut 0 // with BH1/BH2
#define TIME_CUT 1 // in cluster analysis
#define TotCut  1 //for BcOut tracking
#define Chi2Cut  1 //for BcOut tracking

namespace
{
using namespace root;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
auto& gFilter = BH2Filter::GetInstance();
auto& gBH1Mth = BH1Match::GetInstance();
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;

  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  // BH1,2
  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;
  Double_t Btof0Seg;
  Double_t deBtof0;
  Double_t Btof0;
  Double_t CBtof0;

  // BFT
  Int_t    bft_ncl;
  Int_t    bft_ncl_bh1mth;
  Int_t    bft_clsize[NumOfSegBFT];
  Double_t bft_ctime[NumOfSegBFT];
  Double_t bft_clpos[NumOfSegBFT];
  Int_t    bft_bh1mth[NumOfSegBFT];

  // BcOut
  Int_t nlBcOut;
  Int_t ntBcOut;
  Int_t nhBcOut[MaxHits];
  Double_t chisqrBcOut[MaxHits];
  Double_t x0BcOut[MaxHits];
  Double_t y0BcOut[MaxHits];
  Double_t u0BcOut[MaxHits];
  Double_t v0BcOut[MaxHits];

  // K18
  Int_t ntK18;
  Int_t nhK18[MaxHits];
  Double_t p_2nd[MaxHits];
  Double_t p_3rd[MaxHits];
  Double_t delta_2nd[MaxHits];
  Double_t delta_3rd[MaxHits];

  Double_t xin[MaxHits];
  Double_t yin[MaxHits];
  Double_t uin[MaxHits];
  Double_t vin[MaxHits];

  Double_t xout[MaxHits];
  Double_t yout[MaxHits];
  Double_t uout[MaxHits];
  Double_t vout[MaxHits];

  Double_t chisqrK18[MaxHits];
  Double_t xtgtK18[MaxHits];
  Double_t ytgtK18[MaxHits];
  Double_t utgtK18[MaxHits];
  Double_t vtgtK18[MaxHits];
  Double_t xbftK18[MaxHits];
  Double_t ybftK18[MaxHits];
  Double_t ubftK18[MaxHits];
  Double_t vbftK18[MaxHits];

  Double_t theta[MaxHits];
  Double_t phi[MaxHits];

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum          =  0;
  bft_ncl        =  0;
  bft_ncl_bh1mth =  0;
  nlBcOut        =  0;
  ntBcOut        =  0;
  ntK18          =  0;
  Time0Seg       = -1;
  deTime0        = qnan;
  Time0          = qnan;
  CTime0         = qnan;
  Btof0Seg       = qnan;
  deBtof0        = qnan;
  Btof0          = qnan;
  CBtof0         = qnan;

  for(Int_t it=0; it<NumOfSegTrig; it++){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

  for(Int_t it=0; it<NumOfSegBFT; it++){
    bft_clsize[it] = qnan;
    bft_ctime[it]  = qnan;
    bft_clpos[it]  = qnan;
    bft_bh1mth[it] = -1;
  }

  for(Int_t i = 0; i<MaxHits; ++i){
    nhBcOut[i] = 0;
    chisqrBcOut[i] = qnan;
    x0BcOut[i] = qnan;
    y0BcOut[i] = qnan;
    u0BcOut[i] = qnan;
    v0BcOut[i] = qnan;
    p_2nd[i] = qnan;
    p_3rd[i] = qnan;
    delta_2nd[i] = qnan;
    delta_3rd[i] = qnan;
    xin[i] = qnan;
    yin[i] = qnan;
    uin[i] = qnan;
    vin[i] = qnan;
    xout[i] = qnan;
    yout[i] = qnan;
    uout[i] = qnan;
    vout[i] = qnan;
    nhK18[i] = 0;
    chisqrK18[i] = qnan;
    xtgtK18[i] = qnan;
    ytgtK18[i] = qnan;
    utgtK18[i] = qnan;
    vtgtK18[i] = qnan;
    xbftK18[i] = qnan;
    ybftK18[i] = qnan;
    ubftK18[i] = qnan;
    vbftK18[i] = qnan;
    theta[i] = qnan;
    phi[i] = qnan;
  }
}

//_____________________________________________________________________________
namespace root
{
Event event;
TH1   *h[MaxHist];
TTree *tree;
const Int_t BFTHid = 10000;
enum eParticle
{
  kKaon, kPion, nParticle
};
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
  static const Double_t MinDeBH2   = gUser.GetParameter("DeBH2", 0);
  static const Double_t MaxDeBH2   = gUser.GetParameter("DeBH2", 1);
  static const Double_t MinDeBH1   = gUser.GetParameter("DeBH1", 0);
  static const Double_t MaxDeBH1   = gUser.GetParameter("DeBH1", 1);
  static const Double_t MinBeamToF = gUser.GetParameter("BTOF",  1);
  static const Double_t MaxBeamToF = gUser.GetParameter("BTOF",  1);
#endif
#if TIME_CUT
  static const Double_t MinTimeBFT = gUser.GetParameter("TimeBFT", 0);
  static const Double_t MaxTimeBFT = gUser.GetParameter("TimeBFT", 1);
#endif
#if TotCut
  static const Double_t MinTotBcOut = gUser.GetParameter("MinTotBcOut", 0);
#endif
  // static const Double_t MaxMultiHitBcOut = gUser.GetParameter("MaxMultiHitBcOut");

  RawData rawData;
  rawData.DecodeHits();
  HodoAnalyzer hodoAna(rawData);
  DCAnalyzer   DCAna(rawData);

  event.evnum = gUnpacker.get_event_number();

  HF1(1, 0.);

  ///// Trigger Flag
  rawData.DecodeHits("TFlag");
  std::bitset<NumOfSegTrig> trigger_flag;
  {
    for(const auto& hit: rawData.GetHodoRawHitContainer("TFlag")){
      Int_t seg = hit->SegmentId();
      Int_t tdc = hit->GetTdc();
      if(tdc > 0){
	event.trigpat[trigger_flag.count()] = seg;
	event.trigflag[seg]  = tdc;
        trigger_flag.set(seg);
	HF1(10, seg);
	HF1(10+seg, tdc);
      }
    }
  }

  if(trigger_flag[trigger::kSpillOnEnd] || trigger_flag[trigger::kSpillOffEnd])
    return true;

  HF1(1, 1);

  ////////// BH2 time 0
  hodoAna.DecodeHits<BH2Hit>("BH2");
#if HodoCut
  Int_t nhBh2 = hodoAna.GetNHits("BH2");
  if(nhBh2==0) return true;
#endif
  HF1(1, 2);

  //////////////BH2 Analysis
  const auto& cl_time0 = hodoAna.GetTime0BH2Cluster();
  if(cl_time0){
    event.Time0Seg = cl_time0->MeanSeg();
    event.deTime0  = cl_time0->DeltaE();
    event.Time0    = cl_time0->Time0();
    event.CTime0   = cl_time0->CTime0();
  }else{
    return true;
  }
  HF1(1, 3);

  ////////// BH1 Analysis
  hodoAna.DecodeHits("BH1");
#if HodoCut
  Int_t nhBh1 = hodoAna.GetNHits("BH1");
  if(nhBh1==0) return true;
#endif
  HF1(1, 4);

  const auto& cl_btof0 = hodoAna.GetBtof0BH1Cluster();
  if(cl_btof0){
    event.Btof0Seg = cl_btof0->MeanSeg();
    event.deBtof0  = cl_btof0->DeltaE();
    event.Btof0    = cl_btof0->MeanTime() - event.Time0;
    event.CBtof0   = cl_btof0->CMeanTime() - event.CTime0;
  }else{
    return true;
  }

  HF1(1, 5);

  HF1(1, 6);

  std::vector<Double_t> xCand;
  ////////// BFT
  hodoAna.DecodeHits<FiberHit>("BFT");
  {
    // Fiber Cluster
    Int_t ncl_raw = hodoAna.GetNClusters("BFT");
#if TIME_CUT
    hodoAna.TimeCut("BFT", MinTimeBFT, MaxTimeBFT);
#endif
    Int_t ncl = hodoAna.GetNClusters("BFT");
    event.bft_ncl = ncl;
    HF1(BFTHid +100, ncl_raw);
    HF1(BFTHid +101, ncl);
    for(Int_t i=0; i<ncl; ++i){
      const auto& cl = hodoAna.GetCluster("BFT", i);
      if(!cl) continue;
      Double_t clsize = cl->ClusterSize();
      Double_t ctime  = cl->CMeanTime();
      Double_t pos    = cl->MeanPosition();
      // Double_t width  = cl->Width();

      event.bft_clsize[i] = clsize;
      event.bft_ctime[i]  = ctime;
      event.bft_clpos[i]  = pos;

      if(event.Btof0Seg >= 0 && ncl != 1){
	if(gBH1Mth.Judge(pos, event.Btof0Seg)){
	  event.bft_bh1mth[i] = 1;
	  xCand.push_back(pos);
	}
      }else{
	xCand.push_back(pos);
      }

      HF1(BFTHid +102, clsize);
      HF1(BFTHid +103, ctime);
      HF1(BFTHid +104, pos);
    }

    event.bft_ncl_bh1mth = xCand.size();
    HF1(BFTHid + 105, event.bft_ncl_bh1mth);
  }

  HF1(1, 7.);
  if(xCand.size()!=1) return true;
  HF1(1, 10.);

  DCAna.DecodeRawHits();
#if TotCut
  DCAna.TotCutBCOut(MinTotBcOut);
#endif
  ////////////// BC3&4 number of hit in one layer not 0
  Double_t multi_BcOut=0.;
  {
    Int_t nlBcOut = 0;
    for(Int_t l=0; l<NumOfLayersBcOut; ++l){
      const auto& contBcOut = DCAna.GetBcOutHC(l);
      Int_t nhBcOut = contBcOut.size();
      multi_BcOut += Double_t(nhBcOut);
      if(nhBcOut>0) nlBcOut++;
    }
    event.nlBcOut = nlBcOut;
  }

  // if(multi_BcOut/Double_t(NumOfLayersBcOut) > MaxMultiHitBcOut) return true;

  HF1(1, 11.);

  //////////////BCOut tracking
  // BH2Filter::FilterList cands;
  // gFilter.Apply((Int_t)event.Time0Seg-1, *DCAna, cands);
  //DCAna.TrackSearchBcOut(cands, event.Time0Seg-1);
  //  DCAna.TrackSearchBcOut(-1);
  //  DCAna.ChiSqrCutBcOut(10);

  DCAna.TrackSearchBcOut();
 #if Chi2Cut
  DCAna.ChiSqrCutBcOut(20);
 #endif

  Int_t ntBcOut = DCAna.GetNtracksBcOut();
  event.ntBcOut = ntBcOut;
  if(ntBcOut > MaxHits){
    std::cout << "#W too many BcOut tracks : ntBcOut = "
	      << ntBcOut << std::endl;
    ntBcOut = MaxHits;
  }
  HF1(50, Double_t(ntBcOut));
  for(Int_t it=0; it<ntBcOut; ++it){
    const auto& track = DCAna.GetTrackBcOut(it);
    Int_t nh = track->GetNHit();
    Double_t chisqr = track->GetChiSquare();
    Double_t u0 = track->GetU0(),  v0 = track->GetV0();
    Double_t x0 = track->GetX(0.), y0 = track->GetY(0.);

    HF1(51, Double_t(nh));
    HF1(52, chisqr);
    HF1(54, x0); HF1(35, y0);
    HF1(56, u0); HF1(37, v0);
    HF2(58, x0, u0); HF2(39, y0, v0);
    HF2(60, x0, y0);

    event.nhBcOut[it] = nh;
    event.chisqrBcOut[it] = chisqr;
    event.x0BcOut[it] = x0;
    event.y0BcOut[it] = y0;
    event.u0BcOut[it] = u0;
    event.v0BcOut[it] = v0;

    for(Int_t ih=0; ih<nh; ++ih){
      const auto& hit = track->GetHit(ih);
      Int_t layerId=hit->GetLayer()-100;
      HF1(53, layerId);
    }
  }
  if(ntBcOut==0) return true;

  HF1(1, 12.);

  HF1(1, 20.);
  ////////// K18Tracking D2U
  DCAna.TrackSearchK18D2U(xCand);
  Int_t ntK18 = DCAna.GetNTracksK18D2U();
  if(ntK18 > MaxHits){
    std::cout << "#W too many ntK18 "
	      << ntK18 << "/" << MaxHits << std::endl;
    ntK18 = MaxHits;
  }
  event.ntK18 = ntK18;
  HF1(70, Double_t(ntK18));
  for(Int_t i=0; i<ntK18; ++i){
    const auto& k18track = DCAna.GetK18TrackD2U(i);
    if(!k18track) continue;
    const auto& ltrack = k18track->TrackOut();
    std::size_t nh = ltrack->GetNHit();
    Double_t chisqr = ltrack->GetChiSquare();

    Double_t xin=k18track->Xin(), yin=k18track->Yin();
    Double_t uin=k18track->Uin(), vin=k18track->Vin();

    Double_t xt=k18track->Xtgt(), yt=k18track->Ytgt();
    Double_t ut=k18track->Utgt(), vt=k18track->Vtgt();

    Double_t xb=k18track->Xbft(), yb=k18track->Ybft();
    Double_t ub=k18track->Ubft(), vb=k18track->Vbft();

    Double_t xout=k18track->Xout(), yout=k18track->Yout();
    Double_t uout=k18track->Uout(), vout=k18track->Vout();

    Double_t p_2nd=k18track->P();
    Double_t p_3rd=k18track->P3rd();
    Double_t delta_2nd=k18track->Delta();
    Double_t delta_3rd=k18track->Delta3rd();
    Double_t theta = ltrack->GetTheta();
    Double_t phi   = ltrack->GetPhi();

    HF1(74, xt); HF1(75, yt); HF1(76, ut); HF1(77, vt);
    HF2(78, xt, ut); HF2(79, yt, vt); HF2(80, xt, yt);
    HF1(84, xb); HF1(85, yb); HF1(86, ub); HF1(87, vb);
    HF2(88, xb, ub); HF2(89, yb, vb); HF2(90, xb, yb);
    HF1(91, p_3rd); HF1(92, delta_3rd);

    event.p_2nd[i] = p_2nd;
    event.p_3rd[i] = p_3rd;
    event.delta_2nd[i] = delta_2nd;
    event.delta_3rd[i] = delta_3rd;

    event.xin[i] = xin;
    event.yin[i] = yin;
    event.uin[i] = uin;
    event.vin[i] = vin;

    event.xout[i] = xout;
    event.yout[i] = yout;
    event.uout[i] = uout;
    event.vout[i] = vout;

    event.nhK18[i]     = nh;
    event.chisqrK18[i] = chisqr;
    event.xtgtK18[i]   = xt;
    event.ytgtK18[i]   = yt;
    event.utgtK18[i]   = ut;
    event.vtgtK18[i]   = vt;
    event.xbftK18[i]   = xb;
    event.ybftK18[i]   = yb;
    event.ubftK18[i]   = ub;
    event.vbftK18[i]   = vb;
    event.theta[i] = theta;
    event.phi[i]   = phi;
  }

  HF1(1, 22.);

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
  const Int_t    NbinTot =  136;
  const Double_t MinTot  =   -8.;
  const Double_t MaxTot  =  128.;
  const Int_t    NbinTime = 1000;
  const Double_t MinTime  = -500.;
  const Double_t MaxTime  =  500.;

  HB1(1, "Status",  30,   0., 30.);
  HB1(10, "Trigger HitPat", NumOfSegTrig, 0, NumOfSegTrig);
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1(10+i+1, Form("Trigger Flag %d", i+1), 0x1000, 0, 0x1000);
  }

  //BFT
  HB1(BFTHid +21, "BFT CTime U",     NbinTime, MinTime, MaxTime);
  HB1(BFTHid +22, "BFT CTime D",     NbinTime, MinTime, MaxTime);
  HB2(BFTHid +23, "BFT CTime/Tot U", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);
  HB2(BFTHid +24, "BFT CTime/Tot D", NbinTot, MinTot, MaxTot, NbinTime, MinTime, MaxTime);

  HB1(BFTHid +100, "BFT NCluster [Raw]", 100, 0, 100);
  HB1(BFTHid +101, "BFT NCluster [TimeCut]", 10, 0, 10);
  HB1(BFTHid +102, "BFT Cluster Size", 5, 0, 5);
  HB1(BFTHid +103, "BFT CTime (Cluster)", NbinTime, MinTime, MaxTime);
  HB1(BFTHid +104, "BFT Cluster Position",
       NumOfSegBFT, -0.5*(Double_t)NumOfSegBFT, 0.5*(Double_t)NumOfSegBFT);
  HB1(BFTHid +105, "BFT NCluster [TimeCut && BH1Matching]", 10, 0, 10);

  // BcOut
  HB1(50, "#Tracks BcOut", 10, 0., 10.);
  HB1(51, "#Hits of Track BcOut", 20, 0., 20.);
  HB1(52, "Chisqr BcOut", 500, 0., 50.);
  HB1(53, "LayerId BcOut", 15, 12., 27.);
  HB1(54, "X0 BcOut", 400, -100., 100.);
  HB1(55, "Y0 BcOut", 400, -100., 100.);
  HB1(56, "U0 BcOut",  200, -0.20, 0.20);
  HB1(57, "V0 BcOut",  200, -0.20, 0.20);
  HB2(58, "U0%X0 BcOut", 100, -100., 100., 100, -0.20, 0.20);
  HB2(59, "V0%Y0 BcOut", 100, -100., 100., 100, -0.20, 0.20);
  HB2(60, "X0%Y0 BcOut", 100, -100., 100., 100, -100., 100.);

  // K18
  HB1(70, "#Tracks K18", 20, 0., 20.);
  HB1(71, "#Hits of K18Track", 30, 0., 30.);
  HB1(72, "Chisqr K18Track", 500, 0., 100.);
  HB1(73, "LayerId K18Track", 50, 0., 50.);
  HB1(74, "Xtgt K18Track", 200, -100., 100.);
  HB1(75, "Ytgt K18Track", 200, -100., 100.);
  HB1(76, "Utgt K18Track", 300, -0.30, 0.30);
  HB1(77, "Vtgt K18Track", 300, -0.20, 0.20);
  HB2(78, "U%Xtgt K18Track", 100, -100., 100., 100, -0.25, 0.25);
  HB2(79, "V%Ytgt K18Track", 100, -100., 100., 100, -0.10, 0.10);
  HB2(80, "Y%Xtgt K18Track", 100, -100., 100., 100, -100., 100.);
  HB1(84, "Xbft K18Track", 200, -100., 100.);
  HB1(85, "Ybft K18Track", 200, -100., 100.);
  HB1(86, "Ubft K18Track", 300, -0.30, 0.30);
  HB1(87, "Vbft K18Track", 300, -0.20, 0.20);
  HB2(88, "U%Xbft K18Track", 100, -100., 100., 100, -0.25, 0.25);
  HB2(89, "V%Ybft K18Track", 100, -100., 100., 100, -0.10, 0.10);
  HB2(90, "Y%Xbft K18Track", 100, -100., 100., 100, -100., 100.);
  HB1(91, "P K18Track", 500, 0.50, 2.0);
  HB1(92, "dP K18Track", 200, -0.1, 0.1);

  //tree
  HBTree("k18track","Data Summary Table of K18Tracking");
  // Trigger Flag
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  tree->Branch("Time0Seg", &event.Time0Seg,  "Time0Seg/D");
  tree->Branch("deTime0",  &event.deTime0,   "deTime0/D");
  tree->Branch("Time0",    &event.Time0,     "Time0/D");
  tree->Branch("CTime0",   &event.CTime0,    "CTime0/D");
  tree->Branch("Btof0Seg", &event.Btof0Seg,  "Btof0Seg/D");
  tree->Branch("deBtof0",  &event.deBtof0,   "deBtof0/D");
  tree->Branch("Btof0",    &event.Btof0,     "Btof0/D");
  tree->Branch("CBtof0",   &event.CBtof0,    "CBtof0/D");

  //BFT
  tree->Branch("bft_ncl",        &event.bft_ncl,    "bft_ncl/I");
  tree->Branch("bft_ncl_bh1mth", &event.bft_ncl_bh1mth, "bft_ncl_bh1mth/I");
  tree->Branch("bft_clsize",      event.bft_clsize, "bft_clsize[bft_ncl]/I");
  tree->Branch("bft_ctime",       event.bft_ctime,  "bft_ctime[bft_ncl]/D");
  tree->Branch("bft_clpos",       event.bft_clpos,  "bft_clpos[bft_ncl]/D");
  tree->Branch("bft_bh1mth",      event.bft_bh1mth, "bft_bh1mth[bft_ncl]/I");

  // BcOut
  tree->Branch("nlBcOut",   &event.nlBcOut,     "nlBcOut/I");
  tree->Branch("ntBcOut",   &event.ntBcOut,     "ntBcOut/I");
  tree->Branch("nhBcOut",    event.nhBcOut,     "nhBcOut[ntBcOut]/I");
  tree->Branch("chisqrBcOut",event.chisqrBcOut, "chisqrIn[ntBcOut]/D");
  tree->Branch("x0BcOut",    event.x0BcOut,     "x0BcOut[ntBcOut]/D");
  tree->Branch("y0BcOut",    event.y0BcOut,     "y0BcOut[ntBcOut]/D");
  tree->Branch("u0BcOut",    event.u0BcOut,     "u0BcOut[ntBcOut]/D");
  tree->Branch("v0BcOut",    event.v0BcOut,     "v0BcOut[ntBcOut]/D");

  // K18
  tree->Branch("ntK18",      &event.ntK18,     "ntK18/I");
  tree->Branch("nhK18",       event.nhK18,     "nhK18[ntK18]/I");
  tree->Branch("chisqrK18",   event.chisqrK18, "chisqrK18[ntK18]/D");
  tree->Branch("p_2nd",       event.p_2nd,     "p_2nd[ntK18]/D");
  tree->Branch("p_3rd",       event.p_3rd,     "p_3rd[ntK18]/D");
  tree->Branch("delta_2nd",   event.delta_2nd, "delta_2nd[ntK18]/D");
  tree->Branch("delta_3rd",   event.delta_3rd, "delta_3rd[ntK18]/D");
  tree->Branch("xin",    event.xin,   "xin[ntK18]/D");
  tree->Branch("yin",    event.yin,   "yin[ntK18]/D");
  tree->Branch("uin",    event.uin,   "uin[ntK18]/D");
  tree->Branch("vin",    event.vin,   "vin[ntK18]/D");
  tree->Branch("xout",    event.xout,   "xout[ntK18]/D");
  tree->Branch("yout",    event.yout,   "yout[ntK18]/D");
  tree->Branch("uout",    event.uout,   "uout[ntK18]/D");
  tree->Branch("vout",    event.vout,   "vout[ntK18]/D");
  tree->Branch("xtgtK18",    event.xtgtK18,   "xtgtK18[ntK18]/D");
  tree->Branch("ytgtK18",    event.ytgtK18,   "ytgtK18[ntK18]/D");
  tree->Branch("utgtK18",    event.utgtK18,   "utgtK18[ntK18]/D");
  tree->Branch("vtgtK18",    event.vtgtK18,   "vtgtK18[ntK18]/D");
  tree->Branch("xbftK18",    event.xbftK18,   "xbftK18[ntK18]/D");
  tree->Branch("ybftK18",    event.ybftK18,   "ybftK18[ntK18]/D");
  tree->Branch("ubftK18",    event.ubftK18,   "ubftK18[ntK18]/D");
  tree->Branch("vbftK18",    event.vbftK18,   "vbftK18[ntK18]/D");
  tree->Branch("theta",   event.theta,  "theta[ntK18]/D");
  tree->Branch("phi",     event.phi,    "phi[ntK18]/D");

  HPrint();
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
     InitializeParameter<BH2Filter>("BH2FLT")       &&
     InitializeParameter<BH1Match>("BH1MTH")        &&
     InitializeParameter<K18TransMatrix>("K18TM")   &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
