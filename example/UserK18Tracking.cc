// -*- C++ -*-

#include "VEvent.hh"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <TString.h>
#include <TMinuit.h>
#include <TMath.h>
#include <TVector3.h>

#include "ConfMan.hh"
#include "DCAnalyzer.hh"
#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DCLTrackHit.hh"
#include "DCLocalTrack.hh"
#include "DCRawHit.hh"
#include "DCTdcCalibMan.hh"
#include "DetectorID.hh"
#include "EventAnalyzer.hh"
#include "FieldMan.hh"
#include "HSTrack.hh"
#include "HistTools.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "TransferMatrixMan.hh"
#include "UserParamMan.hh"

#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <UnpackerXMLReadDigit.hh>

namespace
{
const auto& gUser = UserParamMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
// non-const ref: CalcBLC1toBLC2 is a non-const method
auto& gTM = TransferMatrixMan::GetInstance();

using seg_t = std::vector<Double_t>;
using tdc_t = std::vector<std::vector<Double_t>>;
TTree* tree;

// number of TMinuit fit parameters (x, u, y, v, dp)
const Int_t    kNParam   = 5;
// transfer matrix order (1 = linear, 2 = with D5Matrix2nd corrections)
const Int_t    kMatrixOrder = 2;
// matrix reference planes (fixed): entrance VI = BLC1a/b pair center
// (BcIn frame origin), exit VO = BLC2a/b pair center (1300 mm upstream
// of FF in the BcOut frame) [mm]
const Double_t kZVI = 0.;
const Double_t kZVO = -1300.;
// max local tracks per side used in the combinatorial fit
const Int_t    kMaxNTrack = 10;

// per-hit data handed to the TMinuit fcn
struct FitHit
{
  Double_t z;     // hit z in the local tracking frame [mm]
  Double_t tilt;  // wire tilt angle [rad]
  Double_t s;     // measured local hit position [mm]
  Double_t reso;  // resolution [mm]
};

// globals for the TMinuit fcn
std::vector<FitHit> g_blc1_hits;  // BLC1 hits, track defined at VI plane
std::vector<FitHit> g_blc2_hits;  // BLC2 hits, track defined at VO plane

// K18 tracks are all exclusive BcIn-BcOut pairs, selected greedily
// from the best chisqrK18 combination; a local track can join only
// one K18 track. All per-track quantities below are indexed by
// the K18 track (best chisqrK18 first).
struct Event
{
UInt_t run_number;
UInt_t event_number;
beam::EBeamFlag beam_flag;
tdc_t trig_flag;
seg_t trig_pat;

Int_t    ntBcIn;
Int_t    ntBcOut;
Int_t    ntK18;
// local-tracking quality of the paired BcIn/BcOut DCLocalTracks
seg_t chisqrBcIn;   // local reduced chisqr
seg_t chisqrBcOut;
tdc_t resBcIn;      // local-track residuals [mm]
tdc_t resBcOut;
// combined K18 momentum fit
seg_t chisqrK18;    // reduced chisqr (chisqr/ndf)
std::vector<Int_t> ndfK18;
tdc_t resK18BcIn;   // hit residuals at the fit params [mm]
tdc_t resK18BcOut;
seg_t delta;        // dp = (p - p0)/p0 fraction
seg_t pk18;         // |p| [GeV/c]
seg_t px;
seg_t py;
seg_t pz;
// fitted state at the entrance plane (VI)
seg_t xin;
seg_t yin;
seg_t uin;
seg_t vin;
// fitted state at VO (matrix output with fitted dp)
seg_t xout;
seg_t yout;
seg_t uout;
seg_t vout;
// measured BcOut local track at the VO plane
seg_t x0BcOut;
seg_t y0BcOut;
seg_t u0BcOut;
seg_t v0BcOut;
  
std::vector<Int_t> rkStatusHS;
tdc_t xvpHS;
tdc_t yvpHS;
tdc_t zvpHS;
tdc_t uvpHS;
tdc_t vvpHS;
tdc_t pvpHS;
tdc_t pxvpHS;
tdc_t pyvpHS;
tdc_t pzvpHS;
};
Event event;

// one BcIn-BcOut combination fitted through the D5 matrix
struct K18Fit
{
  Int_t    ibcin;
  Int_t    ibcout;
  Double_t chisqr;  // reduced
  Int_t    ndf;
  Double_t par[kNParam];  // {xin, uin, yin, vin, delta}
  Double_t vo[4];         // {xout, uout, yout, vout}
  std::vector<Double_t> res1;
  std::vector<Double_t> res2;
};

//_____________________________________________________________________________
// chi2 = sum of BLC1 residuals (at VI) + BLC2 residuals (propagated to VO)
void
K18FitFCN(Int_t& /*npar*/, Double_t* /*gin*/, Double_t& f,
          Double_t* par, Int_t /*flag*/)
{
  const Double_t xin = par[0];
  const Double_t uin = par[1];
  const Double_t yin = par[2];
  const Double_t vin = par[3];
  const Double_t dp  = par[4];

  Double_t chisqr = 0.;

  // BLC1: straight line defined at the VI plane (z = kZVI)
  for(const auto& h : g_blc1_hits){
    Double_t dz   = h.z - kZVI;
    Double_t xcal = xin + uin*dz;
    Double_t ycal = yin + vin*dz;
    Double_t scal = xcal*TMath::Cos(h.tilt) + ycal*TMath::Sin(h.tilt);
    Double_t res  = (h.s - scal)/h.reso;
    chisqr += res*res;
  }

  // propagate VI -> VO through the D5 transfer matrix
  Double_t parblc1[kNParam] = {xin, uin, yin, vin, dp};
  Double_t parblc2[kNParam] = {0., 0., 0., 0., 0.};
  gTM.CalcBLC1toBLC2(parblc1, parblc2, kMatrixOrder);
  const Double_t x_vo = parblc2[0], u_vo = parblc2[1];
  const Double_t y_vo = parblc2[2], v_vo = parblc2[3];

  // BLC2: matrix output is the state at the VO plane (z = kZVO)
  for(const auto& h : g_blc2_hits){
    Double_t dz2  = h.z - kZVO;
    Double_t xcal = x_vo + u_vo*dz2;
    Double_t ycal = y_vo + v_vo*dz2;
    Double_t scal = xcal*TMath::Cos(h.tilt) + ycal*TMath::Sin(h.tilt);
    Double_t res  = (h.s - scal)/h.reso;
    chisqr += res*res;
  }

  f = chisqr;
}

//_____________________________________________________________________________
void
FillFitHits(const DCLocalTrack* track, std::vector<FitHit>& out)
{
  out.clear();
  for(Int_t i=0, n=track->GetNHit(); i<n; ++i){
    const DCLTrackHit* h = track->GetHit(i);
    FitHit fh;
    fh.z    = h->GetZ();
    fh.tilt = h->GetTiltAngle()*TMath::DegToRad();
    fh.s    = h->GetLocalHitPos();
    fh.reso = h->GetResolution();
    out.push_back(fh);
  }
}

//_____________________________________________________________________________
// fit one BcIn-BcOut combination; true if MIGRAD converged
Bool_t
DoK18Fit(const DCLocalTrack* blc1, const DCLocalTrack* blc2, K18Fit& fit)
{
  FillFitHits(blc1, g_blc1_hits);
  FillFitHits(blc2, g_blc2_hits);

  fit.ndf = (Int_t)(g_blc1_hits.size()+g_blc2_hits.size()) - kNParam;
  if(fit.ndf <= 0) return false;

  // initial parameters: BLC1 local track drifted to the VI plane
  const Double_t xin0 = blc1->GetX0() + blc1->GetU0()*kZVI;
  const Double_t yin0 = blc1->GetY0() + blc1->GetV0()*kZVI;

  TMinuit minuit(kNParam);
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(K18FitFCN);
  minuit.DefineParameter(0, "x",  xin0,           0.1,   0., 0.);
  minuit.DefineParameter(1, "u",  blc1->GetU0(),  0.001, 0., 0.);
  minuit.DefineParameter(2, "y",  yin0,           0.1,   0., 0.);
  minuit.DefineParameter(3, "v",  blc1->GetV0(),  0.001, 0., 0.);
  minuit.DefineParameter(4, "dp", 0.,             0.001, 0., 0.);

  Double_t arglist[2] = {1000., 1.};
  Int_t ierflg = 0;
  minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
  if(ierflg != 0) return false;

  Double_t err;
  for(Int_t i=0; i<kNParam; ++i)
    minuit.GetParameter(i, fit.par[i], err);

  Double_t fmin, fedm, errdef;
  Int_t npari, nparx, istat;
  minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);
  fit.chisqr = fmin/fit.ndf;

  // residuals at the fitted parameters
  // (BLC1 at VI, BLC2 propagated to VO)
  const Double_t x = fit.par[0], u = fit.par[1];
  const Double_t y = fit.par[2], v = fit.par[3];
  fit.res1.clear();
  fit.res2.clear();
  for(const auto& h : g_blc1_hits){
    Double_t dz   = h.z - kZVI;
    Double_t scal = (x+u*dz)*TMath::Cos(h.tilt)
                  + (y+v*dz)*TMath::Sin(h.tilt);
    fit.res1.push_back(h.s - scal);
  }
  Double_t p2[kNParam] = {0., 0., 0., 0., 0.};
  gTM.CalcBLC1toBLC2(fit.par, p2, kMatrixOrder);
  for(Int_t i=0; i<4; ++i) fit.vo[i] = p2[i];
  for(const auto& h : g_blc2_hits){
    Double_t dz2  = h.z - kZVO;
    Double_t scal = (fit.vo[0]+fit.vo[1]*dz2)*TMath::Cos(h.tilt)
                  + (fit.vo[2]+fit.vo[3]*dz2)*TMath::Sin(h.tilt);
    fit.res2.push_back(h.s - scal);
  }
  return true;
}


}

//_____________________________________________________________________________
Bool_t
ProcessBegin()
{
  event.run_number = gUnpacker.get_run_number();
  event.event_number = gUnpacker.get_event_number();
  event.beam_flag = beam::kUnknown;
  event.trig_flag.clear();
  event.trig_pat.clear();

  event.ntBcIn  = 0;
  event.ntBcOut = 0;
  event.ntK18   = 0;
  event.chisqrBcIn.clear();
  event.chisqrBcOut.clear();
  event.resBcIn.clear();
  event.resBcOut.clear();
  event.chisqrK18.clear();
  event.ndfK18.clear();
  event.resK18BcIn.clear();
  event.resK18BcOut.clear();
  event.delta.clear();
  event.pk18.clear();
  event.px.clear();
  event.py.clear();
  event.pz.clear();
  event.xin.clear();
  event.yin.clear();
  event.uin.clear();
  event.vin.clear();
  event.xout.clear();
  event.yout.clear();
  event.uout.clear();
  event.vout.clear();
  event.x0BcOut.clear();
  event.y0BcOut.clear();
  event.u0BcOut.clear();
  event.v0BcOut.clear();
  event.rkStatusHS.clear();
  event.xvpHS.clear();
  event.yvpHS.clear();
  event.zvpHS.clear();
  event.uvpHS.clear();
  event.vvpHS.clear();
  event.pvpHS.clear();
  event.pxvpHS.clear();
  event.pyvpHS.clear();
  event.pzvpHS.clear();

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessNormal()
{
  using root::HF1;

  RawData rawData;
  for (const auto& name : DCNameList.at("BcIn"))  rawData.DecodeHits(name);
  for (const auto& name : DCNameList.at("BcOut")) rawData.DecodeHits(name);

  EventAnalyzer evAna;

  HF1("Status", 0);
  rawData.DecodeHits("TriggerFlag");
  evAna.TriggerFlag(rawData);

  HF1("Status", 1);
  rawData.DecodeHits("BAC"); // for beam_flag
  rawData.DecodeHits("BHT"); // for beam_flag
  event.beam_flag = evAna.BeamFlag(rawData);

  HF1("Status", 2);
  DCAnalyzer dcAna(rawData);

  // BcIn (BLC1)
  dcAna.DecodeBcInHits();
  dcAna.TotCut("BLC1a");
  dcAna.TotCut("BLC1b");
  dcAna.DriftTimeCut("BLC1a");
  dcAna.DriftTimeCut("BLC1b");
  dcAna.TrackSearchBcIn();

  // BcOut (BLC2)
  dcAna.DecodeBcOutHits();
  dcAna.TotCut("BLC2a");
  dcAna.TotCut("BLC2b");
  dcAna.DriftTimeCut("BLC2a");
  dcAna.DriftTimeCut("BLC2b");
  dcAna.TrackSearchBcOut();

  const auto& contIn  = dcAna.GetBcInTrackContainer();
  const auto& contOut = dcAna.GetBcOutTrackContainer();
  event.ntBcIn  = contIn.size();
  event.ntBcOut = contOut.size();
  if(contIn.empty() || contOut.empty()) return true;
  if(!gTM.IsReady())                    return true;

  HF1("Status", 3);

  // combinatorial fit of all BcIn-BcOut pairs
  const Int_t nIn  = TMath::Min((Int_t)contIn.size(),  kMaxNTrack);
  const Int_t nOut = TMath::Min((Int_t)contOut.size(), kMaxNTrack);
  std::vector<K18Fit> cands;
  for(Int_t i=0; i<nIn; ++i){
    for(Int_t j=0; j<nOut; ++j){
      K18Fit fit;
      fit.ibcin  = i;
      fit.ibcout = j;
      if(DoK18Fit(contIn[i], contOut[j], fit))
        cands.push_back(std::move(fit));
    }
  }
  std::sort(cands.begin(), cands.end(),
            [](const K18Fit& a, const K18Fit& b){
              return a.chisqr < b.chisqr;
            });

  // greedy exclusive pairing: a local track joins only one K18 track
  std::vector<Bool_t> used1(nIn, false), used2(nOut, false);
  const Double_t p0 = gTM.GetCentralMomentum();
  for(const auto& fit : cands){
    if(used1[fit.ibcin] || used2[fit.ibcout]) continue;
    used1[fit.ibcin]  = true;
    used2[fit.ibcout] = true;

    const DCLocalTrack* blc1 = contIn[fit.ibcin];
    const DCLocalTrack* blc2 = contOut[fit.ibcout];

    event.chisqrBcIn.push_back(blc1->GetChiSquare());
    event.chisqrBcOut.push_back(blc2->GetChiSquare());
    std::vector<Double_t> r1, r2;
    for(Int_t i=0, n=blc1->GetNHit(); i<n; ++i)
      r1.push_back(blc1->GetHit(i)->GetResidual());
    for(Int_t i=0, n=blc2->GetNHit(); i<n; ++i)
      r2.push_back(blc2->GetHit(i)->GetResidual());
    event.resBcIn.push_back(r1);
    event.resBcOut.push_back(r2);

    event.chisqrK18.push_back(fit.chisqr);
    event.ndfK18.push_back(fit.ndf);
    event.resK18BcIn.push_back(fit.res1);
    event.resK18BcOut.push_back(fit.res2);
    event.xin.push_back(fit.par[0]);
    event.uin.push_back(fit.par[1]);
    event.yin.push_back(fit.par[2]);
    event.vin.push_back(fit.par[3]);
    event.delta.push_back(fit.par[4]);
    event.xout.push_back(fit.vo[0]);
    event.uout.push_back(fit.vo[1]);
    event.yout.push_back(fit.vo[2]);
    event.vout.push_back(fit.vo[3]);

    // measured BcOut local track at the VO plane, to be compared
    // with the fitted VO state for frame-consistency checks
    event.x0BcOut.push_back(blc2->GetX0() + blc2->GetU0()*kZVO);
    event.u0BcOut.push_back(blc2->GetU0());
    event.y0BcOut.push_back(blc2->GetY0() + blc2->GetV0()*kZVO);
    event.v0BcOut.push_back(blc2->GetV0());

    // momentum from the fitted delta, directed along the fitted
    // VO-plane (exit) direction
    const Double_t p  = p0*(1. + fit.par[4]);
    const Double_t u  = fit.vo[1];
    const Double_t v  = fit.vo[3];
    const Double_t pz = p/TMath::Sqrt(1. + u*u + v*v);
    event.pk18.push_back(p);
    event.px.push_back(pz*u);
    event.py.push_back(pz*v);
    event.pz.push_back(pz);

    HSTrack hs(fit.vo[0], fit.vo[2], fit.vo[1], fit.vo[3], p);
    hs.Propagate();
    event.rkStatusHS.push_back(hs.IsPassed() ? 1 : 0);
    event.xvpHS.push_back(hs.VPX());
    event.yvpHS.push_back(hs.VPY());
    event.zvpHS.push_back(hs.VPZ());
    event.uvpHS.push_back(hs.VPU());
    event.vvpHS.push_back(hs.VPV());
    event.pvpHS.push_back(hs.VPP());
    std::vector<Double_t> pxvp, pyvp, pzvp;
    const auto& vp_momentum = hs.VPMomentum();
    pxvp.reserve(vp_momentum.size());
    pyvp.reserve(vp_momentum.size());
    pzvp.reserve(vp_momentum.size());
    for(const auto& mom : vp_momentum){
      pxvp.push_back(mom.X());
      pyvp.push_back(mom.Y());
      pzvp.push_back(mom.Z());
    }
    event.pxvpHS.push_back(std::move(pxvp));
    event.pyvpHS.push_back(std::move(pyvp));
    event.pzvpHS.push_back(std::move(pzvp));
  }
  event.ntK18 = event.chisqrK18.size();

  if(event.ntK18 > 0) HF1("Status", 4);

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
  hist::BuildStatus();
  hist::BuildTriggerFlag();

  tree = new TTree("k18", "Data Summary Table of K18Tracking");
  tree->Branch("run_number", &event.run_number);
  tree->Branch("event_number", &event.event_number);
  tree->Branch("beam_flag", &event.beam_flag, "beam_flag/I");
  tree->Branch("trig_flag", &event.trig_flag);
  tree->Branch("trig_pat", &event.trig_pat);

  // BcIn/BcOut local tracking (per K18 track)
  tree->Branch("ntBcIn", &event.ntBcIn);
  tree->Branch("ntBcOut", &event.ntBcOut);
  tree->Branch("chisqrBcIn", &event.chisqrBcIn);
  tree->Branch("chisqrBcOut", &event.chisqrBcOut);
  tree->Branch("resBcIn", &event.resBcIn);
  tree->Branch("resBcOut", &event.resBcOut);
  tree->Branch("x0BcOut", &event.x0BcOut);
  tree->Branch("y0BcOut", &event.y0BcOut);
  tree->Branch("u0BcOut", &event.u0BcOut);
  tree->Branch("v0BcOut", &event.v0BcOut);

  // K18 momentum fit (sorted by chisqrK18, exclusive pairs)
  tree->Branch("ntK18", &event.ntK18);
  tree->Branch("chisqrK18", &event.chisqrK18);
  tree->Branch("ndfK18", &event.ndfK18);
  tree->Branch("resK18BcIn", &event.resK18BcIn);
  tree->Branch("resK18BcOut", &event.resK18BcOut);
  tree->Branch("delta", &event.delta);
  tree->Branch("pk18", &event.pk18);
  tree->Branch("px", &event.px);
  tree->Branch("py", &event.py);
  tree->Branch("pz", &event.pz);
  tree->Branch("xin", &event.xin);
  tree->Branch("yin", &event.yin);
  tree->Branch("uin", &event.uin);
  tree->Branch("vin", &event.vin);
  tree->Branch("xout", &event.xout);
  tree->Branch("yout", &event.yout);
  tree->Branch("uout", &event.uout);
  tree->Branch("vout", &event.vout);
  tree->Branch("rk_statusHS", &event.rkStatusHS);
  tree->Branch("xvpHS", &event.xvpHS);
  tree->Branch("yvpHS", &event.yvpHS);
  tree->Branch("zvpHS", &event.zvpHS);
  tree->Branch("uvpHS", &event.uvpHS);
  tree->Branch("vvpHS", &event.vvpHS);
  tree->Branch("pvpHS", &event.pvpHS);
  tree->Branch("pxvpHS", &event.pxvpHS);
  tree->Branch("pyvpHS", &event.pyvpHS);
  tree->Branch("pzvpHS", &event.pzvpHS);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCTdcCalibMan>("DCTDC")) &&
    (InitializeParameter<DCDriftParamMan>("DCDRFT")) &&
    (InitializeParameter<DCGeomMan>("DCGEO")) &&
    (InitializeParameter<FieldMan>("FLDMAP")) &&
    (InitializeParameter<TransferMatrixMan>("TM")) &&
    (InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
