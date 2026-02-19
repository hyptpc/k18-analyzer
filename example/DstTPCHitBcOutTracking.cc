// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "TPCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DstHelper.hh"
#include "HistTools.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCCluster.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#define RawHit 0
#define RawCluster 1

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gTpcParam = TPCParamMan::GetInstance();
const auto& gCounter  = debug::ObjectCounter::GetInstance();
const Double_t MaxChisqrBcOut = 5.0;
const Double_t MinResidualY = 0.0;
const Double_t MaxResidual = 20.0;
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTpcHit, kBcOut, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[TPCHit]", "[BcOut]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "bcout", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
std::vector<UInt_t> evnumPerFile;
Bool_t SetupReader();
}

//_____________________________________________________________________________
struct Event
{
  Int_t status;
  UInt_t runnum;
  UInt_t evnum;
  std::vector<Double_t> trigpat;
  std::vector<std::vector<Double_t>> trigflag;
  Int_t beamflag;
  std::vector<Double_t> clkTpc;
  std::vector<Double_t> cobo_id;

  Int_t nhTpc;
  std::vector<Double_t> raw_hitpos_x;
  std::vector<Double_t> raw_hitpos_y;
  std::vector<Double_t> raw_hitpos_z;
  std::vector<Double_t> raw_de;
  std::vector<Int_t> raw_padid;
  std::vector<Int_t> raw_layer;
  std::vector<Int_t> raw_row;

  Int_t nclTpc;
  std::vector<Double_t> cluster_x;
  std::vector<Double_t> cluster_y;
  std::vector<Double_t> cluster_z;
  std::vector<Double_t> cluster_de;
  std::vector<Int_t> cluster_size;
  std::vector<Int_t> cluster_layer;
  std::vector<Double_t> cluster_mrow;
  std::vector<Double_t> cluster_de_center;
  std::vector<Double_t> cluster_x_center;
  std::vector<Double_t> cluster_y_center;
  std::vector<Double_t> cluster_z_center;
  std::vector<Int_t> cluster_row_center;

  // BcOut Info
  Int_t ntBcOut;
  Int_t best_bcout_id;
  std::vector<Double_t> chisqrBcOut;
  std::vector<Double_t> x0BcOut;
  std::vector<Double_t> y0BcOut;
  std::vector<Double_t> u0BcOut;
  std::vector<Double_t> v0BcOut;

  // Residual Info (BcOut reference)
  std::vector<Double_t> residual_x_cluster, residual_y_cluster;
  std::vector<Double_t> residual_x_hit, residual_y_hit;

  void clearBasicInfo() {
    runnum   = 0;
    evnum    = 0;
    status   = 0;
    beamflag = beam::kUnknown;
    dst::clear_all(trigpat, trigflag, clkTpc, cobo_id);
  }

  void clearRawHits() {
    nhTpc = 0;
    dst::clear_all(raw_hitpos_x, raw_hitpos_y, raw_hitpos_z,
                       raw_de, raw_padid, raw_layer, raw_row);
  }

  void clearClusters() {
    nclTpc = 0;
    dst::clear_all(cluster_x, cluster_y, cluster_z, cluster_de,
                       cluster_size, cluster_layer, cluster_mrow,
                       cluster_de_center, cluster_x_center, cluster_y_center,
                       cluster_z_center, cluster_row_center);
  }

  void clearBcOut() {
    ntBcOut = 0;
    best_bcout_id = -1;
    dst::clear_all(
      chisqrBcOut, x0BcOut, y0BcOut, u0BcOut, v0BcOut,
      residual_x_cluster, residual_y_cluster,
      residual_x_hit, residual_y_hit
    );
  }

  void clear()
  {
    clearBasicInfo();
    clearRawHits();
    clearClusters();
    clearBcOut();
  }

};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<UInt_t>* runnum;
  TTreeReaderValue<UInt_t>* evnum;
  TTreeReaderValue<std::vector<Double_t>>* trigpat;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* trigflag;
  TTreeReaderValue<Int_t>* beamflag;
  TTreeReaderValue<Int_t>* npadTpc;   // number of pads
  TTreeReaderValue<Int_t>* nhTpc;     // number of hits
  // vector (size=nhTpc)
  TTreeReaderValue<std::vector<Int_t>>* layerTpc;     // layer id
  TTreeReaderValue<std::vector<Int_t>>* rowTpc;       // row id
  TTreeReaderValue<std::vector<Int_t>>* padTpc;       // pad id
  TTreeReaderValue<std::vector<Double_t>>* deTpc;     // dE
  TTreeReaderValue<std::vector<Double_t>>* tTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;    // clock time
  TTreeReaderValue<std::vector<Double_t>>* cobo_id;   // CoBo ID

  TTreeReaderValue<UInt_t>* evnum_bcout;
  TTreeReaderValue<Int_t>* ntBcOut;
  TTreeReaderValue<std::vector<Double_t>>* chisqrBcOut;
  TTreeReaderValue<std::vector<Double_t>>* x0BcOut;
  TTreeReaderValue<std::vector<Double_t>>* y0BcOut;
  TTreeReaderValue<std::vector<Double_t>>* u0BcOut;
  TTreeReaderValue<std::vector<Double_t>>* v0BcOut;
};

namespace root
{
  Event  event;
  Src    src;
  TTree *tree;
}

//_____________________________________________________________________________
int
main(int argc, char **argv)
{
  std::vector<std::string> arg(argv, argv+argc);

  if(!CheckArg(arg))
    return EXIT_FAILURE;
  if(!DstOpen(arg))
    return EXIT_FAILURE;
  if(!gConf.Initialize(arg[kConfFile]))
    return EXIT_FAILURE;
  if(!gConf.InitializeHistograms())
    return EXIT_FAILURE;
  if(!gConf.InitializeUnpacker())
    return EXIT_FAILURE;
  if(!dst::SetupReader())
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if(skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries(TTreeCont);
  if(max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();

  Int_t ievent = skip;
  for(; ievent<nevent && !CatchSignal::Stop(); ++ievent){
    gCounter.check();
    InitializeEvent();
    if(DstRead(ievent)) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6)
            << ievent << std::endl;

  DstClose();

  return EXIT_SUCCESS;
}

//_____________________________________________________________________________
Bool_t
dst::InitializeEvent()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstOpen(std::vector<std::string> arg)
{
  Int_t n_input_files = 0;
  for(const auto& name : TreeName) if(name != "") n_input_files++;

  Int_t open_file = 0;
  Int_t open_tree = 0;
  for(Int_t i=0; i<nArgc; ++i){
    if(TreeName[i] == "") continue;
    open_file += OpenFile(TFileCont[i], arg[i]);
    open_tree += OpenTree(TFileCont[i], TTreeCont[i], TreeName[i]);
  }

  if(open_file!=n_input_files || open_tree!=n_input_files){
    spdlog::error("DstOpen Failed: opened files/trees mismatch based on TreeName definitions."
                  " expected: {}, open_file: {}, open_tree: {}",
                  n_input_files, open_file, open_tree);
    return false;
  }
  if(!CheckEntries(TTreeCont))
    return false;

  TFileCont[kOutFile] = new TFile(arg[kOutFile].c_str(), "recreate");

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstRead(Int_t ievent)
{
  if(ievent%100==0){
    std::cout << "#D Event Number: "
              << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  // Check event numbers
  evnumPerFile = { **src.evnum, **src.evnum_bcout };
  if(!dst::CheckEventNumbers(evnumPerFile, ievent, {TreeName[kTpcHit], TreeName[kBcOut]})){
    return false;
  }

  event.runnum   = **src.runnum;
  event.evnum    = **src.evnum;
  event.trigpat  = **src.trigpat;
  event.trigflag = **src.trigflag;
  event.beamflag = **src.beamflag;
  event.clkTpc   = **src.clkTpc;
  event.cobo_id  = **src.cobo_id;
  HF1("Status", event.status++);

  if(**src.nhTpc == 0)
    return true;
  HF1("Status", event.status++);

  // Check CoBo clock size
  if(event.clkTpc.size() != NumOfSegCOBO){
    spdlog::warn("something is wrong: event.clkTpc.size() != {}", NumOfSegCOBO);
    return true;
  }
  std::vector<Int_t> bad_cobo;
  for(Int_t c = 0; c < NumOfSegCOBO; ++c){
    if(!std::isfinite(event.clkTpc[c])) bad_cobo.push_back(c);
  }
  if(!bad_cobo.empty()){
    std::ostringstream oss;
    for(size_t i = 0; i < bad_cobo.size(); ++i){ oss << (i ? "," : "") << bad_cobo[i]; }
    spdlog::warn("CoBo clock(s) missing (NaN/Inf): cobo={}, skip event", oss.str());
    return true;
  }
  HF1("Status", event.status++);

  // Re-calculate TPC hits
  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, **src.clkTpc);
  HF1("Status", event.status++);

  // check and select best BcOut
  if (**src.ntBcOut <= 0) return true;

  Double_t min_chi2 = 1.0e9;
  Int_t best_itr = -1;

  for (Int_t i_bcout = 0; i_bcout < **src.ntBcOut; ++i_bcout) {
    Double_t chi2 = (**src.chisqrBcOut)[i_bcout];
    if (chi2 < min_chi2) {
      min_chi2 = chi2;
      best_itr = i_bcout;
    }
  }
  event.best_bcout_id = best_itr;
  
  event.ntBcOut     = **src.ntBcOut;
  event.chisqrBcOut = **src.chisqrBcOut;
  event.x0BcOut     = **src.x0BcOut;
  event.y0BcOut     = **src.y0BcOut;
  event.u0BcOut     = **src.u0BcOut;
  event.v0BcOut     = **src.v0BcOut;

  // Tracking parameter
  Double_t u0 = TMath::QuietNaN();
  Double_t v0 = TMath::QuietNaN();
  Double_t x0 = TMath::QuietNaN();
  Double_t y0 = TMath::QuietNaN();

  if(event.best_bcout_id >= 0){
    u0 = (**src.u0BcOut)[event.best_bcout_id];
    v0 = (**src.v0BcOut)[event.best_bcout_id];
    x0 = (**src.x0BcOut)[event.best_bcout_id];
    y0 = (**src.y0BcOut)[event.best_bcout_id];
  }

  // TPC Hit Loop & Residuals
  Int_t nhTpc = 0;
  for (Int_t layer = 0; layer < NumOfLayersTPC; ++layer) {
    const auto& hc_hit = TPCAna.GetTPCHC(layer);
    for (const auto& hit : hc_hit) {
      if (!hit || !hit->IsGood()) continue;
      Int_t row = hit->GetRow();
      const auto& pos = hit->GetPosition();
      ThreeVector localPos(pos.X(), pos.Y(), pos.Z());
      ThreeVector globalPos = gGeom.Local2GlobalPos("HypTPC", localPos);

      // Raw Hit Info
#if RawHit
      event.raw_hitpos_x.push_back(pos.X());
      event.raw_hitpos_y.push_back(pos.Y());
      event.raw_hitpos_z.push_back(pos.Z());
      event.raw_de.push_back(hit->GetCDe());
      event.raw_padid.push_back(hit->GetPad());
      event.raw_layer.push_back(layer);
      event.raw_row.push_back(row);
#endif
      ++nhTpc;

      // Residual Calculation
      Double_t res_x = TMath::QuietNaN();
      Double_t res_y = TMath::QuietNaN();

      if(event.best_bcout_id >= 0){
        Double_t z_g = globalPos.z();
        Double_t x_bc = u0 * z_g + x0;
        Double_t y_bc = v0 * z_g + y0;
        res_x = globalPos.x() - x_bc;
        res_y = globalPos.y() - y_bc;

        HF1(Form("TPCHit_ResX_Layer%02d", layer), res_x);
        HF1("TPCHit_ResX", res_x);
        HF1(Form("TPCHit_ResY_Layer%02d", layer), res_y);
        HF1("TPCHit_ResY", res_y);
        HF2(Form("TPCHit_ResX_vs_X_Layer%02d", layer), x_bc, res_x);
        if (TMath::Abs(res_y) >= MinResidualY) {
          HF1(Form("TPCHit_ResY_Layer%02d_Row%03d", layer, row), res_y);
          HF2(Form("TPCHit_ResY_vs_Y_Layer%02d", layer), y_bc, res_y);
          HF2(Form("TPCHit_ResY_vs_Y_Layer%02d_Row%03d", layer, row), y_bc, res_y);
        }
        HF2("TPCHit_ResX_vs_Layer", layer, res_x);
        HF2("TPCHit_ResY_vs_Layer", layer, res_y);
        if (TMath::Abs(res_x) <= MaxResidual && TMath::Abs(res_y) <= MaxResidual) {
          HF2Poly("TPC_HitPat", hit->GetPad() + 1, 1.);
          HF2("TPCHit_Row_vs_Layer", layer, row);
        }

        Int_t cobo = tpc::GetCoBoId(layer, row);
        Int_t asad = tpc::GetASADId(layer, row);
        Bool_t cobo_valid = (cobo >= 0 && cobo < NumOfSegCOBO);
        Double_t res_y_noclk = TMath::QuietNaN(), res_y_raw = TMath::QuietNaN();
        if (!cobo_valid) {
          spdlog::warn("TPC Hit BcOut ResY vs ClockTime: invalid CoBo id (cobo={}) for layer={} row={}", cobo, layer, row);
        } else if (!std::isfinite((**src.clkTpc)[cobo])) {
          spdlog::warn("TPC Hit BcOut ResY vs ClockTime: non-finite clkTpc[{}]={} for layer={} row={}", cobo, (**src.clkTpc)[cobo], layer, row);
        } else {
          Double_t clk = (**src.clkTpc)[cobo];
          Double_t cclk = 0;
          gTpcParam.GetCClock(layer, row, clk, cclk);
          Double_t ctime_noclk = hit->GetCTime() - cclk;
          Double_t y_noclk = 0, y_raw = 0;
          gTpcParam.GetDriftLength(layer, row, ctime_noclk, y_noclk);
          gTpcParam.GetDriftLength(layer, row, ctime_noclk + clk, y_raw);
          ThreeVector local_noclk(localPos.X(), y_noclk, localPos.Z());
          ThreeVector local_raw(localPos.X(), y_raw, localPos.Z());
          ThreeVector global_noclk = gGeom.Local2GlobalPos("HypTPC", local_noclk);
          ThreeVector global_raw = gGeom.Local2GlobalPos("HypTPC", local_raw);
          res_y_noclk = global_noclk.y() - (v0 * global_noclk.z() + y0);
          res_y_raw   = global_raw.y()  - (v0 * global_raw.z() + y0);

          HF2(Form("TPCHit_ResY_vs_ClockTime_CoBo%d", cobo), clk, res_y);
          HF2(Form("TPCHit_ResY_vs_ClockTime_CoBo%d_RawClock", cobo), clk, res_y_raw);
#ifdef DEBUG_COBO_CLOCK
          HF2(Form("TPCHit_ResY_vs_ClockTime_CoBo%d_NoClock", cobo), clk, res_y_noclk);
#endif
        } // else (cobo_valid)

        Bool_t asad_valid = (asad >= 0 && asad < NumOfAsadTPC);
        if (!asad_valid) {
          spdlog::warn("TPC Hit BcOut ResY vs ClockTime: invalid Asad id (asad={}) for layer={} row={}", asad, layer, row);
        } else if (cobo_valid && std::isfinite((**src.clkTpc)[cobo])) {
          Double_t clk = (**src.clkTpc)[cobo];
          HF2(Form("TPCHit_ResY_vs_ClockTime_Asad%02d", asad), clk, res_y);
          HF2(Form("TPCHit_ResY_vs_ClockTime_Asad%02d_RawClock", asad), clk, res_y_raw);
#ifdef DEBUG_COBO_CLOCK
          HF2(Form("TPCHit_ResY_vs_ClockTime_Asad%02d_NoClock", asad), clk, res_y_noclk);
#endif
        } // else if (cobo_valid && asad_valid)
      } // if(best_bcout_id >= 0)

      event.residual_x_hit.push_back(res_x);
      event.residual_y_hit.push_back(res_y);
    } // for(hit)
  } // for(layer)
  event.nhTpc = nhTpc;

  // TPC Cluster Loop & Residuals
  Int_t nclTpc = 0;
  for (Int_t layer = 0; layer < NumOfLayersTPC; ++layer) {
    const auto& hc_cl = TPCAna.GetTPCClCont(layer);
    for (const auto& cl : hc_cl) {
      if (!cl || !cl->IsGood()) continue;
      Int_t centerRow = cl->GetCenterHit()->GetRow();
      ThreeVector localPos(cl->GetX(), cl->GetY(), cl->GetZ());
      ThreeVector globalPos = gGeom.Local2GlobalPos("HypTPC", localPos);

      // Raw Cluster Info
#if RawCluster
      event.cluster_x.push_back(localPos.X());
      event.cluster_y.push_back(localPos.Y());
      event.cluster_z.push_back(localPos.Z());
      event.cluster_de.push_back(cl->GetDe());
      event.cluster_size.push_back(cl->GetClusterSize());
      event.cluster_layer.push_back(layer);
      event.cluster_mrow.push_back(cl->MeanRow()); // Assuming MeanRow() is correct
      event.cluster_row_center.push_back(centerRow);
      
      TPCHit* centerHit = cl->GetCenterHit();
      const TVector3& centerPos = centerHit->GetPosition();
      event.cluster_de_center.push_back(centerHit->GetCDe());
      event.cluster_x_center.push_back(centerPos.X());
      event.cluster_y_center.push_back(centerPos.Y());
      event.cluster_z_center.push_back(centerPos.Z());
#endif
      ++nclTpc;

      // Residual Calculation
      Double_t res_x = TMath::QuietNaN();
      Double_t res_y = TMath::QuietNaN();

      if(event.best_bcout_id >= 0){
        Double_t z_g = globalPos.z();
        Double_t x_bc = u0 * z_g + x0;
        Double_t y_bc = v0 * z_g + y0;
        res_x = globalPos.x() - x_bc;
        res_y = globalPos.y() - y_bc;

        HF1(Form("TPCCl_ResX_Layer%02d", layer), res_x);
        HF1("TPCCl_ResX", res_x);
        HF1(Form("TPCCl_ResY_Layer%02d", layer), res_y);
        HF1("TPCCl_ResY", res_y);
        HF2(Form("TPCCl_ResX_vs_X_Layer%02d", layer), globalPos.x(), res_x);
        if (TMath::Abs(res_y) >= MinResidualY) {
          HF1(Form("TPCCl_ResY_Layer%02d_Row%03d", layer, centerRow), res_y);
          HF2(Form("TPCCl_ResY_vs_Y_Layer%02d", layer), y_bc, res_y);
          HF2(Form("TPCCl_ResY_vs_Y_Layer%02d_Row%03d", layer, centerRow), y_bc, res_y);
        }
        HF2("TPCCl_ResX_vs_Layer", layer, res_x);
        HF2("TPCCl_ResY_vs_Layer", layer, res_y);
        if (TMath::Abs(res_x) <= MaxResidual && TMath::Abs(res_y) <= MaxResidual) { 
          HF2Poly("TPC_Cluster_HitPat", cl->GetCenterHit()->GetPad() + 1, 1.);
          HF2("TPCCl_Row_vs_Layer", layer, centerRow);
        }

        Int_t cobo = tpc::GetCoBoId(layer, centerRow);
        Int_t asad = tpc::GetASADId(layer, centerRow);
        Bool_t cobo_valid = (cobo >= 0 && cobo < NumOfSegCOBO);
        Double_t res_y_noclk = TMath::QuietNaN(), res_y_raw = TMath::QuietNaN();
        if (!cobo_valid) {
          spdlog::warn("TPC Cluster BcOut ResY vs ClockTime: invalid CoBo id (cobo={}) for layer={} row={}", cobo, layer, centerRow);
        } else if (!std::isfinite((**src.clkTpc)[cobo])) {
          spdlog::warn("TPC Cluster BcOut ResY vs ClockTime: non-finite clkTpc[{}]={} for layer={} row={}", cobo, (**src.clkTpc)[cobo], layer, centerRow);
        } else {
          Double_t clk = (**src.clkTpc)[cobo];
          Double_t cclk = 0;
          gTpcParam.GetCClock(layer, centerRow, clk, cclk);
          Double_t ctime_noclk = cl->GetCenterHit()->GetCTime() - cclk;
          Double_t y_noclk = 0, y_raw = 0;
          gTpcParam.GetDriftLength(layer, centerRow, ctime_noclk, y_noclk);
          gTpcParam.GetDriftLength(layer, centerRow, ctime_noclk + clk, y_raw);
          ThreeVector local_noclk(localPos.X(), y_noclk, localPos.Z());
          ThreeVector local_raw(localPos.X(), y_raw, localPos.Z());
          ThreeVector global_noclk = gGeom.Local2GlobalPos("HypTPC", local_noclk);
          ThreeVector global_raw = gGeom.Local2GlobalPos("HypTPC", local_raw);
          res_y_noclk = global_noclk.y() - (v0 * global_noclk.z() + y0);
          res_y_raw   = global_raw.y()  - (v0 * global_raw.z() + y0);

          HF2(Form("TPCCl_ResY_vs_ClockTime_CoBo%d", cobo), clk, res_y);
          HF2(Form("TPCCl_ResY_vs_ClockTime_CoBo%d_RawClock", cobo), clk, res_y_raw);
#ifdef DEBUG_COBO_CLOCK
          HF2(Form("TPCCl_ResY_vs_ClockTime_CoBo%d_NoClock", cobo), clk, res_y_noclk);
#endif
        } // else (cobo_valid)

        Bool_t asad_valid = (asad >= 0 && asad < NumOfAsadTPC);
        if (!asad_valid) {
          spdlog::warn("TPC Cluster BcOut ResY vs ClockTime: invalid Asad id (asad={}) for layer={} row={}", asad, layer, centerRow);
        } else if (cobo_valid && std::isfinite((**src.clkTpc)[cobo])) {
          Double_t clk = (**src.clkTpc)[cobo];
          HF2(Form("TPCCl_ResY_vs_ClockTime_Asad%02d", asad), clk, res_y);
          HF2(Form("TPCCl_ResY_vs_ClockTime_Asad%02d_RawClock", asad), clk, res_y_raw);
#ifdef DEBUG_COBO_CLOCK
          HF2(Form("TPCCl_ResY_vs_ClockTime_Asad%02d_NoClock", asad), clk, res_y_noclk);
#endif
        } // else if (cobo_valid && asad_valid)
      } // if(best_bcout_id >= 0)

      event.residual_x_cluster.push_back(res_x);
      event.residual_y_cluster.push_back(res_y);
    } // for(cl)
  } // for(layer)
  event.nclTpc = nclTpc;
  HF1("Status", event.status++);

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstClose()
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  const Int_t n = TFileCont.size();
  for(Int_t i=0; i<n; ++i){
    if(TTreeReaderCont[i]) delete TTreeReaderCont[i];
    if(TTreeCont[i])       delete TTreeCont[i];
    if(TFileCont[i])       delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::SetupReader()
{
  // -------------------------------------------------------
  // TPC Hit
  // -------------------------------------------------------
  if (!dst::SetupReader(kTpcHit, "kTpcHit")) return false;

  dst::SetBranch(TTreeReaderCont[kTpcHit], "run_number",   src.runnum);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "event_number", src.evnum);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "trig_pat",     src.trigpat);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "trig_flag",    src.trigflag);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "beam_flag",    src.beamflag);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "npadTpc",      src.npadTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "nhTpc",        src.nhTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "layerTpc",     src.layerTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "rowTpc",       src.rowTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "padTpc",       src.padTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "deTpc",        src.deTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "tTpc",         src.tTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "clkTpc",       src.clkTpc);
  dst::SetBranch(TTreeReaderCont[kTpcHit], "cobo_id",      src.cobo_id);

  // -------------------------------------------------------
  // BcOut
  // -------------------------------------------------------
  if (!dst::SetupReader(kBcOut, "kBcOut")) return false;

  dst::SetBranch(TTreeReaderCont[kBcOut], "event_number", src.evnum_bcout);
  dst::SetBranch(TTreeReaderCont[kBcOut], "ntrack", src.ntBcOut);
  dst::SetBranch(TTreeReaderCont[kBcOut], "chisqr", src.chisqrBcOut);
  dst::SetBranch(TTreeReaderCont[kBcOut], "x0",     src.x0BcOut);
  dst::SetBranch(TTreeReaderCont[kBcOut], "y0",     src.y0BcOut);
  dst::SetBranch(TTreeReaderCont[kBcOut], "u0",     src.u0BcOut);
  dst::SetBranch(TTreeReaderCont[kBcOut], "v0",     src.v0BcOut);

  evnumPerFile.resize(2);
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  hist::BuildStatus();
  hist::BuildTPCHitBcOutTracking();

  tree = new TTree("tpc", "tree of DstTPCTracking");
  tree->Branch("status", &event.status);
  tree->Branch("run_number", &event.runnum);
  tree->Branch("event_number", &event.evnum);
  tree->Branch("trig_pat", &event.trigpat);
  tree->Branch("trig_flag", &event.trigflag);
  tree->Branch("beam_flag", &event.beamflag);
  tree->Branch("clkTpc", &event.clkTpc);
  tree->Branch("cobo_id", &event.cobo_id);

  // BcOut Info
  tree->Branch("ntBcOut", &event.ntBcOut);
  tree->Branch("best_bcout_id", &event.best_bcout_id);
  tree->Branch("chisqrBcOut", &event.chisqrBcOut);
  tree->Branch("x0BcOut", &event.x0BcOut);
  tree->Branch("y0BcOut", &event.y0BcOut);
  tree->Branch("u0BcOut", &event.u0BcOut);
  tree->Branch("v0BcOut", &event.v0BcOut);

  // Residual Info
  tree->Branch("residual_x_cluster", &event.residual_x_cluster);
  tree->Branch("residual_y_cluster", &event.residual_y_cluster);
  tree->Branch("residual_x_hit", &event.residual_x_hit);
  tree->Branch("residual_y_hit", &event.residual_y_hit);

#if RawHit
  tree->Branch("nhTpc", &event.nhTpc);
  tree->Branch("raw_hitpos_x", &event.raw_hitpos_x);
  tree->Branch("raw_hitpos_y", &event.raw_hitpos_y);
  tree->Branch("raw_hitpos_z", &event.raw_hitpos_z);
  tree->Branch("raw_de", &event.raw_de);
  tree->Branch("raw_padid", &event.raw_padid);
  tree->Branch("raw_layer", &event.raw_layer);
  tree->Branch("raw_row", &event.raw_row);
#endif

#if RawCluster
  tree->Branch("nclTpc", &event.nclTpc);
  tree->Branch("cluster_x", &event.cluster_x);
  tree->Branch("cluster_y", &event.cluster_y);
  tree->Branch("cluster_z", &event.cluster_z);
  tree->Branch("cluster_de", &event.cluster_de);
  tree->Branch("cluster_size", &event.cluster_size);
  tree->Branch("cluster_layer", &event.cluster_layer);
  tree->Branch("cluster_row_center", &event.cluster_row_center);
  tree->Branch("cluster_mrow", &event.cluster_mrow);
  tree->Branch("cluster_de_center", &event.cluster_de_center);
  tree->Branch("cluster_x_center", &event.cluster_x_center);
  tree->Branch("cluster_y_center", &event.cluster_y_center);
  tree->Branch("cluster_z_center", &event.cluster_z_center);
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO") &&
     InitializeParameter<TPCParamMan>("TPCPRM", "TPCPHASE") &&
     InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
