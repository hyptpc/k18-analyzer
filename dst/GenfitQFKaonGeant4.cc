// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <TLorentzVector.h>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "FieldMan.hh"
#include "DatabasePDG.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "DCAnalyzer.hh"
#include "TPCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DstHelper.hh"
#include "HodoPHCMan.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCCluster.hh"
#include "TPCVertex.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"
#include "KEKKinematicFit.hh"

#include "HypTPCFitter.hh"
#include "HypTPCTask.hh"

#define SaveRawData 1
#define DebugDisp 0
#define SaveHistograms 0
#define kkevent 0
#define XiRecon 1
#define LLRecon 0

namespace
{
  using namespace root;
  using namespace dst;
  using hddaq::unpacker::GUnpacker;
  const auto qnan = TMath::QuietNaN();
  const auto& gUnpacker = GUnpacker::get_instance();
  const auto& gUser = UserParamMan::GetInstance();
  auto&       gConf = ConfMan::GetInstance();
  const auto& gGeom = DCGeomMan::GetInstance();
  const auto& gCounter = debug::ObjectCounter::GetInstance();
  const auto& gTPC  = TPCParamMan::GetInstance();
  const Int_t MaxTPCHits = 10000;

  static TRandom3 rand_mass(0);
  double xi_res_smear = 0.00345212;  double l_res_smear = 0.00263642;

  //For GenFit Setting
  const bool Const_field = false; //Must be false for linear tracking
  const Int_t verbosity = 0;//0~3;
  //const Int_t verbosity = 1;
  const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");

  const Double_t vtx_scan_range = 150.; //ref
  const Double_t vtx_scan_rangeInsideL = 50.;
  const Double_t vtx_scan_rangeInsidePi = 50.;

  const Double_t xi_masscut = 0.1; const Double_t lambda_masscut = 0.1; //final
  //const Double_t xi_masscut = 0.15; const Double_t lambda_masscut = 0.1; //ref
  const Double_t p_vtx_distcut = 300;
  const Double_t pi_vtx_distcut = 300;
  const Double_t pi2_vtx_distcut = 300;
  const Double_t e_vtx_distcut = 300;
  const Double_t pipi_distcut = 10.; //ref
  const Double_t ppi_distcut = 10.; //ref
  //const Double_t lpi_distcut = 10.;
  const Double_t lpi_distcut = 15.; //ref
  const Double_t phi_kk_distcut = 15.; //ref
  const Double_t xitarget_distcut = 50.; //ref
  const Double_t ltarget_distcut = 25.;

  const Double_t GFppi_distcut = 10.;
  const Double_t GFlpi_distcut = 10.;
  //const Double_t GFlpi_distcut = 15.;
  const Double_t GFphi_kk_distcut = 10.; //ref
  const Double_t GFltarget_distcut = 25.;
  const Double_t GFxitarget_ycut = 20.;
  const Double_t GFltarget_ycut = 20.;

  const Double_t residual_track_distcut = 25.;

  //For gamma reconstruction
  const Double_t gammatarget_distcut = 25.; //temp
  //const Double_t gamma_ecut = 0.5; //for Xi* decay pi0 has momentum < 0.35 GeV/c
  const Double_t gamma_ecut = 2.5; //for Xi* decay pi0 has momentum < 0.35 GeV/c
  const Double_t epem_distcut = 10.; //temp

  //Measured resolutions for multi-track vertexing
  const Double_t duCh2 = 0.001381;
  const Double_t dvCh2 = duCh2;
  const Double_t duDiamond = 0.2796;
  const Double_t dvDiamond = duDiamond;

  Double_t res_uK18 = 0.00288/sqrt(2);
  Double_t res_xK18 = 1.432/sqrt(2);
  Double_t res_vK18 = 0.00334/sqrt(2);
  Double_t res_yK18 = 2.836/sqrt(2);
  Double_t res_uKurama = hypot(res_uK18,duCh2);
  Double_t res_xKurama = res_xK18;
  Double_t res_vKurama = hypot(res_vK18,dvCh2);
  Double_t res_yKurama = res_yK18;
  Double_t res_xXiVtx = 0.6*2;
  Double_t res_yXiVtx = 0.5*2;
  Double_t res_xLdVtx = 0.5*2;
  Double_t res_yLdVtx = 0.5*2;

  const Double_t& HS_field_0 = ConfMan::Get<Double_t>("HSFLDCALIB");
  const Double_t& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
  const Double_t& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");

  const Int_t& CarbonRun = ConfMan::Get<Int_t>("CTARGET");  

  const Double_t HTOF_TIME_RESOLUTION_SIGMA = 0.200;

  const double minMM = 0.38;
  const double maxMM = 0.58;
  // const double minM2Km =  0.01;
  // const double maxM2Km =  0.50;

  const double minThetaKP = 3.5;
  const double maxThetaKP = 4.5;  
  
  const double dedx_smear_nsigma = 1.0;
  const double m2_smear = 0.05;
  // const double dedx_smear_nsigma = 0.0001;
  // const double m2_smear = 0.0001;

  const Double_t sigma_dedx_k[5] = {6.24543, -3.21037, 1.52683, 127.099, -9.1004};
  const Double_t conversion_factor = 12171.3; //HypTPC's ADC to <dE/dx>
  
  int TPCToG4TrackID(std::vector<TVector3>TPCHit, int nhG4,int* tidG4, double* xG4,double* yG4,double* zG4 ,int& nhits){
    std::vector<TVector3> G4Hits;
    for(int ih=0;ih<nhG4;++ih){
      TVector3 G4Hit(xG4[ih],yG4[ih],zG4[ih]);
      G4Hits.push_back(G4Hit);
    }
    int MaxTracks = 1000;
    TH1I Counter("counter","counter",MaxTracks,0,MaxTracks);
    for(auto hit:TPCHit){
      double dl = 5000;
      int G4ID = -1;
      for(int ih=0;ih<nhG4;++ih){
        auto G4Hit = G4Hits.at(ih);
        double dist = (G4Hit - hit).Mag();
        if(dist < dl){
          dl = dist;
          G4ID = tidG4[ih];
        }
      }
      Counter.Fill(G4ID);
    }
    nhits = Counter.GetMaximum();
    int G4id = Counter.GetMaximumBin()-1;
    return G4id;
  }
  TVector3 GetG4Mom(TVector3 TPCHit, vector<TVector3> G4Hits,vector<TVector3>G4Moms){
    int nh = G4Hits.size();
    double dl = 5000;
    TVector3 mom;
    for(int ih=0;ih<nh;++ih){
      auto G4Hit = G4Hits.at(ih);
      double dist = (G4Hit - TPCHit).Mag();
      if(dist < dl){
        dl = dist;
        mom = G4Moms.at(ih);
      }
    }
    return mom;
  }
  int CountHits(int id, double* posx,double* posz,int* tidG4, int nh){
    int count = 0;
    for(int ih=0;ih<nh;++ih){
      double x = posx[ih],z=posz[ih];
      int pad = tpc::findPadID(z,x);
      int layer = tpc::getLayerID(pad);
      int row = tpc::getRowID(pad);
      double val = 0;
      gTPC.GetCDe(layer,row,1,val);
      if(val==0) continue;
      if(tidG4[ih] == id) count++;
    }
    return count;
  }
  TLorentzVector ToHelix(TLorentzVector GlobalLV){
    double E = GlobalLV.E();
    double X = -GlobalLV.X();
    double Y = GlobalLV.Z();
    double Z = GlobalLV.Y();
    return TLorentzVector(X,Y,Z,E);
  }
  TLorentzVector ToGlobal(TLorentzVector HelixLV){
    return ToHelix(HelixLV);
  }

  Double_t HypTPCdEdx1sigmaKaon(Double_t poq){
    Double_t mk = 493.677;
    Double_t par_k[2] = {conversion_factor, mk};
    Double_t dedx_k = Kinematics::HypTPCBethe(&poq, par_k); //P10's <dE/dx>_k
    // 1 sigma of <dE/dx>_k
    Double_t sigma_k = (sigma_dedx_k[0] + sigma_dedx_k[1]*TMath::Abs(poq) +
			sigma_dedx_k[2]*poq*poq + sigma_dedx_k[3]*TMath::Exp(sigma_dedx_k[4]*TMath::Abs(poq)));
    return sigma_k;
  }

  int GetBinIndexBEk(double x){
    double xmin=-0.3;
    double xmax=0.3;
    double nbins=120.;
    
    if(x<xmin||x>xmax) return -1;
    double binwidth = (xmax-xmin)/nbins;    
    int bin = static_cast<int>(std::floor((x-xmin)/binwidth))+1;
    if(bin>nbins) bin = nbins;
    return bin;
  }
  
}

namespace dst
{
  enum kArgc
    {
      kProcess, kConfFile,
      kHelixTrackingGeant4, kOutFile, nArgc
    };
  std::vector<TString> ArgName =
  { "[Process]", "[ConfFile]", "[DstTPCTrackingHelixgeant4]", "[OutFile]" };
  std::vector<TString> TreeName = { "", "", "tpc", "" };
  std::vector<TFile*> TFileCont;
  std::vector<TTree*> TTreeCont;
  std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________________
struct Event
{
  Int_t status;
  Int_t evnum;
  std::vector<Int_t> trigpat;
  std::vector<Int_t> trigflag;

  Int_t nhHtof;
  std::vector<Double_t> HtofSeg;
  std::vector<Double_t> tHtof;
  std::vector<Double_t> dtHtof;
  std::vector<Double_t> deHtof;
  std::vector<Double_t> posHtof;
  std::vector<Int_t> G4tidHtof;

  Int_t nhFtof;
  std::vector<Double_t> FtofSeg;
  std::vector<Double_t> tFtof;
  std::vector<Double_t> deFtof;
  std::vector<Double_t> posFtof;

  Int_t ntK18;
  std::vector<Double_t> pK18;
  std::vector<Double_t> chisqrK18;
  std::vector<Double_t> xtgtK18;
  std::vector<Double_t> ytgtK18;
  std::vector<Double_t> utgtK18;
  std::vector<Double_t> vtgtK18;

  Int_t ntKurama;
  std::vector<Double_t> chisqrKurama;
  std::vector<Double_t> pKurama;
  std::vector<Double_t> qKurama;
  std::vector<Double_t> m2;
  std::vector<Double_t> thetaKurama;
  std::vector<Double_t> xtgtKurama;
  std::vector<Double_t> ytgtKurama;
  std::vector<Double_t> utgtKurama;
  std::vector<Double_t> vtgtKurama;
  std::vector<Double_t> pathwcKurama;
  std::vector<Double_t> xin;
  std::vector<Double_t> yin;
  std::vector<Double_t> zin;
  std::vector<Double_t> pxin;
  std::vector<Double_t> pyin;
  std::vector<Double_t> pzin;
  std::vector<Double_t> xout;
  std::vector<Double_t> yout;
  std::vector<Double_t> zout;
  std::vector<Double_t> pxout;
  std::vector<Double_t> pyout;
  std::vector<Double_t> pzout;

  Int_t nKm;
  Int_t nKp;
  Int_t nKK;
  std::vector<Int_t> inside;
  std::vector<Double_t> vtx;
  std::vector<Double_t> vty;
  std::vector<Double_t> vtz;
  std::vector<Double_t> closeDist;
  std::vector<Double_t> MissMass;
  std::vector<Double_t> MissMassCorr;
  std::vector<Double_t> MissMassCorrDE;
  std::vector<Double_t> MissMassNucl;
  std::vector<Double_t> MissMassNuclCorr;
  std::vector<Double_t> MissMassNuclCorrDE;  
  std::vector<Double_t> pOrg;
  std::vector<Double_t> pCalc;
  std::vector<Double_t> pCorr;
  std::vector<Double_t> pCorrDE;
  std::vector<Double_t> xb;
  std::vector<Double_t> yb;
  std::vector<Double_t> ub;
  std::vector<Double_t> vb;
  std::vector<Double_t> xs;
  std::vector<Double_t> ys;
  std::vector<Double_t> us;
  std::vector<Double_t> vs;
  std::vector<Int_t> Kflag;

  //TPC RK
  std::vector<Int_t> isgoodTPCK18;
  std::vector<Double_t> chisqrTPCK18;
  std::vector<Double_t> pTPCK18;
  std::vector<Double_t> qTPCK18;
  std::vector<Double_t> xtgtTPCK18;
  std::vector<Double_t> ytgtTPCK18;
  std::vector<Double_t> utgtTPCK18;
  std::vector<Double_t> vtgtTPCK18;
  std::vector<Double_t> thetaTPCK18;
  std::vector<Double_t> lhtofTPCK18;
  std::vector<Double_t> xhtofTPCK18;
  std::vector<Double_t> yhtofTPCK18;
  std::vector<std::vector<Double_t>> lvpTPCK18;
  std::vector<std::vector<Double_t>> xvpTPCK18;
  std::vector<std::vector<Double_t>> yvpTPCK18;

  std::vector<Int_t> isgoodTPCKurama;
  std::vector<Int_t> kflagTPCKurama;
  std::vector<Double_t> chisqrTPCKurama;
  std::vector<Double_t> pTPCKurama;
  std::vector<Double_t> qTPCKurama;
  std::vector<Double_t> m2TPCKurama;
  std::vector<Double_t> xtgtTPCKurama;
  std::vector<Double_t> ytgtTPCKurama;
  std::vector<Double_t> utgtTPCKurama;
  std::vector<Double_t> vtgtTPCKurama;
  std::vector<Double_t> thetaTPCKurama;
  std::vector<Double_t> pathTPCKurama;
  std::vector<Double_t> lhtofTPCKurama;
  std::vector<Double_t> xhtofTPCKurama;
  std::vector<Double_t> yhtofTPCKurama;
  std::vector<std::vector<Double_t>> lvpTPCKurama;
  std::vector<std::vector<Double_t>> xvpTPCKurama;
  std::vector<std::vector<Double_t>> yvpTPCKurama;

  std::vector<Int_t> isgoodTPC;
  std::vector<Int_t> insideTPC;
  std::vector<Double_t> vtxTPC;
  std::vector<Double_t> vtyTPC;
  std::vector<Double_t> vtzTPC;
  std::vector<Double_t> closeDistTPC;
  std::vector<Double_t> MissMassTPC;
  std::vector<Double_t> MissMassCorrTPC;
  std::vector<Double_t> MissMassCorrDETPC;
  std::vector<Double_t> pOrgTPC;
  std::vector<Double_t> pCalcTPC;
  std::vector<Double_t> pCorrTPC;
  std::vector<Double_t> pCorrDETPC;
  std::vector<Double_t> thetaTPC;
  std::vector<Double_t> thetaCMTPC;
  std::vector<Double_t> costCMTPC;
  std::vector<Double_t> xbTPC;
  std::vector<Double_t> ybTPC;
  std::vector<Double_t> ubTPC;
  std::vector<Double_t> vbTPC;
  std::vector<Double_t> xsTPC;
  std::vector<Double_t> ysTPC;
  std::vector<Double_t> usTPC;
  std::vector<Double_t> vsTPC;

  std::vector<Double_t> BE;
  std::vector<Double_t> BETPC;
  std::vector<Double_t> BE_LL;
  std::vector<Double_t> BETPC_LL;
  std::vector<Double_t> km_mom_x;
  std::vector<Double_t> km_mom_y;
  std::vector<Double_t> km_mom_z;
  std::vector<Double_t> kp_mom_x;
  std::vector<Double_t> kp_mom_y;
  std::vector<Double_t> kp_mom_z;


  Int_t nclTpc; // Number of clusters
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
  std::vector<Int_t> cluster_houghflag;
  std::vector<Int_t> cluster_G4tid;
  std::vector<Int_t> cluster_G4protonid;

  Int_t remain_nclTpc; // Number of remain clusters not occupied in the tracks
  std::vector<Double_t> remain_cluster_x;
  std::vector<Double_t> remain_cluster_y;
  std::vector<Double_t> remain_cluster_z;
  std::vector<Double_t> remain_cluster_de;
  std::vector<Int_t> remain_cluster_size;
  std::vector<Int_t> remain_cluster_layer;
  std::vector<Double_t> remain_cluster_mrow;
  std::vector<Double_t> remain_cluster_de_center;
  std::vector<Double_t> remain_cluster_x_center;
  std::vector<Double_t> remain_cluster_y_center;
  std::vector<Double_t> remain_cluster_z_center;
  std::vector<Int_t> remain_cluster_row_center;
  std::vector<Int_t> remain_cluster_houghflag;
  std::vector<Int_t> remain_cluster_G4tid;
  std::vector<Int_t> remain_cluster_G4protonid;

  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> trackid; //for Kurama K1.8 tracks
  std::vector<Int_t> isXi;
  std::vector<Int_t> isBeam;
  std::vector<Int_t> isKurama;
  std::vector<Int_t> isK18;
  std::vector<Int_t> isAccidental;
  std::vector<Int_t> isMultiloop;
  std::vector<Int_t> isInTarget;
  std::vector<Int_t> charge; //Helix charge
  std::vector<Int_t> pid;
  std::vector<Double_t> chisqr;
  std::vector<Double_t> pval;
  std::vector<Double_t> purity;
  std::vector<Double_t> efficiency;
  std::vector<Int_t> G4tid;
  std::vector<Int_t> G4pid;
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx; //reference dedx
  std::vector<Double_t> mom0; //Helix momentum at Y = 0
  std::vector<Double_t> path; //Helix path
  std::vector<Int_t> isElectron;
  std::vector<Double_t> nsigma_triton;
  std::vector<Double_t> nsigma_deutron;
  std::vector<Double_t> nsigma_proton;
  std::vector<Double_t> nsigma_kaon;
  std::vector<Double_t> nsigma_pion;
  std::vector<Double_t> nsigma_electron;

  std::vector<std::vector<Double_t>> hitlayer;
  std::vector<std::vector<Double_t>> hitpos_x;
  std::vector<std::vector<Double_t>> hitpos_y;
  std::vector<std::vector<Double_t>> hitpos_z;
  std::vector<std::vector<Double_t>> calpos_x;
  std::vector<std::vector<Double_t>> calpos_y;
  std::vector<std::vector<Double_t>> calpos_z;
  std::vector<std::vector<Double_t>> mom_x;
  std::vector<std::vector<Double_t>> mom_y;
  std::vector<std::vector<Double_t>> mom_z;
  std::vector<std::vector<Double_t>> residual;
  std::vector<std::vector<Double_t>> residual_x;
  std::vector<std::vector<Double_t>> residual_y;
  std::vector<std::vector<Double_t>> residual_z;
  std::vector<std::vector<Double_t>> resolution_x;
  std::vector<std::vector<Double_t>> resolution_y;
  std::vector<std::vector<Double_t>> resolution_z;
  std::vector<std::vector<Double_t>> helix_t;
  std::vector<std::vector<Double_t>> pathhit;
  std::vector<std::vector<Double_t>> alpha;
  std::vector<std::vector<Double_t>> track_cluster_de;
  std::vector<std::vector<Double_t>> track_cluster_size;
  std::vector<std::vector<Double_t>> track_cluster_mrow;
  std::vector<std::vector<Double_t>> track_cluster_de_center;
  std::vector<std::vector<Double_t>> track_cluster_x_center;
  std::vector<std::vector<Double_t>> track_cluster_y_center;
  std::vector<std::vector<Double_t>> track_cluster_z_center;
  std::vector<std::vector<Double_t>> track_cluster_row_center;

  std::vector<Int_t> chargeIndistinguishable;
  std::vector<Double_t> chisqr_inverted;
  std::vector<Double_t> pval_inverted;
  std::vector<Double_t> helix_cx_inverted;
  std::vector<Double_t> helix_cy_inverted;
  std::vector<Double_t> helix_z0_inverted;
  std::vector<Double_t> helix_r_inverted;
  std::vector<Double_t> helix_dz_inverted;
  std::vector<Double_t> mom0_inverted;//Helix momentum at Y = 0
  std::vector<Int_t> pid_inverted;

  Int_t nvtxTpc;
  std::vector<Double_t> vtx_x;
  std::vector<Double_t> vtx_y;
  std::vector<Double_t> vtx_z;
  std::vector<Double_t> vtx_dist;
  std::vector<Double_t> vtx_angle;
  std::vector<std::vector<Double_t>> vtxid;
  std::vector<std::vector<Double_t>> vtxmom_theta;
  std::vector<std::vector<Double_t>> vtxpos_x;
  std::vector<std::vector<Double_t>> vtxpos_y;
  std::vector<std::vector<Double_t>> vtxpos_z;
  std::vector<std::vector<Double_t>> vtxmom_x;
  std::vector<std::vector<Double_t>> vtxmom_y;
  std::vector<std::vector<Double_t>> vtxmom_z;

  std::vector<Int_t> isLambda;
  std::vector<Int_t> ncombiLambda;
  std::vector<Double_t> distLambda;
  std::vector<Double_t> angleLambda;
  std::vector<Double_t> bestmassLambda;
  std::vector<std::vector<Double_t>> massLambda;
  std::vector<std::vector<Double_t>> vtxLambda_x;
  std::vector<std::vector<Double_t>> vtxLambda_y;
  std::vector<std::vector<Double_t>> vtxLambda_z;
  std::vector<std::vector<Double_t>> momLambda;
  std::vector<std::vector<Double_t>> momLambda_x;
  std::vector<std::vector<Double_t>> momLambda_y;
  std::vector<std::vector<Double_t>> momLambda_z;
  std::vector<std::vector<Double_t>> decaysidLambda;
  std::vector<std::vector<Double_t>> decaysmomLambda;
  std::vector<std::vector<Double_t>> decaysmomLambda_x;
  std::vector<std::vector<Double_t>> decaysmomLambda_y;
  std::vector<std::vector<Double_t>> decaysmomLambda_z;

  Int_t nvtxTpcClustered;
  std::vector<Double_t> clusteredVtx_x;
  std::vector<Double_t> clusteredVtx_y;
  std::vector<Double_t> clusteredVtx_z;
  std::vector<std::vector<Double_t>> clusteredVtxid;
  
  // Int_t ncombiLreconfailed;
  // std::vector<Int_t> pidLreconfailed;
  // std::vector<Int_t> piidLreconfailed;
  // std::vector<Double_t> LdecayvtxLreconfailed_x;
  // std::vector<Double_t> LdecayvtxLreconfailed_y;
  // std::vector<Double_t> LdecayvtxLreconfailed_z;
  // std::vector<Double_t> LmassLreconfailed;
  // std::vector<Double_t> LmomLreconfailed;
  // std::vector<Double_t> LmomLreconfailed_x;
  // std::vector<Double_t> LmomLreconfailed_y;
  // std::vector<Double_t> LmomLreconfailed_z;
  // std::vector<Double_t> pmomLreconfailed;
  // std::vector<Double_t> pmomLreconfailed_x;
  // std::vector<Double_t> pmomLreconfailed_y;
  // std::vector<Double_t> pmomLreconfailed_z;
  // std::vector<Double_t> pimomLreconfailed;
  // std::vector<Double_t> pimomLreconfailed_x;
  // std::vector<Double_t> pimomLreconfailed_y;
  // std::vector<Double_t> pimomLreconfailed_z;
  // std::vector<Double_t> ppidistLreconfailed;

  // Int_t ncombiPipair;
  // std::vector<Int_t> pipidPipair;
  // std::vector<Int_t> pimidPipair;
  // std::vector<Double_t> pipmomPipair;
  // std::vector<Double_t> pipmomPipair_x;
  // std::vector<Double_t> pipmomPipair_y;
  // std::vector<Double_t> pipmomPipair_z;
  // std::vector<Double_t> pimmomPipair;
  // std::vector<Double_t> pimmomPipair_x;
  // std::vector<Double_t> pimmomPipair_y;
  // std::vector<Double_t> pimmomPipair_z;
  // std::vector<Double_t> momPipair;
  // std::vector<Double_t> momPipair_x;
  // std::vector<Double_t> momPipair_y;
  // std::vector<Double_t> momPipair_z;
  // std::vector<Double_t> reconLmassPipair;
  // std::vector<Double_t> reconmassPipair;
  // std::vector<Double_t> pipidistPipair;

  std::vector<Int_t> GFfitstatus;
  std::vector<Int_t> GFpdgcode;
  std::vector<Int_t> GFnhtrack;
  std::vector<Double_t> GFcharge;
  std::vector<Double_t> GFchisqr;
  std::vector<Double_t> GFtof;
  //std::vector<Double_t> GFtracklen;
  std::vector<Double_t> GFpval;
  std::vector<std::vector<Double_t>> GFlayer;
  std::vector<std::vector<Double_t>> GFpos_x;
  std::vector<std::vector<Double_t>> GFpos_y;
  std::vector<std::vector<Double_t>> GFpos_z;
  std::vector<std::vector<Double_t>> GFmom;
  std::vector<std::vector<Double_t>> GFmom_x;
  std::vector<std::vector<Double_t>> GFmom_y;
  std::vector<std::vector<Double_t>> GFmom_z;
  std::vector<std::vector<Double_t>> GFresidual_x;
  std::vector<std::vector<Double_t>> GFresidual_y;
  std::vector<std::vector<Double_t>> GFresidual_z;
  std::vector<std::vector<Double_t>> GFresidual_p;
  std::vector<std::vector<Double_t>> GFresidual_px;
  std::vector<std::vector<Double_t>> GFresidual_py;
  std::vector<std::vector<Double_t>> GFresidual_pz;

  std::vector<Int_t> GFinside;

  Int_t GFntTpc_inside;
  Double_t GFprodvtx_x;
  Double_t GFprodvtx_y;
  Double_t GFprodvtx_z;

  std::vector<Double_t> GFtracklen;
  std::vector<Double_t> GFtrack2vtxdist;
  std::vector<Double_t> GFcalctof;
  std::vector<Double_t> GFsegHtof;
  std::vector<Double_t> GFtofHtof;
  std::vector<Double_t> GFtdiffHtof;
  std::vector<Double_t> GFposHtof;
  std::vector<Double_t> GFposx;
  std::vector<Double_t> GFposy;
  std::vector<Double_t> GFposz;
  std::vector<Double_t> GFinvbeta;
  std::vector<Double_t> GFm2;
  std::vector<Double_t> nsigma_tritonHtof;
  std::vector<Double_t> nsigma_deutronHtof;
  std::vector<Double_t> nsigma_protonHtof;
  std::vector<Double_t> nsigma_kaonHtof;
  std::vector<Double_t> nsigma_pionHtof;
  std::vector<Double_t> nsigma_electronHtof;

  Int_t kmflag;  
  Int_t    GFkmid;  
  Double_t GFkmmom;
  Double_t GFkmmom_x;
  Double_t GFkmmom_y;
  Double_t GFkmmom_z;
  Double_t GFkmphi;    
  Double_t GFkmtheta;
  Double_t GFkmtarget_dist;
  Double_t GFkmtargetvtx_x;
  Double_t GFkmtargetvtx_y;
  Double_t GFkmtargetvtx_z;
  Double_t GFkmtargetcenter_x;
  Double_t GFkmtargetcenter_y;
  Double_t GFkmtargetcenter_z;
  Double_t GFkmtargetcenter_dist;
  Double_t GFkmprodvtx_x;
  Double_t GFkmprodvtx_y;
  Double_t GFkmprodvtx_z;
  Double_t GFkmprodvtx_dist;
  Int_t    GFkmhtofid;
  Int_t    GFkmhtofseg;    
  Double_t GFkmtracklen;
  Double_t GFkmtof;
  Double_t GFkmmass2;
  Double_t GFkminvbeta;
  std::vector<Double_t> GFkmposHtof;  

  //Multi-track production vertex
  Double_t GFprodvtx_x_ll;
  Double_t GFprodvtx_y_ll;
  Double_t GFprodvtx_z_ll;
  Double_t GFprodvtx_x_l1;
  Double_t GFprodvtx_y_l1;
  Double_t GFprodvtx_z_l1;
  Double_t GFprodvtx_x_l2;
  Double_t GFprodvtx_y_l2;
  Double_t GFprodvtx_z_l2;
  Double_t GFprodvtx_x_l;
  Double_t GFprodvtx_y_l;
  Double_t GFprodvtx_z_l;

  Bool_t emptyflag;
  Bool_t pimflag;
  Bool_t lpiflag;
  Bool_t lflag;

  std::vector<Int_t> GFdecays_pdgcode;
  std::vector<Int_t> GFdecays_nhtrack;
  std::vector<Double_t> GFdecays_charge;
  std::vector<Double_t> GFdecays_chisqr;
  std::vector<Double_t> GFdecays_pval;
  std::vector<Int_t> GFdecays_htofid;
  std::vector<Double_t> GFdecays_tracklen;
  std::vector<Double_t> GFdecays_tof;
  std::vector<Double_t> GFdecays_mass2;
  std::vector<Double_t> GFdecays_invbeta;
  std::vector<Double_t> GFdecays_mom;
  std::vector<Double_t> GFdecays_mom_x;
  std::vector<Double_t> GFdecays_mom_y;
  std::vector<Double_t> GFdecays_mom_z;
  std::vector<Double_t> GFdecays_CMmom;
  std::vector<Double_t> GFdecays_CMmom_x;
  std::vector<Double_t> GFdecays_CMmom_y;
  std::vector<Double_t> GFdecays_CMmom_z;
  std::vector<Double_t> GFdecays_momloss;
  std::vector<Double_t> GFdecays_eloss;

  std::vector<Int_t> decays_id;
  std::vector<Double_t> decays_mom;
  std::vector<Double_t> decays_mom_x;
  std::vector<Double_t> decays_mom_y;
  std::vector<Double_t> decays_mom_z;
  std::vector<Double_t> decays_CMmom;
  std::vector<Double_t> decays_CMmom_x;
  std::vector<Double_t> decays_CMmom_y;
  std::vector<Double_t> decays_CMmom_z;

  Bool_t pipiflag;

  Int_t accident_multi;
  std::vector<Int_t> accident_id;

  //G4 branches
  Int_t NumberOfTracks;
  vector<int> PIDOfTrack;
  vector<int> ParentIDOfTrack;
  vector<double> VertexOfTrack_x;
  vector<double> VertexOfTrack_y;
  vector<double> VertexOfTrack_z;
  vector<double> MomentumOfTrack;
  vector<double> MomentumOfTrack_x;
  vector<double> MomentumOfTrack_y;
  vector<double> MomentumOfTrack_z;

  int G4kmid;
  int G4kmtid;
  double G4kmvtx_x;// Production vertex, identical to xi vert.
  double G4kmvtx_y;
  double G4kmvtx_z;
  double G4kmmom;
  double G4kmmom_x;
  double G4kmmom_y;
  double G4kmmom_z;

  int G4kpid; // scat Kp for KK reaction or scat P for KP reaction
  int G4kptid;
  double G4kpvtx_x;// Production vertex, identical to xi vert.
  double G4kpvtx_y;
  double G4kpvtx_z;
  double G4kpmom;
  double G4kpmom_x;
  double G4kpmom_y;
  double G4kpmom_z;

  int G4scatkmid;
  int G4scatkmtid;
  double G4scatkmvtx_x;// Production vertex, identical to xi vert.
  double G4scatkmvtx_y;
  double G4scatkmvtx_z;
  double G4scatkmmom;
  double G4scatkmmom_x;
  double G4scatkmmom_y;
  double G4scatkmmom_z;
  
#if XiRecon
  int G4xiid;
  double G4xivtx_x;//Production vertex
  double G4xivtx_y;
  double G4xivtx_z;
  double G4ximom;
  double G4ximom_x;//Momentum at production vtx
  double G4ximom_y;
  double G4ximom_z;

  double xivtx_x;
  double xivtx_y;
  double xivtx_z;

  int G4lid;
  double G4lvtx_x;//Production vertex, identical to xi decay vert
  double G4lvtx_y;
  double G4lvtx_z;
  double G4lmom;
  double G4lmom_x;
  double G4lmom_y;
  double G4lmom_z;
  double lvtx_x;
  double lvtx_y;
  double lvtx_z;

  int G4protonid;//Id of Geant4 Track
  int G4ptid;//G4 Id of TPC Track
  int G4pnh;//Number of G4 Hits
  int G4ptnh;//Number of correct G4 Hits
  double G4pvtx_x;//Production vertex, identical to l decay vert
  double G4pvtx_y;
  double G4pvtx_z;
  double G4pmom;
  double G4pmom_x;
  double G4pmom_y;
  double G4pmom_z;

  int G4extraprotonid;//Id of Geant4 Track
  int G4extraptid;//G4 Id of TPC Track
  int G4extrapnh;//Number of G4 Hits
  int G4extraptnh;//Number of correct G4 Hits
  double G4extrapvtx_x;//Production vertex, identical to l decay vert
  double G4extrapvtx_y;
  double G4extrapvtx_z;
  double G4extrapmom;
  double G4extrapmom_x;
  double G4extrapmom_y;
  double G4extrapmom_z;

  int ptid;
  int pnh;
  double pvtx_x;
  double pvtx_y;
  double pvtx_z;
  double pmom;
  double pmom_x;
  double pmom_y;
  double pmom_z;
  double GFpmom;
  double GFpmom_x;
  double GFpmom_y;
  double GFpmom_z;

  int G4pi1id;
  int G4pi1tid;
  int G4pi1nh;
  int G4pi1tnh;
  double G4pi1vtx_x;
  double G4pi1vtx_y;
  double G4pi1vtx_z;
  double G4pi1mom;
  double G4pi1mom_x;
  double G4pi1mom_y;
  double G4pi1mom_z;

  int pi1tid;
  int pi1nh;
  double pi1vtx_x;
  double pi1vtx_y;
  double pi1vtx_z;
  double pi1mom;
  double pi1mom_x;
  double pi1mom_y;
  double pi1mom_z;
  double GFpi1mom;
  double GFpi1mom_x;
  double GFpi1mom_y;
  double GFpi1mom_z;

  int G4pi2id;
  int G4pi2tid;
  int G4pi2nh;
  int G4pi2tnh;
  double G4pi2vtx_x;
  double G4pi2vtx_y;
  double G4pi2vtx_z;
  double G4pi2mom;
  double G4pi2mom_x;
  double G4pi2mom_y;
  double G4pi2mom_z;

  int pi2tid;
  int pi2nh;
  double pi2vtx_x;
  double pi2vtx_y;
  double pi2vtx_z;
  double pi2mom;
  double pi2mom_x;
  double pi2mom_y;
  double pi2mom_z;
  double GFpi2mom;
  double GFpi2mom_x;
  double GFpi2mom_y;
  double GFpi2mom_z;

  bool lgood,xigood;
  bool p_tracked,pi1_tracked,pi2_tracked;
  bool extrap_tracked;
  double pt_mom0,pi1t_mom0,pi2t_mom0;

  bool lforced,xiforced;
  bool p_k18cut,p_kuramacut,p_isbeamcut,p_accidentalcut,p_pidcut,p_chargecut,p_directioncut,p_distcut;
  bool pi1_k18cut,pi1_kuramacut,pi1_isbeamcut,pi1_accidentalcut,pi1_pidcut,pi1_chargecut,pi1_directioncut,pi1_distcut;
  bool pi2_k18cut,pi2_kuramacut,pi2_isbeamcut,pi2_accidentalcut,pi2_pidcut,pi2_chargecut,pi2_directioncut,pi2_distcut;
  bool lpi_dist_nan,ppi_dist_nan;
  bool l_vertcut,xi_vertcut,xi_targetdistcut;
  double lpi_dist_forced,ppi_dist_forced;
  double xi_mass_forced,l_mass_forced;
#endif
#if LLRecon
  int G4l1id;
  double G4l1vtx_x;
  double G4l1vtx_y;
  double G4l1vtx_z;
  double G4l1mom;
  double G4l1mom_x;
  double G4l1mom_y;
  double G4l1mom_z;

  double l1vtx_x;
  double l1vtx_y;
  double l1vtx_z;

  int G4l2id;
  int G4l2vtx_x;
  int G4l2vtx_y;
  int G4l2vtx_z;
  double G4l2mom;
  double G4l2mom_x;
  double G4l2mom_y;
  double G4l2mom_z;

  double l2vtx_x;
  double l2vtx_y;
  double l2vtx_z;

  double G4llmass;

  int G4p1id;
  int G4p1tid;
  int G4p1nh;
  int G4p1tnh;
  double G4p1vtx_x;
  double G4p1vtx_y;
  double G4p1vtx_z;
  double G4p1mom;
  double G4p1mom_x;
  double G4p1mom_y;
  double G4p1mom_z;

  int p1tid;
  int p1nh;
  double p1vtx_x;
  double p1vtx_y;
  double p1vtx_z;
  double p1mom;
  double p1mom_x;
  double p1mom_y;
  double p1mom_z;
  double GFp1mom;
  double GFp1mom_x;
  double GFp1mom_y;
  double GFp1mom_z;

  int G4p2id;
  int G4p2tid;
  int G4p2nh;
  int G4p2tnh;
  double G4p2vtx_x;
  double G4p2vtx_y;
  double G4p2vtx_z;
  double G4p2mom;
  double G4p2mom_x;
  double G4p2mom_y;
  double G4p2mom_z;

  int p2tid;
  int p2nh;
  double p2vtx_x;
  double p2vtx_y;
  double p2vtx_z;
  double p2mom;
  double p2mom_x;
  double p2mom_y;
  double p2mom_z;
  double GFp2mom;
  double GFp2mom_x;
  double GFp2mom_y;
  double GFp2mom_z;

  int G4pi1id;
  int G4pi1tid;
  int G4pi1nh;
  int G4pi1tnh;
  double G4pi1vtx_x;
  double G4pi1vtx_y;
  double G4pi1vtx_z;
  double G4pi1mom;
  double G4pi1mom_x;
  double G4pi1mom_y;
  double G4pi1mom_z;

  int pi1tid;
  int pi1nh;
  double pi1vtx_x;
  double pi1vtx_y;
  double pi1vtx_z;
  double pi1mom;
  double pi1mom_x;
  double pi1mom_y;
  double pi1mom_z;
  double GFpi1mom;
  double GFpi1mom_x;
  double GFpi1mom_y;
  double GFpi1mom_z;

  int G4pi2id;
  int G4pi2tid;
  int G4pi2nh;
  int G4pi2tnh;
  double G4pi2vtx_x;
  double G4pi2vtx_y;
  double G4pi2vtx_z;
  double G4pi2mom;
  double G4pi2mom_x;
  double G4pi2mom_y;
  double G4pi2mom_z;

  int pi2tid;
  int pi2nh;
  double pi2vtx_x;
  double pi2vtx_y;
  double pi2vtx_z;
  double pi2mom;
  double pi2mom_x;
  double pi2mom_y;
  double pi2mom_z;
  double GFpi2mom;
  double GFpi2mom_x;
  double GFpi2mom_y;
  double GFpi2mom_z;

  bool l1good,l2good,llswap;
  bool p1_tracked,p2_tracked,pi1_tracked,pi2_tracked;
  double p1t_mom0,p2t_mom0,pi1t_mom0,pi2t_mom0;
#endif

  void clear( void )
  {
    evnum = 0;
    status = 0;
    trigpat.clear();
    trigflag.clear();

    nhHtof = 0;
    HtofSeg.clear();
    tHtof.clear();
    dtHtof.clear();
    deHtof.clear();
    posHtof.clear();
    G4tidHtof.clear();

    nhFtof = 0;
    FtofSeg.clear();
    tFtof.clear();
    deFtof.clear();
    posFtof.clear();

    ntK18 = 0;
    pK18.clear();
    chisqrK18.clear();
    xtgtK18.clear();
    ytgtK18.clear();
    utgtK18.clear();
    vtgtK18.clear();

    ntKurama = 0;
    chisqrKurama.clear();
    pKurama.clear();
    qKurama.clear();
    m2.clear();
    thetaKurama.clear();
    xtgtKurama.clear();
    ytgtKurama.clear();
    utgtKurama.clear();
    vtgtKurama.clear();
    pathwcKurama.clear();
    xin.clear();
    yin.clear();
    zin.clear();
    pxin.clear();
    pyin.clear();
    pzin.clear();
    xout.clear();
    yout.clear();
    zout.clear();
    pxout.clear();
    pyout.clear();
    pzout.clear();

    nKm = 0;
    nKp = 0;
    nKK = 0;
    inside.clear();
    vtx.clear();
    vty.clear();
    vtz.clear();
    closeDist.clear();
    MissMass.clear();
    MissMassCorr.clear();
    MissMassCorrDE.clear();
    MissMassNucl.clear();
    MissMassNuclCorr.clear();
    MissMassNuclCorrDE.clear();    
    pOrg.clear();
    pCalc.clear();
    pCorr.clear();
    pCorrDE.clear();
    xb.clear();
    yb.clear();
    ub.clear();
    vb.clear();
    xs.clear();
    ys.clear();
    us.clear();
    vs.clear();
    Kflag.clear();

    isgoodTPCK18.clear();
    chisqrTPCK18.clear();
    qTPCK18.clear();
    pTPCK18.clear();
    xtgtTPCK18.clear();
    ytgtTPCK18.clear();
    utgtTPCK18.clear();
    vtgtTPCK18.clear();
    thetaTPCK18.clear();
    lhtofTPCK18.clear();
    xhtofTPCK18.clear();
    yhtofTPCK18.clear();
    lvpTPCK18.clear();
    xvpTPCK18.clear();
    yvpTPCK18.clear();

    isgoodTPCKurama.clear();
    kflagTPCKurama.clear();
    chisqrTPCKurama.clear();
    pTPCKurama.clear();
    qTPCKurama.clear();
    m2TPCKurama.clear();
    xtgtTPCKurama.clear();
    ytgtTPCKurama.clear();
    utgtTPCKurama.clear();
    vtgtTPCKurama.clear();
    thetaTPCKurama.clear();
    pathTPCKurama.clear();
    lhtofTPCKurama.clear();
    xhtofTPCKurama.clear();
    yhtofTPCKurama.clear();
    lvpTPCKurama.clear();
    xvpTPCKurama.clear();
    yvpTPCKurama.clear();

    isgoodTPC.clear();
    insideTPC.clear();
    vtxTPC.clear();
    vtyTPC.clear();
    vtzTPC.clear();
    closeDistTPC.clear();
    MissMassTPC.clear();
    MissMassCorrTPC.clear();
    MissMassCorrDETPC.clear();
    pOrgTPC.clear();
    pCalcTPC.clear();
    pCorrTPC.clear();
    pCorrDETPC.clear();
    thetaTPC.clear();
    thetaCMTPC.clear();
    costCMTPC.clear();
    xbTPC.clear();
    ybTPC.clear();
    ubTPC.clear();
    vbTPC.clear();
    xsTPC.clear();
    ysTPC.clear();
    usTPC.clear();
    vsTPC.clear();

    BE.clear();
    BETPC.clear();
    BE_LL.clear();
    BETPC_LL.clear();
    km_mom_x.clear();
    km_mom_y.clear();
    km_mom_z.clear();
    kp_mom_x.clear();
    kp_mom_y.clear();
    kp_mom_z.clear();

    nclTpc = 0;
    cluster_x.clear();
    cluster_y.clear();
    cluster_z.clear();
    cluster_de.clear();
    cluster_size.clear();
    cluster_layer.clear();
    cluster_mrow.clear();
    cluster_de_center.clear();
    cluster_x_center.clear();
    cluster_y_center.clear();
    cluster_z_center.clear();
    cluster_row_center.clear();
    cluster_houghflag.clear();
    cluster_G4tid.clear();
    cluster_G4protonid.clear();

    remain_nclTpc = 0;
    remain_cluster_x.clear();
    remain_cluster_y.clear();
    remain_cluster_z.clear();
    remain_cluster_de.clear();
    remain_cluster_size.clear();
    remain_cluster_layer.clear();
    remain_cluster_mrow.clear();
    remain_cluster_de_center.clear();
    remain_cluster_x_center.clear();
    remain_cluster_y_center.clear();
    remain_cluster_z_center.clear();
    remain_cluster_row_center.clear();
    remain_cluster_houghflag.clear();
    remain_cluster_G4tid.clear();
    remain_cluster_G4protonid.clear();

    ntTpc = 0;
    nhtrack.clear();
    trackid.clear();
    isXi.clear();
    isBeam.clear();
    isKurama.clear();
    isK18.clear();
    isAccidental.clear();
    isMultiloop.clear();
    isInTarget.clear();
    charge.clear();
    pid.clear();
    chisqr.clear();
    pval.clear();
    purity.clear();
    efficiency.clear();
    G4tid.clear();
    G4pid.clear();
    helix_cx.clear();
    helix_cy.clear();
    helix_z0.clear();
    helix_r.clear();
    helix_dz.clear();
    dE.clear();
    dEdx.clear();
    mom0.clear();
    path.clear();
    isElectron.clear();
    nsigma_triton.clear();
    nsigma_deutron.clear();
    nsigma_proton.clear();
    nsigma_kaon.clear();
    nsigma_pion.clear();
    nsigma_electron.clear();

    hitlayer.clear();
    hitpos_x.clear();
    hitpos_y.clear();
    hitpos_z.clear();
    calpos_x.clear();
    calpos_y.clear();
    calpos_z.clear();
    mom_x.clear();
    mom_y.clear();
    mom_z.clear();
    residual.clear();
    residual_x.clear();
    residual_y.clear();
    residual_z.clear();
    resolution_x.clear();
    resolution_y.clear();
    resolution_z.clear();
    helix_t.clear();
    pathhit.clear();
    alpha.clear();
    track_cluster_de.clear();
    track_cluster_size.clear();
    track_cluster_mrow.clear();
    track_cluster_de_center.clear();
    track_cluster_x_center.clear();
    track_cluster_y_center.clear();
    track_cluster_z_center.clear();
    track_cluster_row_center.clear();

    chargeIndistinguishable.clear();
    chisqr_inverted.clear();
    pval_inverted.clear();
    helix_cx_inverted.clear();
    helix_cy_inverted.clear();
    helix_z0_inverted.clear();
    helix_r_inverted.clear();
    helix_dz_inverted.clear();
    mom0_inverted.clear();
    pid_inverted.clear();

    nvtxTpc = 0;
    vtx_x.clear();
    vtx_y.clear();
    vtx_z.clear();
    vtx_dist.clear();
    vtx_angle.clear();
    vtxid.clear();
    vtxmom_theta.clear();
    vtxpos_x.clear();
    vtxpos_y.clear();
    vtxpos_z.clear();
    vtxmom_x.clear();
    vtxmom_y.clear();
    vtxmom_z.clear();

    nvtxTpcClustered = 0;
    clusteredVtx_x.clear();
    clusteredVtx_y.clear();
    clusteredVtx_z.clear();
    clusteredVtxid.clear();

    // ncombiLreconfailed = 0;
    // pidLreconfailed.clear();
    // piidLreconfailed.clear();
    // LdecayvtxLreconfailed_x.clear();
    // LdecayvtxLreconfailed_y.clear();
    // LdecayvtxLreconfailed_z.clear();
    // LmassLreconfailed.clear();
    // LmomLreconfailed.clear();
    // LmomLreconfailed_x.clear();
    // LmomLreconfailed_y.clear();
    // LmomLreconfailed_z.clear();
    // pmomLreconfailed.clear();
    // pmomLreconfailed_x.clear();
    // pmomLreconfailed_y.clear();
    // pmomLreconfailed_z.clear();
    // pimomLreconfailed.clear();
    // pimomLreconfailed_x.clear();
    // pimomLreconfailed_y.clear();
    // pimomLreconfailed_z.clear();
    // ppidistLreconfailed.clear();

    // ncombiPipair = 0;
    // pipidPipair.clear();
    // pimidPipair.clear();
    // pipmomPipair.clear();
    // pipmomPipair_x.clear();
    // pipmomPipair_y.clear();
    // pipmomPipair_z.clear();
    // pimmomPipair.clear();
    // pimmomPipair_x.clear();
    // pimmomPipair_y.clear();
    // pimmomPipair_z.clear();
    // momPipair.clear();
    // momPipair_x.clear();
    // momPipair_y.clear();
    // momPipair_z.clear();
    // reconLmassPipair.clear();
    // reconmassPipair.clear();
    // pipidistPipair.clear();

    // isLambda.clear();
    // ncombiLambda.clear();
    // distLambda.clear();
    // angleLambda.clear();
    // bestmassLambda.clear();
    // massLambda.clear();
    // vtxLambda_x.clear();
    // vtxLambda_y.clear();
    // vtxLambda_z.clear();
    // momLambda.clear();
    // momLambda_x.clear();
    // momLambda_y.clear();
    // momLambda_z.clear();
    // decaysidLambda.clear();
    // decaysmomLambda.clear();
    // decaysmomLambda_x.clear();
    // decaysmomLambda_y.clear();
    // decaysmomLambda_z.clear();

    GFcharge.clear();
    GFchisqr.clear();
    GFtof.clear();
    GFtracklen.clear();
    GFpval.clear();
    GFfitstatus.clear();
    GFpdgcode.clear();
    GFnhtrack.clear();
    GFlayer.clear();
    GFpos_x.clear();
    GFpos_y.clear();
    GFpos_z.clear();
    GFmom.clear();
    GFmom_x.clear();
    GFmom_y.clear();
    GFmom_z.clear();
    GFresidual_x.clear();
    GFresidual_y.clear();
    GFresidual_z.clear();
    GFresidual_p.clear();
    GFresidual_px.clear();
    GFresidual_py.clear();
    GFresidual_pz.clear();


    GFinside.clear();

    GFntTpc_inside = 0;
    GFprodvtx_x = qnan;
    GFprodvtx_y = qnan;
    GFprodvtx_z = qnan;

    GFtracklen.clear();
    GFtrack2vtxdist.clear();
    GFcalctof.clear();
    GFsegHtof.clear();
    GFtofHtof.clear();
    GFtdiffHtof.clear();
    GFposHtof.clear();
    GFposx.clear();
    GFposy.clear();
    GFposz.clear();
    GFinvbeta.clear();
    GFm2.clear();
    nsigma_tritonHtof.clear();
    nsigma_deutronHtof.clear();
    nsigma_protonHtof.clear();
    nsigma_kaonHtof.clear();
    nsigma_pionHtof.clear();
    nsigma_electronHtof.clear();

    kmflag = false;
    GFkmid = -1;    
    GFkmmom = qnan;
    GFkmmom_x = qnan;
    GFkmmom_y = qnan;
    GFkmmom_z = qnan;
    GFkmtheta = qnan;
    GFkmphi   = qnan;        
    GFkmtarget_dist = qnan;
    GFkmtargetvtx_x = qnan;
    GFkmtargetvtx_y = qnan;
    GFkmtargetvtx_z = qnan;
    GFkmtargetcenter_dist = qnan;
    GFkmtargetcenter_x = qnan;
    GFkmtargetcenter_y = qnan;
    GFkmtargetcenter_z = qnan;
    GFkmprodvtx_x = qnan;
    GFkmprodvtx_y = qnan;
    GFkmprodvtx_z = qnan;
    GFkmprodvtx_dist = qnan;
    GFkmhtofid = -1;
    GFkmhtofseg = -1;            
    GFkmtracklen = qnan;    
    GFkmtof = qnan;
    GFkmmass2 = qnan;
    GFkminvbeta = qnan;
    GFkmposHtof.clear();    

    GFprodvtx_x_ll = qnan;
    GFprodvtx_y_ll = qnan;
    GFprodvtx_z_ll = qnan;
    GFprodvtx_x_l1 = qnan;
    GFprodvtx_y_l1 = qnan;
    GFprodvtx_z_l1 = qnan;
    GFprodvtx_x_l2 = qnan;
    GFprodvtx_y_l2 = qnan;
    GFprodvtx_z_l2 = qnan;
    GFprodvtx_x_l = qnan;
    GFprodvtx_y_l = qnan;
    GFprodvtx_z_l = qnan;

    GFdecays_pdgcode.clear();
    GFdecays_nhtrack.clear();
    GFdecays_charge.clear();
    GFdecays_chisqr.clear();
    GFdecays_pval.clear();
    GFdecays_htofid.clear();
    GFdecays_tracklen.clear();
    GFdecays_tof.clear();
    GFdecays_mass2.clear();
    GFdecays_invbeta.clear();
    GFdecays_mom.clear();
    GFdecays_mom_x.clear();
    GFdecays_mom_y.clear();
    GFdecays_mom_z.clear();
    GFdecays_CMmom.clear();
    GFdecays_CMmom_x.clear();
    GFdecays_CMmom_y.clear();
    GFdecays_CMmom_z.clear();
    GFdecays_momloss.clear();
    GFdecays_eloss.clear();

    decays_id.clear();
    decays_mom.clear();
    decays_mom_x.clear();
    decays_mom_y.clear();
    decays_mom_z.clear();
    decays_CMmom.clear();
    decays_CMmom_x.clear();
    decays_CMmom_y.clear();
    decays_CMmom_z.clear();
    
    emptyflag = false;
    pimflag = false;
    lflag = false;

    //G4 Initialization
    NumberOfTracks = -1;
    PIDOfTrack.clear();
    ParentIDOfTrack.clear();
    VertexOfTrack_x.clear();
    VertexOfTrack_y.clear();
    VertexOfTrack_z.clear();
    MomentumOfTrack.clear();
    MomentumOfTrack_x.clear();
    MomentumOfTrack_y.clear();
    MomentumOfTrack_z.clear();
    G4kmid = qnan,G4kmtid = qnan;
    G4kmvtx_x = qnan,G4kmvtx_y = qnan,G4kmvtx_z = qnan;
    G4kmmom = qnan,G4kmmom_x = qnan,G4kmmom_y = qnan,G4kmmom_z = qnan;

    G4kpid = qnan,G4kptid = qnan;
    G4kpvtx_x = qnan,G4kpvtx_y = qnan,G4kpvtx_z = qnan;
    G4kpmom = qnan,G4kpmom_x = qnan,G4kpmom_y = qnan,G4kpmom_z = qnan;

    G4scatkmid = qnan,G4scatkmtid = qnan;
    G4scatkmvtx_x = qnan,G4scatkmvtx_y = qnan,G4scatkmvtx_z = qnan;
    G4scatkmmom = qnan,G4scatkmmom_x = qnan,G4scatkmmom_y = qnan,G4scatkmmom_z = qnan;
    
#if XiRecon
    G4xiid = qnan;
    G4xivtx_x = qnan,G4xivtx_y = qnan,G4xivtx_z = qnan;
    G4ximom = qnan,G4ximom_x = qnan,G4ximom_y = qnan,G4ximom_z = qnan;
    xivtx_x = qnan,xivtx_y = qnan,xivtx_z = qnan;

    G4lid = qnan;
    G4lvtx_x = qnan,G4lvtx_y = qnan,G4lvtx_z = qnan;
    G4lmom = qnan,G4lmom_x = qnan,G4lmom_y = qnan,G4lmom_z = qnan;
    lvtx_x = qnan,lvtx_y = qnan,lvtx_z = qnan;

    G4protonid = qnan, G4ptid = qnan, G4pnh = 0, G4ptnh = 0;
    G4pvtx_x = qnan, G4pvtx_y = qnan, G4pvtx_z = qnan;
    G4pmom = qnan, G4pmom_x = qnan, G4pmom_y = qnan, G4pmom_z = qnan;

    G4extraprotonid = -1, G4extraptid = -1, G4extrapnh = 0, G4extraptnh = 0;
    G4extrapvtx_x = qnan, G4extrapvtx_y = qnan, G4extrapvtx_z = qnan;
    G4extrapmom = qnan, G4extrapmom_x = qnan, G4extrapmom_y = qnan, G4extrapmom_z = qnan;

    ptid = qnan,pnh = 0;
    pvtx_x = qnan,pvtx_y = qnan,pvtx_z = qnan;
    pmom = qnan,pmom_x = qnan,pmom_y = qnan,pmom_z = qnan;

    G4pi1id = qnan,G4pi1tid = qnan,G4pi1nh = 0,G4pi1tnh = 0;
    G4pi1vtx_x = qnan,G4pi1vtx_y = qnan,G4pi1vtx_z = qnan;
    G4pi1mom = qnan,G4pi1mom_x = qnan,G4pi1mom_y = qnan,G4pi1mom_z = qnan;
    lgood = false;
    xigood = false;
    p_tracked = false;
    pi1_tracked = false;
    pi2_tracked = false;
    pt_mom0 = qnan;
    pi1t_mom0 = qnan;
    pi2t_mom0 = qnan;

    lforced = false;
    xiforced = false;

    p_accidentalcut = false,p_isbeamcut =false,p_k18cut = false,p_kuramacut = false,p_pidcut = false,p_chargecut = false,p_directioncut = false,p_distcut = false;
    pi1_accidentalcut = false,pi1_isbeamcut = false,pi1_k18cut = false,pi1_kuramacut = false,pi1_pidcut = false,pi1_chargecut = false,pi1_directioncut = false,pi1_distcut = false;
    pi2_accidentalcut = false,pi2_isbeamcut = false,pi2_k18cut = false,pi2_kuramacut = false,pi2_pidcut = false,pi2_chargecut = false,pi2_directioncut = false,pi2_distcut = false;
    lpi_dist_nan = false,ppi_dist_nan = false;
    l_vertcut = false,xi_vertcut = false,xi_targetdistcut;
    lpi_dist_forced = qnan;
    ppi_dist_forced = qnan;
    xi_mass_forced = qnan;
    l_mass_forced = qnan;
#endif
#if LLRecon
    G4l1id =-1;
    G4l1vtx_x = qnan, G4l1vtx_y = qnan, G4l1vtx_z = qnan;
    G4l1mom = qnan, G4l1mom_x = qnan, G4l1mom_y = qnan, G4l1mom_z = qnan;
    l1vtx_x = qnan, l1vtx_y = qnan, l1vtx_z = qnan;

    G4l2id =-1;
    G4l2vtx_x = qnan, G4l2vtx_y = qnan, G4l2vtx_z = qnan;
    G4l2mom = qnan, G4l2mom_x = qnan, G4l2mom_y = qnan, G4l2mom_z = qnan;
    l2vtx_x = qnan, l2vtx_y = qnan, l2vtx_z = qnan;

    G4llmass = qnan;

    G4p1id = -1,G4p1tid = -1,G4p1nh = -1, G4p1tnh = -1;
    G4p1vtx_x = qnan, G4p1vtx_y = qnan, G4p1vtx_z = qnan;
    G4p1mom = qnan, G4p1mom_x = qnan, G4p1mom_y = qnan, G4p1mom_z = qnan;

    G4p2id = -1, G4p2tid = -1, G4p2nh = -1, G4p2tnh = -1;
    G4p2vtx_x = qnan, G4p2vtx_y = qnan, G4p2vtx_z = qnan;
    G4p2mom = qnan, G4p2mom_x = qnan, G4p2mom_y = qnan, G4p2mom_z = qnan;

    G4pi1id = -1, G4pi1tid = -1,G4pi1nh = -1, G4pi1tnh = -1;
    G4pi1vtx_x = qnan, G4pi1vtx_y = qnan, G4pi1vtx_z = qnan;
    G4pi1mom = qnan, G4pi1mom_x = qnan, G4pi1mom_y = qnan, G4pi1mom_z = qnan;

    G4pi2id = -1, G4pi2tid = -1,G4pi2nh = -1, G4pi2tnh = -1;
    G4pi2vtx_x = qnan, G4pi2vtx_y = qnan, G4pi2vtx_z = qnan;
    G4pi2mom = qnan, G4pi2mom_x = qnan, G4pi2mom_y = qnan, G4pi2mom_z = qnan;

    l1good = false; l2good = false; llswap = false;
    p1_tracked = false; p2_tracked = false;
    pi1_tracked = false; pi2_tracked = false;
    p1t_mom0 = qnan; p2t_mom0 = qnan;
    pi1t_mom0 = qnan; pi2t_mom0 = qnan;

#endif
  }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* evnum;

  TTreeReaderValue<Int_t>* nhHtof;
  TTreeReaderValue<std::vector<Double_t>>* HtofSeg;
  TTreeReaderValue<std::vector<Double_t>>* tHtof;
  TTreeReaderValue<std::vector<Double_t>>* dtHtof;
  TTreeReaderValue<std::vector<Double_t>>* deHtof;
  TTreeReaderValue<std::vector<Double_t>>* posHtof;
  TTreeReaderValue<std::vector<Int_t>>* G4tidHtof;

  int nhittpc;
  Int_t ititpc[MaxTPCHits];
  Double_t xtpc[MaxTPCHits];//with resolution
  Double_t ytpc[MaxTPCHits];//with resolution
  Double_t ztpc[MaxTPCHits];//with resolution
  Double_t pxtpc[MaxTPCHits];//with resolution
  Double_t pytpc[MaxTPCHits];//with resolution
  Double_t pztpc[MaxTPCHits];//with resolution

  Int_t NumberOfTracks;
  Int_t PIDOfTrack[1000];
  Int_t ParentIDOfTrack[1000];
  Double_t VertexOfTrack_x[1000];
  Double_t VertexOfTrack_y[1000];
  Double_t VertexOfTrack_z[1000];
  Double_t MomentumOfTrack[1000];
  Double_t MomentumOfTrack_x[1000];
  Double_t MomentumOfTrack_y[1000];
  Double_t MomentumOfTrack_z[1000];

  TTreeReaderValue<Int_t>* nclTpc; // Number of clusters
  TTreeReaderValue<Int_t>* remain_nclTpc; // Number of clusters without tracks
  TTreeReaderValue<std::vector<Double_t>>* cluster_x;
  TTreeReaderValue<std::vector<Double_t>>* cluster_y;
  TTreeReaderValue<std::vector<Double_t>>* cluster_z;
  TTreeReaderValue<std::vector<Double_t>>* cluster_de;
  TTreeReaderValue<std::vector<Int_t>>* cluster_size;
  TTreeReaderValue<std::vector<Int_t>>* cluster_layer;
  TTreeReaderValue<std::vector<Double_t>>* cluster_mrow;
  TTreeReaderValue<std::vector<Double_t>>* cluster_de_center;
  TTreeReaderValue<std::vector<Double_t>>* cluster_x_center;
  TTreeReaderValue<std::vector<Double_t>>* cluster_y_center;
  TTreeReaderValue<std::vector<Double_t>>* cluster_z_center;
  TTreeReaderValue<std::vector<Int_t>>* cluster_row_center;
  TTreeReaderValue<std::vector<Int_t>>* cluster_houghflag;

  TTreeReaderValue<Int_t>* ntTpc; // Number of Tracks
  TTreeReaderValue<std::vector<Int_t>>* nhtrack; // Number of Hits (in 1 tracks)
  TTreeReaderValue<std::vector<Int_t>>* trackid; //for Kurama K1.8 tracks
  TTreeReaderValue<std::vector<Int_t>>* isXi;
  TTreeReaderValue<std::vector<Int_t>>* isBeam;
  TTreeReaderValue<std::vector<Int_t>>* isKurama;
  TTreeReaderValue<std::vector<Int_t>>* isK18;
  TTreeReaderValue<std::vector<Int_t>>* isAccidental;
  TTreeReaderValue<std::vector<Int_t>>* isMultiloop;
  TTreeReaderValue<std::vector<Int_t>>* charge;//Helix charge
  TTreeReaderValue<std::vector<Int_t>>* pid;
  TTreeReaderValue<std::vector<Double_t>>* purity;
  TTreeReaderValue<std::vector<Double_t>>* efficiency;
  TTreeReaderValue<std::vector<Int_t>>* G4tid;
  TTreeReaderValue<std::vector<Double_t>>* chisqr;
  TTreeReaderValue<std::vector<Double_t>>* pval;
  TTreeReaderValue<std::vector<Double_t>>* helix_cx;
  TTreeReaderValue<std::vector<Double_t>>* helix_cy;
  TTreeReaderValue<std::vector<Double_t>>* helix_z0;
  TTreeReaderValue<std::vector<Double_t>>* helix_r;
  TTreeReaderValue<std::vector<Double_t>>* helix_dz;
  TTreeReaderValue<std::vector<Double_t>>* dE;
  TTreeReaderValue<std::vector<Double_t>>* dEdx; //reference dedx
  TTreeReaderValue<std::vector<Double_t>>* mom0;//Helix momentum at Y = 0
  TTreeReaderValue<std::vector<Double_t>>* path;//Helix path
  TTreeReaderValue<std::vector<Int_t>>* isElectron;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_triton;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_deutron;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_proton;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_kaon;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_pion;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_electron;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitlayer;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitpos_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitpos_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitpos_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* calpos_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* calpos_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* calpos_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* mom_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* mom_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* mom_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* helix_t;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* pathhit;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* alpha;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_size;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_mrow;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de_center;

  TTreeReaderValue<std::vector<Int_t>>* chargeIndistinguishable;
  TTreeReaderValue<std::vector<Int_t>>* pid_inverted;
  TTreeReaderValue<std::vector<Double_t>>* chisqr_inverted;
  TTreeReaderValue<std::vector<Double_t>>* pval_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_cx_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_cy_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_z0_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_r_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_dz_inverted;
  TTreeReaderValue<std::vector<Double_t>>* mom0_inverted;//Helix momentum at Y = 0

  TTreeReaderValue<std::vector<Int_t>>* isLambda;
  TTreeReaderValue<std::vector<Int_t>>* ncombiLambda;
  TTreeReaderValue<std::vector<Double_t>>* distLambda;
  TTreeReaderValue<std::vector<Double_t>>* angleLambda;
  TTreeReaderValue<std::vector<Double_t>>* bestmassLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* massLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxLambda_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxLambda_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxLambda_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* momLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* momLambda_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* momLambda_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* momLambda_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* decaysidLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* decaysmomLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* decaysmomLambda_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* decaysmomLambda_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* decaysmomLambda_z;

  TTreeReaderValue<Int_t>* nvtxTpc;
  TTreeReaderValue<std::vector<Double_t>>* vtx_x;
  TTreeReaderValue<std::vector<Double_t>>* vtx_y;
  TTreeReaderValue<std::vector<Double_t>>* vtx_z;
  TTreeReaderValue<std::vector<Double_t>>* vtx_dist;
  TTreeReaderValue<std::vector<Double_t>>* vtx_angle;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxid;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxmom_theta;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxpos_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxpos_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxpos_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxmom_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxmom_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxmom_z;

  TTreeReaderValue<Int_t>* nvtxTpcClustered;
  TTreeReaderValue<std::vector<Double_t>>* clusteredVtx_x;
  TTreeReaderValue<std::vector<Double_t>>* clusteredVtx_y;
  TTreeReaderValue<std::vector<Double_t>>* clusteredVtx_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* clusteredVtxid;

};

namespace root
{
  Event  event;
  Src    src;
  TH1   *h[MaxHist];
  TTree *tree;
  enum eDetHid {
    TPCHid    = 100000,
    genfitHid = 200000
  };
}

//_____________________________________________________________________________
int
main( int argc, char **argv )
{
  std::vector<std::string> arg( argv, argv+argc );

  if( !CheckArg( arg ) )
    return EXIT_FAILURE;
  if( !DstOpen( arg ) )
    return EXIT_FAILURE;
  if( !gConf.Initialize( arg[kConfFile] ) )
    return EXIT_FAILURE;
  if( !gConf.InitializeUnpacker() )
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries( TTreeCont );
  if (max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();

  //Initiallize Geometry, Field, Fitter
  HypTPCFitter* fitter = new HypTPCFitter(tpcGeo.Data(), Const_field);
  //Initiallize the genfit track container
  HypTPCTask& GFTrackCont = HypTPCTask::GetInstance();
  GFTrackCont.SetVerbosity(verbosity);
  std::cout<<"GenFit verbosity = "<<"-1: Silent, 0: Minimum, 1: Errors only, 2: Errors and Warnings, 3: Verbose mode, long term debugging(default)"<<std::endl;
  std::cout<<"Current verbosity = "<<GFTrackCont.GetVerbosity()<<std::endl;

#if 0
  GFTrackCont.DebugMode();
#endif

  Int_t ievent = skip;
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
    gCounter.check();
    InitializeEvent();
    if( DstRead( ievent ) ) tree->Fill();
    GFTrackCont.Clear();
  }

  std::cout << "#D Event Number: " << std::setw(6)
            << ievent << std::endl;

  DstClose();

  delete fitter;
  return EXIT_SUCCESS;
}

//_____________________________________________________________________________
Bool_t
dst::InitializeEvent( void )
{
  event.clear();

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstOpen( std::vector<std::string> arg )
{
  Int_t open_file = 0;
  Int_t open_tree = 0;
  for( Int_t i=0; i<nArgc; ++i ){
    if( i==kProcess || i==kConfFile || i==kOutFile ) continue;
    open_file += OpenFile( TFileCont[i], arg[i] );
    open_tree += OpenTree( TFileCont[i], TTreeCont[i], TreeName[i] );
  }

  if( open_file!=open_tree || open_file!=nArgc-3 )
    return false;
  if( !CheckEntries( TTreeCont ) )
    return false;

  TFileCont[kOutFile] = new TFile( arg[kOutFile].c_str(), "recreate" );

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstRead( Int_t ievent )
{

  static const auto ElectronMass = pdg::ElectronMass();
  static const auto PionMass = pdg::PionMass();
  static const auto KaonMass = pdg::KaonMass();
  static const auto PhiMass = pdg::Mass(333);
  static const auto ProtonMass = pdg::ProtonMass();
  static const auto LambdaMass = pdg::LambdaMass();
  static const auto XiMinusMass = pdg::XiMinusMass();
  static const auto m12C = 11.174864;
  static const auto m11B = 10.252548;
  static const auto m10Be = 9.325504;
  static const auto me = 0.001*0.5109989461;
  static const int XiMinusPdgCode = 3312;
  Double_t pdgmass[3] = {ProtonMass, KaonMass, PionMass};
  TVector3 tgtpos(0, 0, tpc::ZTarget);
  TVector3 qnan_vec = TVector3(qnan, qnan, qnan);

  static const auto mindEdxSigKaon = gUser.GetParameter("MindEdXSigKaon");
  static const auto maxdEdxSigKaon = gUser.GetParameter("MaxdEdXSigKaon");  

  static const auto minM2Km = gUser.GetParameter("MinM2Kaon");
  static const auto maxM2Km = gUser.GetParameter("MaxM2Kaon");
  
  //if( ievent%100000==0 ){
  //if( ievent%1000==0 ){
  //if( ievent%100==0 ){
  //if( ievent%10==0 ){
  if( ievent%1==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  event.evnum = **src.evnum;
  event.NumberOfTracks = src.NumberOfTracks;

  std::vector<TVector3> G4Hits;
  std::vector<TVector3> G4Moms;
  for(int ih=0;ih<src.nhittpc;++ih){
    TVector3 G4Hit(src.xtpc[ih],src.ytpc[ih],src.ztpc[ih]);
    TVector3 G4Mom(src.pxtpc[ih],src.pytpc[ih],src.pztpc[ih]);
    G4Hits.push_back(G4Hit);
    G4Moms.push_back(G4Mom);
  }

  std::vector<Int_t> G4L_trackid;
  for(int it=0;it<=src.NumberOfTracks;++it){
    int pid = src.PIDOfTrack[it];
    if(pid == 3122) G4L_trackid.push_back(it);
  }

  std::vector<Int_t> G4decays_trackid;
  int G4tidScatKm = -1;
  int G4tidScatP = -1;
  int G4tidBeamKm = -1;    
  for(int it=0;it<=src.NumberOfTracks;++it){
    event.PIDOfTrack.push_back(src.PIDOfTrack[it]);
    event.ParentIDOfTrack.push_back(src.ParentIDOfTrack[it]);
    event.VertexOfTrack_x.push_back(src.VertexOfTrack_x[it]);
    event.VertexOfTrack_y.push_back(src.VertexOfTrack_y[it]);
    event.VertexOfTrack_z.push_back(src.VertexOfTrack_z[it]);
    event.MomentumOfTrack.push_back(src.MomentumOfTrack[it]);
    event.MomentumOfTrack_x.push_back(src.MomentumOfTrack_x[it]);
    event.MomentumOfTrack_y.push_back(src.MomentumOfTrack_y[it]);
    event.MomentumOfTrack_z.push_back(src.MomentumOfTrack_z[it]);
    int parent = src.ParentIDOfTrack[it];
    if(it==0) continue;
    if(parent<0) continue;
    int pid = src.PIDOfTrack[it];
    double mom = src.MomentumOfTrack[it]/1000;//MeV to GeV
    double mom_x = src.MomentumOfTrack_x[it]/1000;
    double mom_y = src.MomentumOfTrack_y[it]/1000;
    double mom_z = src.MomentumOfTrack_z[it]/1000;
    double vert_x = src.VertexOfTrack_x[it];
    double vert_y = src.VertexOfTrack_y[it];
    double vert_z = src.VertexOfTrack_z[it];
    if( parent==0){
      G4decays_trackid.push_back(it);
      if(abs(pid)==321){
	if(it==1){
	  event.G4kmid = it;
	  event.G4kmvtx_x = vert_x;
	  event.G4kmvtx_y = vert_y;
	  event.G4kmvtx_z = vert_z;
	  event.G4kmmom = mom;
	  event.G4kmmom_x = -mom_x;
	  event.G4kmmom_y = -mom_y;
	  event.G4kmmom_z = -mom_z;
	  G4tidBeamKm = it;	  	  
	}
	else if(it==3){ // scat Km
	  event.G4scatkmid = it;
	  event.G4scatkmvtx_x = vert_x;
	  event.G4scatkmvtx_y = vert_y;
	  event.G4scatkmvtx_z = vert_z;
	  event.G4scatkmmom = mom;
	  event.G4scatkmmom_x = mom_x;
	  event.G4scatkmmom_y = mom_y;
	  event.G4scatkmmom_z = mom_z;
	  G4tidScatKm = it;	  
	}
      } else if( abs(pid)==2212 ){ // scatP for KP reaction case
	if(it==2){
	  event.G4kpid = it;
	  event.G4kpvtx_x = vert_x;
	  event.G4kpvtx_y = vert_y;
	  event.G4kpvtx_z = vert_z;
	  event.G4kpmom = mom;
	  event.G4kpmom_x = mom_x;
	  event.G4kpmom_y = mom_y;
	  event.G4kpmom_z = mom_z;
	  G4tidScatP = it;	  	  
	}
      }
    }	       
  } 
  vector<int> G4TrackID;
  vector<int> PureHits;

  Double_t pKp = event.G4kpmom;
  Double_t uKp = event.G4kpmom_x/event.G4kpmom_z;
  Double_t vKp = event.G4kpmom_y/event.G4kpmom_z;
  Double_t ptKp = pKp/std::sqrt(1.+uKp*uKp+vKp*vKp);
  TVector3 kp_mom(ptKp*uKp, ptKp*vKp, ptKp);

  Double_t pKm = event.G4kmmom;
  Double_t uKm = event.G4kmmom_x/event.G4kmmom_z;
  Double_t vKm = event.G4kmmom_y/event.G4kmmom_z;
  Double_t ptKm = pKm/std::sqrt(1.+uKm*uKm+vKm*vKm);
  TVector3 km_mom(ptKm*uKm, ptKm*vKm, ptKm);
  
  double thetaKP_rad = kp_mom.Angle(km_mom);
  // Double_t cost = km_mom*kp_mom/(pKm*pKp);
  // Double_t theta = TMath::ACos(cost)*TMath::RadToDeg();
  TLorentzVector LvKm(km_mom, TMath::Hypot(km_mom.Mag(), KaonMass));
  TLorentzVector LvKp(kp_mom, TMath::Hypot(kp_mom.Mag(), ProtonMass));    
  TLorentzVector LvC(0., 0., 0., m12C);
  TLorentzVector LvP(0., 0., 0., ProtonMass);
  TLorentzVector LvRproton = LvKm + LvP - LvKp;
  TLorentzVector LvRc = LvKm + LvC - LvKp;

  // std::cout << " debug " << __FILE__ << " " << __LINE__
  // 	    << " LvKm.M():" << LvKm.M() << " kmmom:" << km_mom.Mag()
  // 	    << " LvKp.M():" << LvKp.M() << " kpmom:" << kp_mom.Mag() 
  // 	    << " LvRproton.M():" << LvRproton.M() << " LvRproton.P():" << LvRproton.P()
  // 	    << " LvRc.M():" << LvRc.M() << " LvRc.P():" << LvRc.P() << std::endl;
  
  Double_t mm_12C = LvRc.Mag(); 
  Double_t binding_energy = m11B + KaonMass - mm_12C; //GeV/c2
  TLorentzVector LvMM = LvKm + LvP - LvKp;
  auto veolcityMM = LvMM.BoostVector();
  TLorentzVector LvKmCM = LvKm;
  TLorentzVector LvKpCM = LvKp;
  LvKmCM.Boost(-veolcityMM);
  LvKpCM.Boost(-veolcityMM);

  event.nKK = 1;
  event.Kflag.push_back(1);
  event.MissMass.push_back((LvKm+LvP-LvKp).M());
  event.MissMassCorr.push_back((LvKm+LvP-LvKp).M());
  event.MissMassCorrDE.push_back((LvKm+LvP-LvKp).M());
  event.MissMassNucl.push_back((LvKm+LvC-LvKp).M());
  event.MissMassNuclCorr.push_back((LvKm+LvC-LvKp).M());
  event.MissMassNuclCorrDE.push_back((LvKm+LvC-LvKp).M());
  event.pOrg.push_back(pKm);
  event.pCalc.push_back(pKm);
  event.pCorr.push_back(pKm);
  event.pCorrDE.push_back(pKm);
  event.ub.push_back(uKm);
  event.vb.push_back(vKm);
  event.us.push_back(uKp);
  event.vs.push_back(vKp);
  event.vtx.push_back(event.G4kpvtx_x);
  event.vty.push_back(event.G4kpvtx_y);
  event.vtz.push_back(event.G4kpvtx_z);
  event.closeDist.push_back(0);
  event.inside.push_back(1);

  event.pKurama.push_back(pKp);
  event.qKurama.push_back(1);
  event.chisqrKurama.push_back(1);
  event.thetaKurama.push_back(kp_mom.Theta()*TMath::RadToDeg());
  event.xtgtKurama.push_back(event.G4kpvtx_x);
  event.ytgtKurama.push_back(event.G4kpvtx_y);
  event.utgtKurama.push_back(uKp);
  event.vtgtKurama.push_back(vKp);
  event.xin.push_back(event.G4kpvtx_x);
  event.yin.push_back(event.G4kpvtx_y);
  event.zin.push_back(event.G4kpvtx_z);
  event.pxin.push_back(event.G4kmmom_x);
  event.pyin.push_back(event.G4kmmom_y);
  event.pzin.push_back(event.G4kmmom_z);
  event.xout.push_back(event.G4kpvtx_x);
  event.yout.push_back(event.G4kpvtx_y);
  event.zout.push_back(event.G4kpvtx_z);
  event.pxout.push_back(event.G4kpmom_x);
  event.pyout.push_back(event.G4kpmom_y);
  event.pzout.push_back(event.G4kpmom_z);

  event.ntK18 = 1;
  event.pK18.push_back(pKm);
  event.chisqrK18.push_back(1);
  event.xtgtK18.push_back(event.G4kmvtx_x);
  event.ytgtK18.push_back(event.G4kmvtx_y);
  event.utgtK18.push_back(uKm);
  event.vtgtK18.push_back(vKm);

  event.ntKurama = 1;
  event.isgoodTPCKurama.push_back(1);
  event.kflagTPCKurama.push_back(1);
  event.chisqrTPCKurama.push_back(1);
  event.pTPCKurama.push_back(pKp);
  event.qTPCKurama.push_back(1);
  event.m2TPCKurama.push_back((LvKp.M2()-KaonMass*KaonMass));
  event.xtgtTPCKurama.push_back(event.G4kpvtx_x);
  event.ytgtTPCKurama.push_back(event.G4kpvtx_y);
  event.utgtTPCKurama.push_back(uKp);
  event.vtgtTPCKurama.push_back(vKp);
  event.thetaTPCKurama.push_back(kp_mom.Theta()*TMath::RadToDeg());
  event.isgoodTPC.push_back(1);
  event.insideTPC.push_back(1);
  event.vtxTPC.push_back(event.G4kpvtx_x);
  event.vtyTPC.push_back(event.G4kpvtx_y);
  event.vtzTPC.push_back(event.G4kpvtx_z-tpc::ZTarget);
  event.closeDistTPC.push_back(0);
  event.MissMassTPC.push_back((LvKm+LvP-LvKp).M());
  event.MissMassCorrTPC.push_back((LvKm+LvP-LvKp).M());
  event.MissMassCorrDETPC.push_back((LvKm+LvP-LvKp).M());
  event.pOrgTPC.push_back(pKm);
  event.pCalcTPC.push_back(pKm);
  event.pCorrTPC.push_back(pKm);
  event.pCorrDETPC.push_back(pKm);
  event.thetaTPC.push_back(thetaKP_rad*TMath::RadToDeg());
  event.thetaCMTPC.push_back(LvKmCM.Theta()*TMath::RadToDeg());
  event.costCMTPC.push_back(cos((LvKmCM.Vect()).Angle(LvKpCM.Vect())));
  event.xbTPC.push_back(event.G4kpvtx_x);
  event.ybTPC.push_back(event.G4kpvtx_y);
  event.ubTPC.push_back(uKm);
  event.vbTPC.push_back(vKm);
  event.xsTPC.push_back(event.G4kpvtx_x);
  event.ysTPC.push_back(event.G4kpvtx_y);
  event.usTPC.push_back(uKp);
  event.vsTPC.push_back(vKp);

  event.BE.resize(event.nKK);
  event.BETPC.resize(event.nKK);
  event.BE_LL.resize(event.nKK);
  event.BETPC_LL.resize(event.nKK);

  event.km_mom_x.resize(event.nKK);
  event.km_mom_y.resize(event.nKK);
  event.km_mom_z.resize(event.nKK);
  event.kp_mom_x.resize(event.nKK);
  event.kp_mom_y.resize(event.nKK);
  event.kp_mom_z.resize(event.nKK);

  HF1( 1, event.status++ ); // debug 0
  //if(std::isnan(event.G4kmmom)) return false;
  if( event.nKK != 1 ) return false;
  if( event.isgoodTPCKurama.size()!=1 ) return false;
  if( event.isgoodTPCKurama[0]!=1 ) return false;
  if( event.insideTPC[0] != 1) return false;
#if kkevent
  if( event.kflagTPCKurama[0]!=1 ) return false;
#endif
  Double_t thetaTPC = event.thetaTPC[0];
  if(CarbonRun){
    if( !(thetaTPC>minThetaKP && thetaTPC<maxThetaKP) ) return true;
  }

  std::cout << __FILE__ << " " << __LINE__ << " " << event.G4kmmom << std::endl;
  
  TVector3 kkvtxTPC(event.vtxTPC[0], event.vtyTPC[0], event.vtzTPC[0] + tpc::ZTarget);
  event.BE[0] = 1000.*binding_energy; //MeV/c2
  Double_t binding_energy_LL = m10Be + 2.*LambdaMass - mm_12C; //GeV/c2
  event.BE_LL[0] = 1000.*binding_energy_LL; //MeV/c2

  event.km_mom_x[0] = km_mom.x();
  event.km_mom_y[0] = km_mom.y();
  event.km_mom_z[0] = km_mom.z();
  event.kp_mom_x[0] = kp_mom.x();
  event.kp_mom_y[0] = kp_mom.y();
  event.kp_mom_z[0] = kp_mom.z();

  TLorentzVector LvRcTPC;
  //TVector3 km_unit = TVector3(event.ubTPC[0], event.vbTPC[0], 1.).Unit();
  //TVector3 km_unit = TVector3(event.utgtK18[0], event.vtgtK18[0], 1.).Unit();
  //TVector3 km_momTPC = km_unit*event.pK18[0];
  TVector3 km_momTPC = km_mom;

  //TVector3 kp_unit = TVector3(event.usTPC[0], event.vsTPC[0], 1.).Unit();
  //TVector3 kp_momTPC = kp_unit*event.pCorrDETPC[0];
  TVector3 kp_momTPC = kp_mom;
  TVector3 miss_momTPC = km_momTPC - kp_momTPC;  

  TLorentzVector LvKmTPC(km_momTPC, TMath::Hypot(km_momTPC.Mag(), KaonMass));
  TLorentzVector LvScatPTPC(kp_momTPC, TMath::Hypot(kp_momTPC.Mag(), ProtonMass));
  TLorentzVector LvCTPC(0., 0., 0., m12C);
  TLorentzVector LvPTPC(0., 0., 0., ProtonMass);
  std::cout << "m12C:" << m12C << " ProtonMass:" << ProtonMass << std::endl;
  TLorentzVector LvScatKmTPC = LvKmTPC + LvPTPC - LvScatPTPC;
  LvRcTPC = LvKmTPC + LvCTPC - LvScatPTPC;

  double mm_12CTPC = LvRcTPC.M();
  //double binding_energyTPC = m11B + KaonMass - (mm_12CTPC - 0.120); //GeV/c2
  double binding_energyTPC = m11B + KaonMass - mm_12CTPC; //GeV/c2
  event.BETPC[0] = binding_energyTPC; //MeV/c2
  
  //event.MissMassCorrDE[0] = event.MissMassCorrDE[0] - 0.120 ;
  event.MissMassCorrDE[0] = event.MissMassCorrDE[0]; 
  double missmass = event.MissMassCorrDE[0];
  if(missmass>minMM&&missmass<maxMM) HF1( 1, event.status++ ); // debug 1
  HF1(3900,event.MissMassCorrDE[0]);
  HF1(13900,-event.BETPC[0]);  
  TVector3 G4SKmMomVec(event.G4scatkmmom_x,event.G4scatkmmom_y,event.G4scatkmmom_z);
  HF1(120, G4SKmMomVec.Phi());
  HF1(121, G4SKmMomVec.Phi()*TMath::RadToDeg());
  
  Double_t pionmip = Kinematics::HypTPCdEdx(3, PionMass*1000., 1.8/TMath::Hypot(PionMass, 1.8)); //MeV/c2
  Double_t HTOF_thr = pionmip*0.1; // 10% of 1.8 GeV/c pi- mip  

  Int_t nhHtof = 0;
  std::vector<Double_t> HtofSeg;
  std::vector<Double_t> tHtof;
  std::vector<Double_t> dtHtof;
  std::vector<Double_t> deHtof;
  std::vector<Double_t> posHtof;
  std::vector<Int_t> G4tidHtof;
  for(int ih=0;ih<**src.nhHtof;++ih){ //decays from Xi-, L and so on
    Double_t seg = ((std::vector<Double_t>)**src.HtofSeg).at(ih);
    Double_t t = ((std::vector<Double_t>)**src.tHtof).at(ih);
    //Double_t dt = ((std::vector<Double_t>)**src.dtHtof).at(ih);
    Double_t de = ((std::vector<Double_t>)**src.deHtof).at(ih);
    Double_t pos = ((std::vector<Double_t>)**src.posHtof).at(ih);
    Int_t G4tid = ((std::vector<Int_t>)**src.G4tidHtof).at(ih);
    auto iter = find(G4decays_trackid.begin(), G4decays_trackid.end(), G4tid);
    if(iter != G4decays_trackid.end()){
      nhHtof++;
      HtofSeg.push_back(seg);
      tHtof.push_back(t);
      //dtHtof.push_back(dt);
      deHtof.push_back(de);
      posHtof.push_back(pos);
      G4tidHtof.push_back(G4tid);
    }
  }
  for(int ih=0;ih<**src.nhHtof;++ih){ //daughters of decays from Xi-, L and so on
    Double_t seg = ((std::vector<Double_t>)**src.HtofSeg).at(ih);
    Double_t t = ((std::vector<Double_t>)**src.tHtof).at(ih);
    //Double_t dt = ((std::vector<Double_t>)**src.dtHtof).at(ih);
    Double_t de = ((std::vector<Double_t>)**src.deHtof).at(ih);
    Double_t pos = ((std::vector<Double_t>)**src.posHtof).at(ih);
    Int_t G4tid = ((std::vector<Int_t>)**src.G4tidHtof).at(ih);
    auto iter = find(G4decays_trackid.begin(), G4decays_trackid.end(), G4tid);
    if(iter != G4decays_trackid.end()) continue;
    if(de < 0.01) continue; //veto neutral particle
    Int_t parentid = event.ParentIDOfTrack[G4tid];
    if(parentid < 0) continue;
    auto iter2 = find(G4decays_trackid.begin(), G4decays_trackid.end(), parentid);
    if(iter2 == G4decays_trackid.end()) continue; //should be a daughter of decays from Xi-, L and so on
    Int_t parentpid = event.ParentIDOfTrack[parentid];
    if(parentpid==2212 || parentpid==321) continue;
    Int_t pid = event.PIDOfTrack[G4tid];

    auto iter3 = find(G4tidHtof.begin(), G4tidHtof.end(), parentid);
    if(iter3 == G4decays_trackid.end()){ //the daughter hits the HTOF
      nhHtof++;
      HtofSeg.push_back(seg);
      tHtof.push_back(t);
      //dtHtof.push_back(dt);
      deHtof.push_back(de);
      posHtof.push_back(pos);
      G4tidHtof.push_back(G4tid);
    }
    else{ //the daughter decays in the HTOF
      for(int id=0;id<nhHtof;++id){
	if(G4tidHtof[id] == parentid){
	  deHtof[id] += de; //add energy mother and daughter's hits

	  nhHtof++;
	  HtofSeg.push_back(seg);
	  tHtof.push_back(t);
	  //dtHtof.push_back(dt);
	  deHtof.push_back(de);
	  posHtof.push_back(pos);
	  G4tidHtof.push_back(G4tid);
	}
      }
    }
  }

  //Filling hits above dE threshold
  for(int id=0;id<nhHtof;++id){
    std::cout << __FILE__ << " " << __LINE__ << " HTOF dE " << deHtof[id] << std::endl;
    if(deHtof[id] < HTOF_thr) continue;
    event.nhHtof++;
    event.HtofSeg.push_back(HtofSeg[id]);
    event.tHtof.push_back(tHtof[id]);
    //event.dtHtof.push_back(dtHtof[id]);
    event.deHtof.push_back(deHtof[id]);
    event.posHtof.push_back(posHtof[id]);
    event.G4tidHtof.push_back(G4tidHtof[id]);
  }

  /*
    event.nhHtof = **src.nhHtof;
    event.HtofSeg = **src.HtofSeg;
    event.tHtof = **src.tHtof;
    //event.dtHtof = **src.dtHtof;
    event.deHtof = **src.deHtof;
    event.posHtof = **src.posHtof;
    event.G4tidHtof = **src.G4tidHtof;
  */
  
  Int_t ntTpc = **src.ntTpc;
  if( ntTpc == 0 ) return true;
  if(missmass>minMM&&missmass<maxMM) HF1( 1, event.status++ ); // debug 2
  
  event.nclTpc = **src.nclTpc;
  event.remain_nclTpc = **src.remain_nclTpc;
#if SaveRawData
  event.cluster_x = **src.cluster_x;
  event.cluster_y = **src.cluster_y;
  event.cluster_z = **src.cluster_z;
  event.cluster_de = **src.cluster_de;
  event.cluster_size = **src.cluster_size;
  event.cluster_layer = **src.cluster_layer;
  event.cluster_mrow = **src.cluster_mrow;
  event.cluster_row_center = **src.cluster_row_center;
  event.cluster_houghflag = **src.cluster_houghflag;
  //  event.cluster_G4tid = **src.cluster_G4tid;

  event.remain_cluster_x.resize(event.remain_nclTpc);
  event.remain_cluster_y.resize(event.remain_nclTpc);
  event.remain_cluster_z.resize(event.remain_nclTpc);
  event.remain_cluster_de.resize(event.remain_nclTpc);
  event.remain_cluster_size.resize(event.remain_nclTpc);
  event.remain_cluster_layer.resize(event.remain_nclTpc);
  event.remain_cluster_mrow.resize(event.remain_nclTpc);
  event.remain_cluster_de_center.resize(event.remain_nclTpc);
  event.remain_cluster_x_center.resize(event.remain_nclTpc);
  event.remain_cluster_y_center.resize(event.remain_nclTpc);
  event.remain_cluster_z_center.resize(event.remain_nclTpc);
  event.remain_cluster_row_center.resize(event.remain_nclTpc);
  event.remain_cluster_houghflag.resize(event.remain_nclTpc);
  Int_t icl_remain = 0;
  for( Int_t icl=0; icl<event.nclTpc; ++icl ){
    if(event.cluster_houghflag[icl]!=0) continue;
    event.remain_cluster_x[icl_remain] = event.cluster_x[icl];
    event.remain_cluster_y[icl_remain] = event.cluster_y[icl];
    event.remain_cluster_z[icl_remain] = event.cluster_z[icl];
    event.remain_cluster_de[icl_remain] = event.cluster_de[icl];
    event.remain_cluster_size[icl_remain] = event.cluster_size[icl];
    event.remain_cluster_layer[icl_remain] = event.cluster_layer[icl];
    event.remain_cluster_mrow[icl_remain] = event.cluster_mrow[icl];
    event.remain_cluster_houghflag[icl_remain] = event.cluster_houghflag[icl];
    icl_remain++;
  }
#endif

  event.ntTpc = ntTpc;
  event.nhtrack = **src.nhtrack;
  event.trackid = **src.trackid;
  event.isXi = **src.isXi;
  event.isBeam = **src.isBeam;
  event.isKurama = **src.isKurama;
  event.isK18 = **src.isK18;
  event.isAccidental = **src.isAccidental;
  event.isMultiloop = **src.isMultiloop;
  event.charge = **src.charge;
  event.pid = **src.pid;
  event.purity = **src.purity;
  event.G4tid = **src.G4tid;
  for(auto t:event.G4tid){
    event.G4pid.push_back(src.PIDOfTrack[t]);
  }

  event.efficiency = **src.efficiency;
  event.chisqr = **src.chisqr;
  event.pval = **src.pval;
  event.helix_cx = **src.helix_cx;
  event.helix_cy = **src.helix_cy;
  event.helix_z0 = **src.helix_z0;
  event.helix_r = **src.helix_r;
  event.helix_dz = **src.helix_dz;
  event.dE = **src.dE;
  event.dEdx = **src.dEdx;  
  event.mom0 = **src.mom0;
  event.path = **src.path;
  event.isElectron = **src.isElectron;
  event.nsigma_triton = **src.nsigma_triton;
  event.nsigma_deutron = **src.nsigma_deutron;
  event.nsigma_proton = **src.nsigma_proton;
  event.nsigma_kaon = **src.nsigma_kaon;
  event.nsigma_pion = **src.nsigma_pion;
  event.nsigma_electron = **src.nsigma_electron;

  event.hitlayer = **src.hitlayer;
  event.hitpos_x = **src.hitpos_x;
  event.hitpos_y = **src.hitpos_y;
  event.hitpos_z = **src.hitpos_z;
  event.calpos_x = **src.calpos_x;
  event.calpos_y = **src.calpos_y;
  event.calpos_z = **src.calpos_z;
  event.mom_x = **src.mom_x;
  event.mom_y = **src.mom_y;
  event.mom_z = **src.mom_z;
  event.residual = **src.residual;
  event.residual_x = **src.residual_x;
  event.residual_y = **src.residual_y;
  event.residual_z = **src.residual_z;
  event.resolution_x = **src.resolution_x;
  event.resolution_y = **src.resolution_y;
  event.resolution_z = **src.resolution_z;
  event.helix_t = **src.helix_t;
  event.pathhit = **src.pathhit;
  event.alpha = **src.alpha;
  event.track_cluster_de = **src.track_cluster_de;
  event.track_cluster_size = **src.track_cluster_size;
  event.track_cluster_mrow = **src.track_cluster_mrow;

  event.chargeIndistinguishable = **src.chargeIndistinguishable;
  event.pid_inverted = **src.pid_inverted;
  event.chisqr_inverted = **src.chisqr_inverted;
  event.pval_inverted = **src.pval_inverted;
  event.helix_cx_inverted = **src.helix_cx_inverted;
  event.helix_cy_inverted = **src.helix_cy_inverted;
  event.helix_z0_inverted = **src.helix_z0_inverted;
  event.helix_r_inverted = **src.helix_r_inverted;
  event.helix_dz_inverted = **src.helix_dz_inverted;
  event.mom0_inverted = **src.mom0_inverted;

  event.nvtxTpc = **src.nvtxTpc;
  event.vtx_x = **src.vtx_x;
  event.vtx_y = **src.vtx_y;
  event.vtx_z = **src.vtx_z;
  event.vtx_dist = **src.vtx_dist;
  event.vtx_angle = **src.vtx_angle;
  event.vtxid = **src.vtxid;
  event.vtxmom_theta = **src.vtxmom_theta;
  event.vtxpos_x = **src.vtxpos_x;
  event.vtxpos_y = **src.vtxpos_y;
  event.vtxpos_z = **src.vtxpos_z;
  event.vtxmom_x = **src.vtxmom_x;
  event.vtxmom_y = **src.vtxmom_y;
  event.vtxmom_z = **src.vtxmom_z;

  event.nvtxTpcClustered = **src.nvtxTpcClustered;
  event.clusteredVtx_x = **src.clusteredVtx_x;
  event.clusteredVtx_y = **src.clusteredVtx_y;
  event.clusteredVtx_z = **src.clusteredVtx_z;
  event.clusteredVtxid = **src.clusteredVtxid;

  Int_t id_kp = -1;
  for(Int_t k=0;k<ntTpc;k++){
    if(event.isKurama[k]==1) id_kp = k;
  }

  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCTracks(**src.ntTpc, **src.isK18, **src.isKurama,
			 **src.charge, **src.nhtrack, **src.helix_cx,
			 **src.helix_cy, **src.helix_z0, **src.helix_r,
			 **src.helix_dz, **src.hitlayer, **src.track_cluster_mrow,
			 **src.helix_t, **src.track_cluster_de, **src.resolution_x,
			 **src.resolution_y, **src.resolution_z, **src.hitpos_x,
			 **src.hitpos_y, **src.hitpos_z);

  HypTPCTask& GFTrackCont = HypTPCTask::GetInstance();
  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix(it);
    if( !tp ) continue;
    // smear dEdx
    double oneSigma = HypTPCdEdx1sigmaKaon(event.mom0[it]);
    event.dEdx[it] = gRandom->Gaus(event.dEdx[it], dedx_smear_nsigma*oneSigma);
    //event.pid[it] = tp -> GetPid();
    event.pid[it] = Kinematics::HypTPCdEdxPID(event.dEdx[it],event.mom0[it]);
    event.isElectron[it] = Kinematics::HypTPCdEdxElectron(event.dEdx[it], event.mom0[it]);
    event.nsigma_triton[it] = Kinematics::HypTPCdEdxNsigmaTriton(event.dEdx[it], event.mom0[it]);
    event.nsigma_deutron[it] = Kinematics::HypTPCdEdxNsigmaDeutron(event.dEdx[it], event.mom0[it]);
    event.nsigma_proton[it] = Kinematics::HypTPCdEdxNsigmaProton(event.dEdx[it], event.mom0[it]);
    event.nsigma_kaon[it]  = Kinematics::HypTPCdEdxNsigmaKaon(event.dEdx[it], event.mom0[it]);
    event.nsigma_pion[it] = Kinematics::HypTPCdEdxNsigmaPion(event.dEdx[it], event.mom0[it]);
    event.nsigma_electron[it] = Kinematics::HypTPCdEdxNsigmaElectron(event.dEdx[it], event.mom0[it]);

    if(event.isKurama[it]==1) GFTrackCont.AddHelixTrack(321, tp);
    else if(event.isElectron[it]==1) GFTrackCont.AddHelixTrack(event.charge[it]*(-11), tp);
    else{
      std::vector<Int_t> pdgcode;
      Kinematics::HypTPCPID_PDGCode(event.charge[it], event.pid[it], pdgcode);
      if((event.pid[it]&4)!=4 && TMath::Abs(event.nsigma_electron[it]) < 3)
	pdgcode.push_back(event.charge[it]*(-11));
      GFTrackCont.AddHelixTrack(pdgcode, tp);
    }

    vector<TVector3> TPCHit;
    for(int ih=0;ih<event.hitpos_x.at(it).size();++ih){
      TPCHit.push_back(TVector3(
				event.hitpos_x.at(it).at(ih),
				event.hitpos_y.at(it).at(ih),
				event.hitpos_z.at(it).at(ih)));
    }
    int nPureHit;
    int G4tid = TPCToG4TrackID(TPCHit,src.nhittpc,src.ititpc,src.xtpc,src.ytpc,src.ztpc,nPureHit);
    G4TrackID.push_back(G4tid);
    PureHits.push_back(nPureHit);
    if(G4tid<0) continue;
#if LLRecon
    if(G4tid == event.G4p1id){
      event.p1_tracked = true;
      event.p1t_mom0 = event.mom0[it];
    }
    if(G4tid == event.G4p2id){
      event.p2_tracked = true;
      event.p2t_mom0 = event.mom0[it];
    }
#elif XiRecon
    if(G4tid == event.G4protonid){
      event.p_tracked = true;
      event.pt_mom0 = event.mom0[it];
    }
    if(event.G4extraprotonid !=-1 && G4tid == event.G4extraprotonid){
      event.extrap_tracked = true;
    }
#endif
    if(G4tid == event.G4pi1id){
      event.pi1_tracked = true;
      event.pi1t_mom0 = event.mom0[it];
    }
    if(G4tid == event.G4pi2id){
      event.pi2_tracked = true;
      event.pi2t_mom0 = event.mom0[it];
    }
    if(G4tid == G4tidScatKm){
      HF2(11, event.charge[it]*event.mom0[it], event.dEdx[it]);
      HF1(24, event.dEdx[it]);      
    }
    if(G4tid == G4tidBeamKm){
      std::cout << " BeamKm isK18=" << event.isBeam[it] << std::endl;
    }
    if(G4tid == G4tidScatP){
      std::cout << " ScatP isKurama=" << event.isKurama[it] << std::endl; 
      event.isKurama[it] = 1;
      // DstTPCTrackingHelixGeant5 shoud be debugged because scat proton is not assigned as isKurama=1
    }
  }

  GFTrackCont.FitTracks();
  if(missmass>minMM&&missmass<maxMM) HF1( 1, event.status++ ); // debug 3

  Int_t GFntTpc = GFTrackCont.GetNTrack();
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    return true;
  }
  if(missmass>minMM&&missmass<maxMM) HF1( 1, event.status++ ); // debug 4
  
  event.GFcharge.resize(GFntTpc);
  event.GFchisqr.resize(GFntTpc);
  event.GFtof.resize(GFntTpc);
  //event.GFtracklen.resize(GFntTpc);
  event.GFpval.resize(GFntTpc);
  event.GFfitstatus.resize(GFntTpc);
  event.GFpdgcode.resize(GFntTpc);
  event.GFnhtrack.resize(GFntTpc);
  event.GFlayer.resize(GFntTpc);
  event.GFpos_x.resize(GFntTpc);
  event.GFpos_y.resize(GFntTpc);
  event.GFpos_z.resize(GFntTpc);
  event.GFmom.resize(GFntTpc);
  event.GFmom_x.resize(GFntTpc);
  event.GFmom_y.resize(GFntTpc);
  event.GFmom_z.resize(GFntTpc);
  event.GFresidual_x.resize(GFntTpc);
  event.GFresidual_y.resize(GFntTpc);
  event.GFresidual_z.resize(GFntTpc);
  event.GFresidual_p.resize(GFntTpc);
  event.GFresidual_px.resize(GFntTpc);
  event.GFresidual_py.resize(GFntTpc);
  event.GFresidual_pz.resize(GFntTpc);
  event.GFinside.resize(GFntTpc);

  event.GFtracklen.resize(GFntTpc);
  event.GFtrack2vtxdist.resize(GFntTpc);
  event.GFcalctof.resize(GFntTpc);
  event.GFposx.resize(GFntTpc);
  event.GFposy.resize(GFntTpc);
  event.GFposz.resize(GFntTpc);
  event.GFsegHtof.resize(GFntTpc);
  event.GFtofHtof.resize(GFntTpc);
  event.GFtdiffHtof.resize(GFntTpc);
  event.GFposHtof.resize(GFntTpc);
  event.GFinvbeta.resize(GFntTpc);
  event.GFm2.resize(GFntTpc);
  event.nsigma_tritonHtof.resize(ntTpc);
  event.nsigma_deutronHtof.resize(ntTpc);
  event.nsigma_protonHtof.resize(ntTpc);
  event.nsigma_kaonHtof.resize(ntTpc);
  event.nsigma_pionHtof.resize(ntTpc);
  event.nsigma_electronHtof.resize(ntTpc);

  Int_t ntrack_intarget = 0;
  Double_t x0[100] = {0};
  Double_t y0[100] = {0};
  Double_t u0[100] = {0};
  Double_t v0[100] = {0};
  for( Int_t igf=0; igf<GFntTpc; ++igf ){
    event.GFfitstatus[igf] = (int)GFTrackCont.TrackCheck(igf);
    if(!GFTrackCont.TrackCheck(igf)) continue;
    int nh = GFTrackCont.GetNHits(igf);
    event.GFlayer[igf].resize(nh);
    event.GFpos_x[igf].resize(nh);
    event.GFpos_y[igf].resize(nh);
    event.GFpos_z[igf].resize(nh);
    event.GFmom[igf].resize(nh);
    event.GFmom_x[igf].resize(nh);
    event.GFmom_y[igf].resize(nh);
    event.GFmom_z[igf].resize(nh);
    event.GFresidual_x[igf].resize(nh);
    event.GFresidual_y[igf].resize(nh);
    event.GFresidual_z[igf].resize(nh);
    event.GFresidual_p[igf].resize(nh);
    event.GFresidual_px[igf].resize(nh);
    event.GFresidual_py[igf].resize(nh);
    event.GFresidual_pz[igf].resize(nh);

    event.GFchisqr[igf] = GFTrackCont.GetChi2NDF(igf);
    event.GFcharge[igf] = GFTrackCont.GetCharge(igf);
    event.GFtof[igf] = GFTrackCont.GetTrackTOF(igf, 0, -1);
    event.GFpval[igf] = GFTrackCont.GetPvalue(igf);
    event.GFnhtrack[igf] = GFTrackCont.GetNHits(igf);
    event.GFpdgcode[igf] = GFTrackCont.GetPDGcode(igf);
    for( Int_t ihit=0; ihit<nh; ++ihit ){
      TVector3 hit = GFTrackCont.GetPos(igf, ihit);
      TVector3 mom = GFTrackCont.GetMom(igf, ihit);
      Int_t layer = (int)event.hitlayer[igf][ihit];
      event.GFlayer[igf][ihit] = layer;
      event.GFmom_x[igf][ihit] = mom.x();
      event.GFmom_y[igf][ihit] = mom.y();
      event.GFmom_z[igf][ihit] = mom.z();
      event.GFmom[igf][ihit] = mom.Mag();
      event.GFpos_x[igf][ihit] = hit.x();
      event.GFpos_y[igf][ihit] = hit.y();
      event.GFpos_z[igf][ihit] = hit.z();

      event.GFresidual_x[igf][ihit] = hit.x() - event.hitpos_x[igf][ihit];
      event.GFresidual_y[igf][ihit] = hit.y() - event.hitpos_y[igf][ihit];
      event.GFresidual_z[igf][ihit] = hit.z() - event.hitpos_z[igf][ihit];

      double chargetest = event.GFcharge[igf]*event.charge[igf];
      event.GFresidual_p[igf][ihit] = mom.Mag() - event.mom0[igf];
      event.GFresidual_px[igf][ihit] = mom.x() - chargetest*event.mom_x[igf][ihit];
      event.GFresidual_py[igf][ihit] = mom.y() - chargetest*event.mom_y[igf][ihit];
      event.GFresidual_pz[igf][ihit] = mom.z() - chargetest*event.mom_z[igf][ihit];
    }//ih
    if(event.isBeam[igf]==1) continue;
    if(event.isK18[igf]==1) continue;
    if(event.isAccidental[igf]==1) continue;
    if(GFTrackCont.IsInsideTarget(igf)){
      event.GFinside[igf] = 1;
      TVector3 posv; TVector3 momv; double len; double tof;
      if(GFTrackCont.ExtrapolateToTargetCenter(igf, posv, momv, len, tof)){
	x0[ntrack_intarget] = posv.x();
	y0[ntrack_intarget] = posv.y();
	u0[ntrack_intarget] = momv.x()/momv.z();
	v0[ntrack_intarget] = momv.y()/momv.z();
	ntrack_intarget++;
      }
    }
    else event.GFinside[igf] = 0;

    TVector3 vertex = Kinematics::MultitrackVertex(ntrack_intarget, x0, y0, u0, v0);
    event.GFntTpc_inside = ntrack_intarget;
    event.GFprodvtx_x = vertex.x();
    event.GFprodvtx_y = vertex.y();
    event.GFprodvtx_z = vertex.z();
    {
      Int_t repid = -1;
      Int_t hitid_htof; Double_t tof; Double_t len;
      TVector3 pos_htof; Double_t track2tgt_dist;
      Bool_t htofextrapolation =
	GFTrackCont.TPCHTOFTrackMatching(igf, repid, vertex,
					 event.HtofSeg, event.posHtof,
					 hitid_htof, tof,
					 len, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	event.GFtracklen[igf] = len;
	event.GFtrack2vtxdist[igf] = track2tgt_dist;
	event.GFcalctof[igf] = tof;
	event.GFposx[igf] = pos_htof.x();
	event.GFposy[igf] = pos_htof.y();
	event.GFposz[igf] = pos_htof.z();
	event.GFsegHtof[igf] = event.HtofSeg[hitid_htof];
	event.GFtofHtof[igf] = event.tHtof[hitid_htof];
	event.GFposHtof[igf] = event.posHtof[hitid_htof];

	Double_t beta = len/event.tHtof[hitid_htof]/MathTools::C();
	event.GFinvbeta[igf] = 1./beta;
	Double_t mass2 = Kinematics::MassSquare(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	event.GFm2[igf] = mass2;

	event.nsigma_tritonHtof[igf] = Kinematics::HypTPCHTOFNsigmaTriton(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	event.nsigma_deutronHtof[igf] = Kinematics::HypTPCHTOFNsigmaDeutron(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	event.nsigma_protonHtof[igf] = Kinematics::HypTPCHTOFNsigmaProton(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	event.nsigma_kaonHtof[igf] = Kinematics::HypTPCHTOFNsigmaKaon(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	event.nsigma_pionHtof[igf] = Kinematics::HypTPCHTOFNsigmaPion(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	event.nsigma_electronHtof[igf] = Kinematics::HypTPCHTOFNsigmaElectron(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
      }
    } //common
  }//igf

  std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;  

  std::vector<Int_t> k_id_container;
  std::vector<Int_t> k_repid_container;
  std::vector<TVector3> k_mom_container, k_vert_container;  
  std::vector<Double_t> L_targetdist_container;
  std::vector<TVector3> L_targetvtx_container;
  std::vector<Int_t>    k_htofextrap_container;
  std::vector<Int_t>    k_htofhitid_container;
  std::vector<Int_t>    k_htofseg_container;
  std::vector<Double_t> k_mass2_container;
  std::vector<Double_t> k_posY_container;
  std::vector<TVector3> k_htofpos_container;
  std::vector<Double_t> k_diffY_container;
  std::vector<Double_t> k_tracklen_container;

  Int_t numK18=0; Int_t numKurama=0; Int_t numBeam=0; Int_t numAcc=0;
  Int_t numKm=0; Int_t numPim=0; Int_t numPip=0; Int_t numP=0; Int_t numPPip=0;
  Int_t numEp=0; Int_t numEm=0;
  Int_t idKm=-1;
  bool kmflag_dedxpid = false;  
  {
    for(Int_t it=0;it<ntTpc;it++){ // kaon
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;        
      if( event.isK18[it]==1 ){
	numK18++;
      } else if( event.isKurama[it]==1 ) {
	numKurama++;
      } else if( event.isBeam[it]==1 ) {
	numBeam++;
      } else if( event.isAccidental[it]==1 ) {
	numAcc++;
      }
      if(event.isElectron[it]==1) continue;
      if(event.isK18[it]==1) continue;
      if(event.isKurama[it]==1) continue;
      if(event.isBeam[it]==1) continue;
      if(event.isAccidental[it]==1) continue;
      if(event.GFinside[it]!=1) continue;
      
      // if((event.pid[it]&2)==2 && event.charge[it]==-1){ //k-
      if(event.charge[it]==-1){ //k-
	if( event.nsigma_kaon[it]>mindEdxSigKaon && event.nsigma_kaon[it]<maxdEdxSigKaon ){       
	  numKm++;
	  idKm = it;
	  continue;
	}
      }
      if(event.isElectron[it]==1){ //e+, e-
	if(event.charge[it]==1){
	  numEp++;
	  continue;
	} else {
	  numEm++;
	  continue;
	}
      }    
      else if((event.pid[it]&4)==4 && (event.pid[it]&1)!=1 && event.charge[it]==1){ //proton
	numP++;
	continue;      
      }
      else if((event.pid[it]&1)==1 && event.charge[it]==-1){ //pi-
	Double_t slope = event.helix_dz[it];
	Double_t helixmom = event.mom0[it];
	Int_t pi_nh = event.helix_t[it].size();      
	TVector3 pi_start = TVector3(event.calpos_x[it][0],
				     event.calpos_y[it][0],
				     event.calpos_z[it][0]);
	TVector3 pi_end = TVector3(event.calpos_x[it][pi_nh-1],
				   event.calpos_y[it][pi_nh-1],
				   event.calpos_z[it][pi_nh-1]);
	double pi_vertex_dist=-999.;
	if(TMath::Abs(slope)<0.05 && TMath::Abs(helixmom)>0.5 &&
	   !(Kinematics::HelixDirection(tgtpos,pi_start,pi_end,pi_vertex_dist))&&
	   (pi_start.x()-pi_end.x())>-10 && (pi_start.x()-pi_end.x())<50.){
	  event.isAccidental[it] = 1;
	  numAcc++;
	  //target_accidental_id_container.push_back(it); //Accidental beam on the target
	  continue; //Accidental K-
	}      
	numPim++;
	continue;
      }
      else if((event.pid[it]&4)!=4 && (event.pid[it]&1)==1 && event.charge[it]==1){ //pi+
	numPip++;
	continue;
      }
      else if(((event.pid[it]&4)==4 || (event.pid[it]&1)==1) && event.charge[it]==1){ //p or pi+ with high-mom
	if(event.isKurama[it]){
	  std::cout << " Proton or pion with high mom measured by Kurama" << std::endl;
	}
	numPPip++;
	continue;
      }
    }
  }

  std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;    

  int l_candidates=1;
  std::vector<Double_t> GFk_targetdist_container(l_candidates, qnan);
  std::vector<Double_t> GFk_tof_container(l_candidates, qnan);
  std::vector<Double_t> GFk_ctau_container(l_candidates, qnan);
  std::vector<TVector3> GFk_targetvtx_container(l_candidates, qnan_vec);
  std::vector<Double_t> GFk_targetcenterdist_container(l_candidates, qnan);
  std::vector<TVector3> GFk_targetcentervtx_container(l_candidates, qnan_vec);
  std::vector<TVector3> GFk_mom_container(l_candidates, qnan_vec);
  std::vector<TVector3> GFk_vert_container(l_candidates, qnan_vec);
  std::vector<Int_t> GFk_id_container(l_candidates, qnan);
  std::vector<Int_t> GFk_repid_container(l_candidates, qnan);
  std::vector<Int_t> GFk_extrap_container(l_candidates, qnan);
  std::vector<Int_t> GFk_htofhitid_container(l_candidates, qnan);
  std::vector<Int_t> GFk_htofseg_container(l_candidates, qnan);
  std::vector<Double_t> GFk_tracklen_container(l_candidates, qnan);
  std::vector<TVector3> GFk_poshtof_container(l_candidates, qnan_vec);
  std::vector<Double_t> GFk_mass2_container(l_candidates, qnan);
  std::vector<Double_t> GFk_invbeta_container(l_candidates, qnan);  

  if(idKm>0){
    std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;      
    std::cout << " total TPC track:"
	      << ntTpc 
	      << " numEm:" << numEm << " numEp:" << numEp
	      << " numPPip:" << numPPip << " numPip:" << numPip
	      << " numP:" << numP << " numKm:" << numKm << std::endl;
  }
  if( numPPip==0&&numPip==0&&numPim==0&&numP==0&&numEm==0&&numEp==0 ){
    std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;
    if(missmass>minMM&&missmass<maxMM) HF1( 1, event.status++ ); // debug 5
    for(int it=0; it<ntTpc; it++){
      if(idKm!=it) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;      
      kmflag_dedxpid=true;
      if ( !event.GFfitstatus[it] ) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;      
      if ( event.isElectron[it]==1 ) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;      
      if ( event.isK18[it]==1 ) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;      
      if ( event.isKurama[it]==1 ) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;      
      if ( event.isBeam[it]==1 ) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;      
      if ( event.isAccidental[it]==1 ) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;      
      Int_t repid_km = -1;
      //if(!GFTrackCont.TrackCheck(it, repid_km)) continue;
      if(event.GFinside[it]!=1) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;      
      Int_t km_nh = event.helix_t[it].size();
      TVector3 km_start = TVector3(event.calpos_x[it][0], event.calpos_y[it][0], event.calpos_z[it][0]);
      TVector3 km_end = TVector3(event.calpos_x[it][km_nh-1], event.calpos_y[it][km_nh-1], event.calpos_z[it][km_nh-1]);
      double vertex_dist = 0.;
      if(!Kinematics::HelixDirection(tgtpos,km_start,km_end,vertex_dist)) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;      
      TVector3 post; TVector3 momKm; double lentgt; double toftgt;
      Int_t extrap=-1;
      if(GFTrackCont.ExtrapolateToTargetCenter(it, post, momKm, lentgt, toftgt)){
      }
      if(momKm.Mag()<0.01) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;      
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;          
      Int_t hitid_htof; Double_t tof_htof; Double_t tracklen_htof; TVector3 pos_htof; Double_t track2tgt_dist; Int_t htofseg;
      Bool_t km_htofextrap =
	GFTrackCont.TPCHTOFTrackMatching(it, repid_km, tgtpos,
					 event.HtofSeg, event.posHtof,
					 hitid_htof, tof_htof,
					 tracklen_htof, pos_htof, track2tgt_dist);
      if(km_htofextrap){
	std::cout << " debug " << __FILE__ << " " << __LINE__ << " HTOF extrap " << std::endl;
	GFk_htofhitid_container[0] = hitid_htof;
	GFk_htofseg_container[0] = event.HtofSeg[hitid_htof];
	GFk_tracklen_container[0] = tracklen_htof;
	GFk_poshtof_container[0] = pos_htof;		
	GFk_tof_container[0] = gRandom->Gaus(event.tHtof[hitid_htof], HTOF_TIME_RESOLUTION_SIGMA);
	GFk_mass2_container[0] =
	  Kinematics::MassSquare(momKm.Mag(), tracklen_htof, GFk_tof_container[0]);
	double kmmass2 = GFk_mass2_container[0];
	GFk_mass2_container[0] = gRandom->Gaus(kmmass2, m2_smear);	
	GFk_invbeta_container[0] =
	  MathTools::C()*event.tHtof[hitid_htof]/tracklen_htof;
      }
      GFk_id_container[0]=idKm;
      GFk_repid_container[0]=repid_km;
      GFk_mom_container[0]=momKm;
      GFk_targetvtx_container[0]=post;
      GFk_targetcentervtx_container[0]=post-tgtpos;
      TVector3 dist = post-tgtpos;
      GFk_targetcenterdist_container[0]=dist.Mag();
      event.kmflag = true;            
    }
  }
  
  if(kmflag_dedxpid){
    if(missmass>minMM&&missmass<maxMM) HF1( 1, event.status++ ); // debug 6
    HF2(10, event.charge[idKm]*event.mom0[idKm], event.dEdx[idKm]);    
    //HF1(3950,event.MissMassCorrDE[0]);
  }
  
  if(event.kmflag){
    std::cout << " debug " << __FILE__ << " " << __LINE__ << " KmFlag " << std::endl;    
    if(missmass>minMM&&missmass<maxMM) HF1( 1, event.status++ ); // debug 7
    int id=0;
    event.GFkmid = GFk_id_container[id];  
    event.GFkmmom = GFk_mom_container[id].Mag();
    event.GFkmmom_x = GFk_mom_container[id].x();
    event.GFkmmom_y = GFk_mom_container[id].y();
    event.GFkmmom_z = GFk_mom_container[id].z();

    TVector3 scat_km_mom(event.GFkmmom_x, event.GFkmmom_y, event.GFkmmom_z);
    TVector3 diff_kmmom = scat_km_mom - miss_momTPC;    
    event.GFkmtheta = scat_km_mom.Theta(); // rad
    event.GFkmphi   = scat_km_mom.Phi();   // rad
    double GFkmCosTheta = scat_km_mom.CosTheta();
    
    event.GFkmtarget_dist = GFk_targetdist_container[id];
    event.GFkmtargetvtx_x = GFk_targetvtx_container[id].x();
    event.GFkmtargetvtx_y = GFk_targetvtx_container[id].y();
    event.GFkmtargetvtx_z = GFk_targetvtx_container[id].z();
    event.GFkmtargetcenter_dist = GFk_targetcenterdist_container[id];
    event.GFkmtargetcenter_x = GFk_targetcentervtx_container[id].x();
    event.GFkmtargetcenter_y = GFk_targetcentervtx_container[id].y();
    event.GFkmtargetcenter_z = GFk_targetcentervtx_container[id].z();
  
    event.GFkmhtofid = GFk_htofhitid_container[id];
    event.GFkmhtofseg = GFk_htofseg_container[id];
    event.GFkmposHtof.push_back(GFk_poshtof_container[id].x());
    event.GFkmposHtof.push_back(GFk_poshtof_container[id].y());
    event.GFkmposHtof.push_back(GFk_poshtof_container[id].z());
    event.GFkmmass2 = GFk_mass2_container[id];
    event.GFkmtracklen = GFk_tracklen_container[id];
    event.GFkminvbeta = GFk_invbeta_container[id];
    event.GFkmtof = GFk_tof_container[id];
    HF1(3950, event.MissMassCorrDE[0]);
    HF1(13950,-event.BETPC[0]);

    if(event.GFkmmom>0.01){
      std::cout << " debug " << __FILE__ << " " << __LINE__ << " GFkmmom>0.01 " << std::endl;
      HF1(20, event.GFkmmass2);
      std::cout << " debug " << __FILE__ << " " << __LINE__
		<< " GFkmmass2 " << event.GFkmmass2  << std::endl;
      HF1(21, event.GFkmtracklen);
      HF1(22, event.GFkmtof);
      HF1(23, event.dEdx[event.GFkmid]);
      HF1(3050, event.GFkmmom);
      HF1(3100, event.GFkmtheta*TMath::RadToDeg());
      HF1(3110, GFkmCosTheta);
      HF1(3120, event.GFkmphi);
      HF1(3121, event.GFkmphi*TMath::RadToDeg());
      HF2(3500, event.GFkmtheta*TMath::RadToDeg(), event.GFkmmom);
      HF2(3510, event.GFkmphi*TMath::RadToDeg(), event.GFkmmom);
      HF2(3520, event.GFkmphi*TMath::RadToDeg(), event.GFkmtheta*TMath::RadToDeg());
      if(event.GFkmmom>0.01&&event.GFkmmass2>minM2Km&&event.GFkmmass2<maxM2Km){
	std::cout << " debug " << __FILE__ << " " << __LINE__ << " KmM2" << std::endl;
	if(missmass>minMM&&missmass<maxMM) HF1( 1, event.status++ ); // debug 8
	std::cout << " debug " << __FILE__ << " " << __LINE__ << " missmass" << std::endl;
	HF1(3951, event.MissMassCorrDE[0]);
	HF1(13951,-event.BETPC[0]);
	// hist
	HF1(4050, event.GFkmmom);
	HF1(4060, diff_kmmom.Mag());
	HF1(4070, miss_momTPC.Mag());
	double bek = -event.BETPC[0];
	for(int i=0; i<8; i++){
	  if( double(i)*0.050-0.100<bek && double(i+1)*0.050-0.100>bek ){
	    HF1(4062+i, diff_kmmom.Mag());
	  }
	}
	int binbek = GetBinIndexBEk(bek);	
	HF1(4062+binbek, diff_kmmom.Mag());	
	if(missmass>minMM&&missmass<maxMM){	
	  HF1(4051, event.GFkmmom);
	  HF1(4061, diff_kmmom.Mag());
	  HF1(4071, miss_momTPC.Mag());	      	    	  
	}
	HF1(4100, event.GFkmtheta*TMath::RadToDeg());
	HF1(4110, GFkmCosTheta);
	HF1(4120, event.GFkmphi);      
	HF1(4121, event.GFkmphi*TMath::RadToDeg());
	HF2(4500, event.GFkmtheta*TMath::RadToDeg(), event.GFkmmom);
	HF2(4510, event.GFkmphi*TMath::RadToDeg(), event.GFkmmom);
	HF2(4520, event.GFkmphi*TMath::RadToDeg(), event.GFkmtheta*TMath::RadToDeg()); 	
      }
    }
    double pmom = kp_momTPC.Mag();
    double ptheta = kp_momTPC.Theta(); //rad
    double pcost = kp_momTPC.CosTheta();
    double pphi = kp_momTPC.Phi(); //rad
    HF1(1050, pmom);
    HF1(1100, ptheta*TMath::RadToDeg());
    HF1(1110, pcost);
    HF1(1120, pphi);      
    HF1(1121, pphi*TMath::RadToDeg());
    HF2(1500, ptheta*TMath::RadToDeg(), pmom);
    HF2(1510, pphi*TMath::RadToDeg(), pmom);
    HF2(1520, pphi*TMath::RadToDeg(), ptheta*TMath::RadToDeg());
    if(event.GFkmmom>0.01&&event.GFkmmass2>minM2Km&&event.GFkmmass2<maxM2Km){
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;      
      HF1(2050, pmom);
      HF1(2100, ptheta*TMath::RadToDeg());
      HF1(2110, pcost);
      HF1(2120, pphi);      
      HF1(2121, pphi*TMath::RadToDeg());
      HF2(2500, ptheta*TMath::RadToDeg(), pmom);
      HF2(2510, pphi*TMath::RadToDeg(), pmom);
      HF2(2520, pphi*TMath::RadToDeg(), ptheta*TMath::RadToDeg());
    }
  }

  std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;    
  //HF1( 1, event.status++ ); // debug 9
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstClose( void )
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  const Int_t n = TFileCont.size();
  for( Int_t i=0; i<n; ++i ){
    if( TTreeReaderCont[i] ) delete TTreeReaderCont[i];
    if( TTreeCont[i] ) delete TTreeCont[i];
    if( TFileCont[i] ) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms( void )
{
  Int_t nbinpoq = 1000;
  Int_t minpoq = -1.5;
  Int_t maxpoq = 1.5;
  Int_t nbindedx = 1000;
  Int_t mindedx = 0;
  Int_t maxdedx = 350;
  Int_t nbinmass2 = 200;
  Int_t minmass2 = -1;
  Int_t maxmass2 =  1;
  
  HB1(1, "Status", 21, 0., 21. );
  HB2(10, "[g4] AnalysisKm <dE/dx>;p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",nbinpoq,minpoq,maxpoq,nbindedx,mindedx,maxdedx);    
  HB2(11, "[g4] GeantKm <dE/dx>;p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",nbinpoq,minpoq,maxpoq,nbindedx,mindedx,maxdedx);

  HB1(20, "[g4] Kaon M2; MassSquare [GeV]; counts ", nbinmass2, minmass2, maxmass2);
  HB1(21, "[g4] Kaon tracklen; tracklen [mm]; counts ", 1000, 0, 1000);
  HB1(22, "[g4] Kaon tof; tof [nsec]; counts ", 500, 0, 10);    
  HB1(23, "[g4] Kaon AnalysisdEdx; dEdx [arb.unit]; counts ", 1000, 0, 350);
  HB1(24, "[g4] Kaon GeantdEdx; dEdx [arb.unit]; counts ", 1000, 0, 350);

  HB1( 50, "[g4] Mom of G4GenScatK- ; #delta momentum [GeV/c]; coutns", 500, 1.5, 2.0);
  HB1(100, "[g4] Theta of G4GenScatK- [deg]; #theta [deg]; counts", 300, 0, 30); 
  HB1(110, "[g4] CosTheta of G4GenScatK- ; Cos(#theta); counts", 200, -1, 1);  
  HB1(120, "[g4] Phi of G4GenScatK- [rad]; #phi [rad]; counts", 500, -TMath::Pi(), TMath::Pi());
  HB1(121, "[g4] Phi of G4GenScatK- [deg]; #phi [deg]; counts", 3600, -180, 180);
  
  HB1(1050, "[g4] Mom of scatP ; momentum [GeV/c]; coutns", 500, 1.5, 2.0);
  HB1(1100, "[g4] Theta of scatP [deg]; #theta [deg]; counts", 300, 0, 30); 
  HB1(1110, "[g4] CosTheta of scatP ; Cos(#theta); counts", 200, -1, 1);
  HB1(1120, "[g4] Phi of scatP [rad]; #phi [rad]; counts", 500, -TMath::Pi(), TMath::Pi());
  HB1(1121, "[g4] Phi of scatP [deg]; #phi [deg]; counts", 3600, -180, 180);
  HB2(1500, "[g4] Mom vs Theta scatP; #theta [deg]; momentum [GeV/c]", 300, 0, 30, 500, 1.5, 2.0);
  HB2(1510, "[g4] Mom vs Phi scatP; #phi [deg]; momentum [GeV/c]", 3600, -180, 180, 500, 1.5, 2.0);
  HB2(1520, "[g4] Theta vs Phi scatP; #phi [deg]; #theta [deg]", 3600, -180, 180, 300, 0, 30);    

  HB1(2050, "[g4][M2] Mom of scatP ; momentum [GeV/c]; coutns", 500, 1.5, 2.0);
  HB1(2100, "[g4][M2] Theta of scatP [deg]; #theta [deg]; counts", 300, 0, 30); 
  HB1(2110, "[g4][M2] CosTheta of scatP ; Cos(#theta); counts", 200, -1, 1);
  HB1(2120, "[g4][M2] Phi of scatP [rad]; #phi [rad]; counts", 500, -TMath::Pi(), TMath::Pi());
  HB1(2121, "[g4][M2] Phi of scatP [deg]; #phi [deg]; counts", 3600, -180, 180);
  HB2(2500, "[g4][M2] Mom vs Theta scatP; #theta [deg]; momentum [GeV/c]", 300, 0, 30, 500, 1.5, 2.0);
  HB2(2510, "[g4][M2] Mom vs Phi scatP; #phi [deg]; momentum [GeV/c]", 3600, -180, 180, 500, 1.5, 2.0);
  HB2(2520, "[g4][M2] Theta vs Phi scatP; #phi [deg]; #theta [deg]", 3600, -180, 180, 300, 0, 30);  

  HB1(3050, "[g4] Mom of scatK- ; momentum [GeV/c]; coutns", 500, 0., 1.0);
  HB1(3100, "[g4] Theta of scatK- [deg]; #theta [deg]; counts", 1800, 0, 180); 
  HB1(3110, "[g4] CosTheta of scatK- ; Cos(#theta); counts", 200, -1, 1);
  HB1(3120, "[g4] Phi of scatK- [rad]; #phi [rad]; counts", 500, -TMath::Pi(), TMath::Pi());
  HB1(3121, "[g4] Phi of scatK- [deg]; #phi [deg]; counts", 3600, -180, 180);
  HB2(3500, "[g4] Mom vs Theta scatK-; #theta [deg]; momentum [GeV/c]", 1800, 0, 180, 500, 0, 1.0);
  HB2(3510, "[g4] Mom vs Phi scatK-; #phi [deg]; momentum [GeV/c]", 3600, -180, 180, 500, 0, 1.0);
  HB2(3520, "[g4] Theta vs Phi scatK-; #phi [deg]; #theta [deg]", 3600, -180, 180, 1800, 0, 180);
  HB1(3900, "[g4] MissingMass KP inclusive; MissingMass [GeV]; Counts", 280, 0., 1.4);
  HB1(3950, "[g4] MissingMass KP exclusive [dEdxPID]; MissingMass [GeV]; Counts", 280, 0., 1.4);
  HB1(3951, "[g4] MissingMass KP exclusive [M2PID]; MissingMass [GeV]; Counts", 280, 0., 1.4);

  HB1(4050, "[g4][M2] Mom of scatK- ; momentum [GeV/c]; coutns", 500, 0., 1.0);
  HB1(4051, "[g4][M2][missmass] Mom of scatK- ; momentum [GeV/c]; coutns", 500, 0., 1.0);
  HB1(4060, "[g4][M2] (MissMom - Mom) of scatK- ; momentum [GeV/c]; coutns", 500, 0., 1.0);
  HB1(4061, "[g4][M2][missmass] (MissMom - Mom) of scatK- ; momentum [GeV/c]; coutns", 500, 0., 1.0);
  for(int i=0; i<8; i++){
    HB1(4062+i, Form("[M2][bek] (MissMom - Mom) of scatK- (%.2f<-BE_K<%.2f); momentum [GeV/c]; coutns",double(i)*0.050-0.100, double(i+1)*0.050-0.100), 500, 0., 1.0);
  }
  
  HB1(4070, "[M2] MissMom of scatK- ; momentum [GeV/c]; coutns", 500, 0., 1.0);
  HB1(4071, "[M2][missmass] MissMom of scatK- ; momentum [GeV/c]; coutns", 500, 0., 1.0);    
  HB1(4100, "[g4][M2] Theta of scatK- [deg]; #theta [deg]; counts", 1800, 0, 180); 
  HB1(4110, "[g4][M2] CosTheta of scatK- ; Cos(#theta); counts", 200, -1, 1);
  HB1(4120, "[g4][M2] Phi of scatK- [rad]; #phi [rad]; counts", 500, -TMath::Pi(), TMath::Pi());
  HB1(4121, "[g4][M2] Phi of scatK- [deg]; #phi [deg]; counts", 3600, -180, 180);
  HB2(4500, "[g4][M2] Mom vs Theta scatK-; #theta [deg]; momentum [GeV/c]", 1800, 0, 180, 1000, 0, 2.0);
  HB2(4510, "[g4][M2] Mom vs Phi scatK-; #phi [deg]; momentum [GeV/c]", 3600, -180, 180, 1000, 0, 2.0);
  HB2(4520, "[g4][M2] Theta vs Phi scatK-; #phi [deg]; #theta [deg]", 3600, -180, 180, 1800, 0, 180);
  
  HB1(13900, "[g4] BE KP inclusive; -BE_{K} [GeV]; Counts", 120, -0.3, 0.3);
  HB1(13950, "[g4] BE KP exclusive [dEdxPID]; -BE_{K} [GeV]; Counts", 120, -0.3, 0.3);
  HB1(13951, "[g4] BE KP exclusive [M2PID]; -BE_{K} [GeV]; Counts", 120, -0.3, 0.3);

  for(int i=0; i<120; i++){
    HB1(14000+i, Form("[M2][bek] (MissMom - Mom) of scatK- (%.3f<-BE_K<%.3f); momentum [GeV/c]; coutns", double(i)*0.005-0.300, double(i+1)*0.005-0.300), 500, 0., 1.0);
  } 

  HBTree( "tpc", "tree of GenfitCarbon" );
  tree->Branch( "status", &event.status );
  tree->Branch( "evnum", &event.evnum );
  tree->Branch( "trigpat", &event.trigpat );
  tree->Branch( "trigflag", &event.trigflag );

  tree->Branch( "nhHtof", &event.nhHtof );
  tree->Branch( "HtofSeg", &event.HtofSeg );
  tree->Branch( "tHtof", &event.tHtof );
  tree->Branch( "dtHtof", &event.dtHtof );
  tree->Branch( "deHtof", &event.deHtof );
  tree->Branch( "posHtof", &event.posHtof );
  tree->Branch( "G4tidHtof", &event.G4tidHtof );

  tree->Branch("NumberOfTracks",&event.NumberOfTracks,"NumberOfTracks/I");
  tree->Branch("PIDOfTrack",&event.PIDOfTrack);
  tree->Branch("ParentIDOfTrack",&event.ParentIDOfTrack);
  tree->Branch("VertexOfTrack_x",&event.VertexOfTrack_x);
  tree->Branch("VertexOfTrack_y",&event.VertexOfTrack_y);
  tree->Branch("VertexOfTrack_z",&event.VertexOfTrack_z);
  tree->Branch("MomentumOfTrack",&event.MomentumOfTrack);
  tree->Branch("MomentumOfTrack_x",&event.MomentumOfTrack_x);
  tree->Branch("MomentumOfTrack_y",&event.MomentumOfTrack_y);
  tree->Branch("MomentumOfTrack_z",&event.MomentumOfTrack_z);

  tree->Branch("G4kmid",&event.G4kmid);
  tree->Branch("G4kmtid",&event.G4kmtid);
  tree->Branch("G4kmvtx_x",&event.G4kmvtx_x);
  tree->Branch("G4kmvtx_y",&event.G4kmvtx_y);
  tree->Branch("G4kmvtx_z",&event.G4kmvtx_z);
  tree->Branch("G4kmmom",&event.G4kmmom);
  tree->Branch("G4kmmom_x",&event.G4kmmom_x);
  tree->Branch("G4kmmom_y",&event.G4kmmom_y);
  tree->Branch("G4kmmom_z",&event.G4kmmom_z);

  tree->Branch("G4kpid",&event.G4kpid);
  tree->Branch("G4kptid",&event.G4kptid);
  tree->Branch("G4kpvtx_x",&event.G4kpvtx_x);
  tree->Branch("G4kpvtx_y",&event.G4kpvtx_y);
  tree->Branch("G4kpvtx_z",&event.G4kpvtx_z);
  tree->Branch("G4kpmom",&event.G4kpmom);
  tree->Branch("G4kpmom_x",&event.G4kpmom_x);
  tree->Branch("G4kpmom_y",&event.G4kpmom_y);
  tree->Branch("G4kpmom_z",&event.G4kpmom_z);

  tree->Branch("G4scatkmid",&event.G4scatkmid);
  tree->Branch("G4scatkmtid",&event.G4scatkmtid);
  tree->Branch("G4scatkmvtx_x",&event.G4scatkmvtx_x);
  tree->Branch("G4scatkmvtx_y",&event.G4scatkmvtx_y);
  tree->Branch("G4scatkmvtx_z",&event.G4scatkmvtx_z);
  tree->Branch("G4scatkmmom",&event.G4scatkmmom);
  tree->Branch("G4scatkmmom_x",&event.G4scatkmmom_x);
  tree->Branch("G4scatkmmom_y",&event.G4scatkmmom_y);
  tree->Branch("G4scatkmmom_z",&event.G4scatkmmom_z);  
  
  // tree->Branch("G4p2id",&event.G4p2id);  
  // tree->Branch("G4p2tid",&event.G4p2tid);
  // tree->Branch("G4p2nh",&event.G4p2nh);
  // tree->Branch("G4p2tnh",&event.G4p2tnh);
  // tree->Branch("p2tid",&event.p2tid);
  // tree->Branch("G4p2vtx_x",&event.G4p2vtx_x);
  // tree->Branch("G4p2vtx_y",&event.G4p2vtx_y);
  // tree->Branch("G4p2vtx_z",&event.G4p2vtx_z);
  // tree->Branch("G4p2mom",&event.G4p2mom);
  // tree->Branch("G4p2mom_x",&event.G4p2mom_x);
  // tree->Branch("G4p2mom_y",&event.G4p2mom_y);
  // tree->Branch("G4p2mom_z",&event.G4p2mom_z);

  tree->Branch( "nclTpc", &event.nclTpc );
  tree->Branch( "remain_nclTpc", &event.nclTpc );

  tree->Branch( "remain_cluster_x", &event.remain_cluster_x );
  tree->Branch( "remain_cluster_y", &event.remain_cluster_y );
  tree->Branch( "remain_cluster_z", &event.remain_cluster_z );
  tree->Branch( "remain_cluster_de", &event.remain_cluster_de );
  tree->Branch( "remain_cluster_size", &event.remain_cluster_size );
  tree->Branch( "remain_cluster_layer", &event.remain_cluster_layer );
#if 0
  tree->Branch( "remain_cluster_row_center", &event.remain_cluster_row_center );
  tree->Branch( "remain_cluster_mrow", &event.remain_cluster_mrow );
  tree->Branch( "remain_cluster_de_center", &event.remain_cluster_de_center );
  tree->Branch( "remain_cluster_x_center", &event.remain_cluster_x_center );
  tree->Branch( "remain_cluster_y_center", &event.remain_cluster_y_center );
  tree->Branch( "remain_cluster_z_center", &event.remain_cluster_z_center );
#endif
  tree->Branch( "remain_cluster_houghflag", &event.remain_cluster_houghflag );
  tree->Branch( "remain_cluster_G4tid" , &event.remain_cluster_G4tid );

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "isInTarget", &event.isInTarget );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "trackid", &event.trackid );
  tree->Branch( "isXi", &event.isXi );
  tree->Branch( "isBeam", &event.isBeam );
  tree->Branch( "isKurama", &event.isKurama );
  tree->Branch( "isK18", &event.isK18 );
  tree->Branch( "isAccidental", &event.isAccidental );
  tree->Branch( "isMultiloop", &event.isMultiloop );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "pid", &event.pid );
  tree->Branch( "purity", &event.purity );
  tree->Branch( "efficiency", &event.efficiency );
  tree->Branch( "G4tid", &event.G4tid );
  tree->Branch( "G4pid", &event.G4pid );
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "pval", &event.pval );
  tree->Branch( "helix_cx", &event.helix_cx );
  tree->Branch( "helix_cy", &event.helix_cy );
  tree->Branch( "helix_z0", &event.helix_z0 );
  tree->Branch( "helix_r", &event.helix_r );
  tree->Branch( "helix_dz", &event.helix_dz );
  tree->Branch( "dE", &event.dE );
  tree->Branch( "dEdx", &event.dEdx );
  tree->Branch( "mom0", &event.mom0 );
  tree->Branch( "path", &event.path );
  tree->Branch( "isElectron", &event.isElectron );
  tree->Branch( "nsigma_triton", &event.nsigma_triton );
  tree->Branch( "nsigma_deutron", &event.nsigma_deutron );
  tree->Branch( "nsigma_proton", &event.nsigma_proton );
  tree->Branch( "nsigma_kaon", &event.nsigma_kaon );
  tree->Branch( "nsigma_pion", &event.nsigma_pion );
  tree->Branch( "nsigma_electron", &event.nsigma_electron );
#if SaveRawData
  tree->Branch( "hitlayer", &event.hitlayer );
  tree->Branch( "hitpos_x", &event.hitpos_x );
  tree->Branch( "hitpos_y", &event.hitpos_y );
  tree->Branch( "hitpos_z", &event.hitpos_z );
  tree->Branch( "calpos_x", &event.calpos_x );
  tree->Branch( "calpos_y", &event.calpos_y );
  tree->Branch( "calpos_z", &event.calpos_z );
  tree->Branch( "residual", &event.residual );
  tree->Branch( "residual_x", &event.residual_x );
  tree->Branch( "residual_y", &event.residual_y );
  tree->Branch( "residual_z", &event.residual_z );
  tree->Branch( "resolution_x", &event.resolution_x);
  tree->Branch( "resolution_y", &event.resolution_y);
  tree->Branch( "resolution_z", &event.resolution_z);
  tree->Branch( "pathhit", &event.pathhit);
  tree->Branch( "alpha", &event.alpha);
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch( "track_cluster_row_center", &event.track_cluster_row_center);
#endif
  tree->Branch( "helix_t", &event.helix_t );
  tree->Branch( "chargeIndistinguishable", &event.chargeIndistinguishable );
  tree->Branch( "chisqr_inverted", &event.chisqr_inverted );
  tree->Branch( "pval_inverted", &event.pval_inverted );
  tree->Branch( "helix_cx_inverted", &event.helix_cx_inverted );
  tree->Branch( "helix_cy_inverted", &event.helix_cy_inverted );
  tree->Branch( "helix_z0_inverted", &event.helix_z0_inverted );
  tree->Branch( "helix_r_inverted", &event.helix_r_inverted );
  tree->Branch( "helix_dz_inverted", &event.helix_dz_inverted );
  tree->Branch( "mom0_inverted", &event.mom0_inverted );
  tree->Branch( "pid_inverted", &event.pid_inverted );

  tree->Branch( "nvtxTpc", &event.nvtxTpc );
  tree->Branch( "vtx_x", &event.vtx_x );
  tree->Branch( "vtx_y", &event.vtx_y );
  tree->Branch( "vtx_z", &event.vtx_z );
  tree->Branch( "vtx_dist", &event.vtx_dist );
  tree->Branch( "vtx_angle", &event.vtx_angle );
  tree->Branch( "vtxid", &event.vtxid );
  tree->Branch( "vtxmom_theta", &event.vtxmom_theta );
  tree->Branch( "vtxpos_x", &event.vtxpos_x );
  tree->Branch( "vtxpos_y", &event.vtxpos_y );
  tree->Branch( "vtxpos_z", &event.vtxpos_z );
  tree->Branch( "vtxmom_x", &event.vtxmom_x );
  tree->Branch( "vtxmom_y", &event.vtxmom_y );
  tree->Branch( "vtxmom_z", &event.vtxmom_z );

  tree->Branch("GFcharge", &event.GFcharge);
  tree->Branch("GFchisqr", &event.GFchisqr);
  tree->Branch("GFtof", &event.GFtof);
  //tree->Branch("GFtracklen", &event.GFtracklen);
  tree->Branch("GFpval", &event.GFpval);
  tree->Branch("GFfitstatus", &event.GFfitstatus);
  tree->Branch("GFpdgcode", &event.GFpdgcode);
  tree->Branch("GFnhtrack", &event.GFnhtrack);
  tree->Branch("GFlayer", &event.GFlayer);
  tree->Branch("GFpos_x", &event.GFpos_x);
  tree->Branch("GFpos_y", &event.GFpos_y);
  tree->Branch("GFpos_z", &event.GFpos_z);
  tree->Branch("GFmom", &event.GFmom);
  tree->Branch("GFmom_x", &event.GFmom_x);
  tree->Branch("GFmom_y", &event.GFmom_y);
  tree->Branch("GFmom_z", &event.GFmom_z);
  tree->Branch("GFresidual_x", &event.GFresidual_x);
  tree->Branch("GFresidual_y", &event.GFresidual_y);
  tree->Branch("GFresidual_z", &event.GFresidual_z);
  tree->Branch("GFresidual_p", &event.GFresidual_p);
  tree->Branch("GFresidual_px", &event.GFresidual_px);
  tree->Branch("GFresidual_py", &event.GFresidual_py);
  tree->Branch("GFresidual_pz", &event.GFresidual_pz);

  tree->Branch("GFntTpc_target", &event.GFntTpc_inside);
  tree->Branch("GFprodvtx_x", &event.GFprodvtx_x);
  tree->Branch("GFprodvtx_y", &event.GFprodvtx_y);
  tree->Branch("GFprodvtx_z", &event.GFprodvtx_z);

  //extrapolation
  tree->Branch("GFinside", &event.GFinside);
  tree->Branch("GFtracklen", &event.GFtracklen);
  tree->Branch("GFtrack2vtxdist", &event.GFtrack2vtxdist);
  tree->Branch("GFcalctof", &event.GFcalctof);
  tree->Branch("GFsegHtof", &event.GFsegHtof);
  tree->Branch("GFtofHtof", &event.GFtofHtof);
  tree->Branch("GFtdiffHtof", &event.GFtdiffHtof);
  tree->Branch("GFposHtof", &event.GFposHtof);
  tree->Branch("GFposx", &event.GFposx);
  tree->Branch("GFposy", &event.GFposy);
  tree->Branch("GFposz", &event.GFposz);
  tree->Branch("GFinvbeta", &event.GFinvbeta);
  tree->Branch("GFm2", &event.GFm2);
  tree->Branch("nsigma_tritonHtof", &event.nsigma_tritonHtof);
  tree->Branch("nsigma_deutronHtof", &event.nsigma_deutronHtof);
  tree->Branch("nsigma_protonHtof", &event.nsigma_protonHtof);
  tree->Branch("nsigma_kaonHtof", &event.nsigma_kaonHtof);
  tree->Branch("nsigma_pionHtof", &event.nsigma_pionHtof);
  tree->Branch("nsigma_electronHtof", &event.nsigma_electronHtof);

  tree->Branch( "nvtxTpcClustered", &event.nvtxTpcClustered );
  tree->Branch( "clusteredVtx_x", &event.clusteredVtx_x );
  tree->Branch( "clusteredVtx_y", &event.clusteredVtx_y );
  tree->Branch( "clusteredVtx_z", &event.clusteredVtx_z );
  tree->Branch( "clusteredVtxid", &event.clusteredVtxid );

  tree->Branch("GFKmTrackId", &event.GFkmid);
  tree->Branch("GFKmMom", &event.GFkmmom);
  tree->Branch("GFKmMom_x", &event.GFkmmom_x);
  tree->Branch("GFKmMom_y", &event.GFkmmom_y);
  tree->Branch("GFKmMom_z", &event.GFkmmom_z);
  tree->Branch("GFKmTheta", &event.GFkmtheta);
  tree->Branch("GFKmPhi", &event.GFkmphi);    
  tree->Branch("GFKmTargetCloseDist", &event.GFkmtarget_dist);
  tree->Branch("GFKmTarget_x", &event.GFkmtargetvtx_x);
  tree->Branch("GFKmTarget_y", &event.GFkmtargetvtx_y);
  tree->Branch("GFKmTarget_z", &event.GFkmtargetvtx_z);
  tree->Branch("GFKmTargetCenter_x", &event.GFkmtargetcenter_x);
  tree->Branch("GFKmTargetCenter_y", &event.GFkmtargetcenter_y);
  tree->Branch("GFKmTargetCenter_z", &event.GFkmtargetcenter_z);
  tree->Branch("GFKmTargetCenterCloseDist", &event.GFkmtargetcenter_dist);
  tree->Branch("GFKmHtofId", &event.GFkmhtofid);
  tree->Branch("GFKmHtofSeg", &event.GFkmhtofseg);
  tree->Branch("GFKmHtofPos", &event.GFkmposHtof);  
  tree->Branch("GFKmMassSquare", &event.GFkmmass2);
  tree->Branch("GFKmInvBeta", &event.GFkminvbeta);         
  tree->Branch("GFKmTrackLen", &event.GFkmtracklen);
  tree->Branch("GFKmTof", &event.GFkmtof);  

  tree->Branch( "ntK18", &event.ntK18);
  tree->Branch( "pK18", &event.pK18);
  tree->Branch( "chisqrK18", &event.chisqrK18);
  tree->Branch( "xtgtK18", &event.xtgtK18);
  tree->Branch( "ytgtK18", &event.ytgtK18);
  tree->Branch( "utgtK18", &event.utgtK18);
  tree->Branch( "vtgtK18", &event.vtgtK18);

  tree->Branch( "isgoodTPCK18", &event.isgoodTPCK18);
  tree->Branch( "chisqrTPCK18", &event.chisqrTPCK18);
  tree->Branch( "qTPCK18", &event.qTPCK18);
  tree->Branch( "pTPCK18", &event.pTPCK18);
  tree->Branch( "xtgtTPCK18", &event.xtgtTPCK18);
  tree->Branch( "ytgtTPCK18", &event.ytgtTPCK18);
  tree->Branch( "utgtTPCK18", &event.utgtTPCK18);
  tree->Branch( "vtgtTPCK18", &event.vtgtTPCK18);
  tree->Branch( "thetaTPCK18", &event.thetaTPCK18);
  tree->Branch( "lhtofTPCK18", &event.lhtofTPCK18);
  tree->Branch( "xhtofTPCK18", &event.xhtofTPCK18);
  tree->Branch( "yhtofTPCK18", &event.yhtofTPCK18);
  tree->Branch( "lvpTPCK18", &event.lvpTPCK18);
  tree->Branch( "xvpTPCK18", &event.xvpTPCK18);
  tree->Branch( "yvpTPCK18", &event.yvpTPCK18);

  tree->Branch( "ntKurama", &event.ntKurama);
  tree->Branch( "pKurama", &event.pKurama);
  tree->Branch( "qKurama", &event.qKurama);
  tree->Branch( "chisqrKurama", &event.chisqrKurama);
  tree->Branch( "m2Kurama", &event.m2);
  tree->Branch( "xtgtKurama", &event.xtgtKurama);
  tree->Branch( "ytgtKurama", &event.ytgtKurama);
  tree->Branch( "utgtKurama", &event.utgtKurama);
  tree->Branch( "vtgtKurama", &event.vtgtKurama);
  tree->Branch( "thetaKurama", &event.thetaKurama);
  tree->Branch( "pathwcKurama", &event.pathwcKurama);
  tree->Branch( "xin", &event.xin);
  tree->Branch( "yin", &event.yin);
  tree->Branch( "zin", &event.zin);
  tree->Branch( "pxin", &event.pxin);
  tree->Branch( "pyin", &event.pyin);
  tree->Branch( "pzin", &event.pzin);
  tree->Branch( "isgoodTPCKurama", &event.isgoodTPCKurama);
  tree->Branch( "kflagTPCKurama", &event.kflagTPCKurama);
  tree->Branch( "chisqrTPCKurama", &event.chisqrTPCKurama);
  tree->Branch( "pTPCKurama", &event.pTPCKurama);
  tree->Branch( "qTPCKurama", &event.qTPCKurama);
  tree->Branch( "m2TPCKurama", &event.m2TPCKurama);
  tree->Branch( "xtgtTPCKurama", &event.xtgtTPCKurama);
  tree->Branch( "ytgtTPCKurama", &event.ytgtTPCKurama);
  tree->Branch( "utgtTPCKurama", &event.utgtTPCKurama);
  tree->Branch( "vtgtTPCKurama", &event.vtgtTPCKurama);
  tree->Branch( "thetaTPCKurama", &event.thetaTPCKurama);
  tree->Branch( "pathTPCKurama", &event.pathTPCKurama);
  tree->Branch( "lhtofTPCKurama", &event.lhtofTPCKurama);
  tree->Branch( "xhtofTPCKurama", &event.xhtofTPCKurama);
  tree->Branch( "yhtofTPCKurama", &event.yhtofTPCKurama);
  tree->Branch( "lvpTPCKurama", &event.lvpTPCKurama);
  tree->Branch( "xvpTPCKurama", &event.xvpTPCKurama);
  tree->Branch( "yvpTPCKurama", &event.yvpTPCKurama);

  tree->Branch( "isgoodTPC", &event.isgoodTPC);
  tree->Branch( "insideTPC", &event.insideTPC);
  tree->Branch( "vtxTPC", &event.vtxTPC);
  tree->Branch( "vtyTPC", &event.vtyTPC);
  tree->Branch( "vtzTPC", &event.vtzTPC);
  tree->Branch( "closeDistTPC", &event.closeDistTPC);
  tree->Branch( "MissMassTPC", &event.MissMassTPC);
  tree->Branch( "MissMassCorrTPC", &event.MissMassCorrTPC);
  tree->Branch( "MissMassCorrDETPC", &event.MissMassCorrDETPC);
  tree->Branch( "pOrgTPC", &event.pOrgTPC);
  tree->Branch( "pCalcTPC", &event.pCalcTPC);
  tree->Branch( "pCorrTPC", &event.pCorrTPC);
  tree->Branch( "pCorrDETPC", &event.pCorrDETPC);
  tree->Branch( "xbTPC", &event.xbTPC);
  tree->Branch( "ybTPC", &event.ybTPC);
  tree->Branch( "ubTPC", &event.ubTPC);
  tree->Branch( "vbTPC", &event.vbTPC);
  tree->Branch( "xsTPC", &event.xsTPC);
  tree->Branch( "ysTPC", &event.ysTPC);
  tree->Branch( "usTPC", &event.usTPC);
  tree->Branch( "vsTPC", &event.vsTPC);

  tree->Branch("nKK", &event.nKK);
  tree->Branch("Kflag", &event.Kflag);
  tree->Branch("MissMass", &event.MissMass);
  tree->Branch("MissMassCorr", &event.MissMassCorr);
  tree->Branch("MissMassCorrDE", &event.MissMassCorrDE);
  tree->Branch("MissMassNucl", &event.MissMassNucl);
  tree->Branch("MissMassNuclCorr", &event.MissMassNuclCorr);
  tree->Branch("MissMassNuclCorrDE", &event.MissMassNuclCorrDE);  
  tree->Branch("vtx", &event.vtx);
  tree->Branch("vty", &event.vty);
  tree->Branch("vtz", &event.vtz);
  tree->Branch("pOrg", &event.pOrg);
  tree->Branch("pCalc", &event.pCalc);
  tree->Branch("pCorr", &event.pCorr);
  tree->Branch("pCorrDE", &event.pCorrDE);
  tree->Branch("xb", &event.xb);
  tree->Branch("yb", &event.yb);
  tree->Branch("ub", &event.ub);
  tree->Branch("vb", &event.vb);
  tree->Branch("xs", &event.xs);
  tree->Branch("ys", &event.ys);
  tree->Branch("us", &event.us);
  tree->Branch("vs", &event.vs);

  tree->Branch("KmMom_x", &event.km_mom_x);
  tree->Branch("KmMom_y", &event.km_mom_y);
  tree->Branch("KmMom_z", &event.km_mom_z);
  tree->Branch("KpMom_x", &event.kp_mom_x);
  tree->Branch("KpMom_y", &event.kp_mom_y);
  tree->Branch("KpMom_z", &event.kp_mom_z);

  tree->Branch("BE", &event.BE);
  tree->Branch("BETPC", &event.BETPC);
  tree->Branch("BE_LL", &event.BE_LL);
  tree->Branch("BETPC_LL", &event.BETPC_LL);

  tree->Branch("LPiflag", &event.lpiflag);
  tree->Branch("PiPiflag", &event.pipiflag);
  tree->Branch("Pimflag", &event.pimflag);
  tree->Branch("Emptyflag", &event.emptyflag);

  tree->Branch("DecaysTrackId", &event.decays_id);
  tree->Branch("DecaysMom", &event.decays_mom);
  tree->Branch("DecaysMom_x", &event.decays_mom_x);
  tree->Branch("DecaysMom_y", &event.decays_mom_y);
  tree->Branch("DecaysMom_z", &event.decays_mom_z);
  tree->Branch("DecaysMomCM", &event.decays_CMmom);
  tree->Branch("DecaysMomCM_x", &event.decays_CMmom_x);
  tree->Branch("DecaysMomCM_y", &event.decays_CMmom_y);
  tree->Branch("DecaysMomCM_z", &event.decays_CMmom_z);

  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("nhittpc",&src.nhittpc);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("ititpc",src.ititpc);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("xtpc",src.xtpc);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("ytpc",src.ytpc);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("ztpc",src.ztpc);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("pxtpc",src.pxtpc);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("pytpc",src.pytpc);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("pztpc",src.pztpc);

  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("NumberOfTracks",&src.NumberOfTracks);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("PIDOfTrack",src.PIDOfTrack);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("ParentIDOfTrack",src.ParentIDOfTrack);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("VertexOfTrack_x",src.VertexOfTrack_x);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("VertexOfTrack_y",src.VertexOfTrack_y);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("VertexOfTrack_z",src.VertexOfTrack_z);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("MomentumOfTrack",src.MomentumOfTrack);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("MomentumOfTrack_x",src.MomentumOfTrack_x);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("MomentumOfTrack_y",src.MomentumOfTrack_y);
  TTreeCont[kHelixTrackingGeant4]->SetBranchAddress("MomentumOfTrack_z",src.MomentumOfTrack_z);

  TTreeReaderCont[kHelixTrackingGeant4] = new TTreeReader( "tpc", TFileCont[kHelixTrackingGeant4] );
  const auto& reader = TTreeReaderCont[kHelixTrackingGeant4];
  src.evnum = new TTreeReaderValue<Int_t>( *reader, "evnum" );

  src.nhHtof = new TTreeReaderValue<Int_t>( *reader, "nhHtof" );
  src.HtofSeg = new TTreeReaderValue<std::vector<Double_t>>( *reader, "HtofSeg" );
  src.tHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "tHtof" );
  //src.dtHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "dtHtof" );
  src.deHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "deHtof" );
  src.posHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "posHtof" );
  src.G4tidHtof = new TTreeReaderValue<std::vector<Int_t>>( *reader, "G4tidHtof" );

  src.nclTpc = new TTreeReaderValue<Int_t>( *reader, "nclTpc" );
  src.remain_nclTpc = new TTreeReaderValue<Int_t>( *reader, "remain_nclTpc" );
#if SaveRawData
  src.cluster_x = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_x" );
  src.cluster_y = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_y" );
  src.cluster_z = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_z" );
  src.cluster_de = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_de" );
  src.cluster_size = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_size" );
  src.cluster_layer = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_layer" );
  src.cluster_mrow = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_mrow" );
  src.cluster_row_center = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_row_center" );
  src.cluster_houghflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_houghflag" );
#endif

  src.ntTpc = new TTreeReaderValue<Int_t>( *reader, "ntTpc" );
  src.nhtrack = new TTreeReaderValue<std::vector<Int_t>>( *reader, "nhtrack" );
  src.trackid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trackid" );
  src.isXi = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isXi" );
  src.isBeam = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isBeam" );
  src.isKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isKurama" );
  src.isK18 = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isK18" );
  src.isAccidental = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isAccidental" );
  src.isMultiloop = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isMultiloop" );
  src.charge = new TTreeReaderValue<std::vector<Int_t>>( *reader, "charge" );
  src.pid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "pid" );
  src.purity = new TTreeReaderValue<std::vector<Double_t>>( *reader, "purity" );
  src.efficiency = new TTreeReaderValue<std::vector<Double_t>>( *reader, "efficiency" );
  src.G4tid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "G4tid" );
  src.chisqr = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqr" );
  src.pval = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pval" );
  src.helix_cx = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cx" );
  src.helix_cy = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cy" );
  src.helix_z0 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_z0" );
  src.helix_r = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_r" );
  src.helix_dz = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_dz" );
  src.dE = new TTreeReaderValue<std::vector<Double_t>>( *reader, "dE" );
  src.dEdx = new TTreeReaderValue<std::vector<Double_t>>( *reader, "dEdx" );
  src.mom0 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "mom0" );
  src.path = new TTreeReaderValue<std::vector<Double_t>>( *reader, "path" );
  src.isElectron = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isElectron" );
  src.nsigma_triton = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_triton" );
  src.nsigma_deutron = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_deutron" );
  src.nsigma_proton = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_proton" );
  src.nsigma_kaon = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_kaon" );
  src.nsigma_pion = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_pion" );
  src.nsigma_electron = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_electron" );
  src.hitlayer = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "hitlayer" );
  src.hitpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "hitpos_x" );
  src.hitpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "hitpos_y" );
  src.hitpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "hitpos_z" );
  src.calpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "calpos_x" );
  src.calpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "calpos_y" );
  src.calpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "calpos_z" );
  src.mom_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "mom_x" );
  src.mom_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "mom_y" );
  src.mom_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "mom_z" );
  src.residual = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual" );
  src.residual_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual_x" );
  src.residual_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual_y" );
  src.residual_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual_z" );
  src.resolution_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "resolution_x" );
  src.resolution_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "resolution_y" );
  src.resolution_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "resolution_z" );
  src.helix_t = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "helix_t" );
  src.pathhit = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "pathhit" );
  src.alpha = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "alpha" );
  src.track_cluster_de = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_de" );
  src.track_cluster_size = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_size" );
  src.track_cluster_mrow = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_mrow" );

  src.chargeIndistinguishable = new TTreeReaderValue<std::vector<Int_t>>( *reader, "chargeIndistinguishable" );
  src.pid_inverted = new TTreeReaderValue<std::vector<Int_t>>( *reader, "pid_inverted" );
  src.chisqr_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqr_inverted" );
  src.pval_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pval_inverted" );
  src.helix_cx_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cx_inverted" );
  src.helix_cy_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cy_inverted" );
  src.helix_z0_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_z0_inverted" );
  src.helix_r_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_r_inverted" );
  src.helix_dz_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_dz_inverted" );
  src.mom0_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "mom0_inverted" );

  src.nvtxTpc = new TTreeReaderValue<Int_t>(*reader,"nvtxTpc");
  src.vtx_x = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx_x" );
  src.vtx_y = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx_y" );
  src.vtx_z = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx_z" );
  src.vtx_dist = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx_dist" );
  src.vtx_angle = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx_angle" );
  src.vtxid = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxid" );
  src.vtxmom_theta = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxmom_theta" );
  src.vtxpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxpos_x" );
  src.vtxpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxpos_y" );
  src.vtxpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxpos_z" );
  src.vtxmom_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxmom_x" );
  src.vtxmom_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxmom_y" );
  src.vtxmom_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxmom_z" );

  src.isLambda = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isLambda" );
  src.ncombiLambda = new TTreeReaderValue<std::vector<Int_t>>( *reader, "ncombiLambda" );
  src.distLambda = new TTreeReaderValue<std::vector<Double_t>>( *reader, "distLambda" );
  src.angleLambda = new TTreeReaderValue<std::vector<Double_t>>( *reader, "angleLambda" );
  src.bestmassLambda = new TTreeReaderValue<std::vector<Double_t>>( *reader, "bestmassLambda" );
  src.massLambda = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "massLambda" );
  src.vtxLambda_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxLambda_x" );
  src.vtxLambda_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxLambda_y" );
  src.vtxLambda_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxLambda_z" );
  src.momLambda = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "momLambda" );
  src.momLambda_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "momLambda_x" );
  src.momLambda_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "momLambda_y" );
  src.momLambda_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "momLambda_z" );
  src.decaysidLambda = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "decaysidLambda" );
  src.decaysmomLambda = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "decaysmomLambda" );
  src.decaysmomLambda_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "decaysmomLambda_x" );
  src.decaysmomLambda_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "decaysmomLambda_y" );
  src.decaysmomLambda_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "decaysmomLambda_z" );

  src.nvtxTpcClustered = new TTreeReaderValue<Int_t>(*reader,"nvtxTpcClustered");
  src.clusteredVtx_x = new TTreeReaderValue<std::vector<Double_t>>( *reader, "clusteredVtx_x" );
  src.clusteredVtx_y = new TTreeReaderValue<std::vector<Double_t>>( *reader, "clusteredVtx_y" );
  src.clusteredVtx_z = new TTreeReaderValue<std::vector<Double_t>>( *reader, "clusteredVtx_z" );
  src.clusteredVtxid = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "clusteredVtxid" );

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")   &&
      InitializeParameter<TPCParamMan>("TPCPRM") &&
      InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
      InitializeParameter<UserParamMan>("USER") &&
      InitializeParameter<FieldMan>("FLDMAP", "HSFLDMAP") &&
      InitializeParameter<HodoPHCMan>("HDPHC") );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  return true;
}
