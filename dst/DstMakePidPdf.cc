// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <TLorentzVector.h>
#include <memory> 

#include <filesystem_util.hh>
#include <UnpackerManager.hh>
#include "TF1.h"

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
#include "PidPdfMan.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCCluster.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCVertex.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

//#include "HypTPCFitter.hh"
//#include "HypTPCTask.hh"

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
auto&       gPidPdf = PidPdfMan::GetInstance();  
const auto& zHSCenter = gGeom.LocalZ("HS");
const Double_t truncatedMean = 0.8; //80%
const Double_t dist_cut = 5;
const Double_t chisqrcut = 100;
const Double_t decutppi = 2.0;
const Double_t gfposycut = 100.0;  

//For GenFit Setting
const Bool_t Const_field = false; //Must be false for linear tracking
const Int_t verbosity = 0;//3;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");

const int kPmin = pidlikeli::kPmin;
const int kPmax = pidlikeli::kPmax;
const int kMomstep = pidlikeli::kDP;
const int kNmom = pidlikeli::kNmom;
const int kNpid = pidlikeli::kNpid; 
const int kNchg = pidlikeli::kNchg;
const int kNtype = pidlikeli::kNtype;
const int kNbe = pidlikeli::kNbe;
const auto& plist = pidlikeli::plist;
const auto& clist = pidlikeli::clist;
const auto& type = pidlikeli::type;  
  
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kGenfitE42, kOutFile, nArgc
};
std::vector<TString> ArgName =
  { "[Process]", "[ConfFile]", "[GenfitE42]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________________
struct Event
{

  Int_t status;
  Int_t runnum;
  Int_t evnum;
  void clear( void )
  {
    runnum = 0;
    evnum = 0;
    status = 0;
  }
};

//_____________________________________________________________________________
struct Src
{
  std::unique_ptr<TTreeReaderValue<Int_t>> runnum;
  std::unique_ptr<TTreeReaderValue<Int_t>> evnum;
};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
  enum eDetHid {
    PadHid    = 100000,
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
  // HypTPCFitter* fitter = new HypTPCFitter(tpcGeo.Data(),Const_field);
  // //Initiallize the genfit track container
  // HypTPCTask& GFtracks = HypTPCTask::GetInstance();
  // GFtracks.SetVerbosity(verbosity);
  // std::cout<<"GenFit verbosity = "<<"-1: Silent, 0: Minimum, 1: Errors only, 2: Errors and Warnings, 3: Verbose mode, long term debugging(default)"<<std::endl;
  // std::cout<<"Current verbosity = "<<GFtracks.GetVerbosity()<<std::endl;

#if 0
  GFtracks.DebugMode();
#endif
  Int_t ievent = skip;
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
    gCounter.check();
    InitializeEvent();
    if( DstRead( ievent ) ) tree->Fill();
  }
  std::cout << "#D Event Number: " << std::setw(6)
            << ievent << std::endl;
  
  gPidPdf.WriteToRootfile(TFileCont[kOutFile]);
  //  gPidPdf.DrawResultsToPdf("anafig/tpc/pidlikeli/PidPdf.pdf");  
  DstClose();

  //  delete fitter;
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
  int open_file = 0;
  int open_tree = 0;
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
dst::DstRead( int ievent )
{
  return true;
  Double_t vtx_scan_range = gUser.GetParameter("VertexScanRange");
  //if( ievent%1000==0 ){
  if( ievent%100==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }
  /*

  GetEntry(ievent);

event.runnum = *(*src.runnum);
  event.evnum = *(*src.evnum);
  event.trigpat = *(*src.trigpat);
  event.trigflag = *(*src.trigflag);

  event.nclTpc = *(*src.nclTpc);
  event.cluster_x = *(*src.cluster_x);
  event.cluster_y = *(*src.cluster_y);
  event.cluster_z = *(*src.cluster_z);
  event.cluster_de = *(*src.cluster_de);
  event.cluster_size = *(*src.cluster_size);
  event.cluster_layer = *(*src.cluster_layer);
  event.cluster_mrow = *(*src.cluster_mrow);
  event.cluster_de_center = *(*src.cluster_de_center);
  event.cluster_x_center = *(*src.cluster_x_center);
  event.cluster_y_center = *(*src.cluster_y_center);
  event.cluster_z_center = *(*src.cluster_z_center);
  event.cluster_row_center = *(*src.cluster_row_center);
  event.cluster_houghflag = *(*src.cluster_houghflag);


  int ntTpc = *(*src.ntTpc);
  event.ntTpc = ntTpc;
  event.nhtrack = *(*src.nhtrack);
  event.isBeam = *(*src.isBeam);
  event.isKurama = *(*src.isKurama);
  event.isK18 = *(*src.isK18);
  event.isAccidental = *(*src.isAccidental);
  event.chisqr = *(*src.chisqr);
  event.helix_cx = *(*src.helix_cx);
  event.helix_cy = *(*src.helix_cy);
  event.helix_z0 = *(*src.helix_z0);
  event.helix_r = *(*src.helix_r);
  event.helix_dz = *(*src.helix_dz);

  event.dE = *(*src.dE);
  event.dEdx = *(*src.dEdx);
  event.mom0 = *(*src.mom0);
  event.charge = *(*src.charge);
  event.path = *(*src.path);
  event.pid = *(*src.pid);

  event.hitlayer = *(*src.hitlayer);
  event.hitpos_x = *(*src.hitpos_x);
  event.hitpos_y = *(*src.hitpos_y);
  event.hitpos_z = *(*src.hitpos_z);
  event.calpos_x = *(*src.calpos_x);
  event.calpos_y = *(*src.calpos_y);
  event.calpos_z = *(*src.calpos_z);
  event.mom_x = *(*src.mom_x);
  event.mom_y = *(*src.mom_y);
  event.mom_z = *(*src.mom_z);
  event.residual = *(*src.residual);
  event.residual_x = *(*src.residual_x);
  event.residual_y = *(*src.residual_y);
  event.residual_z = *(*src.residual_z);
  event.resolution_x = *(*src.resolution_x);
  event.resolution_y = *(*src.resolution_y);
  event.resolution_z = *(*src.resolution_z);
  event.helix_t = *(*src.helix_t);
  event.alpha = *(*src.alpha);
  event.pathhit = *(*src.pathhit);
  event.track_cluster_de = *(*src.track_cluster_de);
  event.track_cluster_size = *(*src.track_cluster_size);
  event.track_cluster_mrow = *(*src.track_cluster_mrow);
  event.track_cluster_de_center = *(*src.track_cluster_de_center);
  event.track_cluster_x_center = *(*src.track_cluster_x_center);
  event.track_cluster_y_center = *(*src.track_cluster_y_center);
  event.track_cluster_z_center = *(*src.track_cluster_z_center);
  event.track_cluster_row_center = *(*src.track_cluster_row_center);

  event.isgoodTPCKurama = *(*src.isgoodTPCKurama);
  event.insideTPC = *(*src.insideTPC);
  event.pTPCKurama = *(*src.pTPCKurama);
  event.qTPCKurama = *(*src.qTPCKurama);
  event.m2TPCKurama = *(*src.m2TPCKurama);
  event.xsTPC = *(*src.xsTPC);
  event.ysTPC = *(*src.ysTPC);
  event.usTPC = *(*src.usTPC);
  event.vsTPC = *(*src.vsTPC);

  event.pK18 = *(*src.pK18);
  event.xbTPC = *(*src.xbTPC);
  event.ybTPC = *(*src.ybTPC);
  event.ubTPC = *(*src.ubTPC);
  event.vbTPC = *(*src.vbTPC);


  event.isElectron.resize(ntTpc);
  event.nsigma_triton.resize(ntTpc);
  event.nsigma_deutron.resize(ntTpc);
  event.nsigma_proton.resize(ntTpc);
  event.nsigma_kaon.resize(ntTpc);
  event.nsigma_pion.resize(ntTpc);
  event.nsigma_electron.resize(ntTpc);
  for(int it=0; it<ntTpc; ++it){
    event.isElectron[it] = Kinematics::HypTPCdEdxElectron(event.dEdx[it], event.mom0[it]);
    event.nsigma_triton[it] = Kinematics::HypTPCdEdxNsigmaTriton(event.dEdx[it], event.mom0[it]);
    event.nsigma_deutron[it] = Kinematics::HypTPCdEdxNsigmaDeutron(event.dEdx[it], event.mom0[it]);
    event.nsigma_proton[it] = Kinematics::HypTPCdEdxNsigmaProton(event.dEdx[it], event.mom0[it]);
    event.nsigma_kaon[it]  = Kinematics::HypTPCdEdxNsigmaKaon(event.dEdx[it], event.mom0[it]);
    event.nsigma_pion[it] = Kinematics::HypTPCdEdxNsigmaPion(event.dEdx[it], event.mom0[it]);
    event.nsigma_electron[it] = Kinematics::HypTPCdEdxNsigmaElectron(event.dEdx[it], event.mom0[it]);
  }

  event.nvtxTpc = *(*src.nvtxTpc);
  event.vtx_x = *(*src.vtx_x);
  event.vtx_y = *(*src.vtx_y);
  event.vtx_z = *(*src.vtx_z);
  event.vtx_dist = *(*src.vtx_dist);
  event.vtx_angle = *(*src.vtx_angle);
  event.vtxid = *(*src.vtxid);
  event.vtxmom_theta = *(*src.vtxmom_theta);
  event.vtxpos_x = *(*src.vtxpos_x);
  event.vtxpos_y = *(*src.vtxpos_y);
  event.vtxpos_z = *(*src.vtxpos_z);
  event.vtxmom_x = *(*src.vtxmom_x);
  event.vtxmom_y = *(*src.vtxmom_y);
  event.vtxmom_z = *(*src.vtxmom_z);

  event.HtofSeg = *(*src.HtofSeg);  
  event.tHtof = *(*src.tHtof);
  event.dtHtof = *(*src.dtHtof);
  event.deHtof = *(*src.deHtof);
  event.posHtof = *(*src.posHtof);
  
  event.GFntTpc = *(*src.GFntTpc);
  event.GFfitstatus = *(*src.GFfitstatus);
  event.GFpdgcode = *(*src.GFpdgcode);
  event.GFnhtrack = *(*src.GFnhtrack);
  event.GFcharge = *(*src.GFcharge);
  event.GFchisqr = *(*src.GFchisqr);
  event.GFtof = *(*src.GFtof);
  event.GFpval = *(*src.GFpval);
  event.GFlayer = *(*src.GFlayer);
  event.GFpos_x = *(*src.GFpos_x);
  event.GFpos_y = *(*src.GFpos_y);
  event.GFpos_z = *(*src.GFpos_z);
  event.GFmom   = *(*src.GFmom);	
  event.GFmom_x = *(*src.GFmom_x);
  event.GFmom_y = *(*src.GFmom_y);
  event.GFmom_z = *(*src.GFmom_z);
  event.GFresidual_x = *(*src.GFresidual_x);
  event.GFresidual_y = *(*src.GFresidual_y);
  event.GFresidual_z = *(*src.GFresidual_z);
  event.GFresidual_p = *(*src.GFresidual_p);
  event.GFresidual_px = *(*src.GFresidual_px);
  event.GFresidual_py = *(*src.GFresidual_py);
  event.GFresidual_pz = *(*src.GFresidual_pz);
  
  event.GFinside = *(*src.GFinside);
  event.GFextrapolationHtof = *(*src.GFextrapolationHtof);

  event.GFntTpc_inside = *(*src.GFntTpc_inside);
  event.GFprodvtx_x = *(*src.GFprodvtx_x);
  event.GFprodvtx_y = *(*src.GFprodvtx_y);
  event.GFprodvtx_z = *(*src.GFprodvtx_z);

  event.GFtracklen          = *(*src.GFtracklen);	    
  event.GFtrack2vtxdist	    = *(*src.GFtrack2vtxdist);    
  event.GFcalctof	    = *(*src.GFcalctof);	    
  event.GFsegHtof	    = *(*src.GFsegHtof);	    
  event.GFtofHtof	    = *(*src.GFtofHtof);	    
  event.GFtdiffHtof	    = *(*src.GFtdiffHtof);	    
  event.GFposHtof	    = *(*src.GFposHtof);	    
  event.GFposx	    	    = *(*src.GFposx);		    
  event.GFposy		    = *(*src.GFposy);		    
  event.GFposz		    = *(*src.GFposz);		    
  event.GFinvbeta	    = *(*src.GFinvbeta);	    
  event.GFm2	    	    = *(*src.GFm2);		    
  event.nsigma_tritonHtof   = *(*src.nsigma_tritonHtof);  
  event.nsigma_deutronHtof  = *(*src.nsigma_deutronHtof); 
  event.nsigma_protonHtof   = *(*src.nsigma_protonHtof);  
  event.nsigma_kaonHtof	    = *(*src.nsigma_kaonHtof);    
  event.nsigma_pionHtof     = *(*src.nsigma_pionHtof);    
  event.nsigma_electronHtof = *(*src.nsigma_electronHtof);

  HF1( 10, ntTpc );
  Int_t GFntTpc = event.GFntTpc;
  if( event.ntTpc == 0 ) return true;

  // for debug
  if(1) return true;
  
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    return true;    
  }
    
  HF1( 1, event.status++ );
  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCTracks(*(*src.ntTpc), *(*src.isK18), *(*src.isKurama),
			 *(*src.charge), *(*src.nhtrack), *(*src.helix_cx),
			 *(*src.helix_cy), *(*src.helix_z0), *(*src.helix_r),
			 *(*src.helix_dz), *(*src.hitlayer), *(*src.track_cluster_mrow),
			 *(*src.helix_t), *(*src.track_cluster_de), *(*src.resolution_x),
			 *(*src.resolution_y), *(*src.resolution_z), *(*src.hitpos_x),
			 *(*src.hitpos_y), *(*src.hitpos_z));

  HypTPCTask& GFTrackCont = HypTPCTask::GetInstance();
  
  std::cout << "debug " << __FILE__ << " " << __LINE__ << " GFntTpc: " << GFntTpc <<  std::endl;    

  for( Int_t igf=0; igf<GFntTpc; ++igf ){
    if(!event.GFfitstatus[igf]) continue;
    if(!event.GFextrapolationHtof[igf]) continue;

    double invb = event.GFinvbeta[igf];
    int charge = event.charge[igf];
    int chg = pidlikeli::ChgToBin(charge);
    double p = event.GFmom[igf][0];
    double dedx = event.dEdx[igf];
    int type = 0;
    int be = 0;
    
    std::cout << "debug " << __FILE__ << " " << __LINE__
	      << "(it,pid,chg,be,mom)=("<< 0 <<",PiKPDE," << chg << ",0," << p
	      << "), (invb,dedx)=(" << invb << "," << dedx <<")"<< std::endl;    
    
    auto pidL = gPidLike.CalcLikeliResultPid(type,chg,be,p,invb,dedx);
    if(pidL.empty()) continue;
    for(int ith=0; ith<kNpid; ith++){
      std::cout << " No." << ith
		<< " Prob: " << pidL.nthPidProb[ith]	
		<< " PId: " << pidL.nthPid[ith] 
		<< " Good: " << pidL.nthPidGood[ith] 
		<< std::endl;
    }
  }
  
  HF1( genfitHid, GFntTpc);
  HF1( 1, event.status++ );

  GFTrackCont.Clear();
  */
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
  /*
  Int_t nbinpoq = 1000;
  Int_t minpoq = -1.5;
  Int_t maxpoq = 1.5;
  Int_t nbindedx = 1000;
  Int_t mindedx = 0;
  Int_t maxdedx = 350;
  Int_t nbinmass2 = 500;
  Int_t minmass2 = -5;
  Int_t maxmass2 = 5;
  Int_t nbininvbeta = 100;
  Int_t mininvbeta = 0;
  Int_t maxinvbeta = 5;
  
  HB1( 1, "Status", 21, 0., 21. );
  HB1( 2, "Genfit Status", 20, 0., 20. );
  HB1( 3, "Genfit Fit Status", 2, 0., 2. );
  HB1( 10, "NTrack TPC", 20, 0., 20. );
  */
  HBTree( "tpc", "tree of DstMakePdf" );
  tree->Branch( "status", &event.status );
  tree->Branch( "runnum", &event.runnum );
  tree->Branch( "evnum", &event.evnum );

  TTreeReaderCont[kGenfitE42] = new TTreeReader( "tpc", TFileCont[kGenfitE42] );
  const auto& reader = TTreeReaderCont[kGenfitE42];
  src.runnum = std::make_unique<TTreeReaderValue<Int_t>>( *reader, "runnum" );
  src.evnum = std::make_unique<TTreeReaderValue<Int_t>>( *reader, "evnum" );
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
      InitializeParameter<PidPdfMan>("PIDPDF") &&
      InitializeParameter<FieldMan>("FLDMAP", "HSFLDMAP") &&
      InitializeParameter<HodoPHCMan>("HDPHC") );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  return true;
}
