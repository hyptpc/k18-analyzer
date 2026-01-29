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
#include "DatabasePDG.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
//#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DstHelper.hh"
#include "HistTools.hh"
#include "HodoPHCMan.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCAnalyzer.hh"
#include "TPCCluster.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#define RawHit 0
#define TrackSearch 1
#define TrackCluster 1
#define TruncatedMean 0

namespace
{
  using namespace root;
  using namespace dst;
  using hddaq::unpacker::GUnpacker;
  const auto& gUnpacker = GUnpacker::get_instance();
  auto&       gConf = ConfMan::GetInstance();
  const auto& gGeom = DCGeomMan::GetInstance();
  const auto& gUser = UserParamMan::GetInstance();
  const auto& gPHC  = HodoPHCMan::GetInstance();
  const auto& gCounter = debug::ObjectCounter::GetInstance();
  const double truncatedMean = 0.8; //80%
}

namespace dst
{
  enum kArgc
    {
      kProcess, kConfFile,
      kTpcHit,  kOutFile, nArgc
    };
  std::vector<TString> ArgName =
    { "[Process]", "[ConfFile]", "[TPCHit]",  "[OutFile]" };
  std::vector<TString> TreeName = { "", "", "tpc", "" };
  std::vector<TFile*> TFileCont;
  std::vector<TTree*> TTreeCont;
  std::vector<TTreeReader*> TTreeReaderCont;
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
  std::vector<Int_t> cluster_houghflag;

  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> isBeam; // isBeam: 1 = Beam, 0 = Scat
  std::vector<Double_t> chisqr;
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx; //reference dedx

  std::vector<Double_t> dEdx_0;
  std::vector<Double_t> dEdx_10;
  std::vector<Double_t> dEdx_20;
  std::vector<Double_t> dEdx_30;
  std::vector<Double_t> dEdx_40;
  std::vector<Double_t> dEdx_50;
  std::vector<Double_t> dEdx_60;
  std::vector<Double_t> dEdx_cor_0;
  std::vector<Double_t> dEdx_cor_10;
  std::vector<Double_t> dEdx_cor_20;
  std::vector<Double_t> dEdx_cor_30;
  std::vector<Double_t> dEdx_cor_40;
  std::vector<Double_t> dEdx_cor_50;
  std::vector<Double_t> dEdx_cor_60;

  std::vector<Double_t> dz_factor;
  std::vector<Double_t> mom0_x;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_y;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_z;//Helix momentum at Y = 0
  std::vector<Double_t> mom0;//Helix momentum at Y = 0

  std::vector<Int_t> charge;//Helix charge
  std::vector<Double_t> path;//Helix path

  std::vector<std::vector<Double_t>> hitlayer;
  std::vector<std::vector<Double_t>> hitpos_x;
  std::vector<std::vector<Double_t>> hitpos_y;
  std::vector<std::vector<Double_t>> hitpos_z;
  std::vector<std::vector<Double_t>> calpos_x;
  std::vector<std::vector<Double_t>> calpos_y;
  std::vector<std::vector<Double_t>> calpos_z;
  std::vector<std::vector<Double_t>> residual;
  std::vector<std::vector<Double_t>> residual_x;
  std::vector<std::vector<Double_t>> residual_y;
  std::vector<std::vector<Double_t>> residual_z;
  std::vector<std::vector<Double_t>> helix_t;
  std::vector<std::vector<Double_t>> pathhit;
  std::vector<std::vector<Double_t>> pathhit_cor;
  std::vector<std::vector<Double_t>> theta_diff;
  std::vector<std::vector<Double_t>> track_cluster_de;
  std::vector<std::vector<Double_t>> track_cluster_size;
  std::vector<std::vector<Double_t>> track_cluster_mrow;
  std::vector<std::vector<Double_t>> track_cluster_de_center;
  std::vector<std::vector<Double_t>> track_cluster_x_center;
  std::vector<std::vector<Double_t>> track_cluster_y_center;
  std::vector<std::vector<Double_t>> track_cluster_z_center;
  std::vector<std::vector<Double_t>> track_cluster_row_center;

  void clearBasicInfo() {
    runnum   = 0;
    evnum    = 0;
    status   = 0;
    beamflag = beam::kUnknown;
    dst::clear_all(trigpat, trigflag, clkTpc);
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
                       cluster_z_center, cluster_row_center, cluster_houghflag);
  }

  void clearHelixTracks() {
    ntTpc = 0;
    dst::clear_all(
      nhtrack, isBeam, chisqr, helix_cx, helix_cy, helix_z0, helix_r, helix_dz, dE, dEdx
      
#if TruncatedMean
      , dEdx_0, dEdx_10, dEdx_20, dEdx_30, dEdx_40, dEdx_50, dEdx_60, dEdx_cor_0, dEdx_cor_10, dEdx_cor_20, dEdx_cor_30, dEdx_cor_40, dEdx_cor_50, dEdx_cor_60
#endif
      , dz_factor, mom0_x, mom0_y, mom0_z, mom0, charge, path,
      
      hitlayer, hitpos_x, hitpos_y, hitpos_z, calpos_x, calpos_y, calpos_z, residual, residual_x, residual_y, residual_z, helix_t,

      pathhit, pathhit_cor, theta_diff, track_cluster_de, track_cluster_size, track_cluster_mrow, track_cluster_de_center, track_cluster_x_center, track_cluster_y_center, track_cluster_z_center, track_cluster_row_center
    );
    
  }
  
  void clear()
  {
    clearBasicInfo();
    clearRawHits();
    clearClusters();
    clearHelixTracks();
  }

  void resizeTracks(Int_t nTracks) {
    dst::resize_all(nTracks,
      nhtrack, isBeam, chisqr, helix_cx, helix_cy, helix_z0, helix_r, helix_dz, dE, dEdx

#if TruncatedMean
      , dEdx_0, dEdx_10, dEdx_20, dEdx_30, dEdx_40, dEdx_50, dEdx_60, dEdx_cor_0, dEdx_cor_10, dEdx_cor_20, dEdx_cor_30, dEdx_cor_40, dEdx_cor_50, dEdx_cor_60
#endif
      , dz_factor, mom0_x, mom0_y, mom0_z, mom0, charge, path,
      
      hitlayer, hitpos_x, hitpos_y, hitpos_z, calpos_x, calpos_y, calpos_z, residual, residual_x, residual_y, residual_z, helix_t,

      pathhit, pathhit_cor, theta_diff, track_cluster_de, track_cluster_size, track_cluster_mrow, track_cluster_de_center, track_cluster_x_center, track_cluster_y_center, track_cluster_z_center, track_cluster_row_center
    );
  }
  void resizeTrackHits(Int_t it, Int_t nh) {
    dst::resize_all(nh,
      hitlayer[it], hitpos_x[it], hitpos_y[it], hitpos_z[it], calpos_x[it], calpos_y[it], calpos_z[it], residual[it], residual_x[it], residual_y[it], residual_z[it], helix_t[it],

      pathhit[it], pathhit_cor[it], theta_diff[it], track_cluster_de[it], track_cluster_size[it], track_cluster_mrow[it], track_cluster_de_center[it], track_cluster_x_center[it], track_cluster_y_center[it], track_cluster_z_center[it], track_cluster_row_center
    );
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
  TTreeReaderValue<std::vector<Double_t>>* pedTpc;    // pedestal
  TTreeReaderValue<std::vector<Double_t>>* rmsTpc;    // rms
  TTreeReaderValue<std::vector<Double_t>>* deTpc;     // dE
  TTreeReaderValue<std::vector<Double_t>>* cdeTpc;    // cdE
  TTreeReaderValue<std::vector<Double_t>>* tTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* ctTpc;     // time
  TTreeReaderValue<std::vector<Double_t>>* chisqrTpc; // chi^2 of signal fitting
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;    // clock time
};

namespace root
{
  Event  event;
  Src    src;
  TTree *tree;

  Double_t
  TranseverseDistance(Double_t x_center, Double_t z_center, Double_t x, Double_t z)
  {
    Double_t dummy = TMath::Sqrt((x-x_center)*(x-x_center) + (z-z_center)*(z-z_center));
    Double_t dist;
    if(x_center-x<0) dist=-1.*dummy;
    else dist=dummy;
    return dist;
  }
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
  if(!gConf.InitializeHistograms())
    return EXIT_FAILURE;
  if( !gConf.InitializeUnpacker() )
    return EXIT_FAILURE;
  if(!dst::SetupReader())
    return EXIT_FAILURE;
  
  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries( TTreeCont );
  if (max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();

  Int_t ievent = skip;
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
    gCounter.check();
    InitializeEvent();
    if( DstRead( ievent ) ) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6)
            << ievent << std::endl;

  DstClose();

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
  if( ievent%100==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  event.runnum   = **src.runnum;
  event.evnum    = **src.evnum;
  event.trigpat  = **src.trigpat;
  event.trigflag = **src.trigflag;
  event.beamflag = **src.beamflag;
  event.clkTpc   = **src.clkTpc;
  HF1("Status", event.status++);

  if( **src.nhTpc == 0 )
    return true;

  HF1("Status", event.status++);
  
  if(event.clkTpc.size() != 8){
    spdlog::warn("something is wrong: event.clkTpc.size() != 8");
    return true;
  }

  HF1("Status", event.status++);
  
  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, **src.clkTpc);
  HF1("Status", event.status++);

  const Double_t VertexScanRange = gUser.GetParameter("VertexScanRange"); //mm

  HF1("Status", event.status++);

#if RawHit
  Int_t nhTpc = 0;
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    auto hc = TPCAna.GetTPCHC( layer );
    for( const auto& hit : hc ){
      if( !hit || !hit->IsGood() )
        continue;
      const auto& pos = hit->GetPosition();
      Double_t x   = pos.X();
      Double_t y   = pos.Y();
      Double_t z   = pos.Z();
      Double_t de  = hit->GetCDe();
      Int_t    pad = hit->GetPad();
      Int_t    row = hit->GetRow();
      event.raw_hitpos_x.push_back(x);
      event.raw_hitpos_y.push_back(y);
      event.raw_hitpos_z.push_back(z);
      event.raw_de.push_back(de);
      event.raw_padid.push_back(pad);
      event.raw_layer.push_back(layer);
      event.raw_row.push_back(row);
      ++nhTpc;
    }
  }
  event.nhTpc = nh_Tpc;
  HF1("Status", event.status++);
#endif

#if RawCluster
  Int_t nclTpc = 0;
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    auto hc = TPCAna.GetTPCClCont( layer );
    for( const auto& cl : hc ){
      if( !cl || !cl->IsGood() )
        continue;
      Double_t x  = cl->GetX();
      Double_t y  = cl->GetY();
      Double_t z  = cl->GetZ();
      Double_t de = cl->GetDe();
      Int_t cl_size = cl->GetClusterSize();
      Double_t mrow = cl->MeanRow();
      Int_t houghFlag = cl->GetHoughFlag();
      TPCHit* centerHit = cl->GetCenterHit();
      const TVector3& centerPos = centerHit->GetPosition();
      Double_t centerDe = centerHit->GetCDe();
      Int_t centerRow = centerHit->GetRow();

      event.cluster_x.push_back(x);
      event.cluster_y.push_back(y);
      event.cluster_z.push_back(z);
      event.cluster_de.push_back(de);
      event.cluster_size.push_back(cl_size);
      event.cluster_layer.push_back(layer);
      event.cluster_mrow.push_back(mrow);
      event.cluster_houghflag.push_back(houghFlag);
      event.cluster_de_center.push_back(centerDe);
      event.cluster_x_center.push_back(centerPos.X());
      event.cluster_y_center.push_back(centerPos.Y());
      event.cluster_z_center.push_back(centerPos.Z());
      event.cluster_row_center.push_back(centerRow);

      ++nclTpc;
    }
  }
  event.nclTpc = nclTpc;
  HF1("Status", event.status++);
#endif
  
#if TrackSearch
  TPCAna.TrackSearchTPCHelix();
#endif

  Int_t ntTpc = TPCAna.GetNTracksTPCHelix();
  event.ntTpc = ntTpc;
  HF1("NTracks_TPC", ntTpc);
  if( event.ntTpc == 0 )
    return true;
  
  HF1("Status", event.status++);
  event.resizeTracks(ntTpc);

  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( it );
    if( !tp ) continue;
    Int_t nh = tp->GetNHit();
    Double_t chisqr = tp->GetChiSquare();
    Double_t helix_cx = tp->Getcx(), helix_cy = tp->Getcy();
    Double_t helix_z0 = tp->Getz0(), helix_r = tp->Getr();
    Double_t helix_dz = tp->Getdz();
    TVector3 Mom0 = tp->GetMom0();
    Int_t isbeam = tp->GetIsBeam();

    event.nhtrack[it] = nh;
    event.isBeam[it] = isbeam;
    event.chisqr[it] = chisqr;
    event.helix_cx[it] = helix_cx;
    event.helix_cy[it] = helix_cy;
    event.helix_z0[it] = helix_z0;
    event.helix_r[it] = helix_r ;
    event.helix_dz[it] = helix_dz;
    event.mom0_x[it] = Mom0.x();
    event.mom0_y[it] = Mom0.y();
    event.mom0_z[it] = Mom0.z();
    event.mom0[it] = Mom0.Mag();
    event.resizeTrackHits(it, nh);

    HF1("NHits_Track_TPC", nh);
    HF1("Chisqr_TPC", chisqr);
    
    Double_t par1[5]={helix_cx, helix_cy, helix_z0,
         helix_r, helix_dz};

    Double_t scantheta1 = VertexScanRange/par1[3]; //mm -> rad.
    Double_t range_theta1[2] = {tp -> GetMint() - scantheta1,
				tp -> GetMaxt() + scantheta1};

    for( Int_t it2=0; it2<ntTpc; ++it2 ){
      if(it2==it) continue;
      TPCLocalTrackHelix *tp2 = TPCAna.GetTrackTPCHelix( it2 );
      if( !tp2 ) continue;
      Double_t helix_cx2 = tp2->Getcx(), helix_cy2 = tp2->Getcy();
      Double_t helix_z02 = tp2->Getz0(), helix_r2 = tp2->Getr();
      Double_t helix_dz2 = tp2->Getdz();

      Double_t par2[5]={helix_cx2, helix_cy2, helix_z02,
			helix_r2, helix_dz2};

      Double_t scantheta2 = VertexScanRange/par2[3]; //mm -> rad.
      Double_t range_theta2[2] = {tp -> GetMint() - scantheta2,
				  tp -> GetMaxt() + scantheta2};

      double closeDistTpc, t1, t2;
      TVector3 vert
	= Kinematics::VertexPointHelix(par1, par2,
				       range_theta1[0], range_theta1[1],
				       range_theta2[0], range_theta2[1],
				       t1, t2, closeDistTpc);

    }

    Double_t min_t = 10000.; Double_t max_t = -10000.;
    Double_t min_layer_t = 0., max_layer_t = 0.;
    Int_t min_layer = 33, max_layer = -1;
    Double_t de=0.;
    std::vector<Double_t> dEdx_vect; std::vector<Double_t> dEdx_cor_vect;
    for( int ih=0; ih<nh; ++ih ){
      TPCLTrackHit *hit = tp->GetHit( ih );
      if( !hit ) continue;
      HF1("HoughDist", hit->GetHoughDist());
      HF1("HoughDistY", hit->GetHoughDistY());

      Int_t layer = hit->GetLayer();
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPosHelix();
      const TVector3& resi_vect = hit->GetResidualVect();

      HF1("LayerId_TPC", layer);

      TPCHit *clhit = hit -> GetHit();
      TPCCluster *cl = clhit -> GetParentCluster();
      Int_t clsize = cl->GetClusterSize();
      Double_t clde = cl->GetDe();
      Double_t mrow = cl->MeanRow(); // same
      TPCHit* centerHit = cl->GetCenterHit();
      const TVector3& centerPos = centerHit->GetPosition();
      Double_t centerDe = centerHit->GetCDe();
      Int_t centerRow = centerHit->GetRow();

#if TrackCluster
      if(ntTpc==1&&nh>=20){
	HF1("Cluster_size", clsize);
	HF1(Form("Cluster_size_layer%d",layer), clsize);
	HF1("Cluster_dE", clde);
	HF1(Form("Cluster_dE_layer%d",layer), clde);
	const TPCHitContainer& hc = cl -> GetHitContainer();
	for(const auto& hits : hc){
	  if(!hits || !hits->IsGood() || !hit->IsGoodForTracking()) continue;
	  const TVector3& pos = hits->GetPosition();
	  Double_t de = hits->GetCDe();
	  Double_t transDist = TranseverseDistance(hitpos.x(), hitpos.z(), pos.x(), pos.z());
	  Double_t ratio = de/clde;
	  HF2("Transverse_diffusion", transDist, ratio);
	  HF2(Form("Transverse_diffusion_layer%d",layer), transDist, ratio);
	}
      }
#endif
      event.track_cluster_de[it][ih] = clde;
      event.track_cluster_size[it][ih] = clsize;
      event.track_cluster_mrow[it][ih] = mrow;
      event.track_cluster_de_center[it][ih] = centerDe;
      event.track_cluster_x_center[it][ih] = centerPos.X();
      event.track_cluster_y_center[it][ih] = centerPos.Y();
      event.track_cluster_z_center[it][ih] = centerPos.Z();
      event.track_cluster_row_center[it][ih] = centerRow;

      Double_t padTheta = tpc::GetTheta(layer, mrow)*acos(-1)/180.;
      Double_t t_cal = hit->GetTheta();
      Double_t thetaDiff = t_cal - padTheta;
      event.theta_diff[it][ih] = thetaDiff;
      Double_t pathHit = tpc::padParameter[layer][5];
      event.pathhit[it][ih] = pathHit;

      //Approximation of pathHit correction
      Double_t pathHit_cor = (pathHit/fabs(cos(thetaDiff)))*sqrt(1.+(pow(helix_dz,2)));
      event.pathhit_cor[it][ih] = pathHit_cor;

      //Charge
      if(min_t>t_cal) min_t = t_cal;
      if(max_t<t_cal) max_t = t_cal;
      if(layer<min_layer){
	min_layer = layer;
	min_layer_t = t_cal;
      }
      if(layer>max_layer){
	max_layer = layer;
	max_layer_t = t_cal;
      }

      de += clde;
      Double_t dEdx = clde/pathHit;
      Double_t dEdx_cor = clde/pathHit_cor;
      dEdx_vect.push_back(dEdx);
      dEdx_cor_vect.push_back(dEdx_cor);

      Double_t residual = hit->GetResidual();
      event.hitlayer[it][ih] = (double)layer;
      event.hitpos_x[it][ih] = hitpos.x();
      event.hitpos_y[it][ih] = hitpos.y();
      event.hitpos_z[it][ih] = hitpos.z();
      event.calpos_x[it][ih] = calpos.x();
      event.calpos_y[it][ih] = calpos.y();
      event.calpos_z[it][ih] = calpos.z();
      event.residual[it][ih] = residual;
      event.residual_x[it][ih] = resi_vect.x();
      event.residual_y[it][ih] = resi_vect.y();
      event.residual_z[it][ih] = resi_vect.z();
      event.helix_t[it][ih] = t_cal;
    }
    if(min_layer_t<max_layer_t) event.charge[it] = 1;
    else event.charge[it] = -1;

    Double_t pathlen = (max_t - min_t)*sqrt(helix_r*helix_r*(1.+helix_dz*helix_dz));
    event.path[it] = pathlen;
    event.dE[it] = de;

    std::sort(dEdx_vect.begin(), dEdx_vect.end());
    std::sort(dEdx_cor_vect.begin(), dEdx_cor_vect.end());
    int n_truncated = (int)(dEdx_vect.size()*truncatedMean);
    for( int ih=0; ih<dEdx_vect.size(); ++ih ){
      if(ih<n_truncated) event.dEdx[it]+=dEdx_cor_vect[ih];
    }
    event.dEdx[it]/=(double)n_truncated;

#if TruncatedMean
    event.dEdx_0[it] = TMath::Mean(dEdx_vect.size(),dEdx_vect.data());
    event.dEdx_cor_0[it] = TMath::Mean(dEdx_cor_vect.size(),dEdx_cor_vect.data());
    event.dEdx_10[it]=0.;
    event.dEdx_20[it]=0.;
    event.dEdx_30[it]=0.;
    event.dEdx_40[it]=0.;
    event.dEdx_50[it]=0.;
    event.dEdx_60[it]=0.;
    event.dEdx_cor_10[it]=0.;
    event.dEdx_cor_20[it]=0.;
    event.dEdx_cor_30[it]=0.;
    event.dEdx_cor_40[it]=0.;
    event.dEdx_cor_50[it]=0.;
    event.dEdx_cor_60[it]=0.;

    int n_10 = (int)(dEdx_vect.size()*0.9);
    int n_20 = (int)(dEdx_vect.size()*0.8);
    int n_30 = (int)(dEdx_vect.size()*0.7);
    int n_40 = (int)(dEdx_vect.size()*0.6);
    int n_50 = (int)(dEdx_vect.size()*0.5);
    int n_60 = (int)(dEdx_vect.size()*0.4);
    for( int ih=0; ih<dEdx_vect.size(); ++ih ){
      if(ih<n_10){
	event.dEdx_10[it]+=dEdx_vect[ih];
	event.dEdx_cor_10[it]+=dEdx_cor_vect[ih];
      }
      if(ih<n_20){
	event.dEdx_20[it]+=dEdx_vect[ih];
	event.dEdx_cor_20[it] +=dEdx_cor_vect[ih];
      }
      if(ih<n_30){
	event.dEdx_30[it]+=dEdx_vect[ih];
	event.dEdx_cor_30[it]+=dEdx_cor_vect[ih];
      }
      if(ih<n_40){
	event.dEdx_40[it]+=dEdx_vect[ih];
	event.dEdx_cor_40[it]+=dEdx_cor_vect[ih];
      }
      if(ih<n_50){
	event.dEdx_50[it]+=dEdx_vect[ih];
	event.dEdx_cor_50[it]+=dEdx_cor_vect[ih];
      }
      if(ih<n_60){
	event.dEdx_60[it]+=dEdx_vect[ih];
	event.dEdx_cor_60[it]+=dEdx_cor_vect[ih];
      }
    }
    event.dEdx_10[it]/=(double)n_10;
    event.dEdx_20[it]/=(double)n_20;
    event.dEdx_30[it]/=(double)n_30;
    event.dEdx_40[it]/=(double)n_40;
    event.dEdx_50[it]/=(double)n_50;
    event.dEdx_60[it]/=(double)n_60;
    event.dEdx_cor_10[it]/=(double)n_10;
    event.dEdx_cor_20[it]/=(double)n_20;
    event.dEdx_cor_30[it]/=(double)n_30;
    event.dEdx_cor_40[it]/=(double)n_40;
    event.dEdx_cor_50[it]/=(double)n_50;
    event.dEdx_cor_60[it]/=(double)n_60;
#endif
    event.dz_factor[it] = sqrt(1.+(pow(helix_dz,2)));
    HF1("mom0", event.mom0[it]);

  }

  HF1("Status", event.status++);
  
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
dst::SetupReader()
{
  if (!dst::SetupReader(kTpcHit, "kTpcHit")) return false;

  const auto& reader = TTreeReaderCont[kTpcHit];

  dst::SetBranch(reader, "run_number",   src.runnum);
  dst::SetBranch(reader, "event_number", src.evnum);
  dst::SetBranch(reader, "trig_pat",     src.trigpat);
  dst::SetBranch(reader, "trig_flag",    src.trigflag);
  dst::SetBranch(reader, "beam_flag",    src.beamflag);
  dst::SetBranch(reader, "npadTpc",      src.npadTpc);
  dst::SetBranch(reader, "nhTpc",        src.nhTpc);
  dst::SetBranch(reader, "layerTpc",     src.layerTpc);
  dst::SetBranch(reader, "rowTpc",       src.rowTpc);
  dst::SetBranch(reader, "padTpc",       src.padTpc);
  dst::SetBranch(reader, "pedTpc",       src.pedTpc);
  dst::SetBranch(reader, "rmsTpc",       src.rmsTpc);
  dst::SetBranch(reader, "deTpc",        src.deTpc);
  dst::SetBranch(reader, "tTpc",         src.tTpc);
  dst::SetBranch(reader, "chisqrTpc",    src.chisqrTpc);
  dst::SetBranch(reader, "clkTpc",       src.clkTpc);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms( void )
{
  hist::BuildStatus();
  hist::BuildTPCHelixTracking();
  

  
  tree = new TTree("tpc", "tree of DstTPCTracking");
  tree->Branch("status", &event.status);
  tree->Branch("run_number", &event.runnum);
  tree->Branch("event_number", &event.evnum);
  tree->Branch("trig_pat", &event.trigpat);
  tree->Branch("trig_flag", &event.trigflag);
  tree->Branch("beam_flag", &event.beamflag);
  tree->Branch("clkTpc", &event.clkTpc);

#if RawHit
  tree->Branch( "nhTpc", &event.nhTpc );
  tree->Branch( "raw_hitpos_x", &event.raw_hitpos_x );
  tree->Branch( "raw_hitpos_y", &event.raw_hitpos_y );
  tree->Branch( "raw_hitpos_z", &event.raw_hitpos_z );
  tree->Branch( "raw_de", &event.raw_de );
  tree->Branch( "raw_padid", &event.raw_padid );
  tree->Branch( "raw_layer", &event.raw_layer );
  tree->Branch( "raw_row", &event.raw_row );
#endif

  tree->Branch( "nclTpc", &event.nclTpc );
  tree->Branch( "cluster_x", &event.cluster_x );
  tree->Branch( "cluster_y", &event.cluster_y );
  tree->Branch( "cluster_z", &event.cluster_z );
  tree->Branch( "cluster_de", &event.cluster_de );
  tree->Branch( "cluster_size", &event.cluster_size );
  tree->Branch( "cluster_layer", &event.cluster_layer );
  tree->Branch( "cluster_row_center", &event.cluster_row_center );
  tree->Branch( "cluster_mrow", &event.cluster_mrow );
  tree->Branch( "cluster_de_center", &event.cluster_de_center );
  tree->Branch( "cluster_x_center", &event.cluster_x_center );
  tree->Branch( "cluster_y_center", &event.cluster_y_center );
  tree->Branch( "cluster_z_center", &event.cluster_z_center );

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "isBeam", &event.isBeam );
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "helix_cx", &event.helix_cx );
  tree->Branch( "helix_cy", &event.helix_cy );
  tree->Branch( "helix_z0", &event.helix_z0 );
  tree->Branch( "helix_r", &event.helix_r );
  tree->Branch( "helix_dz", &event.helix_dz );
  tree->Branch( "mom0_x", &event.mom0_x );
  tree->Branch( "mom0_y", &event.mom0_y );
  tree->Branch( "mom0_z", &event.mom0_z );
  tree->Branch( "mom0", &event.mom0 );
  tree->Branch( "dE", &event.dE );
  tree->Branch( "dEdx", &event.dEdx );

#if TruncatedMean
  tree->Branch( "dEdx_0", &event.dEdx_0 );
  tree->Branch( "dEdx_10", &event.dEdx_10 );
  tree->Branch( "dEdx_20", &event.dEdx_20 );
  tree->Branch( "dEdx_30", &event.dEdx_30 );
  tree->Branch( "dEdx_40", &event.dEdx_40 );
  tree->Branch( "dEdx_50", &event.dEdx_50 );
  tree->Branch( "dEdx_60", &event.dEdx_60 );
  tree->Branch( "dEdx_cor_0", &event.dEdx_cor_0 );
  tree->Branch( "dEdx_cor_10", &event.dEdx_cor_10 );
  tree->Branch( "dEdx_cor_20", &event.dEdx_cor_20 );
  tree->Branch( "dEdx_cor_30", &event.dEdx_cor_30 );
  tree->Branch( "dEdx_cor_40", &event.dEdx_cor_40 );
  tree->Branch( "dEdx_cor_50", &event.dEdx_cor_50 );
  tree->Branch( "dEdx_cor_60", &event.dEdx_cor_60 );
#endif
  tree->Branch( "dz_factor", &event.dz_factor );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "path", &event.path );

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
  tree->Branch( "helix_t", &event.helix_t );
  tree->Branch( "theta_diff", &event.theta_diff);
  tree->Branch( "pathhit", &event.pathhit);
  tree->Branch( "pathhit_cor", &event.pathhit_cor);
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch( "track_cluster_row_center", &event.track_cluster_row_center);

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
      InitializeParameter<HodoPHCMan>("HDPHC") );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  return true;
}

