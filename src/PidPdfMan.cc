// -*- C++ -*-
#include <TMath.h>
#include <iostream>
#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TKey.h>
#include <TROOT.h>
#include <TClass.h>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <fstream>
#include <array>
#include <std_ostream.hh>
#include <numeric>
#include <vector>
#include <iterator>

#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include "TFitResult.h"
#include "TPaveText.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "FuncName.hh"
#include "DeleteUtility.hh"

#include "PidCommon.hh"
#include "PidData.hh"
#include "PidPdfMan.hh"

#define DrawHistToPdf 1
#define ScaleYieldElectronDeutron 1
#define ExtrapolateHighMomProton 1

namespace {

  // for writing objects to rootfile
  template <typename T>
  void WriteClonedObject(T* obj_to_clone, TFile* output_file, const char* name_format, long long id)
  {
    if (!obj_to_clone || !output_file) return;
    auto clone = static_cast<T*>(obj_to_clone->Clone());
    if (!clone) return;
    clone->SetName(Form(name_format, id));
    output_file->WriteObject(clone, clone->GetName());
    delete clone;
  }
  
  template <typename T>
  void WriteClonedObject(const T* obj_to_clone, TFile* output_file, const char* final_name) {
    if (!obj_to_clone || !output_file) return;
    auto clone = static_cast<T*>(obj_to_clone->Clone(final_name));
    clone->SetName(final_name);
    clone->SetTitle(final_name);    
    if (!clone) return;
    output_file->WriteObject(clone, clone->GetName());    
    delete clone;
  }

  // for making good graph
  inline constexpr int ngraph = static_cast<int>(CorrGraph::Graph::COUNT);
  inline constexpr double cutgraph[ngraph]
  = {0.1, 10, 0.01, 1.0, 0.1, 400,
    0.1, 10, 0.01, 1.0, 0.1, 400 };
  inline constexpr double cutgm2 = cutgraph[0];
  inline constexpr double cutgdedx = cutgraph[1]; 

  bool scaleEleDeu = false;
  bool extrapHighP = false;
}

// PidPdfMan
void PidPdfMan::Cleanup() {
    if (s_instance) {
        s_instance->Finalize();
    }
}
FitResultTable PidPdfMan::tbl_{};

//-----------------------------------------------------------------------
PidPdfMan* PidPdfMan::s_instance = nullptr;
PidPdfMan& PidPdfMan::GetInstance() 
{
  if (!s_instance) {
    s_instance = new PidPdfMan();
    std::atexit(PidPdfMan::Cleanup);
  }
  return *s_instance;
}


bool PidPdfMan::HasInstance() {
    return s_instance != nullptr;
}


PidPdfMan::PidPdfMan()
  : m_is_ready(false),
    m_file_name(),
    m_mainfile(nullptr),
    m_paramfile(nullptr)
{
}
PidPdfMan::~PidPdfMan()
{
}
void PidPdfMan::Finalize() {
    if (!m_is_ready) return;
    m_data.reset();
    m_is_ready = false;
    std::cout << "PidPdfMan finalized." << std::endl;
    delete s_instance;
    s_instance = nullptr;
}


bool PidPdfMan::Initialize()
{
#if ScaleYieldElectronDeutron
  scaleEleDeu = true;
#endif
#if ExtrapolateHighMomProton
  extrapHighP = true;
#endif
  
  if(m_is_ready){
    std::cerr << FUNC_NAME
	      << " already initialied" << std::endl;
    return false;
  }
  m_data = std::make_unique<PidData>();

  if(!OpenParamFile()) return false;
  Run();
  DrawResultsToPdf("anafig/tpc/pidlikeli/PidPdfSummary.pdf");
  if(!CloseCurrentFile()) return false;

  m_is_ready = true;

  return true;
}

//_____________________________________________________________________________
Bool_t
PidPdfMan::CdMainFile()
{
  if(m_mainfile){
    m_mainfile->cd();
    return true;
  } else {
    return false;
  }  
}

Bool_t
PidPdfMan::OpenParamFile()
{
  // open pidpdf file
  TString file_name = m_file_name;
  m_file.reset( TFile::Open(file_name,"READ") );
  std::cout << m_file->GetName() << std::endl;
  if(!m_file || m_file->IsZombie()){
    std::cerr << FUNC_NAME << " file open fail : "
	      << m_file_name << std::endl;
    return false;
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
PidPdfMan::ResetParamFile()
{
  if(m_file){
    std::cout << "Close " << m_file->GetName() << std::endl;        
    m_file.reset();
    std::cout << "Done " << std::endl;
    return true;      
  }
}

Bool_t
PidPdfMan::CloseCurrentFile()
{
  if(m_file){
    std::cout << "Close " << m_file->GetName() << std::endl;        
    m_file.reset();
    std::cout << "Done " << std::endl;
    return true;      
  }
  return false;
}

//_____________________________________________________________________________
Bool_t
PidPdfMan::Initialize(const TString& file_name)
{

  m_file_name = file_name;
  return Initialize();
}

// ---------- GetHist by full id ----------
TH2D* PidPdfMan::GetHistById(long long id) const
{  
  TString hname = Form("h%lld",id);
  TH2D* hist_from_file = (TH2D*)m_file->Get(hname);
  if(!hist_from_file) return nullptr;
  return static_cast<TH2D*>(hist_from_file->Clone());
}

//_____________________________________________________________________________
bool PidPdfMan::LoadAllHists()
{
  return true;
  TIter nextKey( m_file->GetListOfKeys() );
  while( TKey* key = static_cast<TKey*>( nextKey() ) ){
    if( std::string(key->GetClassName()) != "TH2D" ) continue;

    TString hname = key->GetName();
    
    if(hname[0] != 'h') continue;
    Long64_t id = std::strtoll(hname.Data()+1, nullptr, 10);
    //TH2D* h = key->ReadObject<TH2D>();
    auto h = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(key->ReadObject<TH2D>()->Clone()));
    if(!h) continue;
    auto clone = static_cast<TH2D*>(h->Clone());    
    h->SetDirectory(nullptr);
    m_histMap[id] = std::move(h);
    
  }
  return !m_histMap.empty();
}

bool PidPdfMan::Run(){
  MakeGoodPdfPiKP(pidlikeli::kDCarbonKP);
  StoreGoodPdfAll();
  StorePriorProbs();
  return true;
}

void PidPdfMan::MakeGoodPdfPiKP(int runmode)
{
  if(runmode==pidlikeli::kDGeneral){ // default
    GetParam(kPion);
    GetParam(kProton);
    MakeGoodPdf(kPion);
    MakeGoodPdf(kProton);
    GetParamKaonUsingOtherParticles();
    MakeGoodPdfKm();
  } else if(runmode==pidlikeli::kDCarbonKK){ // KK on Carbon
    
    return;    
  } else if(runmode==pidlikeli::kDCarbonKP){ // Kp on Carbon
    GetParam(kPion);
    GetParam(kProton);
    MakeGoodPdf(kPion);
    MakeGoodPdf(kProton);
    GetParamKaonUsingOtherParticles();
    MakeGoodPdfKm();    
  } else {
    MakeGoodPdfPiKP(pidlikeli::kDGeneral);
  }

  return;
}

bool PidPdfMan::GetParam(int pid_int, int chgmode)
{
  for(int it_int = 0; it_int < kNtype; ++it_int){
    if(! (it_int==pidlikeli::kTypeGen||it_int==pidlikeli::kTypeLmd||it_int==pidlikeli::kTypeK0)) continue;
    for(int chg_int = 0; chg_int < kNchg; ++chg_int){
      if( (chgmode==kPlus&&chg_int!=kPlus || chgmode==kMinus&&chg_int!=kMinus) ) continue;
      for(int ibe_int = 0; ibe_int < kNbe; ++ibe_int){
        for(int m = 0; m < kNmom; ++m){
          auto type_e = static_cast<pidlikeli::DType>(it_int);
          auto pid_e  = static_cast<pidlikeli::Pid>(pid_int);
          auto chg_e  = static_cast<pidlikeli::Chg>(chg_int);
          auto be_e   = static_cast<pidlikeli::BE>(ibe_int);          
          auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, m);
	  if(!CreateFitFunction(it_int,pid_int,chg_int,ibe_int,m)) continue;
          if(!SetInitFitParsProj(it_int,pid_int,chg_int,ibe_int, m)) continue;
	  if(!ExecFitProj(it_int,pid_int,chg_int,ibe_int,m)) continue;	  
          if(!StoreFitResultProj(it_int, pid_int, chg_int, ibe_int, m)) continue;
          if(!SetGEPoint(it_int, pid_int, chg_int, ibe_int, m)) continue;
          b.valid = true;
        }
        if(!FitGEPoint(it_int, pid_int, chg_int, ibe_int)) continue;

      }
    }
  }
  return true;
}

bool PidPdfMan::GetParamKaonUsingOtherParticles()
{
  int pid_int = -1;
  int fit_pid_int = -1;

  for(int it_int = 0; it_int < kNtype; ++it_int){
    for(int chg_int = 0; chg_int < kNchg; ++chg_int){
      for(int ibe_int = 0; ibe_int < kNbe; ++ibe_int){
	if( !( (it_int==pidlikeli::kTypeGen&&chg_int==kMinus&&ibe_int==0)
	       || (it_int==pidlikeli::kTypeKm&&chg_int==kMinus) )) continue;
	if( it_int==pidlikeli::kTypeGen ){
	  pid_int=kKaon;
	  fit_pid_int=kAllParticles;
	}
	if( it_int==pidlikeli::kTypeKm ){
	  pid_int=kKaon;
	  fit_pid_int=kKaon;
	}
	for(int m = 0; m < kNmom; ++m){
	  
	  if( !CreateFitFunction(pidlikeli::kTypeKm,pid_int,chg_int,ibe_int,m) ) continue;
	  if(!SetInitFitParsDoubleGauss(pidlikeli::kTypeKm,pid_int,chg_int,ibe_int,m)) continue;
	  if(it_int==pidlikeli::kTypeGen){
	    if(!ExecFitDouble2DGaussAllParticles(pidlikeli::kTypeKm,pid_int,chg_int,ibe_int,m) ) continue;

	  }
	  if(it_int==pidlikeli::kTypeKm){
	    if(!ExecFitDouble2DGauss(pidlikeli::kTypeKm,pid_int,chg_int,ibe_int,m) ) continue;	    
	  }

	  if( !StoreFitResultKaon(pidlikeli::kTypeKm,pid_int,chg_int,ibe_int,m) ) continue;
	  if( !SetGEPoint(pidlikeli::kTypeKm,pid_int,chg_int,ibe_int,m) ) continue;
	}
	if( !FitGEPoint(pidlikeli::kTypeKm,pid_int,chg_int,ibe_int) ) continue;	    
      }
    }
  }
  return true;
}
bool PidPdfMan::CreateFitFunction(int it_int, int pid_int, int chg_int, int ibe_int, int m)
{
  auto type_e = static_cast<pidlikeli::DType>(it_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(pid_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(chg_int);
  auto be_e   = static_cast<pidlikeli::BE>(ibe_int);          
  auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, m);
  if(b.functions_created) return true;
  auto h2 = GetHist(it_int, pid_int, chg_int, ibe_int, m);
  if(!h2){
    if(it_int==pidlikeli::kTypeKm&&pid_int==pidlikeli::kKaon&&chg_int==pidlikeli::kMinus&&ibe_int==0){
      h2 = new TH2D(Form("h%dp%dc%dbe%dm%d",it_int,pid_int,chg_int,ibe_int,m),
		    Form("h%dp%dc%dbe%dm%d",it_int,pid_int,chg_int,ibe_int,m),
		    pidlikeli::nbinm2,pidlikeli::minm2,pidlikeli::maxm2,
		    pidlikeli::nbindedx,pidlikeli::mindedx,pidlikeli::maxdedx );
    } else {
      delete h2;
      return false;
    }
  }
  
  h2->SetDirectory(nullptr);
  //fx
  const char* fx_basename = PdfFuncName[static_cast<size_t>(PdfFuncType::FX_PROJ)];
  TString fx_fullname = Form("%s%d_%d_%d_%d_%d",fx_basename,
			     it_int, pid_int, chg_int, ibe_int, m);
  auto fx = new TF1(fx_fullname,
		    pidfunc::RotGauss2DProjXFit,
		    h2->GetXaxis()->GetXmin(),h2->GetXaxis()->GetXmax(),
		    pidfunc::kNparamGauss);	    
  b.set_func(PdfFuncType::FX_PROJ, fx);
  //fy
  const char* fy_basename = PdfFuncName[static_cast<size_t>(PdfFuncType::FY_PROJ)];
  TString fy_fullname = Form("%s%d_%d_%d_%d_%d",fy_basename,
			     it_int, pid_int, chg_int, ibe_int, m);
  auto fy = new TF1(fy_fullname,
		    pidfunc::RotGauss2DProjYFit,
		    h2->GetYaxis()->GetXmin(),h2->GetYaxis()->GetXmax(),
		    pidfunc::kNparamGauss);
  b.set_func(PdfFuncType::FY_PROJ, fy);
  //f2d
  const char* f2d_basename = PdfFuncName[static_cast<size_t>(PdfFuncType::F2D)];
  TString f2d_fullname = Form("%s%d_%d_%d_%d_%d",f2d_basename,
			      it_int, pid_int, chg_int, ibe_int, m);
  auto f2d = new TF2(f2d_fullname,
		     pidfunc::RotGauss2D,
		     h2->GetXaxis()->GetXmin(),h2->GetXaxis()->GetXmax(),
		     h2->GetYaxis()->GetXmin(),h2->GetYaxis()->GetXmax(),
		     pidfunc::kNparamGauss);
  b.set_func(PdfFuncType::F2D, f2d);
  //f2dcorr
  const char* f2dc_basename = PdfFuncName[static_cast<size_t>(PdfFuncType::F2DCORR)];
  TString f2dc_fullname = Form("%s%d_%d_%d_%d_%d",f2dc_basename,
			      it_int, pid_int, chg_int, ibe_int, m);
  auto f2dc = new TF2(f2dc_fullname,
		     pidfunc::RotGauss2D,
		     h2->GetXaxis()->GetXmin(),h2->GetXaxis()->GetXmax(),
		     h2->GetYaxis()->GetXmin(),h2->GetYaxis()->GetXmax(),
		     pidfunc::kNparamGauss);
  b.set_func(PdfFuncType::F2DCORR, f2dc);
  //fd2d
  const char* fd2d_basename = PdfFuncName[static_cast<size_t>(PdfFuncType::FD2D)];
  TString fd2d_fullname = Form("%s%d_%d_%d_%d_%d",fd2d_basename,
			       it_int, pid_int, chg_int, ibe_int, m);
  auto fd2d = new TF2(fd2d_fullname,
		      pidfunc::RotDoubleGauss2D,
		      h2->GetXaxis()->GetXmin(),h2->GetXaxis()->GetXmax(),
		      h2->GetYaxis()->GetXmin(),h2->GetYaxis()->GetXmax(),
		      pidfunc::kNparamDGauss);
  b.set_func(PdfFuncType::FD2D, fd2d);
  //pdf
  const char* pdf_basename = PdfFuncName[static_cast<size_t>(PdfFuncType::PDF)];
  TString pdf_fullname = Form("%s%d_%d_%d_%d_%d",pdf_basename,
			       it_int, pid_int, chg_int, ibe_int, m);
  auto pdf = new TF2(pdf_fullname,
		      pidfunc::RotGauss2D,
		      h2->GetXaxis()->GetXmin(),h2->GetXaxis()->GetXmax(),
		      h2->GetYaxis()->GetXmin(),h2->GetYaxis()->GetXmax(),
		      pidfunc::kNparamGauss);
  b.set_func(PdfFuncType::PDF, pdf);
  
  
  b.functions_created = true;
  delete h2;
  return true;
}

bool PidPdfMan::CreateFitFunctionFive2DGauss(int chg_int, int ibe_int, int m)
{
  int it_int = pidlikeli::kTypeGen;
  int ip_int = pidlikeli::kAllParticles;
  auto type_e = static_cast<pidlikeli::DType>(it_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(ip_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(chg_int);
  auto be_e   = static_cast<pidlikeli::BE>(ibe_int);          
  auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, m);
  if(b.five2dgauss_created) return true;
  auto h2 = GetHist(it_int, ip_int, chg_int, ibe_int, m);
  if(!h2){
      return false;
  }  
  h2->SetDirectory(nullptr);
  //ff2d
  const char* ff2d_basename = PdfFuncName[static_cast<size_t>(PdfFuncType::FFIVE2D)];
  TString ff2d_fullname = Form("%s%d_%d_%d_%d_%d",ff2d_basename,
			       it_int, ip_int, chg_int, ibe_int, m);
  auto ff2d = new TF2(ff2d_fullname,
		      pidfunc::RotFiveGauss2D,
		      h2->GetXaxis()->GetXmin(),h2->GetXaxis()->GetXmax(),
		      h2->GetYaxis()->GetXmin(),h2->GetYaxis()->GetXmax(),
		      pidfunc::kNparamGauss*pidlikeli::kNpid);
  b.set_func(PdfFuncType::FFIVE2D, ff2d);  
  // b.ff2d = new TF2(Form("ff2dt%dp%dc%dbe%dm%d",it_int,ip_int,chg_int,ibe_int,m),
  // 		   pidfunc::RotFiveGauss2D,
  // 		   h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax(),
  // 		   h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax(),
  // 		   pidfunc::kNparamGauss*pidlikeli::kNpid);
  b.five2dgauss_created = true;
  delete h2;
  return true;
}


bool PidPdfMan::ExecFitProj(int it_int,int pid_int,int chg_int,int ibe_int,int m)
{  
  auto type_e = static_cast<pidlikeli::DType>(it_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(pid_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(chg_int);
  auto be_e   = static_cast<pidlikeli::BE>(ibe_int);          
  auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, m);
  
  auto h2 = GetHist(it_int, pid_int, chg_int, ibe_int, m);
  if(!h2 || h2->GetEntries() < 1){
    delete h2;
    return false;
  }
  h2->SetDirectory(nullptr);
  TH1D *h1x = h2->ProjectionX("h1x");
  TH1D *h1y = h2->ProjectionY("h1y"); 
  std::cout << " h1x->GetEntries(): " << h1x->GetEntries()
	    << " h1x->GetMean(): " << h1x->GetMean()
	    << " h1x->GetStdDev(): " << h1x->GetStdDev() << std::endl;
  auto fx = b.get_func<TF1>(PdfFuncType::FX_PROJ);
  auto fy = b.get_func<TF1>(PdfFuncType::FY_PROJ);  
  //h1x->Fit(b.fx, "Q");
  h1x->Fit(fx, "Q");	    
  std::cout << " h1y->GetEntries(): " << h1y->GetEntries()
	    << " h1y->GetMean(): "    << h1y->GetMean()
	    << " h1y->GetStdDev(): "  << h1y->GetStdDev() << std::endl;	  
  //  h1y->Fit(b.fy, "Q");
  h1y->Fit(fy, "Q");
  delete h1x;
  delete h1y;  
  delete h2;
  
  return true;
}

bool PidPdfMan::ExecFitDouble2DGauss(int it_int,int pid_int,int chg_int,int ibe_int,int m)
{
  auto h = GetHist(it_int, pid_int, chg_int, ibe_int, m);
  if (!h || h->GetEntries() < 1) return false;
  h->SetDirectory(nullptr);

  auto type_e = static_cast<pidlikeli::DType>(it_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(pid_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(chg_int);
  auto be_e   = static_cast<pidlikeli::BE>(ibe_int);          
  auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, m);
  auto fd2d = b.get_func<TF2>(PdfFuncType::FD2D);
  //if(!b.fd2d) return false;
  if(!fd2d) return false;
  //h->Fit(b.fd2d,"Q");
  //h->Fit(b.fd2d,"NR");
  h->Fit(fd2d,"NR");
  delete h;
  return true;
}

bool PidPdfMan::ExecFitDouble2DGaussAllParticles(int it_int,int pid_int,int chg_int,int ibe_int,int m)
{
  int it_all_int = pidlikeli::kTypeGen;
  int pid_all_int = pidlikeli::kAllParticles;
  auto h = GetHist(it_all_int, pid_all_int, chg_int, ibe_int, m);
  if (!h || h->GetEntries() < 1) return false;
  h->SetDirectory(nullptr);

  auto type_e = static_cast<pidlikeli::DType>(it_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(pid_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(chg_int);
  auto be_e   = static_cast<pidlikeli::BE>(ibe_int);          
  auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, m);
  auto fd2d = b.get_func<TF2>(PdfFuncType::FD2D);
  //  if(!b.fd2d) return false;
  if(!fd2d) return false;
  //h->Fit(b.fd2d,"Q");
  //h->Fit(b.fd2d,"NR");
  h->Fit(fd2d,"NR");  
  delete h;
  return true;
}

bool PidPdfMan::SetInitFitParsDoubleGauss(int type, int ip, int chg, int be, int momid)
{
  int be_pi=0;
  auto h = GetHist(type, ip, chg, be, momid);
  if (!h || h->GetEntries() < 1){
    if( !(type==pidlikeli::kTypeKm&&ip==pidlikeli::kKaon
	  &&chg==pidlikeli::kMinus&&be==0) ){    
      return false;
    }
  }

  auto type_e = static_cast<pidlikeli::DType>(type);
  auto type_pi_e = static_cast<pidlikeli::DType>(pidlikeli::kTypeK0);
  auto pid_e = static_cast<pidlikeli::Pid>(ip);
  auto pid_pi_e = static_cast<pidlikeli::Pid>(kPion);
  auto chg_e = static_cast<pidlikeli::Chg>(chg);
  auto be_e = static_cast<pidlikeli::BE>(be);
  auto be_pi_e = static_cast<pidlikeli::BE>(be_pi);
  
  auto &bkm= m_data->at_bin(type_e, pid_e, chg_e, be_e, momid); // K-
  auto &bpm= m_data->at_bin(type_pi_e, pid_pi_e, chg_e, be_pi_e, momid); // pi-
  auto fkpi = bkm.get_func<TF2>(PdfFuncType::FD2D);
  auto fpi  = bpm.get_func<TF2>(PdfFuncType::F2D);
  // TF2*fkpi = bkm.fd2d;
  // TF2*fpi = bpm.f2d;
  if(!fkpi||!fpi) return false;
  for(int i=0; i<pidfunc::kNparamGauss-1; i++){
    fkpi->FixParameter(i+pidfunc::kNparamGauss, fpi->GetParameter(i));
  }
  // Fix/Set param
  int lastparamid=pidfunc::kParamIdYield;
  fkpi->SetParameter(lastparamid+pidfunc::kNparamGauss, fpi->GetParameter(lastparamid));
  if(!SetInitFitPars2D(type,kKaon,chg,be,momid,h,fkpi)) return false;
  delete h;    
  return true;  
}

bool PidPdfMan::SetInitFitParsProj(int type, int pid, int chg, int be, int momid)
{
  auto h = GetHist(type, pid, chg, be, momid);
  if (!h || h->GetEntries() < 1) return false;
  auto type_e = static_cast<pidlikeli::DType>(type);
  auto pid_e = static_cast<pidlikeli::Pid>(pid);
  auto chg_e = static_cast<pidlikeli::Chg>(chg);
  auto be_e = static_cast<pidlikeli::BE>(be);
  auto &b= m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
  for(int i=0; i<2; i++){
    TF1* f; 
    if(i==0) f = b.get_func<TF1>(PdfFuncType::FX_PROJ);
    if(i==1) f = b.get_func<TF1>(PdfFuncType::FY_PROJ);
  // set initial fitting parameter condition
    if(!h) return false;  
    if(!f) return false;  
    //f->SetParNames("x0","y0","sigx","sigy","theta","N");
    double meanx     = -999.;
    double meany     = -999.;
    double sigx      = -999.;
    double sigy      = -999.;
    double rotangle  = -999.;
    double count     = -999.;
    double mom = pidlikeli::BinToMom(momid);

    double cm2min = pidlikeli::cutm2min[pid] - 0.01;
    double cm2max = pidlikeli::cutm2max[pid] + 0.01;
    double cdemin = pidlikeli::cutdedxmin[pid] - 1.;
    double cdemax = pidlikeli::cutdedxmax[pid] + 1.;

    double massvalue = pidlikeli::mass[pid];
    double m2 = massvalue * massvalue;  
    meanx = m2;
    if(chg==1) meanx = m2*(-1.);
    meany = pidlikeli::dedx_pid(mom,pid);
    sigx=pidfunc::CalcSigM2(mom,pid);  
    sigy=pidfunc::CalcSigdEdx(mom,pid);  
    rotangle = 0.;
    int yield = h->GetEntries();
    int tot = yield;  
    f->SetParameters(meanx,meany,sigx,sigy,rotangle,yield);
    f->FixParameter(4,rotangle);
    if(chg==0){
      f->SetParLimits(0, cm2min, cm2max);
    } else {
      f->SetParLimits(0,-cm2max,-cm2min);
    }
  
    if(pid==kPion){ // pion
      f->SetParLimits(1,meany-10,meany+10);
      f->SetParLimits(2,0.001,1.0);
      f->SetParLimits(3,0,50.);
      //f->SetParLimits(4,0,TMath::Pi()/2.);
      f->SetParLimits(5,0.1,10*tot);
      
      f->FixParameter(0,meanx);
      f->FixParameter(1,meany);      
      f->FixParameter(2,sigx);
      f->FixParameter(3,sigy);                           
    }
    if(pid==kKaon){
      f->SetParLimits(1,meany-10,meany+10);
      f->SetParLimits(2,0.01,0.5);
      f->SetParLimits(3,0,50.);
      f->SetParLimits(5,0.1,10*tot);      
    }
    if(pid==kProton){
      f->SetParLimits(1,meany-10,meany+10);        
      f->SetParLimits(2,0,0.5);
      f->SetParLimits(3,0,50.);
      //f->SetParLimits(4,0,TMath::Pi()/2.);
      f->SetParLimits(5,0.1,10*tot);

      f->FixParameter(0,meanx);
      f->FixParameter(1,meany);      
      f->FixParameter(2,sigx);
      f->FixParameter(3,sigy);                           
    }
  }
  delete h;  
  return true;
}

bool PidPdfMan::ExecFitFive2DGauss(int chg, int be, int momid)
{
  if(momid<25) return false;
  if(!CreateFitFunctionFive2DGauss(chg,be,momid)) return false;;    
  int it_int = pidlikeli::kTypeGen;
  int ip_int = pidlikeli::kAllParticles;
  auto h = GetHist(it_int, ip_int, chg, be, momid);
  if (!h || h->GetEntries() < 1){
    return false;
  }
  
  auto type_e = static_cast<pidlikeli::DType>(it_int);
  auto pid_e = static_cast<pidlikeli::Pid>(ip_int);
  auto chg_e = static_cast<pidlikeli::Chg>(chg);
  auto be_e = static_cast<pidlikeli::BE>(be);
  auto &b= m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
  if(!b.five2dgauss_created) return false;
  //TF2 * fitmodel = b.ff2d;
  auto fitmodel  = b.get_func<TF2>(PdfFuncType::FFIVE2D);

  for(const auto& config: pidlikeli::gFitConfigs) {
    int pid = static_cast<int>(config.pid);
    const TString& particle_name = pidlikeli::plist[pid];
    int param_offset = pid * pidfunc::kNparamGauss;
    for(int iparam=0; iparam<pidfunc::kNparamGauss; iparam++){
      TString param_name = CorrFunc::FuncNames[iparam];
      param_name.Remove(0,1);
      int global_param_index = pid*pidfunc::kNparamGauss+iparam;
      fitmodel->SetParName(global_param_index,
			Form("%s_%s", param_name.Data(),particle_name.Data()));      
    }
    
    // set initial parameter
    double param[pidfunc::kNparamGauss];
    std::array<double, pidfunc::kNparamGauss> params;
    if(pid<=pidlikeli::kProton) { // pi,K,p case
      const auto& stored_params = GetStoredFitParamsOfRefPdf(config, momid);
      if(stored_params.valid){
	for(int iparam=0; iparam<pidfunc::kNparamGauss; iparam++){
	  params[iparam] = stored_params.p[iparam];
	}
      } else {
	std::cerr << "Warning: No valid stored params for " << config.name << ". Using calculated values." << std::endl;
	params = GetCalcParametersRef(pid,momid);
      }      
    } else { // d,e case
      params = GetCalcParametersRef(pid,momid);
    }    
    if(chg==pidlikeli::kMinus) {
      params[pidfunc::kParamIdM2] = -1.0 * std::abs(params[pidfunc::kParamIdM2]);
    } else {
      params[pidfunc::kParamIdM2] = std::abs(params[pidfunc::kParamIdM2]);
    }
    
    int n = pidfunc::kNparamGauss-pidfunc::kNparamYield;
    for(int iparam=0; iparam<n; iparam++){
      // fitmodel->FixParameter(param_offset+iparam,params[iparam]);      
      fitmodel->SetParameter(param_offset+iparam,params[iparam]);
      double prange = params[iparam]*0.5 + 1e-6;
      double parammin = params[iparam] - prange;
      double parammax = params[iparam] + prange;      
      fitmodel->SetParLimits(param_offset+iparam,parammin,parammax);
      if(iparam!=pidfunc::kParamIdM2&&parammin<0.0){
	fitmodel->SetParLimits(param_offset+iparam,0.0,parammax);	
      }
      fitmodel->FixParameter(param_offset+pidfunc::kParamIdRotAngle, 0.0);
    }
    
    // Integral in ROI
    double m2min_roi = params[pidfunc::kParamIdM2]-2.0*params[pidfunc::kParamIdSigM2];
    double m2max_roi = params[pidfunc::kParamIdM2]+2.0*params[pidfunc::kParamIdSigM2];
    double dedxmin_roi = params[pidfunc::kParamIddEdx]-2.0*params[pidfunc::kParamIdSigdEdx];
    double dedxmax_roi = params[pidfunc::kParamIddEdx]+2.0*params[pidfunc::kParamIdSigdEdx];
    int binxmin = h->GetXaxis()->FindBin(m2min_roi);
    int binxmax = h->GetXaxis()->FindBin(m2max_roi);
    int binymin = h->GetYaxis()->FindBin(dedxmin_roi);
    int binymax = h->GetYaxis()->FindBin(dedxmax_roi);
    double yield_init = h->Integral(binxmin,binxmax,binymin,binymax);
    if(yield_init<1.0) yield_init = 0.0;
    // Set yield parameter
    fitmodel->SetParameter(param_offset+n, yield_init); 
    fitmodel->SetParLimits(param_offset+n, 0.0, h->GetEntries() * 10);
    
    if(pid==pidlikeli::kDeutron
       ||(pid==pidlikeli::kElectron)
       ||(pid==pidlikeli::kProton)
       ||(pid==pidlikeli::kKaon)
       ||(pid==pidlikeli::kProton&&chg==1)
       ||(pid==pidlikeli::kKaon&&chg==0)
       ){
      //    if(pid!=pidlikeli::kPion){
      fitmodel->FixParameter(param_offset+n, 0.0);
    }
  }
  //  TFitResultPtr fit_result = h->Fit(fitmodel, "SMLI");
  TFitResultPtr fit_result = h->Fit(fitmodel, "SMLQ");
  if (fit_result && fit_result->Status() == 0) {
    auto& cell = PidPdfMan::at(pidlikeli::kTypeGen, pidlikeli::kAllParticles, chg, be, momid);
    
    for (int i = 0; i < fitmodel->GetNpar(); ++i) {
      cell.p[i]    = fit_result->Parameter(i);
      cell.perr[i] = fit_result->ParError(i);
    }
    cell.chi2ndf = fit_result->Chi2() / fit_result->Ndf();
    cell.status  = fit_result->Status();
    cell.valid   = true;
 
    delete h;
    delete fitmodel;
    return true;
  }
  if(!fit_result) std::cout << " Null fit result " << std::endl;
  if(fit_result && fit_result->Status() != 0) std::cout << " fit result: "
					  << fit_result->Status() << std::endl;  
  std::cerr << "Warning: Fit failed for chg=" << chg << ", mom=" << momid << std::endl;
  delete h;
  delete fitmodel;
  return false;
}

bool
PidPdfMan::ExecFitProj1DGaussForYield(int chg, int be, int momid, 
				      std::vector<double>& params, std::vector<double>& errors)
{
  int it_int = pidlikeli::kTypeGen;
  int ip_int = pidlikeli::kAllParticles;
  auto type_e = static_cast<pidlikeli::DType>(it_int);
  auto pid_e = static_cast<pidlikeli::Pid>(ip_int);
  auto chg_e = static_cast<pidlikeli::Chg>(chg);
  auto be_e = static_cast<pidlikeli::BE>(be);
  auto &b= m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
  
  TH2* h2d = GetHist(pidlikeli::kTypeGen, pidlikeli::kAllParticles, chg, be, momid);
  if (!h2d || h2d->GetEntries() < 200) {
    if(h2d) delete h2d;
    return false;
  }
  TH1D* h_m2 = h2d->ProjectionX("h_m2_proj");
  TH1D* h_dedx = h2d->ProjectionY("h_dedx_proj");

  const int n_params_total = pidlikeli::kNpid * pidfunc::kNparamGauss;
  TF1* fit_m2 = new TF1("fit_m2",
			pidfunc::Penta1DGaussForM2,
			h_m2->GetXaxis()->GetXmin(),
			h_m2->GetXaxis()->GetXmax(),
			n_params_total);
  TF1* fit_dedx = new TF1("fit_dedx",
			  pidfunc::Penta1DGaussFordEdx,
			  h_dedx->GetXaxis()->GetXmin(),
			  h_dedx->GetXaxis()->GetXmax(),
			  n_params_total);  
  
  for (const auto& config : pidlikeli::gFitConfigs) {
    int pid = static_cast<int>(config.pid);
    std::array<double, pidfunc::kNparamGauss> params = GetCalcParametersRef(pid, momid);
    int offset = pid * pidfunc::kNparamGauss;
    const TString& particle_name = pidlikeli::plist[pid];
    for(int iparam=0; iparam<pidfunc::kNparamGauss; iparam++){
      TString param_name = CorrFunc::FuncNames[iparam];
      param_name.Remove(0,1);
      int global_param_index = pid*pidfunc::kNparamGauss+iparam;
      fit_m2->SetParName(global_param_index,
			 Form("%s_%s", param_name.Data(),particle_name.Data()));
      fit_dedx->SetParName(global_param_index,
			   Form("%s_%s", param_name.Data(),particle_name.Data()));
    }
    for(int iparam=0; iparam<pidfunc::kNparamGauss; iparam++){
      fit_m2->SetParameter(offset+iparam, params[iparam]);
      fit_dedx->SetParameter(offset+iparam, params[iparam]);
      double m2range = 0.1;
      double derange = 15;	
      double pmin = params[iparam] - m2range;
      double pmax = params[iparam] + m2range;
      if(iparam==pidfunc::kParamIdM2||iparam==pidfunc::kParamIddEdx){
	fit_m2->SetParLimits(offset+iparam,params[iparam]-m2range,params[iparam]+m2range);
	fit_dedx->SetParLimits(offset+iparam,params[iparam]-derange,params[iparam]+derange);
	if(pmin<0){
	  fit_m2->SetParLimits(offset+iparam,0.0,params[iparam]+m2range);
	  fit_dedx->SetParLimits(offset+iparam,0.0,params[iparam]+derange);
	}
      }
      if(iparam==pidfunc::kParamIdSigM2){
	fit_m2->SetParLimits(offset+iparam,0.0,5.0);
	fit_dedx->FixParameter(offset+iparam,0.0);	  
      }
      if(iparam==pidfunc::kParamIdSigdEdx){
	fit_m2->FixParameter(offset+iparam,0.0);	  
	fit_dedx->SetParLimits(offset+iparam,0.0,40);
      }	
      fit_m2->SetParLimits(offset+pidfunc::kParamIdYield, 0., h2d->GetEntries()*10);
      fit_dedx->SetParLimits(offset+pidfunc::kParamIdYield, 0., h2d->GetEntries()*10);
      if(!( (pid==pidlikeli::kPion)
	    ||(pid==pidlikeli::kProton&&chg==0)
	    // ||(pid==pidlikeli::kDeutron&&chg==0)
	    // ||(pid==pidlikeli::kElectron)	    	    
	    ||(pid==pidlikeli::kKaon&&chg==1) ) )
	{
	  //if(pid!=pidlikeli::kPion&&(){
	  fit_dedx->FixParameter(offset+iparam, 0.);
	  fit_m2->FixParameter(offset+iparam, 0.);
	}
    }
    fit_dedx->FixParameter(offset+pidfunc::kParamIdRotAngle, 0.);
  }
  TFitResultPtr fit_res_m2 = h_m2->Fit(fit_m2, "SLQ");
  TFitResultPtr fit_res_dedx = h_dedx->Fit(fit_dedx, "SLQ");
  
  bool success
    = fit_res_m2.Get() && fit_res_m2->IsValid()
    && fit_res_dedx.Get() && fit_res_dedx->IsValid();
  if (success) {
    auto& storage_cell = m_data->at_bin(type_e,pid_e,chg_e,be_e,momid);
    storage_cell.set_func(PdfFuncType::F1D_M2, static_cast<TF1*>(fit_m2->Clone()));
    storage_cell.set_func(PdfFuncType::F1D_DEDX, static_cast<TF1*>(fit_dedx->Clone()));
    storage_cell.valid = true;

    params.clear();
    errors.clear();
    params.reserve(n_params_total);
    errors.reserve(n_params_total);
    
    auto& cell = PidPdfMan::at(pidlikeli::kTypeGen, pidlikeli::kAllParticles, chg, be, momid);    
    for(int i = 0; i < n_params_total; ++i) {
      if(fit_res_m2->Parameter(pidfunc::kParamIdYield)<10) return false;
      int param_type = i % pidfunc::kNparamGauss;
      if (param_type == pidfunc::kParamIddEdx
	  || param_type == pidfunc::kParamIdSigdEdx) {
	params.push_back(fit_res_dedx->Parameter(i));
	errors.push_back(fit_res_dedx->ParError(i));
      } else {
	params.push_back(fit_res_m2->Parameter(i));
	errors.push_back(fit_res_m2->ParError(i));
      }
    }
    cell.valid = true;
    cell.status = fit_res_m2->Status();
  }  
  
  delete h2d;
  delete h_m2;
  delete h_dedx;
  delete fit_m2;
  delete fit_dedx;
  
  return success;
}

void PidPdfMan::ProcessFitResultsAndSetGEPointForYield(int chg, int be)
{
    int it = pidlikeli::kTypeGen;
    for (int im = 0; im < kNmom; ++im) {
      const auto& cell_all = PidPdfMan::at(it, pidlikeli::kAllParticles, chg, be, im);
      
        if (!cell_all.valid) continue;
        for (int ip = 0; ip < pidlikeli::kNpid; ++ip) {
            int offset = ip * pidfunc::kNparamGauss;            
            auto& cell_pid = PidPdfMan::at(it, ip, chg, be, im);
            for (int i_param = 0; i_param < pidfunc::kNparamGauss; ++i_param) {
                cell_pid.p[i_param] = cell_all.p[offset + i_param];
                cell_pid.perr[i_param] = cell_all.perr[offset + i_param];
            }
            cell_pid.valid = true;
            SetGEPoint(it, ip, chg, be, im);
        }
    }
}

bool PidPdfMan::SetInitFitPars2D(int type, int pid, int chg, int be, int momid, TH2D*h, TF2*f)
{
  if(!f) return false;  
  //f->SetParNames("x0","y0","sigx","sigy","theta","N");
  double meanx     = -999.;
  double meany     = -999.;
  double sigx      = -999.;
  double sigy      = -999.;
  double rotangle  = -999.;
  double count     = -999.;
  double mom = pidlikeli::BinToMom(momid);

  
  double cm2min = pidlikeli::cutm2min[pid] - 0.01;
  double cm2max = pidlikeli::cutm2max[pid] + 0.01;
  double cdemin = pidlikeli::cutdedxmin[pid] - 1.;
  double cdemax = pidlikeli::cutdedxmax[pid] + 1.;

  double massvalue = pidlikeli::mass[pid];
  double m2 = massvalue * massvalue;  
  meanx = m2;
  if(chg==1) meanx = (-1.)*m2;
  meany = pidlikeli::dedx_pid(mom,pid);
  sigx=pidfunc::CalcSigM2(mom,pid);  
  sigy=pidfunc::CalcSigdEdx(mom,pid);  
  rotangle = 0.;
  int yield = 0.;
  if(h){
    yield = h->GetEntries();
  } else {
    yield = 100;
  }
  int tot = yield;  
  f->SetParameters(meanx,meany,sigx,sigy,rotangle,yield);
  f->FixParameter(4,rotangle);
  if(chg==0){
    f->SetParLimits(0, cm2min, cm2max);
  } else {
    f->SetParLimits(0,-cm2max,-cm2min);
  }
  if(pid==kPion){ // pion
    f->SetParLimits(1,meany-10,meany+10);
    f->SetParLimits(2,0.001,1.0);
    f->SetParLimits(3,0,50.);
    //f->SetParLimits(4,0,TMath::Pi()/2.);
    f->SetParLimits(5,0.1,10*tot);      
    //f->FixParameter(0,meanx);                
  }
  if(pid==kKaon){
    f->SetParLimits(1,meany-10,meany+10);
    f->SetParLimits(2,0.01,0.5);
    f->SetParLimits(3,0,50.);
    //f->SetParLimits(4,0,TMath::Pi()/2.);
    f->SetParLimits(5,0.1,10*tot);      
    //f->FixParameter(0,meanx);
    //f->FixParameter(2,sigx);
    //f->FixParameter(0,meanx);
    //f->FixParameter(1,meany);      
    f->FixParameter(2,sigx);
    f->FixParameter(3,sigy);                               
  }
  if(pid==kProton){
    f->SetParLimits(1,meany-10,meany+10);        
    f->SetParLimits(2,0,0.5);
    f->SetParLimits(3,0,50.);
    //f->SetParLimits(4,0,TMath::Pi()/2.);
    f->SetParLimits(5,0.1,10*tot);
    //f->FixParameter(0,meanx);   
  }
  return true;
}

// ---------- GetHist by 5 components ----------
TH2D* PidPdfMan::GetHist(int type,int pid,int charge,int be,int mom) const
{
  //if( !OpenParamFile() ) return nullptr;
  //  if(!m_is_ready) return nullptr;
  if (   type<0   || type>=kNtype
	 || pid<0    || pid>=kNpidAll // pi,k,p,d,e+all
	 || charge<0 || charge>=kNchg
	 || be<0     || be>=kNbe
	 || mom<0    || mom>=kNmom ) return nullptr; 
  long long id = fac_t*(type+1)+fac_p*pid+fac_c*charge+fac_b*be+fac_m*mom;
  return GetHistById(id);
}

bool PidPdfMan::CheckFitFunc(int type_int, int pid_int, int chg_int, int be_int, int mom_int) const
{
    auto type_e = static_cast<pidlikeli::DType>(type_int);
    auto pid_e  = static_cast<pidlikeli::Pid>(pid_int);
    auto chg_e  = static_cast<pidlikeli::Chg>(chg_int);
    auto be_e   = static_cast<pidlikeli::BE>(be_int);

    const auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, mom_int);
    auto f2d = b.get_func<TF2>(PdfFuncType::F2D);
    return b.valid && f2d;
}

bool PidPdfMan::CheckFitFunc() const
{
  for(int it_int=0; it_int<kNtype; it_int++){
    for(int ip_int=0; ip_int<kNpidAll; ip_int++){
      for(int ic_int=0; ic_int<kNchg; ic_int++){
	for(int ibe_int=0; ibe_int<kNbe; ibe_int++){
	  for(int im=0; im<kNmom; im++){
	    bool check = CheckFitFunc(it_int,ip_int,ic_int,ibe_int,im);
	    //const auto &b = m_bin[it][ip][ic][ibe][im];
	    auto type_e = static_cast<pidlikeli::DType>(it_int);
	    auto pid_e  = static_cast<pidlikeli::Pid>(ip_int);
	    auto chg_e  = static_cast<pidlikeli::Chg>(ic_int);
	    auto be_e   = static_cast<pidlikeli::BE>(ibe_int);	    
	    const auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, im);
	    if(!check){
	    }
	  }
	}
      }
    }
  }
  return true;
}
// ----------- FitResult -----------------------------------------------
bool
PidPdfMan::StoreFitResultKaon(int it_int,int ip_int,int ic_int,int ibe_int,int momid)
{
  auto type_e = static_cast<pidlikeli::DType>(it_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(ip_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(ic_int);
  auto be_e   = static_cast<pidlikeli::BE>(ibe_int);
  const auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
  auto f = b.get_func<TF2>(PdfFuncType::FD2D);
  if(!f) return false;
  auto &cell = PidPdfMan::at(it_int,ip_int,ic_int,ibe_int,momid);
  if(!cell.valid){
    // if(1)
    //   throw std::runtime_error("Input value cannot be negative.");
    for(int i=0; i<pidfunc::kNparamGauss; ++i){
      double param = f->GetParameter(i);
      double perror = f->GetParError(i);    
      cell.p[i]    = param;
      cell.perr[i] = perror;
      std::cout << "param,error: " << param << "," << perror << std::endl;    
    }  
    cell.chi2ndf = f->GetChisquare()/f->GetNDF();
    cell.valid   = true;
    //    TF2* f2d = b.f2d;
    auto f2d  = b.get_func<TF2>(PdfFuncType::F2D);
    for(int i=0; i<pidfunc::kNparamGauss; i++){
      f2d->SetParameter(i,f->GetParameter(i));
    }
  }
  return true;
}

bool
PidPdfMan::StoreFitResultProj(int it,int pid,int ic,int be,int momid)
{
  auto type_e = static_cast<pidlikeli::DType>(it);
  auto pid_e  = static_cast<pidlikeli::Pid>(pid);
  auto chg_e  = static_cast<pidlikeli::Chg>(ic);
  auto be_e   = static_cast<pidlikeli::BE>(be);	      
  const auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
  auto fx = b.get_func<TF2>(PdfFuncType::FX_PROJ);
  auto fy = b.get_func<TF2>(PdfFuncType::FY_PROJ);  
  if(!fx||!fy) return false;
  auto &cell = PidPdfMan::at(it,pid,ic,be,momid);
  for(int i=0;i<2;++i){
    double paramx = fx->GetParameter(i*2);
    double perrorx = fx->GetParError(i*2);
    double paramy = fy->GetParameter(i*2+1);
    double perrory = fy->GetParError(i*2+1);        
    cell.p[i*2]    = paramx;
    cell.perr[i*2] = perrorx;
    cell.p[i*2+1]    = paramy;
    cell.perr[i*2+1] = perrory;    
    std::cout << "paramx,y / errorx,y: "
	      << paramx << ", " << paramy
	      << " / " << perrorx << ", " << perrory << std::endl;
  }
  for(int i=0; i<2; i++){
    double paramx = fx->GetParameter(i+4);
    double perrorx = fx->GetParError(i+4);
    cell.p[i+4] = paramx;
    cell.perr[i+4] = perrorx;
  }
  std::cout << " fx yield: " << fx->GetParameter(5)
	    << " fy yield: " << fy->GetParameter(5) << std::endl;
  cell.chi2ndf = fx->GetChisquare()/fx->GetNDF();
  //cell.chi2ndf = fy->GetChisquare()/fy->GetNDF();  
  cell.valid   = true;
  
  // set param to f2d
  auto f2d = b.get_func<TF2>(PdfFuncType::F2D);
  f2d->SetParameters(cell.p);
  f2d->SetParErrors(cell.perr);
  
  return true;
}


bool PidPdfMan::CheckGoodFitdEdx(int t_int, int p_int, int c_int, int be_int, int momid){
  auto type_e = static_cast<pidlikeli::DType>(t_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(p_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(c_int);
  auto be_e   = static_cast<pidlikeli::BE>(be_int);
  auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
  return b.is_gooddEdx;
}

bool PidPdfMan::CheckGoodFitM2(int t_int, int p_int, int c_int, int be_int, int momid){
  auto type_e = static_cast<pidlikeli::DType>(t_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(p_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(c_int);
  auto be_e   = static_cast<pidlikeli::BE>(be_int);
  auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
  return b.is_goodM2;
}

bool PidPdfMan::CheckGoodFitM2dEdx(int t_int, int p_int, int c_int, int be_int, int momid) const
{
  auto type_e = static_cast<pidlikeli::DType>(t_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(p_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(c_int);
  auto be_e   = static_cast<pidlikeli::BE>(be_int);
  auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
  return b.is_goodM2dEdx;
}

bool PidPdfMan::MakeGoodPdf(int pid_int)
{
  for(int it=0; it < kNtype; it++){
    for(int ic=0; ic < kNchg; ic++){
      for(int ibe=0; ibe < kNbe; ibe++){
        for(int im=0; im < kNmom; im++){
	  auto type_e = static_cast<pidlikeli::DType>(it);
	  auto pid_e  = static_cast<pidlikeli::Pid>(pid_int);
	  auto chg_e  = static_cast<pidlikeli::Chg>(ic);
	  auto be_e   = static_cast<pidlikeli::BE>(ibe);
	  auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, im);
	  auto f2d = b.get_func<TF2>(PdfFuncType::F2D);
	  //if(!b.f2d||!b.f2dc) continue;
	  if(!f2d) continue;
	  double mom = pidlikeli::BinToMom(im);
	  auto f2dc = b.get_func<TF2>(PdfFuncType::F2DCORR);	  
          if(CheckGoodFitM2dEdx(it, pid_int, ic, ibe, im)&&mom<0.7){
            f2dc->SetParameters(f2d->GetParameters());
	    for(int i=0; i<pidfunc::kNparamGauss; i++){
	      std::cout << "f2d->GetParameter(" << i << "): " << f2d->GetParameter(i) << std::endl;
	      std::cout << "f2dc->GetParameter(" << i << "): " << f2dc->GetParameter(i) << std::endl;	      
	    }
          } else {
            std::array<double,pidfunc::kNparamGauss> params = GetCalcParameters(it,pid_int,ic,ibe,im);
            f2dc->SetParameters(params.data());
          }
        }
      }
    }
  }
  return true;
}

bool PidPdfMan::MakeGoodPdfKm()
{
  int pid_int=kKaon;
  for(int it=0; it < kNtype; it++){
    for(int ic=0; ic < kNchg; ic++){
      for(int ibe=0; ibe < kNbe; ibe++){
        for(int im=0; im < kNmom; im++){
	  auto type_e = static_cast<pidlikeli::DType>(it);
	  auto pid_e  = static_cast<pidlikeli::Pid>(pid_int);
	  auto chg_e  = static_cast<pidlikeli::Chg>(ic);
	  auto be_e   = static_cast<pidlikeli::BE>(ibe);
	  auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, im);
	  auto f2d = b.get_func<TF2>(PdfFuncType::F2D);
	  if(!f2d) continue;
	  auto f2dc = b.get_func<TF2>(PdfFuncType::F2DCORR);	  
          if(CheckGoodFitM2dEdx(it, pid_int, ic, ibe, im)){
            f2dc->SetParameters(f2d->GetParameters());	
          } else {
            std::array<double, pidfunc::kNparamGauss> params = GetCalcParameters(it, pid_int, ic, ibe, im);
            f2dc->SetParameters(params.data());
          }
        }
      }
    }
  }
  return true;
}

void PidPdfMan::StoreGoodPdfAll()
{

  std::cout << "Storing good PDFs with normalized yield (1.0)..." << std::endl;
  if (!m_data) {
    std::cerr << "Error in StoreGoodPdfAll: m_data is not initialized." << std::endl;
    return;
  }
  
  for (int it = 0; it < kNtype; ++it) {
    for (int ip = 0; ip < kNpidAll; ++ip) {
      for (int ic = 0; ic < kNchg; ++ic) {
	for (int ib = 0; ib < kNbe; ++ib) {
	  for (int im = 0; im < kNmom; ++im) {                        
	    auto& bin = m_data->at_bin(
				       static_cast<pidlikeli::DType>(it),
				       static_cast<pidlikeli::Pid>(ip),
				       static_cast<pidlikeli::Chg>(ic),
				       static_cast<pidlikeli::BE>(ib),
				       im
				       );
	    //TF2* source_func = bin.f2dc;
	    auto source_func  = bin.get_func<TF2>(PdfFuncType::F2DCORR);
	    if (!source_func) {
	      continue;
	    }
	    auto pdf = bin.get_func<TF2>(PdfFuncType::PDF);
	    if (!pdf) {
	      pdf = static_cast<TF2*>(source_func->Clone());
	      const char* pdf_basename
		= PdfFuncName[static_cast<size_t>(PdfFuncType::PDF)];
	      pdf->SetName(Form("%s_%s", pdf_basename,source_func->GetName()));
	    }
	    for (int i = 0; i < pidfunc::kParamIdYield; ++i) {
	      pdf->SetParameter(i, source_func->GetParameter(i));
	    }
	    pdf->SetParameter(pidfunc::kParamIdYield, 1.0);
	  }
	}
      }
    }
  }
  std::cout << "Finished storing good PDFs." << std::endl;
}


double PidPdfMan::GetYield(int it_int, int ip_int, int ic_int, int ib_int, int im)
{
  double yield = 0.;
  if(CheckGoodFitM2dEdx(it_int, ip_int, ic_int, ib_int, im)){
    yield = getFitYield(it_int, ip_int, ic_int, ib_int, im);
  } else {
    yield = InterpYield(it_int, ip_int, ic_int, ib_int, im);
  }
  return yield;
}

double PidPdfMan::GetTotalYield(int it_int, int ic_int, int ib_int, int im)
{
  double totalyield = 0.;
  for(int ip=0; ip<kNpid; ip++){
    totalyield += GetYield(it_int, ip, ic_int, ib_int, im);
  }
  return totalyield;
}

void PidPdfMan::FillPriorHist()
{
  for(int it=0; it<kNtype; it++){
    for(int ip=0; ip<kNpid; ip++){    
      for(int ic=0; ic<kNchg; ic++){
	for(int ib=0; ib<kNbe; ib++){
	  for(int im=0; im<kNmom; im++){
	    FillPriorHist(it, ip, ic, ib, im);
	  }
	}
      }
    }
  }
}

void PidPdfMan::FillPriorHist(int ip)
{
  for(int it=0; it<kNtype; it++){
    for(int ic=0; ic<kNchg; ic++){
      for(int ib=0; ib<kNbe; ib++){
	for(int im=0; im<kNmom; im++){
	  FillPriorHist(it, ip, ic, ib, im);
	}
      }
    }
  }
}

bool PidPdfMan::FillPriorHist(int it_int, int ip_int, int ic_int, int ib_int, int im)
{
  if( !m_data || !m_data->m_h5_priors) {
    std::cerr << "Error: Prior hist (THn) is already initialized. " << std::endl;
    return false;
  }
  double priorprob = GetPriorProb(it_int,ip_int,ic_int,ib_int,im);
  if(GetTotalYield(it_int,ic_int,ib_int,im)>100){
    int bin_coords[5] = {it_int+1, ip_int+1, ic_int+1, ib_int+1, im+1};
    long long global_bin = m_data->m_h5_priors->GetBin(bin_coords);
    m_data->m_h5_priors->SetBinContent(global_bin, priorprob);
  }

  return true;
}

std::vector<double> PidPdfMan::CalculatePriorProb(const TFitResultPtr& fit_result) const
{
  if (!fit_result) {
    std::cerr << "Error: Fit failed and returned a null pointer. Cannot calculate ratios." << std::endl;
    return {};
  } else if (fit_result->Status() != 0) {    
    std::cerr << "Warning: Fit did not converge properly (status=" << fit_result->Status() 
	      << "). Cannot calculate ratios." << std::endl;
    return {};
  }
  std::vector<double> yields;
  yields.reserve(pidlikeli::kNpid);

  for (int ip = 0; ip < pidlikeli::kNpid; ++ip) {
    int yield_param_index = ip * pidfunc::kNparamGauss + pidfunc::kParamIdYield;
    double current_yield = fit_result->Parameter(yield_param_index);
      
    if(current_yield<0 || std::isnan(current_yield) || std::isinf(current_yield)) {
      std::cerr << "Warning: Invalid yield value (" << current_yield
		<< ") found for PID " << ip << ". Cannot calculate ratios." << std::endl;
      return {};
    }
    yields.push_back(current_yield);
  }
  double total_yield = std::accumulate(yields.begin(), yields.end(), 0.0);

  if (total_yield < 1e-9) {
    std::cerr << "Warning: Total yield is nearly zero. Cannot calculate ratios." << std::endl;
    return {};
  }
  std::vector<double> ratios;
  ratios.reserve(pidlikeli::kNpid);
  for (double yield : yields) {
    double priprob = yield / total_yield;
    ratios.push_back(priprob);
    std::cout << " Prior Prob: " << priprob << std::endl;
  }

  return ratios;
}

void PidPdfMan::StorePriorProbs()
{
  std::cout << "Starting generation of prior probabilities..." << std::endl;

  int it_gen = pidlikeli::kTypeGen;
  if (!m_data || !m_data->m_h5_priors) {
    std::cerr << "Error: Cannot generate priors, PidPdfMan not ready." << std::endl;
    return;
  }
  
  for (int ic = 0; ic < pidlikeli::kNchg; ++ic) {
    for (int ib = 0; ib < pidlikeli::kNbe; ++ib) {
      if(ib!=0) continue;
      std::map<int, std::vector<double>> fit_params_map;
      std::map<int, std::vector<double>> fit_errors_map;
      std::map<int, bool> fit_success_map;

      // make good TG
      for (int im = 0; im < pidlikeli::kNmom; ++im) {
	std::vector<double> temp_params, temp_errors;
	bool fit_succeeded = ExecFitProj1DGaussForYield(ic, ib, im, temp_params, temp_errors);
	
	if (fit_succeeded) {
	  for (int ip = 0; ip < kNpid; ++ip) {
	    auto& cell_pid = PidPdfMan::at(it_gen, ip, ic, ib, im);
	    int offset = ip * pidfunc::kNparamGauss;
	    for (int i_param = 0; i_param < pidfunc::kNparamGauss; ++i_param) {
	      cell_pid.p[i_param] = temp_params[offset + i_param];
	      cell_pid.perr[i_param] = temp_errors[offset + i_param];
	    }
	    cell_pid.valid = true;
	  }	  
	  if(scaleEleDeu){
	    double yieldfac_e  = pidlikeli::yieldfac_e;
	    double yieldfac_d  = pidlikeli::yieldfac_d;	  
	    auto& cell_pi = PidPdfMan::at(it_gen, pidlikeli::kPion, ic, ib, im);
	    auto& cell_p  = PidPdfMan::at(it_gen, pidlikeli::kProton, ic, ib, im);
	    auto& cell_e = PidPdfMan::at(it_gen, pidlikeli::kElectron, ic, ib, im);
	    auto& cell_d = PidPdfMan::at(it_gen, pidlikeli::kDeutron, ic, ib, im);  
	    cell_e.p[pidfunc::kParamIdYield] = cell_pi.p[pidfunc::kParamIdYield]*yieldfac_e;
	    cell_e.perr[pidfunc::kParamIdYield] = cell_pi.perr[pidfunc::kParamIdYield]*yieldfac_e;
	    if (ic == pidlikeli::kPlus) {
	      cell_d.p[pidfunc::kParamIdYield] = cell_p.p[pidfunc::kParamIdYield]*yieldfac_d;
	      cell_d.perr[pidfunc::kParamIdYield] = cell_p.perr[pidfunc::kParamIdYield]*yieldfac_d;
	    } else {
	      cell_d.p[pidfunc::kParamIdYield] = 0.0;
	      cell_d.perr[pidfunc::kParamIdYield] = 0.0;
	    }
	  }
	  for (int ip=0; ip<pidlikeli::kNpid; ip++){
	    SetYieldGEPoint(it_gen, ip, ic, ib, im);
	  }
	}
      }
      for (int im = 0; im < pidlikeli::kNmom; ++im) {
	printf("\n--- Processing: Chg=%d, BE=%d, Mom=%.2f GeV/c ---\n", ic, ib, pidlikeli::BinToMom(im));
	double total_yield_this_mom = 0;
	std::vector<double> yields(kNpid);

	for (int ip = 0; ip < kNpid; ++ip) {
	  auto& cell_pid = PidPdfMan::at(it_gen, ip, ic, ib, im);
	  double yield_val = 0.0;
	  bool was_good_fit = cell_pid.valid && CheckGoodFitM2dEdx(it_gen, ip, ic, ib, im);
	  double mom = pidlikeli::BinToMom(im);	  
	  double mom_lth = pidlikeli::extrapolateLowMomProton;
	  double mom_hth = pidlikeli::extrapolateHighMomProton;	  
	  if (extrapHighP && ip==pidlikeli::kProton && ic==pidlikeli::kPlus && (mom < mom_lth || mom > mom_hth) ) {
            was_good_fit = false;
	  }
	  if(cell_pid.valid) std::cout << "cell valid " << std::endl;
	  if(CheckGoodFitM2dEdx(it_gen, ip, ic, ib, im)) std::cout << "goodfit m2 dedx  " << std::endl;
	  if (was_good_fit) {
	    yield_val = cell_pid.p[pidfunc::kParamIdYield];
	  } else {
	    yield_val = InterpYield(it_gen, ip, ic, ib, im);
	  }
	  if ( extrapHighP && ip==pidlikeli::kProton && mom>mom_hth) {
	    yield_val = InterpYield(it_gen, ip, ic, ib, im);
	    double yield_err = yield_val*0.1;
	    cell_pid.p[pidfunc::kParamIdYield] = yield_val;
	    cell_pid.perr[pidfunc::kParamIdYield] = yield_err;
	    cell_pid.valid = true;
	    if(scaleEleDeu){
	      double yieldfac_d  = pidlikeli::yieldfac_d;	  
	      auto& cell_d = PidPdfMan::at(it_gen, pidlikeli::kDeutron, ic, ib, im);	    
	      if (ic == pidlikeli::kPlus) {
		cell_d.p[pidfunc::kParamIdYield] = yieldfac_d*yield_val;
		cell_d.perr[pidfunc::kParamIdYield] = yieldfac_d*yield_err;
	      }
	    }	    
	  }
	  if (yield_val < 10 || std::isnan(yield_val) || mom<0.1 ) yield_val = 0.1;
	  //if (yield_val < 10 || std::isnan(yield_val)) yield_val = InterpYield(it_gen, ip, ic, ib, im);
	  yields[ip] = yield_val;
	  total_yield_this_mom += yields[ip];
	  printf("  PID: %-8s | Yield: %10.2f (%s)\n",
		 pidlikeli::plist[ip].Data(), yield_val, was_good_fit ? "Fit" : "Interpolated");
	}
	if (total_yield_this_mom > 1e-9) {
	  printf("  --------------------------------------------------\n");
	  for (int ip = 0; ip < kNpid; ++ip) {
	    double prior = yields[ip] / total_yield_this_mom;
	    printf("  PID: %-8s | Prior Probability: %.4f\n", pidlikeli::plist[ip].Data(), prior);
	    int bin_coords[5] = {pidlikeli::kTypeForPrior + 1, ip + 1, ic + 1, ib + 1, im + 1};
	    m_data->m_h5_priors->SetBinContent(m_data->m_h5_priors->GetBin(bin_coords), prior);
	  }
	}
      }
    }
    std::cout << "Finished generation of prior probabilities." << std::endl;
  }
}

void PidPdfMan::CalculateAndStorePriorProbs()
{
  std::cout << "Starting generation of prior probabilities..." << std::endl;
  if (!m_data || !m_data->m_h5_priors) {
    std::cerr << "Error: Cannot generate priors, PidPdfMan not ready." << std::endl;
    return;
  }

  for (int ic = 0; ic < pidlikeli::kNchg; ++ic) {
    for (int ib = 0; ib < pidlikeli::kNbe; ++ib) {
      if(ib!=0) continue;
      for (int im = 0; im < pidlikeli::kNmom; ++im) {
	// // Fitting
	if(!CreateFitFunctionFive2DGauss(ic,ib,im)) continue;  
	Bool_t fit_res = ExecFitFive2DGauss(ic, ib, im);
	//Bool_t fit_res = ExecFitSimultaneous1DGauss(ic, ib, im);
	if(fit_res){
	  std::vector<double> yields;
	  yields.reserve(pidlikeli::kNpid);
	  for (int ip = 0; ip < kNpid; ++ip) {
	    const auto& cell = PidPdfMan::at(pidlikeli::kTypeGen, ip, ic, ib, im);
	    double yield_val = cell.p[pidfunc::kParamIdYield];
	    yields.push_back(yield_val > 0 ? yield_val : 0.0);
	  }
	  double total_yield = std::accumulate(yields.begin(), yields.end(), 0.0);
	  if (total_yield > 1e-9) {
	    int type_idx = pidlikeli::kTypeForPrior;
	    for (int ip = 0; ip < kNpid; ++ip) {
	      double prior = yields[ip] / total_yield;
	      int bin_coords[5] = {type_idx + 1, ip + 1, ic + 1, ib + 1, im + 1};
	      long long global_bin = m_data->m_h5_priors->GetBin(bin_coords);
	      m_data->m_h5_priors->SetBinContent(global_bin, prior);
	    }
	  }
	} else {
	  std::cout << "Skipping prior calculation for chg=" << ic << ", mom=" << im 
		    << " due to RooFit failure." << std::endl;
	}	
      }
    }
  }
  std::cout << "Finished generation of prior probabilities." << std::endl;
}


void
PidPdfMan::GetPriorProb()
{
  double nInt[kNpid];
  double nTotal = 0.;
  // --- prior ---
  for(int it=0; it<kNtype; it++){
    for(int ip=0; ip<kNpid; ip++){    
      for(int ic=0; ic<kNchg; ic++){
	for(int ib=0; ib<kNbe; ib++){
	  for(int im=0; im<kNmom; im++){
	    double priorprob = GetPriorProb(it,ip,ic,ib,im);
	  }
	}
      }
    }
  }
}
double
PidPdfMan::GetPriorProb(int it_int, int ip_int, int chg_int, int ibe_int, int im_int)
{
  double nInt[kNpid];
  double nTotal = 0.;
  double priorprob = 0.;
  // pid
  for(int pid_int=0; pid_int<kNpid; ++pid_int){
    nInt[pid_int] = GetYield(it_int,pid_int,chg_int,ibe_int,im_int);
    nTotal += nInt[pid_int];
  }


  auto type_e = static_cast<pidlikeli::DType>(it_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(ip_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(chg_int);
  auto be_e   = static_cast<pidlikeli::BE>(ibe_int);          
  //  auto &prior = m_data->at_prior(type_e, pid_e, chg_e, be_e, im_int);
  double prior;
  prior = nInt[ip_int]/nTotal;
  priorprob = prior;
  if(prior<1e-10){
    prior = 1e-500;
  }
  return priorprob;
}

double PidPdfMan::getFitMuM2(int t,int p,int c,int b,int momid)
{
  
  auto &cell = PidPdfMan::at(t,p,c,b,momid);
  return cell.valid ? cell.p[0] : std::numeric_limits<double>::quiet_NaN();
}
double PidPdfMan::getFitErrMuM2(int t,int p,int c,int b,int momid)
{
  auto &cell = PidPdfMan::at(t,p,c,b,momid);
  return cell.valid ? cell.perr[0] : std::numeric_limits<double>::quiet_NaN();
}
double PidPdfMan::getFitMudEdx(int t,int p,int c,int b,int momid)
{
  auto &cell = PidPdfMan::at(t,p,c,b,momid);
  return cell.valid ? cell.p[1] : std::numeric_limits<double>::quiet_NaN();
}
double PidPdfMan::getFitErrMudEdx(int t,int p,int c,int b,int momid)
{
  auto &cell = PidPdfMan::at(t,p,c,b,momid);
  return cell.valid ? cell.perr[1] : std::numeric_limits<double>::quiet_NaN();
}
double PidPdfMan::getFitSigM2(int t,int p,int c,int b,int momid)
{
  auto &cell = PidPdfMan::at(t,p,c,b,momid);
  return cell.valid ? cell.p[2] : std::numeric_limits<double>::quiet_NaN();
}
double PidPdfMan::getFitErrSigM2(int t,int p,int c,int b,int momid)
{
  auto &cell = PidPdfMan::at(t,p,c,b,momid);
  return cell.valid ? cell.perr[2] : std::numeric_limits<double>::quiet_NaN();
}
double PidPdfMan::getFitSigdEdx(int t,int p,int c,int b,int momid)
{
  auto &cell = PidPdfMan::at(t,p,c,b,momid);
  return cell.valid ? cell.p[3] : std::numeric_limits<double>::quiet_NaN();
}
double PidPdfMan::getFitErrSigdEdx(int t,int p,int c,int b,int momid)
{
  auto &cell = PidPdfMan::at(t,p,c,b,momid);
  return cell.valid ? cell.perr[3] : std::numeric_limits<double>::quiet_NaN();
}
double PidPdfMan::getFitAngle(int t,int p,int c,int b,int momid)
{
  auto &cell = PidPdfMan::at(t,p,c,b,momid);
  return cell.valid ? cell.p[4] : std::numeric_limits<double>::quiet_NaN();
}
double PidPdfMan::getFitErrAngle(int t,int p,int c,int b,int momid)
{
  auto &cell = PidPdfMan::at(t,p,c,b,momid);
  return cell.valid ? cell.perr[4] : std::numeric_limits<double>::quiet_NaN();
}
double PidPdfMan::getFitYield(int t,int p,int c,int b,int momid)
{
  auto &cell = PidPdfMan::at(t,p,c,b,momid);
  return cell.valid ? cell.p[5] : std::numeric_limits<double>::quiet_NaN();
}
double PidPdfMan::getFitErrYield(int t,int p,int c,int b,int momid)
{
  auto &cell = PidPdfMan::at(t,p,c,b,momid);
  return cell.valid ? cell.perr[5] : std::numeric_limits<double>::quiet_NaN();
}
double PidPdfMan::getFitValue(int t,int p,int c,int b,int momid, int paramid)
{
  paramid = paramid%pidfunc::kNparamGauss;
  auto &cell = PidPdfMan::at(t,p,c,b,momid);
  return cell.valid ? cell.p[paramid] : std::numeric_limits<double>::quiet_NaN();
}
double PidPdfMan::getFitErrValue(int t,int p,int c,int b,int momid, int paramid)
{
  auto &cell = PidPdfMan::at(t,p,c,b,momid);
  return cell.valid ? cell.perr[paramid] : std::numeric_limits<double>::quiet_NaN();
}

bool PidPdfMan::SetGEPoint(int t_int, int p_int, int c_int, int b_int, int momid)
{
    auto type_e = static_cast<pidlikeli::DType>(t_int);
    auto pid_e  = static_cast<pidlikeli::Pid>(p_int);
    auto chg_e  = static_cast<pidlikeli::Chg>(c_int);
    auto be_e   = static_cast<pidlikeli::BE>(b_int);
    //    auto& b = m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
    const auto& cell = PidPdfMan::at(t_int, p_int, c_int, b_int, momid);
    if (!cell.valid) return false;
    std::map<CorrGraph::Graph, int> param_map = {
        {CorrGraph::Graph::MeanM2,      pidfunc::kParamIdM2},
        {CorrGraph::Graph::MeandEdx,    pidfunc::kParamIddEdx},
        {CorrGraph::Graph::SigM2,       pidfunc::kParamIdSigM2},
        {CorrGraph::Graph::SigdEdx,     pidfunc::kParamIdSigdEdx},
        {CorrGraph::Graph::RotAngle,    pidfunc::kParamIdRotAngle},
        {CorrGraph::Graph::Yield,       pidfunc::kParamIdYield},
        {CorrGraph::Graph::MeanM2Mom,   pidfunc::kParamIdM2},
        {CorrGraph::Graph::MeandEdxMom, pidfunc::kParamIddEdx},
        {CorrGraph::Graph::SigM2Mom,    pidfunc::kParamIdSigM2},
        {CorrGraph::Graph::SigdEdxMom,  pidfunc::kParamIdSigdEdx},
        {CorrGraph::Graph::RotAngleMom, pidfunc::kParamIdRotAngle},
        {CorrGraph::Graph::YieldMom,    pidfunc::kParamIdYield}
    };

    CorrGraph& cgi = m_data->at_ge_init(type_e, pid_e, chg_e, be_e);
    double mom = pidlikeli::BinToMom(momid);
    double beta = pidlikeli::MomToBetaPid(mom,p_int);

    using Gr = CorrGraph::Graph;
    for (const auto& pair : param_map) {
        CorrGraph::Graph graph_type = pair.first;
        int param_id = pair.second;
        auto cgigraph = cgi.graphs[pidlikeli::scast(graph_type)];
        int n = cgigraph->GetN();
        TString graph_name = CorrGraph::GraphNames[pidlikeli::scast(graph_type)];        
        double p_val = cell.p[param_id];
        double p_err = cell.perr[param_id];
        double x_value = graph_name.EndsWith("Mom") ? mom : beta;
        double x_err = graph_name.EndsWith("Mom") ? pidlikeli::kDP / 2.0 : 0.01;        
        cgigraph->SetPoint(n, x_value, p_val);
        cgigraph->SetPointError(n, x_err, p_err);
    }
    
    // good
    double mm2 = getFitMuM2(t_int, p_int, c_int, b_int, momid);
    double mdedx = getFitMudEdx(t_int, p_int, c_int, b_int, momid);
    double emm2 = getFitErrMuM2(t_int, p_int, c_int, b_int, momid);    
    double emdedx = getFitErrMudEdx(t_int, p_int, c_int, b_int, momid);


    std::cout << " mm2,emm2,cut: " << mm2 << "," << emm2 << "," << cutgm2
	      << " mdedx,emdedx,cut: " << mdedx << "," << emdedx << "," << cutgdedx << std::endl;
    CorrGraph& cg = m_data->at_ge(type_e, pid_e, chg_e, be_e);
    double cm2min = pidlikeli::cutm2min[p_int];
    double cm2max = pidlikeli::cutm2max[p_int];
    double cdemin = pidlikeli::cutdedxmin[p_int];
    double cdemax = pidlikeli::cutdedxmax[p_int];
    if(p_int==kAllParticles){
      cm2min = pidlikeli::cutm2min[kKaon];
      cm2max = pidlikeli::cutm2max[kKaon];
      cdemin = pidlikeli::cutdedxmin[kKaon];
      cdemax = pidlikeli::cutdedxmax[kKaon];      
    }
    bool m2cut
      = ( (c_int==0&&mm2>cm2min&&mm2<cm2max)||(c_int==1&&mm2>-cm2max&&mm2<-cm2min) )
      && (cutgm2>emm2);
    bool dedxcut
      =  (mdedx>cdemin&&mdedx<cdemax)
      && (cutgdedx>emdedx);
    
    auto& b = m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
    if(m2cut){
      b.is_goodM2 = true;
    }
    if(dedxcut){
      b.is_gooddEdx = true;
    }
    if(m2cut && dedxcut){
      b.is_goodM2dEdx = true;
      std::cout << "[GOOD] mm2,emm2,cut: " << mm2 << "," << emm2 << "," << cutgm2
		<< " mdedx,emdedx,cut: " << mdedx << "," << emdedx << "," << cutgdedx << std::endl;    
      for (size_t i = 0; i < static_cast<size_t>(Gr::COUNT); ++i) {
	auto& cggraph = cg.graphs[i];
	int n = cggraph->GetN();
	double p = getFitValue(t_int, p_int, c_int, b_int, momid, i);
	double perr = getFitErrValue(t_int, p_int, c_int, b_int, momid, i);
	TString graph_name = CorrGraph::GraphNames[i];
	std::cout << "p,error: " << p << "," << perr << std::endl;
	double x_value = graph_name.EndsWith("Mom") ? mom : beta;
	double x_err = graph_name.EndsWith("Mom") ? kDP : 0.01;
	//if(perr<piddata::cutgraph[i]){
	cggraph->SetPoint(n,x_value,p);
	cggraph->SetPointError(n,x_err,perr);
	  //}
      }
    }
    b.valid = true;
    return true;
}

bool PidPdfMan::FitGEPoint(int t_int, int p_int, int c_int, int b_int)
{
  auto type_e = static_cast<pidlikeli::DType>(t_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(p_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(c_int);
  auto be_e   = static_cast<pidlikeli::BE>(b_int);

  CorrGraph& cg = m_data->at_ge(type_e, pid_e, chg_e, be_e);  
  CorrFunc& cf = m_data->at_f_good(type_e, pid_e, chg_e, be_e);
  
  using Gr = CorrGraph::Graph;
  using Fn = CorrFunc::Func;
    
  for (size_t i = 0; i < static_cast<size_t>(Fn::COUNT); ++i) {
    auto graph_type = static_cast<Gr>(i);
    auto func_type = static_cast<Fn>(i);
    TF1* func_to_fit = cf.functions[pidlikeli::scast(func_type)];
    auto& g = cg.graphs[pidlikeli::scast(graph_type)];
    TString graph_name = CorrGraph::GraphNames[i];
    if(graph_name.EndsWith("Mom")) continue;
    if (func_to_fit){
      if(g->GetN()>0){
	g->Fit(func_to_fit, "R");
      } else {
	std::cout << "Info: Graph for " << func_to_fit->GetName() << " is empty. Skipping fit." << std::endl;
      }
    } else {
      std::cerr << "Error in FitGEPoint: The fit function (TF1) at index "
		<< i << " is null." << std::endl;
    }
  }
  return true;
}

bool PidPdfMan::SetYieldGEPoint(int t_int, int p_int, int c_int, int b_int, int momid)
{
    auto type_e = static_cast<pidlikeli::DType>(t_int);
    auto pid_e  = static_cast<pidlikeli::Pid>(p_int);
    auto chg_e  = static_cast<pidlikeli::Chg>(c_int);
    auto be_e   = static_cast<pidlikeli::BE>(b_int);

    double mm2 = getFitValue(t_int, p_int, c_int, b_int, momid, pidfunc::kParamIdM2);
    double mdedx = getFitValue(t_int, p_int, c_int, b_int, momid, pidfunc::kParamIddEdx);
    double emm2 = getFitErrValue(t_int, p_int, c_int, b_int, momid, pidfunc::kParamIdSigM2);
    double emdedx = getFitErrValue(t_int, p_int, c_int, b_int, momid, pidfunc::kParamIdSigdEdx);

    double cm2min = pidlikeli::cutm2min[p_int];
    double cm2max = pidlikeli::cutm2max[p_int];
    bool m2cut = ((c_int == 0 && mm2 > cm2min && mm2 < cm2max) || (c_int == 1 && mm2 > -cm2max && mm2 < -cm2min))
                 && (cutgm2 > emm2);
    bool dedxcut = (mdedx > pidlikeli::cutdedxmin[p_int] && mdedx < pidlikeli::cutdedxmax[p_int])
                   && (cutgdedx > emdedx);

    auto& b = m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
    
    bool is_good_fit = false;
    bool is_special_particle = ( (p_int == pidlikeli::kElectron)
				 || (p_int==pidlikeli::kDeutron&&c_int==pidlikeli::kPlus) );
    
    if ( scaleEleDeu&&is_special_particle ){
      is_good_fit = true;
    } else {
      is_good_fit = m2cut && dedxcut;
    }
    
    if ( is_good_fit
	 || ( extrapHighP&&t_int==pidlikeli::kTypeForPrior&&p_int==pidlikeli::kProton&&c_int==pidlikeli::kPlus ) ) {
      b.is_goodM2dEdx = true;
	
      CorrGraph& cg = m_data->at_ge(type_e, pid_e, chg_e, be_e);
        
      double yield_val = getFitValue(t_int, p_int, c_int, b_int, momid, pidfunc::kParamIdYield);
      double yield_err = getFitErrValue(t_int, p_int, c_int, b_int, momid, pidfunc::kParamIdYield);

      double mom = pidlikeli::BinToMom(momid);
      double mom_err = pidlikeli::kDP / 2.0;
        
      auto graph_yield_mom = cg.graphs[pidlikeli::scast(CorrGraph::Graph::YieldMom)];
      int n_yield_mom = graph_yield_mom->GetN();
      graph_yield_mom->SetPoint(n_yield_mom, mom, yield_val);
      graph_yield_mom->SetPointError(n_yield_mom, mom_err, yield_err);

      return true;
    }

    return false;
}

std::array<double, pidfunc::kNparamGauss> PidPdfMan::GetCalcParameters(int t, int pid, int c, int b, int m) const
{
  std::array<double, pidfunc::kNparamGauss> calcparam;
  
  double mom = pidlikeli::BinToMom(m);
  double massvalue = pidlikeli::mass[pid];
  double m2 = massvalue * massvalue;
  if(c==1) m2 = (-1.)*m2;
  double dedx = pidlikeli::dedx_pid(mom, pid);
  double m2sigma = pidfunc::CalcSigM2(pid, mom);
  double dedxsigma = pidfunc::CalcSigdEdx(pid, mom);
  double rotangle = 0.;
  double yield = InterpYield(t, pid, c, b, m);

  calcparam = {m2, dedx, m2sigma, dedxsigma, rotangle, yield};
 
  return calcparam;
}

std::array<double, pidfunc::kNparamGauss> 
PidPdfMan::GetCalcParametersRef(int pid, int momid) const
{
  const pidlikeli::ParticleFitConfig* current_config = nullptr;
  for (const auto& cfg : pidlikeli::gFitConfigs) {
    if (static_cast<int>(cfg.pid) == pid) {
      current_config = &cfg;
      break;
    }
  }

  if (current_config) {
    return GetCalcParametersRef(*current_config, momid);
  } else {
    std::cerr <<"Error: No fit configuration found for PID"<< pid << std::endl;
    return {};
  }
}

std::array<double, pidfunc::kNparamGauss> 
PidPdfMan::GetCalcParametersRef(const pidlikeli::ParticleFitConfig& config, int momid) const
{
  int pid = static_cast<int>(config.pid);
  int chg = static_cast<int>(config.source_chg);
    
  double momentum = pidlikeli::BinToMom(momid);

  //m2 calc
  double m2 = pidlikeli::mass[pid] * pidlikeli::mass[pid];
  if (chg == pidlikeli::kMinus) {
    m2 *= -1.0;
  }

  //dedx calc
  double dedx = pidlikeli::dedx_pid(momentum, pid);

  // sigma calc
  int sigma_pid = static_cast<int>(config.sigma_source_pid);
  double sigma_m2 = pidfunc::CalcSigM2(momentum, sigma_pid);
  double sigma_dedx = pidfunc::CalcSigdEdx(momentum, sigma_pid);

  // fix rotangle
  double rot_angle = 0.0;
    
  // calc yield
  int type_idx = static_cast<int>(config.source_type);
  int be_idx = static_cast<int>(config.source_be);
  double yield = InterpYield(type_idx, pid, chg, be_idx, momid);
  if (std::isnan(yield)) {
    yield = 10.0; // default value
  }

  return {m2, dedx, sigma_m2, sigma_dedx, rot_angle, yield};
}

double PidPdfMan::InterpYield(int t_int, int p_int, int c_int, int b_int, int m) const
{
  double yield = 0.;
  auto type_e = static_cast<pidlikeli::DType>(t_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(p_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(c_int);
  auto be_e   = static_cast<pidlikeli::BE>(b_int);
  const auto &ge = m_data->at_ge(type_e, pid_e, chg_e, be_e).graphs[pidlikeli::scast(CorrGraph::Graph::YieldMom)];
  double mom = pidlikeli::BinToMom(m);

  if ( extrapHighP && t_int == pidlikeli::kTypeForPrior
       && p_int == pidlikeli::kProton && c_int == pidlikeli::kPlus) {
    const double mom_hth = pidlikeli::extrapolateHighMomProton;// GeV/c
    const double mom_lth = pidlikeli::extrapolateLowMomProton;// GeV/c    

    if (mom <= mom_lth) {
      double sum_of_yields = 0.0;
      int num_points_for_avg = 0;
      const int n_avg_points = 5;
      int start_mom_bin = pidlikeli::MomToBin(mom_lth) + (n_avg_points - 1);
      if (start_mom_bin < 0) start_mom_bin = 0;

      for (int i = 0; i < n_avg_points; ++i) {
	int current_bin = start_mom_bin - i;
	if (current_bin <= 0) continue;

	//if (CheckGoodFitM2dEdx(t_int, p_int, c_int, b_int, current_bin)) {
	const auto& cell = PidPdfMan::at(t_int, p_int, c_int, b_int, current_bin);
	sum_of_yields += cell.p[pidfunc::kParamIdYield];
	num_points_for_avg++;
	//getchar();
	//}
      }
      if (num_points_for_avg > 0) {	
	double val = sum_of_yields / num_points_for_avg;
	return val;
      }      
    }
    if (mom >= mom_hth) {
      double sum_of_yields = 0.0;
      int num_points_for_avg = 0;
      const int n_avg_points = 5;
      
      int start_mom_bin = pidlikeli::MomToBin(mom_hth) - (n_avg_points - 1);
      if (start_mom_bin < 0) start_mom_bin = 0;

      for (int i = 0; i < n_avg_points; ++i) {
	int current_bin = start_mom_bin + i;
	if (current_bin >= kNmom) continue;

	//if (CheckGoodFitM2dEdx(t_int, p_int, c_int, b_int, current_bin)) {
	const auto& cell = PidPdfMan::at(t_int, p_int, c_int, b_int, current_bin);
	sum_of_yields += cell.p[pidfunc::kParamIdYield];
	num_points_for_avg++;
	//getchar();
	//}
      }
      if (num_points_for_avg > 0) {	
	double val = sum_of_yields / num_points_for_avg;
	//const auto& cell_th = PidPdfMan::at(t_int, p_int, c_int, b_int, pidlikeli::MomToBin(mom_th));	
	//val = cell_th.p[pidfunc::kParamIdYield];
	// 	  << " charge: " << c_int <<"return yield: " << val << std::endl;		
	return val;
      }
    }
  }
  

  if (ge->GetN() == 0) {
    // throw std::invalid_argument("ge is empty in InterpYield");
    return std::numeric_limits<double>::quiet_NaN();
  }
    
  int N = ge->GetN();
  std::vector<std::pair<double,double>> pts;
  pts.reserve(N);
  for(int i=0;i<N;++i){
    double x, y;
    ge->GetPoint(i, x, y);
    pts.emplace_back(x,y);
  }
  std::sort(pts.begin(), pts.end(),
	    [](auto &a, auto &b){ return a.first < b.first; });
  
  if (mom <= pts.front().first)  return pts.front().second;
  if (mom >= pts.back().first)   return pts.back().second;
  
  int lo = 0, hi = N-2;
  while(lo <= hi){
    int mid = (lo+hi)/2;
    if (pts[mid].first < mom) {
      if (pts[mid+1].first >= mom) {
	lo = mid;
	break;
      }
      lo = mid+1;
    }
    else {
      hi = mid-1;
    }
  }
  int i0 = lo;
  double x0 = pts[i0].first,   y0 = pts[i0].second;
  double x1n= pts[i0+1].first, y1n= pts[i0+1].second;
  double slp = (mom - x0)/(x1n - x0);
  yield = y0 + slp*(y1n - y0);
  
  return yield;
}

const FitResultPars&
PidPdfMan::GetStoredFitParamsOfRefPdf(const pidlikeli::ParticleFitConfig& config, int momid) const
{
  // get index from source config
  int type_idx = static_cast<int>(config.source_type);
  int pid_idx  = static_cast<int>(config.pid);
  int chg_idx  = static_cast<int>(config.source_chg);
  int be_idx   = static_cast<int>(config.source_be);

  return PidPdfMan::at(type_idx, pid_idx, chg_idx, be_idx, momid);
}

Long64_t PidPdfMan::WriteToRootfile(TFile* fout)
{
  if(!fout || !fout->IsOpen()) return -1;
  fout->cd();

  //  TTree  
  constexpr int be = 0;
  Long64_t nWritten =0;
  for(int it_int=0; it_int<kNtype; it_int++){
    for(int ip_int=0; ip_int<kNpidAll; ip_int++){
      for(int ic_int=0; ic_int<kNchg; ic_int++){
	for(int ibe_int=0; ibe_int<kNbe; ibe_int++){
	  auto type_e = static_cast<pidlikeli::DType>(it_int);
	  auto pid_e  = static_cast<pidlikeli::Pid>(ip_int);
	  auto chg_e  = static_cast<pidlikeli::Chg>(ic_int);
	  auto be_e   = static_cast<pidlikeli::BE>(ibe_int);	    
	  for(int im=0; im<kNmom; im++){
	    auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, im);
	    auto f2d = b.get_func<TF2>(PdfFuncType::F2D);
	    if (!f2d) continue;
	    long long id
	      = fac_t *(it_int+1)
	      + fac_p * ip_int
	      + fac_c * ic_int
	      + fac_b * ibe_int
	      + fac_m * im;
	    if(f2d){
	    }
	    auto f2dc = b.get_func<TF2>(PdfFuncType::F2DCORR);
	    if (f2dc) {
	      const char* f2dc_basename
		= PdfFuncName[static_cast<size_t>(PdfFuncType::F2DCORR)];
	      TString name = TString(f2dc_basename) + "%lld";
	      WriteClonedObject(f2dc, fout, name.Data(), id);
	      ++nWritten;
	    }
	    auto fd2d = b.get_func<TF2>(PdfFuncType::FD2D);	    
	    if (fd2d) {
	      const char* fd2d_basename
		= PdfFuncName[static_cast<size_t>(PdfFuncType::FD2D)];
	      TString name = TString(fd2d_basename) + "%lld";	      
	      WriteClonedObject(fd2d, fout, name.Data(), id);
	      ++nWritten;
	    }
	    auto fx = b.get_func<TF1>(PdfFuncType::FX_PROJ);    
	    if(fx){
	      const char* fx_basename
		= PdfFuncName[static_cast<size_t>(PdfFuncType::FX_PROJ)];
	      TString name = TString(fx_basename) + "%lld";	      	      
	      WriteClonedObject(fx, fout, name.Data(), id);
	      ++nWritten;	      
	    }
	    auto fy = b.get_func<TF1>(PdfFuncType::FY_PROJ);    	    
	    if(fy){
	      const char* fy_basename
		= PdfFuncName[static_cast<size_t>(PdfFuncType::FX_PROJ)];
	      TString name = TString(fy_basename) + "%lld";	      	      	  
	      WriteClonedObject(fy, fout,name.Data(), id);
	      ++nWritten;	      
	    }
	    auto pdf = b.get_func<TF2>(PdfFuncType::PDF);
	    if(pdf){
	      const char* pdf_basename
		= PdfFuncName[static_cast<size_t>(PdfFuncType::PDF)];
	      TString name = TString(pdf_basename) + "%lld";	      
	      WriteClonedObject(pdf, fout, name.Data(), id);
	      ++nWritten;	      
	    }
	  }
	  long long idfunc
	    = fac_t *(it_int+1)
	    + fac_p * ip_int
	    + fac_c * ic_int
	    + fac_b * ibe_int;
	  const auto& cg = m_data->at_ge(type_e, pid_e, chg_e, be_e);
	  using Gr = CorrGraph::Graph;
	  for (size_t i = 0; i < static_cast<size_t>(Gr::COUNT); ++i) {
	    const auto& graph_to_write = cg.graphs[i];
	    if(graph_to_write->GetN()>0){
	      TString graph_enum_name = CorrGraph::GraphNames[i];
	      TString final_name = Form("%s%lld",graph_enum_name.Data(),idfunc);
	      WriteClonedObject(graph_to_write, fout, final_name.Data());
	      ++nWritten;		
	    }
	  }
	  const auto& cf = m_data->at_f_good(type_e, pid_e, chg_e, be_e);
	  using Fn = CorrFunc::Func;
	  for (size_t i = 0; i < static_cast<size_t>(Fn::COUNT); ++i) {
	    const auto& func_to_write = cf.functions[i];
	    if(func_to_write){
	      TString func_enum_name = CorrFunc::FuncNames[i];
	      TString final_name = Form("%s%lld",func_enum_name.Data(),idfunc);
	      WriteClonedObject(func_to_write, fout, final_name.Data());
	      ++nWritten;		
	    }
	  }			
	}
      }
    }
  }
  if (m_data && m_data->m_h5_priors) {
    m_data->m_h5_priors->Write("priprob");
    nWritten++;
  }
  std::cout << "nWritten: " << nWritten << std::endl;
  return nWritten;
    //  }
}

Bool_t PidPdfMan::WritePDFToRootfile()
{
  if(!m_mainfile) return false;
  Long64_t n = WriteToRootfile(m_mainfile.get());
  if(n>0) return true;
  if(n==0 or n==-1) return false;
}

void PidPdfMan::DrawResultsToPdf()
{
  const std::string& output_filename = "";
  DrawResultsToPdf(output_filename);
}

void PidPdfMan::DrawResultsToPdf(const std::string& output_filename)
{
  FileStat_t file_info;
  if (gSystem->GetPathInfo(m_file_name.Data(), file_info) != 0) {
    std::cerr << "Error in DrawResultsToPdf: Input file does not exist: "
	      << m_file_name.Data() << std::endl;
    return;
  }

#if DrawHistToPdf
    std::cout << "Drawing results to " << output_filename << " ..." << std::endl;
    
    gROOT->SetBatch(true);
    //gStyle->SetOptStat(1);
    gStyle->SetPalette(kRainBow);

    TCanvas dummy_canvas("dummy", "dummy", 1, 1);
    dummy_canvas.Print((output_filename + "[").c_str());

    const int plots_per_canvas = 20;
    const int mom_step = kNmom / plots_per_canvas;
    static int canvas_counter = 0;

    for (int it = 0; it < kNtype; ++it) {
      for (int ic = 0; ic < kNchg; ++ic) {
	for (int ib = 0; ib < kNbe; ++ib) {
	  for (int ip = 0; ip < kNpidAll; ++ip) {	    
	    auto type_e = static_cast<pidlikeli::DType>(it);
	    auto pid_e  = static_cast<pidlikeli::Pid>(ip);
	    auto chg_e  = static_cast<pidlikeli::Chg>(ic);
	    auto be_e   = static_cast<pidlikeli::BE>(ib);
	    {
	      TString c_name = Form("c_hist_%d", canvas_counter);
	      TString c_title = Form("HISTs: Type=%d, PID=%d, Chg=%d, BE=%d", it, ip, ic, ib);
	      TCanvas* c_th2 = new TCanvas(c_name, c_title, 1600, 1600);	    
	      c_th2->Divide(5,5);

	      std::vector<TH2D*> hist_clones;
	      for (int i_plot = 0; i_plot < plots_per_canvas; ++i_plot) {
		int im = i_plot * mom_step;
		if (im >= kNmom) continue;		      
		TH2D* hist_clone = GetHist(it,ip,ic,ib,im);
		hist_clones.push_back(hist_clone);
	      }
	      bool drawn = false;
	      for (size_t i_plot = 0; i_plot < hist_clones.size(); ++i_plot) {
		c_th2->cd(i_plot + 1)->SetLogz(1);
		//c_th2->cd(i_plot + 1);			
		TH2D* hist = hist_clones[i_plot];
		if (hist) {
		  int im = i_plot * mom_step;
		  hist->SetTitle(Form("Index: %d, Mom = %.2f GeV/c", im, pidlikeli::BinToMom(im)));
		  hist->Draw("COLZ");
		  //hist->Draw("SURF2");		
		  drawn = true;
		}	      
	      }
	      if(drawn){
		c_th2->Print(output_filename.c_str());
	      }
	      for (auto h : hist_clones) { if (h) delete h; }
	      delete c_th2;
	    }
	    
	    {	    
	      TCanvas* c_fd2d = new TCanvas(Form("c_fd2d_%d", canvas_counter++), Form("fd2d: type=%d,pid=%d,chg=%d,be=%d", it, ip, ic, ib), 1600, 1600);
	      c_fd2d->Divide(5,5);
	      bool fd2d_drawn = false;
	      for (int i_plot = 0; i_plot < plots_per_canvas; ++i_plot) {
		int momid = i_plot * mom_step;
		if (momid >= kNmom) continue;
		c_fd2d->cd(i_plot + 1);
		const auto& bin = m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
		if (auto f = bin.get_func<TF2>(PdfFuncType::FD2D)) {
		  f->SetTitle(Form("%s, Mom = %.2f GeV/c", f->GetName(),pidlikeli::BinToMom(momid)));
		  f->Draw("SURF2");
		  fd2d_drawn = true;
		}	      
	      }
	      if (fd2d_drawn) c_fd2d->Print(output_filename.c_str());
	      delete c_fd2d;


	      TCanvas* c_f2dc = new TCanvas(Form("c_f2dc_%d", canvas_counter++), Form("f2dc: type=%d,pid=%d,chg=%d,be=%d", it, ip, ic, ib), 1600, 1600);
	      c_f2dc->Divide(5, 5);
	      bool f2dc_drawn = false;
	      for (int i_plot = 0; i_plot < plots_per_canvas; ++i_plot) {
		int momid = i_plot * mom_step;
		if (momid >= kNmom) continue;
		c_f2dc->cd(i_plot + 1);
		const auto& bin = m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
		if (auto f = bin.get_func<TF2>(PdfFuncType::F2DCORR)) {
		  f->SetTitle(Form("%s, Mom = %.2f GeV/c", f->GetName(),pidlikeli::BinToMom(momid)));
		  f->Draw("SURF2");
		  f2dc_drawn = true;
		}
	      }
	      if (f2dc_drawn) c_f2dc->Print(output_filename.c_str());
	      delete c_f2dc;
	    }

	    {
	      TString c_name = Form("c_graph_%d", canvas_counter++);
	      TString c_title = Form("Graphs: type=%d, pid=%d, chg=%d, be=%d", it, ip, ic, ib);
	      TCanvas* c_graph = new TCanvas(c_name, c_title, 1600, 1600);
	      c_graph->Divide(6,4);
	      const auto& corr_graph = m_data->at_ge(type_e, pid_e, chg_e, be_e);
	      const auto& corr_graphinit = m_data->at_ge_init(type_e, pid_e, chg_e, be_e);	    
	      bool graph_drawn = false;
	      for (size_t ig = 0; ig < static_cast<size_t>(CorrGraph::Graph::COUNT); ++ig) {
		c_graph->cd(ig + 1);
		TGraphErrors* graph = corr_graph.graphs[ig];
		if (graph && graph->GetN() > 0) {
		  graph->SetTitle(CorrGraph::GraphNames[ig]);
		  graph->SetMarkerSize(0.7);			    
		  graph->SetMarkerStyle(20);
		  graph->Draw("AP");
		  graph_drawn = true;
		}
		c_graph->cd(ig + 1 + 12);	      
		TGraphErrors* graphinit = corr_graphinit.graphs[ig];
		if (graphinit && graphinit->GetN() > 0) {
		  graphinit->SetTitle(CorrGraph::GraphNames[ig]+"-init");
		  graphinit->SetMarkerSize(0.7);			    
		  graphinit->SetMarkerStyle(20);
		  graphinit->Draw("AP");
		}	      
	      }

	      if (graph_drawn) {
		c_graph->cd(0);
		TPaveText* title_box = new TPaveText(0.1, 0.96, 0.9, 0.99, "NDC");
		title_box->SetFillColor(kWhite);
		title_box->AddText(c_title);
		title_box->Draw();
		c_graph->Print(output_filename.c_str());
	      }
	      //if (graph_drawn) c_graph->Print(output_filename.c_str());
	      delete c_graph;
	    }
	  }
	}
      }
    }
    {
      int it = pidlikeli::kTypeGen;
      int ip = pidlikeli::kAllParticles;
      const int plots_per_page_1d = 10;
      const int mom_step_1d = kNmom / plots_per_page_1d;
      std::vector<TF1*> components_m2, components_dedx;      
      for (int ic = 0; ic < kNchg; ++ic) {
	for (int ib = 0; ib < kNbe; ++ib) {
	  if(ib!=0) continue;
	  TString c_name = Form("c_1dfit_%d_%d_%d", ic, ib, canvas_counter++);
	  TString c_title = Form("1D Projection Fits: Type=%d, Chg=%d, BE=%d", it, ic, ib);
	  TCanvas* c_1dfit = new TCanvas(c_name, c_title, 1600, 1200);
	  c_1dfit->Divide(5,4);

	  bool drawn_anything = false;
	  for (int i_plot = 0; i_plot < plots_per_page_1d; ++i_plot) {
	    int im = i_plot * mom_step_1d;
	    if (im >= kNmom) continue;
	    TH2* h2 = GetHist(it, ip, ic, ib, im);
	    if (!h2 || h2->GetEntries() < 1) {
	      if (h2) delete h2;
	      continue;
	    }
	    
	    const auto& cell = m_data->at_bin(
					      static_cast<pidlikeli::DType>(it),
					      static_cast<pidlikeli::Pid>(ip),
					      static_cast<pidlikeli::Chg>(ic),
					      static_cast<pidlikeli::BE>(ib),
					      im);
	    
	    //if (!cell.valid) continue;
	    TF1* func_m2 = cell.get_func<TF1>(PdfFuncType::F1D_M2);
	    TF1* func_dedx = cell.get_func<TF1>(PdfFuncType::F1D_DEDX);
	    
	    if (!func_m2 || !func_dedx){
	      delete h2;
	      continue;
	    }
	    std::cout << "--- Drawing mom=" << im << " ---" << std::endl;
	    std::cout << "m2 Function Parameters:" << std::endl;
	    for (int i=0; i < func_m2->GetNpar(); ++i) {
	      if (std::abs(func_m2->GetParameter(i)) > 1e-9) {
		std::cout << "  p" << i << " (" << func_m2->GetParName(i) << ") = " 
			  << func_m2->GetParameter(i) << std::endl;
	      }
	    }
          
	    TH1D* h_m2 = h2->ProjectionX(Form("h_m2_draw_%d_%d_%d", ic, ib, im));
	    TH1D* h_dedx = h2->ProjectionY(Form("h_dedx_draw_%d_%d_%d", ic, ib, im));
	    h_m2->SetTitle(Form("m2 (mom=%.2f);m^{2} (GeV^{2}/c^{4});Entries", pidlikeli::BinToMom(im)));
	    h_dedx->SetTitle(Form("dE/dx (mom=%.2f);dE/dx;Entries", pidlikeli::BinToMom(im)));
	    drawn_anything = true;

	    bool is_good_fit = CheckGoodFitM2dEdx(it, ip, ic, ib, im);
	    if (is_good_fit) {
	      func_m2->SetLineColor(kRed);
	      func_m2->SetLineStyle(1);
	      func_dedx->SetLineColor(kRed);
	      func_dedx->SetLineStyle(1);
	    } else {
	      func_m2->SetLineColor(kGray + 1);
	      //func_m2->SetLineStyle(2);
	      func_m2->SetLineStyle(1);
	      func_dedx->SetLineColor(kGray + 1);
	      func_dedx->SetLineStyle(2);
	    }	   

	    c_1dfit->cd(i_plot+1);
	    std::cout << "i_plot: " << i_plot << std::endl;
	    //gPad->SetLogy();
	    if(h_m2){
	      //h_m2->Draw("E");
	      h_m2->DrawClone();
	      func_m2->Draw("SAME");
	    }
	    c_1dfit->cd(i_plot + 11);
	    //gPad->SetLogy();
	    if(h_dedx){
	      //h_dedx->Draw("E");
	      h_dedx->DrawClone();
	      func_dedx->Draw("SAME");
	    }

	    for (int i_pid = 0; i_pid < pidlikeli::kNpid; ++i_pid) {
	      int offset = i_pid * pidfunc::kNparamGauss;
	      //M2
	      TF1* comp_m2 = new TF1(Form("comp_m2_%d_%d", im, i_pid),
				     pidfunc::Single1DGaussForM2,
				     h_m2->GetXaxis()->GetXmin(),
				     h_m2->GetXaxis()->GetXmax(),
				     pidfunc::kNparamGauss);
	      for(int i=0; i<pidfunc::kNparamGauss; i++)
		{
		  comp_m2->SetParameter(i, func_m2->GetParameter(offset+i));
		}
	      
	      comp_m2->SetLineColor(kGray + 2);
	      comp_m2->SetLineStyle(2);
	      c_1dfit->cd(i_plot + 1);
	      comp_m2->Draw("SAME");
	      components_m2.push_back(comp_m2);

	      // dEdx
	      TF1* comp_dedx = new TF1(Form("comp_dedx_%d_%d", im, i_pid),
				       pidfunc::Single1DGaussFordEdx,
				       h_dedx->GetXaxis()->GetXmin(),
				       h_dedx->GetXaxis()->GetXmax(),
				       pidfunc::kNparamGauss);
	      for(int i=0; i<pidfunc::kNparamGauss; i++){
		comp_dedx->SetParameter(i, func_dedx->GetParameter(offset+i));
	      }
	
	      comp_dedx->SetLineColor(kGray + 2);
	      comp_dedx->SetLineStyle(2);
	      c_1dfit->cd(i_plot + 11);
	      comp_dedx->Draw("SAME");
	      components_dedx.push_back(comp_dedx);
	    }
	    delete h2;
	    delete h_m2;
	    delete h_dedx;
	  }
	  if (drawn_anything) {
	    c_1dfit->Print(output_filename.c_str());
	  }
	  delete c_1dfit;
	}
      }      
    }
    {
      int it = pidlikeli::kTypeForPrior;
      int ib = 0;
      int colors[] = {kRed, kBlue, kGreen + 2, kMagenta, kOrange + 7};

      for (int ic = 0; ic < kNchg; ++ic) {
	TString charge_str = (ic == pidlikeli::kPlus) ? "Positive" : "Negative";

	{
	  TString c_name = Form("c_yield_summary_chg%d", ic);
	  TString c_title = Form("Yield vs Momentum (%s Charge)", charge_str.Data());
	  TCanvas* c_yield_summary = new TCanvas(c_name, c_title, 900, 700);
	  gPad->SetGridx(); gPad->SetGridy();

	  //gPad->SetLogy();
	  TMultiGraph* mg_yield = new TMultiGraph();
	  mg_yield->SetTitle(c_title + ";Momentum (GeV/c);Yield");

	  TLegend* leg_yield = new TLegend(0.7, 0.65, 0.9, 0.9);
	  leg_yield->SetHeader("Particle Type");

	  for (int ip = 0; ip < kNpid; ++ip) {
	    const auto& corr_graph = m_data->at_ge(
						   static_cast<pidlikeli::DType>(it), static_cast<pidlikeli::Pid>(ip),
						   static_cast<pidlikeli::Chg>(ic), static_cast<pidlikeli::BE>(ib));
                
	    TGraphErrors* g = corr_graph.graphs[pidlikeli::scast(CorrGraph::Graph::YieldMom)];
	    if (g && g->GetN() > 0) {
	      g->SetLineColor(colors[ip]);
	      g->SetMarkerColor(colors[ip]);
	      g->SetMarkerStyle(20 + ip);
	      mg_yield->Add(g, "P");
	      leg_yield->AddEntry(g, pidlikeli::plist[ip].Data(), "lp");
	    }
	  }
            
	  if (mg_yield->GetListOfGraphs() && mg_yield->GetListOfGraphs()->GetEntries() > 0) {
	    mg_yield->Draw("A");
	    leg_yield->Draw();
	  }
	  c_yield_summary->Print(output_filename.c_str());
	  // delete c_yield_summary;
	  // delete mg_yield;
	  // delete leg_yield;
	}

	{
	  TString c_name = Form("c_prior_summary_chg%d", ic);
	  TString c_title = Form("Prior Probability vs Momentum (%s Charge)", charge_str.Data());
	  TCanvas* c_prior_summary = new TCanvas(c_name, c_title, 900, 700);
	  gPad->SetGridx(); gPad->SetGridy();

	  TMultiGraph* mg_prior = new TMultiGraph();
	  mg_prior->SetTitle(c_title + ";Momentum (GeV/c);Prior Probability");
	  mg_prior->SetMinimum(0.0);
	  mg_prior->SetMaximum(1.05);

	  TLegend* leg_prior = new TLegend(0.7, 0.65, 0.9, 0.9);
	  leg_prior->SetHeader("Particle Type");
            
	  std::vector<TGraph*> temp_graphs; 

	  for (int ip = 0; ip < kNpid; ++ip) {
	    TGraph* g_prior = new TGraph();
	    temp_graphs.push_back(g_prior);

	    for (int im = 0; im < kNmom; ++im) {
	      int bin_coords[] = {pidlikeli::kTypeForPrior + 1, ip + 1, ic + 1, ib + 1, im + 1};
	      double prior = m_data->m_h5_priors->GetBinContent(m_data->m_h5_priors->GetBin(bin_coords));
	      if (prior > 0) {
		double mom = pidlikeli::BinToMom(im);
		g_prior->SetPoint(g_prior->GetN(), mom, prior);
	      }
	    }
                
	    if (g_prior->GetN() > 0) {
	      g_prior->SetLineColor(colors[ip]);
	      g_prior->SetLineWidth(2);
	      mg_prior->Add(g_prior, "L");
	      leg_prior->AddEntry(g_prior, pidlikeli::plist[ip].Data(), "l");
	    }
	  }
            
	  if (mg_prior->GetListOfGraphs() && mg_prior->GetListOfGraphs()->GetEntries() > 0) {
	    mg_prior->Draw("A");
	    leg_prior->Draw();
	  }
	  c_prior_summary->Print(output_filename.c_str());
            
	  // delete c_prior_summary;
	  // delete mg_prior;
	  // delete leg_prior;
	  // for (auto g : temp_graphs) delete g;
	}
      }
    }
    dummy_canvas.Print((output_filename + "]").c_str());
    std::cout << "Drawing results finished." << std::endl;
#endif    
}
