// -*- C++ -*-
#include "PidPdfMan.hh"
#include "PidCommon.hh"
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

#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <vector>
#include <iterator>
#include "TFitResult.h"

#include "FuncName.hh"
#include "DeleteUtility.hh"
#include <TSystem.h>

#include "PidCommon.hh"
#include "PidData.hh"
#include "PidPdfMan.hh"

#define DrawHistToPdf 1

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

  // for setting fit config
  const int kNparamYield = 1;
  const int n = pidfunc::kNparamGauss-kNparamYield;
  const std::array<ParticleFitConfig, n> gFitConfigs = {{
      // { PID, name, Parameter configuration source (type, chg, be), sigma source (pid) }
      { pidlikeli::Pid::Pi, "Pi", pidlikeli::DType::K0, pidlikeli::Chg::Minus, pidlikeli::BE::All, pidlikeli::Pid::Pi },
      { pidlikeli::Pid::K,  "K",  pidlikeli::DType::Km, pidlikeli::Chg::Minus, pidlikeli::BE::RsK, pidlikeli::Pid::K  },
      { pidlikeli::Pid::P,  "P",  pidlikeli::DType::Lmd,pidlikeli::Chg::Plus,  pidlikeli::BE::All, pidlikeli::Pid::P  },
      { pidlikeli::Pid::D,  "D",  pidlikeli::DType::General, pidlikeli::Chg::Plus, pidlikeli::BE::All, pidlikeli::Pid::P },
      { pidlikeli::Pid::E,  "E",  pidlikeli::DType::General, pidlikeli::Chg::Minus, pidlikeli::BE::All, pidlikeli::Pid::Pi}
    }};
  
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

  std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;
  
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
    std::cout << "debug " << __FILE__ << __LINE__
	      << " " << __func__ << " hname: " << hname << std::endl;
    
    if(hname[0] != 'h') continue;
    Long64_t id = std::strtoll(hname.Data()+1, nullptr, 10);
    //TH2D* h = key->ReadObject<TH2D>();
    auto h = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(key->ReadObject<TH2D>()->Clone()));
    if(!h) continue;
    auto clone = static_cast<TH2D*>(h->Clone());    
    h->SetDirectory(nullptr);
    m_histMap[id] = std::move(h);
    
    std::cout << "debug " << __FILE__ << __LINE__
	      << " " << __func__ << std::endl;    
    
  }
  return !m_histMap.empty();
}

bool PidPdfMan::Run(){
  std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;      
  // GetParam(kPion);
  // GetParam(kProton);
  // MakeGoodPdf(kPion);
  // MakeGoodPdf(kProton);
  // GetParamKaonUsingOtherParticles();
  // MakeGoodPdfKm();
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
  } else if(runmode==runmode==pidlikeli::kDCarbonKP){ // Kp on Carbon
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
  std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__ << std::endl;
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
	  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;	  
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
  std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__ << std::endl;
  for(int it_int = 0; it_int < kNtype; ++it_int){
    for(int chg_int = 0; chg_int < kNchg; ++chg_int){
      for(int ibe_int = 0; ibe_int < kNbe; ++ibe_int){
	if( !( (it_int==pidlikeli::kTypeGen&&chg_int==kMinus&&ibe_int==0)
	       || (it_int==pidlikeli::kTypeKm&&chg_int==kMinus) )) continue;
	//if(!(it_int==pidlikeli::kTypeGen&&chg_int==kMinus&&ibe_int==0)) continue;
	if( it_int==pidlikeli::kTypeGen ){
	  pid_int=kKaon;
	  fit_pid_int=kAllParticles;
	}
	if( it_int==pidlikeli::kTypeKm ){
	  pid_int=kKaon;
	  fit_pid_int=kKaon;
	}
	for(int m = 0; m < kNmom; ++m){
	    std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__
		      << " it,ip,ic,ib,im: " << it_int << "," << pid_int << "," << chg_int << "," << ibe_int << "," << m << std::endl;
	  
	  if( !CreateFitFunction(pidlikeli::kTypeKm,pid_int,chg_int,ibe_int,m) ) continue;
	  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;
	  if(!SetInitFitParsDoubleGauss(pidlikeli::kTypeKm,pid_int,chg_int,ibe_int,m)) continue;
	  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;	  
	  if(it_int==pidlikeli::kTypeGen){
	    if(!ExecFitDouble2DGaussAllParticles(pidlikeli::kTypeKm,pid_int,chg_int,ibe_int,m) ) continue;
	    std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__
		      << " it,ip,ic,ib,im: " << it_int << "," << pid_int << "," << chg_int << "," << ibe_int << "," << m << std::endl;
	  }
	  if(it_int==pidlikeli::kTypeKm){
	    if(!ExecFitDouble2DGauss(pidlikeli::kTypeKm,pid_int,chg_int,ibe_int,m) ) continue;	    
	  }
	  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;
	  if( !StoreFitResultKaon(pidlikeli::kTypeKm,pid_int,chg_int,ibe_int,m) ) continue;
	  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;	  	  
	  if( !SetGEPoint(pidlikeli::kTypeKm,pid_int,chg_int,ibe_int,m) ) continue;
	  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;	  	  	  
	}
	if( !FitGEPoint(pidlikeli::kTypeKm,pid_int,chg_int,ibe_int) ) continue;	    
      }
    }
  }
  return true;
}
bool PidPdfMan::CreateFitFunction(int it_int, int pid_int, int chg_int, int ibe_int, int m)
{
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;
  auto type_e = static_cast<pidlikeli::DType>(it_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(pid_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(chg_int);
  auto be_e   = static_cast<pidlikeli::BE>(ibe_int);          
  auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, m);
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ 
	    << " t,p,c,b,m: " << it_int+1 << "," << pid_int << ","
	    << chg_int << "," << ibe_int << "," << m << std::endl;
  if(b.functions_created) return true;
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;
  auto h2 = GetHist(it_int, pid_int, chg_int, ibe_int, m);
  if(!h2){
    if(it_int==pidlikeli::kTypeKm&&pid_int==pidlikeli::kKaon&&chg_int==pidlikeli::kMinus&&ibe_int==0){
      h2 = new TH2D(Form("h%dp%dc%dbe%dm%d",it_int,pid_int,chg_int,ibe_int,m),
		    Form("h%dp%dc%dbe%dm%d",it_int,pid_int,chg_int,ibe_int,m),
		    pidlikeli::nbinm2,pidlikeli::minm2,pidlikeli::maxm2,
		    pidlikeli::nbindedx,pidlikeli::mindedx,pidlikeli::maxdedx );
    } else {
      std::cout<<"debug "<<__FILE__<<":"<<__LINE__<<" "<<__func__<<std::endl;
      delete h2;
      return false;
    }
  }
  
  h2->SetDirectory(nullptr);
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;	        
  b.fx = new TF1(Form("fx%dp%dc%dbe%dm%d",it_int,pid_int,chg_int,ibe_int,m),
		 pidfunc::RotGauss2DProjXFit,
		 h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax(), pidfunc::kNparamGauss);
  b.fy = new TF1(Form("fy%dp%dc%dbe%dm%d",it_int,pid_int,chg_int,ibe_int,m),
		 pidfunc::RotGauss2DProjYFit,
		 h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax(), pidfunc::kNparamGauss);
  b.f2d = new TF2(Form("f2dt%dp%dc%dbe%dm%d",it_int,pid_int,chg_int,ibe_int,m),
		  pidfunc::RotGauss2D,
		  h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax(),
		  h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax(), pidfunc::kNparamGauss);
  b.f2dc = new TF2(Form("f2dct%dp%dc%dbe%dm%d",it_int,pid_int,chg_int,ibe_int,m),
		      pidfunc::RotGauss2D,
		      h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax(),
		      h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax(), pidfunc::kNparamGauss);
  b.fd2d = new TF2(Form("fd2dt%dp%dc%dbe%dm%d",it_int,pid_int,chg_int,ibe_int,m),
		   pidfunc::RotDoubleGauss2D,
		   h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax(),
		   h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax(), pidfunc::kNparamDGauss);  
  b.functions_created = true;
  delete h2;

  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__
	    << " functions created "<< std::endl;	   
  return true;
}

bool PidPdfMan::CreateFitFunctionFive2DGauss(int chg_int, int ibe_int, int m)
{
  int it_int = pidlikeli::kTypeGen;
  int ip_int = pidlikeli::kAllParticles;
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;
  auto type_e = static_cast<pidlikeli::DType>(it_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(ip_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(chg_int);
  auto be_e   = static_cast<pidlikeli::BE>(ibe_int);          
  auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, m);
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ 
	    << " t,p,c,b,m: " << it_int+1 << "," << ip_int << ","
	    << chg_int << "," << ibe_int << "," << m << std::endl;
  if(b.five2dgauss_created) return true;
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;
  auto h2 = GetHist(it_int, ip_int, chg_int, ibe_int, m);
  if(!h2){
      return false;
  }  
  h2->SetDirectory(nullptr);
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;	  
  b.ff2d = new TF2(Form("ff2dt%dp%dc%dbe%dm%d",it_int,ip_int,chg_int,ibe_int,m),
		   pidfunc::RotFiveGauss2D,
		   h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax(),
		   h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax(),
		   pidfunc::kNparamGauss*pidlikeli::kNpid);
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
  h1x->Fit(b.fx, "Q");	  
  std::cout << " h1y->GetEntries(): " << h1y->GetEntries()
	    << " h1y->GetMean(): "    << h1y->GetMean()
	    << " h1y->GetStdDev(): "  << h1y->GetStdDev() << std::endl;	  
  h1y->Fit(b.fy, "Q");
  //std::cout << "Fit 1D hist with 1D gaussian " << std::endl;
  std::cout << "debug " << __FILE__ << __LINE__ << " ip: " << pid_int
	    << " sigma1 " << b.fx->GetParameter(2) << " sigma2 " << b.fy->GetParameter(3) << std::endl;
  delete h1x;
  delete h1y;  
  delete h2;
  
  return true;
}

bool PidPdfMan::ExecFitDouble2DGauss(int it_int,int pid_int,int chg_int,int ibe_int,int m)
{
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;	       
  auto h = GetHist(it_int, pid_int, chg_int, ibe_int, m);
  if (!h || h->GetEntries() < 1) return false;
  h->SetDirectory(nullptr);

  auto type_e = static_cast<pidlikeli::DType>(it_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(pid_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(chg_int);
  auto be_e   = static_cast<pidlikeli::BE>(ibe_int);          
  auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, m);
  if(!b.fd2d) return false;
  //h->Fit(b.fd2d,"Q");
  h->Fit(b.fd2d,"NR");
  delete h;
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__
	    << " Execute fitting double 2D gaussian "<< std::endl;	       
  return true;
}

bool PidPdfMan::ExecFitDouble2DGaussAllParticles(int it_int,int pid_int,int chg_int,int ibe_int,int m)
{
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;
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
  if(!b.fd2d) return false;
  //h->Fit(b.fd2d,"Q");
  h->Fit(b.fd2d,"NR");
  delete h;
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__
	    << " Execute fitting double 2D gaussian "<< std::endl;	       
  return true;
}

bool PidPdfMan::SetInitFitParsDoubleGauss(int type, int ip, int chg, int be, int momid)
{
  std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;  
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
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__
	    << " Fix with pion parameters"<< std::endl;
  TF2*fkpi = bkm.fd2d;
  TF2*fpi = bpm.f2d;
  if(!fkpi||!fpi) return false;
  for(int i=0; i<pidfunc::kNparamGauss-1; i++){
    std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__
	      << "FixParam" << i << " " << bpm.f2d->GetParameter(i) << std::endl;
    fkpi->FixParameter(i+pidfunc::kNparamGauss, fpi->GetParameter(i));
    //fkpi->SetParameter(i+6, fpi->GetParameter(i));
  }
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__
	    << " Fix with pion parameters"<< std::endl;	       
  // Fix/Set param
  int lastparamid=5;
  fkpi->SetParameter(lastparamid+pidfunc::kNparamGauss, fpi->GetParameter(lastparamid));
  if(!SetInitFitPars2D(type,kKaon,chg,be,momid,h,bkm.fd2d)) return false;
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__
	    << " Set initial fitting parameters"<< std::endl;
  delete h;    
  return true;  
}

bool PidPdfMan::SetInitFitParsProj(int type, int pid, int chg, int be, int momid)
{
  std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;

  auto h = GetHist(type, pid, chg, be, momid);
  if (!h || h->GetEntries() < 1) return false;
  auto type_e = static_cast<pidlikeli::DType>(type);
  auto pid_e = static_cast<pidlikeli::Pid>(pid);
  auto chg_e = static_cast<pidlikeli::Chg>(chg);
  auto be_e = static_cast<pidlikeli::BE>(be);
  auto &b= m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
  for(int i=0; i<2; i++){
    TF1* f; 
    if(i==0) f = b.fx;
    if(i==1) f = b.fy;
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

    std::cout << "debug " << __FILE__ << __LINE__
	      << " " << __func__ << std::endl;
  
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
      //f->SetParLimits(4,0,TMath::Pi()/2.);
      //f->SetParLimits(5,tot*0.5,2*tot);
      //f->FixParameter(0,meanx);
      // f->FixParameter(0,meanx);
      // f->FixParameter(2,sigx);                                  
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

TFitResultPtr PidPdfMan::ExecFitFive2DGauss(int chg, int be, int momid)
{
  std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__ << std::endl;  
  int it_int = pidlikeli::kTypeGen;
  int ip_int = pidlikeli::kAllParticles;
  auto h = GetHist(it_int, ip_int, chg, be, momid);
  if (!h || h->GetEntries() < 1){
    return false;
  }
  std::cout<<"debug "<<__FILE__<<__LINE__<<" "<<__func__<<std::endl;
  
  auto type_e = static_cast<pidlikeli::DType>(it_int);
  auto pid_e = static_cast<pidlikeli::Pid>(ip_int);
  auto chg_e = static_cast<pidlikeli::Chg>(chg);
  auto be_e = static_cast<pidlikeli::BE>(be);
  auto &b= m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
  if(!b.five2dgauss_created) return false;
  TF2 * fitmodel = b.ff2d;

  std::cout<<"debug "<<__FILE__<<__LINE__<<" "<<__func__<<std::endl;
  
  for(const auto& config: gFitConfigs) {
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
    
    std::cout<<"debug "<<__FILE__<<__LINE__<<" "<<__func__<<std::endl;
    
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
    
    std::cout<<"debug "<<__FILE__<<__LINE__<<" "<<__func__<<std::endl;
    
    if(chg==pidlikeli::kMinus) {
      params[pidfunc::kParamIdM2] = -1.0 * std::abs(params[pidfunc::kParamIdM2]);
    } else {
      params[pidfunc::kParamIdM2] = std::abs(params[pidfunc::kParamIdM2]);
    }
    
    int n = pidfunc::kNparamGauss-kNparamYield;
    for(int iparam=0; iparam<n; iparam++){
      std::cout<<"debug "<<__FILE__<<__LINE__<<" "<<__func__
	       << " SetParameter " <<std::endl;	          
      fitmodel->FixParameter(param_offset+iparam,params[iparam]);
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
    if(yield_init<1.0) yield_init = 1.0;
    // Set yield parameter
    fitmodel->SetParameter(param_offset+n, yield_init); 
    fitmodel->SetParLimits(param_offset+n, 0.0, h->GetEntries() * 2);
  }
  std::cout<<"debug "<<__FILE__<<__LINE__<<" "<<__func__<<std::endl;	    
  TFitResultPtr fit_result = h->Fit(fitmodel, "LS");
  //h->Fit(fitmodel, "LS");
  fit_result->Print();
  //TFitResultPtr fit_result = nullptr;
  delete h;  
  return fit_result;  
}

bool PidPdfMan::SetInitFitPars2D(int type, int pid, int chg, int be, int momid, TH2D*h, TF2*f)
{
  std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;
  //auto h = GetHist(type, pid, chg, be, momid);
  
  // if(!h) return false;    
  // if (h->GetEntries()<1){
  //   std::cout << "No entry in " << h->GetName() << std::endl;
  //   return false;
  // }
  // set initial fitting parameter condition
  if(!f) return false;  
  //f->SetParNames("x0","y0","sigx","sigy","theta","N");
  double meanx     = -999.;
  double meany     = -999.;
  double sigx      = -999.;
  double sigy      = -999.;
  double rotangle  = -999.;
  double count     = -999.;
  double mom = pidlikeli::BinToMom(momid);

  std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__ << std::endl;
  
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
    return b.valid && b.f2d;
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
	    if(check){
	      std::cout << "debug " << __FILE__ << __LINE__ << " "
			<< b.f2d->GetName() << " created " << std::endl;
	    } else {
	      std::cout << "debug " << __FILE__ << __LINE__
			<< " not created " << std::endl;
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
  std::cout<<"debug "<<__FILE__<<":"<<__LINE__<<" "<<__func__<<std::endl;
  auto type_e = static_cast<pidlikeli::DType>(it_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(ip_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(ic_int);
  auto be_e   = static_cast<pidlikeli::BE>(ibe_int);
  const auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
  if(!b.fd2d) return false;
  TF2*f = b.fd2d;
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
    TF2* f2d = b.f2d;
    for(int i=0; i<pidfunc::kNparamGauss; i++){
      f2d->SetParameter(i,f->GetParameter(i));
    }
  }
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__
	    << " store fitting result "<< std::endl;	         
  return true;
}

bool
PidPdfMan::StoreFitResultProj(int it,int pid,int ic,int be,int momid)
{
  std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;
  auto type_e = static_cast<pidlikeli::DType>(it);
  auto pid_e  = static_cast<pidlikeli::Pid>(pid);
  auto chg_e  = static_cast<pidlikeli::Chg>(ic);
  auto be_e   = static_cast<pidlikeli::BE>(be);	      
  const auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
  if(!b.fx||!b.fy) return false;
  TF1*fx = b.fx;
  TF1*fy = b.fy;  
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
  b.f2d->SetParameters(cell.p);
  b.f2d->SetParErrors(cell.perr);
  
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

bool PidPdfMan::CheckGoodFitM2dEdx(int t_int, int p_int, int c_int, int be_int, int momid){
  auto type_e = static_cast<pidlikeli::DType>(t_int);
  auto pid_e  = static_cast<pidlikeli::Pid>(p_int);
  auto chg_e  = static_cast<pidlikeli::Chg>(c_int);
  auto be_e   = static_cast<pidlikeli::BE>(be_int);
  auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
  return b.is_goodM2dEdx;
}

bool PidPdfMan::MakeGoodPdf(int pid_int)
{
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;
  for(int it=0; it < kNtype; it++){
    for(int ic=0; ic < kNchg; ic++){
      for(int ibe=0; ibe < kNbe; ibe++){
        for(int im=0; im < kNmom; im++){
	  auto type_e = static_cast<pidlikeli::DType>(it);
	  auto pid_e  = static_cast<pidlikeli::Pid>(pid_int);
	  auto chg_e  = static_cast<pidlikeli::Chg>(ic);
	  auto be_e   = static_cast<pidlikeli::BE>(ibe);
	  auto &b = m_data->at_bin(type_e, pid_e, chg_e, be_e, im);
	  //if(!b.f2d||!b.f2dc) continue;
	  if(!b.f2d) continue;
	  double mom = pidlikeli::BinToMom(im);
          if(CheckGoodFitM2dEdx(it, pid_int, ic, ibe, im)&&mom<0.7){
            b.f2dc->SetParameters(b.f2d->GetParameters());
	    for(int i=0; i<pidfunc::kNparamGauss; i++){
	      std::cout << "f2d->GetParameter(" << i << "): " << b.f2d->GetParameter(i) << std::endl;
	      std::cout << "f2dc->GetParameter(" << i << "): " << b.f2dc->GetParameter(i) << std::endl;	      
	    }
          } else {
            std::array<double,pidfunc::kNparamGauss> params = GetCalcParameters(it,pid_int,ic,ibe,im);
            b.f2dc->SetParameters(params.data());
          }
        }
      }
    }
  }
  return true;
}

bool PidPdfMan::MakeGoodPdfKm()
{
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl;
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
	  if(!b.f2d) continue;
	  //if(!b.f2d||!b.f2dc) continue;
          if(CheckGoodFitM2dEdx(it, pid_int, ic, ibe, im)){
            b.f2dc->SetParameters(b.f2d->GetParameters());	
	    std::cout << "debug " << __FILE__ << ":" << __LINE__
		      << " " << __func__ << std::endl;	    		      
          } else {
            std::array<double, pidfunc::kNparamGauss> params = GetCalcParameters(it, pid_int, ic, ibe, im);
            b.f2dc->SetParameters(params.data());
	    std::cout << "debug " << __FILE__ << ":" << __LINE__
		      << " " << __func__ << std::endl;	    
          }
        }
      }
    }
  }
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__
	    << " Make good PDF Kminus" << std::endl;	      
  return true;
}

void PidPdfMan::StoreGoodPdfAll()
{

  std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__
	    << " Generating PDF ..." << std::endl;  
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
	    TF2* source_func = bin.f2dc;	    
	    if (!source_func) {
	      continue;
	    }
	    if (!bin.pdf) {
	      bin.pdf = static_cast<TF2*>(source_func->Clone());
	      bin.pdf->SetName(Form("pdf_%s", source_func->GetName()));
	    }
	    for (int i = 0; i < pidfunc::kParamIdYield; ++i) {
	      bin.pdf->SetParameter(i, source_func->GetParameter(i));
	    }
	    bin.pdf->SetParameter(pidfunc::kParamIdYield, 1.0);
	  }
	}
      }
    }
  }
  std::cout << "Finished storing good PDFs." << std::endl;
}


double PidPdfMan::GetYield(int it_int, int ip_int, int ic_int, int ib_int, int im)
{
  std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;
  
  double yield = 0.;
  if(CheckGoodFitM2dEdx(it_int, ip_int, ic_int, ib_int, im)){
    yield = getFitYield(it_int, ip_int, ic_int, ib_int, im);
  } else {
    yield = InterpYield(it_int, ip_int, ic_int, ib_int, im);
  }
  std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;  
  return yield;
}

double PidPdfMan::GetTotalYield(int it_int, int ic_int, int ib_int, int im)
{
  std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;
  
  double totalyield = 0.;
  for(int ip=0; ip<kNpid; ip++){
    totalyield += GetYield(it_int, ip, ic_int, ib_int, im);
  }
  return totalyield;
}

void PidPdfMan::FillPriorHist()
{
  std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;  
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
  std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__ << std::endl;
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
  std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;
  
  if( !m_data || !m_data->m_h5_priors) {
    std::cerr << "Error: Prior hist (THn) is already initialized. " << std::endl;
    return false;
  }
  std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;  
  double priorprob = GetPriorProb(it_int,ip_int,ic_int,ib_int,im);
  if(GetTotalYield(it_int,ic_int,ib_int,im)>100){
    int bin_coords[5] = {it_int+1, ip_int+1, ic_int+1, ib_int+1, im+1};
    long long global_bin = m_data->m_h5_priors->GetBin(bin_coords);
    std::cout << "debug " << __FILE__ << __LINE__
	      << " " << __func__ << std::endl;    
    m_data->m_h5_priors->SetBinContent(global_bin, priorprob);
  }

  std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;
  
  return true;
}

std::vector<double> PidPdfMan::CalculatePriorProb(const TFitResultPtr& fit_result) const
{
  std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__ << std::endl;  
  if (!fit_result || fit_result->Status() != 0) {
    if (fit_result) {
      std::cerr << "Warning: Fit did not converge properly (status=" << fit_result->Status() << "). Cannot calculate ratios." << std::endl;
    } else {
      std::cerr << "Null fitting result" << std::endl;
    }
    return {}; // null vector
  }
  std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__ << std::endl;
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
  std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__ << std::endl;    
  double total_yield = std::accumulate(yields.begin(), yields.end(), 0.0);

  if (total_yield < 1e-9) {
    std::cerr << "Warning: Total yield is nearly zero. Cannot calculate ratios." << std::endl;
    return {};
  }
  std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__ << std::endl;
  std::vector<double> ratios;
  ratios.reserve(pidlikeli::kNpid);
  for (double yield : yields) {
    double priprob = yield / total_yield;
    ratios.push_back(priprob);
    std::cout << " Prior Prob: " << priprob << std::endl;
  }
  std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__ << std::endl;  
  return ratios;
}

void PidPdfMan::StorePriorProbs()
{
  std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__ << std::endl;  
  std::cout << "Starting generation of prior probabilities..." << std::endl;
  if (!m_data || !m_data->m_h5_priors) {
    if(!m_data) std::cout << 2 << std::endl;
    if(!m_data->m_h5_priors) std::cout << 3 << std::endl;
    std::cerr << "Error: Cannot generate priors, PidPdfMan not ready." << std::endl;
    return;
  }

  for (int ic = 0; ic < pidlikeli::kNchg; ++ic) {
    for (int ib = 0; ib < pidlikeli::kNbe; ++ib) {
      if(ib!=0) continue;
      for (int im = 0; im < pidlikeli::kNmom; ++im) {
	std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__ <<
	  " Processing bin: chg=" << ic << ", be=" << ib << ", mom=" << im << std::endl;
	// // Fitting
	if(!CreateFitFunctionFive2DGauss(ic,ib,im)) continue;  
	TFitResultPtr fit_res = ExecFitFive2DGauss(ic, ib, im);
	std::vector<double> priors = CalculatePriorProb(fit_res);
	
	if (!priors.empty()) {
	  int type_idx = pidlikeli::kTypeGen; 
	  for (int ip = 0; ip < kNpid; ++ip) {
	    int bin_coords[5] = {type_idx + 1, ip + 1, ic + 1, ib + 1, im + 1};
	    long long global_bin = m_data->m_h5_priors->GetBin(bin_coords);
	    m_data->m_h5_priors->SetBinContent(global_bin, priors[ip]);
	  }
	} else {
	  std::cout<<"debug "<<__FILE__<<__LINE__<<" "<<__func__
		   << " Empty vector " <<std::endl;
	}
      }
    }
  }
  std::cout << "Finished generation of prior probabilities." << std::endl;
}

void
PidPdfMan::GetPriorProb()
{
  std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__ << std::endl;
  
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
  std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;
  
  double nInt[kNpid];
  double nTotal = 0.;
  double priorprob = 0.;
  // pid
  for(int pid_int=0; pid_int<kNpid; ++pid_int){
    nInt[pid_int] = GetYield(it_int,pid_int,chg_int,ibe_int,im_int);
    nTotal += nInt[pid_int];
  }


  std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;
  
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
    std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;

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
    std::cout << "debug " << __FILE__ << __LINE__
	    << " " << __func__ << std::endl;
    
    auto type_e = static_cast<pidlikeli::DType>(t_int);
    auto pid_e  = static_cast<pidlikeli::Pid>(p_int);
    auto chg_e  = static_cast<pidlikeli::Chg>(c_int);
    auto be_e   = static_cast<pidlikeli::BE>(b_int);

    CorrGraph& cgi = m_data->at_ge_init(type_e, pid_e, chg_e, be_e);
    double mom = pidlikeli::BinToMom(momid);
    double beta = pidlikeli::MomToBetaPid(mom,p_int);

    using Gr = CorrGraph::Graph;
    for (size_t i = 0; i < static_cast<size_t>(Gr::COUNT); ++i) {
      auto graph_type = static_cast<Gr>(i);
      auto cgigraph = cgi.graphs[pidlikeli::scast(graph_type)];
      int n = cgigraph->GetN();
      TString graph_name = CorrGraph::GraphNames[i];      
      double p = getFitValue(t_int, p_int, c_int, b_int, momid, i);
      double perr = getFitErrValue(t_int, p_int, c_int, b_int, momid, i);
      double x_value = graph_name.EndsWith("Mom") ? mom : beta;
      double x_err = graph_name.EndsWith("Mom") ? kDP : 0.01;
      cgigraph->SetPoint(n,x_value,p);
      cgigraph->SetPointError(n,x_err,perr);
     // std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__
     // 	<< " " << cgigraph->GetName() << " SetPoint " << std::endl;      
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
    
    // std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__
    // 	      << " M2:" << mm2 << " cut[" << cm2min << "," << cm2max << "] "
    // 	      << " dEdx:" << mdedx << " cut[" << cdemin << "," << cdemax << "]"
    // 	      << std::endl;
    
    auto& b = m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
    if(m2cut){
      b.is_goodM2 = true;
      std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__
		<< " GoodM2" << std::endl;      
    }
    if(dedxcut){
      b.is_gooddEdx = true;
      std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__
		<< " GoodDedx" << std::endl;            
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
	// std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__
	// 	  << " " << cggraph->GetName() << " SetPointGood" << std::endl;	
      }
    }
    b.valid = true;
    std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__
	      << " Set GraphErrors points " << std::endl;	        
    return true;
}

bool PidPdfMan::FitGEPoint(int t_int, int p_int, int c_int, int b_int)
{
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__ << std::endl; 
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
  std::cout << "debug " << __FILE__ << ":" << __LINE__ << " " << __func__
	    << " Fit GraphErrors points " << std::endl;	           
  return true;
}

std::array<double, pidfunc::kNparamGauss> PidPdfMan::GetCalcParameters(int t, int pid, int c, int b, int m) const
{
  std::cout << "debug " << __FILE__ << __LINE__ << " " << __func__ << std::endl;  
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
  const ParticleFitConfig* current_config = nullptr;
  for (const auto& cfg : gFitConfigs) {
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
PidPdfMan::GetCalcParametersRef(const ParticleFitConfig& config, int momid) const
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
PidPdfMan::GetStoredFitParamsOfRefPdf(const ParticleFitConfig& config, int momid) const
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
	    if (!b.f2d) continue;
	    long long id
	      = fac_t *(it_int+1)
	      + fac_p * ip_int
	      + fac_c * ic_int
	      + fac_b * ibe_int
	      + fac_m * im;
	    if(b.f2d){
	    }
	      /*
	      fout->WriteObject(b.f2d.get(), Form("f2d%lld", id));
	      ++nWritten;
	      */
	    if (b.f2dc) {
	      WriteClonedObject(b.f2dc, fout, "f2dc%lld", id);
	      ++nWritten;
	    }
	    if (b.fd2d) {
	      WriteClonedObject(b.fd2d, fout, "fd2d%lld", id);
	      ++nWritten;
	    }
	    if(b.fx){
	      WriteClonedObject(b.fx, fout, "fx%lld", id);
	      ++nWritten;	      
	    }
	    if(b.fy){
	      WriteClonedObject(b.fy, fout, "fy%lld", id);
	      ++nWritten;	      
	    }
	    if(b.pdf){
	      WriteClonedObject(b.pdf, fout, "pdf%lld", id);
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
      for (int ip = 0; ip < kNpidAll; ++ip) {
	for (int ic = 0; ic < kNchg; ++ic) {
	  for (int ib = 0; ib < kNbe; ++ib) {
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
		    
	    TCanvas* c_fd2d = new TCanvas(Form("c_fd2d_%d", canvas_counter++), Form("fd2d: type=%d,pid=%d,chg=%d,be=%d", it, ip, ic, ib), 1600, 1600);
	    c_fd2d->Divide(5,5);
	    bool fd2d_drawn = false;
	    for (int i_plot = 0; i_plot < plots_per_canvas; ++i_plot) {
	      int momid = i_plot * mom_step;
	      if (momid >= kNmom) continue;
	      c_fd2d->cd(i_plot + 1);
	      const auto& bin = m_data->at_bin(type_e, pid_e, chg_e, be_e, momid);
	      if (bin.fd2d) {
		bin.fd2d->SetTitle(Form("%s, Mom = %.2f GeV/c", bin.fd2d->GetName(),pidlikeli::BinToMom(momid)));
		bin.fd2d->Draw("SURF2");
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
	      if (bin.f2dc) {
		bin.f2dc->SetTitle(Form("%s, Mom = %.2f GeV/c", bin.f2dc->GetName(),pidlikeli::BinToMom(momid)));
		bin.f2dc->Draw("SURF2");
		f2dc_drawn = true;
		// TPaveText* pt = new TPaveText(0.05, 0.9, 0.25, 1.0, "blNDC");
		// pt->SetFillColor(0);
		// pt->SetBorderSize(0);
		// pt->SetTextAlign(12);
		// pt->AddText(Form("Index: %d", i_plot));
		// pt->Draw();
	      }
	    }
	    if (f2dc_drawn) c_f2dc->Print(output_filename.c_str());
	    delete c_f2dc;

		    
	    TCanvas* c_graph = new TCanvas(Form("c_graph_%d", canvas_counter++), Form("Graphs: type=%d,pid=%d,chg=%d,be=%d", it, ip, ic, ib), 1600, 1600);
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
	    if (graph_drawn) c_graph->Print(output_filename.c_str());
	    delete c_graph;	    
	  }
	}
      }
    }
    dummy_canvas.Print((output_filename + "]").c_str());
    std::cout << "Drawing results finished." << std::endl;
#endif    
}
