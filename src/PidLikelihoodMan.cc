// // -*- C++ -*-
// src/PidLikelihoodMan.cc

#include "PidLikelihoodMan.hh"
#include <iostream>
#include <TFile.h>
#include <numeric>
#include "PidCommon.hh"
#include "PidData.hh"

// === Singleton, Constructor, etc. ===
PidLikelihoodMan* PidLikelihoodMan::s_instance = nullptr;
PidLikelihoodMan& PidLikelihoodMan::GetInstance() 
{
  if (!s_instance) {
    s_instance = new PidLikelihoodMan();
    std::atexit(PidLikelihoodMan::Cleanup);
  }
  return *s_instance;
}

void PidLikelihoodMan::Cleanup() {
    if (s_instance) {
        s_instance->Finalize();
    }
}

bool PidLikelihoodMan::HasInstance() {
    return s_instance != nullptr;
}


PidLikelihoodMan::PidLikelihoodMan()
  : m_is_ready(false),
    m_file_name()
{  
}
PidLikelihoodMan::~PidLikelihoodMan()
{ 
}
void PidLikelihoodMan::Finalize() {
    if (!m_is_ready) return;
    m_is_ready = false;
    std::cout << "PidLikelihoodMan finalized." << std::endl;
    delete s_instance;
    s_instance = nullptr;
}

Bool_t PidLikelihoodMan::Initialize()
{
  if (m_is_ready) return true;
  
  // save current dir
  TDirectory* savedDir = gDirectory;
  TFile*      savedFile = gFile;

  // open pidpdf file
  TString file_name = m_file_name;

  TFile * file = TFile::Open(file_name,"READ");
  if(!file || file->IsZombie()){
    std::cerr << __func__ << " file open fail : "
	      << m_file_name << std::endl;
    return false;
  }

  m_prior_hist.reset(dynamic_cast<THnD*>(file->Get("priprob")));
  if (!m_prior_hist) {
    std::cerr << "Error: PriorProbabilities THnD not found in " << file_name << std::endl;
    file->Close();
    if (savedFile) savedFile->cd();
    else if (savedDir) savedDir->cd();
    return false;
  }

  for (int ip = 0; ip < pidlikeli::kNpid; ++ip) {  
    const auto& config = pidlikeli::gFitConfigs[ip];
        
    int ref_type_idx = static_cast<int>(config.source_type);
    int ref_chg_idx  = static_cast<int>(config.source_chg);
    int ref_be_idx   = static_cast<int>(config.source_be);
        
    for (int im = 0; im < pidlikeli::kNmom; ++im) {
      long long id = pidlikeli::fac_t * (ref_type_idx + 1)
	+ pidlikeli::fac_p * ip
	+ pidlikeli::fac_c * ref_chg_idx
	+ pidlikeli::fac_b * ref_be_idx
	+ pidlikeli::fac_m * im;            
      TString pdf_name = Form("pdf%lld", id);
      TF2* pdf_from_file = dynamic_cast<TF2*>(file->Get(pdf_name));
      if (pdf_from_file) {
	m_pdf_table[ip][im].pdf = std::make_unique<TF2>(*static_cast<TF2*>(pdf_from_file->Clone()));	
      }
    }
  }
  file->Close();
  if (savedFile) savedFile->cd();
  else if (savedDir) savedDir->cd();
  m_is_ready = true;
  std::cout << "PidLikelihoodMan initialized successfully." << std::endl;
  return true;
}

Bool_t PidLikelihoodMan::Initialize(const TString& file_name)
{
  m_file_name = file_name;
  return Initialize();
}

// === Posterior prob calculation ===
std::vector<double> PidLikelihoodMan::CalculatePosterior(int charge, double momentum, double m2, double dedx) const
{
  if ( (std::isnan(m2) || std::isnan(dedx)) ) return {};
  if( !m_is_ready ) return {};
  std::cout << " debug " << __FILE__ << " " << __LINE__ << " " << __func__ << std::endl;
  int mom_bin = pidlikeli::MomToBin(momentum);
  int chg_bin = pidlikeli::ChgToBin(charge);
  if( mom_bin<0||chg_bin<0 ) return {};
  std::cout << " debug " << __FILE__ << " " << __LINE__ << " " << __func__ << std::endl;
  std::vector<double> log_posteriors;
  log_posteriors.reserve(pidlikeli::kNpid);
  double total_posterior_sum = 0.0;

  for (int ip = 0; ip < pidlikeli::kNpid; ++ip) {
    int bin_coords[] =
      {pidlikeli::kTypeForPrior + 1, ip + 1, chg_bin + 1, pidlikeli::kBEForPrior + 1, mom_bin + 1};
    double prior = m_prior_hist->GetBinContent(m_prior_hist->GetBin(bin_coords)); 
    double log_prior = (prior > 1e-10) ? TMath::Log(prior) : -100.0;
    
    const auto& pdf_bin = m_pdf_table[ip][mom_bin];
    if (!pdf_bin.pdf) {
      log_posteriors.push_back(-1e10);
      continue;
    }
    TF2* pdf = pdf_bin.pdf.get();

    double params[pidfunc::kNparamGauss];
    pdf->GetParameters(params);

    double temp_params[pidfunc::kNparamGauss];
    std::copy(params, params + pidfunc::kNparamGauss, temp_params);
        
    temp_params[0] = std::abs(temp_params[0]);
    if (charge<1) temp_params[0] *= -1.0;
    temp_params[5] = 1.0;

    double xy_vals[] = {m2, dedx};
    double log_likelihood = pidfunc::LogRotGauss2D(xy_vals, temp_params);
    log_posteriors.push_back(log_likelihood + log_prior);
  }
  std::cout << " debug " << __FILE__ << " " << __LINE__ << " " << __func__ << std::endl;
  double max_log_post = -1e10;
  for (double lp : log_posteriors) {
    if (lp > max_log_post) max_log_post = lp;
  }

  if (max_log_post <= -1e9) {
        return std::vector<double>(pidlikeli::kNpid, 0.0);
  }
  std::cout << " debug " << __FILE__ << " " << __LINE__ << " " << __func__ << std::endl;
  double posterior_sum_exp = 0.0;
  for (double lp : log_posteriors) {
    posterior_sum_exp += TMath::Exp(lp - max_log_post);
  }
  
  if (total_posterior_sum > 1e-9) {
    for (double& post : log_posteriors) {
      post /= total_posterior_sum;
    }
  }

  std::vector<double> final_posteriors;
  final_posteriors.reserve(pidlikeli::kNpid);
  if (posterior_sum_exp > 1e-9) {
    for (double lp : log_posteriors) {
      final_posteriors.push_back(TMath::Exp(lp - max_log_post) / posterior_sum_exp);
    }
  } else {
    std::cout << " debug " << __FILE__ << " " << __LINE__ << " " << __func__ << std::endl;
    return std::vector<double>(pidlikeli::kNpid, 0.0);
  }
  std::cout << " debug " << __FILE__ << " " << __LINE__ << " " << __func__ << std::endl;
  return final_posteriors;
}

std::vector<PidResult> 
PidLikelihoodMan::GetRankedPid(int charge, double momentum, double m2, double dedx) const
{
  std::vector<double> posteriors = CalculatePosterior(charge, momentum, m2, dedx);

    if (posteriors.empty()) {
        return {};
    }

    std::vector<PidResult> results;
    results.reserve(pidlikeli::kNpid);
    
    for (int ip = 0; ip < pidlikeli::kNpid; ++ip) {
      results.push_back({static_cast<pidlikeli::Pid>(ip), posteriors[ip]});
    }
    
    std::sort(results.begin(), results.end(), 
              [](const PidResult& a, const PidResult& b) {
		return a.postprob > b.postprob;
              });

    return results;
}
