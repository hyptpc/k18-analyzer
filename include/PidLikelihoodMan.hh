// -*- C++ -*-
// // include/PidLikelihoodMan.hh

#ifndef PIDLIKELIHOODMAN_HH
#define PIDLIKELIHOODMAN_HH

#include <vector>
#include <string>
#include <memory>
#include "THn.h"
#include "TF2.h"
#include "PidCommon.hh"
#include "PidData.hh" // for ParticleFitConfig

struct PidResult {
  pidlikeli::Pid pid;
  double postprob;
};

struct LHPdfBin {
  std::unique_ptr<TF2> pdf;
};

using LHPdfTable = 
  std::array< // pid
    std::array<LHPdfBin, pidlikeli::kNmom>,
  pidlikeli::kNpid>;


class PidLikelihoodMan
{
public:
  static PidLikelihoodMan& GetInstance();
  ~PidLikelihoodMan();
  Bool_t Initialize();
  Bool_t Initialize(const TString& file_name);
  void SetFileName(const TString& file_name) { m_file_name = file_name; }  
  Bool_t IsReady() const { return m_is_ready; }
  static const TString& ClassName();  
  static bool HasInstance();
  std::vector<double> CalculatePosterior(int charge, double momentum, double m2, double dedx) const;
  std::vector<PidResult> GetRankedPid(int charge, double momentum, double m2, double dedx) const;
  void DrawResultsToPdf(const std::string& output_filename);
  void DrawResultsToPdf();    

private:
  PidLikelihoodMan(); // Singleton
  PidLikelihoodMan(const PidLikelihoodMan&) = delete;
  PidLikelihoodMan& operator=(const PidLikelihoodMan&) = delete;
  void Finalize();
  static void Cleanup();
  
  //member
  static PidLikelihoodMan* s_instance;
  bool m_is_ready = false;
  TString  m_file_name;      
  std::unique_ptr<THnD> m_prior_hist;
  //std::array<std::unique_ptr<TF2>, pidlikeli::kNpid> m_ref_pdfs;
  LHPdfTable m_pdf_table;
  
};


inline const TString&
PidLikelihoodMan::ClassName()
{
  static TString s_name("PidLikelihoodMan");
  return s_name;
}
#endif
