#if TKIN
#ifndef KINEMATICAL_FIT_HH
#define KINEMATICAL_FIT_HH

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TString.h>
#include <TMath.h>

#include <KinFitter/TKinFitter.h>
#include <KinFitter/TFitParticlePxPyPz.h>
#include <KinFitter/TFitConstraintM.h>
#include <KinFitter/TFitConstraintEp.h>
//#include <Math/ProbFuncMathCore.h>

class KinematicalFit
{
 public:
  static KinematicalFit&      GetInstance( void );
  static const std::string& ClassName( void );
  ~KinematicalFit( void );

 private:
  KinematicalFit( void );
  KinematicalFit( const KinematicalFit& );
  KinematicalFit& operator =( const KinematicalFit& );
  
  bool        m_is_ready;
  //  TString     target_name;
  TKinFitter* kinfitter;
  TLorentzVector TL_kfit[6];

 public:
  bool Initialize(bool yamaga=true);
  void Ldn(TVector3 l[6]);
  void Lpn(TVector3 l[6]);
  double get_chi2()      const { return kinfitter->getS(); } 
  double get_NDF()       const { return kinfitter->getNDF(); } 
  double get_status()    const { return kinfitter->getStatus(); } 
  //  double get_prob()      const { return ROOT::Math::chisquared_cdf_c(get_chi2(), get_NDF()); }
  double get_prob()      const { return TMath::Prob(get_chi2(), get_NDF()); }
  TLorentzVector get_lv(const int i) { return TL_kfit[i]; }
};
//______________________________________________________________________________
inline KinematicalFit&
KinematicalFit::GetInstance( void )
{
  static KinematicalFit g_instance;
  return g_instance;
}
//______________________________________________________________________________
inline const std::string&
KinematicalFit::ClassName( void )
{
  static std::string g_name("KinematicalFit");
  return g_name;
}
#endif
#endif
