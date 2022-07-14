//  Authors: Wooseung Jung

//GenKEK
#include "HypTPCFitter.hh"
#include "HypTPCField.hh"

//GenFit
#include <KalmanFitterRefTrack.h>
#include <KalmanFitter.h>
#include <KalmanFitterRefTrack.h>
#include <DAF.h>
#include <FieldManager.h>
#include <FitStatus.h>
#include <MaterialEffects.h>
#include <TGeoMaterialInterface.h>
#include <Track.h>

//k18-analyzer
#include <UserParamMan.hh>

#include <cstddef>
#include <iostream>

#define LogDEBUG(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogERROR(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogWARNING(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl

namespace{ const auto& gUser = UserParamMan::GetInstance(); }

ClassImp(HypTPCFitter)

genfit::AbsKalmanFitter* HypTPCFitter::_fitter = nullptr;
int HypTPCFitter::GenFitFitter = -1;
TGeoManager *HypTPCFitter::_tgeo_manager =  nullptr;

HypTPCFitter::HypTPCFitter(const std::string& tgeo_file_name, const bool m_is_const)
{
  _tgeo_manager = new TGeoManager("Geometry", "HypTPC geometry");
  _tgeo_manager = TGeoManager::Import(tgeo_file_name.data());

  genfit::FieldManager::getInstance()->init(new HypTPCField(m_is_const));
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  GenFitFitter = gUser.GetParameter("Fitter");
  unsigned int MinIter = gUser.GetParameter("nIteration",0);
  unsigned int MaxIter = gUser.GetParameter("nIteration",1);

  // init fitter
  if(GenFitFitter==0) _fitter = new genfit::KalmanFitterRefTrack();
  else if(GenFitFitter==1) _fitter = new genfit::KalmanFitter();
  else if(GenFitFitter==2) _fitter = new genfit::DAF(true);
  else if(GenFitFitter==3) _fitter = new genfit::DAF(false);
  else _fitter = new genfit::KalmanFitterRefTrack();
  _fitter -> setMaxIterations(MaxIter);
  _fitter -> setMinIterations(MinIter);

  genfit::Exception::quiet(true);
}

genfit::AbsKalmanFitter* HypTPCFitter::GetFitter(){
  return _fitter;
}
