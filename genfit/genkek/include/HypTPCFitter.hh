//  Authors: Wooseung Jung

#ifndef HYPTPCFITTER_HH
#define HYPTPCFITTER_HH

//GenKEK
#include "HypTPCField.hh"
#include "HypTPCTrack.hh"

//ROOT
#include <TGeoManager.h>

//GenFit
#include <FieldManager.h>
#include <FitStatus.h>
#include <Exception.h>
#include <AbsKalmanFitter.h>

//STL
#include <string>

class HypTPCFitter{

public:

  HypTPCFitter(const std::string& tgeo_file_name, const bool m_is_const=false);
  ~HypTPCFitter(){ delete _fitter; }
  /*!
   * Fitter control:
   * default & 0: KalmanFitterRefTrack()
   * 1: KalmanFitter()
   * 2: genfit::DAF(true) : DafRef
   * 3: genfit::DAF(false) : DafSimple
   */
  static genfit::AbsKalmanFitter* GetFitter();
  static int GenFitFitter;
  static genfit::AbsKalmanFitter* _fitter;

private:

  TGeoManager* _tgeo_manager;

  //ClassDef(HypTPCFitter, 1);

};  //class HypTPCFitter

#endif // HypTPCFitter_hh
