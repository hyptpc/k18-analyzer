#include "HSTrack.hh"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <TMinuit.h>

#include <std_ostream.hh>

#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"
#include "TrackHit.hh"
#include "K18TrackD2U.hh"

namespace {
  const auto &gGeom = DCGeomMan::GetInstance();

}

HSTrack::HSTrack(Double_t x, Double_t y,
		 Double_t u, Double_t v, Double_t p)
  :m_status(kInit),xInit(x), yInit(y), uInit(u), vInit(v),m_initial_momentum(p),
   m_polarity(0.),m_is_good(true)
{
  s_status[kInit]                = "Initialized";
  s_status[kPassed]              = "Passed";
  s_status[kExceedMaxPathLength] = "Exceed Max Path Length";
  s_status[kExceedMaxStep]       = "Exceed Max Step";
  s_status[kFailedSave]          = "Failed to Save";
  s_status[kFatal]               = "Fatal";

}

HSTrack::~HSTrack(){}

Bool_t HSTrack::Propagate() {
  m_status = kInit;
  
  if(m_initial_momentum <= 0){
    hddaq::cout << FUNC_NAME << " initial momentum must be exist"
		<< m_initial_momentum << std::endl;
    m_status = kFatal;
    return false;
  }

  static const auto gPosBcOut = gGeom.GetGlobalPosition("BC3-X1").Z();
  static const auto localBcOut = gGeom.GetLocalZ("BC3-X1");

  const Double_t xIn = xInit;
  const Double_t yIn = yInit;
  const Double_t uIn = uInit;
  const Double_t vIn = vInit;
  const Double_t pz = m_initial_momentum/std::sqrt(1.+uIn*uIn+vIn*vIn);
  const ThreeVector posIn(xIn-50.,yIn,gPosBcOut-localBcOut);
  const ThreeVector momIn(pz*uIn , pz*vIn, pz);
  //  std::cout << xIn << "\t" <<yIn << "\t" << uIn << "\t" << vIn << std::endl;
  RKCordParameter iniCord(posIn,momIn);
  RKCordParameter prevCord;
  RKHitPointContainer preHPntCont;

  m_HitPointCont = RK::MakeHSHPContainer();
  RKHitPointContainer prevHPntCont;
  m_status = (RKstatus)RK::Extrap(iniCord,m_HitPointCont);
  //  std::cout << m_status << std::endl;
  if (m_status != kPassed || !SaveTrack()) return false;

  return true;
}

    
//_____________________________________________________________________________
Bool_t
HSTrack::SaveTrack()
{
  //  const RKcalcHitPoint& hpTgt  = m_HitPointCont.rbegin()->second;

  const Int_t TGTid = gGeom.DetectorId("K18Target");
  const RKcalcHitPoint& hpTgt  = m_HitPointCont.HitPointOfLayer(TGTid);
  const ThreeVector& posTgt = hpTgt.PositionInGlobal();
  const ThreeVector& momTgt = hpTgt.MomentumInGlobal();
  m_tgt_position = gGeom.Global2LocalPos(TGTid, posTgt);
  m_tgt_momentum = gGeom.Global2LocalDir(TGTid, momTgt);

  const Int_t IdBH2 =gGeom.DetectorId("BH2");
  const RKcalcHitPoint& hpBH2 = m_HitPointCont.HitPointOfLayer(IdBH2);
  const ThreeVector &posBH2 = hpBH2.PositionInGlobal();
  const ThreeVector &momBH2 = hpBH2.MomentumInGlobal();  
  m_bh2_position = gGeom.Global2LocalPos(IdBH2,posBH2);
  m_bh2_momentum = gGeom.Global2LocalDir(IdBH2,momBH2);
  
  const Int_t IdVP1 = gGeom.DetectorId("VPHS1");
  const RKcalcHitPoint& hpVP1 = m_HitPointCont.HitPointOfLayer(IdVP1);
  const ThreeVector &posVP1 = hpVP1.PositionInGlobal();
  const ThreeVector &momVP1 = hpVP1.MomentumInGlobal();
  m_vp1_position = gGeom.Global2LocalPos(IdVP1,posVP1);
  m_vp1_momentum = gGeom.Global2LocalDir(IdVP1,momVP1);

  const Int_t IdVP3 = gGeom.DetectorId("VPHS3");
  const RKcalcHitPoint& hpVP3 = m_HitPointCont.HitPointOfLayer(IdVP3);
  const ThreeVector &posVP3 = hpVP3.PositionInGlobal();
  const ThreeVector &momVP3 = hpVP3.MomentumInGlobal();
  m_vp3_position = gGeom.Global2LocalPos(IdVP3,posVP3);
  m_vp3_momentum = gGeom.Global2LocalDir(IdVP3,momVP3);
  
  const Int_t IdHtof = m_HitPointCont.rbegin()->first;
  const RKcalcHitPoint &hpHtof = m_HitPointCont.rbegin()->second;
  const ThreeVector &posHtof = hpHtof.PositionInGlobal();
  const ThreeVector &momHtof = hpHtof.MomentumInGlobal();
  m_htof_position = gGeom.Global2LocalPos(IdHtof,posHtof);
  m_htof_momentum = gGeom.Global2LocalDir(IdHtof,momHtof);
  

  m_path_length_total = std::abs(hpBH2.PathLength()-hpTgt.PathLength());

  return true;
}
