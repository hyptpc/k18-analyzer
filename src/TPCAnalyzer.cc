// -*- C++ -*-

#include "TPCAnalyzer.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include "DCGeomMan.hh"
#include "DCGeomRecord.hh"
#include "DebugCounter.hh"
#include "DebugTimer.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "TPCRawData.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"
#include "TPCPadHelper.hh"
#include "TPCRawHit.hh"
#include "TPCHit.hh"
#include "TRandom3.h"


namespace
{
TRandom3 RandGen;

const auto& gUser   = UserParamMan::GetInstance();
const auto& gGeom   = DCGeomMan::GetInstance();

const double ztgt = tpc::ZTarget;
}
//_____________________________________________________________________________
TPCAnalyzer::TPCAnalyzer()
  : m_is_decoded(n_type),
    m_TPCHitCont(NumOfLayersTPC+1)
{
  for(Int_t i=0; i<n_type; ++i){
    m_is_decoded[i] = false;
  }
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCAnalyzer::~TPCAnalyzer()
{
  ClearTPCHits();

  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::DecodeTPCHits(TPCRawData &TPCrawData,
			   const std::vector<Double_t> clock)
{
  if(m_is_decoded[kTPC]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  ClearTPCHits();
  //ClearTPCClusters();

  for(Int_t layer=0; layer<=NumOfLayersTPC; ++layer){
    for(const auto& rhit: TPCrawData.GetTPCCorHitContainer(layer)){
      auto hit = new TPCHit(rhit);
      Int_t row  = hit -> GetRow();
      Int_t cobo_id = tpc::GetCoBoId(layer, row);
      if(cobo_id >= clock.size()){
	delete hit;
	return false; //No cobo input
      }

      if(hit->DoFit() && hit->Calculate(clock[cobo_id])){
        m_TPCHitCont[layer].push_back(hit);
      }else{
        delete hit;
      }
    }
  }

// #if 0 // Cluster analysis will be done by RecalcTPCHits() in Dst.
//   static const Double_t MaxYDif = gUser.GetParameter("MaxYDifClusterTPC");
//   for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
//     MakeUpTPCClusters(m_TPCHitCont[layer], m_TPCClCont[layer], MaxYDif);
//   }
// #endif

  m_is_decoded[kTPC] = true;
  return true;
}


//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCHits()
{
  del::ClearContainerAll(m_TPCHitCont);
}
