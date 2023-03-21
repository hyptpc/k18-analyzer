#ifndef HYPTPC_VERTEX_HH
#define HYPTPC_VERTEX_HH

//GenKEK
#include "HypTPCSpacepointMeasurement.hh"
#include "HypTPCHit.hh"
#include "HypTPCTask.hh"
#include "HypTPCTrack.hh"

//k18-analyzer
#include "TPCLocalTrackHelix.hh"
#include "TPCLocalTrack.hh"
#include "TPCLTrackHit.hh"
#include "TPCRiemannTrack.hh"
#include "TPCHitCluster.hh"
#include "TPCTrack.hh"
#include "TPCVertex.hh"

//GenFit
#include <AbsMeasurement.h>
#include <Track.h>
#include <MeasurementFactory.h>
#include <MeasurementProducer.h>

#include <GFRaveVertexFactory.h>
#include <GFRaveVertex.h>
#include <GFRaveTrackParameters.h>

//ROOT
#include <TClonesArray.h>

//STL
#include <vector>

class HypTPCVertex {

public:

  HypTPCVertex();
  virtual ~HypTPCVertex();
  
  Bool_t Execute(std::vector<TPCTrack *> &tpcTrackArray, std::vector<TPCVertex*> &gfVertexArray);

private:

  Bool_t isInitialized;

  TString fGFRaveVertexMethod;

  std::vector<genfit::Track*> genfitTrackArray;
  std::vector<genfit::Track*> gfTrackToVertexArray;

protected:

  genfit::GFRaveVertexFactory *fVertexFactory = nullptr;

  ClassDef(HypTPCVertex, 1)
  
};

#endif

  
