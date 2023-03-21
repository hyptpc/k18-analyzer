#include "TPCVertex.hh"

#include <GFRaveVertex.h>

TPCVertex::TPCVertex() {
  Clear();
}

TPCVertex::TPCVertex(const genfit::GFRaveVertex &vertex) {

  Clear();

  fVertexID = vertex.getId();
  fNumTrack = vertex.getNTracks();

  fPosition = 10. * vertex.getPos();

  fCovariance.ResizeTo(vertex.getCov());
  fCovariance = 100. * vertex.getCov();

  fChisqr = vertex.getChi2()/vertex.getNdf();
}

TPCVertex::~TPCVertex() {}

void TPCVertex::Clear() {

  fVertexID = -1;
  fNumTrack = -1;

  fPosition.SetXYZ(0,0,0);
  fCovariance(3,3);
  fChisqr = -9999;

  isCollisionVertex = false;
  isDecayVertex = false;
}

