#include "HypTPCVertex.hh"
#include "HypTPCTask.hh"

#include "TPCRiemannTrack.hh"
#include "TPCTrack.hh"

ClassImp(HypTPCVertex);

namespace {
  HypTPCTask &fHypTPCTask = HypTPCTask::GetInstance();
}

HypTPCVertex::HypTPCVertex() {
  
  isInitialized = false;
  
  fVertexFactory = new genfit::GFRaveVertexFactory(3, false);
  //  fGFRaveVertexMethod = "avf-smoothing:1-Tini:256-ratio:0.25-sigmacut:5";
  //  fGFRaveVertexMethod = "kalman-smoothing:1";
  fGFRaveVertexMethod = "avr-smoothing:1-minweight:0.5-primcut:2-seccut:4";
  fVertexFactory -> setMethod(fGFRaveVertexMethod.Data());
}

HypTPCVertex::~HypTPCVertex() {

  genfitTrackArray.clear();
  gfTrackToVertexArray.clear();

  delete fVertexFactory;
  
}

Bool_t HypTPCVertex::Execute(std::vector<TPCTrack *> &tpcTrackArray,std::vector<TPCVertex *> &gfVertexArray) {
  
  for (Int_t nt = 0 ; nt < tpcTrackArray.size() ; nt++) {
    auto tpcTrack = tpcTrackArray.at(nt);
    genfitTrackArray.push_back(tpcTrack->GetGenFitTrack());
    
    if (tpcTrack->GetStatus() == 0) continue;
    gfTrackToVertexArray.push_back(tpcTrack->GetGenFitTrack());
  }

  auto numTracks = gfTrackToVertexArray.size();

  if (numTracks > 1)
    isInitialized = true;    
  
  if (isInitialized == false) {
    std::cout << "Number of Tracks is :  " << numTracks  << std::endl;
    return false;
  }

  //  fVertexFactory -> setBeamspot(fVertex, fCovVertex);
  
  std::vector<genfit::GFRaveVertex *> vertices;
	std::cout << "here"<< std::endl;
  
  try {
    fVertexFactory->findVertices(&vertices, gfTrackToVertexArray);
  } catch (...) {
  }

  Int_t numVertices = vertices.size();

  for (Int_t iVtx = 0 ; iVtx < numVertices ; iVtx++) {

    genfit::GFRaveVertex *vertex = static_cast<genfit::GFRaveVertex*>(vertices[iVtx]);
    TPCVertex *tpcVertex = new TPCVertex(*vertex);
					 
    if (numVertices > 0) {
      auto pos = 10*vertex->getPos();
      std::cout << pos.X() << "\t" << pos.Y() << "\t" << pos.Z() << std::endl;
    }
    
    gfVertexArray.push_back(tpcVertex);
    
    Int_t numTracks = vertex->getNTracks();

    std::cout << numTracks << std::endl;

    for (Int_t iTrack = 0 ; iTrack < numTracks ; iTrack++) {
      genfit::GFRaveTrackParameters *par = vertex -> getParameters(iTrack);
      const genfit::Track *track = par->getTrack();

      std::vector<genfit::Track *>::iterator iter = std::find(gfTrackToVertexArray.begin(), gfTrackToVertexArray.end(), track);

      if (iter != gfTrackToVertexArray.end()) {
	
	Int_t index = std::distance(gfTrackToVertexArray.begin(), iter);
	genfit::Track *gfTrack = gfTrackToVertexArray.at(index);
	std::vector<genfit::Track *>::iterator iter_org = std::find(genfitTrackArray.begin(), genfitTrackArray.end(), gfTrack);
	Int_t index_org = std::distance(genfitTrackArray.begin(), iter_org);

	std::cout << index << "\t" << index_org << std::endl;
	
	auto tpcTrack = tpcTrackArray.at(index_org);

	if (tpcTrack->GetVertexID() < 0)
	  tpcTrack->SetVertexID(iVtx);
	
	TVector3 momVertex, posVertex;
	fHypTPCTask.GetMomentumWithVertex(gfTrack, 10*vertex->getPos(), momVertex, posVertex);

	tpcTrack->SetMomentum(momVertex);
	tpcTrack->SetVertex(posVertex);
	
	std::cout << "MomVtx" << std::endl;
	std::cout << momVertex.X() << "\t" << momVertex.Y() << "\t" << momVertex.Z() << std::endl;
	std::cout << posVertex.X() << "\t" << posVertex.Y() << "\t" << posVertex.Z() << std::endl;
	  
	//	Double_t effCurv1, effCurv2, effCurv3;
	//	Double_t charge;
	//	charge = fHypTPCTask.DetermineCharge(gfTrack, vertex->getPos(), effCurv1, effCurv2, effCurv3);
      }
      
    }
  }
  
  return true;
}

