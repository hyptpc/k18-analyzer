//  Authors: Wooseung Jung

#ifndef HYPTPCTRACK_HH
#define HYPTPCTRACK_HH

//GenKEK
#include "HypTPCSpacepointMeasurement.hh"
#include "HypTPCHit.hh"

//k18-analyzer
#include "TPCLocalTrackHelix.hh"
#include "TPCLocalTrack.hh"

//GenFit
#include <AbsMeasurement.h>
#include <Track.h>
#include <MeasurementFactory.h>
#include <MeasurementProducer.h>

//ROOT
#include <TClonesArray.h>

//STL
#include <vector>

class HypTPCTrack{

public:

  HypTPCTrack();
  virtual ~HypTPCTrack(){};
  void Clear();
  genfit::Track* GetTrack(int ith) const;
  void AddReps(int ith, int pdg);
  //single pid hypothesis
  void AddHelixTrack(int pdg, TPCLocalTrackHelix *tp);
  void AddLinearTrack(int pdg, TPCLocalTrack *tp, double momentum);
  void AddReconstructedTrack(int pdg, TVector3 posSeed, TVector3 momSeed); //add recontructed track for extrapolation
  //multiple pid hypotheses
  void AddHelixTrack(std::vector<int> pdg, TPCLocalTrackHelix *tp);
  void AddLinearTrack(std::vector<int> pdg, TPCLocalTrack *tp, double momentum);
  int GetNTrack() const { return _genfitTrackArray -> GetEntriesFast(); }

protected:

  static TClonesArray *_hitClusterArray;
  static TClonesArray *_genfitTrackArray;

private:

  int TPCDetID=0;
  genfit::MeasurementFactory<genfit::AbsMeasurement> *_measurementFactory;
  genfit::MeasurementProducer<HypTPCHit, genfit::HypTPCSpacepointMeasurement> *_measurementProducer;

  ClassDef(HypTPCTrack, 1)

}; //class HypTPCTrack.hh

#endif // HypTPCTrack_hh
