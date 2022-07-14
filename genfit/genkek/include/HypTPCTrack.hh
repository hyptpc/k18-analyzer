//  Authors: Wooseung Jung

#ifndef HYPTPCTRACK_HH
#define HYPTPCTRACK_HH

//GenKEK
#include "HypTPCSpacepointMeasurement.hh"
#include "HypTPCHit.hh"

//k18-analyzer
#include "TPCLocalTrackHelix.hh"
#include "TPCLocalTrack.hh"
#include "TPCLTrackHit.hh"

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
  virtual ~HypTPCTrack(){}
  void Init();
  genfit::Track* GetTrack(int ith) const;
  void AddReps(int ith, int pdg);
  void AddHelixTrack(int pdg, TPCLocalTrackHelix *tp);
  void AddHelixTrack(std::vector<int> pdg, TPCLocalTrackHelix *tp);
  int GetNTrack() const { return _genfitTrackArray -> GetEntriesFast(); }
  int GetNHit(int ith) const;

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
