// HistTools.h

#ifndef HistTools_h
#define HistTools_h

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMacro.h>
#include <TString.h>
#include <TSystem.h>
#include <TVector3.h>

namespace hist
{
void BuildStatus();
void BuildTriggerFlag();
void BuildDAQ();

void BuildHodoRaw(Bool_t flag_beam_particle=false);
void BuildHodoHit(Bool_t flag_beam_particle=false);
void BuildHodoCluster(Bool_t flag_beam_particle=false);

void BuildDCRaw(const TString& dcname, Bool_t flag_beam_particle=false);
void BuildDCHit(const TString& dcname, Bool_t flag_beam_particle=false);
void BuildDCTrack(const TString& dcname, Bool_t flag_beam_particle=false);

void BuildTPCHit();
void BuildTPCTrackingCommon();
void BuildTPCLineTrackParam();
void BuildTPCTracking(Bool_t calib_flag = false);

void BuildTPCBcOutTracking(Bool_t calib_flag = false);
void BuildTPCHitBcOutTracking(Bool_t calib_flag = false);
void BuildTPCHelixTracking(Bool_t calib_flag = false);
void BuildTPCHSHelixResidual();
void BuildTPCHelixLambda();
void BuildTPCLambdaEta();
void BuildTPCHelixK0Short();

void BuildD5Tracking(Bool_t flag_beam_particle=false);
void BuildD5WireResidual();

enum CoBoClockTimeFlags : UInt_t {
  kCoBoClockTime_Hit     = 1u << 0,
  kCoBoClockTime_Cluster = 1u << 1,
  kCoBoClockTime_Track   = 1u << 2,
};
void BuildCoBoClockTime(UInt_t flags = 0);
}

#endif
