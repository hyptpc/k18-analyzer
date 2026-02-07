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
void BuildTPCBasic();
void BuildTPCTracking();
void BuildTPCBcOutTracking();
void BuildTPCHitBcOutTracking();
void BuildTPCHelixTracking();
}

#endif
