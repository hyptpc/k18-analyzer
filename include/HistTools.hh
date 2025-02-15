// HistTools.h

#ifndef HistTools_h
#define HistTools_h 1

#include <string>
#include <vector>
#include <iostream>
#include <cmath>

#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMacro.h"
#include "TSystem.h"
#include "TString.h"

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
}

#endif
