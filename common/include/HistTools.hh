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

void BuildHodoRaw(Bool_t flag_beam_particle=false);
void BuildHodoHit(Bool_t flag_beam_particle=false);

void BuildDCRaw(Bool_t flag_beam_particle=false);
void BuildDCHit(Bool_t flag_beam_particle=false);

void BuildDAQ();
}

#endif
