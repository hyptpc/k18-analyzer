// -*- C++ -*-
/**
 *  file: BeamRunCondition.hh
 *  Resolve beam species / momentum from run_num.
 *
 *  Priority (particle and momentum independently, strongest first):
 *    1. UserParam BeamParticle / BeamMom (explicit override)
 *    2. run table
 *  A run without a table entry has no beam condition unless UserParam gives it (no guessing).
 *
 *  Mass is derived from species (no separate mass field).
 */

#ifndef BEAM_RUN_CONDITION_HH
#define BEAM_RUN_CONDITION_HH

#include <cmath>
#include <optional>

#include <Rtypes.h>

enum class BeamSpecies : Int_t {
  Pion    = 0,
  Kaon    = 1,
  Proton  = 2,
  Unknown = 3,
};

enum class BeamConditionSource : Int_t {
  None = 0,
  RunTable,
  UserParam,
};

struct BeamRunCondition {
  UInt_t run_num = 0;
  BeamSpecies particle = BeamSpecies::Unknown;
  std::optional<Double_t> momentum_gev_c; // positive magnitude [GeV/c]
};

struct ResolvedBeamCondition {
  UInt_t run_num = 0;
  BeamSpecies particle = BeamSpecies::Unknown;
  std::optional<Double_t> momentum_gev_c;
  BeamConditionSource particle_source = BeamConditionSource::None;
  BeamConditionSource momentum_source = BeamConditionSource::None;

  Bool_t ValidParticle() const
  {
    return particle == BeamSpecies::Pion
        || particle == BeamSpecies::Kaon
        || particle == BeamSpecies::Proton;
  }
  Bool_t ValidMomentum() const
  {
    return momentum_gev_c.has_value()
        && *momentum_gev_c > 0.
        && std::isfinite(*momentum_gev_c);
  }
};

const char* ToString(BeamSpecies species);
const char* ToString(BeamConditionSource source);

/// Table lookup only (no UserParam). Missing run_num → nullopt.
std::optional<BeamRunCondition> FindBeamRunCondition(UInt_t run_num);

/// Table, then UserParam overrides. Logs once per run_num via spdlog.
/// The result for the last run_num is kept, so calling it for every track / event is cheap.
ResolvedBeamCondition ResolveBeamCondition(UInt_t run_num);

/// Pion/Kaon/Proton → mass [GeV/c^2]; Unknown → nullopt.
std::optional<Double_t> BeamMassGeV(BeamSpecies species);

/// Expand ranges, check overlaps / invalid entries. Returns false on failure.
Bool_t ValidateBeamRunConditionTable();

#endif
