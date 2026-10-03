// -*- C++ -*-
/**
 *  file: BeamRunCondition.hh
 *  Resolve beam species / momentum / charge from run_num.
 *
 *  Priority (particle, momentum, and charge independently, strongest first):
 *    1. UserParam BeamParticle / BeamMom / BeamCharge (explicit override)
 *    2. run table
 *  A run without a table entry has no beam condition unless UserParam gives it
 *  (no guessing from BeamSpecies). Normal data jobs should not need BeamCharge:
 *  every current table row carries charge. BeamCharge is for Geant4 / special
 *  runs / intentional polarity overrides.
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
  std::optional<Int_t> charge;            // ±1 when set
};

struct ResolvedBeamCondition {
  UInt_t run_num = 0;
  BeamSpecies particle = BeamSpecies::Unknown;
  std::optional<Double_t> momentum_gev_c;
  std::optional<Int_t> charge;
  BeamConditionSource particle_source = BeamConditionSource::None;
  BeamConditionSource momentum_source = BeamConditionSource::None;
  BeamConditionSource charge_source = BeamConditionSource::None;

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
  Bool_t ValidCharge() const
  {
    return charge.has_value()
        && (*charge == -1 || *charge == 1);
  }
};

const char* ToString(BeamSpecies species);
const char* ToString(BeamConditionSource source);

/// Table lookup only (no UserParam). Missing run_num → nullopt.
std::optional<BeamRunCondition> FindBeamRunCondition(UInt_t run_num);

/// Table, then UserParam overrides. Logs once per run_num via spdlog.
/// The result for the last run_num is kept, so calling it for every track / event is cheap.
ResolvedBeamCondition ResolveBeamCondition(UInt_t run_num);

/// Beam charge for run_num (table / BeamCharge). Unresolved: warn once per run, return fallback.
Int_t BeamChargeWithFallback(UInt_t run_num, Int_t fallback, const char* caller);

/// Pion/Kaon/Proton → mass [GeV/c^2]; Unknown → nullopt.
std::optional<Double_t> BeamMassGeV(BeamSpecies species);

/// Expand ranges, check overlaps / invalid entries. Returns false on failure.
Bool_t ValidateBeamRunConditionTable();

#endif
