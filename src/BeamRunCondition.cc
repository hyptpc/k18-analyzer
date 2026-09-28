// -*- C++ -*-

#include "BeamRunCondition.hh"

#include <cmath>
#include <unordered_map>
#include <unordered_set>

#include "DatabasePDG.hh"
#include "UserParamMan.hh"

#include <spdlog/spdlog.h>

namespace
{

struct BeamRunRange {
  UInt_t first_run_num;
  UInt_t last_run_num;
  BeamSpecies particle;
  Double_t momentum_gev_c;  // positive magnitude [GeV/c]; always set
};

// Beam condition per run, from the run summaries. Each entry covers consecutive runs with the
// same particle and momentum; runs absent from the table (dummy, junk, clock, beam tune,
// no beam, runs without a known particle, ...) have no entry, i.e. no beam condition unless
// UserParam gives one.
// Ordered by run_num; the section comments give the beam momentum of the entries below them.
const BeamRunRange kBeamRunRanges[] = {
  // --- 0.735 GeV/c ---
  {1916, 1933, BeamSpecies::Pion,    0.735},
  {1934, 1947, BeamSpecies::Kaon,    0.735},
  {1971, 2051, BeamSpecies::Kaon,    0.735},
  {2053, 2057, BeamSpecies::Kaon,    0.735},

  // --- 1.144 GeV/c ---
  {2058, 2062, BeamSpecies::Pion,    1.144},

  // --- 0.645 GeV/c ---
  {2105, 2121, BeamSpecies::Kaon,    0.645},
  {2123, 2123, BeamSpecies::Kaon,    0.645},
  {2125, 2125, BeamSpecies::Kaon,    0.645},
  {2127, 2132, BeamSpecies::Kaon,    0.645},
  {2134, 2140, BeamSpecies::Kaon,    0.645},
  {2142, 2143, BeamSpecies::Kaon,    0.645},

  // --- 0.735 GeV/c ---
  {2144, 2147, BeamSpecies::Kaon,    0.735},
  {2150, 2150, BeamSpecies::Kaon,    0.735},
  {2153, 2168, BeamSpecies::Kaon,    0.735},
  {2170, 2171, BeamSpecies::Kaon,    0.735},
  {2173, 2174, BeamSpecies::Kaon,    0.735},
  {2177, 2180, BeamSpecies::Kaon,    0.735},
  {2208, 2213, BeamSpecies::Kaon,    0.735},
  {2215, 2216, BeamSpecies::Kaon,    0.735},
  {2260, 2260, BeamSpecies::Kaon,    0.735},
  {2263, 2264, BeamSpecies::Kaon,    0.735},
  {2266, 2272, BeamSpecies::Kaon,    0.735},
  {2274, 2287, BeamSpecies::Kaon,    0.735},
  {2290, 2300, BeamSpecies::Kaon,    0.735},
  {2304, 2306, BeamSpecies::Kaon,    0.735},
  {2309, 2309, BeamSpecies::Kaon,    0.735},
  {2311, 2318, BeamSpecies::Kaon,    0.735},
  {2321, 2324, BeamSpecies::Kaon,    0.735},
  {2327, 2327, BeamSpecies::Kaon,    0.735},
  {2329, 2329, BeamSpecies::Kaon,    0.735},
  {2331, 2331, BeamSpecies::Kaon,    0.735},
  {2333, 2333, BeamSpecies::Kaon,    0.735},
  {2336, 2337, BeamSpecies::Kaon,    0.735},
  {2339, 2339, BeamSpecies::Kaon,    0.735},
  {2341, 2341, BeamSpecies::Kaon,    0.735},
  {2343, 2344, BeamSpecies::Kaon,    0.735},
  {2346, 2347, BeamSpecies::Kaon,    0.735},
  {2349, 2349, BeamSpecies::Kaon,    0.735},
  {2351, 2352, BeamSpecies::Kaon,    0.735},
  {2354, 2357, BeamSpecies::Kaon,    0.735},
  {2359, 2366, BeamSpecies::Kaon,    0.735},
  {2368, 2370, BeamSpecies::Kaon,    0.735},
  {2372, 2374, BeamSpecies::Kaon,    0.735},
  {2376, 2380, BeamSpecies::Kaon,    0.735},
  {2382, 2393, BeamSpecies::Kaon,    0.735},
  {2395, 2395, BeamSpecies::Kaon,    0.735},
  {2408, 2421, BeamSpecies::Kaon,    0.735},
  {2423, 2428, BeamSpecies::Kaon,    0.735},
  {2431, 2432, BeamSpecies::Kaon,    0.735},
  {2434, 2443, BeamSpecies::Kaon,    0.735},
  {2447, 2454, BeamSpecies::Kaon,    0.735},
  {2456, 2460, BeamSpecies::Kaon,    0.735},
  {2462, 2463, BeamSpecies::Kaon,    0.735},
  {2465, 2466, BeamSpecies::Kaon,    0.735},
  {2468, 2471, BeamSpecies::Kaon,    0.735},
  {2485, 2487, BeamSpecies::Kaon,    0.735},

  // --- 1.000 GeV/c ---
  {2489, 2489, BeamSpecies::Pion,    1.000},
  {2491, 2492, BeamSpecies::Proton,  1.000},

  // --- 0.814 GeV/c ---
  {2494, 2494, BeamSpecies::Proton,  0.814},
  {2496, 2498, BeamSpecies::Proton,  0.814},
  {2500, 2500, BeamSpecies::Kaon,    0.814},
  {2502, 2502, BeamSpecies::Pion,    0.814},

  // --- 0.645 GeV/c ---
  {2504, 2504, BeamSpecies::Proton,  0.645},
  {2506, 2506, BeamSpecies::Kaon,    0.645},
  {2508, 2509, BeamSpecies::Pion,    0.645},

  // --- 0.400 GeV/c ---
  {2511, 2512, BeamSpecies::Pion,    0.400},

  // --- 0.300 GeV/c ---
  {2514, 2514, BeamSpecies::Pion,    0.300},

  // --- 1.000 GeV/c ---
  {2516, 2516, BeamSpecies::Pion,    1.000},
  {2518, 2518, BeamSpecies::Proton,  1.000},

  // --- 0.814 GeV/c ---
  {2520, 2520, BeamSpecies::Pion,    0.814},
  {2522, 2522, BeamSpecies::Proton,  0.814},

  // --- 0.645 GeV/c ---
  {2524, 2524, BeamSpecies::Pion,    0.645},
  {2527, 2527, BeamSpecies::Proton,  0.645},

  // --- 0.400 GeV/c ---
  {2529, 2529, BeamSpecies::Pion,    0.400},
  {2531, 2531, BeamSpecies::Proton,  0.400},
  {2534, 2535, BeamSpecies::Pion,    0.400},

  // --- 0.300 GeV/c ---
  {2537, 2537, BeamSpecies::Pion,    0.300},
  {2542, 2542, BeamSpecies::Pion,    0.300},
  {2544, 2544, BeamSpecies::Proton,  0.300},
  {2546, 2546, BeamSpecies::Proton,  0.300},

  // --- 0.735 GeV/c ---
  {2560, 2562, BeamSpecies::Pion,    0.735},
  {2568, 2570, BeamSpecies::Pion,    0.735},
  {2580, 2581, BeamSpecies::Proton,  0.735},

  // --- 0.933 GeV/c ---
  {2585, 2585, BeamSpecies::Proton,  0.933},

  // --- 0.755 GeV/c ---
  {2587, 2587, BeamSpecies::Proton,  0.755},

  // --- 0.715 GeV/c ---
  {2589, 2590, BeamSpecies::Proton,  0.715},

  // --- 0.645 GeV/c ---
  {2592, 2592, BeamSpecies::Proton,  0.645},

  // --- 1.000 GeV/c ---
  {2594, 2596, BeamSpecies::Pion,    1.000},
  {2599, 2599, BeamSpecies::Pion,    1.000},
  {2601, 2604, BeamSpecies::Pion,    1.000},
  {2606, 2607, BeamSpecies::Pion,    1.000},

  // --- 0.715 GeV/c ---
  {2609, 2612, BeamSpecies::Kaon,    0.715},
  {2615, 2616, BeamSpecies::Kaon,    0.715},
  {2618, 2619, BeamSpecies::Kaon,    0.715},
  {2622, 2636, BeamSpecies::Kaon,    0.715},
  {2647, 2648, BeamSpecies::Kaon,    0.715},
  {2650, 2653, BeamSpecies::Kaon,    0.715},
  {2655, 2655, BeamSpecies::Kaon,    0.715},
  {2680, 2680, BeamSpecies::Kaon,    0.715},
  {2682, 2684, BeamSpecies::Kaon,    0.715},
  {2686, 2687, BeamSpecies::Kaon,    0.715},
  {2689, 2693, BeamSpecies::Kaon,    0.715},
  {2695, 2695, BeamSpecies::Kaon,    0.715},
  {2697, 2697, BeamSpecies::Kaon,    0.715},
  {2699, 2700, BeamSpecies::Kaon,    0.715},
  {2702, 2706, BeamSpecies::Kaon,    0.715},

  // --- 0.755 GeV/c ---
  {2708, 2711, BeamSpecies::Kaon,    0.755},
  {2713, 2713, BeamSpecies::Kaon,    0.755},
  {2715, 2715, BeamSpecies::Kaon,    0.755},
  {2718, 2719, BeamSpecies::Kaon,    0.755},
  {2721, 2726, BeamSpecies::Kaon,    0.755},
  {2729, 2729, BeamSpecies::Kaon,    0.755},
  {2731, 2732, BeamSpecies::Kaon,    0.755},
  {2734, 2736, BeamSpecies::Kaon,    0.755},
  {2738, 2739, BeamSpecies::Kaon,    0.755},
  {2741, 2742, BeamSpecies::Kaon,    0.755},
  {2747, 2755, BeamSpecies::Kaon,    0.755},
  {2757, 2757, BeamSpecies::Kaon,    0.755},
  {2759, 2764, BeamSpecies::Kaon,    0.755},
  {2767, 2772, BeamSpecies::Kaon,    0.755},

  // --- 0.715 GeV/c ---
  {2774, 2785, BeamSpecies::Kaon,    0.715},
  {2787, 2788, BeamSpecies::Kaon,    0.715},
  {2790, 2794, BeamSpecies::Kaon,    0.715},
  {2796, 2797, BeamSpecies::Kaon,    0.715},

  // --- 0.790 GeV/c ---
  {2799, 2800, BeamSpecies::Kaon,    0.790},
  {2803, 2806, BeamSpecies::Kaon,    0.790},
  {2815, 2818, BeamSpecies::Kaon,    0.790},
  {2820, 2820, BeamSpecies::Kaon,    0.790},
  {2822, 2828, BeamSpecies::Kaon,    0.790},
  {2833, 2834, BeamSpecies::Kaon,    0.790},
  {2836, 2837, BeamSpecies::Kaon,    0.790},
  {2841, 2842, BeamSpecies::Kaon,    0.790},
  {2844, 2844, BeamSpecies::Kaon,    0.790},

  // --- 0.600 GeV/c ---
  {2846, 2848, BeamSpecies::Kaon,    0.600},

  // --- 0.933 GeV/c ---
  {2851, 2851, BeamSpecies::Pion,    0.933},
  {2853, 2854, BeamSpecies::Pion,    0.933},
  {2856, 2861, BeamSpecies::Kaon,    0.933},
  {2863, 2868, BeamSpecies::Kaon,    0.933},
  {2870, 2875, BeamSpecies::Kaon,    0.933},

  // --- 0.685 GeV/c ---
  {2883, 2884, BeamSpecies::Kaon,    0.685},
  {2887, 2889, BeamSpecies::Kaon,    0.685},
  {2891, 2896, BeamSpecies::Kaon,    0.685},
  {2898, 2898, BeamSpecies::Kaon,    0.685},
  {2901, 2901, BeamSpecies::Kaon,    0.685},
  {2903, 2906, BeamSpecies::Kaon,    0.685},
  {2909, 2910, BeamSpecies::Kaon,    0.685},

  // --- 0.814 GeV/c ---
  {2912, 2914, BeamSpecies::Kaon,    0.814},
  {2916, 2918, BeamSpecies::Kaon,    0.814},
  {2921, 2924, BeamSpecies::Kaon,    0.814},
  {2926, 2928, BeamSpecies::Kaon,    0.814},
  {2930, 2935, BeamSpecies::Kaon,    0.814},

  // --- 0.842 GeV/c ---
  {2939, 2940, BeamSpecies::Kaon,    0.842},
  {2942, 2943, BeamSpecies::Kaon,    0.842},

  // --- 0.645 GeV/c ---
  {2977, 2979, BeamSpecies::Kaon,    0.645},
  {2981, 2981, BeamSpecies::Kaon,    0.645},
  {2983, 2983, BeamSpecies::Kaon,    0.645},
  {2985, 2990, BeamSpecies::Kaon,    0.645},

  // --- 0.870 GeV/c ---
  {2992, 2993, BeamSpecies::Kaon,    0.870},
  {2995, 2997, BeamSpecies::Kaon,    0.870},
  {2999, 3000, BeamSpecies::Kaon,    0.870},

  // --- 0.665 GeV/c ---
  {3005, 3010, BeamSpecies::Kaon,    0.665},
  {3012, 3012, BeamSpecies::Kaon,    0.665},
  {3014, 3014, BeamSpecies::Kaon,    0.665},
  {3016, 3016, BeamSpecies::Kaon,    0.665},
  {3018, 3019, BeamSpecies::Kaon,    0.665},
  {3021, 3021, BeamSpecies::Kaon,    0.665},
  {3023, 3023, BeamSpecies::Kaon,    0.665},
  {3025, 3025, BeamSpecies::Kaon,    0.665},
  {3027, 3027, BeamSpecies::Kaon,    0.665},

  // --- 0.842 GeV/c ---
  {3030, 3032, BeamSpecies::Kaon,    0.842},
  {3034, 3036, BeamSpecies::Kaon,    0.842},
  {3038, 3038, BeamSpecies::Kaon,    0.842},
  {3040, 3040, BeamSpecies::Kaon,    0.842},
  {3042, 3043, BeamSpecies::Kaon,    0.842},

  // --- 0.735 GeV/c ---
  {3642, 3644, BeamSpecies::Kaon,    0.735},
  {3646, 3648, BeamSpecies::Kaon,    0.735},

  // --- 0.964 GeV/c ---
  {3650, 3652, BeamSpecies::Pion,    0.964},
  {3654, 3666, BeamSpecies::Pion,    0.964},
  {3668, 3669, BeamSpecies::Pion,    0.964},

  // --- 1.144 GeV/c ---
  {3688, 3703, BeamSpecies::Proton,  1.144},

  // --- 0.964 GeV/c ---
  {3704, 3712, BeamSpecies::Proton,  0.964},
  {3730, 3737, BeamSpecies::Pion,    0.964},
  {3739, 3742, BeamSpecies::Pion,    0.964},

  // --- 1.144 GeV/c ---
  {3743, 3749, BeamSpecies::Pion,    1.144},

  // --- 0.735 GeV/c ---
  {3750, 3759, BeamSpecies::Pion,    0.735},
  {3761, 3761, BeamSpecies::Pion,    0.735},
  {3762, 3779, BeamSpecies::Kaon,    0.735},
  {3781, 3781, BeamSpecies::Kaon,    0.735},

  // --- 0.902 GeV/c ---
  {3782, 3782, BeamSpecies::Kaon,    0.902},
  {3784, 3784, BeamSpecies::Kaon,    0.902},

  // --- 0.685 GeV/c ---
  {3813, 3813, BeamSpecies::Kaon,    0.685},

  // --- 0.902 GeV/c ---
  {3830, 3844, BeamSpecies::Kaon,    0.902},
  {3846, 3847, BeamSpecies::Kaon,    0.902},

  // --- 0.933 GeV/c ---
  {3848, 3853, BeamSpecies::Kaon,    0.933},
  {3855, 3855, BeamSpecies::Kaon,    0.933},
  {3857, 3861, BeamSpecies::Kaon,    0.933},
  {3863, 3867, BeamSpecies::Kaon,    0.933},

  // --- 0.814 GeV/c ---
  {3868, 3876, BeamSpecies::Kaon,    0.814},
};

Bool_t
IsKnownParticle(BeamSpecies s)
{
  return s == BeamSpecies::Pion
      || s == BeamSpecies::Kaon
      || s == BeamSpecies::Proton;
}

std::unordered_map<UInt_t, BeamRunCondition>
BuildTable(Bool_t* ok)
{
  std::unordered_map<UInt_t, BeamRunCondition> map;
  *ok = true;

  for (const auto& r : kBeamRunRanges) {
    if (r.first_run_num > r.last_run_num) {
      spdlog::error("BeamRunCondition: invalid range {}-{}",
                    r.first_run_num, r.last_run_num);
      *ok = false;
      continue;
    }
    if (!IsKnownParticle(r.particle)) {
      spdlog::error("BeamRunCondition: invalid particle in range {}-{}",
                    r.first_run_num, r.last_run_num);
      *ok = false;
      continue;
    }
    if (!(r.momentum_gev_c > 0.) || !std::isfinite(r.momentum_gev_c)) {
      spdlog::error("BeamRunCondition: invalid momentum in range {}-{} ({})",
                    r.first_run_num, r.last_run_num, r.momentum_gev_c);
      *ok = false;
      continue;
    }
    for (UInt_t n = r.first_run_num; n <= r.last_run_num; ++n) {
      if (map.find(n) != map.end()) {
        spdlog::error("BeamRunCondition: overlapping run_num={}", n);
        *ok = false;
        continue;
      }
      BeamRunCondition c;
      c.run_num = n;
      c.particle = r.particle;
      c.momentum_gev_c = r.momentum_gev_c;
      map.emplace(n, c);
    }
  }
  return map;
}

const std::unordered_map<UInt_t, BeamRunCondition>&
Table()
{
  static Bool_t ok = true;
  static const auto table = BuildTable(&ok);
  static const Bool_t logged = [&]() {
    if (!ok)
      spdlog::error("BeamRunCondition: table validation failed at init");
    else
      spdlog::info("BeamRunCondition: loaded {} run entries", table.size());
    return true;
  }();
  (void)logged;
  return table;
}

void
ApplyUserParamOverrides(ResolvedBeamCondition& out)
{
  const auto& gUser = UserParamMan::GetInstance();
  if (gUser.Has("BeamMom")) {
    const Double_t mom = gUser.Get("BeamMom");
    if (mom > 0. && std::isfinite(mom)) {
      out.momentum_gev_c = mom;
      out.momentum_source = BeamConditionSource::UserParam;
    } else {
      spdlog::error("BeamRunCondition: run_num={} invalid BeamMom={}",
                    out.run_num, mom);
    }
  }
  if (gUser.Has("BeamParticle")) {
    const Double_t raw = gUser.Get("BeamParticle");
    if (!std::isfinite(raw) || raw != std::floor(raw)) {
      spdlog::error("BeamRunCondition: run_num={} BeamParticle={} "
                    "must be integer 0/1/2",
                    out.run_num, raw);
    } else {
      const Int_t code = static_cast<Int_t>(raw);
      if (code == 0 || code == 1 || code == 2) {
        out.particle = static_cast<BeamSpecies>(code);
        out.particle_source = BeamConditionSource::UserParam;
      } else {
        spdlog::error("BeamRunCondition: run_num={} BeamParticle={} "
                      "out of range (use 0=Pion, 1=Kaon, 2=Proton)",
                      out.run_num, code);
      }
    }
  }
}

// Result of the last ResolveBeamCondition() call. The run table and UserParam do not change
// during a job, so a repeated run_num (every track / event of a run) returns this copy.
struct LastResolved
{
  Bool_t valid = false;
  ResolvedBeamCondition value;
};
LastResolved g_last_resolved;

void
LogResolvedOnce(const ResolvedBeamCondition& r)
{
  static std::unordered_set<UInt_t> logged;
  if (!logged.insert(r.run_num).second)
    return;

  const Bool_t any = r.ValidParticle() || r.ValidMomentum()
      || r.particle_source != BeamConditionSource::None
      || r.momentum_source != BeamConditionSource::None;
  if (!any) {
    spdlog::warn("BeamRunCondition: run_num={} unresolved "
                 "(no run-table entry, no BeamParticle / BeamMom)",
                 r.run_num);
    return;
  }
  if (r.ValidMomentum()) {
    spdlog::info("BeamRunCondition: run_num={} particle={} momentum={:.3f} GeV/c "
                 "particle_source={} momentum_source={}",
                 r.run_num, ToString(r.particle), *r.momentum_gev_c,
                 ToString(r.particle_source), ToString(r.momentum_source));
  } else {
    spdlog::info("BeamRunCondition: run_num={} particle={} momentum=none "
                 "particle_source={} momentum_source={}",
                 r.run_num, ToString(r.particle),
                 ToString(r.particle_source), ToString(r.momentum_source));
  }
}

} // namespace

//_____________________________________________________________________________
const char*
ToString(BeamSpecies species)
{
  switch (species) {
    case BeamSpecies::Pion:    return "Pion";
    case BeamSpecies::Kaon:    return "Kaon";
    case BeamSpecies::Proton:  return "Proton";
    case BeamSpecies::Unknown: return "Unknown";
  }
  return "Invalid";
}

//_____________________________________________________________________________
const char*
ToString(BeamConditionSource source)
{
  switch (source) {
    case BeamConditionSource::None:              return "None";
    case BeamConditionSource::RunTable:          return "RunTable";
    case BeamConditionSource::UserParam:         return "UserParam";
  }
  return "Invalid";
}

//_____________________________________________________________________________
std::optional<BeamRunCondition>
FindBeamRunCondition(UInt_t run_num)
{
  const auto& table = Table();
  const auto it = table.find(run_num);
  if (it == table.end())
    return std::nullopt;
  return it->second;
}

//_____________________________________________________________________________
std::optional<Double_t>
BeamMassGeV(BeamSpecies species)
{
  switch (species) {
    case BeamSpecies::Pion:    return pdg::PionMass();
    case BeamSpecies::Kaon:    return pdg::KaonMass();
    case BeamSpecies::Proton:  return pdg::ProtonMass();
    case BeamSpecies::Unknown: return std::nullopt;
  }
  return std::nullopt;
}

//_____________________________________________________________________________
// Priority (strongest first): 1) BeamMom / BeamParticle (UserParam)  2) run table.
// A run without a table entry stays unresolved unless UserParam gives the missing item.
ResolvedBeamCondition
ResolveBeamCondition(UInt_t run_num)
{
  if (g_last_resolved.valid && g_last_resolved.value.run_num == run_num)
    return g_last_resolved.value;

  ResolvedBeamCondition out;
  out.run_num = run_num;
  out.particle = BeamSpecies::Unknown;
  out.particle_source = BeamConditionSource::None;
  out.momentum_source = BeamConditionSource::None;

  if (const auto hit = FindBeamRunCondition(run_num)) {
    out.particle = hit->particle;
    out.particle_source = BeamConditionSource::RunTable;
    out.momentum_gev_c = hit->momentum_gev_c;
    out.momentum_source = BeamConditionSource::RunTable;
  }

  ApplyUserParamOverrides(out);
  LogResolvedOnce(out);
  g_last_resolved.valid = true;
  g_last_resolved.value = out;
  return out;
}

//_____________________________________________________________________________
Bool_t
ValidateBeamRunConditionTable()
{
  Bool_t ok = true;
  const auto table = BuildTable(&ok);
  if (!ok)
    return false;
  for (const auto& kv : table) {
    if (kv.first != kv.second.run_num) {
      spdlog::error("BeamRunCondition: key/run_num mismatch {} vs {}",
                    kv.first, kv.second.run_num);
      ok = false;
    }
    if (!IsKnownParticle(kv.second.particle)) {
      spdlog::error("BeamRunCondition: bad particle at run_num={}", kv.first);
      ok = false;
    }
    if (!kv.second.momentum_gev_c
        || !(*kv.second.momentum_gev_c > 0.)
        || !std::isfinite(*kv.second.momentum_gev_c)) {
      spdlog::error("BeamRunCondition: bad momentum at run_num={}", kv.first);
      ok = false;
    }
  }
  return ok;
}
