// -*- C++ -*-

#include "Lh2ElossMan.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include <TMath.h>

#include "FuncName.hh"

#include <spdlog/spdlog.h>

namespace
{
constexpr Double_t kEps = 1.e-12;
constexpr Double_t kInf = std::numeric_limits<Double_t>::infinity();

const char* const kParticleName[] = {"proton", "pion", "kaon"};
const char* const kRegionName[] = {"vessel", "gap", "gap_cross"};
const char* const kExitName[] = {"side", "cap"};
const char* const kHolderName[] = {"window", "g10"};

static_assert(sizeof(kParticleName) / sizeof(kParticleName[0])
                == Lh2ElossMan::kNParticle,
              "kParticleName size");
static_assert(sizeof(kRegionName) / sizeof(kRegionName[0])
                == Lh2ElossMan::kNRegion,
              "kRegionName size");
static_assert(sizeof(kExitName) / sizeof(kExitName[0])
                == Lh2ElossMan::kNExit,
              "kExitName size");
static_assert(sizeof(kHolderName) / sizeof(kHolderName[0])
                == Lh2ElossMan::kNHolder,
              "kHolderName size");

// Ten grids required per particle in the map (gap has side only; L = 0).
const Int_t kRequiredGrid[][3] = {
  {Lh2ElossMan::kVessel, Lh2ElossMan::kSide, Lh2ElossMan::kWindow},
  {Lh2ElossMan::kVessel, Lh2ElossMan::kSide, Lh2ElossMan::kG10},
  {Lh2ElossMan::kVessel, Lh2ElossMan::kCap, Lh2ElossMan::kWindow},
  {Lh2ElossMan::kVessel, Lh2ElossMan::kCap, Lh2ElossMan::kG10},
  {Lh2ElossMan::kGap, Lh2ElossMan::kSide, Lh2ElossMan::kWindow},
  {Lh2ElossMan::kGap, Lh2ElossMan::kSide, Lh2ElossMan::kG10},
  {Lh2ElossMan::kGapCross, Lh2ElossMan::kSide, Lh2ElossMan::kWindow},
  {Lh2ElossMan::kGapCross, Lh2ElossMan::kSide, Lh2ElossMan::kG10},
  {Lh2ElossMan::kGapCross, Lh2ElossMan::kCap, Lh2ElossMan::kWindow},
  {Lh2ElossMan::kGapCross, Lh2ElossMan::kCap, Lh2ElossMan::kG10},
};
constexpr Int_t kNRequiredGrid
  = static_cast<Int_t>(sizeof(kRequiredGrid) / sizeof(kRequiredGrid[0]));

struct GeoKey
{
  const char* name;
  Double_t Lh2ElossMan::Geometry::* member;
  Bool_t positive;
};

const GeoKey kGeoKeys[] = {
  {"TargetRadius", &Lh2ElossMan::Geometry::target_radius, true},
  {"TargetHalfHeight", &Lh2ElossMan::Geometry::target_half_height, true},
  {"TargetCenterX", &Lh2ElossMan::Geometry::target_center_x, false},
  {"TargetCenterY", &Lh2ElossMan::Geometry::target_center_y, false},
  {"TargetCenterZ", &Lh2ElossMan::Geometry::target_center_z, false},
  {"VesselRadius", &Lh2ElossMan::Geometry::vessel_radius, true},
  {"VesselHalfHeight", &Lh2ElossMan::Geometry::vessel_half_height, true},
  {"HolderInnerRadius", &Lh2ElossMan::Geometry::holder_inner_radius, true},
  {"HolderRadius", &Lh2ElossMan::Geometry::holder_radius, true},
  {"HolderWindowHalfPhi", &Lh2ElossMan::Geometry::holder_window_half_phi, true},
  {"HolderWindowHalfHeight", &Lh2ElossMan::Geometry::holder_window_half_height, true},
  {"HolderBottomY", &Lh2ElossMan::Geometry::holder_bottom_y, false},
};
constexpr Int_t kNGeoKey
  = static_cast<Int_t>(sizeof(kGeoKeys) / sizeof(kGeoKeys[0]));

//_____________________________________________________________________________
Int_t
IndexOf(
  const char* const* names,
  Int_t n,
  const std::string& s
)
{
  for (Int_t i = 0; i < n; ++i) {
    if (s == names[i])
      return i;
  }
  return -1;
}

//_____________________________________________________________________________
Bool_t
StrictlyIncreasing(
  const std::vector<Double_t>& v
)
{
  for (std::size_t i = 1; i < v.size(); ++i) {
    if (!(v[i] > v[i - 1]))
      return false;
  }
  return true;
}

//_____________________________________________________________________________
// Linear weights on a grid axis; clamps outside the range (single-node -> i0=i1).
void
AxisWeight(
  const std::vector<Double_t>& edges,
  Double_t x,
  Int_t& i0,
  Int_t& i1,
  Double_t& w
)
{
  const Int_t n = static_cast<Int_t>(edges.size());
  if (n <= 1 || x <= edges.front()) {
    i0 = i1 = 0;
    w = 0.;
    return;
  }
  if (x >= edges.back()) {
    i0 = i1 = n - 1;
    w = 0.;
    return;
  }
  const auto it = std::upper_bound(edges.begin(), edges.end(), x);
  i0 = static_cast<Int_t>(it - edges.begin()) - 1;
  i1 = i0 + 1;
  w = (x - edges[i0]) / (edges[i1] - edges[i0]);
}

}  // namespace

//_____________________________________________________________________________
struct Lh2ElossMan::MapReader
{
  std::ifstream ifs;
  TString file;
  Int_t nline = 0;

  Bool_t Next(std::istringstream& iss);
  // Error at the current line; always returns false.
  Bool_t Error(const TString& msg) const;
  // Error about the whole file (no line number); always returns false.
  Bool_t FileError(const TString& msg) const;
};

//_____________________________________________________________________________
// Skip blank / comment-only lines; set iss to the next line that has a token.
Bool_t
Lh2ElossMan::MapReader::Next(
  std::istringstream& iss
)
{
  std::string raw;
  while (std::getline(ifs, raw)) {
    ++nline;
    const auto hash = raw.find('#');
    if (hash != std::string::npos)
      raw.erase(hash);
    iss.clear();
    iss.str(raw);
    std::string tok;
    if (!(iss >> tok))
      continue;
    // Rewind so the caller can parse the full line from the first token.
    iss.clear();
    iss.str(raw);
    return true;
  }
  return false;
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::MapReader::Error(
  const TString& msg
) const
{
  spdlog::error("{}:{}: {}", file.Data(), nline, msg.Data());
  return false;
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::MapReader::FileError(
  const TString& msg
) const
{
  spdlog::error("{}: {}", file.Data(), msg.Data());
  return false;
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::ReadAxis(
  MapReader& in,
  const char* name,
  std::vector<Double_t>& axis
)
{
  std::istringstream iss;
  if (!in.Next(iss))
    return in.Error(Form("unexpected EOF reading %s axis", name));
  std::string label;
  if (!(iss >> label) || label != name)
    return in.Error(Form("expected '%s' axis", name));
  Double_t v = 0.;
  while (iss >> v) {
    if (!std::isfinite(v))
      return in.Error(Form("non-finite %s axis value", name));
    axis.push_back(v);
  }
  if (axis.empty() || !StrictlyIncreasing(axis))
    return in.Error(Form("invalid %s axis", name));
  return true;
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::ReadRow(
  MapReader& in,
  Double_t p,
  Double_t L,
  Int_t nc,
  Double_t* loss
)
{
  std::istringstream iss;
  if (!in.Next(iss))
    return in.Error("unexpected EOF in grid body");
  Double_t pv = 0.;
  Double_t Lv = 0.;
  if (!(iss >> pv >> Lv))
    return in.Error("bad grid row");
  if (pv != p || Lv != L)
    return in.Error(Form("p/L mismatch (got %g %g, expected %g %g)",
                         pv, Lv, p, L));
  for (Int_t ic = 0; ic < nc; ++ic) {
    Double_t val = 0.;
    if (!(iss >> val) || !std::isfinite(val) || val < 0.)
      return in.Error("bad loss value");
    loss[ic] = val;
  }
  Double_t junk = 0.;
  if (iss >> junk)
    return in.Error("extra tokens in grid row");
  return true;
}

//_____________________________________________________________________________
Lh2ElossMan::Lh2ElossMan()
{
  Clear();
}

//_____________________________________________________________________________
Lh2ElossMan::~Lh2ElossMan()
{
}

//_____________________________________________________________________________
void
Lh2ElossMan::Clear()
{
  m_is_ready = false;
  m_geom = Geometry{};
  m_grid.assign(static_cast<std::size_t>(NGridSlots()), Grid{});
}

//_____________________________________________________________________________
Int_t
Lh2ElossMan::GridIndex(
  Int_t particle,
  Int_t region,
  Int_t exit,
  Int_t holder
)
{
  return ((particle * kNRegion + region) * kNExit + exit) * kNHolder + holder;
}

//_____________________________________________________________________________
std::string
Lh2ElossMan::MissingGrids(
  Int_t particle
) const
{
  if (particle < 0 || particle >= kNParticle)
    return "invalid particle";
  std::string miss;
  for (Int_t i = 0; i < kNRequiredGrid; ++i) {
    const Int_t rg = kRequiredGrid[i][0];
    const Int_t ex = kRequiredGrid[i][1];
    const Int_t ho = kRequiredGrid[i][2];
    if (!GridAt(particle, rg, ex, ho).Empty())
      continue;
    if (!miss.empty())
      miss += ", ";
    miss += kRegionName[rg];
    miss += ' ';
    miss += kExitName[ex];
    miss += ' ';
    miss += kHolderName[ho];
  }
  return miss;
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::HasAnyGrid(
  Int_t particle
) const
{
  for (Int_t i = 0; i < kNRequiredGrid; ++i) {
    const Int_t* g = kRequiredGrid[i];
    if (!GridAt(particle, g[0], g[1], g[2]).Empty())
      return true;
  }
  return false;
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::HasParticle(
  Int_t particle
) const
{
  return m_is_ready && MissingGrids(particle).empty();
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::Initialize(
  const TString& file_name
)
{
  m_file_name = file_name;
  return Initialize();
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::Initialize()
{
  Clear();
  if (m_file_name.IsNull()) {
    spdlog::error("{} empty file name", FUNC_NAME.Data());
    return false;
  }

  MapReader in;
  in.ifs.open(m_file_name.Data());
  in.file = m_file_name;
  if (!in.ifs.is_open())
    return in.Error("cannot open");

  std::istringstream iss;
  while (in.Next(iss)) {
    std::string key;
    iss >> key;
    const Bool_t ok
      = (key == "grid") ? ReadGrid(in, iss) : ReadGeoKey(in, key, iss);
    if (!ok)
      return false;
  }
  if (!CheckGeometry(in) || !CheckGrids(in))
    return false;

  m_is_ready = true;
  spdlog::info("{} loaded {}", FUNC_NAME.Data(), m_file_name.Data());
  return true;
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::ReadGeoKey(
  MapReader& in,
  const std::string& key,
  std::istringstream& iss
)
{
  const GeoKey* gkey = nullptr;
  for (Int_t i = 0; i < kNGeoKey; ++i) {
    if (key == kGeoKeys[i].name) {
      gkey = &kGeoKeys[i];
      break;
    }
  }
  if (!gkey)
    return in.Error(Form("unknown key '%s'", key.c_str()));
  if (!std::isnan(m_geom.*(gkey->member)))
    return in.Error(Form("duplicate key '%s'", key.c_str()));

  Double_t val = 0.;
  if (!(iss >> val) || !std::isfinite(val))
    return in.Error(Form("invalid value for '%s'", key.c_str()));
  std::string extra;
  if (iss >> extra)
    return in.Error(Form("extra tokens for '%s'", key.c_str()));
  if (gkey->positive && !(val > 0.))
    return in.Error(Form("'%s' must be > 0", key.c_str()));

  m_geom.*(gkey->member) = val;
  return true;
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::ReadGrid(
  MapReader& in,
  std::istringstream& iss
)
{
  std::string s_pa, s_rg, s_ex, s_ho;
  if (!(iss >> s_pa >> s_rg >> s_ex >> s_ho))
    return in.Error("bad grid header");
  std::string extra;
  if (iss >> extra)
    return in.Error("extra tokens in grid header");

  const Int_t pa = IndexOf(kParticleName, kNParticle, s_pa);
  const Int_t rg = IndexOf(kRegionName, kNRegion, s_rg);
  const Int_t ex = IndexOf(kExitName, kNExit, s_ex);
  const Int_t ho = IndexOf(kHolderName, kNHolder, s_ho);
  if (pa < 0 || rg < 0 || ex < 0 || ho < 0)
    return in.Error(Form("unknown grid name (%s %s %s %s)",
                         s_pa.c_str(), s_rg.c_str(), s_ex.c_str(),
                         s_ho.c_str()));
  // gap vertices never exit through a liquid surface; only side grids exist.
  if (rg == kGap && ex != kSide)
    return in.Error("gap grid must use exit side");

  Grid& slot = GridAt(pa, rg, ex, ho);
  if (!slot.Empty())
    return in.Error("duplicate grid");

  Grid gr;
  if (!ReadAxis(in, "p", gr.p) || !ReadAxis(in, "L", gr.L)
      || !ReadAxis(in, "cos", gr.cos))
    return false;

  const Int_t np = static_cast<Int_t>(gr.p.size());
  const Int_t nL = static_cast<Int_t>(gr.L.size());
  const Int_t nc = static_cast<Int_t>(gr.cos.size());
  gr.loss.assign(static_cast<std::size_t>(np * nL * nc), 0.);

  // Row order is p outer, L inner; each row must match the axis values exactly.
  for (Int_t ip = 0; ip < np; ++ip) {
    for (Int_t iL = 0; iL < nL; ++iL) {
      const std::size_t base
        = static_cast<std::size_t>((ip * nL + iL) * nc);
      if (!ReadRow(in, gr.p[static_cast<std::size_t>(ip)],
                   gr.L[static_cast<std::size_t>(iL)], nc,
                   &gr.loss[base]))
        return false;
    }
  }
  slot = std::move(gr);
  return true;
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::CheckGeometry(
  MapReader& in
) const
{
  for (Int_t i = 0; i < kNGeoKey; ++i) {
    if (std::isnan(m_geom.*(kGeoKeys[i].member)))
      return in.FileError(Form("missing key '%s'", kGeoKeys[i].name));
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::CheckGrids(
  MapReader& in
) const
{
  // Proton is mandatory; other particles only if at least one of their grids appears.
  for (Int_t pa = 0; pa < kNParticle; ++pa) {
    if (pa != kProton && !HasAnyGrid(pa))
      continue;
    const std::string miss = MissingGrids(pa);
    if (!miss.empty())
      return in.FileError(Form("%s grids incomplete: missing %s",
                               kParticleName[pa], miss.c_str()));
  }
  return true;
}

//_____________________________________________________________________________
Lh2ElossMan::Cylinder
Lh2ElossMan::Liquid() const
{
  return {m_geom.target_radius, -m_geom.target_half_height,
          m_geom.target_half_height};
}

//_____________________________________________________________________________
Lh2ElossMan::Cylinder
Lh2ElossMan::Vessel() const
{
  return {m_geom.vessel_radius, -m_geom.vessel_half_height,
          m_geom.vessel_half_height};
}

//_____________________________________________________________________________
// Open top (yhi = +inf): rays can leave through the free upper end.
Lh2ElossMan::Cylinder
Lh2ElossMan::Holder() const
{
  return {m_geom.holder_radius, m_geom.holder_bottom_y, kInf};
}

//_____________________________________________________________________________
TVector3
Lh2ElossMan::TargetCenter() const
{
  return {m_geom.target_center_x, m_geom.target_center_y,
          m_geom.target_center_z};
}

//_____________________________________________________________________________
// Overlap of the ray x+t*u (t>=0) with r_xz < R and ylo < y < yhi.
Lh2ElossMan::Segment
Lh2ElossMan::RayInCylinder(
  const TVector3& x,
  const TVector3& u,
  const Cylinder& cyl
)
{
  Segment out;
  // Radial interval from the quadratic in the XZ plane.
  const Double_t a = u.X() * u.X() + u.Z() * u.Z();
  const Double_t c = x.X() * x.X() + x.Z() * x.Z() - cyl.R * cyl.R;
  Double_t r_lo = 0.;
  Double_t r_hi = 0.;
  if (a > kEps) {
    const Double_t b = x.X() * u.X() + x.Z() * u.Z();
    const Double_t disc = b * b - a * c;
    if (disc <= 0.)
      return out;
    const Double_t s = std::sqrt(disc);
    r_lo = (-b - s) / a;
    r_hi = (-b + s) / a;
  } else {
    // Ray parallel to the Y axis: either forever inside or never.
    if (c >= 0.)
      return out;
    r_lo = -kInf;
    r_hi = kInf;
  }

  Double_t y_lo = 0.;
  Double_t y_hi = 0.;
  if (std::abs(u.Y()) > kEps) {
    const Double_t ta = (cyl.ylo - x.Y()) / u.Y();
    const Double_t tb = (cyl.yhi - x.Y()) / u.Y();
    y_lo = std::min(ta, tb);
    y_hi = std::max(ta, tb);
  } else {
    if (!(cyl.ylo < x.Y() && x.Y() < cyl.yhi))
      return out;
    y_lo = -kInf;
    y_hi = kInf;
  }

  out.t0 = std::max({r_lo, y_lo, 0.});
  out.t1 = std::min(r_hi, y_hi);
  if (!(out.t1 > out.t0 + kEps))
    return out;
  // Which bound sets t1: side wall vs end plane (cap / bottom).
  out.side_ends = (r_hi <= y_hi);
  out.hit = true;
  return out;
}

//_____________________________________________________________________________
// Strict inequalities on the boundaries (same as the map / holderref2).
Int_t
Lh2ElossMan::FindRegion(
  const TVector3& x,
  const TVector3& u
) const
{
  const Double_t r = std::hypot(x.X(), x.Z());
  const Bool_t in_vessel
    = (r < m_geom.vessel_radius
       && std::abs(x.Y()) < m_geom.vessel_half_height);
  if (in_vessel)
    return kVessel;

  const Bool_t in_holder
    = (r < m_geom.holder_inner_radius
       && x.Y() > m_geom.holder_bottom_y);
  if (!in_holder)
    return kOutside;

  // Between vessel and holder: gap_cross if the ray later enters the vessel.
  const Segment s = RayInCylinder(x, u, Vessel());
  return s.hit ? kGapCross : kGap;
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::InLiquid(
  const TVector3& x
) const
{
  const Double_t r = std::hypot(x.X(), x.Z());
  return r < m_geom.target_radius
    && std::abs(x.Y()) < m_geom.target_half_height;
}

//_____________________________________________________________________________
void
Lh2ElossMan::SetLiquidPath(
  const TVector3& x,
  const TVector3& u,
  Track& trk
) const
{
  const Segment s = RayInCylinder(x, u, Liquid());
  if (!s.hit)
    return;
  trk.L = s.t1 - s.t0;
  trk.exit = s.side_ends ? kSide : kCap;
}

//_____________________________________________________________________________
// Where the ray leaves the holder (vertex is already inside; t0 = 0).
void
Lh2ElossMan::SetHolderExit(
  const TVector3& x,
  const TVector3& u,
  Track& trk
) const
{
  const Segment s = RayInCylinder(x, u, Holder());
  // Straight up the open top: treat as G10 with cos = 0.
  if (!s.hit || !std::isfinite(s.t1)) {
    trk.holder = kG10;
    trk.cos_alpha = 0.;
    return;
  }
  if (!s.side_ends) {
    trk.holder = kG10;
    trk.cos_alpha = std::abs(u.Y());
    return;
  }

  // Side wall: window if |phi| is within HalfPhi of 0 or 180 deg (beam axis).
  const TVector3 h = x + s.t1 * u;
  const Double_t phi = TMath::ATan2(h.X(), h.Z()) * TMath::RadToDeg();
  const Double_t hp = m_geom.holder_window_half_phi;
  const Bool_t win
    = ((std::abs(phi) < hp || std::abs(phi) > 180. - hp)
       && std::abs(h.Y()) < m_geom.holder_window_half_height);
  trk.holder = win ? kWindow : kG10;
  // Incidence cosine on the cylindrical wall, clamped to [0, 1].
  const Double_t c
    = (h.X() * u.X() + h.Z() * u.Z()) / m_geom.holder_radius;
  trk.cos_alpha = std::min(std::max(c, 0.), 1.);
}

//_____________________________________________________________________________
Lh2ElossMan::Track
Lh2ElossMan::Classify(
  const TVector3& vtx,
  const TVector3& dir
) const
{
  Track trk;
  if (!(dir.Mag() > 0.))
    return trk;
  const TVector3 x = vtx - TargetCenter();
  const TVector3 u = dir.Unit();

  trk.region = FindRegion(x, u);
  if (trk.region == kOutside)
    return trk;
  trk.in_liquid = InLiquid(x);
  // gap never crosses liquid; L stays 0 and exit stays kSide.
  if (trk.region != kGap)
    SetLiquidPath(x, u, trk);
  SetHolderExit(x, u, trk);
  return trk;
}

//_____________________________________________________________________________
Double_t
Lh2ElossMan::Grid::Eval(
  Double_t p_MeV,
  Double_t L_mm,
  Double_t c
) const
{
  Int_t ip0 = 0, ip1 = 0, iL0 = 0, iL1 = 0, ic0 = 0, ic1 = 0;
  Double_t wp = 0., wL = 0., wc = 0.;
  AxisWeight(p, p_MeV, ip0, ip1, wp);
  AxisWeight(L, L_mm, iL0, iL1, wL);
  AxisWeight(cos, c, ic0, ic1, wc);

  const Int_t nL = static_cast<Int_t>(L.size());
  const Int_t nc = static_cast<Int_t>(cos.size());
  auto at = [&](Int_t ip, Int_t iL, Int_t ic) -> Double_t {
    return loss[static_cast<std::size_t>((ip * nL + iL) * nc + ic)];
  };

  Double_t s = 0.;
  const Double_t ap[2] = {1. - wp, wp};
  const Double_t aL[2] = {1. - wL, wL};
  const Double_t ac[2] = {1. - wc, wc};
  const Int_t ips[2] = {ip0, ip1};
  const Int_t iLs[2] = {iL0, iL1};
  const Int_t ics[2] = {ic0, ic1};
  for (Int_t i = 0; i < 2; ++i) {
    for (Int_t j = 0; j < 2; ++j) {
      for (Int_t k = 0; k < 2; ++k)
        s += ap[i] * aL[j] * ac[k] * at(ips[i], iLs[j], ics[k]);
    }
  }
  return std::max(s, 0.);
}

//_____________________________________________________________________________
Double_t
Lh2ElossMan::Loss(
  const Track& trk,
  Double_t p_MeV,
  Int_t particle
) const
{
  if (!m_is_ready)
    return 0.;
  if (trk.region == kOutside)
    return 0.;
  if (!(p_MeV > 0.) || particle < 0 || particle >= kNParticle
      || trk.holder < 0 || trk.holder >= kNHolder
      || trk.exit < 0 || trk.exit >= kNExit
      || trk.region < 0 || trk.region >= kNRegion)
    return TMath::QuietNaN();

  const Grid& gr = GridAt(particle, trk.region, trk.exit, trk.holder);
  if (gr.Empty())
    return TMath::QuietNaN();
  return gr.Eval(p_MeV, trk.L, trk.cos_alpha);
}
