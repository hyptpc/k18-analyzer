// -*- C++ -*-

#include "Lh2ElossMan.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <TMath.h>

#include "FuncName.hh"

#include <spdlog/spdlog.h>

namespace
{
constexpr Double_t kEps = 1.e-12;
constexpr Double_t kInf = std::numeric_limits<Double_t>::infinity();

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

// Geometry keys ("<key> <value>" lines, any order, each exactly once).
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
const GeoKey*
FindGeoKey(
  const std::string& key
)
{
  for (const GeoKey& g : kGeoKeys) {
    if (key == g.name)
      return &g;
  }
  return nullptr;
}

//_____________________________________________________________________________
// Integer enum value in [0, n); -1 otherwise.
Int_t
ToTag(
  Double_t v,
  Int_t n
)
{
  if (!(v >= 0. && v < n) || v != std::floor(v))
    return -1;
  return static_cast<Int_t>(v);
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
Bool_t
AllFinite(
  const std::vector<Double_t>& v
)
{
  return std::all_of(v.begin(), v.end(),
                     [](Double_t x) { return std::isfinite(x); });
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
// Reads the map line by line. A line is an optional leading key followed by
// numbers; '#' starts a comment.
struct Lh2ElossMan::MapLine
{
  std::string key;  // empty if the line starts with a number
  std::vector<Double_t> v;
};

//_____________________________________________________________________________
struct Lh2ElossMan::MapReader
{
  enum EStatus { kLine, kEof, kBadToken };
  using Line = MapLine;

  std::ifstream ifs;
  TString file;
  Int_t nline = 0;

  EStatus Next(Line& ln);
  // Next line that must exist, have no key and be numeric.
  Bool_t ExpectNumbers(std::vector<Double_t>& v, const char* what);
  // Next line that must be "<key> <one value>".
  Bool_t ExpectKey(const char* key, Double_t& value);
  // Error at the current line; always returns false.
  Bool_t Error(const TString& msg) const;
  // Error about the whole file (no line number); always returns false.
  Bool_t FileError(const TString& msg) const;
};

//_____________________________________________________________________________
Lh2ElossMan::MapReader::EStatus
Lh2ElossMan::MapReader::Next(
  Line& ln
)
{
  std::string raw;
  while (std::getline(ifs, raw)) {
    ++nline;
    raw.erase(std::min(raw.find('#'), raw.size()));
    std::istringstream iss(raw);
    std::string first;
    if (!(iss >> first))
      continue;
    ln.key.clear();
    ln.v.clear();
    std::istringstream num(first);
    Double_t x = 0.;
    if (num >> x && num.eof())
      ln.v.push_back(x);
    else
      ln.key = first;
    while (iss >> x)
      ln.v.push_back(x);
    return iss.eof() ? kLine : kBadToken;
  }
  return kEof;
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::MapReader::ExpectNumbers(
  std::vector<Double_t>& v,
  const char* what
)
{
  Line ln;
  const EStatus st = Next(ln);
  if (st == kEof)
    return Error(Form("unexpected end of file reading %s", what));
  if (st == kBadToken || !ln.key.empty())
    return Error(Form("non-numeric token in %s", what));
  v = std::move(ln.v);
  return true;
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::MapReader::ExpectKey(
  const char* key,
  Double_t& value
)
{
  Line ln;
  if (Next(ln) != kLine || ln.key != key || ln.v.size() != 1)
    return Error(Form("expected '%s <value>'", key));
  value = ln.v[0];
  return true;
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
// Missing required grids as "region exit holder" triplets; empty if complete.
std::string
Lh2ElossMan::MissingGrids(
  Int_t particle
) const
{
  if (particle < 0 || particle >= kNParticle)
    return "invalid particle";
  std::string miss;
  for (Int_t i = 0; i < kNRequiredGrid; ++i) {
    const Int_t* g = kRequiredGrid[i];
    if (!GridAt(particle, g[0], g[1], g[2]).Empty())
      continue;
    if (!miss.empty())
      miss += ", ";
    miss += Form("%d %d %d", g[0], g[1], g[2]);
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
// Top level: geometry key lines and grid blocks (each starts with "particle").
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

  MapLine ln;
  for (;;) {
    const auto st = in.Next(ln);
    if (st == MapReader::kEof)
      break;
    if (st == MapReader::kBadToken)
      return in.Error("non-numeric value");
    const Bool_t ok
      = (ln.key == "particle") ? ReadGrid(in, ln) : ReadGeoKey(in, ln);
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
  const MapLine& ln
)
{
  if (ln.key.empty())
    return in.Error("number line outside a grid block");
  const GeoKey* g = FindGeoKey(ln.key);
  if (!g)
    return in.Error(Form("unknown key '%s'", ln.key.c_str()));
  if (!std::isnan(m_geom.*(g->member)))
    return in.Error(Form("duplicate key '%s'", g->name));
  if (ln.v.size() != 1 || !std::isfinite(ln.v[0]))
    return in.Error(Form("'%s' needs one finite value", g->name));
  if (g->positive && !(ln.v[0] > 0.))
    return in.Error(Form("'%s' must be > 0", g->name));
  m_geom.*(g->member) = ln.v[0];
  return true;
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::CheckGeometry(
  MapReader& in
) const
{
  for (const GeoKey& g : kGeoKeys) {
    if (std::isnan(m_geom.*(g.member)))
      return in.FileError(Form("missing key '%s'", g.name));
  }
  // Nested cylinders: liquid inside vessel inside holder.
  const Geometry& g = m_geom;
  if (!(g.target_radius <= g.vessel_radius
        && g.vessel_radius < g.holder_inner_radius
        && g.holder_inner_radius < g.holder_radius
        && g.target_half_height <= g.vessel_half_height))
    return in.FileError("geometry is not nested (liquid <= vessel < holder)");
  return true;
}

//_____________________________________________________________________________
// Block: particle / region / exit / holder lines, then axes and rows.
Bool_t
Lh2ElossMan::ReadGrid(
  MapReader& in,
  const MapLine& first
)
{
  if (first.v.size() != 1)
    return in.Error("expected 'particle <value>'");
  Double_t t[4] = {first.v[0], 0., 0., 0.};
  if (!in.ExpectKey("region", t[1]) || !in.ExpectKey("exit", t[2])
      || !in.ExpectKey("holder", t[3]))
    return false;
  const Int_t pa = ToTag(t[0], kNParticle);
  const Int_t rg = ToTag(t[1], kNRegion);
  const Int_t ex = ToTag(t[2], kNExit);
  const Int_t ho = ToTag(t[3], kNHolder);
  if (pa < 0 || rg < 0 || ex < 0 || ho < 0)
    return in.Error("grid label out of range"
                    " (particle 0-2 region 0-2 exit 0-1 holder 0-1)");
  // gap vertices never exit through a liquid surface; only side grids exist.
  if (rg == kGap && ex != kSide)
    return in.Error("gap grid must use exit 0 (side)");
  Grid& slot = GridAt(pa, rg, ex, ho);
  if (!slot.Empty())
    return in.Error("duplicate grid");

  Grid gr;
  if (!ReadAxis(in, "p axis", gr.p) || !ReadAxis(in, "L axis", gr.L)
      || !ReadAxis(in, "cos axis", gr.cos))
    return false;
  if (!ReadBody(in, gr))
    return false;
  slot = std::move(gr);
  return true;
}

//_____________________________________________________________________________
Bool_t
Lh2ElossMan::ReadAxis(
  MapReader& in,
  const char* what,
  std::vector<Double_t>& axis
)
{
  if (!in.ExpectNumbers(axis, what))
    return false;
  if (!AllFinite(axis) || !StrictlyIncreasing(axis))
    return in.Error(Form("%s must be finite and strictly increasing", what));
  return true;
}

//_____________________________________________________________________________
// Rows "p L loss(cos_0) ... loss(cos_n)", p outer and L inner; p and L must
// repeat the axis values exactly.
Bool_t
Lh2ElossMan::ReadBody(
  MapReader& in,
  Grid& gr
)
{
  const std::size_t nL = gr.L.size();
  const std::size_t nc = gr.cos.size();
  gr.loss.clear();
  gr.loss.reserve(gr.p.size() * nL * nc);
  std::vector<Double_t> row;
  for (const Double_t p : gr.p) {
    for (const Double_t L : gr.L) {
      if (!in.ExpectNumbers(row, "grid row"))
        return false;
      if (row.size() != 2 + nc)
        return in.Error(Form("grid row needs %zu values, got %zu",
                             2 + nc, row.size()));
      if (row[0] != p || row[1] != L)
        return in.Error(Form("p/L mismatch (got %g %g, expected %g %g)",
                             row[0], row[1], p, L));
      for (std::size_t ic = 0; ic < nc; ++ic) {
        const Double_t val = row[2 + ic];
        if (!std::isfinite(val) || val < 0.)
          return in.Error("loss must be finite and >= 0");
        gr.loss.push_back(val);
      }
    }
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
      return in.FileError(Form("particle %d grids incomplete,"
                               " missing (region exit holder): %s",
                               pa, miss.c_str()));
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
