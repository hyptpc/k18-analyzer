// -*- C++ -*-

#ifndef LH2_ELOSS_MAN_HH
#define LH2_ELOSS_MAN_HH

#include <string>
#include <vector>

#include <TMath.h>
#include <TString.h>
#include <TVector3.h>

//_____________________________________________________________________________
// LH2 momentum-loss map. Values are numbers only; keys label them
// ('#' starts a comment).
//   - Geometry: one "<key> <value>" line per Geometry field, any order, all
//     required (TargetRadius ... HolderBottomY) [mm; HolderWindowHalfPhi deg].
//   - Grid blocks, each
//       particle <EParticle>
//       region   <ERegion>
//       exit     <EExit>
//       holder   <EHolder>
//       p axis [MeV/c]
//       L axis [mm]
//       cos axis
//       np*nL rows: p L loss(cos_0) ... loss(cos_n) [MeV/c]  (p outer, L inner)
// Cylinder axis is Y in the TPC frame (target centre from the geometry keys).
class Lh2ElossMan
{
public:
  enum EParticle { kProton = 0, kPion = 1, kKaon = 2, kNParticle = 3 };
  // Map regions 0..2; kOutside has no grid.
  enum ERegion {
    kVessel = 0,
    kGap = 1,
    kGapCross = 2,
    kNRegion = 3,
    kOutside = kNRegion
  };
  enum EExit { kSide = 0, kCap = 1, kNExit = 2 };
  enum EHolder { kWindow = 0, kG10 = 1, kNHolder = 2 };

  // Straight ray from a vertex: inputs needed by the momentum-loss map.
  struct Track
  {
    Int_t region = kOutside;
    Bool_t in_liquid = false;  // vertex inside the LH2 cylinder
    Double_t L = 0.;           // path length in LH2 along the ray [mm]
    Int_t exit = kSide;        // surface where the LH2 path ends (kSide if L=0)
    Int_t holder = -1;         // EHolder; -1 if region == kOutside
    Double_t cos_alpha = TMath::QuietNaN();  // incidence cosine on the holder
  };

  struct Geometry
  {
    Double_t target_radius = TMath::QuietNaN();
    Double_t target_half_height = TMath::QuietNaN();
    Double_t target_center_x = TMath::QuietNaN();
    Double_t target_center_y = TMath::QuietNaN();
    Double_t target_center_z = TMath::QuietNaN();
    Double_t vessel_radius = TMath::QuietNaN();
    Double_t vessel_half_height = TMath::QuietNaN();
    Double_t holder_inner_radius = TMath::QuietNaN();
    Double_t holder_radius = TMath::QuietNaN();
    Double_t holder_window_half_phi = TMath::QuietNaN();  // deg
    Double_t holder_window_half_height = TMath::QuietNaN();
    Double_t holder_bottom_y = TMath::QuietNaN();  // relative to target centre
  };

  static const TString& ClassName();
  static Lh2ElossMan& GetInstance();
  ~Lh2ElossMan();

private:
  Lh2ElossMan();
  Lh2ElossMan(const Lh2ElossMan&);
  Lh2ElossMan& operator=(const Lh2ElossMan&);

  struct Grid
  {
    std::vector<Double_t> p;    // MeV/c, strictly increasing
    std::vector<Double_t> L;    // mm
    std::vector<Double_t> cos;  // [0,1]
    // index = ip * (nL * ncos) + iL * ncos + ic  (cos fastest)
    std::vector<Double_t> loss;

    Bool_t Empty() const { return p.empty() || loss.empty(); }
    Double_t Eval(Double_t p_MeV, Double_t L_mm, Double_t c) const;
  };

  // Finite Y-axis cylinder in the target-centre frame.
  struct Cylinder
  {
    Double_t R;
    Double_t ylo;
    Double_t yhi;
  };

  // Parameter interval of a ray inside a cylinder (t >= 0).
  struct Segment
  {
    Bool_t hit = false;
    Double_t t0 = 0.;
    Double_t t1 = 0.;
    Bool_t side_ends = false;  // true if t1 is set by the side wall
  };

public:
  Bool_t Initialize();
  Bool_t Initialize(const TString& file_name);
  Bool_t IsReady() const { return m_is_ready; }
  Bool_t HasParticle(Int_t particle) const;
  const TString& FileName() const { return m_file_name; }

  // Geometry only (no map lookup). dir need not be normalized.
  Track Classify(const TVector3& vtx, const TVector3& dir) const;
  // Momentum loss [MeV/c] to add at the vertex; p_MeV = |p| measured by the TPC.
  // 0 for kOutside; NaN if the particle is not in the map or p_MeV is not > 0.
  Double_t Loss(const Track& trk, Double_t p_MeV, Int_t particle) const;

private:
  // Line reader for the map file (defined in .cc).
  struct MapReader;
  struct MapLine;

  static Int_t GridIndex(Int_t particle, Int_t region, Int_t exit, Int_t holder);
  static Int_t NGridSlots()
  {
    return kNParticle * kNRegion * kNExit * kNHolder;
  }
  Grid& GridAt(Int_t particle, Int_t region, Int_t exit, Int_t holder)
  {
    return m_grid[static_cast<std::size_t>(GridIndex(particle, region, exit, holder))];
  }
  const Grid& GridAt(Int_t particle, Int_t region, Int_t exit, Int_t holder) const
  {
    return m_grid[static_cast<std::size_t>(GridIndex(particle, region, exit, holder))];
  }

  static Segment RayInCylinder(const TVector3& x, const TVector3& u,
                               const Cylinder& c);
  static Bool_t ReadAxis(MapReader& in, const char* what,
                         std::vector<Double_t>& axis);
  static Bool_t ReadBody(MapReader& in, Grid& gr);

  Cylinder Liquid() const;
  Cylinder Vessel() const;
  Cylinder Holder() const;
  TVector3 TargetCenter() const;

  Int_t FindRegion(const TVector3& x, const TVector3& u) const;
  Bool_t InLiquid(const TVector3& x) const;
  void SetLiquidPath(const TVector3& x, const TVector3& u, Track& trk) const;
  void SetHolderExit(const TVector3& x, const TVector3& u, Track& trk) const;

  Bool_t ReadGeoKey(MapReader& in, const MapLine& ln);
  Bool_t ReadGrid(MapReader& in, const MapLine& first);
  Bool_t CheckGeometry(MapReader& in) const;
  Bool_t CheckGrids(MapReader& in) const;
  std::string MissingGrids(Int_t particle) const;
  Bool_t HasAnyGrid(Int_t particle) const;

  void Clear();

  Bool_t m_is_ready = false;
  TString m_file_name;
  Geometry m_geom;
  std::vector<Grid> m_grid;  // size NGridSlots(); empty Grid = absent
};

//_____________________________________________________________________________
inline const TString&
Lh2ElossMan::ClassName()
{
  static TString s_name("Lh2ElossMan");
  return s_name;
}

//_____________________________________________________________________________
inline Lh2ElossMan&
Lh2ElossMan::GetInstance()
{
  static Lh2ElossMan s_instance;
  return s_instance;
}

#endif
