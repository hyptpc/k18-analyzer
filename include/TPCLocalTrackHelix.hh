// -*- C++ -*-

#ifndef TPC_LOCAL_TRACK_HELIX_HH
#define TPC_LOCAL_TRACK_HELIX_HH

#include <vector>
#include <functional>

#include <std_ostream.hh>

#include "TVector3.h"
#include "ThreeVector.hh"
#include "DetectorID.hh"
#include "TPCHit.hh"
#include "TPCCluster.hh"
#include "TPCLTrackHit.hh"

class TPCHit;
class TPCCluster;

//______________________________________________________________________________
class TPCLocalTrackHelix
{
public:
  static const TString ClassName();
  explicit TPCLocalTrackHelix();
  ~TPCLocalTrackHelix();
  TPCLocalTrackHelix(TPCLocalTrackHelix *init); //deep copy

private:
  TPCLocalTrackHelix & operator =(const TPCLocalTrackHelix &init);

private:
  bool   m_is_fitted;     // flag of DoFit()
  bool   m_is_calculated; // flag of Calculate()
  int    m_hough_flag;
  std::vector<TPCLTrackHit*> m_hit_array;
  std::vector<int> m_hit_order;
  std::vector<double> m_hit_t;

  //equation of Helix
  //x = -X;
  //y = Z - Tgtz;
  //z = Y;
  //x = p[0] + p[3]*cos(t+theta0);
  //y = p[1] + p[3]*sin(t+theta0);
  //z = p[2] + p[4]*p[3]*(t+theta0);
  // track coordinate origin is target, ***NOT TPC center***
  //Hough params
  double m_Acx;
  double m_Acy;
  double m_Az0;
  double m_Ar;
  double m_Adz;
  //Track params
  double m_cx;
  double m_cy;
  double m_z0;
  double m_r;
  double m_dz;

  Double_t m_closedist; //closest distance from the track to the target
  double m_chisqr;
  double m_n_iteration;
  TVector3 m_mom0;
  TVector3 m_mom0_corP;
  TVector3 m_mom0_corN;
  double m_median_t;
  double m_min_t;
  double m_max_t;
  double m_path;
  double m_transverse_path;
  int m_charge;
  int m_fitflag;
  int m_vtxflag; //for track seperation around the target
  int m_isBeam;
  int m_isK18;
  int m_isKurama;
  int m_isAccidental;

  int m_searchtime; //millisec
  int m_fittime; //millisec

  std::vector<double> m_cx_exclusive;
  std::vector<double> m_cy_exclusive;
  std::vector<double> m_z0_exclusive;
  std::vector<double> m_r_exclusive;
  std::vector<double> m_dz_exclusive;
  std::vector<double> m_chisqr_exclusive;
  std::vector<TVector3> m_vp; //RK virtual plane

public:
  void         AddTPCHit(TPCLTrackHit *hit);
  void         EraseHits(std::vector<Int_t> delete_hits);
  void         ClearHits();
  void         Calculate();
  void         DeleteNullHit();
  //bool         DoHelixFit(int MinHits);
  bool         DoHelixFit();
  bool         DoFit(int MinHits=0);
  int          FinalizeTrack(int &delete_hit);
  int          FinalizeTrack();
  void         SortHitOrder(); //Sort hits by theta
  bool         ResidualCheck(TVector3 pos, TVector3 Res, double &resi);
  bool         ResidualCheck(TVector3 pos, double xzwindow, double ywindow, double &resi);
  bool         ResidualCheck(TVector3 pos, double xzwindow, double ywindow);
  int          Side(TVector3 pos);
  void         InvertChargeCheck();
  void         AddVPHit(TVector3 vp);
  bool         DoVPFit();
  void         Print(const TString& arg="", bool print_allhits = false) const;

  void     CalcClosestDist();
  TVector3 CalcHelixMom(double par[5], double y) const;
  TVector3 CalcHelixMom_t(double par[5], double t) const;
  TVector3 CalcHelixMom_corP(double par[5], double y) const;
  TVector3 CalcHelixMom_t_corP(double par[5], double t) const;
  TVector3 CalcHelixMom_corN(double par[5], double y) const;
  TVector3 CalcHelixMom_t_corN(double par[5], double t) const;

  TPCLTrackHit* GetHit(std::size_t nth) const;
  TPCLTrackHit* GetHitInOrder(std::size_t nth) const;
  int          GetNHit() const { return m_hit_array.size(); }
  int          GetNPad() const;
  int          GetOrder(int i) const { return m_hit_order[i]; }

  double   GetChiSquare() const { return m_chisqr; }
  Double_t GetClosestDist() const { return m_closedist; }
  int      GetNIteration() const { return m_n_iteration; }
  int      GetHoughFlag(void) const {return m_hough_flag; }
  int      GetFitFlag(void) const {return m_fitflag; }
  int      GetVtxFlag(void) const {return m_vtxflag;} //vtx is in the target or not
  int      GetIsBeam(void) const { return m_isBeam; }
  int      GetIsK18(void) const { return m_isK18; }
  int      GetIsKurama(void) const { return m_isKurama; }
  int      GetIsAccidental(void) const { return m_isAccidental; }
  bool     IsFitted() const { return m_is_fitted; }
  bool     IsCalculated() const { return m_is_calculated; }
  int      GetSearchTime() const { return m_searchtime; }
  int      GetFitTime() const { return m_fittime; }
  TVector3 GetMom0() const { return m_mom0; }// Momentum at Y = 0
  TVector3 GetMom0_corP() const { return m_mom0_corP; }// Momentum at Y = 0
  TVector3 GetMom0_corN() const { return m_mom0_corN; }// Momentum at Y = 0
  TVector3 GetPosition(double par[5], double t) const;
  int      GetNDF() const;
  double Getcx() const { return m_cx; }
  double Getcy() const { return m_cy; }
  double Getz0() const { return m_z0; }
  double Getr() const { return m_r; }
  double Getdz() const { return m_dz; }
  double GetAcx() const { return m_Acx; }
  double GetAcy() const { return m_Acy; }
  double GetAz0() const { return m_Az0; }
  double GetAr() const { return m_Ar; }
  double GetAdz() const { return m_Adz; }
  double GetTheta(int i) const { return m_hit_t[i]; }
  double GetTcal(TVector3 pos);
  double GetMint(void) const { return m_min_t; }
  double GetMaxt(void) const { return m_max_t; }
  double GetPath(void) const { return m_path; }
  double GetTransversePath(void) const { return m_transverse_path; }
  int    GetCharge(void) const { return m_charge; }
  double GetTrackdE();
  double GetdEdx(double truncatedMean = 1.0);

  void SetHoughFlag(int hough_flag);
  void SetFlag(int flag);
  void SetFitFlag(int flag) { m_fitflag = flag; }
  void SetVtxFlag(int flag) { m_vtxflag = flag; }
  void SetIsBeam(int flag=1) { m_isBeam = flag; }
  void SetIsK18(int flag=1) { m_isK18 = flag; }
  void SetIsKurama(int flag=1) { m_isKurama = flag; }
  void SetIsAccidental(int flag=1) { m_isAccidental = flag; }
  void SetSearchTime(int time) { m_searchtime = time; }
  void SetFitTime(int time) { m_fittime = time; }
  void SetAcx(double Acx) { m_Acx = Acx; }
  void SetAcy(double Acy) { m_Acy = Acy; }
  void SetAz0(double Az0) { m_Az0 = Az0; }
  void SetAr(double Ar) { m_Ar = Ar; }
  void SetAdz(double Adz) { m_Adz = Adz; }
  void SetHoughParam(double *par){ m_Acx = par[0]; m_Acy = par[1]; m_Az0 = par[2]; m_Ar = par[3]; m_Adz = par[4]; };
  void SetParamUsingHoughParam(void);

  //exclusive
  void CalculateExclusive();
  void DoFitExclusive();

};

//_____________________________________________________________________________
inline const TString
TPCLocalTrackHelix::ClassName()
{
  static TString s_name("TPCLocalTrackHelix");
  return s_name;
}

#endif
