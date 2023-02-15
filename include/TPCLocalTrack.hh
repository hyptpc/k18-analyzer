// -*- C++ -*-

#ifndef TPC_LOCAL_TRACK_HH
#define TPC_LOCAL_TRACK_HH

#include <vector>
#include <functional>

#include <TH2D.h>
#include <TMinuit.h>
#include <TString.h>
#include <TVector3.h>

#include <std_ostream.hh>

#include "DetectorID.hh"
#include "TPCCluster.hh"
#include "TPCHit.hh"
#include "TPCLTrackHit.hh"

class TPCHit;
class TPCCluster;

//_____________________________________________________________________________
class TPCLocalTrack
{
public:
  static const TString& ClassName();
  TPCLocalTrack();
  ~TPCLocalTrack();

private:
  TPCLocalTrack(const TPCLocalTrack &);
  TPCLocalTrack & operator =(const TPCLocalTrack &);

private:
  Bool_t   m_is_fitted;     // flag of DoFit()
  Bool_t   m_is_calculated; // flag of Calculate()
  int      m_fitflag;
  std::vector<TPCLTrackHit*> m_hit_array;

  // track coordinate origin is target, ***NOT TPC center***
  //Hough params
  Double_t m_Ax;
  Double_t m_Ay;
  Double_t m_Au;
  Double_t m_Av;
  //Track params
  Double_t m_x0;
  Double_t m_y0;
  Double_t m_u0;
  Double_t m_v0;

  Double_t m_chisqr;
  Double_t m_n_iteration;
  int m_searchtime; //millisec
  int m_fittime; //millisec

  std::vector<Double_t> m_x0_exclusive;
  std::vector<Double_t> m_y0_exclusive;
  std::vector<Double_t> m_u0_exclusive;
  std::vector<Double_t> m_v0_exclusive;

public:

  void          AddTPCHit(TPCLTrackHit *hit);
  void          ClearHits();
  void          Calculate();
  void          CalcChisquare();
  void          DeleteNullHit();
  Bool_t        DoFit();
  Bool_t        DoFit(Int_t min_hits); // obsolete
  Bool_t        DoLinearFit(Int_t min_hits); // obsolete
  Int_t         GetNDF() const { return m_hit_array.size() - 4; }
  Int_t         GetNHit() const { return m_hit_array.size(); }
  TPCLTrackHit* GetHit(Int_t i) const { return m_hit_array.at(i); }
  Bool_t        IsFitted() const { return m_is_fitted; }
  Bool_t        IsCalculated() const { return m_is_calculated; }
  Bool_t        ResidualIsWithinResolution(const TVector3& position,
                                           const TVector3& resolution);

  void SetAx(Double_t Ax) { m_Ax = Ax; }
  void SetAy(Double_t Ay) { m_Ay = Ay; }
  void SetAu(Double_t Au) { m_Au = Au; }
  void SetAv(Double_t Av) { m_Av = Av; }

  void SetParamUsingHoughParam(void);
  void SetHoughFlag(int hough_flag);
  void SetFitFlag(int flag) { m_fitflag = flag; }
  void SetSearchTime(int time) { m_searchtime = time; }
  void SetFitTime(int time) { m_fittime = time; }

  int GetFitFlag(void) const {return m_fitflag;}
  int GetSearchTime() const { return m_searchtime; }
  int GetFitTime() const { return m_fittime; }

  Double_t GetX0() const { return m_x0; }
  Double_t GetY0() const { return m_y0; }
  Double_t GetU0() const { return m_u0; }
  Double_t GetV0() const { return m_v0; }

  Double_t GetAx() const { return m_Ax; }
  Double_t GetAy() const { return m_Ay; }
  Double_t GetAu() const { return m_Au; }
  Double_t GetAv() const { return m_Av; }

  Double_t GetChiSquare() const { return m_chisqr; }
  Double_t GetX(Double_t z) const { return m_x0+m_u0*z; }
  Double_t GetY(Double_t z) const { return m_y0+m_v0*z; }
  Double_t GetS(Double_t z, Double_t tilt) const { return GetX(z)*std::cos(tilt)+GetY(z)*std::sin(tilt); }
  Int_t    GetNIteration() const { return m_n_iteration; }
  Double_t GetTheta() const;
  Double_t GetTrackdE();
  void Print(const TString& arg="", bool print_allhits = false) const;

  //exclusive
  void  CalculateExclusive();
  void  DoLinearFitExclusive();

};

//_____________________________________________________________________________
inline const TString&
TPCLocalTrack::ClassName()
{
  static TString s_name("TPCLocalTrack");
  return s_name;
}

/*
//_____________________________________________________________________________
struct TPCLTrackComp_Nhit
  : public std::binary_function <TPCLocalTrack *, TPCLocalTrack *, Bool_t>
{
  Bool_t operator()(const TPCLocalTrack * const p1,
                    const TPCLocalTrack * const p2) const
    {
      Int_t n1=p1->GetNHit(), n2=p2->GetNHit();

      if(n1>=n2)
        return true;
      else
        return false;
    }
};

//_____________________________________________________________________________
struct TPCLTrackComp_Chisqr
  : public std::binary_function <TPCLocalTrack *, TPCLocalTrack *, Bool_t>
{
  Bool_t operator()(const TPCLocalTrack * const p1,
                    const TPCLocalTrack * const p2) const
    {
      Double_t chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();

      if (chi1 <= chi2)
        return true;
      else
        return false;
    }
};

//_____________________________________________________________________________
struct TPCLTrackComp // TODO
  : public std::binary_function <TPCLocalTrack *, TPCLocalTrack *, Bool_t>
{
  Bool_t operator()(const TPCLocalTrack * const p1,
                    const TPCLocalTrack * const p2) const
    {
      Int_t n1=p1->GetNHit(), n2=p2->GetNHit();
      Double_t chi1=p1->GetChiSquare(), chi2=p2->GetChiSquare();
      if(n1 > n2) return true;
      if(n2 > n1) return false;
      return (chi1 <= chi2);
    }
};
*/
#endif
