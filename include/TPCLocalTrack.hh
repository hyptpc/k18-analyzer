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
class DCAnalyzer;

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
  std::vector<TPCLTrackHit*> m_hit_array;
  // std::vector<TPCCluster*> m_cluster_array;

  // track coordinate origin is target, ***NOT TPC center***
  Double_t m_Ax;
  Double_t m_Ay;
  Double_t m_Au;
  Double_t m_Av;

  Double_t m_Az;
  Double_t m_Bz;

  Double_t m_x0;
  Double_t m_y0;
  Double_t m_u0;
  Double_t m_v0;
  Double_t m_a;
  Double_t m_b;
  Double_t m_chisqr;
  Bool_t   m_good_for_tracking;
  Double_t m_n_iteration;

  Double_t m_de;
  TMinuit *m_minuit;

public:
  static const Int_t NumOfParam = 4;

  void          AddTPCHit(TPCLTrackHit *hit);
  // void          AddTPCCluster(TPCCluster *cluster);//not supported
  void          ClearHits();
  void          Calculate();
  void          CalcChisquare();
  void          DeleteNullHit();
  Bool_t        DoFit();
  Bool_t        DoFit(Int_t min_hits); // obsolete
  Bool_t        DoFitLinear(Int_t min_hits); // obsolete
  Int_t         GetNDF() const { return m_hit_array.size() - NumOfParam;  }
  Int_t         GetNHit() const { return m_hit_array.size();  }
  TPCLTrackHit* GetHit(Int_t i) const { return m_hit_array.at(i); }
  Bool_t        IsFitted() const { return m_is_fitted; }
  Bool_t        IsCalculated() const { return m_is_calculated; }
  Bool_t        Residual_check(TVector3 pos, TVector3 Res);

  void SetAx(Double_t Ax) { m_Ax = Ax; }
  void SetAy(Double_t Ay) { m_Ay = Ay; }
  void SetAu(Double_t Au) { m_Au = Au; }
  void SetAv(Double_t Av) { m_Av = Av; }
  void SetAz(Double_t Az){  m_Az = Az; }
  void SetBz(Double_t Bz){  m_Bz = Bz; }
  void SetDe(Double_t de) { m_de = de; }

  Double_t GetX0() const { return m_x0; }
  Double_t GetY0() const { return m_y0; }
  Double_t GetU0() const { return m_u0; }
  Double_t GetV0() const { return m_v0; }

  Double_t GetAx() const { return m_Ax; }
  Double_t GetAy() const { return m_Ay; }
  Double_t GetAu() const { return m_Au; }
  Double_t GetAv() const { return m_Av; }

  Double_t GetAz () const { return m_Az;  }
  Double_t GetBz () const { return m_Bz;  }


  Double_t GetChiSquare() const { return m_chisqr; }
  // Double_t GetChiX() const { return m_Chix; }
  // Double_t GetChiY() const { return m_Chiy; }
  // Double_t GetChiU() const { return m_Chiu; }
  // Double_t GetChiV() const { return m_Chiv; }
  Double_t GetX(Double_t z) const { return m_x0+m_u0*z; }
  Double_t GetY(Double_t z) const { return m_y0+m_v0*z; }
  Double_t GetS(Double_t z, Double_t tilt) const { return GetX(z)*std::cos(tilt)+GetY(z)*std::sin(tilt); }
  Int_t    GetNIteration() const { return m_n_iteration; }
  Double_t GetTheta() const;
  Bool_t   GoodForTracking() const { return m_good_for_tracking; }
  Bool_t   GoodForTracking(Bool_t status)
    { Bool_t ret = m_good_for_tracking; m_good_for_tracking = status; return ret; }
  Double_t GetDe() const { return m_de; }
  void   Print(const TString& arg="", std::ostream& ost=hddaq::cout) const;

};

//_____________________________________________________________________________
inline const TString&
TPCLocalTrack::ClassName()
{
  static TString s_name("TPCLocalTrack");
  return s_name;
}

//_____________________________________________________________________________
inline std::ostream&
operator <<(std::ostream& ost,
            const TPCLocalTrack& track)
{
  track.Print("", ost);
  return ost;
}

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

#endif
