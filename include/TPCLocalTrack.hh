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
  bool   m_is_fitted;     // flag of DoFit()
  bool   m_is_calculated; // flag of Calculate()
  std::vector<TPCLTrackHit*> m_hit_array;
  std::vector<TPCCluster*> m_cluster_array;
  // TH2D *hist;

  double m_Ax;
  double m_Ay;
  double m_Au;
  double m_Av;
  // double m_Chix;
  // double m_Chiy;
  // double m_Chiu;
  // double m_Chiv;

  double m_Az;
  double m_Bz;

  double m_x0;
  double m_y0;
  double m_u0;
  double m_v0;
  double m_a;
  double m_b;
  double m_chisqr;
  bool   m_good_for_tracking;
  double m_n_iteration;
  // for SSD
  double m_de;
  TMinuit *m_minuit;

public:
  void         AddTPCHit(TPCLTrackHit *hit);
  void         AddTPCCluster(TPCCluster *cluster);//not supported
  void         ClearHits();
  void         Calculate();
  void         DeleteNullHit();
  bool         DoLinearFit(Int_t MinHits);
  //  bool         DoHelixFit();
  bool         DoFit(Int_t MinHits);
  Int_t          GetNDF() const;
  Int_t          GetNHit() const { return m_hit_array.size();  }
  TPCLTrackHit* GetHit(std::size_t nth) const;
  bool         IsFitted() const { return m_is_fitted; }
  bool         IsCalculated() const { return m_is_calculated; }
  bool         Residual_check(TVector3 pos, TVector3 Res);
  void         CalcChi2(void);


  void SetAx(double Ax) { m_Ax = Ax; }
  void SetAy(double Ay) { m_Ay = Ay; }
  void SetAu(double Au) { m_Au = Au; }
  void SetAv(double Av) { m_Av = Av; }
  void SetAz(double Az){  m_Az = Az; }
  void SetBz(double Bz){  m_Bz = Bz; }
  // void SetChix(double Chix) { m_Chix = Chix; }
  // void SetChiy(double Chiy) { m_Chiy = Chiy; }
  // void SetChiu(double Chiu) { m_Chiu = Chiu; }
  // void SetChiv(double Chiv) { m_Chiv = Chiv; }
  void SetDe(double de) { m_de = de; }

  double GetX0() const { return m_x0; }
  double GetY0() const { return m_y0; }
  double GetU0() const { return m_u0; }
  double GetV0() const { return m_v0; }

  double GetAx() const { return m_Ax; }
  double GetAy() const { return m_Ay; }
  double GetAu() const { return m_Au; }
  double GetAv() const { return m_Av; }

  double GetAz () const { return m_Az;  }
  double GetBz () const { return m_Bz;  }


  double GetChiSquare() const { return m_chisqr; }
  // double GetChiX() const { return m_Chix; }
  // double GetChiY() const { return m_Chiy; }
  // double GetChiU() const { return m_Chiu; }
  // double GetChiV() const { return m_Chiv; }
  double GetX(double z) const { return m_x0+m_u0*z; }
  double GetY(double z) const { return m_y0+m_v0*z; }
  double GetS(double z, double tilt) const { return GetX(z)*std::cos(tilt)+GetY(z)*std::sin(tilt); }
  Int_t    GetNIteration() const { return m_n_iteration; }
  double GetTheta() const;
  bool   GoodForTracking() const { return m_good_for_tracking; }
  bool   GoodForTracking(bool status)
  { bool ret = m_good_for_tracking; m_good_for_tracking = status; return ret; }
  double GetDe() const { return m_de; }
  void   Print(const std::string& arg="", std::ostream& ost=hddaq::cout) const;

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
  : public std::binary_function <TPCLocalTrack *, TPCLocalTrack *, bool>
{
  bool operator()(const TPCLocalTrack * const p1,
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
  : public std::binary_function <TPCLocalTrack *, TPCLocalTrack *, bool>
{
  bool operator()(const TPCLocalTrack * const p1,
                   const TPCLocalTrack * const p2) const
  {
    double chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();

    if (chi1 <= chi2)
      return true;
    else
      return false;
  }
};

//_____________________________________________________________________________
struct TPCLTrackComp // TODO
  : public std::binary_function <TPCLocalTrack *, TPCLocalTrack *, bool>
{
  bool operator()(const TPCLocalTrack * const p1,
		   const TPCLocalTrack * const p2) const
  {
    Int_t n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(), chi2=p2->GetChiSquare();
    if(n1 > n2) return true;
    if(n2 > n1) return false;
    return (chi1 <= chi2);

  }
};

#endif
