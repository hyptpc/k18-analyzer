// -*- C++ -*-

#ifndef KURAMA_TRACK_HH
#define KURAMA_TRACK_HH

#include <vector>
#include <iosfwd>
#include <iostream>
#include <functional>

#include <TString.h>

#include <std_ostream.hh>

#include "RungeKuttaUtilities.hh"
#include "ThreeVector.hh"

class DCLocalTrack;
class TrackHit;
class DCAnalyzer;
class Hodo2Hit;

//_____________________________________________________________________________
class KuramaTrack
{
public:
  static TString ClassName();
  KuramaTrack(DCLocalTrack *track_in, DCLocalTrack *track_out);
  ~KuramaTrack();

private:
  KuramaTrack(const KuramaTrack&);
  KuramaTrack& operator =(const KuramaTrack&);

public:
  enum RKstatus { kInit,
		  kPassed,
		  kExceedMaxPathLength,
		  kExceedMaxStep,
		  kFailedTraceLast,
		  kFailedGuess,
		  kFailedSave,
		  kFatal,
		  nRKstatus };

private:
  TString                 s_status[nRKstatus];
  RKstatus                m_status;
  DCLocalTrack*           m_track_in;
  DCLocalTrack*           m_track_out;
  Double_t                m_tof_seg;
  Double_t                m_initial_momentum;
  std::vector<TrackHit*>  m_hit_array;
  RKHitPointContainer     m_HitPointCont;
  Int_t                   m_n_iteration;
  Int_t                   m_nef_iteration;
  Double_t                m_chisqr;
  Double_t                m_polarity;
  ThreeVector             m_primary_position;
  ThreeVector             m_primary_momentum;
  Double_t                m_path_length_tof;
  Double_t                m_path_length_total;
  ThreeVector             m_tof_pos;
  ThreeVector             m_tof_mom;
  RKCordParameter         m_cord_param;
  Bool_t                  m_is_good;

public:
  DCLocalTrack*      GetLocalTrackIn() { return m_track_in;}
  DCLocalTrack*      GetLocalTrackOut() { return m_track_out; }
  Bool_t             DoFit();
  Bool_t             DoFit(RKCordParameter iniCord);
  Bool_t             DoFitMinuit();
  Bool_t             Status() const { return m_status; }
  Int_t              Niteration() const { return m_n_iteration; }
  void               SetInitialMomentum(Double_t p) { m_initial_momentum = p; }
  const ThreeVector& PrimaryPosition() const { return m_primary_position; }
  const ThreeVector& PrimaryMomentum() const { return m_primary_momentum; }
  Double_t           PrimaryMomMag() const { return m_primary_momentum.Mag();}
  Double_t           PathLengthToTOF() const { return m_path_length_tof; }
  Double_t           PathLengthTotal() const { return m_path_length_total; }
  Double_t           TofSeg() const { return m_tof_seg; }
  const ThreeVector& TofPos() const { return m_tof_pos; }
  const ThreeVector& TofMom() const { return m_tof_mom; }
  Double_t           chisqr() const { return m_chisqr; }
  Double_t           Polarity() const { return m_polarity; }
  std::size_t        GetNHits() const { return m_hit_array.size(); }
  TrackHit*          GetHit(std::size_t nth) const;
  TrackHit*          GetHitOfLayerNumber(Int_t lnum) const;
  Bool_t             GoodForAnalysis() const { return m_is_good; }
  Bool_t             GoodForAnalysis(Bool_t status)
    { m_is_good = status; return status; }
  Double_t           GetInitialMomentum() const { return m_initial_momentum; }
  Bool_t             GetTrajectoryLocalPosition(Int_t layer, Double_t& x, Double_t& y) const;
  void               Print(const TString& arg="", std::ostream& ost=hddaq::cout);
  Bool_t             ReCalc(Bool_t applyRecursively=false);

private:
  void     ClearHitArray();
  Double_t CalcChiSqr(const RKHitPointContainer &hpCont) const;
  void     FillHitArray();
  Bool_t   GuessNextParameters(const RKHitPointContainer &hpCont,
                               RKCordParameter &Cord,
                               Double_t &estDeltaChisqr,
                               Double_t &lambdaCri, Double_t dmp=0.) const;
  void     PrintCalcHits(const RKHitPointContainer &hpCont,
                         std::ostream &ost = std::cout) const;
  Bool_t   SaveCalcPosition(const RKHitPointContainer &hpCont);
  Bool_t   SaveTrackParameters(const RKCordParameter &cp);

};

//_____________________________________________________________________________
inline TString
KuramaTrack::ClassName()
{
  static TString s_name("KuramaTrack");
  return s_name;
}

//_____________________________________________________________________________
struct KuramaTrackComp
  : public std::binary_function<KuramaTrack*, KuramaTrack*, Bool_t>
{
  Bool_t operator()(const KuramaTrack* const p1,
                    const KuramaTrack* const p2) const
  {
    Int_t n1=p1->GetNHits(), n2=p2->GetNHits();
    if(n1>n2+1) return true;
    if(n2>n1+1) return false;
    return (p1->chisqr())<(p2->chisqr());
  }
};

#endif
