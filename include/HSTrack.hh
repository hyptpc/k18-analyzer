// -*- C++ -*-
#ifndef HSTRACK_HH
#define HSTRACK_HH

#include <vector>
#include <iosfwd>
#include <iostream>
#include <functional>

#include <TString.h>
#include <TVector3.h>

#include <std_ostream.hh>

#include "RungeKuttaUtilities.hh"

class DCLocalTrack;
class TrackHit;
class DCAnalyzer;
class Hodo2Hit;
class K18TrackD2U;

//_____________________________________________________________________________
class HSTrack
{
public:
  static const TString& ClassName();
  HSTrack(Double_t x, Double_t y,
	  Double_t u, Double_t v, Double_t p);
  ~HSTrack();
  
private:
  HSTrack(const HSTrack&);
  HSTrack& operator =(const HSTrack&);

public:
  enum RKstatus { kInit,
		  kPassed,
		  kExceedMaxPathLength,
		  kExceedMaxStep,
		  kFailedSave,
		  kFatal,
		  nRKstatus };

private:

  TString                 s_status[nRKstatus];
  RKstatus                m_status;
  Double_t xInit, yInit;
  Double_t uInit, vInit;
  Double_t                m_initial_momentum;
  RKHitPointContainer     m_HitPointCont;
  Double_t                m_polarity;
  TVector3                m_tgt_position;
  TVector3                m_tgt_momentum;

  TVector3                m_bh2_position;
  TVector3                m_bh2_momentum;

  TVector3                m_vp1_position;
  TVector3                m_vp1_momentum;

  TVector3                m_vp3_position;
  TVector3                m_vp3_momentum;

  TVector3                m_htof_position;
  TVector3                m_htof_momentum;

  Double_t                m_path_length_total;
  RKCordParameter         m_cord_param;
  Bool_t                  m_is_good;
  
public:
  Bool_t          Propagate();
  Bool_t          Propagate(RKCordParameter iniCord);
  Double_t        GetInitialMomentum() const { return m_initial_momentum; }
  Bool_t          GoodForAnalysis() const { return m_is_good; }
  Bool_t          GoodForAnalysis(Bool_t status)
  { m_is_good = status; return status; }
  Double_t        PathLength() const { return m_path_length_total; }
  Double_t        Polarity() const { return m_polarity; }

  const TVector3& TgtMomentum() const { return m_tgt_momentum; }
  Double_t        TgtMomMag() const { return m_tgt_momentum.Mag();}
  const TVector3& TgtPosition() const { return m_tgt_position; }
  
  const TVector3& BH2Momentum() const { return m_bh2_momentum; }
  const TVector3& BH2Position() const { return m_bh2_position; }

  const TVector3& VP1Momentum() const { return m_vp1_momentum; }
  const TVector3& VP1Position() const { return m_vp1_position; }

  const TVector3& VP3Momentum() const { return m_vp3_momentum; }
  const TVector3& VP3Position() const { return m_vp3_position; }

  const TVector3& HtofMomentum() const { return m_htof_momentum; }
  const TVector3& HtofPosition() const { return m_htof_position; }

  void            Print(const TString& arg="", std::ostream& ost=hddaq::cout);
  void            SetInitialMomentum(Double_t p) { m_initial_momentum = p; }
  Bool_t          ExtrapStatus() const { return m_status; }
  
private:
  void     ClearHitArray();
  Bool_t   SaveTrack();

};

//_____________________________________________________________________________
inline const TString&
HSTrack::ClassName()
{
  static TString s_name("HSTrack");
  return s_name;
}

//_____________________________________________________________________________
#endif
