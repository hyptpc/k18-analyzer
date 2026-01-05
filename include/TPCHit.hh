// -*- C++ -*-

#ifndef TPC_HIT_HH
#define TPC_HIT_HH

#include "DCHit.hh"

#include <cmath>
#include <vector>
#include <deque>
#include <string>
#include <numeric>

#include <std_ostream.hh>

#include <TVector3.h>
#include <ThreeVector.hh>


class TPCRawHit;

//_____________________________________________________________________________
class TPCHit : public DCHit
{
public:
  static TString ClassName();
  TPCHit(TPCRawHit* rhit);
  ~TPCHit();

private:
  TPCHit();
  TPCHit(const TPCHit&);
  TPCHit& operator =(const TPCHit&);

protected:
  TPCRawHit*            m_rhit;
  Int_t                 m_layer;
  Int_t                 m_row;
  Double_t              m_padtheta;
  Double_t              m_padlength;
  Double_t              m_mrow;
  Int_t                 m_pad;
  Double_t              m_pedestal;
  Double_t              m_rms;
  Double_t 		m_raw_rms;
  std::vector<Double_t> m_de;
  std::vector<Double_t> m_sigma;
  std::vector<Double_t> m_time;
  std::vector<Double_t> m_chisqr;
  std::vector<Double_t> m_cde;
  std::vector<Double_t> m_ctime;
  std::vector<Double_t> m_drift_length; // this means Y (beam height = 0)
  std::vector<TVector3> m_position;
  Bool_t                m_is_good;
  Int_t                 m_is_calculated;

public:
  void            AddHit(Double_t de, Double_t time, Double_t sigma=0.,
                         Double_t chisqr=0.);
  Bool_t          Calculate(Double_t clock=0.);
  Bool_t          DoFit();
  Double_t        GetCDe(Int_t i=0) const { return m_cde.at(i); }
  Int_t           GetCDeSize() const { return m_cde.size(); }
  Double_t        GetChisqr(Int_t i=0) const { return m_chisqr.at(i); }
  Int_t           GetChisqrSize() const { return m_chisqr.size(); }
  Double_t        GetCTime(Int_t i=0) const { return m_ctime.at(i); }
  Int_t           GetCTimeSize() const { return m_ctime.size(); }
  Double_t        GetDe(Int_t i=0) const { return m_de.at(i); }
  Double_t        GetSigma(Int_t i=0) const { return m_sigma.at(i); }
  Int_t           GetDeSize() const { return m_de.size(); }
  Double_t        GetDriftLength(Int_t i=0) const
  { return m_drift_length.at(i); }
  Int_t           GetDriftLengthSize() const
  { return m_drift_length.size(); }
  Int_t           GetNHits() const { return m_de.size(); }
  Int_t           GetPad() const { return m_pad; }
  TPCRawHit*      GetRawHit() const { return m_rhit; }
  Int_t           GetLayer() const { return m_layer; }
  Int_t           GetRow() const { return m_row; }
  Double_t        GetPedestal() const { return m_pedestal; }
  Double_t        GetRMS() const { return m_rms; }
  Double_t	  GetRawRMS()const{return m_raw_rms;}
  Double_t        GetX(Int_t i=0) const { return m_position.at(i).X(); }
  Double_t        GetY(Int_t i=0) const { return m_position.at(i).Y(); }
  Double_t        GetZ(Int_t i=0) const { return m_position.at(i).Z(); }
  Double_t        GetPadTheta() { return m_padtheta; }
  Double_t        GetPadLength() const { return m_padlength; }
  Double_t        GetMRow() const { return m_mrow; }

  Double_t        GetTime(Int_t i=0) const { return m_time.at(i); }
  Int_t           GetTimeSize() const { return m_time.size(); }
  Bool_t          IsGood() const;
  void            Print(const std::string& arg="", std::ostream& ost=hddaq::cout) const;

protected:
  void ClearRegisteredHits();
};

//_____________________________________________________________________________
inline TString
TPCHit::ClassName()
{
  static TString s_name("TPCHit");
  return s_name;
}

//_____________________________________________________________________________
inline std::ostream&
operator <<(std::ostream& ost, const TPCHit& hit)
{
  hit.Print("", ost);
  return ost;
}

#endif
