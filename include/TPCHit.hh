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

class TPCLTrackHit;
class TPCRawHit;

//_____________________________________________________________________________
class TPCHit : public DCHit
{
public:
  static TString ClassName( void );
  TPCHit( TPCRawHit* rhit );
  TPCHit( Int_t padid, TVector3 pos, Double_t charge); // for cluster hit
  ~TPCHit( void );

private:
  TPCHit( void );
  TPCHit( const TPCHit& );
  TPCHit& operator =( const TPCHit& );

protected:
  TPCRawHit*            m_rhit;
  Int_t                 m_row;
  Int_t                 m_pad;
  Double_t              m_pedestal;
  Double_t              m_rms;
  std::vector<Double_t> m_de;
  std::vector<Double_t> m_time;
  std::vector<Double_t> m_drift_length; // this means Y (beam height = 0)
  std::vector<Double_t> m_chisqr;
  Double_t              m_charge;
  TVector3              m_pos;
  Int_t                 m_is_good;
  Int_t                 m_is_calculated;

  ///// for TPC(MWPC)
  Double_t m_mrow;
  Bool_t   m_tpc_flag;
  Int_t m_clsize;

  Double_t m_resx;
  Double_t m_resy;
  Double_t m_resz;
  Double_t m_res;

  Bool_t m_belong_track;

  DCHit* m_hit_xz;
  DCHit* m_hit_yz;

  std::vector<TPCLTrackHit*> m_register_container;

public:
  Bool_t          BelongToTrack( void ) const { return m_belong_track; }
  Bool_t          CalcTPCObservables( void );
  Double_t        GetChisqr( Int_t i ) const { return m_chisqr.at(i); }
  Int_t           GetChisqrSize( void ) const { return m_chisqr.size(); }
  Double_t        GetDe( Int_t i ) const { return m_de.at(i); }
  Int_t           GetDeSize( void ) const { return m_de.size(); }
  Double_t        GetDriftLength( Int_t i ) const
    { return m_drift_length.at(i); }
  Int_t           GetDriftLengthSize( void ) const
    { return m_drift_length.size(); }
  Int_t           GetNHits( void ) const { return m_de.size(); }
  Int_t           GetPad( void ) const { return m_pad; }
  TPCRawHit*      GetRawHit( void ) const { return m_rhit; }
  Int_t           GetRow( void ) const { return m_row; }
  Double_t        GetPedestal( void ) const { return m_pedestal; }
  Double_t        GetRMS( void ) const { return m_rms; }
  Double_t        GetCharge( void ) const { return m_charge; }
  const TVector3& GetPos( void ) const { return m_pos; }
  Double_t        GetX( void ) const { return m_pos.X(); }
  Double_t        GetY( void ) const { return m_pos.Y(); }
  Double_t        GetZ( void ) const { return m_pos.Z(); }
  DCHit*          GetHitXZ( void ) const { return m_hit_xz; }
  DCHit*          GetHitYZ( void ) const { return m_hit_yz; }
  Double_t        GetMRow( void ) const { return m_mrow; }
  Bool_t          GetTPCFlag( void ) const { return m_tpc_flag; }
  Int_t           GetClusterSize( void ) const { return m_clsize; }
  Double_t        GetResolutionX( void );
  Double_t        GetResolutionY( void );
  Double_t        GetResolutionZ( void );
  Double_t        GetResolution( void );
  Double_t        GetTime( Int_t i ) const { return m_time.at(i); }
  Int_t           GetTimeSize( void ) const { return m_time.size(); }
  Bool_t          IsGood( void ) const { return m_is_good; }
  // Bool_t          IsWithinRange( void ) const
  //   { return m_pair_cont.at(nh).dl_range; }
  void            JoinTrack( void ) { m_belong_track = true; }
  void            Print( const std::string& arg="",
                         std::ostream& ost=hddaq::cout ) const;
  void            QuitTrack( void ) { m_belong_track = false;}
  Bool_t          ReCalculate( Bool_t recursive=false );
  void            RegisterHits( TPCLTrackHit *hit )
    { m_register_container.push_back( hit ); }
  void            SetPad( Int_t pad ) { m_pad = pad; }
  void            SetRow( Int_t row ) { m_row  = row; }
  void            SetCharge( Double_t charge ) { m_charge = charge; }
  void            SetPos( TVector3 pos ) { m_pos = pos; }
  void            SetWirePosition( Double_t wpos ) { m_wpos = wpos; }
  void            SetMRow( Double_t mrow ) { m_mrow = mrow; }
  void            SetClusterSize( Double_t clsize ) { m_clsize = clsize; }
  void            SetTPCFlag( Bool_t flag ) { m_tpc_flag = flag; }
  void            SetResX( Double_t resx ) { m_resx = resx; }
  void            SetResY( Double_t resy ) { m_resy = resy; }
  void            SetResZ( Double_t resz ) { m_resz = resz; }
  void            SetRes( Double_t res ) { m_res = res; }

protected:
  void ClearRegisteredHits( void );
};

//_____________________________________________________________________________
inline TString
TPCHit::ClassName( void )
{
  static TString s_name( "TPCHit" );
  return s_name;
}

//_____________________________________________________________________________
inline std::ostream&
operator <<( std::ostream& ost, const TPCHit& hit )
{
  hit.Print( "", ost );
  return ost;
}

#endif
