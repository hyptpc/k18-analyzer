// -*- C++ -*-

#ifndef TPC_RAW_HIT_HH
#define TPC_RAW_HIT_HH

#include <vector>

#include <TMath.h>
#include <TString.h>

typedef std::vector<Int_t> FADC_t;

//_____________________________________________________________________________
class TPCRawHit
{
public:
  static TString ClassName( void );
  TPCRawHit( Int_t layer, Int_t row );
  ~TPCRawHit( void );

private:
  Int_t  m_pad_id;
  Int_t  m_layer_id;
  Int_t  m_row_id;
  FADC_t m_fadc;

public:
  void          AddFadc( Int_t adc );
  const FADC_t& Fadc( void ) const { return m_fadc; }
  Int_t         LayerId( void ) const { return m_layer_id; }
  Int_t         PadId( void ) const { return m_pad_id; }
  void          Print( Option_t* opt="" ) const;
  Int_t         RowId( void ) const { return m_row_id; }
  Double_t      LocMax( void ) const
    { return TMath::LocMax( m_fadc.size(), m_fadc.data() ); }
  Double_t      MaxAdc( void ) const
    { return TMath::MaxElement( m_fadc.size(), m_fadc.data() ); }
  Double_t      Mean( void ) const
    { return TMath::Mean( m_fadc.size(), m_fadc.data() ); }
  Double_t      RMS( void ) const
    { return TMath::RMS( m_fadc.size(), m_fadc.data() ); }
};

//_____________________________________________________________________________
inline TString
TPCRawHit::ClassName( void )
{
  static TString s_name( "TPCRawHit" );
  return s_name;
}

#endif
