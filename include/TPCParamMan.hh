// -*- C++ -*-

#ifndef TPC_PARAM_MAN_HH
#define TPC_PARAM_MAN_HH

#include <map>
#include <TString.h>

//_____________________________________________________________________________
class TPCAParam
{
public:
  TPCAParam( Double_t offset, Double_t gain )
    : m_offset( offset ),
      m_gain( gain )
    {}
  ~TPCAParam( void )
    {}

private:
  TPCAParam( void );
  TPCAParam( const TPCAParam& );
  TPCAParam& operator =( const TPCAParam& );

private:
  Double_t m_offset;
  Double_t m_gain;

public:
  Double_t Offset( void ) const { return m_offset; }
  Double_t Gain( void ) const { return m_gain; }
  Double_t CDeltaE( Double_t de ) const { return ( de - m_offset ) * m_gain; }
  Double_t Offset( Double_t cde ) const { return cde / m_gain + m_offset; }
  Double_t P0( void ) const { return m_offset; }
  Double_t P1( void ) const { return m_gain; }
};

//_____________________________________________________________________________
class TPCTParam
{
public:
  TPCTParam( Double_t offset, Double_t gain )
    : m_offset( offset ),
      m_gain( gain )
    {}
  ~TPCTParam( void )
    {}

private:
  TPCTParam( void );
  TPCTParam( const TPCTParam& );
  TPCTParam& operator =( const TPCTParam& );

private:
  Double_t m_offset;
  Double_t m_gain;

public:
  Double_t Offset( void ) const { return m_offset; }
  Double_t Gain( void ) const { return m_gain; }
  Double_t P0( void ) const { return m_offset; }
  Double_t P1( void ) const { return m_gain; }
  Double_t Time( Double_t tdc ) const { return ( tdc - m_offset ) * m_gain; }
  Double_t Tdc( Double_t time ) const { return time / m_gain + m_offset; }
};

//_____________________________________________________________________________
class TPCYParam
{
public:
  TPCYParam( Double_t offset, Double_t drift_velocity )
    : m_offset( offset ),
      m_drift_velocity( drift_velocity )
    {}
  ~TPCYParam( void )
    {}

private:
  TPCYParam( void );
  TPCYParam( const TPCYParam& );
  TPCYParam& operator =( const TPCYParam& );

private:
  Double_t m_offset;
  Double_t m_drift_velocity;

public:
  Double_t Offset( void ) const { return m_offset; }
  Double_t DriftLength( Double_t time ) const { return Y( time ); }
  Double_t DriftVelocity( void ) const { return m_drift_velocity; }
  Double_t P0( void ) const { return m_offset; }
  Double_t P1( void ) const { return m_drift_velocity; }
  Double_t Time( Double_t dl ) const
    { return dl / m_drift_velocity + m_offset; }
  Double_t Y( Double_t time ) const
    { return ( time - m_offset ) * m_drift_velocity; }
};

//_____________________________________________________________________________
class TPCParamMan
{
public:
  static const TString& ClassName( void );
  static TPCParamMan&   GetInstance( void );
  ~TPCParamMan( void );

private:
  TPCParamMan( void );
  TPCParamMan( const TPCParamMan& );
  TPCParamMan& operator =( const TPCParamMan& );

private:
  enum eAorT { kAdc, kTdc, kY, kATY };
  typedef std::map<Int_t, TPCAParam*> AContainer;
  typedef std::map<Int_t, TPCTParam*> TContainer;
  typedef std::map<Int_t, TPCYParam*> YContainer;
  typedef AContainer::const_iterator AIterator;
  typedef TContainer::const_iterator TIterator;
  typedef YContainer::const_iterator YIterator;
  Bool_t     m_is_ready;
  TString    m_file_name;
  AContainer m_APContainer;
  TContainer m_TPContainer;
  YContainer m_YPContainer;

public:
  Bool_t GetCDe( Int_t layer, Int_t row, Double_t de,
                 Double_t &cde ) const;
  Bool_t GetCTime( Int_t layer, Int_t row, Double_t time,
                   Double_t &ctime ) const;
  Double_t GetC_Clock( Int_t layer, Int_t row, Double_t time ) const;
  Bool_t GetDriftLength( Int_t layer, Int_t row, Double_t time,
                         Double_t& y ) const
    { return GetY( layer, row, time, y ); }
  Bool_t GetY( Int_t layer, Int_t row, Double_t time, Double_t& y ) const;
  Bool_t Initialize( void );
  Bool_t Initialize( const TString& file_name );
  Bool_t IsReady( void ) const { return m_is_ready; }
  void   SetFileName( const TString& file_name ) { m_file_name = file_name; }

private:
  void       ClearACont( void );
  void       ClearTCont( void );
  void       ClearYCont( void );
  TPCAParam* GetAmap( Int_t layer, Int_t row ) const;
  TPCTParam* GetTmap( Int_t layer, Int_t row ) const;
  TPCYParam* GetYmap( Int_t layer, Int_t row ) const;
};

//_____________________________________________________________________________
inline const TString&
TPCParamMan::ClassName( void )
{
  static TString g_name( "TPCParamMan" );
  return g_name;
}

//_____________________________________________________________________________
inline TPCParamMan&
TPCParamMan::GetInstance( void )
{
  static TPCParamMan s_instance;
  return s_instance;
}

#endif
