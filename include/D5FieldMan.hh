// -*- C++ -*-

#ifndef D5_FIELD_MAN_HH
#define D5_FIELD_MAN_HH

#include <map>

#include <TString.h>

//_____________________________________________________________________________
// Run-by-run D5 field (param/D5): mean D5 field [T] and its rms during each run.
// The K18 central momentum follows from the D5 field:
//   p0 [GeV/c] = 0.299792458 * |B| [T] * R [m] * FieldFactor
class D5FieldMan
{
public:
  static const TString& ClassName( void );
  static D5FieldMan&    GetInstance( void );
  ~D5FieldMan( void );

  static constexpr Double_t kRadius      = 1.8984; // D5 bending radius [m]
  static constexpr Double_t kFieldFactor = 1.1;
  // A run whose D5 rms exceeds this changed D5 during the run [T].
  static constexpr Double_t kMaxRms      = 1.e-3;

private:
  D5FieldMan( void );
  D5FieldMan( const D5FieldMan& );
  D5FieldMan& operator =( const D5FieldMan& );

private:
  struct Entry { Double_t field; Double_t rms; };
  std::map<UInt_t, Entry> m_table;
  Bool_t                  m_is_ready;
  TString                 m_file_name;

public:
  Bool_t   Initialize( void );
  Bool_t   Initialize( const TString& file_name );
  Bool_t   IsReady( void ) const { return m_is_ready; }
  // false if the run is not in the table
  Bool_t   GetField( UInt_t run, Double_t& field, Double_t& rms ) const;
  static Double_t CentralMomentum( Double_t field );
};

//_____________________________________________________________________________
inline const TString&
D5FieldMan::ClassName( void )
{
  static TString s_name( "D5FieldMan" );
  return s_name;
}

//_____________________________________________________________________________
inline D5FieldMan&
D5FieldMan::GetInstance( void )
{
  static D5FieldMan s_instance;
  return s_instance;
}

#endif
