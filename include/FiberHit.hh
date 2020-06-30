// -*- C++ -*-

#ifndef FIBER_HIT_HH
#define FIBER_HIT_HH

#include <string>

#include <std_ostream.hh>

#include "Hodo1Hit.hh"

class FLHit;

//______________________________________________________________________________
class FiberHit : public Hodo1Hit
{
public:
  explicit  FiberHit( HodoRawHit *object, const char* name );
  virtual  ~FiberHit( void );

private:
  FiberHit( void );
  FiberHit( const FiberHit& object );
  FiberHit& operator =( const FiberHit& object );

protected:
  std::string m_detector_name;
  int         m_segment;
  int         m_ud;
  double      m_position;
  double      m_offset;
  int         m_pair_id;
  bool        m_status;
  std::vector<FLHit*> m_hit_container;
  double      m_adc_hg;
  double      m_adc_lg;
  double      m_pedcor_hg;
  double      m_pedcor_lg;
  double      m_mip_hg;
  double      m_mip_lg;
  double      m_dE_hg;
  double      m_dE_lg;

  struct data_pair{
    double time_l;
    double time_t;
    double ctime_l;
    double tot;
    int    index_t;
  };

  std::vector<data_pair> m_pair_cont;

public:
  void   SetDetectorName( const char* name ){ m_detector_name = name; }
  void   SetDetectorName( const std::string& name ){ m_detector_name = name; }
  bool   Calculate( void );
  // Call super class method
  int    GetNLeading( void )  const { return Hodo1Hit::GetNumOfHit(0);    }
  int    GetNTrailing( void ) const { return Hodo1Hit::GetNumOfHit(1);    }
  double GetLeading( int n=0 )  const { return m_ud==0? m_raw->GetTdc1(n)  : m_raw->GetTdc2(n);}
  double GetTrailing( int n=0 ) const { return m_ud==0? m_raw->GetTdcT1(n) : m_raw->GetTdcT2(n);}

  // Call member in this class
  int    GetNPair( void )       const { return m_pair_cont.size();        }
  double GetTime( int n=0 )     const { return m_pair_cont.at(n).time_l;  } // Leading
  double GetCTime( int n=0 )    const { return m_pair_cont.at(n).ctime_l; } // Leading
  double GetTimeT( int n=0 )    const { return m_pair_cont.at(n).time_t;  } // Trailing
  double GetWidth( int n=0 )    const { return m_pair_cont.at(n).tot;     }
  double GetTot( int n=0 )      const { return m_pair_cont.at(n).tot;     }

  double GetPosition( void )  const { return m_position + m_offset;       }

  int    PairId( void )       const { return m_pair_id;                   }
  //  virtual double SegmentId( void )    const { return m_segment;                   }

  double GetAdcHG( void ) const { return m_adc_hg; }
  double GetAdcLG( void ) const { return m_adc_lg; }
  double GetMipHG( void ) const { return m_mip_hg; }
  double GetMipLG( void ) const { return m_mip_lg; }
  double GetDeHG( void ) const { return m_dE_hg; }
  double GetDeLG( void ) const { return m_dE_lg; }
  void SetPedestalCor( double deltaHG, double deltaLG ){ m_pedcor_hg = deltaHG; m_pedcor_lg = deltaLG; }

  void   Print( const std::string& arg="", std::ostream& ost=hddaq::cout ) const;
  void   RegisterHits( FLHit* hit )   { m_hit_container.push_back(hit);   }

  virtual bool ReCalc( bool allpyRecursively=false )
  { return FiberHit::Calculate(); }

  static bool CompFiberHit( const FiberHit* left, const FiberHit* right );

};

//______________________________________________________________________________
inline bool
FiberHit::CompFiberHit( const FiberHit* left, const FiberHit* right )
{
  return left->PairId() < right->PairId();
}

#endif
