/**
 *  file: DCLocalTrack.hh
 *  date: 2017.04.10
 *
 */

#ifndef DC_LOCAL_TRACK_HH
#define DC_LOCAL_TRACK_HH

#include <vector>
#include <functional>

#include <std_ostream.hh>

#include "ThreeVector.hh"
#include "DCLTrackHit.hh"
#include "DetectorID.hh"

class DCLTrackHit;
class DCAnalyzer;

const double MaxChi2CFT  = 150.;

//______________________________________________________________________________
class DCLocalTrack
{
public:
  explicit DCLocalTrack( void );
  ~DCLocalTrack( void );

private:
  DCLocalTrack( const DCLocalTrack & );
  DCLocalTrack & operator =( const DCLocalTrack & );

private:
  bool   m_is_fitted;     // flag of DoFit()
  bool   m_is_calculated; // flag of Calculate()
  std::vector<DCLTrackHit*> m_hit_array;
  std::vector<DCLTrackHit*> m_hit_arrayUV;
  double m_Ax;
  double m_Ay;
  double m_Au;
  double m_Av;
  double m_Chix;
  double m_Chiy;
  double m_Chiu;
  double m_Chiv;

  double m_x0;
  double m_y0;
  double m_u0;
  double m_v0;
  double m_a;
  double m_b;
  double m_chisqr;
  bool   m_good_for_tracking;
  // for SSD
  double m_de;
  // for Honeycomb
  double m_chisqr1st; // 1st iteration for honeycomb
  double m_n_iteration;
  // for CFT
  int m_xyFitFlag; // 0:y=ax+b, 1:x=ay+b
  int m_zTrackFlag;
  double m_Axy, m_Bxy;
  double m_Az , m_Bz ;
  double m_chisqrXY;
  double m_chisqrZ;
  double m_vtx_z;
  double m_theta;
  ThreeVector m_Dir, m_Pos0;
  double m_meanseg[NumOfPlaneCFT*2];
  // phi before pos. calib.
  double m_phi_ini_before[NumOfPlaneCFT*2];  //initial phi from segNo
  double m_phi_track_before[NumOfPlaneCFT*2];// phi from tracking
  double m_dphi_before[NumOfPlaneCFT*2];     // delta phi (ini - calc)
  double m_dist_phi_before[NumOfPlaneCFT*2]; // residual of XY plane
  // phi after 
  double m_phi_ini[NumOfPlaneCFT*2];  //initial phi from segNo
  double m_phi_track[NumOfPlaneCFT*2];// phi from tracking
  double m_dphi[NumOfPlaneCFT*2];     // delta phi (ini - calc)
  double m_dist_phi[NumOfPlaneCFT*2]; // residual of XY plane
  // z before pos. calib.
  double m_z_ini_before[NumOfPlaneCFT*2];    //initial z from segNo&phi
  double m_z_track_before[NumOfPlaneCFT*2];  // z from tracking
  double m_dz_before[NumOfPlaneCFT*2];     // delta z (ini - calc)
  double m_dist_uv_before[NumOfPlaneCFT*2]; // residual of Z_XY plane
  // z after 
  double m_z_ini[NumOfPlaneCFT*2];    //initial z from segNo&phi
  double m_z_track[NumOfPlaneCFT*2];  // z from tracking
  double m_dz[NumOfPlaneCFT*2];     // delta z (ini - calc)
  double m_dist_uv[NumOfPlaneCFT*2]; // residual of Z_XY plane
  // dr
  double m_dr[NumOfPlaneCFT*2];     // delta r (org - tracked r)
  double m_r[NumOfPlaneCFT*2];     // cariblated r 

  // dE value
  double m_sum_adc[NumOfPlaneCFT*2]; double m_max_adc[NumOfPlaneCFT*2];
  double m_sum_mip[NumOfPlaneCFT*2]; double m_max_mip[NumOfPlaneCFT*2];
  double m_sum_dE[NumOfPlaneCFT*2];  double m_max_dE[NumOfPlaneCFT*2]; 
  double m_total_adc   , m_total_max_adc;
  double m_total_mip   , m_total_max_mip;
  double m_total_dE    , m_total_max_dE;
  double m_total_dE_phi, m_total_max_dE_phi; // only phi sum
  double m_total_dE_uv , m_total_max_dE_uv;  // only uv  sum


public:
  void         AddHit( DCLTrackHit *hitp );
  void         AddHitUV( DCLTrackHit *hitp );
  void         Calculate( void );
  void         CalculateCFT( void );
  void         DeleteNullHit( void );
  bool         DoFit( void );
  bool         DoFitBcSdc( void );
  bool         FindLayer( int layer ) const;
  int          GetNDF( void ) const;
  int          GetNHit( void ) const { return m_hit_array.size();  }
  int          GetNHitUV(void )const { return m_hit_arrayUV.size();}
  int          GetNHitSFT( void ) const;
  int          GetNHitY( void ) const;
  DCLTrackHit* GetHit( std::size_t nth ) const;
  DCLTrackHit* GetHitUV( std::size_t nth ) const;
  DCLTrackHit* GetHitOfLayerNumber( int lnum ) const;
  double       GetWire( int layer ) const;
  bool         HasHoneycomb( void ) const;
  bool         IsFitted( void ) const { return m_is_fitted; }
  bool         IsCalculated( void ) const { return m_is_calculated; }

  void SetAx( double Ax ) { m_Ax = Ax; }
  void SetAy( double Ay ) { m_Ay = Ay; }
  void SetAu( double Au ) { m_Au = Au; }
  void SetAv( double Av ) { m_Av = Av; }
  void SetAz( double Az ){  m_Az = Az; }
  void SetBz( double Bz ){  m_Bz = Bz; }
  void SetChix( double Chix ) { m_Chix = Chix; }
  void SetChiy( double Chiy ) { m_Chiy = Chiy; }
  void SetChiu( double Chiu ) { m_Chiu = Chiu; }
  void SetChiv( double Chiv ) { m_Chiv = Chiv; }
  void SetDe( double de ) { m_de = de; }

  double GetX0( void ) const { return m_x0; }
  double GetY0( void ) const { return m_y0; }
  double GetU0( void ) const { return m_u0; }
  double GetV0( void ) const { return m_v0; }

  //For XUV Tracking
  bool DoFitVXU( void );

  double GetVXU_A( void ) const { return m_a; }
  double GetVXU_B( void ) const { return m_b; }
  double GetVXU( double z ) const { return m_a*z+m_b; }
  double GetAx( void ) const { return m_Ax; }
  double GetAy( void ) const { return m_Ay; }
  double GetAu( void ) const { return m_Au; }
  double GetAv( void ) const { return m_Av; }

  double GetAz ( void ) const { return m_Az;  }
  double GetBz ( void ) const { return m_Bz;  }

  double GetDifVXU( void ) const ;
  double GetDifVXUSDC34( void ) const;
  double GetChiSquare( void ) const { return m_chisqr; }
  double GetChiSquare1st( void ) const { return m_chisqr1st; }
  double GetChiX( void ) const { return m_Chix; }
  double GetChiY( void ) const { return m_Chiy; }
  double GetChiU( void ) const { return m_Chiu; }
  double GetChiV( void ) const { return m_Chiv; }
  double GetX( double z ) const { return m_x0+m_u0*z; }
  double GetY( double z ) const { return m_y0+m_v0*z; }
  double GetS( double z, double tilt ) const { return GetX(z)*std::cos(tilt)+GetY(z)*std::sin(tilt); }
  int    GetNIteration( void ) const { return m_n_iteration; }
  double GetTheta( void ) const;
  bool   GoodForTracking( void ) const { return m_good_for_tracking; }
  bool   GoodForTracking( bool status )
  { bool ret = m_good_for_tracking; m_good_for_tracking = status; return ret; }
  bool   ReCalc( bool ApplyRecursively=false );
  double GetDe( void ) const { return m_de; }
  void   Print( const std::string& arg="", std::ostream& ost=hddaq::cout ) const;
  void   PrintVXU( const std::string& arg="" ) const;

  // for CFT
  bool   DoFitPhi ( void );
  bool   DoFitPhi2nd( void );
  bool   DoFitPhi16pp( void );
  bool   DoFitUV ( void );
  bool   DoFitUV2nd( int calib ); //0:w/o, 1:w/ UV calib.
  bool   SetCalculatedValueCFT();
  bool   GetCrossPointR(double r, double *phi, int layer);
  void   SetAxy( double Axy){ m_Axy= Axy;}
  void   SetBxy( double Bxy){ m_Bxy= Bxy;}
  double GetAxy( void ) const { return m_Axy; }
  double GetBxy( void ) const { return m_Bxy; }
  double GetChiSquareXY( void ) const { return m_chisqrXY; }
  double GetChiSquareZ ( void ) const { return m_chisqrZ; }
  double calcPhi(double x, double y);
  double CalculateZpos(double phi, DCLTrackHit *cl);
  double GetVtxZ( void ) const { return m_vtx_z; }
  double GetThetaCFT( void ) const { return m_theta; }
  int    GetCFTxyFlag( void ) const { return m_xyFitFlag;  }
  int    GetCFTzFlag ( void ) const { return m_zTrackFlag; }

  double GetPhiIniBefore(  int layer) const {return m_phi_ini_before[layer]  ;}
  double GetPhiTrackBefore(int layer) const {return m_phi_track_before[layer];}
  double GetdPhiBefore(    int layer) const {return m_dphi_before[layer]     ;}
  double GetPhiIni(  int layer) const {return m_phi_ini[layer]  ;}
  double GetPhiTrack(int layer) const {return m_phi_track[layer];}
  double GetdPhi(    int layer) const {return m_dphi[layer]     ;}

  double GetZIniBefore(    int layer) const {return m_z_ini_before[layer]    ;}
  double GetZTrackBefore(  int layer) const {return m_z_track_before[layer]  ;}
  double GetdZBefore(      int layer) const {return m_dz_before[layer]       ;}
  double GetZIni(    int layer) const {return m_z_ini[layer]    ;}
  double GetZTrack(  int layer) const {return m_z_track[layer]  ;}
  double GetdZ(      int layer) const {return m_dz[layer]       ;}
  double GetdR(      int layer) const {return m_dr[layer]     ; }
  double GetR(       int layer) const {return m_r[layer]     ; }

  ThreeVector GetPos0( void ) const { return m_Pos0; }
  ThreeVector GetDir( void )  const { return m_Dir; }
  double GetSumAdc(  int layer) const {return m_sum_adc[layer];   }
  double GetSumMIP(  int layer) const {return m_sum_mip[layer];   }
  double GetSumdE (  int layer) const {return m_sum_dE[layer];   }
  double GetMaxAdc(  int layer) const {return m_max_adc[layer];   }
  double GetMaxMIP(  int layer) const {return m_max_mip[layer];   }
  double GetMaxdE (  int layer) const {return m_max_dE[layer];   }

  double GetTotalSumAdc( void ) const {return m_total_adc;   }
  double GetTotalSumMIP( void ) const {return m_total_mip;   }
  double GetTotalSumdE ( void ) const {return m_total_dE ;   }
  double GetTotalSumdEphi ( void ) const {return m_total_dE_phi ;   }
  double GetTotalSumdEuv  ( void ) const {return m_total_dE_uv ;   }
  double GetTotalMaxAdc( void ) const {return m_total_max_adc;   }
  double GetTotalMaxMIP( void ) const {return m_total_max_mip;   }
  double GetTotalMaxdE ( void ) const {return m_total_max_dE ;   }
  double GetTotalMaxdEphi ( void ) const {return m_total_max_dE_phi ;   }
  double GetTotalMaxdEuv  ( void ) const {return m_total_max_dE_uv ;   }

};


//______________________________________________________________________________
inline
std::ostream&
operator <<( std::ostream& ost,
	     const DCLocalTrack& track )
{
  track.Print( "", ost );
  return ost;
}

//______________________________________________________________________________
struct DCLTrackComp
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1,
		   const DCLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(), chi2=p2->GetChiSquare();
    if( n1>n2+1 )
      return true;
    if( n2>n1+1 )
      return false;
    if( n1<=4 || n2<=4 )
      return ( n1 >= n2 );
    if( n1==n2 )
      return ( chi1 <= chi2 );

    return ( chi1-chi2 <= 3./(n1-4) );// 3-sigma
  }
};

//______________________________________________________________________________
struct DCLTrackComp1
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1,
		   const DCLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
    if(n1>n2) return true;
    if(n2>n1) return false;
    return ( chi1<=chi2 );
  }

};

//______________________________________________________________________________
struct DCLTrackComp2
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1,
		   const DCLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
    if(n1<n2) return true;
    if(n2<n1) return false;
    return ( chi1<=chi2 );
  }

};

//______________________________________________________________________________
struct DCLTrackComp3
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1,
		   const DCLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
    double a1= std::abs(1.-chi1), a2=std::abs(1.-chi2);
    if(a1<a2) return true;
    if(a2<a1) return false;
    return (n1<=n2);
  }

};

//______________________________________________________________________________
struct DCLTrackComp4
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1,
		   const DCLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
    if( (n1>n2+1) && (std::abs(chi1-chi2)<2.) )
      return true;
    if( (n2>n1+1) && ( std::abs(chi1-chi2)<2.) )
      return false;

    return (chi1<=chi2);
  }
};

//______________________________________________________________________________

struct DCLTrackComp_Nhit 
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1,
                   const DCLocalTrack * const p2 ) const
  {
    const int minimum=8;

    int n1=p1->GetNHit(), n2=p2->GetNHit();

    if( n1>=n2 )
      return true;
    else
      return false;
  }
};

//______________________________________________________________________________

struct DCLTrackComp_Chisqr
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1,
                   const DCLocalTrack * const p2 ) const
  {
    double chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();

    if (chi1 <= chi2)
      return true;
    else
      return false;
  }
};


//______________________________________________________________________________
struct DCLTrackCompSdcInFiber
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1,
		   const DCLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
    int NofFiberHit1 = 0;
    int NofFiberHit2 = 0;
    for(int ii=0;ii<n1;ii++){
      int layer = p1->GetHit(ii)->GetLayer();
      if( layer > 6 ) NofFiberHit1++;
    }
    for(int ii=0;ii<n2;ii++){
      int layer = p2->GetHit(ii)->GetLayer();
      if( layer > 6 ) NofFiberHit2++;
    }

    if( (n1>n2+1) ){
      return true;
    }
    else if( (n2>n1+1)  ){
      return false;
    }
    else if( NofFiberHit1 > NofFiberHit2 ){
      return true;
    }
    else if( NofFiberHit2 > NofFiberHit1 ){
      return false;
    }
    else{
      return (chi1<=chi2);
    }

  }
};

//______________________________________________________________________________
struct DCLTrackCompSdcOut
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1,
		   const DCLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(), chi2=p2->GetChiSquare();
    if( n1 > n2 ) return true;
    if( n2 > n1 ) return false;
    return ( chi1 <= chi2 );
  }
};

//______________________________________________________________________________
struct DCLTrackCompCFT
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1,
		   const DCLocalTrack * const p2 ) const
  {
    int nphi1=p1->GetNHit(), nphi2=p2->GetNHit();
    int nuv1=p1->GetNHitUV(), nuv2=p2->GetNHitUV();
    int n1=nphi1+nuv1, n2=nphi2+nuv2;
    
    double chiXY1=p1->GetChiSquareXY(),chiZ1=p2->GetChiSquareZ();
    double chiXY2=p2->GetChiSquareXY(),chiZ2=p2->GetChiSquareZ();
    double chi1=sqrt(chiXY1*chiXY1+chiZ1*chiZ1);
    double chi2=sqrt(chiXY2*chiXY2+chiZ2*chiZ2);
    
    if( (n1==n2)&&(chi1<MaxChi2CFT)) return (chi1<=chi2);
    if( (n1>n2) &&(chi1<MaxChi2CFT)) return true;
    if( (n2>n1) &&(chi2<MaxChi2CFT)) return false;    
    return (chi1<=chi2);


  }
};

//______________________________________________________________________________
struct DCLTrackCompCFT16pp
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1,
		   const DCLocalTrack * const p2 ) const
  {
    
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(), chi2=p2->GetChiSquare();
    if( n1 > n2 ) return true;
    if( n2 > n1 ) return false;
    return ( chi1 <= chi2 );
    
  }
};

#endif
