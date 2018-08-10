/**
 *  file: DCLocalTrack.cc
 *  date: 2017.04.10
 *
 */

#include "DCLocalTrack.hh"

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <stdexcept>
#include <sstream>

#include <std_ostream.hh>

#include "DCAnalyzer.hh"
#include "DCLTrackHit.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"
#include "HodoParamMan.hh"

namespace
{
  const std::string& class_name("DCLocalTrack");
  const DCGeomMan& gGeom = DCGeomMan::GetInstance();
  const double& zK18tgt = gGeom.LocalZ("K18Target");
  const double& zTgt    = gGeom.LocalZ("Target");

  const HodoParamMan& gHodo = HodoParamMan::GetInstance();

  const int ReservedNumOfHits  = 16;
  const int DCLocalMinNHits    =  4;
  const int DCLocalMinNHitsVXU =  2;// for SSD
  const int MaxIteration       = 100;// for Honeycomb
  const double MaxChisqrDiff   = 1.0e-3;
  const int CFTLocalMinNHits   =  3;
}

//______________________________________________________________________________
DCLocalTrack::DCLocalTrack( void )
  : m_is_fitted(false),
    m_is_calculated(false),
    m_Ax(0.), m_Ay(0.), m_Au(0.), m_Av(0.),
    m_Chix(0.), m_Chiy(0.), m_Chiu(0.), m_Chiv(0.),
    m_x0(0.), m_y0(0.),
    m_u0(0.), m_v0(0.),
    m_a(0.),  m_b(0.),
    m_chisqr(1.e+10),
    m_good_for_tracking(true),
    m_de(0.),
    m_chisqr1st(1.e+10),
    m_Axy(0.), m_Bxy(0.), m_Az(0.), m_Bz(0.),
    m_chisqrXY(1.e+10),
    m_chisqrZ(1.e+10),
    m_vtx_z(-999.),
    m_xyFitFlag(-1),
    m_zTrackFlag(-1),
    m_theta(-999),
    m_n_iteration(0)
{
  m_hit_array.reserve( ReservedNumOfHits );
  m_hit_arrayUV.reserve( ReservedNumOfHits );
  debug::ObjectCounter::increase(class_name);

  for(int i=0; i<NumOfPlaneCFT*2; i++){
    m_meanseg[i]   = -999;
    m_phi_ini[i]   = -999;  
    m_phi_track[i] = -999;
    m_dphi[i]      = -999;        
    m_z_ini[i]     = -999;    
    m_z_track[i]   = -999;  
    m_dz[i]        = -999;           
    m_dist_phi[i]  = -999;           
  }

}

//______________________________________________________________________________
DCLocalTrack::~DCLocalTrack( void )
{
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void
DCLocalTrack::AddHit( DCLTrackHit *hitp )
{
  if( hitp )
    m_hit_array.push_back( hitp );
}
//______________________________________________________________________________
void
DCLocalTrack::AddHitUV( DCLTrackHit *hitp )
{
  if( hitp )
    m_hit_arrayUV.push_back( hitp );
}

//______________________________________________________________________________
void
DCLocalTrack::Calculate( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( IsCalculated() ){
    hddaq::cerr << "#W " << func_name << " "
		<< "already called" << std::endl;
    return;
  }

  const std::size_t n = m_hit_array.size();
  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = m_hit_array[i];
    double z0 = hitp->GetZ();
    hitp->SetCalPosition( GetX(z0), GetY(z0) );
    hitp->SetCalUV( m_u0, m_v0 );
    // for Honeycomb
    if( !hitp->IsHoneycomb() )
      continue;
    double scal = hitp->GetLocalCalPos();
    double wp   = hitp->GetWirePosition();
    double dl   = hitp->GetDriftLength();
    double ss   = scal-wp>0 ? wp+dl : wp-dl;
    hitp->SetLocalHitPos( ss );
  }
  m_is_calculated = true;
}

//______________________________________________________________________________
void
DCLocalTrack::CalculateCFT( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( IsCalculated() ){
    hddaq::cerr << "#W " << func_name << " "
		<< "already called" << std::endl;
    return;
  }

  const std::size_t n = m_hit_array.size();
  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = m_hit_array[i];
    double z0 = hitp->GetZ();
    hitp->SetCalPosition( GetX(z0), GetY(z0) );
    hitp->SetCalUV( m_u0, m_v0 );    
    // for Honeycomb
    if( !hitp->IsHoneycomb() )
      continue;
    /*
    double scal = hitp->GetLocalCalPos();
    double wp   = hitp->GetWirePosition();
    double dl   = hitp->GetDriftLength();
    double ss   = scal-wp>0 ? wp+dl : wp-dl;
    */
    //hitp->SetLocalHitPos( ss );
    hitp->SetLocalHitPos( 0 );    
  }

  const std::size_t nUV = m_hit_arrayUV.size();
  for( std::size_t i=0; i<nUV; ++i ){
    DCLTrackHit *hitp = m_hit_arrayUV[i];
    double z0 = hitp->GetZ();
    hitp->SetCalPosition( GetX(z0), GetY(z0) );
    hitp->SetCalUV( m_u0, m_v0 );    
    // for Honeycomb
    if( !hitp->IsHoneycomb() )
      continue;
    hitp->SetLocalHitPos( 0 );    
  }

  m_is_calculated = true;
}

//______________________________________________________________________________
int
DCLocalTrack::GetNDF( void ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const std::size_t n = m_hit_array.size();
  int ndf = 0;
  for( std::size_t i=0; i<n; ++i ){
    if( m_hit_array[i] ) ++ndf;
  }
  return ndf-4;
}

//______________________________________________________________________________
int
DCLocalTrack::GetNHitSFT(void) const
{
  int n_sft=0;
  for(const auto& hit : m_hit_array){
    if(hit->GetLayer() > 6) ++n_sft;
  }

  return n_sft;
}

//______________________________________________________________________________
int
DCLocalTrack::GetNHitY(void) const
{
  int n_y=0;
  for(const auto& hit : m_hit_array){
    if(hit->GetTiltAngle() > 1) ++n_y;
  }

  return n_y;
}

//______________________________________________________________________________
DCLTrackHit*
DCLocalTrack::GetHit( std::size_t nth ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( nth<m_hit_array.size() )
    return m_hit_array[nth];
  else
    return 0;
}

//______________________________________________________________________________
DCLTrackHit*
DCLocalTrack::GetHitUV( std::size_t nth ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( nth<m_hit_arrayUV.size() )
    return m_hit_arrayUV[nth];
  else
    return 0;
}

//______________________________________________________________________________
DCLTrackHit*
DCLocalTrack::GetHitOfLayerNumber( int lnum ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const std::size_t n = m_hit_array.size();
  for( std::size_t i=0; i<n; ++i ){
    if( m_hit_array[i]->GetLayer()==lnum )
      return m_hit_array[i];
  }
  return 0;
}

//______________________________________________________________________________
void
DCLocalTrack::DeleteNullHit( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( std::size_t i=0; i<m_hit_array.size(); ++i ){
    DCLTrackHit *hitp = m_hit_array[i];
    if( !hitp ){
      hddaq::cout << func_name << " "
		  << "null hit is deleted" << std::endl;
      m_hit_array.erase( m_hit_array.begin()+i );
      --i;
    }
  }
}

//______________________________________________________________________________
bool
DCLocalTrack::DoFit( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( IsFitted() ){
    hddaq::cerr << "#W " << func_name << " "
		<< "already called" << std::endl;
    return false;
  }

  DeleteNullHit();

  const std::size_t n = m_hit_array.size();
  if( n < DCLocalMinNHits ){
    hddaq::cout << "#D " << func_name << " "
		<< "the number of layers is too small : " << n << std::endl;
    return false;
  }

  const int nItr = HasHoneycomb() ? MaxIteration : 1;

  double prev_chisqr = m_chisqr;
  std::vector <double> z0(n), z(n), wp(n),
    w(n), s(n), ct(n), st(n), coss(n);
  std::vector<bool> honeycomb(n);
  for( int iItr=0; iItr<nItr; ++iItr ){
    for( std::size_t i=0; i<n; ++i ){
      DCLTrackHit *hitp = m_hit_array[i];
      int    lnum = hitp->GetLayer();
      honeycomb[i] = hitp->IsHoneycomb();
      wp[i] = hitp->GetWirePosition();
      z0[i] = hitp->GetZ();
      double ww = gGeom.GetResolution( lnum );
      w[i] = 1./(ww*ww);
      double aa = hitp->GetTiltAngle()*math::Deg2Rad();
      ct[i] = std::cos(aa); st[i] = std::sin(aa);
      double ss = hitp->GetLocalHitPos();
      double dl = hitp->GetDriftLength();
      double dsdz = m_u0*std::cos(aa)+m_v0*std::sin(aa);
      double dcos = std::cos( std::atan(dsdz) );
      coss[i] = dcos;
      double dsin = std::sin( std::atan(dsdz) );
      double ds = dl * dcos;
      double dz = dl * dsin;
      double scal = iItr==0 ? ss : GetS(z[i],aa);
      if( honeycomb[i] ){
	s[i] = scal-wp[i]>0 ? wp[i]+ds : wp[i]-ds;
	z[i] = scal-wp[i]>0 ? z0[i]-dz : z0[i]+dz;
      }else{
	s[i] = ss;
	z[i] = z0[i];
      }
    }

    double x0, u0, y0, v0;
    if( !math::SolveGaussJordan( z, w, s, ct, st,
				 x0, u0, y0, v0 ) ){
      hddaq::cerr << func_name << " Fitting failed" << std::endl;
      return false;
    }

    double chisqr = 0.;
    double de     = 0.;
    for( std::size_t i=0; i<n; ++i ){
      double scal = (x0+u0*z0[i])*ct[i]+(y0+v0*z0[i])*st[i];
      double ss   = wp[i]+(s[i]-wp[i])/coss[i];
      double res  = honeycomb[i] ? (ss-scal)*coss[i] : s[i]-scal;
      chisqr += w[i]*res*res;
    }
    chisqr /= GetNDF();

    if( iItr==0 ) m_chisqr1st = chisqr;

    // if worse, not update
    if( prev_chisqr-chisqr>0. ){
      m_x0 = x0;
      m_y0 = y0;
      m_u0 = u0;
      m_v0 = v0;
      m_chisqr = chisqr;
      m_de     = de;
    }

    // judge convergence
    if( prev_chisqr-chisqr<MaxChisqrDiff ){
#if 0
      // if( chisqr<200. && GetTheta()>4. )
      if( chisqr<20. )
      {
	if( iItr==0 ) hddaq::cout << "=============" << std::endl;
	hddaq::cout << func_name << " NIteration : " << iItr << " "
		    << "chisqr = " << std::setw(10) << std::setprecision(4)
		    << m_chisqr << " "
		    << "diff = " << std::setw(20) << std::left
		    << m_chisqr-m_chisqr1st << " ndf = " << GetNDF() << std::endl;
      }
#endif
      m_n_iteration = iItr;
      break;
    }

    prev_chisqr = chisqr;
  }

  m_is_fitted = true;
  return true;
}

//______________________________________________________________________________
bool
DCLocalTrack::DoFitBcSdc( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( IsFitted() ){
    hddaq::cerr << "#W " << func_name << " "
		<< "already called" << std::endl;
    return false;
  }

  DeleteNullHit();

  std::size_t n = m_hit_array.size();
  if( n < DCLocalMinNHits ) return false;

  std::vector <double> z, w, s, ct, st;
  z.reserve(n); w.reserve(n); s.reserve(n);
  ct.reserve(n); st.reserve(n);

  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = m_hit_array[i];
    int lnum = hitp->GetLayer();
    double ww = gGeom.GetResolution( lnum );
    double zz = hitp->GetZ();
    if ( lnum>=113&&lnum<=124 ){ // BcOut
      zz -= zK18tgt - zTgt;
    }
    if ( lnum>=11&&lnum<=18 ){ // SsdIn/Out
      zz -= zK18tgt - zTgt;
    }
    double aa = hitp->GetTiltAngle()*math::Deg2Rad();
    z.push_back( zz ); w.push_back( 1./(ww*ww) );
    s.push_back( hitp->GetLocalHitPos() );
    ct.push_back( cos(aa) ); st.push_back( sin(aa) );
  }

  if( !math::SolveGaussJordan( z, w, s, ct, st,
			       m_x0, m_u0, m_y0, m_v0 ) ){
    hddaq::cerr << func_name << " Fitting fails" << std::endl;
    return false;
  }

  double chisqr = 0.;
  for( std::size_t i=0; i<n; ++i ){
    double ww=w[i], zz=z[i];
    double scal=GetX(zz)*ct[i]+GetY(zz)*st[i];
    chisqr += ww*(s[i]-scal)*(s[i]-scal);
  }
  chisqr /= GetNDF();
  m_chisqr = chisqr;

  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = m_hit_array[i];
    if( hitp ){
      int lnum = hitp->GetLayer();
      double zz = hitp->GetZ();
      if ( (lnum>=113&&lnum<=124) ){ // BcOut
        zz -= zK18tgt - zTgt;
      }
      if ( (lnum>=11&&lnum<=18) ){ // SsdIn/Out
        zz -= zK18tgt - zTgt;
      }
      hitp->SetCalPosition( GetX(zz), GetY(zz) );
    }
  }

  m_is_fitted = true;
  return true;
}

//______________________________________________________________________________
bool
DCLocalTrack::DoFitVXU( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( IsFitted() ){
    hddaq::cerr << "#W " << func_name << " "
		<< "already called" << std::endl;
    return false;
  }

  DeleteNullHit();

  const std::size_t n = m_hit_array.size();
  if( n<DCLocalMinNHitsVXU ) return false;

  double w[n+1],z[n+1],x[n+1];

  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = m_hit_array[i];
    int lnum = hitp->GetLayer();
    w[i] = gGeom.GetResolution( lnum );
    z[i] = hitp->GetZ();
    x[i] = hitp->GetLocalHitPos();
#if 0
    hddaq::cout << "" << std::endl;
    hddaq::cout << "**********" << std::endl;
    hddaq::cout << std::setw(10) << "layer = " << lnum
		<< std::setw(10) << "wire  = " << hitp->GetWire() << " "
		<< std::setw(20) << "WirePosition = "<<hitp->GetWirePosition() << " "
		<< std::setw(20) << "DriftLength = "<<hitp->GetDriftLength() << " "
		<< std::setw(20) << "hit position = "<<hitp->GetLocalHitPos()<< " "
		<< std::endl;
    hddaq::cout << "**********" << std::endl;
    hddaq::cout << "" << std::endl;
#endif
  }

  double A=0, B=0, C=0, D=0, E=0;// <-Add!!
  for( std::size_t i=0; i<n; ++i ){
    A += z[i]/(w[i]*w[i]);
    B += 1/(w[i]*w[i]);
    C += x[i]/(w[i]*w[i]);
    D += z[i]*z[i]/(w[i]*w[i]);
    E += x[i]*z[i]/(w[i]*w[i]);
  }

  m_a = (E*B-C*A)/(D*B-A*A);
  m_b = (D*C-E*A)/(D*B-A*A);

  double chisqr = 0.;
  for( std::size_t i=0; i<n; ++i ){
    chisqr += (x[i]-m_a*z[i]-m_b)*(x[i]-m_a*z[i]-m_b)/(w[i]*w[i]);
  }

  if(n==2) chisqr  = 0.;
  else     chisqr /= n-2.;
  m_chisqr = chisqr;

  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = m_hit_array[i];
    if( hitp ){
      double zz = hitp->GetZ();
      hitp->SetLocalCalPosVXU( m_a*zz+m_b );
    }
  }

  m_is_fitted = true;
  return true;
}

// for CFT
//______________________________________________________________________________
bool
DCLocalTrack::DoFitPhi( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( IsFitted() ){
    hddaq::cerr << "#W " << func_name << " "
		<< "already called" << std::endl;
    return false;
  }
  DeleteNullHit();

  const std::size_t n = m_hit_array.size();
  if( n < CFTLocalMinNHits ){
    hddaq::cout << "#D " << func_name << " "
    		<< "the number of layers is too small : " << n << std::endl;
    return false;
  }

  std::vector <double> x, y, s;
  std::vector <int> l;
  x.reserve(n); y.reserve(n); s.reserve(n);

  // initial Phi 
  double iniPhi[n]; 
  double delta_Phi[n]; 
  double Layer[n]; 
  double sumPhi = 0.;
  double Phi0 = -999.;

  for(int jj=0; jj<n; ++jj ){   
    DCLTrackHit *ini_hitp = m_hit_array[jj];
    Layer[jj] = ini_hitp->GetLayer();    
    iniPhi[jj] = ini_hitp->GetPositionPhi();
    sumPhi += iniPhi[jj];    

    m_meanseg[(int)Layer[jj]] = ini_hitp->GetMeanSeg();
    m_phi_ini[(int)Layer[jj]] = iniPhi[jj];
  }
  Phi0 = sumPhi/(double)n;
  double Phi0_ = Phi0;

  if(Phi0<=40||Phi0>=320){ 
    Phi0+=180.; 
    if(Phi0>360){Phi0-=360.;}
  }

  if(n<=4){// normal tracking
    // initial Phi cut
    for(int k=0; k<n; ++k ){      
      if(Phi0_<=40||Phi0_>=320){ 
	iniPhi[k]+=180.;	
	if(iniPhi[k]>360){iniPhi[k]-=360.;}
      }
      delta_Phi[k] = iniPhi[k] -Phi0;
      if( fabs(delta_Phi[k])>70 && fabs(delta_Phi[k])<270 ){ 	  
	return false;
      }      
    }
  }

  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = m_hit_array[i];
    if( hitp ){
      int lnum = hitp->GetLayer();

      double rr = hitp->GetPositionR();
      double pphi = hitp->GetPositionPhi()*math::Deg2Rad();      
      //double rr = hitp->GetPositionRMax();
      //double pphi = hitp->GetPositionPhiMax()*math::Deg2Rad();
      double xx = rr*cos(pphi);
      double yy = rr*sin(pphi);
      double ss = 0.75/sqrt(12);  

      x.push_back( xx ); y.push_back( yy ); s.push_back(ss);
      l.push_back( lnum );

#if 0
      double phi = hitp->GetPositionPhi();
      std::cout << std::setw(10) << "layer = " << lnum 
		<< std::setw(10) << "x  = " << xx << ", y  =  " << yy
		<< ", phi = " << pphi*math::Rad2Deg()
		<< std::endl;
#endif
    }
  }

  double A=0., B=0., C=0., D=0., E=0., F=0.;

  // y = ax + b
  for (int i=0; i<n; i++) {
    A += x[i]/(s[i]*s[i]);
    B += 1./(s[i]*s[i]);
    C += y[i]/(s[i]*s[i]);
    D += x[i]*x[i]/(s[i]*s[i]);
    E += x[i]*y[i]/(s[i]*s[i]);
    F += y[i]*y[i]/(s[i]*s[i]);
  }

  m_Axy=(E*B-C*A)/(D*B-A*A);
  m_Bxy=(D*C-E*A)/(D*B-A*A);
  
#if 0
  std::cout << "Axy = " << m_Axy
	    << "Bxy = " << m_Bxy << ", phi = " << atan(m_Axy)*math::Rad2Deg()
	    << std::endl;  
#endif  

  if (fabs(m_Axy) <= 1.) {  //y  = ax + b
    m_xyFitFlag=0;
    m_chisqrXY = 0.;
    for (int i=0; i<n; i++) {
      double ycal = m_Axy*x[i] + m_Bxy;
      m_chisqrXY += (y[i]-ycal)*(y[i]-ycal)/(s[i]*s[i]);

      // +phi or -phi             
      if(m_Axy>=0){
	if(x[i]>=0){
	  m_dist_phi[l[i]] = (y[i]-m_Axy*x[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}else{
	  m_dist_phi[l[i]] = -(y[i]-m_Axy*x[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}
      }else{
	if(x[i]<=0){
	  m_dist_phi[l[i]] = -(y[i]-m_Axy*x[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}else{
	  m_dist_phi[l[i]] = (y[i]-m_Axy*x[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}
      }      

    }
    m_chisqrXY /= n;

  }else{ // x = ay + b    
    m_xyFitFlag=1;

    A=0.; B=0.; C=0.; D=0.; E=0.; F=0.;
    for (int i=0; i<n; i++) {
      A += y[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += x[i]/(s[i]*s[i]);
      D += y[i]*y[i]/(s[i]*s[i]);
      E += x[i]*y[i]/(s[i]*s[i]);
      F += x[i]*x[i]/(s[i]*s[i]);
    }
    m_Axy=(E*B-C*A)/(D*B-A*A);
    m_Bxy=(D*C-E*A)/(D*B-A*A);

    m_chisqrXY = 0.;
    for (int i=0; i<n; i++) {
      double xcal = m_Axy*y[i] + m_Bxy;
      m_chisqrXY += (x[i]-xcal)*(x[i]-xcal)/(s[i]*s[i]);

      // +phi or -phi 
      if(m_Axy>=0){
	if(y[i]>=0){
	  m_dist_phi[l[i]] = -(x[i]-m_Axy*y[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}else{
	  m_dist_phi[l[i]] = (x[i]-m_Axy*y[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}
      }else{
	if(y[i]<=0){
	  m_dist_phi[l[i]] = (x[i]-m_Axy*y[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}else{
	  m_dist_phi[l[i]] = -(x[i]-m_Axy*y[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}
      }

    }
    m_chisqrXY /= n;

  }

  return true;
}

// 1st spiral layer tracking
bool DCLocalTrack::DoFitUV()
{
  const std::string funcname = "[CFTLocalTrack::DoFitUV()]";
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t n = m_hit_arrayUV.size();

  if(n < CFTLocalMinNHits ) return false;


  // initial z 
  double iniZ[n]; 
  double delta_z[n]; 
  double delta_dz[n]; 
  double sum_dz=0; 
  double c_dz=0; 
  double delta_n=0; 
  int Layer[n]; 
  double rrr[n]; 
  double phi;
  
  for(int jj=0; jj<n; ++jj ){
    DCLTrackHit *ini_hitp = m_hit_arrayUV[jj];
    Layer[jj] = ini_hitp->GetLayer();
    m_meanseg[(int)Layer[jj]] = ini_hitp->GetMeanSeg();
    int l_geo = Layer[jj] +301;
    if(Layer[jj]>7){l_geo -= 8;}
    //double r = ini_hitp->GetRMax();
    double r = ini_hitp->GetPositionR();

    bool ret = GetCrossPointR(r, &phi, Layer[jj]);
    if (!ret) {
      //      std::cout << funcname << " GetCrossPointR false r=" << r
      //		<< ", phi=" << phi << ", phi0=" << Phi0_ << std::endl;
      return false;
    }    
    
    iniZ[jj] = CalculateZpos(phi, ini_hitp);
    m_phi_ini[(int)Layer[jj]] = phi;
    
  }

#if 0  
  // z slope cut
  for(int k=0; k<n-1; ++k ){
    delta_z[k] = (iniZ[k+1]-iniZ[k])/(Layer[k+1]-Layer[k]) *2;     
    std::cout << "delta_Z=" << delta_z[k] 
    	      <<", Layer " << Layer[k] 
    	      << " and " << Layer[k+1] << std::endl ;    
    sum_dz+=delta_z[k];
  }
  c_dz = sum_dz/((double)n-1);
  std::cout << " c_dz=" << c_dz << std::endl ;
#endif  
  
  std::vector <double> z, xy, s;
  z.reserve(n); xy.reserve(n); s.reserve(n);
  

  for( std::size_t i=0; i<n; ++i ){
    
    DCLTrackHit *hitp = m_hit_arrayUV[i];
    if( hitp ){
      int lnum = hitp->GetLayer();
      int l_geo = lnum +301;
      if(lnum>7){l_geo -= 8;}
      double r = hitp->GetPositionR();
      
      double phi;
      bool ret = GetCrossPointR(r, &phi, lnum);
      if (!ret) {
	//std::cout << funcname << " return at GetCrossPointR" << std::endl;
	return false;
      }

      double z1 = CalculateZpos(phi, hitp);
      m_z_ini[lnum] = z1; 

      double xy1=-999.;
      if (m_xyFitFlag==0) { // x
	xy1 = r*cos(phi*math::Deg2Rad());
      }else{ // y
	xy1 = r*sin(phi*math::Deg2Rad());
      }
      
      double ss = geomMan.GetResolution( l_geo );
      z.push_back( z1 ); xy.push_back( xy1 ); s.push_back(ss);
      rrr[i] = r;
    
    
#if 0
    std::cout << std::setw(10) << "layer = " << lnum 
	      << std::setw(10) << "z  = " << z1 << ", xy  =  " << xy1
	      << ", phi = " << phi
	      << std::endl;
#endif    
    }
  }
  
  
  double slope = (xy[n-1]-xy[0])/(z[n-1]-z[0]);
  if (fabs(slope)<=1) {
    // x = az + b
    m_zTrackFlag = 0;
    
    double A=0., B=0., C=0., D=0., E=0., F=0.;

    for (int i=0; i<n; i++) {
      A += z[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += xy[i]/(s[i]*s[i]);
      D += z[i]*z[i]/(s[i]*s[i]);
      E += z[i]*xy[i]/(s[i]*s[i]);
      F += xy[i]*xy[i]/(s[i]*s[i]);
    }
    
    m_Az=(E*B-C*A)/(D*B-A*A);
    m_Bz=(D*C-E*A)/(D*B-A*A);
    
    if (m_xyFitFlag==0) {
      m_u0 = m_Az;
      m_x0 = m_Bz;
      m_v0 = m_Axy*m_Az;
      m_y0 = m_Axy*m_Bz+m_Bxy;
    } else {
      m_u0 = m_Axy*m_Az;
      m_x0 = m_Axy*m_Bz+m_Bxy;
      m_v0 = m_Az;
      m_y0 = m_Bz;
    }
    
    m_chisqrZ = 0.;
    for (int i=0; i<n; i++) {
      double xycal = m_Az*z[i] + m_Bz;
      double xcal = m_u0*z[i] + m_x0;
      double ycal = m_v0*z[i] + m_y0;

      double zcal = (xy[i]-m_Bz)/m_Az;      	      
      m_z_track[Layer[i]] = zcal;
      
      m_chisqrZ += (xy[i]-xycal)*(xy[i]-xycal)/(s[i]*s[i]);
      
      // inner(-) or outer(+)       
      if(m_xyFitFlag==0){// x = az + b
	if(xy[i]>=0){
	  m_dist_uv[Layer[i]] = (xy[i]-m_Az*z[i]-m_Bz)/sqrt(1+m_Az*m_Az);
	}else{
	  m_dist_uv[Layer[i]] = -(xy[i]-m_Az*z[i]-m_Bz)/sqrt(1+m_Az*m_Az);
	}
      }else{ // y = az + b
	if(xy[i]<=0){
	  m_dist_uv[Layer[i]] = -(xy[i]-m_Az*z[i]-m_Bz)/sqrt(1+m_Az*m_Az);
	}else{
	  m_dist_uv[Layer[i]] = (xy[i]-m_Az*z[i]-m_Bz)/sqrt(1+m_Az*m_Az);
	}
      }
      
    }
    m_chisqrZ /= n;   
    m_vtx_z = -m_Bz/m_Az;
    


  } else {
    // z = ax + b or z = ay + b
    m_zTrackFlag = 1;    
    
    double A=0., B=0., C=0., D=0., E=0., F=0.;
    for (int i=0; i<n; i++) {
      A += xy[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += z[i]/(s[i]*s[i]);
      D += xy[i]*xy[i]/(s[i]*s[i]);
      E += z[i]*xy[i]/(s[i]*s[i]);
      F += z[i]*z[i]/(s[i]*s[i]);
    }

    m_Az=(E*B-C*A)/(D*B-A*A);
    m_Bz=(D*C-E*A)/(D*B-A*A);

    if (m_xyFitFlag==0) {
      m_u0 = 1./m_Az;
      m_x0 = -m_Bz/m_Az;
      m_v0 = m_Axy/m_Az;
      m_y0 = -m_Axy*m_Bz/m_Az+m_Bxy;
    } else {
      m_u0 = m_Axy/m_Az;
      m_x0 = -m_Axy*m_Bz/m_Az+m_Bxy;
      m_v0 = 1./m_Az;
      m_y0 = -m_Bz/m_Az;
    }

    m_chisqrZ = 0.;
    for (int i=0; i<n; i++) {
      double zcal = m_Az*xy[i] + m_Bz;
      m_z_track[Layer[i]] = zcal;
      m_chisqrZ += (z[i]-zcal)*(z[i]-zcal)/(s[i]*s[i]);

      // realistic       
      m_dist_uv[Layer[i]] = -(z[i]-m_Az*xy[i]-m_Bz)/sqrt(1+m_Az*m_Az);

    }
    m_chisqrZ /= n;    
    m_vtx_z = m_Bz;
    
  }  
  
  return true;
}

// for CFT 2nd tracking (position correction)
//______________________________________________________________________________
bool
DCLocalTrack::DoFitPhi2nd( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( IsFitted() ){
    hddaq::cerr << "#W " << func_name << " "
		<< "already called" << std::endl;
    return false;
  }
  DeleteNullHit();

  const std::size_t n = m_hit_array.size();
  if( n < CFTLocalMinNHits ){
    hddaq::cout << "#D " << func_name << " "
    		<< "the number of layers is too small : " << n << std::endl;
    return false;
  }

  std::vector <double> x, y, s;
  std::vector <int> l;
  x.reserve(n); y.reserve(n); s.reserve(n);

  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = m_hit_array[i];
    if( hitp ){
      int lnum = hitp->GetLayer();
      int seg  = (int)hitp->GetMeanSeg();
      m_meanseg[lnum] = seg;

      double rr = hitp->GetPositionR();
      //double rr = hitp->GetPositionRMax();
      double xx = rr*cos(m_phi_ini[lnum]*math::Deg2Rad());
      double yy = rr*sin(m_phi_ini[lnum]*math::Deg2Rad());
 
      
      double ss = 0.75/sqrt(12);  

      double zz=-1500.;
      // z dependence
      if(m_zTrackFlag==0){
	//(x or y) = az + b
	if(m_xyFitFlag==0){
	  zz = (xx - m_Bz)/m_Az;
	}else if(m_xyFitFlag==1){
	  zz = (yy - m_Bz)/m_Az;
	}
      }else if(m_zTrackFlag==1){
	//z = a*(x or y) + b
	if(m_xyFitFlag==0){
	  zz = m_Az*xx + m_Bz;
	}else if(m_xyFitFlag==1){
	  zz = m_Az*yy + m_Bz;
	}
      }      
      m_z_track[lnum] = zz;

      // parameter for position correction
      int ll = lnum;
      if(lnum>7){lnum -= 8;}
      double par[6], par2nd[6];
      for(int ip=0; ip<6; ip++){
	par[ip] = gHodo.GetPar(113,lnum,seg,0,ip);	
	//par2nd[ip] = gHodo.GetPar(113,lnum,seg,1,ip);
      }      



#if 1 // correction depending on z
      if(zz>=0&&zz<400){
	double d_phi = par[0]*pow(zz,0) +par[1]*pow(zz,1) 
	  +par[2]*pow(zz,2) +par[3]*pow(zz,3)
	  +par[4]*pow(zz,4) +par[5]*pow(zz,5) ;
	if(fabs(d_phi)>6){
	  std::cout << "layer=" << lnum <<  ",seg=" << m_meanseg[lnum]
		    <<  ", z=" << zz <<  ", d_phi = " <<d_phi  << std::endl;
	}else{
	  m_phi_ini[lnum] += d_phi; //axis test
	}
      }else{}
#endif
      
      // re:calc
      xx = rr*cos(m_phi_ini[lnum]*math::Deg2Rad() );
      yy = rr*sin(m_phi_ini[lnum]*math::Deg2Rad() );
            
      x.push_back( xx ); y.push_back( yy ); s.push_back(ss);
      l.push_back( lnum );

#if 0
      double phi = hitp->GetPositionPhi();
      std::cout << std::setw(10) << "layer = " << lnum 
		<< std::setw(10) << "x  = " << xx << ", y  =  " << yy
		<< ", phi = " << pphi*math::Rad2Deg()
		<< std::endl;
      std::cout << "param0  layer" <<lnum << ", seg=" <<seg 
		<< " par0=" <<par[0] << ", par1=" <<par[1] 
		<< ", par2=" <<par[2] << ", par3=" <<par[3]
		<< ", par4=" <<par[4] << ", par5=" <<par[5] << std::endl;
#endif
    }
  }

  double A=0., B=0., C=0., D=0., E=0., F=0.;

  // y = ax + b
  for (int i=0; i<n; i++) {
    A += x[i]/(s[i]*s[i]);
    B += 1./(s[i]*s[i]);
    C += y[i]/(s[i]*s[i]);
    D += x[i]*x[i]/(s[i]*s[i]);
    E += x[i]*y[i]/(s[i]*s[i]);
    F += y[i]*y[i]/(s[i]*s[i]);
  }

  m_Axy=(E*B-C*A)/(D*B-A*A);
  m_Bxy=(D*C-E*A)/(D*B-A*A);
  
#if 0
  std::cout << "Axy = " << m_Axy
	    << "Bxy = " << m_Bxy << ", phi = " << atan(m_Axy)*math::Rad2Deg()
	    << std::endl;  
#endif  

  if (fabs(m_Axy) <= 1.) {  //y  = ax + b
    m_xyFitFlag=0;
    m_chisqrXY = 0.;
    for (int i=0; i<n; i++) {
      double ycal = m_Axy*x[i] + m_Bxy;
      m_chisqrXY += (y[i]-ycal)*(y[i]-ycal)/(s[i]*s[i]);

      // +phi or -phi             
      if(m_Axy>=0){
	if(x[i]>=0){
	  m_dist_phi[l[i]] = (y[i]-m_Axy*x[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}else{
	  m_dist_phi[l[i]] = -(y[i]-m_Axy*x[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}
      }else{
	if(x[i]<=0){
	  m_dist_phi[l[i]] = -(y[i]-m_Axy*x[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}else{
	  m_dist_phi[l[i]] = (y[i]-m_Axy*x[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}
      }      

    }
    m_chisqrXY /= n;

  }else{ // x = ay + b    
    m_xyFitFlag=1;

    A=0.; B=0.; C=0.; D=0.; E=0.; F=0.;
    for (int i=0; i<n; i++) {
      A += y[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += x[i]/(s[i]*s[i]);
      D += y[i]*y[i]/(s[i]*s[i]);
      E += x[i]*y[i]/(s[i]*s[i]);
      F += x[i]*x[i]/(s[i]*s[i]);
    }
    m_Axy=(E*B-C*A)/(D*B-A*A);
    m_Bxy=(D*C-E*A)/(D*B-A*A);

    m_chisqrXY = 0.;
    for (int i=0; i<n; i++) {
      double xcal = m_Axy*y[i] + m_Bxy;
      m_chisqrXY += (x[i]-xcal)*(x[i]-xcal)/(s[i]*s[i]);

      // +phi or -phi 
      if(m_Axy>=0){
	if(y[i]>=0){
	  m_dist_phi[l[i]] = -(x[i]-m_Axy*y[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}else{
	  m_dist_phi[l[i]] = (x[i]-m_Axy*y[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}
      }else{
	if(y[i]<=0){
	  m_dist_phi[l[i]] = (x[i]-m_Axy*y[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}else{
	  m_dist_phi[l[i]] = -(x[i]-m_Axy*y[i]-m_Bxy)/sqrt(1+m_Axy*m_Axy);
	}
      }

    }
    m_chisqrXY /= n;

  }

  return true;
}

// 2nd spiral layer tracking
bool DCLocalTrack::DoFitUV2nd()
{
  const std::string funcname = "[CFTLocalTrack::DoFitUV2nd()]";
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t n = m_hit_arrayUV.size();
  if(n < CFTLocalMinNHits ) return false;

  // initial z 
  int Layer[n];   
  std::vector <double> z, xy, s;
  z.reserve(n); xy.reserve(n); s.reserve(n);
  

  for( std::size_t i=0; i<n; ++i ){
    
    DCLTrackHit *hitp = m_hit_arrayUV[i];
    if( hitp ){
      int lnum = hitp->GetLayer();
      Layer[i] = lnum;
      m_meanseg[lnum] = hitp->GetMeanSeg();
      int seg = (int)m_meanseg[lnum];
      int l_geo = lnum +301;
      if(lnum>7){l_geo -= 8;}
      double r = hitp->GetPositionR();
      
      double phi;
      bool ret = GetCrossPointR(r, &phi, lnum);
      if (!ret) {
	//std::cout << funcname << " return at GetCrossPointR" << std::endl;
	return false;
      }

      double z1 = CalculateZpos(phi, hitp);
      m_z_ini[lnum] = z1;      
      m_phi_track[lnum] = phi;    

      // parameter for position correction
      int ll = lnum;
      if(lnum>7){lnum -= 8;}
      double par[6], par2nd[6];;
      for(int ip=0; ip<6; ip++){
	par[ip] = gHodo.GetPar(113,lnum,seg,0,ip);	
	par2nd[ip] = gHodo.GetPar(113,lnum,seg,1,ip);
      }      

#if 1 // position correction
      // divide <180 or >=180
      //double phis=phi;
      double phis=-(phi-90.); // axis test
      if(phis<0){phis+=360.;}
      else if(phis>360){phis-=360.;}
      double d_z = 0;
      //if(phi<180){// use 1st par
      if(phis<180){// use 1st par
	d_z = par[0]*pow(phis,0) +par[1]*pow(phis,1) 
	    +par[2]*pow(phis,2) +par[3]*pow(phis,3)
	  +par[4]*pow(phis,4) +par[5]*pow(phis,5) ;	  
      }else{// use 2nd par
	d_z = par2nd[0]*pow(phis,0) +par2nd[1]*pow(phis,1) 
	  +par2nd[2]*pow(phis,2) +par2nd[3]*pow(phis,3)
	  +par2nd[4]*pow(phis,4) +par2nd[5]*pow(phis,5) ;
      }
      if(fabs(d_z)>6){z1 += (d_z/2.);}
      else{z1 += d_z;}            
      //if(fabs(d_z)>6){z1 -= (d_z/2.);} // axis test
      //else{z1 -= d_z;} // axis test
#endif      
      m_z_ini[lnum] = z1;      

      double xy1=-999.;
      if (m_xyFitFlag==0) { // x
	xy1 = r*cos(m_phi_track[lnum]*math::Deg2Rad());
      }else{ // y
	xy1 = r*sin(m_phi_track[lnum]*math::Deg2Rad());
      }      
      double ss = geomMan.GetResolution( l_geo );
      z.push_back( z1 ); xy.push_back( xy1 ); s.push_back(ss);
    
    
#if 0
    std::cout << std::setw(10) << "layer = " << lnum 
	      << std::setw(10) << "z  = " << z1 << ", xy  =  " << xy1
	      << ", phi = " << phi
	      << std::endl;
#endif    
    }
  }
  
  
  double slope = (xy[n-1]-xy[0])/(z[n-1]-z[0]);
  if (fabs(slope)<=1) {
    // x = az + b
    m_zTrackFlag = 0;
    
    double A=0., B=0., C=0., D=0., E=0., F=0.;

    for (int i=0; i<n; i++) {
      A += z[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += xy[i]/(s[i]*s[i]);
      D += z[i]*z[i]/(s[i]*s[i]);
      E += z[i]*xy[i]/(s[i]*s[i]);
      F += xy[i]*xy[i]/(s[i]*s[i]);
    }
    
    m_Az=(E*B-C*A)/(D*B-A*A);
    m_Bz=(D*C-E*A)/(D*B-A*A);
    
    if (m_xyFitFlag==0) {
      m_u0 = m_Az;
      m_x0 = m_Bz;
      m_v0 = m_Axy*m_Az;
      m_y0 = m_Axy*m_Bz+m_Bxy;
    } else {
      m_u0 = m_Axy*m_Az;
      m_x0 = m_Axy*m_Bz+m_Bxy;
      m_v0 = m_Az;
      m_y0 = m_Bz;
    }
    
    m_chisqrZ = 0.;
    for (int i=0; i<n; i++) {
      double xycal = m_Az*z[i] + m_Bz;
      double xcal = m_u0*z[i] + m_x0;
      double ycal = m_v0*z[i] + m_y0;

      double zcal = (xy[i]-m_Bz)/m_Az;      	      
      m_z_track[Layer[i]] = zcal;
      
      m_chisqrZ += (xy[i]-xycal)*(xy[i]-xycal)/(s[i]*s[i]);
      
      // inner(-) or outer(+)       
      if(m_xyFitFlag==0){// x = az + b
	if(xy[i]>=0){
	  m_dist_uv[Layer[i]] = (xy[i]-m_Az*z[i]-m_Bz)/sqrt(1+m_Az*m_Az);
	}else{
	  m_dist_uv[Layer[i]] = -(xy[i]-m_Az*z[i]-m_Bz)/sqrt(1+m_Az*m_Az);
	}
      }else{ // y = az + b
	if(xy[i]<=0){
	  m_dist_uv[Layer[i]] = -(xy[i]-m_Az*z[i]-m_Bz)/sqrt(1+m_Az*m_Az);
	}else{
	  m_dist_uv[Layer[i]] = (xy[i]-m_Az*z[i]-m_Bz)/sqrt(1+m_Az*m_Az);
	}
      }
      
    }
    m_chisqrZ /= n;   
    m_vtx_z = -m_Bz/m_Az;
    


  } else {
    // z = ax + b or z = ay + b
    m_zTrackFlag = 1;    
    
    double A=0., B=0., C=0., D=0., E=0., F=0.;
    for (int i=0; i<n; i++) {
      A += xy[i]/(s[i]*s[i]);
      B += 1./(s[i]*s[i]);
      C += z[i]/(s[i]*s[i]);
      D += xy[i]*xy[i]/(s[i]*s[i]);
      E += z[i]*xy[i]/(s[i]*s[i]);
      F += z[i]*z[i]/(s[i]*s[i]);
    }

    m_Az=(E*B-C*A)/(D*B-A*A);
    m_Bz=(D*C-E*A)/(D*B-A*A);

    if (m_xyFitFlag==0) {
      m_u0 = 1./m_Az;
      m_x0 = -m_Bz/m_Az;
      m_v0 = m_Axy/m_Az;
      m_y0 = -m_Axy*m_Bz/m_Az+m_Bxy;
    } else {
      m_u0 = m_Axy/m_Az;
      m_x0 = -m_Axy*m_Bz/m_Az+m_Bxy;
      m_v0 = 1./m_Az;
      m_y0 = -m_Bz/m_Az;
    }

    m_chisqrZ = 0.;
    for (int i=0; i<n; i++) {
      double zcal = m_Az*xy[i] + m_Bz;
      m_z_track[Layer[i]] = zcal;
      m_chisqrZ += (z[i]-zcal)*(z[i]-zcal)/(s[i]*s[i]);

      // realistic       
      m_dist_uv[Layer[i]] = -(z[i]-m_Az*xy[i]-m_Bz)/sqrt(1+m_Az*m_Az);

    }
    m_chisqrZ /= n;    
    m_vtx_z = m_Bz;
    
  }  
  
  return true;
}


//______________________________________________________________________________
bool
DCLocalTrack::FindLayer( int layer ) const
{
  const std::size_t n = m_hit_array.size();
  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = m_hit_array[i];
    if( !hitp )
      continue;
    if( layer == hitp->GetLayer() )
      return true;
  }
  return false;
}

//______________________________________________________________________________
double
DCLocalTrack::GetWire( int layer ) const
{
  const std::size_t n = m_hit_array.size();
  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = m_hit_array[i];
    if( !hitp )
      continue;
    if( layer == hitp->GetLayer() )
      return hitp->GetWire();
  }
  return math::nan();
}

//______________________________________________________________________________
double
DCLocalTrack::GetDifVXU( void ) const
{
  static const double Cu = cos(  15.*math::Deg2Rad() );
  static const double Cv = cos( -15.*math::Deg2Rad() );
  static const double Cx = cos(   0.*math::Deg2Rad() );

  return
    pow( m_Av/Cv - m_Ax/Cx, 2 ) +
    pow( m_Ax/Cx - m_Au/Cu, 2 ) +
    pow( m_Au/Cu - m_Av/Cv, 2 );
}

//______________________________________________________________________________
double
DCLocalTrack::GetDifVXUSDC34( void ) const
{
  static const double Cu = cos(  30.*math::Deg2Rad() );
  static const double Cv = cos( -30.*math::Deg2Rad() );
  static const double Cx = cos(   0.*math::Deg2Rad() );

  return
    pow( m_Av/Cv - m_Ax/Cx, 2 ) +
    pow( m_Ax/Cx - m_Au/Cu, 2 ) +
    pow( m_Au/Cu - m_Av/Cv, 2 );
}

//______________________________________________________________________________
double
DCLocalTrack::GetTheta( void ) const
{
  double cost = 1./std::sqrt(1.+m_u0*m_u0+m_v0*m_v0);
  return std::acos(cost)*math::Rad2Deg();
}

//______________________________________________________________________________
bool
DCLocalTrack::HasHoneycomb( void ) const
{
  for( std::size_t i=0, n = m_hit_array.size(); i<n; ++i ){
    DCLTrackHit *hitp = m_hit_array[i];
    if( !hitp ) continue;
    if( hitp->IsHoneycomb() ) return true;
  }
  return false;
}

//______________________________________________________________________________
bool
DCLocalTrack::ReCalc( bool applyRecursively )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  std::size_t n = m_hit_array.size();
  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = m_hit_array[i];
    if( hitp ) hitp->ReCalc( applyRecursively );
  }

  bool status = DoFit();
  if( !status ){
    hddaq::cerr << "#W " << func_name << " "
		<< "Recalculation fails" << std::endl;
  }

  return status;
}

//______________________________________________________________________________
void
DCLocalTrack::Print( const std::string& arg, std::ostream& ost ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  PrintHelper helper( 3, std::ios::fixed, ost );

  const int w = 8;
  ost << func_name << " " << arg << std::endl
      << " X0 : " << std::setw(w) << std::left << m_x0
      << " Y0 : " << std::setw(w) << std::left << m_y0
      << " U0 : " << std::setw(w) << std::left << m_u0
      << " V0 : " << std::setw(w) << std::left << m_v0;
  helper.setf( std::ios::scientific );
  ost << " Chisqr : " << std::setw(w) << m_chisqr << std::endl;
  helper.setf( std::ios::fixed );
  const std::size_t n = m_hit_array.size();
  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = m_hit_array[i];
    if( !hitp ) continue;
    int lnum = hitp->GetLayer();
    double zz = hitp->GetZ();
    double s  = hitp->GetLocalHitPos();
    double res = hitp->GetResidual();
    // double aa = hitp->GetTiltAngle()*math::Deg2Rad();
    // double scal=GetX(zz)*cos(aa)+GetY(zz)*sin(aa);
    const std::string& h = hitp->IsHoneycomb() ? "+" : "-";
    ost << "[" << std::setw(2) << i << "]"
	<< " #"  << std::setw(2) << lnum << h
	<< " S " << std::setw(w) << s
	<< " ( " << std::setw(w) << GetX(zz)
	<< ", "  << std::setw(w) << GetY(zz)
	<< ", "  << std::setw(w) << zz
	<< " )"
	<< " " << std::setw(w) << s
	<< " -> " << std::setw(w) << res << std::endl;
	// << " -> " << std::setw(w) << s-scal << std::endl;
  }
  ost << std::endl;
}

//______________________________________________________________________________
void
DCLocalTrack::PrintVXU( const std::string& arg ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  PrintHelper helper( 3, std::ios::fixed );

  const int w = 10;
  hddaq::cout << func_name << " " << arg << std::endl
	      << "(Local X = A*z + B) "
	      << " A : " << std::setw(w) << m_a
	      << " B : " << std::setw(w) << m_b;
  helper.setf( std::ios::scientific );
  hddaq::cout << " Chisqr : " << std::setw(w) << m_chisqr << std::endl;
  helper.setf( std::ios::fixed );

  const std::size_t n = m_hit_array.size();
  for( std::size_t i=0; i<n; ++i ){
    const DCLTrackHit * const hitp = m_hit_array[i];
    if( !hitp ) continue;
    int    lnum = hitp->GetLayer();
    double zz   = hitp->GetZ();
    double s    = hitp->GetLocalHitPos();
    double res  = hitp->GetResidualVXU();
    const std::string& h = hitp->IsHoneycomb() ? "+" : "-";
    hddaq::cout << "[" << std::setw(2) << i << "] "
		<< " #"  << std::setw(2) << lnum << h
		<< " S " << std::setw(w) << s
		<< " ( " << std::setw(w) << (m_a*zz+m_b)
		<< ", "  << std::setw(w) << zz
		<< " )"
		<< " " << std::setw(w) << s
		<< " -> " << std::setw(w) << res << std::endl;
  }
  hddaq::cout << std::endl;
}

// for CFT tracking
bool DCLocalTrack::GetCrossPointR(double r, double *phi, int layer)
{
  double a = m_Axy;
  double b = m_Bxy;
  double x1, y1, x2, y2;

  if (m_xyFitFlag==0) {
    double hanbetsu = a*a*b*b-(1.+a*a)*(b*b-r*r);
    if (hanbetsu<0){
      return false;
    }else {
      x1 = (-a*b+sqrt(hanbetsu))/(1+a*a);
      x2 = (-a*b-sqrt(hanbetsu))/(1+a*a);
      y1 = m_Axy*x1+m_Bxy;
      y2 = m_Axy*x2+m_Bxy;
    }
  }else if (m_xyFitFlag==1) {
    double hanbetsu = a*a*b*b-(1.+a*a)*(b*b-r*r);
    if (hanbetsu<0){
      return false;
    }else {
      y1 = (-a*b+sqrt(hanbetsu))/(1+a*a);
      y2 = (-a*b-sqrt(hanbetsu))/(1+a*a);
      x1 = m_Axy*y1+m_Bxy;
      x2 = m_Axy*y2+m_Bxy;
    }
  }


  
  double phi1 = calcPhi(x1, y1);
  double phi2 = calcPhi(x2, y2);
  if(phi1 < 0){phi1 += 360.;}
  if(phi2 < 0){phi2 += 360.;}

  int nPhiHit = GetNHit();
  int llayer[8];
  for (int i=0; i<8; i++){llayer[i] = -1;}

  bool flagZeroCross = false;
  bool flagZone1 = false;
  bool flagZone4 = false;
  for (int i=0; i<nPhiHit; i++) {
    double phiVal = m_hit_array[i]->GetPositionPhi();
    llayer[i] = m_hit_array[i]->GetLayer();
    if (phiVal>=0 && phiVal<=90){
      flagZone1 = true;
    }else if (phiVal>=270 && phiVal<=360){
      flagZone4 = true;
    }
  }

  if (flagZone1 || flagZone4){
    flagZeroCross = true;
  }

  double phi__[4] = {-100,-100,-100,-100};

  double meanPhi1=0, meanPhi2=0.;// 1:layer0~7, 2:layer8~15
  double n1=0., n2=0.;
  for (int i=0; i<nPhiHit; i++) {
    double phiVal = m_hit_array[i]->GetPositionPhi();
    
    phi__[i]=phiVal;
    //std::cout << i << ":" <<phiVal << std::endl;
    if(llayer[i]<8 && phi){meanPhi1 += phiVal;n1+=1.;}
    else{meanPhi2 += phiVal;n2+=1.;}
  }
  meanPhi1 /= n1;
  meanPhi2 /= n2;
  //std::cout << std::endl;

  double meanPhi;
  if(layer<8){meanPhi=meanPhi1;}
  else{meanPhi=meanPhi2;}

  if( fabs(phi1-meanPhi)<fabs(phi2-meanPhi) ){
    *phi = phi1;
  }else if(fabs(fabs(phi1-(meanPhi))-360)<fabs(phi2-meanPhi)){
    *phi = phi1;
  }else{*phi = phi2;}

  return true;
}

double DCLocalTrack::calcPhi(double x, double y)
{
  if (x>=0 && y>=0)
    return atan(y/x)*math::Rad2Deg();
  else if (x<0 && y>=0)
    return 180. + atan(y/x)*math::Rad2Deg();
  else if (x<0 && y<0)
    return 180. + atan(y/x)*math::Rad2Deg();
  else if (x>=0 && y<0)
    return 360. + atan(y/x)*math::Rad2Deg();

}

double DCLocalTrack::CalculateZpos(double phi, DCLTrackHit *cl)
{
  double z0 = -1000.;
  int layer =cl->GetLayer();
  if(layer>7){layer -= 8;}
  double phi0 =cl->GetPositionPhi();//phi@z=0
  double seg =cl->GetMeanSeg();
  double slope = 0;

  if(layer==0||layer==4){
    slope = 400. /360.;
    z0 =  slope * (-1)*phi0;
  }else if(layer==2||layer==6){
    slope = -400. /360.;
    z0 = 400 + slope * (-1)*phi0;
  }

  if(z0>400)z0-=400.;
  if(z0<0)z0+=400.;

  double a1 = z0+slope*phi;
  double a2 = (z0-400.)+slope*phi;
  if (slope<0)
    a2 = (z0+400.)+slope*phi;

  if(a1>400)a1-=400;
  if(a1<0)a1+=400;
  if(a2>400)a2-=400;
  if(a2<0)a2+=400;
  
  if (a1>=0 && a1<=400){
      if(a1<0)
      std::cout << "a1   layer=" << layer << ", seg=" << seg
		<<", z0=" << z0 << ", slope" << slope
		<< ", phi0=" << phi0 << ", a1=" << a1 << std::endl;  
    return (a1);
  }
  else if (a2>=0 && a2<=400){
    if(a2<0)
      std::cout << "a2  layer=" << layer << ", seg=" << seg
		<< ", z0="<< z0 << ", slope" << slope
		<< ", phi0=" << phi0 << ", a2=" << a2 << std::endl;
    return (a2);
  }
  else{
    return -1000.;
  }
}

bool DCLocalTrack::SetCalculatedValueCFT()
{
  const std::string funcname = "[DCLocalTrack::SetCalculatedValueCFT()]";
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  std::size_t n = m_hit_array.size();

  ThreeVector pos1, pos2;

  // phi on track for not hit layer
  int Nlayer = 8;
  if(n>4){Nlayer=16;}
  for( std::size_t i=0; i<Nlayer; ++i ){
    double phi;
    int lnum=i, seg=0;
    int l_geo = i +301;
    if(i>7){l_geo-=8;}
    double r = gGeom.CalcCFTPositionR(l_geo,seg);
    bool ret = GetCrossPointR(r, &phi, lnum);
    if (!ret) {
      //std::cout << funcname << " return at GetCrossPointR" << std::endl;
      return false;
    }
    m_phi_track[lnum] = phi;    
  }

  for( std::size_t i=0; i<n; ++i ){
    DCLTrackHit *hitp = m_hit_array[i];
    if( hitp ){
      int lnum = hitp->GetLayer();
      double phi;

      int seg=0;
      int l_geo = lnum +301;
      if(lnum>7){l_geo -= 8;}
      //double r = hitp->GetRMax();
      double r = hitp->GetPositionR();
      bool ret = GetCrossPointR(r, &phi, lnum);
      if (!ret) {
	std::cout << funcname << " return at GetCrossPointR" << std::endl;
	return false;
      }
      
      m_dphi[lnum] = m_phi_ini[lnum] - phi;
      m_phi_track[lnum] = phi;

      double xcal = r*cos(phi*math::Deg2Rad());
      double ycal = r*sin(phi*math::Deg2Rad());
      double zcal;
      if (m_zTrackFlag==0) {

	if (m_xyFitFlag==0) {
	  zcal = (xcal- m_Bz)/m_Az;
	} else {
	  zcal = (ycal- m_Bz)/m_Az;
	}

	if (i==0) {
	  pos1 = ThreeVector(m_u0*zcal+m_x0, m_v0*zcal+m_y0, zcal);
	} else if (i==n-1) {
	  pos2 = ThreeVector(m_u0*zcal+m_x0, m_v0*zcal+m_y0, zcal);
	}

      } else if (m_zTrackFlag==1){

	if (m_xyFitFlag==0) {
	  zcal = m_Az *xcal + m_Bz;
	  ycal = m_Axy*xcal + m_Bxy; 
	} else {
	  zcal = m_Az *ycal + m_Bz;
	  xcal = m_Axy*ycal + m_Bxy; 
	}

	if (i==0) {
	  pos1 = ThreeVector(xcal, ycal, zcal);
	} else if (i==n-1) {
	  pos2 = ThreeVector(xcal, ycal, zcal);
	}

      }

    }
  }

  m_Pos0 = pos1;
  m_Dir = pos2-pos1;
  
  double D=(m_Dir.x()*m_Dir.x()+m_Dir.y()*m_Dir.y()+m_Dir.z()*m_Dir.z());
  m_theta = acos(m_Dir.z()/sqrt(D))*180./acos(-1.);

  
  std::size_t nUV = m_hit_arrayUV.size();
  for( std::size_t i=0; i<nUV; ++i ){
    DCLTrackHit *hitp = m_hit_arrayUV[i];
    if( hitp ){
      int lnum = hitp->GetLayer();
      int l_geo = lnum +301;
      if(lnum>7){l_geo -= 8;}
      //double r = hitp->GetRMax();
      double r = hitp->GetPositionR();
      double phi;      
      bool ret = GetCrossPointR(r, &phi, lnum);
      if (!ret) {
	//std::cout << funcname << " return at GetCrossPointR" << std::endl;
	return false;
      }
      m_phi_track[lnum] = phi;

      double z1 = CalculateZpos(phi, hitp);
      double xycal = m_Az*z1 + m_Bz;
      double xcal = m_u0*z1 + m_x0;
      double ycal = m_v0*z1 + m_y0;

      m_dz[lnum] = m_z_ini[lnum] - m_z_track[lnum];

    }
  }
  
  return true;

}
