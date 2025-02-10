// -*- C++ -*-

#ifndef LOCAL_TRACK_HH
#define LOCAL_TRACK_HH

#include <iomanip>
#include <iostream>
#include <vector>
#include <map>

#include "DCHit.hh"
#include "DetectorID.hh"
#include "MathTools.hh"
#include "TVector3.h"

struct TrackHit
{
  int cid;
  int layer;
  int wire;

  TVector3 wpos;
  TVector3 wdir;
  double dt;
  double dl;
  double tilt;
  double rotation;
  TVector3 hitpos;
  TVector3 trackpos;
  double residual;

  TrackHit( DCHit* hit, int nh )
  {
    cid=hit->GetDetId();
    wpos=hit->GetWirePosition();
    wdir=hit->GetWireDirection();
    layer=hit->GetLayer();
    wire=hit->GetWire();
    dt=hit->GetDriftTime(nh);
    dl=hit->GetDriftLength(nh);
    tilt=hit->GetTiltAngle();
    rotation=hit->GetRotationAngle();
  }
  inline void Print( void ) const
  {
    std::cout << std::setw(10)  << "TrackHit"
	      << std::setw(5)  << cid
	      << std::setw(5)  << layer
	      << std::setw(5)  << wire
	      << std::setw(10)  << wpos.Z()
	      << std::setw(10)  << wpos.X()
	      << std::setw(10) << dt << std::endl;
  }
};

class LocalTrack
{
public:

private:
  std::vector<TrackHit> m_hit[2];

public:
  LocalTrack();
  //  LocalTrack(const LocalTrack &right);
  LocalTrack(LocalTrack* track);
  LocalTrack(const int &detid);
  ~LocalTrack();

  int AddHit(DCHit *hit, int nh, int xy){
    TrackHit newhit(hit,nh);
    if(CID==DetIdBLC1&&newhit.cid==DetIdBLC1b&&newhit.layer<8) newhit.layer+=8;
    if(CID==DetIdBLC2&&newhit.cid==DetIdBLC2b&&newhit.layer<8) newhit.layer+=8;
    m_hit[xy].push_back(newhit);
    //    CID=hit->GetDetId();
    //    newhit.Print();
    return nhit(xy);
  }
  int AddHit(TrackHit hit,int xy){
    if(CID==DetIdBLC1&&hit.cid==DetIdBLC1b&&hit.layer<8) hit.layer+=8;
    if(CID==DetIdBLC2&&hit.cid==DetIdBLC2b&&hit.layer<8) hit.layer+=8;
    m_hit[xy].push_back(hit);
    return nhit(xy);
  }
  TrackHit hit( const int &xy, const int &i) { return m_hit[xy][i]; }
  int nhit( const int &xy )  { return m_hit[xy].size(); }
  double resid( const int &xy, int i )  { return m_hit[xy][i].residual; }
  int layer( const int &xy, int i )  { return m_hit[xy][i].layer; }
  bool DeleteHit( const int &layer,const int &i );
  bool CompareTrackHit(LocalTrack *tr);

 private:
  std::vector<double> clusterTimes;

 public:
  void AddClusterTime(const double& t) { clusterTimes.push_back(t); }
  int    nclustertimes() const { return clusterTimes.size(); }
  double clustertime(const int &i) const { return clusterTimes.at(i); }
  double GetTrackTime()    const { return math::mean(clusterTimes); }
  double GetTrackTimeRMS() const { return math::max(clusterTimes)-math::min(clusterTimes); }
  bool CheckRange(const double &ll,const double &ul){
    return (GetTrackTime()>ll&&GetTrackTime()<ul);
  }

 private:
  int CID;
  double A, B, C;		/* a track in XZ plane : Ax + Bz + C = 0 */
  double D, E, F;		/* a track in YZ plane : Dy + Ez + F = 0 */
  double GA, GB, GC;		/* Parameters for a global track */
  double GD, GE, GF;		/* Parameters for a global track */
  double GZ;

  int xzDof, yzDof;
  double xzChi, yzChi;
  bool good;

 public:
  int cid()  const { return CID; }
  double a() const { return A; }
  double b() const { return B; }
  double c() const { return C; }
  double d() const { return D; }
  double e() const { return E; }
  double f() const { return F; }
  double x() const  { return -C/A; }
  double dx() const { return -B/A; }
  double y() const  { return -F/D; }
  double dy() const { return -E/D; }
  double ga() const { return GA; }
  double gb() const { return GB; }
  double gc() const { return GC; }
  double gd() const { return GD; }
  double ge() const { return GE; }
  double gf() const { return GF; }
  double gx() const  { return -GC/GA; }
  double gdx() const { return -GB/GA; }
  double gy() const  { return -GF/GD; }
  double gdy() const { return -GE/GD; }
  double gz() const { return GZ; }
  void labc( double &a, double &b, double &c ) { a=A; b=B; c=C; }
  void ldef( double &d, double &e, double &f ) { d=D; e=E; f=F; }
  void abc( double &a, double &b, double &c ) const { a=A; b=B; c=C; }
  void def( double &d, double &e, double &f ) const { d=D; e=E; f=F; }
  void gabc( double &a, double &b, double &c ) const { a=GA; b=GB; c=GC; }
  void gdef( double &d, double &e, double &f ) const { d=GD; e=GE; f=GF; }
  double chi2xz() const { return xzChi; }
  double chi2yz() const { return yzChi; }
  double chi2all() const { return (xzChi*xzDof + yzChi*yzDof)/(double)(xzDof+yzDof); }
  int dofxz() const { return xzDof; }
  int dofyz() const { return yzDof; }

  bool isgood() const { return good; }
  void SetGood() { good=true; }
  void SetBad()  { good=false; }

  void SetABC( const double &a, const double &b, const double &c ){ A=a; B=b; C=c; }
  void SetDEF( const double &d, const double &e, const double &f ){ D=d; E=e; F=f; }
  void SetChisqr( const int &xy, const double &val ) {
    if(xy==0) xzChi=val;
    else if(xy==1) yzChi=val;
  }
  void SetDof( const int &xy, const int &val ) {
    if(xy==0) xzDof=val;
    else if(xy==1) yzDof=val;
  }
  double chi2(const int &xy) const { return xy ? yzChi : xzChi; }
  int dof(const int &xy) const { return xy ? yzDof : xzDof ; }

 public:
  bool XYLocalPosatZ( const double &z, double &x, double &y);
  bool XYPosatZ( const double &z, double &x, double &y );
  bool ZXPosatY( const double &y, double &z, double &x );
  bool ZYPosatX( const double &x, double &z, double &y );
  TVector3 GetLocalPosatZ(const double &z);
  TVector3 GetPosatZ(const double &z);
  TVector3 GetPosatX(const double &z);
  TVector3 GetPosatY(const double &z);
  TVector3 GetMomDir();
  TVector3 GetLocalMomDir();
  bool DoFit();
  bool LeastSquareFit( const int &xy, TString option="" );
  bool LinearFit( TString option="" );
  bool dRLineFit( TString option="" );

public:
  // void Calculate();
  void Clear();
  void Print();
  void LineToLine( const TVector3 &x1, const TVector3 &a1,
		   const TVector3 &x2, const TVector3 &a2,
		   const double &dl,
		   double &dist,
		   TVector3 &xest, TVector3 &next );
  bool ConvLocalToGlobal();

};

inline void LocalTrack::LineToLine( const TVector3 &x1, const TVector3 &a1,
				    const TVector3 &x2, const TVector3 &a2,
				    const double &dl,
                                   double &dist,
                                   TVector3 &xest, TVector3 &next )
{
  TVector3 x = x2-x1;
  double a =  a1.Dot(a1);
  double b = -a1.Dot(a2);
  double c =  a2.Dot(a1);
  double d = -a2.Dot(a2);
  double A1 = a1.Dot(x);
  double A2 = a2.Dot(x);

  double D = a*d-b*c;

  TVector3 x2p;
  if( fabs(D)<0.00000000000001 ){
    dist = sqrt(x.Mag2()-A1*A1);
  }
  else{
    double s = (a*A2-c*A1)/D;
    double t = (d*A1-b*A2)/D;
    xest  = x1 + t*a1;
    x2p   = x2 + s*a2;
    dist = (xest-x2p).Mag();
    next = x2p+(xest-x2p)*(dl/dist);
  }
}

#endif
