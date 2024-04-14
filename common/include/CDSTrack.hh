#ifndef CDSTrack_h
#define CDSTrack_h 1

#include <vector>
#include <iostream>
#include <iomanip>

#include "TVector3.h"
#include "DCHit.hh"
#include "Hodo2Hit.hh"
#include "MathTools.hh"
#include "Event.hh"

class DCCluster;

#define TMPPARAM 0

class CDSTrack
{
public:
  struct TrackHit
  {
    int layer;
    int wire;    
    int nth;    
    TVector3 wpos;
    TVector3 wdir;
    double dt;
    double cdt;
    double dl;
    double cdl;
    double tot;
    TVector3 hitpos;
    TVector3 trackpos;
    double resolution;
    double residual;
    TrackHit( const DCHit* hit, int nh )
    {
      nth=nh;
      wpos=hit->GetWirePosition();
      wdir=hit->GetWireDirection();
      layer=hit->GetLayer();
      wire=hit->GetWire();
      dt=hit->GetDriftTime(nh);
      cdt=dt;
      dl=hit->GetDriftLength(nh);
      cdl=dl;
      tot=hit->GetTOT(nh);
      resolution=-1;
    }
    TrackHit( const TVector3 vtx, int rz, double resol )
    {
      nth=rz;
      wpos=vtx;
      if(rz==0){
	wdir=TVector3(0,1,1);
      }else{
	wdir=TVector3(1,0,1);
      }
      layer=0;
      wire=0;
      dt=0;
      cdt=0;
      dl=0;
      cdl=0;
      tot=0;
      resolution=resol;
    }
    inline void Print( void ) const
    {
      std::cout << std::setw(10) << "CDCTrackHit"
       		<< std::setw(5)  << layer
       		<< std::setw(5)  << wire
		<< std::setw(10) << wpos.X()
		<< std::setw(10) << wpos.Y()
		<< std::setw(10) << dt
		<< std::setw(10) << cdt
		<< std::setw(10) << dl
		<< std::setw(10) << cdl
		<< std::setw(10) << tot << std::endl;
    }
  };

 private:
  std::vector<TrackHit> m_hit[15];
  double ChiSquare;
  int Dof;
  double Param[5];
  // param[0]: y0 (x=0)
  // param[1]: dy/dx
  // param[0]: x0
  // param[1]: y0
  // param[2]: phi(dy/dx)  (-pi/2, 3*pi/2)
  // param[3]: z0
  // param[4]: theta (dr/dz)
  // for line track param[0:2]: position closest to R=0, param[3]:theta (dr/dz), param[4]:phi (dy/dx)
  int fittinglevel;
  int trackID;
  bool goodflag;
  bool lineflag;

  typedef std::vector<double> parContainer;
  std::map<int,parContainer> parCont;
  std::map<int,TVector3> vtxContainer;
  std::map<int,double> tofContainer;
  std::map<int,double> flContainer;

 public:
  CDSTrack();
  CDSTrack( CDSTrack *track);
  virtual ~CDSTrack() {};
  void Clear();  

 private:
  void CheckCharge();

 public:
  //  void CalcELoss(double mass=-1);
  bool AddCluster(DCCluster *cl);
  int AddHit(const DCHit *hit, int nh){
    TrackHit newhit(hit,nh);    
    m_hit[newhit.layer].push_back(newhit);
    return nhits(newhit.layer);
  }
  int AddVertex(TVector3 vtx, double resol){
    TrackHit newhit(vtx,0,resol); // r    
    m_hit[newhit.layer].push_back(newhit);
    TrackHit newhit2(vtx,1,resol); // r    
    m_hit[newhit2.layer].push_back(newhit2);
    return nhits(newhit.layer);
  }

  int FittingLevel() const { return fittinglevel; }
  int trackid()      const { return trackID; }
  int nhits() const;
  double totsum()  const;
  double totave()  const { return totsum()/nhits(); }
  double totave2() const { return totsum()/nhits()/sqrt(1+Param[4]*Param[4]); }
  int nhits( const int &layer ) const { 
    if(layer<15&&layer>=0) return m_hit[layer].size(); 
    return 0;
  }
  TrackHit *hit( const int &layer, const int &i ) { return &m_hit[layer][i]; }
  int wire( const int &layer, const int &i ) const { return m_hit[layer][i].wire; }
  bool DeleteHit( const int &layer,const int &i );

  CDCTrackHits get_trackhits();

  int naxialhits();
  int nstereohits();
  int nhitlayers();
  int naxialhitlayers();
  int nstereohitlayers();
  int nsuperlayers();
  int naxialsuperlayers();
  int nstereosuperlayers();

  bool FirstHelixFitting();
  bool HelixFitting();
  bool CircleFitting();

  bool FirstLineFitting();
  bool LineFitting();
  bool LineAxialFitting();
  //calc? sethitpos? calcchi2? vertex? get position ?

  bool Calc();
  bool IsGood()  { return goodflag; }
  bool IsLine()  { return lineflag; }
  void SetGoodFlag( const bool &flag ) { goodflag = flag; }
  void SetLineFlag( const bool &flag=true ) { lineflag = flag; }
  void SetTrackID(  const int &id ) { trackID = id; }
  void SetHitPos();

  TVector3 GetMomentumVector(const TVector3 &pos);
  TVector3 GetMomentumVector(const TVector3 &pos,const double *aparam);
  TVector3 GetPositionatR(const double &r);

  //  bool IsLine() {return Param[2]==0 ? true : false; }
  double param(const int &i) { return Param[i]; }
  double CircleR() const { return fabs(1./Param[2]); }
  double CircleX() const { return (Param[0]+1./Param[2])*cos(Param[1]); }
  double CircleY() const { return (Param[0]+1./Param[2])*sin(Param[1]); }
  int  nParamSets() const { return (int)parCont.size(); }

  void AddParameters( const int &id, const double *par,const TVector3 &vtx, const double &tof,const double &totl);
  void GetParameters( double *aparam);
  bool GetGParameters(double *aparam);
  bool GetParameters( const int &id, double *aparam, TVector3 &vtx, double &tof, double &length);
  bool GetParameters( const int &id, double *aparam, TVector3 &vtx);
  bool GetNthParameters( const int &n, int &id, double *aparam, TVector3 &vtx);

  void SetParameters(const double *aparam);
  void SetChi2( const double &chi ) { ChiSquare = chi; }
  void SetDof( const double &adof ) { Dof = adof; }

  void CalcChi2();
  double chi2()		const { return ChiSquare; }
  double dof()		const { return Dof; }
  double pt( )		const { return pt(Param); }
  double mom( )		const { return mom(Param); }
  int    charge()	const { return pt()>0 ? 1 : -1; }               
  double pt(  const double *par )	const;
  double mom( const double *par )	const;
  int pid_tot();

  TVector3 GetPos0() const{ return TVector3(Param[0],Param[1],Param[3]); }
  TVector3 GetDir0() const{ TVector3 tmp; tmp.SetMagThetaPhi(1,Param[4],Param[2]); return tmp; }
  bool GetMomentum(const TVector3 &pos, const double &mass, TVector3 &p,double &tof,double &fl,
		   bool ELOSS=false, bool GLOBAL=true); 
  bool GetMomentum(const TVector3 &pos, const double &mass, TVector3 &p,
		   bool ELOSS=false, bool GLOBAL=true); 
  bool GetVertex(const TVector3 &pos, const TVector3 &dir, TVector3 &lpos, TVector3 &hpos);
  bool CalcVertexTimeLength(const TVector3 &pos,const TVector3 &dir,const double &mass,
			    TVector3 &lpos, TVector3 &hpos,double &time, double &length, 
			    bool ADDPAR=false);
  bool CalcEnergyLoss(const TVector3 &pos, const double &mass, double *par1, TVector3& tmppos);
  bool Retiming(double cdhtime, double beta, bool SLEW=false);
  bool RemoveBadHits(double threshold_chi2=5,int threshold_nhits=13);

  void Print();
};

#endif
