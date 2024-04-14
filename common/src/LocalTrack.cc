// -*- C++ -*-

#include <map>
#include <iomanip>

#include "LocalTrack.hh"
#include "LineFit.hh"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "BLDCWireMapMan.hh"
#include "DebugCounter.hh"
#include "RootData.hh"

#define DEBUG 0
namespace
{
  const double mm=0.1;
  const double cm=10*mm;
  const std::string& class_name("LocalTrack");
  const BLDCWireMapMan&       gBLDC  = BLDCWireMapMan::GetInstance();
  // const UserParamMan& gUser = UserParamMan::GetInstance();
  // const XTMapMan& gXt = XTMapMan::GetInstance();
  const bool SelectTDC1st  = true;
  const double SpatialResolutionOfBLDC=0.2*mm; // [mm]
  TVector3 DEFV(-999,-999,-999);
}

LocalTrack::LocalTrack() :
  good(true)
{
  m_hit[0].clear();
  m_hit[1].clear();
  A=B=C=D=E=F=-999.;
  GA=GB=GC=GD=GE=GF=-999.;
  xzDof=yzDof=0;
  xzChi=yzChi=-999.;
  //  std::cout<<class_name<<" increased"<<std::endl;
  //  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
LocalTrack::LocalTrack( LocalTrack *tra )
{
  for(int i=0;i<2;i++) {
    for(int j=0;j<tra->nhit(i);j++){
      m_hit[i].push_back((tra->hit(i,j)));
    }
  }
  for(int j=0;j<tra->nclustertimes();j++){
    AddClusterTime(tra->clustertime(j));
  }
  CID=tra->cid();
  A=tra->a();
  B=tra->b();
  C=tra->c();
  D=tra->d();
  E=tra->e();
  F=tra->f();
  GA=tra->ga();
  GB=tra->gb();
  GC=tra->gc();
  GD=tra->gd();
  GE=tra->ge();
  GF=tra->gf();
  GZ=tra->gz();

  xzDof=tra->dofxz();
  yzDof=tra->dofyz();
  xzChi=tra->chi2xz();
  yzChi=tra->chi2yz();
  good=tra->isgood();
}
//______________________________________________________________________________
LocalTrack::LocalTrack(const int &detid) :
  CID(detid),good(true)
{
  //  std::cout<<class_name<<" increased"<<std::endl;
  //  debug::ObjectCounter::increase(class_name);
}

// LocalTrack::LocalTrack(const LocalTrack &right)
// {
//   //  std::cout<<class_name<<" increased"<<std::endl;
//   //  debug::ObjectCounter::increase(class_name);
//}

LocalTrack::~LocalTrack()
{
  //  std::cout<<class_name<<" decreased"<<std::endl;
  //  debug::ObjectCounter::decrease(class_name);
}
void LocalTrack::Clear()
{
  clusterTimes.clear();
  m_hit[0].clear();
  m_hit[1].clear();
}

//______________________________________________________________________________
bool LocalTrack::DeleteHit( const int &xy ,const int &i)
{
  std::vector<TrackHit>::iterator it = m_hit[xy].begin();
  for( int j=0; j<i; j++ ) ++it;
  m_hit[xy].erase( it ) ;
  return true;
}

bool LocalTrack::CompareTrackHit( LocalTrack *tr )
{
  bool check=false;
  for(int xy=0;xy<2;xy++){
    for(int i1=0;i1<this->nhit(xy);i1++){
      for(int i2=0;i2<tr->nhit(xy);i2++){
        TrackHit h1=this->hit(xy,i1);
        TrackHit h2=tr->hit(xy,i2);
	if(h1.cid  !=h2.cid  ) continue;
	if(h1.layer!=h2.layer) continue;
	if(h1.wire !=h2.wire ) continue;
	check= true;
	return true;
      }
    }
  }
  return check;
}

bool LocalTrack::XYLocalPosatZ( const double &z, double &x, double &y)
{
  if( fabs(A) < 1.0e-9 ) return false;
  if( fabs(D) < 1.0e-9 ) return false;

  double tmpB=B,tmpE=E;
  x = -1.*( tmpB*z + C )/A;
  y = -1.*( tmpE*z + F )/D;
  return true;
}

bool LocalTrack::XYPosatZ( const double &z, double &x, double &y )
{
  if( fabs(GA) < 1.0e-9 ) return false;
  if( fabs(GD) < 1.0e-9 ) return false;

  x = -1.*( GB*(z-GZ) + GC )/GA;
  y = -1.*( GE*(z-GZ) + GF )/GD;

  return true;
}

bool LocalTrack::ZYPosatX( const double &x, double &z, double &y )
{
  if( fabs(GB) < 1.0e-9 ) return false;
  if( fabs(GD) < 1.0e-9 ) return false;

  z = -1.*( GA*x +GC ) / GB + GZ;
  y = -1.*( GE*(z-GZ) + GF )/GD;
  return true;
}

bool LocalTrack::ZXPosatY( const double &y, double &z, double &x )
{
  if( fabs(GA) < 1.0e-9 ) return false;
  if( fabs(GE) < 1.0e-9 ) return false;

  z = -1.*( GD*y +GF ) / GE + GZ;
  x = -1.*( GB*(z-GZ) + GC )/GA;
  return true;
}

TVector3 LocalTrack::GetLocalPosatZ(const double &z)
{
  double tmpx,tmpy;
  XYLocalPosatZ(z,tmpx,tmpy);
  return TVector3(tmpx,tmpy,z);
}
TVector3 LocalTrack::GetPosatZ(const double &z)
{
  double tmpx,tmpy;
  XYPosatZ(z,tmpx,tmpy);
  return TVector3(tmpx,tmpy,z);
}
TVector3 LocalTrack::GetPosatY(const double &y)
{
  double tmpx,tmpz;
  ZXPosatY(y,tmpz,tmpx);
  return TVector3(tmpx,y,tmpz);
}

TVector3 LocalTrack::GetPosatX(const double &x)
{
  double tmpy,tmpz;
  ZXPosatY(x,tmpz,tmpy);
  return TVector3(x,tmpy,tmpz);
}
TVector3 LocalTrack::GetMomDir()
{
  TVector3 tmpdir(-GB/GA,-GE/GD,1);
  return tmpdir;
}
TVector3 LocalTrack::GetLocalMomDir()
{
  TVector3 tmpdir(-B/A,-E/D,1);
  return tmpdir;
}


static const unsigned int LRMASK=0x01;
inline static unsigned int LR( int hid, unsigned int key )
{
  key = key >> hid;
  key = key & LRMASK;
  return key;
}
bool LocalTrack::DoFit()
{
  if(CID==DetIdFDC){
    //    return 0;
    return LinearFit()&&dRLineFit() && ConvLocalToGlobal();
  }
  bool tmp1=LeastSquareFit(0);
  bool tmp2=LeastSquareFit(1);
  bool tmp3=ConvLocalToGlobal();
  bool tmp=tmp1&&tmp2&&tmp3;
#if DEBUG
  std::cout<<__PRETTY_FUNCTION__<<"  "<<tmp1<<"  "<<tmp2<<"  "<<tmp3<<std::endl;
#endif
  return tmp;
}
bool LocalTrack::LeastSquareFit( const int &xy, TString option )
{
#if DEBUG
  std::cout << "!!! LocalTrack::LeastSquareFit()  "<<option << std::endl;
#endif
  int np = nhit(xy);
  const int MAXHIT=32;
  if(np>MAXHIT){
    return false;
  }
  if( np < 3 ){
    return false;
  }
  double sigma=SpatialResolutionOfBLDC;
  double minchi = 1.0e+9;
  double canda=0., candb=0., candc=0.;
  unsigned int candkey=0x0;
  unsigned int key = 0x0;
  int nconb = (int)pow(2,np);
  if(option.Contains("mwpc")) nconb=1;

  double x[MAXHIT],y[MAXHIT],w[MAXHIT];
  double tmp,tmp2;
  for( int ic=0; ic<nconb; ic++ ){
    int j=0;
    for( int i=0; i<np;i++){
      x[j] = m_hit[xy][i].wpos.Z();
      y[j] = m_hit[xy][i].wpos.X();
      double dl = m_hit[xy][i].dl;
      if(!option.Contains("mwpc")){
	w[j]=sigma;
	unsigned int lr = LR(j,key);
	if( lr==0 ) y[j] += dl;
	else y[j] -= dl;
      }else{
	w[j]=sigma;
      }
      // std::cout<<i<<"  "<<j<<"  "<<x[j]<<"  "<<y[j]<<"  "<<dl<<"  "<<w[j]<<std::endl;
      j++;
    }
    double par[2],chi2;
    math::LineFit(x,y,w,np,par,chi2);
    if( (chi2 < minchi) || (chi2==minchi && candb*candb<par[1]*par[1]) ){
      canda = 1.;
      candb = -1.*par[1];
      candc = -1.*par[0];
      candkey=key;
      minchi=chi2;
    }
    key++;
  }
  if(minchi==1.0e+9){
#if DEBUG
    std::cout<<"LocalTrack::LeastSquareFit():  no track candidate "<<std::endl;
#endif
    return false;
  }
  if(xy)
    SetDEF(canda,candb,candc);
  else
    SetABC(canda,candb,candc);
  SetDof(xy,np-2);
  SetChisqr(xy,minchi);
  if(option.Contains("mwpc")) return true;
  int j=0;
  //  std::cout<<"-----------------------------__"<<std::endl;
  for( int i=0; i<np; i++ ){
    double dl = m_hit[xy][i].dl;
    double zpos= m_hit[xy][i].wpos.Z();
    double wpos = m_hit[xy][i].wpos.X();
    double tmpx,tmpy;
    unsigned int lr = LR(i,candkey);
    double xpos=wpos;
    if( lr==0 ) xpos+=dl;
    else        xpos-=dl;
    TVector3 tmppos=GetLocalPosatZ(zpos);
    if(xy==0){
      m_hit[xy][i].hitpos=TVector3(xpos,0,zpos);
      tmppos.SetY(0);
      m_hit[xy][i].trackpos=tmppos;
      m_hit[xy][i].residual=xpos-tmppos.X();
    }else{
      m_hit[xy][i].hitpos=TVector3(0,xpos,zpos);
      tmppos.SetX(0);
      m_hit[xy][i].trackpos=tmppos;
      m_hit[xy][i].residual=xpos-tmppos.Y();
    }
    //    std::cout<<xy<<"  "<<i<<" / "<<np<<"  z:"<<zpos<<"  w:"<<wpos<<" lr:"<<lr<<" dl: "<<dl<<" trapos: "<<tmppos[xy]<<" hitpos:"<<m_hit[xy][i].hitpos<<"  resid:"<<m_hit[xy][i].residual<<std::endl;
  }
  return true;
}

bool LocalTrack::LinearFit( TString option )
{
#if 0
  std::cout << "!!! LinearTrack::LinearFit()" << std::endl;
#endif
  int MaxNumOfHitsInTrack=64;
  int np = nhit(0);
  if( np > MaxNumOfHitsInTrack ){
    return false;
  }
  if( np < 2 ) return false;
  double x[MaxNumOfHitsInTrack], y[MaxNumOfHitsInTrack], theta[MaxNumOfHitsInTrack], w[MaxNumOfHitsInTrack];
  int np_org = np;
  int nlr=0;
  int ihitlr[MaxNumOfHitsInTrack];
  double minchi = 1.0e+9;
  double canda=0., candb=0., candc=0.,candd=0.;
  unsigned int candkey=0x0;
  unsigned int key = 0x0;
  nlr=np;

  double sigma=SpatialResolutionOfBLDC;
  int nconb = (int)pow(2,nlr);
  for( int ic=0; ic<nconb; ic++ ){
    for( int i=0; i<np;i++){
      x[i] = m_hit[0][i].wpos.Z();
      y[i] = m_hit[0][i].wpos.X();
      w[i]=sigma;
      double dl = m_hit[0][i].dl;
      unsigned int lr = LR(i,key);
      if( lr==0 ) y[i] += dl;
      else y[i] -= dl;
      theta[i]=m_hit[0][i].rotation*TMath::DegToRad();

    }
    double right[4];
    double c2=0.,s2=0.,cs=0.,zc2=0.,zs2=0.,zcs=0.,z2c2=0.,z2s2=0.,z2cs=0;
    for(int i=0;i<4;i++){
      right[i]=0.;
    }
    for( int i=0; i<np; i++ ){
      //std::cout<<"x: "<<x[i]<<" ,y: "<<y[i]<<" ,w: "<<w[i]<<" ,theta: "<<theta[i]<<std::endl;
      c2 += pow(TMath::Cos(theta[i]),2);
      s2 += pow(TMath::Sin(theta[i]),2);
      cs += TMath::Cos(theta[i])*TMath::Sin(theta[i]);
      zc2 += x[i]*pow(TMath::Cos(theta[i]),2);
      zs2 += x[i]*pow(TMath::Sin(theta[i]),2);
      zcs += x[i]*TMath::Cos(theta[i])*TMath::Sin(theta[i]);
      z2c2 += x[i]*x[i]*pow(TMath::Cos(theta[i]),2);
      z2s2 += x[i]*x[i]*pow(TMath::Sin(theta[i]),2);
      z2cs += x[i]*x[i]*TMath::Cos(theta[i])*TMath::Sin(theta[i]);
      right[0]+=y[i]*TMath::Cos(theta[i]);
      right[1]+=y[i]*x[i]*TMath::Cos(theta[i]);
      right[2]+=y[i]*TMath::Sin(theta[i]);
      right[3]+=y[i]*x[i]*TMath::Sin(theta[i]);
    }
    TMatrixD mat(4,4);
    mat[0][0]= c2; mat[0][1]=zc2; mat[0][2]=cs; mat[0][3]=zcs;
    mat[1][0]= zc2; mat[1][1]=z2c2; mat[1][2]=zcs; mat[1][3]=z2cs;
    mat[2][0]= cs; mat[2][1]=zcs; mat[2][2]=s2; mat[2][3]=zs2;
    mat[3][0]= zcs; mat[3][1]=z2cs; mat[3][2]=zs2; mat[3][3]=z2s2;

    TMatrixD aa;
    aa.Use(4,1,right);
    mat.InvertFast();
    TMatrixD par(4,1);
    par.Mult(mat,aa);

    double chi = 0;
    double tempx,tempy,trackpos;
    for( int i=0; i<np; i++ ){
      tempx = par[0][0]+par[1][0]*x[i];
      tempy = par[2][0]+par[3][0]*x[i];
      trackpos = tempx*TMath::Cos(theta[i])+tempy*TMath::Sin(theta[i]);
      chi += pow((y[i]-trackpos)/w[i],2);
    }
    if(np<5) chi=-1;
    else chi /= (double)(np-4);
    if( (chi < minchi) ){
      canda = par[0][0];
      candb = par[1][0];
      candc = par[2][0];
      candd = par[3][0];
      candkey = key;
      minchi=chi;
    }
    key++;
  }
  SetABC(1.,-candb,-canda);
  SetDEF(1.,-candd,-candc);
  SetDof(0,np-4);
  SetChisqr(0,minchi);
  return true;
}
bool LocalTrack::dRLineFit( TString option )
{
  int MaxNumOfHitsInTrack=20;
  if( nhit(0) > MaxNumOfHitsInTrack )
    return false;

  TVector3 wirepos[20];
  TVector3 wiredir[20];
  double weight[MaxNumOfHitsInTrack];
  double drift[MaxNumOfHitsInTrack];

  int allhit=0;
  double sigma=SpatialResolutionOfBLDC;
  for(int n=0;n<(int)nhit(0);n++)
    {
      wirepos[n]=m_hit[0][n].wpos;
      wiredir[n]=m_hit[0][n].wdir;
      weight[n]=sigma;
      drift[n]=m_hit[0][n].dl;
#if DEBUG
      std::cout<<"---- "<<n<<"  / "<<nhit(0)<<std::endl;
      wirepos[n].Print();
      wiredir[n].Print();
      std::cout<<drift[n]<<std::endl;
#endif
      allhit++;
    }

  double tmppar[6]={x(),y(),0,
		    dx(),dy(),1};
  LineFit *fit=new LineFit(tmppar,wirepos,wiredir,weight,drift,allhit);
  fit->GetParameters(tmppar);
  xzChi=fit->chisquare();
  xzDof=fit->dof();
  delete fit;
  TVector3 pos(tmppar[0],tmppar[1],tmppar[2]);
  TVector3 dir(tmppar[3],tmppar[4],tmppar[5]);
  dir=1./dir.Z()*dir;
  pos=pos-dir*pos.Z();
  A=1;
  B=-dir.X();
  C=-pos.X();
  D=1;
  E=-dir.Y();
  F=-pos.Y();

#if DEBUG
  std::cout<<"------calc residuals for dRLinearFit ----"<<std::endl;
#endif
  int xy=0;
  TVector3 tmppos=GetLocalPosatZ(0);
  TVector3 tmpdir=GetLocalMomDir();
  for( int i=0; i<nhit(0); i++ ){
    double dl = m_hit[xy][i].dl;
    double dist;
    TVector3 fittmp1,fittmp2;
    LineToLine(m_hit[xy][i].wpos,m_hit[xy][i].wdir,tmppos,tmpdir,0,dist,fittmp1,fittmp2);
    m_hit[xy][i].residual=dist-dl;
#if DEBUG
    std::cout<<i<<"  "<<dist<<"  "<<dl<<"  "<<dist-dl<<std::endl;
#endif
  }
  return true;
}

bool LocalTrack::ConvLocalToGlobal()
{
  TVector3 gpos;
  TVector3 grot;
  if(!gBLDC.GetGParam(m_hit[0][0].cid,gpos,grot)) return false;
  // std::cout<<"================="<<std::endl;
  // std::cout<<"Global position: "<<m_hit[0][0].cid; gpos.Print();
  // std::cout<<"Global rotation: "<<m_hit[0][0].cid; grot.Print();
  TVector3 pos(-C/A,-F/D,0);
  pos.RotateZ(grot.Z()*TMath::DegToRad());
  pos.RotateX(grot.X()*TMath::DegToRad());
  pos.RotateY(grot.Y()*TMath::DegToRad());
  TVector3 dir(-B/A,-E/D,1);
  dir.RotateZ(grot.Z()*TMath::DegToRad());
  dir.RotateX(grot.X()*TMath::DegToRad());
  dir.RotateY(grot.Y()*TMath::DegToRad());
  GA = 1.;
  GB = -dir.X()/dir.Z();
  GC = -(pos.X()-dir.X()/dir.Z()*pos.Z()+gpos.X());
  GD = 1.;
  GE = -dir.Y()/dir.Z();
  GF = -(pos.Y()-dir.Y()/dir.Z()*pos.Z()+gpos.Y());
  GZ = gpos.Z();
  return true;
}

void LocalTrack::Print(){
  std::cout<< "DC_ID:" << this->cid() <<std::endl;
  for(int xy =0;xy<2;xy++){
   for(int i=0;i<this->nhit(xy);i++){
     TrackHit hit=this->hit(xy,i);
     std::cout<<"layer: "<<hit.layer
	      <<", wire: "<<hit.wire
	      <<", hitpos: "<<hit.hitpos.X()
	      <<", "<<hit.hitpos.Y()
	      <<", "<<hit.hitpos.Z()
	      <<", wpos: "<<hit.wpos.X()
	      <<", "<<hit.wpos.Y()
	      <<", "<<hit.wpos.Z()
	      <<std::endl;
   }
  }
  std::cout<< "x: Chi2,Dof,a,b,c :" << this->chi2xz() << "\t" << this->dofxz()
   	   << "\t" << this->a() << "\t" << this->b() << "\t" << this->c() << std::endl;
  std::cout<< "           ga,b,c :" << this->chi2xz() << "\t" << this->dofxz()
  	   << "\t" << this->ga() << "\t" << this->gb() << "\t" << this->gc() << std::endl;
  std::cout<< "y: Chi2,Dof,a,b,c :" << this->chi2yz() << "\t" << this->dofyz()
  	   << "\t" << this->d() << "\t" << this->e() << "\t" << this->f() << std::endl;
  std::cout<< "           ga,b,c :" << this->chi2yz() << "\t" << this->dofyz()
  	   << "\t" << this->gd() << "\t" << this->ge() << "\t" << this->gf() << std::endl;
  std::cout<<"------------------------------------------------------"<<std::endl;
}
