#include "CDSTrack.hh"
#include "MathTools.hh"
#include "ELossTools.hh"
#include "GeomTools.hh"
#include "LineFit.hh"
#include "HelixFit.hh"
#include "DCCluster.hh"
#include "UserParamMan.hh"
#include "CDCWireMapMan.hh"
#include "HodoPHCMan.hh"
#include "XTMapMan.hh"
#include "DetectorID.hh"
#include "DatabasePDG.hh"
#define DEBUG 0
//______________________________________________________________________________
namespace
{
  const double mm=0.1;
  const double cm=10*mm;
  const double m=1e3*mm;
  const double ResolutionOfCDC=0.18*mm;

  const double cdh_rin=54.4*cm;

  void print_TV3(TVector3 vec){
    std::cout<<"( "
	     <<std::setw(8)<<vec.X()<<","
	     <<std::setw(8)<<vec.Y()<<","
	     <<std::setw(8)<<vec.Z()<<")";
  }
  const int super[7][3]={{0,1,2},
			 {3,4,-1},
			 {5,6,-1},
			 {7,8,-1},
			 {9,10,-1},
			 {11,12,-1},
			 {13,14,-1}
  };
  const int axsuper[3][3]={{0,1,2},
			   {7,8,-1},
			   {13,14,-1}
  };
  const int stsuper[4][3]={
    {3,4,-1},
    {5,6,-1},
    {9,10,-1},
    {11,12,-1}
  };
  const int axlayer[7]={0,1,2,7,8,13,14};
  const int stereo[8]={3,4,5,6,9,10,11,12};

  const UserParamMan& gUser = UserParamMan::GetInstance();
  const CDCWireMapMan& gCDC = CDCWireMapMan::GetInstance();
  const XTMapMan& gXT       = XTMapMan::GetInstance();
  const HodoPHCMan& gSlew   = HodoPHCMan::GetInstance();
  double magneticfield;
  double tot_pid_param1;
  double tot_pid_param2;
}
//______________________________________________________________________________
CDSTrack::CDSTrack()
{  
  Clear();
  magneticfield=gUser.GetParameter("CDSFIELD");
  tot_pid_param1=gUser.GetParameter("CDCTOTPID",0);
  tot_pid_param2=gUser.GetParameter("CDCTOTPID",1);
}
//______________________________________________________________________________
CDSTrack::CDSTrack( CDSTrack *tra )
{
  for(int i=0;i<15;i++) {
    for(int j=0;j<tra->nhits(i);j++){
      m_hit[i].push_back(*(tra->hit(i,j)));
    }
  }
  fittinglevel=tra->FittingLevel();
  lineflag=tra->IsLine();
  for(int i=0;i<5;i++) Param[i]=tra->param(i);
}
//______________________________________________________________________________
void CDSTrack::Clear()
{
  for(int i=0;i<15;i++) { m_hit[i].clear(); }
  lineflag=false;
  goodflag=true;
  fittinglevel=0;
  for(int i=0;i<5;i++) Param[i]=0;
}
//______________________________________________________________________________
void CDSTrack::GetParameters(double *aparam)
{
  for(int i=0;i<5;i++) aparam[i]=Param[i];
}
//______________________________________________________________________________
bool CDSTrack::GetGParameters(double *aparam)
{
  TVector3 tmp;
  double tmp1,tmp2;
  return GetParameters(DetIdCDC,aparam,tmp,tmp1,tmp2);
}
//______________________________________________________________________________
bool CDSTrack::GetParameters(const int &id, double *aparam, TVector3 &vtx)
{
  double tmp1,tmp2;
  return GetParameters(id,aparam,vtx,tmp1,tmp2);
}
//______________________________________________________________________________
bool CDSTrack::GetParameters(const int &id, double *aparam, TVector3 &vtx, double &tof, double &length)
{
  vtx.SetXYZ(999,999,999);
  if(parCont.find(id)!=parCont.end()){
    for(int i=0;i<5;i++) aparam[i]=parCont[id][i];
    if(vtxContainer.find(id)!=vtxContainer.end()){
      vtx=vtxContainer[id];
      tof=tofContainer[id];
      length=flContainer[id];
      return true;
    }
  }
  return false;
}
//______________________________________________________________________________
bool CDSTrack::GetNthParameters(const int &n, int &id, double *aparam, TVector3 &vtx)
{
  std::map<int,parContainer>::iterator it = parCont.begin();
  for( int j=0; j<n; j++ ) ++it;
  id=it->first;
  for(int i=0;i<5;i++) aparam[i]=(it->second)[i];
  vtx=vtxContainer[id];

  return true;
}
//______________________________________________________________________________
void CDSTrack::SetParameters(const double *aparam)
{
  for(int i=0;i<5;i++) Param[i]=aparam[i];
}
//______________________________________________________________________________
bool CDSTrack::DeleteHit( const int &layer ,const int &i) 
{
  //  std::cout<<"[CDSTrack::DeleteHit] begin"<<std::endl;
  if(layer<0 || 14<layer) return false;
  std::vector<TrackHit>::iterator it = m_hit[layer].begin();
  for( int j=0; j<i; j++ ) ++it;
  m_hit[layer].erase( it ) ;
  //  std::cout<<"[CDSTrack::DeleteHit] end"<<std::endl;
  return true;
}
//______________________________________________________________________________
CDCTrackHits CDSTrack::get_trackhits()
{
  CDCTrackHits tra;
  tra.set_param(Param);
  for(int l=0;l<15;l++){
    for(int i=0;i<nhits(l);i++){
      tra.addhit(hit(l,i)->layer,hit(l,i)->wire,hit(l,i)->nth);
    }
  }     
  return tra;
}
//______________________________________________________________________________
bool CDSTrack::AddCluster(DCCluster *cl)
{
  for(int ih=0;ih<cl->nhit();ih++){
    AddHit(cl->hit(ih),cl->nth(ih));
  }
  return true;
}
//______________________________________________________________________________
int CDSTrack::nhits() const
{
  int tmpnhits=0;
  for(int i=0 ;i<15;i++)
    {
      tmpnhits+=nhits(i);
    }
  return tmpnhits;
}
//______________________________________________________________________________
double CDSTrack::totsum() const
{
  double tmpsum=0;
  for(int layer=0 ;layer<15;layer++)
    {
      for(int n=0;n<(int)nhits(layer);n++)
	{
	  tmpsum+=m_hit[layer][n].tot;	  
	}
    }
  return tmpsum;
}
//______________________________________________________________________________
int CDSTrack::nsuperlayers()
{
  int tmpnsupers=0;
  for(int i=0 ;i<7;i++) {
    for(int j=0 ;j<3;j++) {
      int layer=super[i][j];
      if(nhits(layer)>0) {
	tmpnsupers++;
	break;
      }
    }
  }
  return tmpnsupers;
}
//______________________________________________________________________________
int CDSTrack::naxialsuperlayers()
{
  int tmpnsupers=0;
  for(int i=0 ;i<3;i++) {
    for(int j=0 ;j<3;j++){
      int layer=axsuper[i][j];
      if(nhits(layer)>0){
	tmpnsupers++;
	break;
      }
    }
  }
  return tmpnsupers;
}
//______________________________________________________________________________
int CDSTrack::nstereosuperlayers()
{
  int tmpnsupers=0;
  for(int i=0 ;i<4;i++) {
    for(int j=0 ;j<3;j++){
      int layer=stsuper[i][j];
      if(nhits(layer)>0){
	tmpnsupers++;
	break;
      }
    }
  }
  return tmpnsupers;
}
//______________________________________________________________________________
int CDSTrack::nhitlayers()
{
  int tmpnlayers=0;
  for(int i=0 ;i<15;i++)
    {
      if(nhits(i)>0) tmpnlayers++;
    }
  return tmpnlayers;
}
//______________________________________________________________________________
int CDSTrack::naxialhits()
{
  int tmpnhits=0;
  for(int i=0 ;i<7;i++)
    {
      int layer=axlayer[i];
      tmpnhits+=nhits(layer);
    }
  return tmpnhits;
}
//______________________________________________________________________________
int CDSTrack::naxialhitlayers()
{
  int tmpnlayers=0;
  for(int i=0 ;i<7;i++)
    {
      int layer=axlayer[i];
      if(nhits(layer)>0) tmpnlayers++;
    }
  return tmpnlayers;
}
//______________________________________________________________________________
int CDSTrack::nstereohits()
{

  int tmpnlayers=0;
  for(int i=0 ;i<8;i++)
    {
      int layer=stereo[i];
      if(nhits(layer)>0) tmpnlayers++;
    }
  return tmpnlayers;
}
//______________________________________________________________________________
int CDSTrack::nstereohitlayers()
{
  int tmpnhits=0;
  for(int i=0 ;i<8;i++)
    {
      int layer=stereo[i];
      tmpnhits+=nhits(layer);
    }
  return tmpnhits;
}
//______________________________________________________________________________
bool CDSTrack::FirstHelixFitting()
{
  // std::cout<<__PRETTY_FUNCTION__<<std::endl;
  // set nearest point on wire to the Circle for each stereo layer in XY plane
  double rc=CircleR();
  double xc=CircleX();
  double yc=CircleY();
  for(int stlayer=0 ;stlayer<8;stlayer++)
    {
      int layer=stereo[stlayer];
      for(int n=0;n<(int)nhits(layer);n++)
	{
	  TrackHit* tmphit=hit(layer,n);	  
	  TVector3 circlepos=math::NearestPointCircleToLine(tmphit->wpos,tmphit->wdir,
							    rc,xc,yc);
	  double dis=(tmphit->wpos-circlepos).Perp(); // distance in R direction
	  double tmpz=tmphit->wpos.z()+dis/tmphit->wdir.Perp()*tmphit->wdir.Z();
	  circlepos.SetZ(tmpz);
	  tmphit->hitpos=circlepos; 
	}
    }  
  // phi-Z plane tracking
  // z = - 1/kappa * dipangle * phi + dz
  // kappa = param[2] in CircleFit
  // dipangle = param[4]
  // dz = param[3]
  int allhit=0;
  double s_x=0,s_y=0,s_xx=0,s_xy=0;
  for(int stlayer=0 ;stlayer<8;stlayer++)
    {
      int layer=stereo[stlayer];
      for(int n=0;n<(int)nhits(layer);n++)
	{
	  TrackHit *tmphit=hit(layer,n);
	  double z=tmphit->hitpos.Z();
	  double phi=math::CalcHelixPhi(tmphit->hitpos.X(),tmphit->hitpos.Y(),Param);
	  phi/=Param[2];
	  s_x+=phi;
	  s_y+=z;
	  s_xx+=phi*phi;
	  s_xy+=phi*z;
	  allhit++;		  
	  //	  std::cout<<layer<<"  "<<n<<"  "<<z<<"  "<<phi<<std::endl;
	}
    }
  double aa,cc;
  aa=(allhit*s_xy-s_x*s_y)/(allhit*s_xx-s_x*s_x);
  cc=(s_xx*s_y-s_xy*s_x)/(allhit*s_xx-s_x*s_x);
  //  bb=-1;
  Param[3]=cc;Param[4]=-aa;
  fittinglevel=3;
  CalcChi2();
  //  SetHitPos();       
  return true;
}
//______________________________________________________________________________
void CDSTrack::CheckCharge()
{
  double initphi=-9999,finiphi=-9999;
  double xc=CircleX();
  double yc=CircleY();
  for(int l=0;l<7;l++)
    {
      int layer=axlayer[l];
      for(int n=0;n<(int)nhits(layer);n++)
	{
	  TrackHit* tmphit=hit(layer,n);
	  TVector3 wirepos=tmphit->wpos+tmphit->wdir*0.5;
	  double tmpphi=math::CalcDeg(wirepos.X()-xc,wirepos.Y()-yc);
	  if(initphi<-1000) initphi=tmpphi;
	  finiphi=tmpphi;
	}
    }
  double phi_diff=finiphi-initphi;
  if( phi_diff>180 ) phi_diff-=360;
  if( phi_diff<-180 ) phi_diff+=360;
  if( phi_diff>0 && Param[2]>0) {
    Param[0]=-Param[0];
    Param[2]=-Param[2];
    Param[4]=-Param[4];
    Param[1]+=TMath::Pi();
  }
  else if( phi_diff<0 && Param[2]<0){
    Param[0]=-Param[0];
    Param[2]=-Param[2];
    Param[4]=-Param[4];
    Param[1]-=TMath::Pi();
  }
  if(Param[1]<0) Param[1]+=2*TMath::Pi();
  else if(Param[1]>2*TMath::Pi()) Param[1]-=2*TMath::Pi();
}
//______________________________________________________________________________
bool CDSTrack::CircleFitting()
{   
  double x[40],y[40],weight[40];
  int nall=0;
  // FirstFitting with MWPC mode
  for(int il=0;il<7;il++)
    {
      int layer=axlayer[il];
      for(int n=0;n<(int)nhits(layer);n++)
	{
	  TrackHit* cdc=hit(layer,n);
	  x[nall]=cdc->wpos.X();
	  y[nall]=cdc->wpos.Y();
	  weight[nall]=3.0*mm;
	  nall++;
	}
    }
#if DEBUG
  for(int i=0;i<nall;i++) std::cout<<i<<" / "<<nall<<", x= "<<x[i]<<" , y= "<<y[i]<<std::endl;
#endif
  double chilr=99999;double chilrtmp=99999;
  double paramlr[5];
  if(!math::CircleFit(x,y,weight,nall,paramlr,chilrtmp)){
    //    std::cout<<"circlefitting: mwpc fitting failed!!!"<<std::endl;
    return false;
  }
#if DEBUG
  std::cout<<"chisqare of circle fit with MWDC mode: "<<chilrtmp<<std::endl;
#endif
  double xc=paramlr[0];
  double yc=paramlr[1];
  //fitting all LR combinations
  int ncombi=pow(2,nall);
  for(int ilr=0;ilr<ncombi;ilr++)
    {
      nall=0;	
      for(int il=0;il<7;il++)
	{
	  int layer=axlayer[il];
	  for(int n=0;n<(int)nhits(layer);n++)
	    {
	      TrackHit* cdc=hit(layer,n);
	      x[nall]=cdc->wpos.X();
	      y[nall]=cdc->wpos.Y();
	      double dl=cdc->cdl;
	      if(dl<0) dl=0.0001*mm;	      
	      double rtmp=sqrt( (x[nall]-xc)*(x[nall]-xc)+(y[nall]-yc)*(y[nall]-yc) );
	      double sintmp= (y[nall]-yc)/rtmp;
	      double costmp= (x[nall]-xc)/rtmp;
	      double sign=0.0;
	      int check=ilr;	      
	      if(check<pow(2,nall) ) sign=1;
	      else 
		{
		  int cflag;
		  for(int n=0;n<nall;n++)
		    {
		      cflag=check%2;
		      check=(check-cflag)/2;
		    }
		  cflag=check%2;
		  if(cflag==0)  sign=1;
		  else if(cflag==1) sign=-1;
		}
	      x[nall]+= sign*dl*costmp;
	      y[nall]+= sign*dl*sintmp;
	      weight[nall]=ResolutionOfCDC;
	      nall++;
	    }
	}
      if(!math::CircleFit(x,y,weight,nall,Param,chilrtmp))
	{
	  continue;
	}
      if(chilrtmp<chilr)
	{
	  chilr=chilrtmp;
	  for(int n=0;n<5;n++) paramlr[n]=Param[n];
	}
    }//LR
  Param[0] = sqrt(paramlr[0]*paramlr[0]+paramlr[1]*paramlr[1]) - paramlr[2] ;
  Param[1] = atan2(paramlr[1],paramlr[0]);
  Param[2] = 1./paramlr[2];
  if(Param[1]<0) Param[1]+=2*TMath::Pi();
  else if(Param[1]>2*TMath::Pi()) Param[1]-=2*TMath::Pi();
  SetChi2(chilr);
#if DEBUG
  std::cout<<"chisquare of circlefitting: "<<chilr<<std::endl;
#endif
  if( chi2() >99999 )     return false;
  CheckCharge();  
  fittinglevel=2;
  //  SetHitPos(cdsMan);
  return true;
}
//______________________________________________________________________________
bool CDSTrack::HelixFitting()
{
  int allhit=0;
  TVector3 wirepos[50],wiredir[50];
  double weight[50],drift[50];

  for(int layer=0;layer<15;layer++)
    {
      for(int n=0;n<(int)nhits(layer);n++)
	{
	  TrackHit* tmphit=hit(layer,n);
	  wirepos[allhit]=tmphit->wpos+tmphit->wdir*0.5;
	  //	  std::cout<<layer<<"  "<<n<<"  "; wirepos[allhit].Print();
	  if( tmphit->wdir.Perp()<0.01 ){ // axial layer
	    wirepos[allhit].SetZ(-9999);
	  }
	  wiredir[allhit]=tmphit->wdir;
	  if(tmphit->resolution<0) tmphit->resolution=ResolutionOfCDC;
	  weight[allhit]=tmphit->resolution;
	  drift[allhit]=tmphit->cdl;
	  allhit++;
	  if(allhit>=50) return false;
	}      
    }    
  HelixFit *helixfit=new HelixFit(Param,wirepos,wiredir,weight,drift,allhit);
  helixfit->GetParameters(Param); 
  Dof= helixfit->dof();
  ChiSquare=helixfit->chisquare()/Dof;
  // if(!helixfit->stat()) return false;
  delete helixfit;
  if( chi2() >9999 )      return false;
  SetHitPos();
  fittinglevel=4;
  return true;
}
//______________________________________________________________________________
bool CDSTrack::FirstLineFitting()
{
  // std::cout<<__PRETTY_FUNCTION__<<std::endl;
  // set nearest point on wire to the Circle for each stereo layer in XY plane
  for(int stlayer=0 ;stlayer<8;stlayer++)
    {
      int layer=stereo[stlayer];
      for(int n=0;n<(int)nhits(layer);n++)
	{
	  TrackHit* tmphit=hit(layer,n);	  
	  TVector3 tmppos=math::NearestPointLineToLine(GetPos0(),GetDir0(),tmphit->wpos,tmphit->wdir,0);
	  double dis=(tmphit->wpos-tmppos).Perp(); // distance in R direction
	  double tmpz=tmphit->wpos.z()+dis/tmphit->wdir.Perp()*tmphit->wdir.Z();
	  tmppos.SetZ(tmpz);
	  tmphit->hitpos=tmppos; //
	}
    }  
  // phi-Z plane tracking
  // z = - 1/kappa * dipangle * phi + dz
  // dipangle = param[4]
  // dz = param[3]
  int allhit=0;
  double s_x=0,s_y=0,s_xx=0,s_xy=0;
  for(int stlayer=0 ;stlayer<8;stlayer++)
    {
      int layer=stereo[stlayer];
      for(int n=0;n<(int)nhits(layer);n++)
	{
	  TrackHit *tmphit=hit(layer,n);
	  double z=tmphit->hitpos.Z();
	  double tmpr=tmphit->hitpos.Perp();
	  s_x+=tmpr;
	  s_y+=z;
	  s_xx+=tmpr*tmpr;
	  s_xy+=z*tmpr;
	  allhit++;		  
	}
    }
  double aa,cc; //z = aa*r+cc
  aa=(allhit*s_xy-s_x*s_y)/(allhit*s_xx-s_x*s_x);
  cc=(s_xx*s_y-s_xy*s_x)/(allhit*s_xx-s_x*s_x);
  //  bb=-1;
  Param[3]=aa*GetPos0().Perp()+cc; 
  Param[4]=TMath::Pi()/2.-TMath::ATan(aa);
  fittinglevel=3;
  CalcChi2();
  return true;
}
//______________________________________________________________________________
bool CDSTrack::LineFitting()
{
  int allhit=0;
  TVector3 wirepos[50],wiredir[50];
  double weight[50],drift[50];

  for(int layer=0;layer<15;layer++)
    {
      for(int n=0;n<(int)nhits(layer);n++)
	{
	  TrackHit* tmphit=hit(layer,n);
	  wirepos[allhit]=tmphit->wpos+tmphit->wdir*0.5;
	  wiredir[allhit]=tmphit->wdir;
	  if(tmphit->resolution<0) tmphit->resolution=ResolutionOfCDC;
	  weight[allhit]=tmphit->resolution;
	  drift[allhit]=tmphit->cdl;
	  allhit++;
	  if(allhit>=50) return false;
	}      
    }    
  
  LineFit *linefit=new LineFit(Param,wirepos,wiredir,weight,drift,allhit);
  linefit->GetParameters(Param); 
  Dof= linefit->dof();
  ChiSquare=linefit->chisquare()/Dof;
  // if(!helixfit->stat()) return false;
  delete linefit;
  if( chi2() >9999 )      return false;
  SetHitPos();
  fittinglevel=4;
  return true;
}
//______________________________________________________________________________
bool CDSTrack::LineAxialFitting()
{
  double x[40],y[40],weight[40];
  int nall=0;
  // FirstFitting with MWPC mode
  double tmpphi=Param[2];
  for(int il=0;il<7;il++)
    {
      int layer=axlayer[il];
      for(int n=0;n<(int)nhits(layer);n++)
	{
	  TrackHit* cdc=hit(layer,n);
	  TVector3 tmppos=cdc->wpos;
	  tmppos.RotateZ(-tmpphi);
	  x[nall]=tmppos.X();
	  y[nall]=tmppos.Y();
	  weight[nall]=5.0*mm;
	  nall++;
	}
    }
  double tmpchi2=99999;
  double param[2];
  if(!math::LineFit(x,y,weight,nall,param,tmpchi2)){
    return false;
  }
  // y = p0 + p1 *x
  //fitting all LR combinations
  double sintmp=1/sqrt(1+param[1]*param[1]);
  double costmp=-sintmp*param[1];
#if DEBUG
  std::cout<<"chisqare of line fit with MWDC mode: "<<tmpchi2<<std::endl;
  std::cout<<"   p0: "<<param[0]
	   <<",  p1: "<<param[1]<<std::endl;
  for( int i=0; i<nall; i++ ){
    double tmpy=param[0]+param[1]*x[i];
    double tmp=pow((y[i]-tmpy)/weight[i],2);    
    std::cout<<x[i]<<"  "<<y[i]<<"  "<<tmpy<<"  "<<tmp<<std::endl;
  } 
#endif
  int ncombi=pow(2,nall);
  double tmpminchi2=99999;
  double tmpparam[2];  
  for(int ilr=0;ilr<ncombi;ilr++)
    {
      nall=0;	
      for(int il=0;il<7;il++)
	{
	  int layer=axlayer[il];
	  for(int n=0;n<(int)nhits(layer);n++)
	    {
	      TrackHit* cdc=hit(layer,n);
	      TVector3 tmppos=cdc->wpos;
	      tmppos.RotateZ(-tmpphi);
	      x[nall]=tmppos.X();
	      y[nall]=tmppos.Y();
	      double dl=cdc->cdl;
	      if(dl<0) dl=0.0001*mm;	      
	      double sign=0.0;
	      int check=ilr;	      
	      if(check<pow(2,nall) ) sign=1;
	      else 
		{
		  int cflag;
		  for(int n=0;n<nall;n++)
		    {
		      cflag=check%2;
		      check=(check-cflag)/2;
		    }
		  cflag=check%2;
		  if(cflag==0)  sign=1;
		  else if(cflag==1) sign=-1;
		}
	      y[nall]+= sign*dl;
	      weight[nall]=ResolutionOfCDC;
	      nall++;
	    }
	}
      if(!math::LineFit(x,y,weight,nall,tmpparam,tmpchi2))
	{
	  continue;
	}
      if(tmpchi2<tmpminchi2)
	{
	  tmpminchi2=tmpchi2;
	  for(int n=0;n<2;n++) param[n]=tmpparam[n];
	}
    }//LR
  // y = p0 + p1 *x
  // -p1*x + y - p0 =0;
  // pos (0,p0), dir (1,p1)
  double k=-param[0]*param[1]/(1+param[1]*param[1]);
  TVector3 tmppos0(k,param[0]+k*param[1],0);
  tmppos0.RotateZ(tmpphi);
  Param[0] = tmppos0.X();
  Param[1] = tmppos0.Y();
  Param[2] = TMath::ATan(param[1])+tmpphi;
  Param[3]=0;
  Param[4]=TMath::Pi()/2.;
  SetChi2(tmpminchi2);
  if( chi2() >99999 )     return false;
  fittinglevel=2;
  return true;
}
//______________________________________________________________________________
void CDSTrack::SetHitPos()
{    
  for(int layer=0;layer<15;layer++)
    {
      for(int n=0;n<(int)nhits(layer);n++)
	{
	  TrackHit *tmphit=hit(layer,n);	  
	  double dis;
	  TVector3 lnest,hnest;
	  double dl=tmphit->cdl;
	  if(dl<0) dl=0.002*mm;
	  if(IsLine()){
	    math::LineToLine( GetPos0(), GetDir0(), tmphit->wpos, tmphit->wdir, dl, dis, hnest, lnest);
	  }else{	    
	    if(!math::LineToHelix(tmphit->wpos,tmphit->wdir,Param,lnest,hnest,dis) )
	      { 		  
		tmphit->residual=999*cm;
		tmphit->hitpos=TVector3(999,999,999);
		continue;
	      }
	    TVector3 tmp=(hnest-lnest).Unit();
	    lnest+=dl*tmp;
	  }
	  tmphit->residual	= dl-dis;
	  tmphit->hitpos	= lnest;
	  tmphit->trackpos	= hnest;
	}
    }
}
//______________________________________________________________________________
TVector3 CDSTrack::GetPositionatR(const double &r)
{
  if(IsLine()){
    return math::CalcLinePosatR(GetPos0(),GetDir0(),r);
  }
  return math::CalcHelixPosatR(Param,r);
}
//______________________________________________________________________________
TVector3 CDSTrack::GetMomentumVector(const TVector3 &pos)
{
  double tmppar[5];
  GetParameters(tmppar);
  return GetMomentumVector(pos,tmppar);
}
//______________________________________________________________________________
TVector3 CDSTrack::GetMomentumVector(const TVector3 &pos, const double *aparam)
{
  double phi=math::CalcHelixPhi(pos.x(),pos.y(),aparam);  
  double pttmp = fabs(pt(aparam));
  double px = pttmp*(-1*sin(aparam[1]+phi));
  double py = pttmp*(cos(aparam[1]+phi));
  double pz = pttmp*(aparam[4]);
  TVector3 p; p.SetXYZ(px,py,pz);
  return p;
}
//______________________________________________________________________________
bool CDSTrack::Calc()
{
  double global[5];
  gCDC.HelixLocalToGlobal(Param,global,charge());
  double tmpr=15.3*cm+(53.0-15.3)*cm/2.; 
  TVector3 tmppos=math::CalcHelixPosatR(global,tmpr);
  TVector3 tmppos2=math::CalcHelixPosatR(Param,tmpr);
#if 0
  for(int i=0;i<5;i++){
    std::cout<<std::setw(10)<<Param[i]<<" -> "<<std::setw(10)<<global[i]<<std::endl;
  }
  std::cout<<"local : ";tmppos2.Print();
  std::cout<<"global: ";tmppos.Print();
#endif
  parCont.clear();
  vtxContainer.clear();
  AddParameters(DetIdCDC,global,tmppos,0,0);
  return true;
}
//______________________________________________________________________________
void CDSTrack::AddParameters(const int &id,const double *aparam, const TVector3 &vtx,const double &tof,const double &length)
{
  parContainer tmp;
  for(int i=0;i<5;i++) tmp.push_back(aparam[i]);
  parCont[id]=tmp;
  vtxContainer[id]=vtx;
  tofContainer[id]=tof;
  flContainer[id]=length;
}
//______________________________________________________________________________
void CDSTrack::CalcChi2()
{  
  Double_t chisq=0.;
  Int_t dof = 0;
  TVector3 fittmp,fittmp2;
  Double_t fitx, fity, dist;
  for(int layer=0;layer<15;layer++)
    {
      for(int n=0;n<(int)nhits(layer);n++)
	{
	  TrackHit* tmphit=hit(layer,n);
	  TVector3 wirepos=tmphit->wpos+tmphit->wdir*0.5;
	  TVector3 wiredir=tmphit->wdir;
	  double weight=ResolutionOfCDC;
	  double drift=tmphit->cdl;
	  if(IsLine()){		
	    math::LineToLine( wirepos, wiredir, GetPos0(), GetDir0(), 0, dist, fittmp, fittmp2);
	  }else{
	    if( tmphit->wdir.Perp()<0.01 ){ // axial layer
	      math::PointToCircle(Param, wirepos.x(), wirepos.y(), fitx, fity,dist);
	    }else{
	      math::LineToHelix( wirepos, wiredir, Param, fittmp, fittmp2, dist);
	    }
	  }
	  chisq += pow( (dist-drift)/weight , 2 );
	  dof++;
	}      
    }    
  Dof=  dof - 5;
  ChiSquare = chisq;
}
//______________________________________________________________________________
double CDSTrack::pt(const double *par) const 
{
  return magneticfield/m/m*math::C()*mm/par[2]; // T
}
//______________________________________________________________________________
double CDSTrack::mom(const double *par) const 
{
  return pt(par)*sqrt(1+par[4]*par[4]); 
}
//______________________________________________________________________________
int CDSTrack::pid_tot(){
  if(pt()>0){
    if(mom()>(tot_pid_param1+tot_pid_param2*totave2() )){
      return kProton; //proton
    }else{
      return kPiPlus; //pi plus
    }
  }else{
    return kPiMinus; //pi minus
  }
}
//______________________________________________________________________________
bool CDSTrack::GetMomentum(const TVector3 &pos, const double &mass, TVector3 &p,bool ELOSS, bool GLOBAL)
{
  double tmp1,tmp2; 
  return GetMomentum(pos,mass,p,tmp1,tmp2,ELOSS,GLOBAL);
}
bool CDSTrack::GetMomentum(const TVector3 &pos, const double &mass, TVector3 &p, double &tof,double &length,bool ELOSS, bool GLOBAL)
{
  //  std::cout<<"========= CDSTrack::GetMomentum"<<std::endl;
  // calculate momentum at given pos.
  // if you would like to consider dE correction, you should set parameters before hand
  // by calling CalcVertexTimeLength()
  TVector3 initpos;
  double tmppar[5];
  if(ELOSS){
    int volid=geom::GetID(pos);
    if(!GetParameters(volid,tmppar,initpos,tof,length)){
      // std::cout<<"[CDSTrack::GetMomentm] no parameters found for detid = "<<volid<<std::endl;
      if(!CalcEnergyLoss(pos,mass,tmppar,initpos) ) return false;      
      if(!GetParameters(volid,tmppar,initpos,tof,length)) return false;
    }
  }
  else if(GLOBAL){
    if(!GetParameters(DetIdCDC,tmppar,initpos,tof,length)) return false;
  }
  else GetParameters(tmppar); // CDC local without energy loss
  
  if(ELOSS){
    double momout;
    double tmptof,tmpl;
    if(!eloss::CalcHelixElossPointToPoint(tmppar,initpos,pos,mom(tmppar),mass,momout,tmptof,tmpl)){
      std::cout<<"[CDSTrack::GetMomentm] Failed in CalcHelixElossPtoP"<<std::endl;
      return false;
    }
    tof+=tmptof;
    length+=tmpl;
    double newpar[5];
    math::ChangePivot(TVector3(0,0,0),pos,tmppar,newpar,charge());
    // T = V*s/m2
    newpar[2]=magneticfield/m/m*math::C()*mm/(momout/sqrt(1+newpar[4]*newpar[4]));
    math::ChangePivot(pos,TVector3(0,0,0),newpar,tmppar,charge());
  }
  p=GetMomentumVector(pos,tmppar);
  return true;

}
//______________________________________________________________________________
bool CDSTrack::GetVertex(const TVector3 &pos, const TVector3 &dir, TVector3 &lpos, TVector3 &hpos)
{
  // vertex search considering curvature change by dE
  // return true if the vertex
  TVector3 tmpvtx;
  TVector3 tmp1,tmp2;
  double tmp[5];
  int id;
  double dis, tmpdis=999;
  //  std::cout<<"---------------"<<std::endl;
  for(int i=0;i<nParamSets();i++){
    GetNthParameters(i,id,tmp,tmpvtx);
    if(!math::LineToHelix(pos,dir,tmp,tmp1,tmp2,dis) ){
#if DEBUG
      std::cout<<"math::LineToHelix failed in CDSTrack::GetVertex"<<std::endl;
#endif
      continue;
    }
#if DEBUG
    std::cout<<"=== id,getid,dis "<<id<<"  "<<geom::GetID(tmp2)<<"  "<<dis<<std::endl;
    tmpvtx.Print();
    tmp1.Print();
    tmp2.Print();
#endif
    if(dis<tmpdis||geom::GetID(tmp2)==id){
      lpos=tmp1;
      hpos=tmp2;
      tmpdis=dis;
      if(geom::GetID(tmp2)==id) return true;
    }
  }
  if(tmpdis<100) return true;
  //  if(lpos.Mag()<0.01)  Print();
  return false;
}
//______________________________________________________________________________
bool CDSTrack::CalcVertexTimeLength(const TVector3 &beam_pos,const TVector3 &beam_dir,const double &mass,TVector3 &lpos, TVector3 &hpos,double &tof, double &length, bool ADDPAR)
{
  //  std::cout<<"[CDSTrack::CalcVertexTimeLength] start !!! "<<mass<<"  "<<mom()<<std::endl;
  if(mass<=0 || TMath::Abs(mom())>10){
    std::cout<<"[CDSTrack::CalcVertexTimeLength] failed mass mom !!! "<<mass<<"  "<<mom()<<std::endl;
    return false;
  }
  int sign=1; // add energy loss
  double momout=mom();
  double momin=momout;
  double tmppar[5],newpar[5],newpar2[5];
  int id=DetIdCDC;
  TVector3 pos_in,pos_out;
  GetParameters(id,tmppar,pos_in,tof,length);
  TString mat=geom::GetMaterial(pos_in);
  length=0; tof=0;

  int count=0;
  TString newmat;
  double tmptof,tmpl,dis;
  while(count<30){
    //
    if(!math::LineToHelix(beam_pos,beam_dir,tmppar,lpos,hpos,dis) ){
      std::cout<<"[CDSTrack::CalcVertexTimeLength] faile in LineToHelix"<<std::endl;
      return false;
    }
    if(!geom::HelixStepToNextVolume(tmppar,pos_in,pos_out,tmpl,newmat,id)){
      std::cout<<"[CDSTrack::CalcVertexTimeLength] faile in HelixStepToNextVolume"<<std::endl;
      return false; 
    }
#if 0
    std::cout<<std::setw(15)<<geom::GetVolName(pos_in)
	     <<" ("<<std::setprecision(5)<<std::setw(8)<<pos_in.X()
	     <<"," <<std::setprecision(5)<<std::setw(8)<<pos_in.Y()
	     <<"," <<std::setprecision(5)<<std::setw(8)<<pos_in.Z()<<")"
	     <<"->"<<std::setw(15)<<geom::GetVolName(pos_out)
	     <<" ("<<std::setprecision(5)<<std::setw(8)<<pos_out.X()
	     <<"," <<std::setprecision(5)<<std::setw(8)<<pos_out.Y()
	     <<"," <<std::setprecision(5)<<std::setw(8)<<pos_out.Z()<<")"
	     <<std::endl;
#endif
    double tmpl2=math::CalcHelixArc(tmppar,pos_in,hpos);
    if(!eloss::CalcdE(momin,mass,tmpl,sign,mat,momout,tmptof) ){
      std::cout<<"[CDSTrack::CalcVertexTimeLength] faile in CalcdE"<<std::endl;
      return false;
    }
    math::ChangePivot(TVector3(0,0,0),pos_out,tmppar,newpar,charge());
    newpar[2]=magneticfield/m/m*math::C()*mm/(momout/sqrt(1+newpar[4]*newpar[4])); // check later
    math::ChangePivot(pos_out,TVector3(0,0,0),newpar,newpar2,charge());
    if(tmpl>tmpl2)   break; // cross the vertex point
    if(id==DetIdCDCCFRP||id==DetIdCDC){
      // the length and tof are calcurated for between vertex and the CDC inner wall
      length=0;
      tof=0;
    }else{
      length+=tmpl;
      tof+=tmptof;
    }
    if(ADDPAR) AddParameters((int)id,newpar2,pos_out,tof,length);
    for(int i=0;i<5;i++) tmppar[i]=newpar2[i];
    pos_in=pos_out;
    momin=momout;
    mat=newmat;
    count++;
  }
  mat=newmat;
  tmpl=math::CalcHelixArc(tmppar,pos_in,hpos);
  eloss::CalcdE(momin,mass,tmpl,sign,mat,momout,tmptof);
  if(id!=DetIdCDCCFRP&&id!=DetIdCDC){
    tof+=tmptof;
    length+=tmpl;
  }
  return true;
}
//______________________________________________________________________________
bool CDSTrack::CalcEnergyLoss(const TVector3 &pos, const double &mass, double *par, TVector3& pos_out)
{
  if(mass<=0 || TMath::Abs(mom())>10){
    std::cout<<"#W too large momentum or negative mass"<<std::endl;
    return false;
  }
  int sign=1; // add energy loss
  double tmppar[5],newpar[5];
  int id=DetIdCDC;
  int id_goal=geom::GetID(pos);
  double tof,length;
  GetParameters(id,par,pos_out,tof,length);
  int count=0;
  TString newmat;
  double tmptof,tmpl,momout;
  double dis=(pos_out-pos).Mag();
  while(count++<30){
    TVector3 pos_in=pos_out;
    for(int i=0;i<5;i++) tmppar[i]=par[i];
    TString mat=geom::GetMaterial(pos_in);
    if(!geom::HelixStepToNextVolume(tmppar,pos_in,pos_out,tmpl,newmat,id)){
      std::cout<<"[CDSTrack::CalcEnergyLoss] failed in HelixStepToNextVolume"<<std::endl;
      return false; 
    }
    //    std::cout<<id<<"  "<<mat<<" ->  "<<newmat<<std::endl;
    if(GetParameters(id,par,pos_in,tof,length)){
      if((pos_in-pos).Mag()<(pos_out-pos).Mag())  continue;      
    }
    double momin=mom(tmppar);
    if(!eloss::CalcdE(momin,mass,tmpl,sign,mat,momout,tmptof) ){
      std::cout<<"[CDSTrack::CalcEnergyLoss] failed in CalcdE"<<std::endl;
      return false;
    }
    math::ChangePivot(TVector3(0,0,0),pos_out,tmppar,newpar,charge());
    newpar[2]=magneticfield/m/m*math::C()*mm/(momout/sqrt(1+newpar[4]*newpar[4])); // check later
    math::ChangePivot(pos_out,TVector3(0,0,0),newpar,par,charge());
    if(id==DetIdCDCCFRP||id==DetIdCDC){
      // the length and tof are calcurated for between vertex and the CDC inner wall
      length=0;
      tof=0;
    }else{
      length+=tmpl;
      tof+=tmptof;
    }
    AddParameters((int)id,par,pos_out,tof,length);
    if(id==id_goal) return true;
    double tmpdis=(pos_out-pos).Mag();
    if(tmpdis<dis) dis=tmpdis;
    else return true;
  }
  return false;
}
//______________________________________________________________________________
bool CDSTrack::Retiming(double cdhtime, double beta, bool SLEW)
{
  // cdh time with respect to the beta=1 cds particle
  for(int layer=0;layer<15;layer++){
    for(int n=0;n<(int)nhits(layer);n++)
      {
	TrackHit *cdc= hit(layer,n);
	double ctime=0;
	if(0){
	  double dis=math::CalcHelixArc(Param,cdc->hitpos.Perp(),cdh_rin);
	  double time_cdc=dis/beta/(TMath::C()*1e-6*mm);
	  ctime = cdhtime-time_cdc;
	  cdc->cdt=cdc->dt-ctime+cdc->hitpos.Z()/30.; // 30*cm/ns
	}
	if(SLEW){
#if UNIDAQ
	  double obeta2=1./pow(beta,2);
	  gSlew.DoCorrection(DetIdCDC,cdc->layer,cdc->wire,0,cdc->dt,obeta2,ctime);
#else
	  gSlew.DoCorrection(DetIdCDC,cdc->layer,cdc->wire,0,cdc->dt,cdc->tot,ctime);
#endif
	  cdc->cdt=ctime;
	}
	cdc->cdl=gXT.CalcDriftLength( DetIdCDC , cdc->layer , cdc->wire, cdc->cdt );
      }
  }  
  HelixFitting();
  return true;
}
//______________________________________________________________________________
bool CDSTrack::RemoveBadHits(double threshold_sigma,int threshold_nhits)
{
  while(nhits()>threshold_nhits&&nsuperlayers()==7){
    int tmpl=-1,tmpn=-1;
    double tmpresid=-1;
    for(int layer=0;layer<15;layer++){
      for(int n=0;n<(int)nhits(layer);n++)
	{
	  TrackHit *cdc= hit(layer,n);
	  if(fabs(cdc->residual)>tmpresid){
	    tmpl=layer;
	    tmpn=n;
	    tmpresid=fabs(cdc->residual);
	  }
	}
    }    
#if DEBUG
    std::cout<<nhits()<<"  "<<nsuperlayers()
	     <<"  "<<tmpl<<"   "<<tmpn<<"  "<<tmpresid<<std::endl;
#endif
    if(tmpresid/ResolutionOfCDC<threshold_sigma) break;
    DeleteHit(tmpl,tmpn);
    if(IsLine()){
      LineFitting();
    }else{
      HelixFitting();
    }
  }
  return true;
}
//______________________________________________________________________________
void CDSTrack::Print()
{
  std::cout<<"=================================="<<std::endl;
  std::cout<<"TrackID      : "<<trackid()<<std::endl;
  std::cout<<"Fitting level: "<<fittinglevel<<std::endl;
  std::cout<<"Chi2      : "<<chi2()<<std::endl;
  std::cout<<"Params[5]:"<<std::endl;
  for(int i=0;i<5;i++) std::cout<<"  "<<Param[i];
  std::cout<<std::endl;
  if(IsLine()){
    std::cout<<"Pos0: "; GetPos0().Print();
    std::cout<<"Dir0: "; GetDir0().Print();
  }else{
    std::cout<<"CirCenterX: "<<CircleX()<<std::endl;
    std::cout<<"CirCenterY: "<<CircleY()<<std::endl;
    std::cout<<"CirRho    : "<<CircleR()<<std::endl;
  }
  for(int lay=0;lay<15;lay++){
    for(int i=0;i<nhits(lay);i++){
      if(i>100) break;
      std::cout<<"layer: "<<std::setw(5)<<lay
	       <<std::setw(3)<<i<<" / "
	       <<std::setw(3)<<nhits(lay)
	       <<std::setw(5)<<hit(lay,i)->wire
	       <<std::setprecision(3)<<std::setw(12)<<hit(lay,i)->residual;
      print_TV3(hit(lay,i)->wpos);
      print_TV3(hit(lay,i)->trackpos);
      std::cout<<std::endl;
    }
  }
  
  for(int i=0;i<nParamSets();i++){
    int id=-1;
    double tmp[5]={0,0,0,0,0};
    TVector3 pos;
    GetNthParameters(i,id,tmp,pos);
    std::cout<<"id  "<<id<<"  ";
    for(int i=0;i<5;i++) std::cout<<tmp[i]<<"  ";
    std::cout<<pos.Perp()<<"  ";
    pos.Print();   
  }
}
