//#include "Parameters.hh"
#include "Utility_Helix.hh"
#include "TMath.h"
#include <TRandom.h>
#include <TVector3.h>
#include <iostream>
#include <sstream>
#include <string>


// // TVector3 Utility_Helix::GaussPosition( TVector3 tv )
// {
//   double x = 0. , y = 0., z = 0.;//target center
//   double tvx = tv.x();
//   double tvy = tv.y();
//   double tvz = tv.z();
//   double sigx = TPCResolutionXPara;
//   double sigy = TPCResolutionYPara;
//   double sigz = TPCResolutionZPara;

//   static TRandom *random = new TRandom();
//   //random->SetSeed((int) tvx);

//   tvx += random->Gaus(0.0, sigx);
//   tvy += random->Gaus(0.0, sigy);
//   tvz += random->Gaus(0.0, sigz);

//   return TVector3 (tvx, tvy, tvz);

// }

double Utility_Helix::CalcRad( const double x, const double y )
{
  double theta=atan( y/x );
  if(x<0 && y>=0) theta=TMath::Pi()+theta;
  else if(x<0 && y<0) theta=TMath::Pi()+theta;
  else if(x>=0 && y<0) theta=2.0*TMath::Pi()+theta;
  return theta;
}

double Utility_Helix::CalcHelixPhi(const double &x, const double &y, const double *par)
{
  // return value between -pi and pi                                                                        
  double sin_c, cos_c;
  double cx=(par[0]+1./par[2])*cos(par[1]);
  double cy=(par[0]+1./par[2])*sin(par[1]);

  cos_c=x-cx;
  sin_c=y-cy;
  double phi;
  phi=Utility_Helix::CalcRad(cos_c, sin_c);
  // phi=atan2(sin_c, cos_c);

  if(phi>TMath::Pi()) phi -= 2*TMath::Pi();
  if(par[2] < 0 ) phi -= par[1];
  else phi -= (par[1]+TMath::Pi() );
  while(fabs(phi)>TMath::Pi())
    {
      if(phi>0) phi -= 2*TMath::Pi();
      else phi += 2.*TMath::Pi();
    }

  return phi;
}

bool Utility_Helix::HelixToHelix(const double *par1, const double  *par2,TVector3 &fitpos1,TVector3 &fitpos2, double &dis )
{
  double posphi=Utility_Helix::CalcHelixPhi(0,0,par1);
  double distmp=999;
  TVector3 pos,tmp;
  //for(int scale=0;scale<3;scale++)
  for(int scale=0;scale<5;scale++)
    {
      double tmpphi=posphi;
      for(int n=0;n<150;n++)
	//for(int n=0;n<30;n++)
        {
          double phi=tmpphi+(n-10)*pow(10,-scale)*3*par1[2];
          pos=Utility_Helix::GetPosition(phi,par1);
          if(!Utility_Helix::PointToHelix(pos,par2,tmp,distmp) ){distmp=999;}
          if(n==0) { dis=distmp; posphi=phi; fitpos2=tmp;}
          else if(distmp<dis) { dis=distmp; posphi=phi; fitpos2=tmp;}
        }
    }
  fitpos1=Utility_Helix::GetPosition(posphi,par1);
  //  double phi2=Utillity::CalcHelixPhi(fitpos2,par2);
  //  std::cout<<"posphi = "<<posphi<<"  "<<phi2<<std::endl;
  //  if(dis<10) return true;
  //  else return false;
  return true;
}

TVector3 Utility_Helix::GetPosition(const double &helixphi, const double *par)
{
  TVector3 pos;
  if(par[2]!=0){
    // This is the equation for Helix formula
    pos.SetX( par[0]*cos(par[1])+1./par[2]*( cos(par[1])-cos(par[1]+helixphi) ) );
    pos.SetY( par[0]*sin(par[1])+1./par[2]*( sin(par[1])-sin(par[1]+helixphi) ) );
    pos.SetZ( par[3]-1./par[2]*par[4]*helixphi );
  }
  else if(par[2]==0)
    {
      pos.SetXYZ(par[0]*cos(par[1])-helixphi*(-sin(par[1]) ),
                 par[0]*sin(par[1])-helixphi*( cos(par[1]) ),
                 par[3]+par[4]*helixphi);
    }

  return pos;
}

bool Utility_Helix::PointToHelix(const TVector3 &hitpos, TVector3 &fitpos, const double *par)
{
  double phi=0., phi_b=0., phi_a = 0.;
  phi=Utility_Helix::CalcHelixPhi(hitpos.x(), hitpos.y(), par); //helix phi, between -pi and pi

  double cx=( par[0]+1./par[2] )*cos(par[1]);
  double cy=( par[0]+1./par[2] )*sin(par[1]);
  double dz=par[3];
  double rho=fabs(1./par[2]);

  double philen = 1./par[2]*sqrt(1. + par[4]*par[4]);
  //double philen = sqrt( (cx-hitpos.x())*(cx-hitpos.x()) + (cy-hitpos.y())*(cy-hitpos.y()) + (dz-hitpos.z())*(dz-hitpos.z()) )  -  rho;
  double dist = 1.;
  int trial = 14;

  while(dist < 128*128)
    {
      phi_b = phi-dist/philen;
      phi_a = phi+dist/philen;
      double dlen_b = dfunc_PTH(hitpos, phi_b, par);
      double dlen_a = dfunc_PTH(hitpos, phi_a, par);
      if(dlen_b*dlen_a<=0.) break;
      else 
	{ 
	  dist *= 2.;
	  trial++;
	}
    }
  
  if(dist >= 128*128)
    {
      std::cout<<"[HelixFit] Can not find PTH initial param!!"<<std::endl;
      //std::cout<<"[HelixFit] philen: "<<philen<<std::endl;
      fitpos.SetXYZ(9999, 9999, 9999);
      return false;
    }

  //Bisection Method
  for( int i = 0; i< trial; i++)
    {
      phi_b=phi-dist/philen;
      phi_a=phi+dist/philen;
      double dlen_b=dfunc_PTH(hitpos, phi_b, par);
      double dlen=dfunc_PTH(hitpos, phi, par);
      if(dlen*dlen_b <= 0) 
	{
	  phi=(phi_b+phi)/2.0;
	  dist=dist/2.0;
	}
      else 
	{
	  phi=(phi_a+phi)/2.0;
	  dist=dist/2.0;
	}
    }
  fitpos=GetPosition(phi, par);
  return true;
}

double Utility_Helix::dfunc_PTH( const TVector3 &pos, const double &helixphi, const double *par)
{
  double inv=1./par[2];
  double x  = par[0]*cos(par[1]) + inv*(cos(par[1]) - cos(par[1]+helixphi));
  double y  = par[0]*sin(par[1]) + inv*(sin(par[1]) - sin(par[1]+helixphi));
  double z  = par[3] - inv*par[4]*helixphi;
  double dx = inv*sin(par[1]+helixphi);
  double dy = -inv*cos(par[1]+helixphi);
  double dz = -inv*par[4];

  double S = 2.*( (x - pos.x())*dx + (y - pos.y())*dy + (z - pos.z())*dz);
  return S;
}

bool Utility_Helix::PointToHelix(const TVector3 &hitpos, const double *par,TVector3 &fitpos ,double &dis)
{
  if(!PointToHelix(hitpos,fitpos,par)) return false;
  TVector3 tmp;
  tmp=fitpos-hitpos;
  dis=tmp.Mag();
  return true;
}

//GeV/c unit
TVector3 Utility_Helix::CalcHelixMom(const double par[5], double z)
{
  const double Const = 0.299792458; // =c/10^9                                                                             
  const double dMagneticField = 1*FieldPara; //T, "-1" is needed.                                                                   
  double phi2 = (par[3]-z)*par[2]/par[4];
  double pt = fabs(1/par[2])*(Const*dMagneticField)/1000.; // GeV/c
  double px = pt*(-1*sin(par[1]+phi2));
  double py = pt*(cos(par[1]+phi2));
  double pz = pt*(par[4]);
  return TVector3(px,py,pz);
}

bool Utility_Helix::EnergyLossBethe( const double distance, const TLorentzVector lv, TLorentzVector &lv_corr)
{
  if( distance == 0.)
    {
      lv_corr = lv;
      return true;
    }

  double dE_dx; /*MeV/cm*/
  double eloss; /*MeV*/

  double beta;
  double thickness = distance/10.0; /*cm*/
  double m = lv.M()*1000.0;  /*mass of incident particle(Mev)*/

  double E_0;
  double E;
  TVector3 v3_p(lv.X(), lv.Y(), lv.Z());
  TVector3 p_dir (v3_p.x()/v3_p.Mag(), v3_p.y()/v3_p.Mag(), v3_p.z()/v3_p.Mag());

  double p = lv.P()*1000.0;
  //double delta=0.01; /*cm*/                         
  double delta=0.1; /*cm*/
  int i;
  double length=0.0;
  double total_eloss=0.0;

  double p_corr;
  double E_corr;

  //E_0=MyGamma(beta)*m; 
  //p=MyGamma(beta)*beta*m;

  E_0= sqrt(m*m + p*p);
  E=E_0;
  beta = Utility_Helix::MyBeta(E,p);
  //printf("beta=%f, E=%f, p=%f\n",beta, E, p);

  //std::cout<<"[Bethe] mass: "<< m <<std::endl;

  if (beta<=0.0) 
    {
      p_corr = p/1000.0;
      E_corr = E/1000.0;
      return true;
    }
  
  dE_dx=0.0;
  eloss=0.0;
  for(i=0;i<=thickness/delta;i++)
    {
      dE_dx=calc_dE_dx(beta);
      eloss=dE_dx*delta;
      E=E+eloss;
      if(E<m)
	{
	  //fprintf(stderr,"particle stops in material at %5.3fcm\n",length);
	  p_corr = p/1000.0;
	  E_corr = E/1000.0;
	  return false;
	  //break;                                                                                                                 
	}
      p=sqrt(pow(E,2.0)-pow(m,2.0));
      beta=Utility_Helix::MyBeta(E,p);
      length=length+delta;
      total_eloss=total_eloss+eloss;
#if 0
      printf("beta:%5.3f\n",beta);
      printf("dE_dx:%5.3f\teloss:%5.3f\n",dE_dx,eloss);
      printf("E:%5.3f(MeV)\tp:%5.3f(MeV/c)\n",E,p);
      printf("length:%5.3f(cm)\n",length);
      printf("total energy loss:%5.3f(MeV)\n",total_eloss);                                        
#endif
      //getchar();                                                                                     
    }
  p_corr = p/1000.0;
  E_corr = E/1000.0;
  
  TLorentzVector LV_corr(p_dir.x()*p_corr, p_dir.y()*p_corr, p_dir.z()*p_corr, E_corr);
  lv_corr = LV_corr;
  
  #if 0
  std::cout<<"[Bethe formula] momentum change: "<< p - lv.P()*1000.0  <<" MeV/c"<<std::endl;
  std::cout<<"[Bethe formula] distance: "<< distance  <<" mm"<<std::endl;
  #endif


  return true;

}

double Utility_Helix::calc_dE_dx(double beta)
{
  double value;
  const double C=0.1535; /*MeVcm^2/g*/
  const double m_e=0.511;
  double logterm;
  double gamma_2;
  double W_max;
  double gamma;
  double X;
  double delta;

  //double rho=0.0709;   /*g/cm^3 (C)*/
  double rho=TargetDensity;   /*g/cm^3 (C)*/
  double I=21.8;     /*eV*/
  double Z_A=0.99216;
  int z=1;
  double C0=-3.2632;
  double X0=0.4759;
  double X1=1.9215;
  double a=0.13483;
  double M=5.6249;

  gamma = Utility_Helix::MyGamma(beta);
  X = log10(beta*gamma);
  if (X<=X0)
    delta=0.0;
  else if (X0<X && X<X1)
    delta=4.6052*X+C0+a*pow((X1-X),M);
  else if (X>=X1)
    delta=4.6052*X+C0;

  gamma_2=1/(1-pow(beta,2.0));
  //printf("sqrt(gamma_2):%f\n",sqrt(gamma_2));                                                                                

  W_max=2.0*m_e*pow(beta,2.0)*gamma_2;

  logterm=log(2.0*m_e*gamma_2*pow(beta,2.0)*W_max*pow(10.0,12.0)/pow(I,2.0))-2.0*pow(beta,2.0)-delta/2.0;
  //logterm=log(2.0*m_e*gamma_2*pow(beta,2.0)*W_max*pow(10.0,12.0)/pow(I,2.0))-2.0*pow(beta,2.0);

  value=C*rho*Z_A*pow((double)z,2.0)*logterm/pow(beta,2.0);

  return value;

}

double Utility_Helix::MyGamma(double beta)
{
  double value;
  value=1.0/sqrt(1.0-pow(beta,2.0));
  return value;
}

double Utility_Helix::MyBeta(double energy,double mormentum)
{
  double value;
  value=mormentum/energy;
  return value;
}

double Utility_Helix::CalcFlightLengthInTarget( const double par[5], const TVector3 vertex )
{
  const double Const = 0.299792458; // =c/10^9                                                        
  const double dMagneticField = 1*FieldPara; //T, "-1" is needed.                                      

  double vertex_phi = Utility_Helix::CalcHelixPhi(vertex.x(), vertex.y(), par);
  double radius = fabs(1/par[2]);
  double phi_unit = 1.0/radius; // mimimum length 1.0 mm
  double phi = vertex_phi;

  if(par[2] > 0.) phi_unit = -1.*phi_unit;

  double length_inside_target = 0.0;
  double total_length = 0.0;
  
  double phia = phi;
  double phib = phi + phi_unit;
  TVector3 posa = GetPosition(phia, par);
  TVector3 posb = GetPosition(phib, par);
  
  int num_try = 0;
  bool check_targeta = CheckInTarget(posa);
  bool check_targetb = CheckInTarget(posb);

  TVector3 target_pos(0., -143., 0.);
  double helix_phi = CalcHelixPhi(0., -143., par);
  TVector3 fit_pos;
  double dist;
  /*
  Utility_Helix::PointToHelix(target_pos, par, fit_pos, dist);

  if( fabs(sqrt(dist*dist - fit_pos.z()*fit_pos.z())) > 30.)
    {
      return 0.0;
    }

  if( fabs(fit_pos.z()) > 50.)
    {
      return 0.0;
    }
  */
  while( num_try < 200 )
    {
      check_targeta = CheckInTarget(posa);
      check_targetb = CheckInTarget(posb);
      if(check_targeta == true && check_targetb ==false)
	{
	  break;
	}
      double arc = Utility_Helix::CalcHelixArc(par, posa, posb);
      total_length += arc;
      if(CheckInTarget(posa) == true && CheckInTarget(posb) == true )
	{
	  length_inside_target += arc;
	}
      
      GetPosition( phi, par );
      phia += phi_unit;
      phib += phi_unit;
      posa = GetPosition(phia, par);
      posb = GetPosition(phib, par);
      
#if 0
      std::cout<<"[HelixArcTarget] try: "<<num_try << std::endl;
      std::cout<<"pos x: "<<posa.x()<< " pos y: "<<posa.y() << " pos z: "<< posa.z() << std::endl;
      std::cout<<"check inside target: "<< CheckInTarget(posa) << std::endl;
      std::cout<<" arc: "<< arc<<" Length in target: " 
	       << length_inside_target << " Total length: "<< total_length <<std::endl;
#endif 
      num_try++;
    }

#if 0
  std::cout<<"[HelixArcTarget] try: "<<num_try << std::endl;
  std::cout<<"closest distance: "<< dist << std::endl;
  std::cout<<"pos x: "<<posa.x()<< " pos y: "<<posa.y() << " pos z: "<< posa.z() << std::endl;
  std::cout<<"check inside target: "<< CheckInTarget(posb) << std::endl;
  std::cout<<"Length in target: " 
	   << length_inside_target << " Total length: "<< total_length <<std::endl;
#endif 

  return length_inside_target;

}

double Utility_Helix::CalcHelixArc(const double par[5], const TVector3 &pos1,const TVector3 &pos2)
{
  /* double phi1 = par[2]/par[4]*(par[3]-pos1.Z()); */
  /* double phi2 = par[2]/par[4]*(par[3]-pos2.Z()); */
  double phi1=CalcHelixPhi(pos1, par);
  double phi2=CalcHelixPhi(pos2, par);
  double rho= fabs(1./par[2]);
  double larc=TMath::Abs(rho*(phi1-phi2))*sqrt(1+par[4]*par[4]);
  return larc;
}

bool Utility_Helix::CheckInTarget(TVector3 pos)
{
  double distt = sqrt(pos.x()*pos.x() + (pos.y()-TargetPositionZ)*(pos.y()-TargetPositionZ));
  if( fabs(distt) < TargetRadius && fabs(pos.z()) < TargetHeight)
    {
      return true;
    }
  else return false;
}

bool Utility_Helix::CheckInTargetC(TVector3 pos)
{
  if( fabs(pos.x()) < TargetXC && fabs(pos.z()) < TargetYC && fabs(pos.y()-TargetPositionZ) < TargetZC)
    {
      return true;
    }
  else return false;
}



double Utility_Helix::CalcHelixPhi(const TVector3 &pos, const double *par)
{
  return Utility_Helix::CalcHelixPhi(pos.X(),pos.Y(),par);
}

TVector3 Utility_Helix::GetVertex(TVector3 v1, TVector3 v2)
{
  TVector3 vertex((v1.x()+v2.x())/2., (v1.y()+v2.y())/2., (v1.z()+v2.z())/2.);
  return vertex;
}

TVector3 Utility_Helix::BeamPMod( TVector3 tv )
{

  const double sigma = 0.002;
  double p = tv.Mag();

  static TRandom *random = new TRandom();
  double p_mod = random->Gaus(p, p*sigma);

  TVector3 tv_mod(tv.x()/tv.Mag()*p_mod, tv.y()/tv.Mag()*p_mod, tv.z()/tv.Mag()*p_mod);
  return tv_mod;
}

double Utility_Helix::TofEMod( double edep )
{

  const double sigma = 0.2;

  static TRandom *random = new TRandom();
  double edep_mod = random->Gaus(edep, edep*sigma);

  return edep_mod;
}

void Utility_Helix::LineToLine(const TVector3 &x1, const TVector3 &p1, const TVector3 &x2, const TVector3 &p2, TVector3 &vp, double &dist)
{
  // X1 = x1 + t*a1                                                                  
  // X2 = x2 + s*a2                                                                  
  // (X1-X2)*a1=0                                                                    
  // (X1-X2)*a2=0                                                                    
  //    ||                                                                           
  //    \/                                                                           
  // a*t + b*s = A1                                                                  
  // c*t + d*s = A2                                                                  
  // xest : on x1,a1                                                                 
  TVector3 x = x2-x1;
  double a =  p1.Dot(p1);
  double b = -p1.Dot(p2);
  double c =  p2.Dot(p1);
  double d = -p2.Dot(p2);
  double A1 = p1.Dot(x);
  double A2 = p2.Dot(x);

  double D = a*d-b*c;

  TVector3 x2p;
  TVector3 xest;

  if( fabs(D)<0.00000000000000001 )
    {
      dist = sqrt(x.Mag2()-A1*A1);
      TVector3 vertex_point( (x1.x()+x2.x())/2., (x1.y()+x2.y())/2., (x1.z()+x2.z())/2.);
      vp = vertex_point;
    }
  else
    {
      double s = (a*A2-c*A1)/D;
      double t = (d*A1-b*A2)/D;
      xest  = x1 + t*p1;
      x2p   = x2 + s*p2;
      TVector3 vertex_point( (xest.x()+x2p.x())/2., (xest.y()+x2p.y())/2., (xest.z()+x2p.z())/2.);
      vp = vertex_point;
      dist = (xest-x2p).Mag();
    }

}

double Utility_Helix::CalcFlightLengthInTargetC( const double par[5], const TVector3 vertex )
{
  const double Const = 0.299792458; // =c/10^9                                                        
  const double dMagneticField = 1*FieldPara; //T, "-1" is needed.                                      

  double vertex_phi = Utility_Helix::CalcHelixPhi(vertex.x(), vertex.y(), par);
  double radius = fabs(1/par[2]);
  double phi_unit = 1.0/radius; // mimimum length 1.0 mm
  double phi = vertex_phi;

  if(par[2] > 0.) phi_unit = -1.*phi_unit;

  double length_inside_target = 0.0;
  double total_length = 0.0;
  
  double phia = phi;
  double phib = phi + phi_unit;
  TVector3 posa = GetPosition(phia, par);
  TVector3 posb = GetPosition(phib, par);
  
  int num_try = 0;
  bool check_targeta = CheckInTargetC(posa);
  bool check_targetb = CheckInTargetC(posb);

  TVector3 target_pos(0., -143., 0.);
  double helix_phi = CalcHelixPhi(0., -143., par);
  TVector3 fit_pos;
  double dist;
  /*
  Utility_Helix::PointToHelix(target_pos, par, fit_pos, dist);

  if( fabs(sqrt(dist*dist - fit_pos.z()*fit_pos.z())) > 30.)
    {
      return 0.0;
    }

  if( fabs(fit_pos.z()) > 50.)
    {
      return 0.0;
    }
  */
  while( num_try < 200 )
    {
      check_targeta = CheckInTargetC(posa);
      check_targetb = CheckInTargetC(posb);
      if(check_targeta == true && check_targetb ==false)
	{
	  break;
	}
      double arc = Utility_Helix::CalcHelixArc(par, posa, posb);
      total_length += arc;
      if(CheckInTargetC(posa) == true && CheckInTargetC(posb) == true )
	{
	  length_inside_target += arc;
	}
      
      GetPosition( phi, par );
      phia += phi_unit;
      phib += phi_unit;
      posa = GetPosition(phia, par);
      posb = GetPosition(phib, par);
      
#if 0
      std::cout<<"[HelixArcTarget] try: "<<num_try << std::endl;
      std::cout<<"pos x: "<<posa.x()<< " pos y: "<<posa.y() << " pos z: "<< posa.z() << std::endl;
      std::cout<<"check inside target: "<< CheckInTargetC(posa) << std::endl;
      std::cout<<" arc: "<< arc<<" Length in target: " 
	       << length_inside_target << " Total length: "<< total_length <<std::endl;
#endif 
      num_try++;
    }

#if 0
  std::cout<<"[HelixArcTarget] try: "<<num_try << std::endl;
  std::cout<<"closest distance: "<< dist << std::endl;
  std::cout<<"pos x: "<<posa.x()<< " pos y: "<<posa.y() << " pos z: "<< posa.z() << std::endl;
  std::cout<<"check inside target: "<< CheckInTargetC(posb) << std::endl;
  std::cout<<"Length in target: " 
	   << length_inside_target << " Total length: "<< total_length <<std::endl;
#endif 

  return length_inside_target;

}



bool Utility_Helix::LineToHelix(const TVector3 &aa, const TVector3 &dlinea, const double *par, TVector3 &lnest, TVector3 &hnest, double &dis)
{
  // aa: origin of the line                                                                                   
  // dline: direction of the line                                                                            
  // par: helix parameter                                                                                    
  // lnest: nearest point on the line                                                                        
  // hnest: nearest point on the helix                                                                       
  TVector3 a(aa.x(), aa.y(), aa.z());
  TVector3 dline(dlinea.x(), dlinea.y(), dlinea.z());

  TVector3 dlineu=dline.Unit();
  double phi=0,phi_b=0,phi_a=0;
  phi=Utility_Helix::CalcHelixPhi(a.x(),a.y(),par);

  double philen=1./par[2]*sqrt(1+par[4]*par[4]);
  double dist=1.;
  int trial=14; //sigma_position<1.2 micron                                                                  

  // find LTH initial param                                                                                  
  while(dist<128*128)
    {
      phi_b=phi-dist/philen;
      phi_a=phi+dist/philen;
      double dlen_b=dfunc_LTH(a,dlineu,phi_b,par);
      double dlen_a=dfunc_LTH(a,dlineu,phi_a,par);
      if(dlen_b*dlen_a<=0) break;
      else {dist*=2;trial++;}
    }

  if(dist>=128*128)
    {
#if 1
      std::cout<<"[Utility_Helix]LineToHelix Can not find LTH inital param!!"<<std::endl;
#endif
      dis=-9999.;
      return false;
    }
  //Bisection Method                                                                                         
  for(int i=0;i<trial;i++)
    {
      phi_b=phi-dist/philen;
      phi_a=phi+dist/philen;
      double dlen_b=dfunc_LTH(a,dlineu,phi_b,par);
      //      double dlen_a=dfunc_LTH(a,dlineu,phi_a,par);                                                   
      double dlen=dfunc_LTH(a,dlineu,phi,par);
      if(dlen*dlen_b<=0 ) {phi=(phi_b+phi)/2.0;dist=dist/2.0;}
      else  {phi=(phi_a+phi)/2.0;dist=dist/2.0;}
    }
  //                                                                                                         
  // helix nearest point                                                                                    

 
  hnest=Utility_Helix::GetPosition(phi,par);
  double k=(hnest-a)*dlineu;
  // line nearest point                                                                                      
  lnest=a+k*dlineu;
  TVector3 tmp;
  tmp=hnest-lnest;
  dis=tmp.Mag();



  return true;
}



double Utility_Helix::dfunc_LTH(const TVector3 &lpos,const TVector3 &dline,const double &helixphi,const double *par)
{
  // function to return the differential (dphi) of                                                   
  // the distance between Line and a point on Helix which corresponds to helixphi                   
  //                                                                                                  
  // lpos  : hit position on the line                                                                      
  // dline : direction of the line                                                                  
  // helixphi :phi correspondig to hit position                                                        
  // par : parameter of helix                                                                   
  // drho=par[0], phi0=par[1], rho=par[2], dz=par[3], tlam=par[4];                                          
  TVector3 hpos,dhpos;
  hpos=Utility_Helix::GetPosition(helixphi,par);
  // dhpos: momentum direction on helix                                                                     
  dhpos.SetXYZ( (1./par[2])*sin(par[1]+helixphi)  ,
                (-1./par[2])*cos(par[1]+helixphi) ,
                (-1./par[2])*par[4]                );
  double k=hpos*dline-lpos*dline; // minimum distance point on the line                               
  double dk=dhpos*dline; // dk/dphi                                                              
  double S = 2.*(lpos+k*dline-hpos)*(dk*dline-dhpos); // d( (lpos+k*dline) -hpos)/d phi          
  return S;
}
