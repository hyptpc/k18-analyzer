#include "ELossTools.hh"
#include "GeomTools.hh"
#include "MathTools.hh"
#include <iomanip>
#define DEBUG 0
namespace{
  const double mm=0.1;
  const double cm=10*mm;
  const double um=1e-3*mm;
  const double default_cdc_step=1*cm;
}
//______________________________________________________________________________
bool eloss::CalcHelixElossToVertex(const double param[5],const TVector3 &vertex, const double &momin, const double &mass, double &momout,double &tof,const double &rstart){
  //  int sign=1; // add energy loss
  TVector3 pos1=math::CalcHelixPosatR(param,rstart);
  double tmp;
  return eloss::CalcHelixElossPointToPoint(param,pos1,vertex,momin,mass,momout,tof,tmp);
}
//______________________________________________________________________________
bool eloss::CalcHelixElossToNextBoundary(const double param[5],TVector3 &pos1,const double &momin, const double &mass,
					 double &momout,  double &tof, double &length,int &id, int sign)
{
  if(pos1.Mag()>700*mm) return false;// out of CDC
  TString mat;
  TVector3 pos2;
  if(!geom::HelixStepToNextVolume(param,pos1,pos2,length,mat,id)) return false;
  geom::GetIDMat(pos1,mat);
  pos1=pos2;
  bool flag=eloss::CalcdE(momin,mass,length,sign,mat,momout,tof);
  return flag;
}
//______________________________________________________________________________
bool eloss::CalcElossBeam(const TVector3 &in, const TVector3 &out, 
			  const double &momin, const double &mass, 
			  double &momout,double &tof){
  bool flag=eloss::CalcElossBeam(in,out,momin,mass,-1,momout,tof);
  return flag;
}
//______________________________________________________________________________
bool eloss::CalcElossBeam(const TVector3 &in, const TVector3 &out, 
			  const double &momin, const double &mass,const int &sign, 
			  double &momout,double &tof){
  double tmpmom=momin;
  TString mat=geom::GetMaterial(in);
  TVector3 tmppos=in; 
  TString newmat;
  double length;
  double sumlength=0;
  tof=0;
  double tmptof=-999;
  int count=0;
  while(count<1000&&
	geom::StepToNextVolume(tmppos,out,length,newmat)){
    if(!eloss::CalcdE(tmpmom,mass,length,sign,mat,momout,tmptof,false))  return false;
    tof+=tmptof;
    tmpmom=momout;
    mat=newmat;
    sumlength+=length;
    tmppos+=(out-tmppos).Unit()*length;
    count++;
  }
  // if(newmat.IsNull()) return false; // I do not understand why this line is, 20170817. if in and out positions are in the same volume, then newmat would be not defined.
  length=(out-tmppos).Mag();//-sumlength;
  if(!eloss::CalcdE(tmpmom,mass,length,sign,mat,momout,tmptof,false)) return false;
  tof+=tmptof;
#if DEBUG
  //#if 1
  std::cout<<"  "<<mass<<"  "<<mat<<"  "<<length<<"  "<<(momout-tmpmom)*1000<<"  "<<tmptof<<"  "<<tof<<std::endl;
#endif
  //  exit(-1);
  return true;
}
//______________________________________________________________________________
bool eloss::CalcHelixElossPointToPoint(const double param[5],const TVector3 &posin, const TVector3 &posout, const double &momin, const double &mass, double &momout,double &tof,double &totl,double offs){
  //offs???
  if(mass<=0){
    std::cout<<"#W mass 0 or negative"<<std::endl;
    return false;
  }
  // definde step  
  int sign=1; // add energy loss
  double phiin=math::CalcHelixPhi(posin,param);
  double phiout=math::CalcHelixPhi(posout,param);
  if(phiin*param[2]>phiout*param[2]) sign=-1; // 1: backto vertex, -1: outgoing

  double step=default_cdc_step;
  double tmptof,length=0;
  double tmpmomin=momin;
  int count=0;
  totl=0;
  tof=0;
  momout=momin;  

  TVector3 pos1=posin;
  TVector3 pos2=math::CalcHelixStep(param,pos1,sign*step);
  double phi1=math::CalcHelixPhi(pos1,param);
  double phi2=math::CalcHelixPhi(pos2,param);
  int id1=geom::GetID(pos1);
  int id2=geom::GetID(pos2);

  while(count<50){
    if(pos1.Mag()>700*mm){      
      std::cout<<"#W track too far away from the target"<<std::endl;
#if DEBUG
      pos1.Print();
      posin.Print();
      posout.Print();
#endif
      return false; //out of CDC  
    }
    // calculate the length within the same volume
    while(step>50*um&&(phi1-phiout)*(phi2-phiout)>0){
      if(geom::IsSameVolume(pos1,pos2)){
	length+=step;
	pos1=pos2;
	pos2=math::CalcHelixStep(param,pos1,sign*step);
	phi1=math::CalcHelixPhi(pos1,param);
	phi2=math::CalcHelixPhi(pos2,param);
      }else{
	step/=10.;
	pos2=math::CalcHelixStep(param,pos1,sign*step);
	phi1=math::CalcHelixPhi(pos1,param);
	phi2=math::CalcHelixPhi(pos2,param);
      }    
      if(length>100*cm){
	std::cout<<"#W track too long"<<std::endl;
	return false; // too long
      }
    }
    // caculate dE
    if((phi1-phiout)*(phi2-phiout)<0){ // reached to the final position
      length+=math::CalcHelixArc(param,pos1,posout);
      //length-=offs; // ??
      if(length<0) length=10*um;
      TString mat=geom::GetMaterial(pos1);
      eloss::CalcdE(tmpmomin,mass,length,sign,mat,momout,tmptof);
      tof+=tmptof;
      totl+=length;
      return true;
    }else{ // cross to the next volume
      step*=10;
      TString mat=geom::GetMaterial(pos1);
      pos2=math::CalcHelixStep(param,pos1,sign*step);
      //      step*=sign; // to be positive??
      //   geom::CrossBoundary(pos1,pos2,step,tmpmat,tmpstep,tmp);
      //      length+=step;
      double tmpstep=geom::CrossBoundary(pos1,pos2,step);
      length+=tmpstep; // step or tempstep ???
      eloss::CalcdE(tmpmomin,mass,length,sign,mat,momout,tmptof);
      tof+=tmptof;
      totl+=length;
      // set parameters for the next step
      length=step-tmpstep;
      tmpmomin=momout;
      step=default_cdc_step;
      pos1=math::CalcHelixStep(param,pos1,sign*(tmpstep+1*um));
      pos2=math::CalcHelixStep(param,pos1,sign*step);
      id1=geom::GetID(pos1);
      id2=geom::GetID(pos2);
      phi1=math::CalcHelixPhi(pos1,param);
      phi2=math::CalcHelixPhi(pos2,param);
      count++;
    }
  }
  std::cout<<"#W too many loops"<<std::endl;
  return false;
}
//______________________________________________________________________________
bool eloss::CalcdE(const double param[5], const TVector3 &vertex, const double &rmax, const double &rmin,
		   const double &momin, const double &mass, const TString &material, 
		   double &momout, double &tof)
{
  // for helix rmax->rmin or rmax->vertex
  int sign=1;
  double rvtx=vertex.Perp();
  if(rvtx>rmax) return false;
  double length=0;
  if(rvtx>rmin){
    length=math::CalcHelixArc(param,rmax,rvtx);
  }else{
    length=math::CalcHelixArc(param,rmax,rmin);
  }
  return eloss::CalcdE(momin,mass,length,sign,material,momout,tof);
}
//______________________________________________________________________________
bool eloss::CalcdE(const double &momin,const double &mass, const double &length,const int &sign,
		   const TString &matname,
		   double &momout, double &tof,const bool &CORR)
{
  // length in cm
  momout=momin;
  double tmpmass=mass;
  if(mass<=0) tmpmass=511e-6; // electron mass
  double Ein=sqrt(tmpmass*tmpmass+momin*momin);
  double pm= momin>0 ? 1. : -1;
  double beta=TMath::Abs(momin)/Ein;

  double step=500*um; //cm
  if(tmpmass>0.5&&momin<0.4) step=100*um; 
  if(tmpmass<0.5&&momin<0.2) step=100*um; 
  if(matname.Contains("Air")||matname.Contains("Gas")||matname.Contains("Vacuum")) step*=100;
  double toteloss=0.;
  double totlength=0.;
  tof=0;
  double betamid=0.;
  for(int i=1;i<=length/step;i++){
    double dE=eloss::CalcdEdX(beta,matname,CORR); // GeV
    if(dE<0)  dE=0;//return false;
    dE*=step*cm;
    double Eout=Ein+sign*dE; //step
    if(Eout<tmpmass||beta<=0.0){ // stop
      momout*=pm;
      return false;
    }    
    momout=sqrt(Eout*Eout-tmpmass*tmpmass);
    betamid=(beta+momout/Eout)*0.5;
    tof+=step/(betamid*TMath::C()*1e-6*mm); // m/s -> mm/ns
    beta=momout/Eout;
    Ein=Eout;
    toteloss+=dE;
    totlength+=step;
  }
  double dE=eloss::CalcdEdX(beta,matname,CORR); // per 1cm
  if(dE<0 ) dE=0;
  dE*=(length-totlength)*cm;
  //  totlength+=(length-step*i);
  double Eout=Ein+sign*dE;
  if(Eout<tmpmass||beta<=0.0){
    momout*=pm;
    return false;
  }
  momout=sqrt(Eout*Eout-tmpmass*tmpmass);
  betamid=(beta+momout/Eout)*0.5;
  tof+=(length-totlength)/(betamid*TMath::C()*1e-6*mm);
  toteloss+=dE;
  momout*=pm;
  
#if DEBUG
  //#if 1
  double tmpein=sqrt(mass*mass+momin*momin)-mass;
  double tmpeout=sqrt(mass*mass+momout*momout)-mass;
  std::cout  <<std::setw(10)<<matname
	     <<"  mass "<<std::setw(10)<<mass
	     <<"  totlength/length "<<std::setw(10)<<totlength
	     <<" / "<<std::setw(10)<<length
    //	     <<"   "<<std::setw(10)<<momin<<" -> "<<std::setw(10)<<momout
	     <<"   "<<std::setw(10)<<tmpein<<" -> "<<std::setw(10)<<tmpeout
	     <<" eloss (MeV) = "<<std::setw(15)<<toteloss*1e3
	     <<" time "<<std::setw(10)<<tof<<std::endl;
#endif
  return true;  
}
//______________________________________________________________________________
double eloss::CalcdEdX(const double &beta,const TString &matname, bool CORR){
  // calculate dE per 1cm
  double rho,I,Z_A;
  double Z,C0,a,M,X1,X0;
  // rho, A,Z from knucl3_ag: KnuclMaterials
  // I,C0,a,M,X1,X0 from Atomic Data and Nuclear Data Tables 30, 261-271 (1984)
  if( !matname.CompareTo("Aluminum") ){
    rho=2.7;//g/cm3
    I=166.;// eV
    Z_A=0.48181;
    C0=-4.2395;
    a=0.08024;
    M=3.6345;
    X1=3.0127;
    X0=0.1708;
    Z=13.;
  }
  else if(!matname.CompareTo("CDCGas") ){
    // Ar-Ethane 50-50, simply take avarage
    rho=(1.356*0.5+1.782*0.5)/1000.;
    I=(188.+45.4)/2.;
    Z_A=(0.45059+0.59861)/2.;  
    C0=-(11.9480+9.1043)/2.;
    a=(0.19714+0.09627)/2.;
    M=(2.9618+3.6095)/2.;
    X1=(4.4855+3.8743)/2.;
    X0=(1.7635+1.5107)/2.;
    Z=(18.+18/8.)/2.;
  }
  else if(!matname.CompareTo("CFRP") ){
    rho=1.700;
    I=78;
    Z_A=0.49954;  
    C0=-3.1550;
    a=0.20762;
    M=2.9532;
    X1=2.5387;
    X0=0.0480;
    Z=6.;
  }
  else if(!matname.CompareTo("Scinti") || !matname.CompareTo("Plastic")){
    rho=1.032;
    I=64.7;
    Z_A=0.54141;  
    C0=-3.1997;
    a=0.16101;
    M=3.2393;
    X1=2.4855;
    X0=0.1464;
    Z=3.5;
  } 
  else if(!matname.CompareTo("Aerogel") ){
    //SiO2, take medium of Si & O
    rho=0.2;
    I=(173+95)/2.;
    Z_A=0.5;  
    C0=-(10.7+4.435)/2.;
    a=(0.15+0.118)/2.;
    M=(3.255+3.29)/2.;
    X1=(2.8715+4.3213)/2.;
    X0=(0.2014+1.7541)/2.;
    Z=11;
  } 
  else if(!matname.CompareTo("Mylar") ){
    rho=1.39;
    I=78.7;
    Z_A=0.52037;
    C0=-3.3262;
    a=0.12679;
    M=3.3076;
    X1=2.6507;
    X0=0.1562;
    Z=50./96.*11.;
  }
  else if(!matname.CompareTo("Beryllium") ){
    rho=1.848;
    I=63.7;
    Z_A=0.44384;  
    C0=-2.7847;
    a=0.80392;
    M=2.4339;
    X1=1.6922;
    X0=0.0592;
    Z=4.;
  }
  else if(!matname.CompareTo("AlBeMet") ){
    //Al 38% Be 62%
    rho=2.071;
    I	=63.7	*0.62+166.   *0.38;
    Z_A	=0.44384*0.62+0.48181*0.38;  
    C0	=-2.7847*0.62-4.2395 *0.38;
    a	=0.80392*0.62+0.08024*0.38;
    M	=2.4339	*0.62+3.6345 *0.38;
    X1	=1.6922	*0.62+3.0127 *0.38;
    X0	=0.0592	*0.62+0.1708 *0.38;
    Z	=4.	*0.62+13     *0.38;
  }
  else if(!matname.CompareTo("Iron") ){
    rho=7.874;
    I=286.;
    Z_A=0.46556;  
    C0=-4.2911;
    a=0.14680;
    M=2.9632;
    X1=3.1531;
    X0=-0.0012;
    Z=26.;
  }
  else if(!matname.CompareTo("LHelium-3") ){
#if E73
    rho=0.07;
#else
    rho=0.08;
#endif
    I=41.8;
    Z_A=2./3.;  
    C0=-11.1393;
    a=0.13443;
    M=5.8347;
    X1=3.6122;
    X0=2.2017;
    Z=2.;
  }
  else if(!matname.CompareTo("LHelium-4") ){
    // temporary set the same values as lHe3, except for rho,Z_A
    rho=0.145;
    I=41.8;
    Z_A=2./4.;  
    C0=-11.1393;
    a=0.13443;
    M=5.8347;
    X1=3.6122;
    X0=2.2017;
    Z=2.;
  }
  else if(!matname.CompareTo("LHydrogen") ){
    rho=0.071;
    I=21.8;
    Z_A=1./1.;  
    C0=-3.2632;
    a=0.13483;
    M=5.6249;
    X1=1.9215;
    X0=0.4759;
    Z=1.;
  }
  else if(!matname.CompareTo("LDeuterium") ){
    rho=0.168;
    I=21.8;
    Z_A=1./2.;  
    C0=-3.2632;
    a=0.13483;
    M=5.6249;
    X1=1.9215;
    X0=0.4759;
    Z=1.;
  }
  else if(!matname.CompareTo("Air") || !matname.CompareTo("BLDCGas") ){
    rho=1.2048/1000.;
    I=85.7;
    Z_A=0.49919;  
    C0=-10.5961;
    a=0.10914;
    M=3.3994;
    X1=4.2759;
    X0=1.7418;
    Z=6.8;
  }
  else if(!matname.CompareTo("Vacuum") ){
    rho=1.2048/100000.;
    I=85.7;
    Z_A=0.49919;  
    C0=-10.5961;
    a=0.10914;
    M=3.3994;
    X1=4.2759;
    X0=1.7418;
    Z=6.8;
  }else if(!matname.CompareTo("PbF2") || !matname.CompareTo("PbG") ){
    return 0;
  }else{
    std::cout<<"Material name "<<matname<<" is not defined"<<std::endl;
    //    exit(1);
    return 0;
  }
  double val2=eloss::CalcdEdX(beta,rho,I,Z_A,Z,C0,a,M,X0,X1,CORR); // per 1cm
#if 0 //DEBUG
  std::cout<<matname<<"  beta "<<beta<<"  rho "<<rho<<"  I "<<I<<" eloss = "<<val2<<std::endl;
#endif
  return val2;
}
//______________________________________________________________________________
double eloss::CalcdEdX(const double &beta,const double &rho, const double &I, const double &Z_A, const double &Z,
		       const double &C0, const double &a, const double &M, const double &X0, const double &X1, bool CORR){
  // calculate dE per 1cm
  const double C=0.1535; /*MeVcm^2/g*/
  const double m_e=0.511; // MeV/c^2

  double gamma_2=1/(1-pow(beta,2.0));
  double gamma=sqrt(gamma_2);
  double W_max=2.0*m_e*pow(beta,2.0)*gamma_2;
  
  double X=0.0;
  double C_shell=0.0;
  double delta=0.0;
  double corr=0.0;
  if(CORR){
    double eta=beta*gamma;
    C_shell=(0.422377/pow(eta,2)+0.0304043*pow(eta,-4)-0.00038106*pow(eta,-6))*1e-6*pow(I,2)
      +(3.850190*pow(eta,-2)-0.1667989*pow(eta,-4)+0.00157955*pow(eta,-6))*1e-9*pow(I,3);
    X = log10(eta);
    if (X<=X0)
      delta=0.0;
    else if (X0<X && X<X1) 
      delta=4.6052*X+C0+a*pow((X1-X),M);
    else if (X>=X1)  
      delta=4.6052*X+C0;
    corr=-delta-2*C_shell/Z;
  }
  
  //  std::cout<<"delta,C_shell   "<<delta<<"  "<<C_shell<<std::endl;
  int z=1;
  double logterm=log(2.0*m_e*gamma_2*pow(beta,2.0)*W_max*pow(10.0,12.0)/pow(I,2.0))-2.0*pow(beta,2.0)+corr;
  double val=C*rho*Z_A*pow((double)z,2.0)*logterm/pow(beta,2.0);
  
  return val/1000.; // MeV->GeV
}
