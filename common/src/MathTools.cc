// MathTools.cc
#include "MathTools.hh"

namespace math
{
  bool ChangePivot(const TVector3 &oldpivot, const TVector3 &newpivot,
		   const double oldpar[5], double newpar[5], const int &charge)
  {
    // par[0]=d_rho, par[1]= phi_0, par[2]= kappa/alpha, par[3]=dz, par[4]=tanL
    TVector3 diff=oldpivot-newpivot;
    double tmppar[5];
    for(int i=0;i<5;i++)  tmppar[i]=oldpar[i];
    double tmpx=diff.X() + (tmppar[0]+1./tmppar[2])*TMath::Cos(tmppar[1]);
    double tmpy=diff.Y() + (tmppar[0]+1./tmppar[2])*TMath::Sin(tmppar[1]);
    newpar[1]=TMath::ATan2(tmpy,tmpx) + TMath::Pi()/2.*( 1 - tmppar[2]/TMath::Abs(tmppar[2]) );

    if(newpar[1]-oldpar[1]>TMath::Pi()){
      newpar[1]-=TMath::TwoPi();
    }
    if(newpar[1]-oldpar[1]<-TMath::Pi()){
      newpar[1]+=TMath::TwoPi();
    }

    newpar[0]=tmpx*TMath::Cos(newpar[1])+tmpy*TMath::Sin(newpar[1])-1./tmppar[2];
    newpar[2]=tmppar[2];
    newpar[3]=diff.Z()+tmppar[3]-1./tmppar[2]*(newpar[1]-tmppar[1])*tmppar[4];
    newpar[4]=tmppar[4];

    return true;
  }

  bool HelixToHelix(const double *par1, const double  *par2,TVector3 &fitpos1,TVector3 &fitpos2, double &dis )
  { 
    double posphi=math::CalcHelixPhi(0,0,par1); 
    double distmp=9999;
    TVector3 pos,tmp;
    for(int scale=0;scale<3;scale++)
      {
	double tmpphi=posphi;
	for(int n=0;n<20;n++)
	  {	  
	    double phi=tmpphi+(n-10)*pow(10,-scale)*3*par1[2];
	    pos=math::GetPosition(phi,par1);
	    if(!math::PointToHelix(pos,par2,tmp,distmp) ){distmp=9999;}
	    if(n==0) { dis=distmp; posphi=phi; fitpos2=tmp;}
	    else if(distmp<dis) { dis=distmp; posphi=phi; fitpos2=tmp;}
	  }
      }
    fitpos1=math::GetPosition(posphi,par1);
    return true;
  }

  bool HelixToHelixWresl(const double *par1, const double  *par2,TVector3 &fitpos1,TVector3 &fitpos2, double &dis )
  { 
    double posphi=math::CalcHelixPhi(0,0,par1); 
    double distmp=999;
    double rdistmp=999;
    double rdis=999;
    TVector3 pos,tmp;
    for(int scale=0;scale<3;scale++)
      {
	double tmpphi=posphi;
	for(int n=0;n<20;n++)
	  {	  
	    double phi=tmpphi+(n-10)*pow(10,-scale)*3*par1[2];
	    pos=math::GetPosition(phi,par1);
	    if(!math::PointToHelix(pos,par2,tmp,distmp) ){distmp=999;}
	    TVector3 rpos=pos-tmp;
	    rdistmp=sqrt( (rpos.X()*rpos.X()+rpos.Y()*rpos.Y())/0.3/0.3+rpos.Z()*rpos.Z()/1.0);
 
	    if(n==0) { dis=distmp; posphi=phi; fitpos2=tmp; rdis=rdistmp;}
	    else if(rdistmp<rdis) { dis=distmp; posphi=phi; fitpos2=tmp; rdis=rdistmp;}
	  }
      }
    fitpos1=math::GetPosition(posphi,par1);
    return true;
  }
  
  TVector3 NearestPointCircleToLine(const TVector3 &pos, const TVector3 &dir,
				    const double &rho ,const double &xc, const double &yc)
  {
    double aw,bw,cw;
    // aw*x+bw*y+cw=0	 	    
    if(TMath::Abs( dir.X() )>0.01)
      {
	// y = dy/dx *x + C
	aw=dir.Y()/dir.X();
	bw=-1;
	cw= -(aw*pos.X()+bw*pos.Y() );
      }
    else
      {
	// check !!
	bw=dir.X()/dir.Y();
	aw= -1;
	cw= -(aw*pos.X()+bw*pos.Y() ); 
      }
    double x_p,y_p,x_n,y_n;
    //    std::cout<<__func__<<"aw,bw,cw,rho,xc,yc: "<<aw<<"  "<<bw<<"  "<<cw<<"  "<<rho<<"  "<<xc<<"  "<<yc<<std::endl;

    LineToCircle(aw,bw,cw,rho,xc,yc,x_p,y_p,x_n,y_n);
    //    std::cout<<__func__<<"xp,yp,xn,yn: "<<x_p<<"  "<<y_p<<"  "<<x_n<<"  "<<y_n<<std::endl;
    TVector3 Pos_p(x_p,y_p,0);
    TVector3 Pos_n(x_n,y_n,0);
    if( (Pos_p-pos).Perp() < (Pos_n-pos).Perp() )
      return Pos_p;
    return Pos_n;
  }
  TVector3 NearestPointLineToLine(const TVector3 &pos, const TVector3 &dir,
				  const TVector3 &x2, const TVector3 &a2,
				  const double &dl)
  {
    double dist;
    TVector3 xest, next;
    TVector3 tmpx2(x2);
    TVector3 tmpa2(a2);
    tmpx2.SetZ(0);
    tmpa2.SetZ(0);
    LineToLine(pos,dir,tmpx2,tmpa2,dl,dist,xest,next);
    return xest;
  }

  bool LineToCircle(const double &a,const double &b, const double &c, 
		    const double &rho ,const double &xc, const double &yc, 
		    double &x_p,double &y_p,double &x_n,double &y_n)
  {
    // in XY 2 dim 
    // line (wire) : a * x + b * y + c =0
    // circle : (x-xc)^2 + (y-yc)^2 = rho^2 
    // x_p,x_n : 2 soultions
  
    //   double dis=fabs(a*xc+b*yc+c)/sqrt(a*a+b*b);
    //   if(dis>rho){
    // #if DEBUG
    //     std::cout<<"no crossing point! rho dis "<<rho<<" "<<dis<<" "<<c<<std::endl;
    // #endif     
    //     // should be modified for short track analysis
    //     // if the hit is at the edge of the wire, it may be killed.
    //     return false;
    //   }

    double p[3];
    if(TMath::Abs(b)>0.001)
      {
	p[0]=a*a/(b*b)+1;
	p[1]=2*a*c/(b*b)-2*xc+2*yc*a/b;
	p[2]=xc*xc+c*c/(b*b)+2*yc*c/b+yc*yc-rho*rho;
      
	x_p=(-p[1]+sqrt(p[1]*p[1]-4*p[0]*p[2]) )/(2*p[0]);
	x_n=(-p[1]-sqrt(p[1]*p[1]-4*p[0]*p[2]) )/(2*p[0]);
	y_p= -(a*x_p+c)/b;
	y_n= -(a*x_n+c)/b;
      }
    else if(TMath::Abs(a)>0.001)
      {
	p[0]=b*b/(a*a)+1;
	p[1]=2*b*c/(a*a)-2*yc+2*xc*b/a;
	p[2]=yc*yc+c*c/(a*a)+2*xc*c/a+xc*xc-rho*rho;
      
	y_p=(-p[1]+sqrt(p[1]*p[1]-4*p[0]*p[2]) )/(2*p[0]);
	y_n=(-p[1]-sqrt(p[1]*p[1]-4*p[0]*p[2]) )/(2*p[0]);
	x_p= -(b*y_p+c)/a;
	x_n= -(b*y_n+c)/a;

      }
    double dis2=sqrt( (x_p-xc)*(x_p-xc)+(y_p-yc)*(y_p-yc) );
    if(fabs(dis2-rho)>100) {
      // std::cout<<"not macth rho! missing calc!"<<std::endl; 
      // std::cout<<a<<"  "<<b<<"  "<<c<<std::endl;
      // std::cout<<rho<<"  "<<xc<<"  "<<yc<<"  "<<x_p<<"  "<<y_p<<std::endl;
      return false;
    }
    return true;
  }


  double CalcHelixDCA( const double par1[5], const double par2[5], TVector3 &vtx1, TVector3 &vtx2, TVector3 &vtx )
  {
    const int NUM = 2;
    const int npoint = 100;
    const double initz = 0; //mm
    double region = 600; //+/-mm
    TVector3 pos1[npoint+1];
    TVector3 pos2[npoint+1];
    int num = 1;
    TVector3 now_vtx1;
    TVector3 now_vtx2;
    TVector3 now_vtx;
    double nowz = initz;
    double minl = 999.0;
    while(num<=NUM){
      //cerr<<"---"<<num<<endl;
      for(int i=0; i<npoint+1; i++){
	double z=nowz+(2*region/npoint)*(-npoint/2+i);
	//cerr<<z<<endl;
	pos1[i] = CalcHelixPosatZ(par1,z);
	pos2[i] = CalcHelixPosatZ(par2,z);
      }
      for(int i=0; i<npoint+1; i++){
	for(int j=0; j<npoint+1; j++){
	  TVector3 diff = pos1[i]-pos2[j];
	  double l = diff.Mag();
	  if(l<minl){
	    minl = l;
	    now_vtx1 = pos1[i];
	    now_vtx2 = pos2[j];
	    now_vtx = now_vtx1+now_vtx2;
	    now_vtx *=0.5;
	    nowz = now_vtx.z();
	  }
	}
      }
      num++;
      region = 2*region/10;
    }
    vtx1 = now_vtx1;
    vtx2 = now_vtx2;
    vtx = now_vtx;

    // fine search by k.t.
    double step1 = 0.5;
    double step2 = 0.5;
    double dl = minl;
    double dldiff = 9999;
    TVector3 v1[4],v2[4];
    double dlp[4];
    int counter = 0;
    while( 0.0001<dldiff && counter<10000 ){
      double z1 = vtx1.Z();
      double z2 = vtx2.Z();
      v1[0] = CalcHelixPosatZ( par1, z1+step1 );   v2[0] = CalcHelixPosatZ( par2, z2+step2 );  dlp[0] = (v1[0]-v2[0]).Mag();
      v1[1] = CalcHelixPosatZ( par1, z1+step1 );   v2[1] = CalcHelixPosatZ( par2, z2-step2 );  dlp[1] = (v1[1]-v2[1]).Mag();
      v1[2] = CalcHelixPosatZ( par1, z1-step1 );   v2[2] = CalcHelixPosatZ( par2, z2+step2 );  dlp[2] = (v1[2]-v2[2]).Mag();
      v1[3] = CalcHelixPosatZ( par1, z1-step1 );   v2[3] = CalcHelixPosatZ( par2, z2-step2 );  dlp[3] = (v1[3]-v2[3]).Mag();
      bool valnewed=false;
      for( int i=0; i<4; i++ ){
	if( dlp[i]<dl ){
	  vtx1 = v1[i]; vtx2 = v2[i];
	  vtx = (vtx1+vtx2)*0.5;
	  dldiff = fabs(dl-dlp[i]); dl = dlp[i]; 
	  valnewed=true;
	}
      }
      if( !valnewed ){
	step1 *= 0.5; step2 *= 0.5;
      }
      counter++;
    }
    if( counter==10000 ){
      std::cout << "!!! too many loops !!! " << std::endl;
    }
    return dl;
  }

  double dfunc_LTH(const TVector3 &lpos,const TVector3 &dline,const double &helixphi,const double *par)
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
    hpos=math::GetPosition(helixphi,par);

    // dhpos: momentum direction on helix
    dhpos.SetXYZ( (1./par[2])*sin(par[1]+helixphi)  ,
		  (-1./par[2])*cos(par[1]+helixphi) ,
		  (-1./par[2])*par[4]                );

    double k=hpos*dline-lpos*dline; // minimum distance point on the line
    double dk=dhpos*dline; // dk/dphi

    double S = 2.*(lpos+k*dline-hpos)*(dk*dline-dhpos); // d( (lpos+k*dline) -hpos)/d phi
    return S;
  }
  
  
  bool LineToHelix(const TVector3 &a, const  TVector3 &dline ,
		   const double *par, TVector3 &lnest,TVector3 &hnest,
		   double &dis)
  {
    // a: origin of the line
    // dline: direction of the line
    // par: helix parameter
    // lnest: nearest point on the line
    // hnest: nearest point on the helix
    TVector3 dlineu=dline.Unit();
    double phi=0,phi_b=0,phi_a=0;
    phi=math::CalcHelixPhi(a.x(),a.y(),par);

    double philen=1./par[2]*sqrt(1+par[4]*par[4]);
    double dist=1.;
    int trial=14; //sigma_position<1.2 micron 

    // find LTH initial param
    while(dist<1024)
      {
	phi_b=phi-dist/philen;
	phi_a=phi+dist/philen;
	double dlen_b=dfunc_LTH(a,dlineu,phi_b,par);
	double dlen_a=dfunc_LTH(a,dlineu,phi_a,par);
	if(dlen_b*dlen_a<=0) break;
	else {dist*=2;trial++;}
      }
  
    if(dist>=1024) 
      {
#if 0
	std::cout<<"Can not find LTH inital param!!"<<std::endl;
#endif
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
    hnest=math::GetPosition(phi,par);
    double k=(hnest-a)*dlineu;
    // line nearest point
    lnest=a+k*dlineu;
    TVector3 tmp;
    tmp=hnest-lnest;
    dis=tmp.Mag();
    return true;
  }
}

