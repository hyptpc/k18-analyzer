// CDCWireMapMan.cpp

#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <new>
#include <cctype>
#include <map>
#include <utility>

#include "MathTools.hh"
#include "DetectorID.hh"
#include "CDCWireMapMan.hh"
#define OFFSET 1
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
namespace{  
  const double DEFVECT[3]={-999.,-999.,-999.};
  const int MAXCHAR = 144;  
  const int WIREMASK =0x00FF;  /* 8 (0-255) */
  const int LAYERMASK=0x001F;  /* 5 (0- 31) */
  const int WIRESHIFT = 0;
  const int LAYERSHIFT=10;
  const int ASDMASK  = 0x00FF;
  const int ASDSHIFT = 16;
}
inline int 
KEY(const int &layer,const int &wire)
{ 
  return ( ( (layer&LAYERMASK )<<LAYERSHIFT ) |
	   ( (wire&WIREMASK)<<WIRESHIFT ) );
}
inline int 
ASDKEY(const int &x, const int &y)
{ 
  return ( ( (x&ASDMASK)<<ASDSHIFT ) |
	   ( (y&ASDMASK) ) );
}
CDCWireMap::CDCWireMap()
{
  Layer = Wire = -1;
  Radius = Phi = Tilt = 0.0;
  for(int i=0;i<3;i++){
    Pos[i]= -999.;
    Posp[i]= -999.;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDCWireMap::CDCWireMap( int layer, int wire )
{
  Layer = layer; Wire = wire;
  Radius = Phi = Tilt = 0.0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDCWireMap::CDCWireMap( int layer, int wire,
			double rad, double phi, double tilt, double zlen )

{
  Layer = layer; Wire = wire;
  Radius = rad; Phi = phi; Tilt = tilt;

  Pos[0]=Radius*cos(Phi*TMath::DegToRad());
  Pos[1]=Radius*sin(Phi*TMath::DegToRad());
  Pos[2]=zlen/2.;

#if OFFSET
  Posp[0]=Radius*cos((Phi+Tilt)*TMath::DegToRad());
  Posp[1]=Radius*sin((Phi+Tilt)*TMath::DegToRad());
#else
  double S = zlen*tan(Tilt*TMath::DegToRad());
  double rad2 = Radius*Radius;
  double cost = 1.-S*S/(2.*rad2);
  double sint = S/(2.*rad2)*sqrt(4.*rad2-S*S);
  if( 0<TMath::Abs(tilt) ){
    Posp[0]=Pos[0]*cost-Pos[1]*sint;
    Posp[1]=Pos[0]*sint+Pos[1]*cost;  
  }else{
    Posp[0]=Pos[0];
    Posp[1]=Pos[1];
  }
#endif
  Posp[2]=-zlen/2.;
  //  std::cout<<layer<<"  "<<wire<<"  "<<Pos[0]<<"  "<<Pos[1]<<"  "<<Pos[2]<<std::endl;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDCWireMapMan::CDCWireMapMan()
  :m_isready(false)
{
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void CDCWireMapMan::Clear()
{
  GX = GY = GZ = dGX = dGY = dGZ = -999.;
  for( int i=0; i<NumOfCDCLayers; i++ ){
    Radius[i] = Phi0[i] = dPhi[i] = Tilt[i] = ZLengthOfWire[i] = 0.;
    NWires[i] = 0;
  }
  InnerRadius = OuterRadius = 0.;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::Initialize( const TString &filename1, const TString &filename2 )
{
  WireMapFileName=filename1;
  ASDMapFileName=filename2;
  return Initialize();
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::Initialize( const std::string &filename1, const std::string &filename2 )
{
  WireMapFileName=filename1;
  ASDMapFileName=filename2;
  return Initialize();
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::Initialize()
{
  static const TString funcname = "CDCWireMapMan::Initialize";
  // std::cout << "[" << funcname << "] Initialization start ..."<<std::endl;

  FILE *fp0, *fp1;
  if( (fp0=fopen(WireMapFileName.Data(), "r"))==0 ){
    std::cerr << funcname << " File open fail : " << WireMapFileName << std::endl;
    exit(-1);
  }
  if( (fp1=fopen(ASDMapFileName.Data(), "r"))==0 ){
    std::cerr << funcname << " File open fail : " << ASDMapFileName << std::endl;
    exit(-1);
  }
  
  int cid,n;
  double x,y,z,dx,dy,dz;
  double r1, r2,phi,zlen,dphi;
  int layer, nwires, slayer;
  double radius,tilt;
  char str[MAXCHAR];
  
  while( fgets(str, MAXCHAR, fp0)!=0 ){
    if( str[0]=='#' ) continue;
    if( (n=sscanf(str, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &cid, &layer, &x, &y, &z, &dx, &dy, &dz, &r1, &r2, &phi, &zlen, &dphi ))==13 ){
      if(cid!=DetIdCDC||layer!=-1) continue;
      GX=x; GY=y; GZ=z; 
      dGX=dx; dGY=dy; dGZ=dz; 
      InnerRadius=r1; OuterRadius=r2; 
      ZLengthOfWire[0]=zlen; 
      MotherVolume=(int)dphi;
      RotationAngle = dz; 
    }
    else if( (n=sscanf(str, "%d %d %d %lf %lf %lf %lf %lf %d", &cid, &layer, &nwires, &radius, &phi, &dphi, &tilt, &zlen, &slayer ))==9 ){
      if( layer<=NumOfCDCLayers ){
	NWires[layer] = nwires;
	Radius[layer] = radius;
	Phi0[layer] = phi;
	dPhi[layer] = dphi;
	Tilt[layer] = tilt;
	ZLengthOfWire[layer]=zlen; 
      }
    }else if( (n=sscanf(str, "%d %d %d %lf %lf %lf %lf", &cid, &layer, &nwires, &radius, &phi, &dphi, &tilt ))==7 ){
      if( layer<=NumOfCDCLayers ){
	NWires[layer] = nwires;
	Radius[layer] = radius;
	Phi0[layer] = phi;
	dPhi[layer] = dphi;
	Tilt[layer] = tilt;
	ZLengthOfWire[layer]=zlen; 
      }      
    }
  }
  fclose(fp0);
  
  CDCWireMap  *fMap = 0;
  wContainer.clear();
  asdContainer.clear();
  
  int wire,asdnum,asdch;
  while( fgets(str, MAXCHAR, fp1)!=0 ){
    if( str[0]=='#' ) continue;
    if( (n=sscanf(str, "%d %d %d %d", &layer, &wire, &asdnum, &asdch ))==4 ){
      asdContainer[ASDKEY(asdnum,asdch)] = ASDKEY(layer,wire);
    }
  }
  fclose(fp1);

  for( asdnum=0; asdnum<118; asdnum++)
    {
      for( asdch=0; asdch<16; asdch++ )
	{
	  if( GetLayerWire(asdnum,asdch,layer,wire) )
	    {    
	      radius = Radius[layer];
	      phi = Phi0[layer]+dPhi[layer]*wire;
	      tilt = Tilt[layer];	      
	      if( (fMap = new CDCWireMap(layer,wire,radius,phi,tilt,ZLengthOfWire[layer])) )
	      	{
		  unsigned int key= KEY(layer,wire);
	      	  wContainer[key] = *fMap;
	      	  delete fMap;
	      	}
	    }
	}
    }
  m_isready=true;
  return true;
}
  
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int CDCWireMapMan::GetSuperLayer( const int &layer ) const
{
  if(layer<0) return -1;
  else if(layer<3) return 1;
  else if(layer<5) return 2;
  else if(layer<7) return 3;
  else if(layer<9) return 4;
  else if(layer<11) return 5;
  else if(layer<13) return 6;
  else if(layer<15) return 7;
  return -1;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetWire( const int &layer, const int &wire,
			     int &slayer, double &rad, double &phi, double &tilt) const
{
  unsigned int key = KEY(layer,wire);
  std::map <unsigned int, CDCWireMap>::const_iterator wi = wContainer.find(key);
  const CDCWireMap *map=0;
  if( wi != wContainer.end() ){
    map = &(wi->second);
  }
  if(map){
    slayer = GetSuperLayer(layer);
    rad    = map->radius();
    phi    = map->phi();
    tilt   = map->tilt();
    return true;
  }
  else{
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TVector3 CDCWireMapMan::GetWirePos( const int &layer, const int &wire) const
{
  unsigned int key = KEY(layer,wire);
  std::map <unsigned int, CDCWireMap>::const_iterator wi = wContainer.find(key);
  const CDCWireMap *map=0;
  if( wi != wContainer.end() ){
    map = &(wi->second);
  }
  if(map){
    return map->pos();
  }
  else{
    return TVector3(DEFVECT);
  }
}

TVector3 CDCWireMapMan::GetWireDir( const int &layer, const int &wire) const
{
  unsigned int key = KEY(layer,wire);
  std::map <unsigned int, CDCWireMap>::const_iterator wi = wContainer.find(key);
  const CDCWireMap *map=0;
  if( wi != wContainer.end() ){
    map = &(wi->second);
  }
  if(map){
    return map->posp()-map->pos();
  }
  else{
    return TVector3(DEFVECT);
  }
}

TVector3 CDCWireMapMan::GetWirePosp( const int &layer, const int &wire) const
{
  unsigned int key = KEY(layer,wire);
  std::map <unsigned int, CDCWireMap>::const_iterator wi = wContainer.find(key);
  const CDCWireMap *map=0;
  if( wi != wContainer.end() ){
    map = &(wi->second);
  }
  if(map){
    return map->posp();
  }
  else{
    return TVector3(DEFVECT);
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetLayerWire( const int &asdnum, const int &asdch, int &layer, int &wire ) const
{
  unsigned int key = ASDKEY(asdnum,asdch);
  std::map <unsigned int, unsigned int>::const_iterator asdi = asdContainer.find(key);
  if( asdi != asdContainer.end() ){
    unsigned int lw = asdi->second;
    layer = (((lw)>>ASDSHIFT)&ASDMASK);
    wire  = ( (lw)           &ASDMASK);
    return true;
  }
  return false;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetWirePosDir( const int &layer, const int &wire , TVector3 &pos, TVector3 &dir) const
{
  unsigned int key = KEY(layer,wire);
  std::map <unsigned int, CDCWireMap>::const_iterator wi = wContainer.find(key);
  const CDCWireMap *map=0;
  if( wi != wContainer.end() ){
    map = &(wi->second);
  }

  if(map){
    pos=map->pos();
    dir=map->posp()-map->pos();
    return true;
  }
  else{
    return false;
  }
}    
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetGWirePosDir(  const int &layer, const int &wire , TVector3 &pos, TVector3 &dir) const
{
  unsigned int key = KEY(layer,wire);
  std::map <unsigned int, CDCWireMap>::const_iterator wi = wContainer.find(key);
  const CDCWireMap *map=0;
  if( wi != wContainer.end() ){
    map = &(wi->second);
  }
  if(map){
    TVector3 tmppos=map->posp();
    TVector3 tmpposp=map->pos();
    TVector3 posp;
    LocalToGlobal(tmpposp,pos);
    LocalToGlobal(tmppos ,posp);
    dir=posp-pos;
    return true;
  }
  else{
    return false;
  }
}    
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::PrintSimpleWireMap( std::ostream &p_out )
{
  int slayer;
  double rad, phi, tilt;
  for( int layer=0; layer<NumOfCDCLayers; layer++ ){
    for( int wire=0; wire<200; wire++ ){
      if( GetWire( layer, wire, slayer, rad, phi, tilt ) ){
	p_out.setf(std::ios::showpoint);
	p_out << std::setw(5) << DetIdCDC
	      << std::setw(5) << layer
	      << std::setw(5) << wire
	      << std::setw(5) << rad
	      << std::setw(5) << phi
	      << std::setw(5) << tilt
	      << std::endl;
      }
    }
  }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::PrintASDMap( std::ostream &p_out )
{
  int layer,wire;
  for( int asdnum=0; asdnum<118; asdnum++ ){
    for( int asdch=0; asdch<16; asdch++ ){
      if( GetLayerWire( asdnum,asdch,layer,wire ) ){
	p_out.setf(std::ios::showpoint);
	p_out << std::setw(5) << DetIdCDC
	      << std::setw(5) << layer
	      << std::setw(5) << wire
	      << std::setw(5) << asdnum
	      << std::setw(5) << asdch
	      << std::endl;
      }
    }
  }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetGeom( const int &layer, double &radius, double &phi0, double &dphi, double &tilt ) const
{
  if( 0<=layer && layer<NumOfCDCLayers ){
    radius = Radius[layer]; phi0 = Phi0[layer]; dphi = dPhi[layer]; tilt = Tilt[layer];
    return true;
  }
  else{
    radius = -999.; phi0 = -999.; dphi = -999.; tilt = -999.;
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetGPOS( double &gx, double &gy, double &gz,
			     double &dgx, double &dgy, double &dgz ) const
{
  gx = GX; gy = GY; gz = GZ; dgx = dGX; dgy = dGY; dgz = dGZ;
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetGParam( double *par ) const
{
  par[0] = GX; par[1] = GY; par[2] = GZ; par[3] = dGX; par[4] = dGY; par[5] = dGZ;
  par[6] = InnerRadius; par[7]=OuterRadius; par[8]=360; par[9]=ZLengthOfWire[0]; par[10]=MotherVolume; 
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GetFrame( double &zlen, double &rin, double &rout ) const
{
  zlen = ZLengthOfWire[0]; rin = InnerRadius; rout = OuterRadius;
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::HelixLocalToGlobal( double *local , double *global, double charge) const
{
  double tmp[5];
  for(int i=0;i<5;i++) tmp[i]=local[i];
  if(tmp[2]>0)         tmp[1]+=RotationAngle*TMath::DegToRad();
  else if(tmp[2]<0)    tmp[1]-=RotationAngle*TMath::DegToRad();
  for(int i=0;i<5;i++) global[i]=tmp[i];
  //  std::cout<<"GX,GY,GZ=  "<<GX<<"  "<<GY<<"  "<<GZ<<std::endl;
  math::ChangePivot(TVector3(0,0,0),TVector3(-GX,-GY,-GZ),tmp,global,charge);
  return true;
}  
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::LocalToGlobal( const TVector3 &local , TVector3 &global) const
{
  global=local;
  global.RotateZ(RotationAngle*TMath::DegToRad());
  global+=TVector3(GX,GY,GZ);
  return true;
}  
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCWireMapMan::GlobalToLocal( const TVector3 &global , TVector3 &local) const
{
  local=global-TVector3(GX,GY,GZ);;
  local.RotateZ(-RotationAngle*TMath::DegToRad());
  return true;
}  
