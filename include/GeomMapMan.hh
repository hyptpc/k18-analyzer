// GeomMapMan.h
#ifndef GeomMapMan_h
#define GeomMapMan_h 1

#include <map>
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"
#include "TPolyLine.h"

class GeomMap
{
 public:
  GeomMap();
  ~GeomMap() {};
 private:
  static const int npar=11;
  double R,Theta;
  double Pos[3];
  double Rot[3];
  double Size[4];
  double LightVelocity;
  
 public:
  inline void SetParam(const double *par,const double &scale=10);
  inline void GetParam(double *par) const;

  void SetPos( const TVector3 &pos ) { Pos[0]=pos.X(),Pos[1]=pos.Y(),Pos[2]=pos.Z(); }          
  void SetX( const double &x ) { Pos[0]=x; }
  void SetY( const double &y ) { Pos[1]=y; }
  void SetZ( const double &z ) { Pos[2]=z; }
  void SetRTheta( const double &r, const double &theta)
  { R=r; Theta=theta; Pos[0]=r*TMath::Cos(theta*TMath::DegToRad());  Pos[1]=r*TMath::Sin(theta*TMath::DegToRad()); }
  void SetRot( const TVector3 &rot ) { Rot[0]=rot.X(),Rot[1]=rot.Y(),Rot[2]=rot.Z(); }          
  void SetRotX( const double &x ) { Rot[0]=x; }
  void SetRotY( const double &y ) { Rot[1]=y; }
  void SetRotZ( const double &z ) { Rot[2]=z; }
  void SetSize( const double *size ) { for( int i =0; i<4;i++ ) Size[i]=size[i]; }
  void SetLength( const double &len ) { Size[1] = len; }
  void SetWidth( const double &w ) { Size[0] = w; }
  void SetThick( const double &t ) { Size[2] = t; }
  void SetLightVelocity( const double &lv ) { LightVelocity = lv; }

  inline TVector3 GetPos() const; 
  inline TVector3 GetPos(const double &ctsub) const;
  TVector3 GetRot() const { return TVector3(Rot); }
  const double* GetSize() const { return Size; }

  double GetX() const { return Pos[0]; }
  double GetY() const { return Pos[1]; }
  double GetZ() const { return Pos[2]; }
  double GetR() const { return R; }
  double GetTheta() const { return Theta; }

  double GetRotX() const { return Rot[0]; }
  double GetRotY() const { return Rot[1]; }
  double GetRotZ() const { return Rot[2]; }

  double GetLength() const { return IsBox() ? Size[1] : Size[3]; }
  double GetWidth() const { return IsBox() ? Size[0] : Size[0]*Size[2]*TMath::DegToRad(); }
  double GetThick() const { return IsBox() ? Size[2] : Size[1]-Size[0]; }

  double GetRmin() const { return Size[0]; }
  double GetRmax() const { return Size[1]; }
  double GetRmean() const { return (Size[0]+Size[0])/2.; }
  double GetdPhi() { return Size[2]; }

  bool IsBox() const { return Size[3]>0 ? false : true; }
  bool IsTube() const { return Size[3]>0 ? true : false; }
  bool IsCartesian() const { return Size[3]<0 ? false : true; }
  bool IsCyl() const { return Size[3]<0 ? true : false; }

  double GetLightVelocity() const { return LightVelocity; }
  void PrintMap( std::ostream &p_out = std::cout ) const;
  //  ClassDef( GeomMap, 1 );
};
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
inline void GeomMap::PrintMap( std::ostream &p_out )const {
  double tmppar[npar];
  GetParam(tmppar);
  for(int i=0;i<npar;i++)  
    p_out<<std::setw(10)<<"   "<<tmppar[i];
  p_out<<std::endl;
}

inline void GeomMap::SetParam( const double *param, const double &scale){
  for(int i=0;i<3;i++)
    Pos[i]=param[i]*scale;
  for(int i=0;i<3;i++)
    Rot[i]=param[i+3];
  for(int i=0;i<4;i++)
    Size[i]=param[i+6]*scale;
  if(IsCyl()) SetRTheta(param[0]*scale,param[1]);

  if(IsTube()){
    Size[2]=param[8];
    if(Size[2]<90){
      double tmpr=(param[0]+param[1])/2.*scale; double tmptheta=param[5];
      SetRTheta(tmpr,tmptheta);    
    }
  }
  LightVelocity=param[10];
}

inline TVector3 GeomMap::GetPos() const {
  if(IsCartesian())  return TVector3(Pos); 
  TVector3 tmppos(R,0,Pos[2]);
  tmppos.RotateZ(Theta*TMath::DegToRad());
  return tmppos;
}
inline TVector3 GeomMap::GetPos(const double &ctsub) const {
  if(IsCartesian()){    
    TVector3 dis(0,ctsub*LightVelocity,0);
    return GetPos()+dis;
  }else{
    TVector3 dis(0,0,ctsub*LightVelocity);
    // dis.RotateX(Rot[0]);
    // dis.RotateY(Rot[1]);
    // dis.RotateZ(Rot[2]);
    return GetPos()+dis;
  }
  return GetPos();
}

inline void GeomMap::GetParam( double *param ) const {
  for(int i=0;i<3;i++)
    param[i]=Pos[i];
  for(int i=0;i<3;i++)
    param[i+3]=Rot[i];
  for(int i=0;i<4;i++)
    param[i+6]=Size[i];
  if(IsCyl()){
    //  if(IsCyl()||(IsTube()&&Size[2]<90)){
    //    std::cout<<"cylindrical cordinate !!!"<<std::endl;
    param[0]=R; param[1]=Theta; 
  }
  param[10]=LightVelocity;
}

class GeomMapMan
{
 public:
  static GeomMapMan& GetInstance(void);
  static const std::string& ClassName( void );
  ~GeomMapMan();

  void SetFileNameCDS( const TString & filenameCDS );
  void SetFileNameBL( const TString & filenameBL );
  void SetFileNameHall( const TString & filenameHall );
  bool Initialize( const char *filename1, const char *filename2 );
  bool Initialize( const std::string& filename1, const std::string& filename2);
  bool Initialize( const char *filename1, const char *filename2, const char* filename3 );
  bool Initialize( const std::string& filename1, const std::string& filename2, const std::string& filename3 );
  bool Initialize();

 private:
  GeomMapMan();
  GeomMapMan(const GeomMapMan &right);

  static const int npar=11;
  TString FileNameBL;
  TString FileNameHall;
  TString FileNameCDS;

  bool m_isready;
  
  typedef std::map <unsigned int, GeomMap> GeomMapContainer;
  GeomMapContainer geomContainer;
  bool ReadFile(const TString filename);
  const GeomMap *GetMap(const int &cid, const int &seg) const;  
 
 public:
  bool GetParam( const int &cid, const int &seg, double *par) const;
  bool GetPos( const int &cid, const int &seg, TVector3 &pos);
  bool GetPos( const int &cid, const int &seg, const double &lv, TVector3 &pos);
  bool GetRot( const int &cid, const int &seg, TVector3 &rot);
  bool GetSize( const int &cid, const int &seg, const double *par );
  bool GetLightVelocity( const int &cid, const int &seg, double &lv);  
  bool GetGPos( const int &cid, const int &seg, TVector3 &pos) const;
  bool GetGSurface( const int &cid, const int &seg, TVector3 &pos);
  double GetGSurfaceZ( const int &cid, const int &seg);
  bool GetGPos( const int &cid, const int &seg, const double &lv, TVector3 &pos);
  //  bool GetGRot( const int &cid, const int &seg, TVector3 &rot);

  bool GetXYCDS( const int &cid, const int &seg, TPolyLine &pline);
  //  bool GetYZSlice( const int &cid, const int &seg, TPolyLine &pline);
  double GetZ( const int &cid, const int &seg);
  bool GetZX( const int &cid, const int &seg, TPolyLine &pline, const bool &GLOBAL);

  bool SetParam( const int &cid, const int &seg, const double *par);
  bool SetPos( const int &cid, const int &seg, const TVector3 &pos);
  bool SetRot( const int &cid, const int &seg, const TVector3 &rot);
  bool SetSize( const int &cid, const int &seg, const double *par);
  bool SetLightVelocity( const int &cid, const int &seg,const double &lv);

  TString GetFileNameCDS() { return FileNameCDS; }
  TString GetFileNameBL() { return FileNameBL; }
  TString GetFileNameHall() { return FileNameHall; }

  void PrintMap( const int &id,std::ostream &p_out = std::cout );
  void PrintMapBL(std::ostream &p_out = std::cout );
  void PrintMapCDS(std::ostream &p_out = std::cout );
};
inline GeomMapMan&
GeomMapMan::GetInstance( void )
{
  static GeomMapMan g_instance;
  return g_instance;
}
inline const std::string&
GeomMapMan::ClassName( void )
{
  static std::string g_name("GeomMapMan");
  return g_name;
}

#endif
