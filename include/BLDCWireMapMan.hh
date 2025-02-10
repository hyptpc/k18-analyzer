// BLDCWireMapMan.h

#ifndef BLDCWireMapMan_h
#define BLDCWireMapMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "TVector3.h"
#include "TMath.h"
#include "GeomMapMan.hh"

class BLDCWireMap
{
public:
  BLDCWireMap();
  ~BLDCWireMap() {};

 private:
  int nWire;
  double Z;
  int XY;
  double XY0, dXY, WireLength, TiltAngle, RotationAngle;
  double WirePhi;

 public:
  void SetParam( const int &nw, const double &z, const int &xy, const double &xy0, const double &dxy,
		 const double &wl, const double &tilt, const double &ra );
  void SetParam(double *param);

  void SetNWire( const int &nw ) { nWire = nw; }
  void SetZ(     const double &z ) { Z = z; }
  void SetXY( const int &xy ) { XY = xy; }
  void SetXY0( const double &xy0 ) { XY0 = xy0; }
  void SetdXY( const double &dxy ) { dXY = dxy; }
  void SetWireLength( const double &wl ) { WireLength = wl; }
  void SetTiltAngle( const double &tilt ) { TiltAngle = tilt; }
  void SetRotationAngle( const double &ra ) { RotationAngle = ra; }
  void SetWirePhi( const double &phi ) { WirePhi = phi; }

  int GetNWire()const  { return nWire; }
  double GetZ() const { return Z; }
  int GetXY()   const  { return XY; }
  double GetXY0() const { return XY0; } 
  double GetdXY() const { return dXY; }
  double GetWireLength() const { return WireLength; }
  double GetTiltAngle()const  { return TiltAngle; }
  double GetRotationAngle()const { return RotationAngle; }
  double GetWirePhi() const{ return WirePhi; }
};

class BLDCWireMapMan
{
public:
  static BLDCWireMapMan& GetInstance(void);
  static const std::string& ClassName( void );
  ~BLDCWireMapMan(void) {}
  void SetFileName( const TString & filename ) { FileName = filename; }
  bool Initialize();
  bool Initialize( const char *file_name );
  bool Initialize( const std::string& file_name );
  
private:
  BLDCWireMapMan( void );

 private:
  TString FileName;

  typedef std::map < unsigned int, BLDCWireMap > BLDCWireMapContainer;
  BLDCWireMapContainer bldcContainer;
  typedef std::map < unsigned int, GeomMap > GeomMapContainer;
  GeomMapContainer geomContainer;

  static const unsigned int KEYMASK = 0x000F;
  static const unsigned int CMASK   = 0x00FF;
  static const unsigned int LMASK   = 0x00FF;
  static const int          CSHIFT  = 4;
  static const int          LSHIFT  = 16;
  static const unsigned int KEYFLAG = 0x0003; 
  inline int KEY(const int &cid, const int &layer) const
  {  return ((((cid)&CMASK)<<CSHIFT) | (((layer)&LMASK)<<LSHIFT) | KEYFLAG ); }
  static const int MAXCHAR = 256;
  bool m_isready;

  
 public:
  bool IsReady( void ) const { return m_isready; }
  const GeomMap * GetGMap( const int &cid, const int &layer) const;
  const BLDCWireMap * GetWireMap( const int &cid, const int &layer) const;

  TVector3 CalcWirePosition( const int &cid, const int &layer, const int &wire) const;
  TVector3 CalcWireDirection( const int &cid, const int &layer, const int &wire) const;
  double GetTiltAngle( const int &cid, const int & layer) const;
  double GetRotationAngle( const int &cid, const int & layer) const;
  double GetLocalZ( const int &cid, const int &layer) const ;
  bool GetParam( const int &cid, const int &layer,
		 int &nw, double &z, int &xy, double &xy0, double &dxy, 
		 double &wl, double &tilt, double &ra ) const;
  bool GetXY0( const int &cid, const int &layer, double &xy0) const;
  int GetXY( const int &cid, const int &layer) const;
  // bool SetXY0( const int &cid, const int &layer,const double &xy0);
  // bool SetParam( const int &cid, const int &layer,
  // 		 const int &nw,const double &z,const int &xy,const double &xy0,const double &dxy, 
  // 		 const double &wl,const double &tilt,const double &ra );
  int  GetNWire( const int &cid, const int &layer ) const;
  bool GetGParam( const int &cid, TVector3 &pos, TVector3 &rot ) const;

  bool GetGParam( const int &cid, const int &seg, double *par) const;
  //  bool SetGParam( const int &cid, const int &seg, const double *par);

  TString GetFileName() { return FileName; }

  void PrintMap( const int &id,std::ostream &p_out = std::cout ) const ;
  void PrintMapBL(std::ostream &p_out = std::cout )const ;
};
//______________________________________________________________________________
inline BLDCWireMapMan&
BLDCWireMapMan::GetInstance( void )
{
  static BLDCWireMapMan g_instance;
  return g_instance;
}
inline const std::string&
BLDCWireMapMan::ClassName( void )
{
  static std::string g_name("BLDCWireMan");
  return g_name;
}
#endif
