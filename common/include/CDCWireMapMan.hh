// CDCWireMapMan.h

#ifndef CDCWireMapMan_h
#define CDCWireMapMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "TVector3.h"
#include "TString.h"
#include "TMath.h"

class CDCWireMap
{
 private:
  int Layer, Wire;
  int SuperLayer;
  double Radius, Phi, Tilt;
  double Pos[3],Posp[3];

 public:
  CDCWireMap();
  CDCWireMap( int layer , int wire );
  CDCWireMap( int layer , int wire, 
	      double rad, double phi, double tilt, double zlen );
  ~CDCWireMap() {}

 public:
  void SetLayer( const int &layer ) { Layer = layer; }
  void SetWire( const int &wire ) { Wire = wire; }
  void SetRadius( const double &rad ) { Radius = rad; }
  void SetPhi( const double &phi ) { Phi = phi; }
  void SetTilt( const double &tilt ) { Tilt = tilt; }

  int layer()       const { return Layer; }
  int wire()        const { return Wire; }
  double radius()   const { return Radius; }
  double phi()      const { return Phi; }
  double tilt()     const { return Tilt; }
  TVector3 pos()    const { return TVector3(Pos); }
  TVector3 posp()   const { return TVector3(Posp); }
};

class CDCWireMapMan
{
 public:
  static CDCWireMapMan& GetInstance( void );
  static const std::string& ClassName( void );
  ~CDCWireMapMan( void ){};

 private:
  CDCWireMapMan( void );
  TString WireMapFileName;
  TString ASDMapFileName;
  bool m_isready;

 public:
  bool IsReady( void ) const { return m_isready; }
  void SetWireMapFileName( const TString & filename ) { WireMapFileName = filename; }
  void SetASDMapFileName( const TString & filename ) { ASDMapFileName = filename; }
  bool Initialize( const TString &filename1,
		   const TString &filename2);
  bool Initialize();
  bool Initialize( const std::string& filename1,
		   const std::string& filename2
		   );

  void Clear();

 private:
  static const int NumOfCDCLayers=15;
  typedef std::map <unsigned int,  CDCWireMap> WireMapContainer;
  typedef std::map <unsigned int, unsigned int> ASDMapContainer;
  WireMapContainer wContainer;
  ASDMapContainer asdContainer;

  double GX,GY,GZ,dGX,dGY,dGZ;
  int    NWires[NumOfCDCLayers];
  double Radius[NumOfCDCLayers], Phi0[NumOfCDCLayers];
  double dPhi[NumOfCDCLayers], Tilt[NumOfCDCLayers];
  double ZLengthOfWire[NumOfCDCLayers];

  double RotationAngle;
  double InnerRadius, OuterRadius;
  int MotherVolume;

 public:
  TString GetFileName1() const { return WireMapFileName; }
  TString GetFileName2() const { return ASDMapFileName; }
  TString GetWireMapFileName() const { return WireMapFileName; }
  TString GetASDMapFileName() const { return ASDMapFileName; }

  int GetSuperLayer( const int &layer) const;
  bool GetWire( const int &layer, const int &wire,
		int &slayer, double &rad, double &phi, double &tilt) const;
  TVector3 GetWirePos( const int &layer, const int &wire) const;
  TVector3 GetWireDir( const int &layer, const int &wire) const;
  TVector3 GetWirePosp( const int &layer, const int &wire) const;

  bool GetLayerWire( const int &asdnum, const int &asdch, int &layer, int &wire) const;

  bool GetWirePosDir( const int &layer, const int &wire , TVector3 &pos, TVector3 &dir) const;
  bool GetGWirePosDir( const int &layer, const int &wire , TVector3 &pos, TVector3 &dir) const;

  bool GetGeom( const int &layer, double &radius, 
		double &phi0, double &dphi, double &tilt ) const ;
  bool GetGPOS( double &gx, double &gy, double &gz, 
		double &dgx, double &dgy, double &dgz ) const ;
  bool GetGParam( double *par ) const;
  bool GetFrame( double &zlen, double &rin, double &rout )const;

  bool HelixLocalToGlobal( double *local , double *global, double charge) const;  
  bool LocalToGlobal( const TVector3 &local , TVector3 &global) const;  
  bool GlobalToLocal( const TVector3 &global, TVector3 &local) const;  

  double gx() const { return GX; }
  double gy() const { return GY; }
  double gz() const { return GZ; }
  double dgx() const { return dGX; }
  double dgy() const { return dGY; }
  double dgz() const { return dGZ; }
  void gparam(double *par) const { par[0]=GX; par[1]=GY; par[2]=GZ; par[3]=RotationAngle*TMath::DegToRad(); }
/*   double thetaoffset() const { return Thetaoffset; } */
/*   double phioffset() const { return Phioffset; } */
  int    nwires( const int &layer ) const { return NWires[layer]; }
  double radius( const int &layer ) const { return Radius[layer]; }
  double phi0( const int &layer ) const { return Phi0[layer]; }
  double dphi( const int &layer ) const { return dPhi[layer]; }
  double tilt( const int &layer ) const { return Tilt[layer]; }
  int    nw( const int &layer ) const { return (int)(360/dPhi[layer]+0.2); }

  double zlen() const { return ZLengthOfWire[0]; }
  double rin()  const { return InnerRadius; }
  double rout() const { return OuterRadius; }
  double rot()  const { return RotationAngle; }

  bool PrintSimpleWireMap( std::ostream &p_out = std::cout );
  bool PrintASDMap( std::ostream &p_out = std::cout );
};
//______________________________________________________________________________
inline CDCWireMapMan&
CDCWireMapMan::GetInstance( void )
{
  static CDCWireMapMan g_instance;
  return g_instance;
}
inline const std::string&
CDCWireMapMan::ClassName( void )
{
  static std::string g_name("CDCWireMan");
  return g_name;
}

#endif
