#ifndef PARTICLE_hh
#define PARTICLE_hh 1

#include <vector>
#include <map>
#include "TVector3.h"
#include "TLorentzVector.h"
// #include "GlobalVariables.h"
// #include "HodoscopeLikeHit.h"
// #include "FieldMapMan.h"
// #include "ConfMan.h"
//#include "TObject.h"
#include "DatabasePDG.hh"

class pBeam //: public TObject
{
 public:
  pBeam();
  virtual ~pBeam(){};

 private:
  int PID; //geant particle code
  //  int TrigID[20];
  int T0seg;
  int BHDseg;
  double T0time;
  double T0ctime;
  double BHDX;
  double tofbhdt0;
  double Momentum;

  bool BPCTRACK, BLCTRACK, MOM;

  TVector3 Vertex;
  TVector3 T0Pos;
  TVector3 BPCPos;//at FF
  TVector3 BPCDir;
  TVector3 BLCPos;//at FF
  TVector3 BLCDir;

 public:
  void SetPID(		const int &pid)		{ PID		= pid; }
  void SetT0seg(	const int &seg)		{ T0seg		= seg; }
  void SetT0Time(	const double &time)	{ T0time	= time; }
  void SetT0CTime(	const double &time)	{ T0ctime	= time; }
  void SetT0Pos(	const TVector3 &pos)	{ T0Pos		= pos; }
  void SetBHDseg(	const int &seg)		{ BHDseg	= seg; }
  void SetBHDT0TOF(	const double &tof)	{ tofbhdt0	= tof; }
  void SetBHDX(		const double &x )	{ BHDX		= x; }

  void SetVertex(TVector3 vec) { Vertex = vec; }  
  void SetBPCPos(TVector3 vec) { BPCPos = vec; BPCTRACK = true; }
  void SetBPCDir(TVector3 vec) { BPCDir = vec; BPCTRACK = true; }
  void SetBLCPos(TVector3 vec) { BLCPos = vec; BLCTRACK = true; }
  void SetBLCDir(TVector3 vec) { BLCDir = vec; BLCTRACK = true; }
    
  void SetMomentum(const double &mom) { Momentum = mom; MOM = true; } 

  int           pid()      const { return PID; }
  int		t0seg()    const { return T0seg; }
  double	t0time()   const { return T0time; }
  double	t0ctime()  const { return T0ctime; }
  TVector3	t0pos()    const { return T0Pos; }
  int		bhdseg()   const { return BHDseg; }
  double	bhdt0tof() const { return tofbhdt0; }
  double	bhdx()     const { return BHDX; }

  TVector3 beammom() const { return BPCDir.Unit()*Momentum; } 
  TVector3 bpcpos()  const { return BPCPos; }
  TVector3 bpcdir()  const { return BPCDir; }
  TVector3 blcpos()  const { return BLCPos; }
  TVector3 blcdir()  const { return BLCDir; }
  TVector3 vertex()  const { return Vertex; }

  TVector3 bpcpos(const double &z) const {
    return BPCPos +  BPCDir * (z / BPCDir.Z());
  }
  TVector3 blcpos(const double &z) const {
    return BLCPos +  BLCDir * (z / BLCDir.Z());
  }

  double mom()  const { return Momentum; }
  double mass() const { return pdg::Mass(pid()); }
  double beta() const { return mom()/sqrt(mass()*mass()+mom()*mom()); }
  bool isbpc()  const { return BPCTRACK; }
  bool isblc()  const { return BLCTRACK; }
  bool ismom()  const { return MOM; }

  double CalcVertexTime(const TVector3 &vertex);
  double CalcVertexMom( const TVector3 &vertex,const int &sign=-1);

  TLorentzVector GetLorentzVector() 
    { 
    TLorentzVector lv;
    lv.SetVectM(BPCDir.Unit()*Momentum,mass());
    return lv;
    }

  TLorentzVector GetLorentzVector(const TVector3 &vertex, const int &sign=-1) 
  { 
    TLorentzVector lv;
    lv.SetVectM(BPCDir.Unit()*CalcVertexMom(vertex,sign),mass());
    return lv;
  }
  
  //  ClassDef( pBeam, 1 );
};

class pCDS// : public TObject
{
 public:
  pCDS();
  virtual ~pCDS(){};

 private:
  int trackID;
  int daughterID1;
  int daughterID2;

  int PID;
  int CombID;

  double Momentum;
  double RawMomentum;
  double Beta;
  double Gamma;
  double TOF;
  double Mass;
  double Mass2;
  //  double PDGMass;
  double VertexDistance;
  double VertexBeamDistance;
  double ProductBeamDCA;
  double Param[5];

  double OpenAngle;
  double AngleLab;
  double FlightLength;
  double FlightTime;
  double Dt;
  double CDt;
  double Chi;

  double DeCDC;
  double DeCDH;

  TVector3 Vertex;
  TVector3 VertexCDC;
  TVector3 VertexBeam;
  TVector3 MomDir;

  std::vector<int> CDHseg;
  std::vector<int> IHseg;

  TLorentzVector lmom_daughter1;
  TLorentzVector lmom_daughter2;

 public:
  void SetTrackID(	const int &id)		{ trackID	= id; }
  void SetDaughterID1(	const int &id)		{ daughterID1	= id; }
  void SetDaughterID2(	const int &id)		{ daughterID2	= id; }

  void SetDaughterLMom1(const TLorentzVector &l){ lmom_daughter1	= l; }
  void SetDaughterLMom2(const TLorentzVector &l){ lmom_daughter2	= l; }

  void SetCombID(	const int &id)		{ CombID	= id; }
  void SetPID(		const int &pid)		{ PID		= pid; }
  void SetMomentum(	const double &mom)	{ Momentum	= mom; }
  void SetRawMomentum(	const double &mom)	{ RawMomentum	= mom; }
  void SetMass(		const double &mass)	{ Mass		= mass; }
  void SetMass2(	const double &mass)	{ Mass2		= mass; }
  void SetBeta(		const double &beta)	{ Beta		= beta; }
  void SetGamma(	const double &beta)	{ Gamma		= beta; }
  void SetTOF(		const double &tof)	{ TOF		= tof; }
  void SetFL(		const double &fl)	{ FlightLength	= fl; }
  void SetFT(		const double &ft)	{ FlightTime	= ft; }
  void SetDt(		const double &dt)	{ Dt		= dt; }
  void SetCDt(		const double &dt)	{ CDt		= dt; }
  void SetDt1(		const double &dt)	{ Dt		= dt; }
  void SetDt2(		const double &dt)	{ CDt		= dt; }
  void SetChi2(		const double &chi)	{ Chi		= chi; }
  void SetDeCDH(	const double &de)	{ DeCDH		= de; }
  void SetDeCDC(	const double &de)	{ DeCDC		= de; }
  void AddDeCDH(	const double &de)	{ DeCDH		+= de; }

  void SetVertex(     TVector3 vec) { Vertex     =vec; }
  void SetVertexCDC(  TVector3 vec) { VertexCDC  =vec; }
  void SetVertexBeam( TVector3 vec) { VertexBeam =vec; }
  void SetMomDir(     TVector3 vec) { MomDir     =vec; }
  void SetVDis(  const double &vdis ){ VertexDistance     =vdis; }
  void SetVBDis( const double &vdis ){ VertexBeamDistance =vdis; }
  void SetPBDCA( const double &vdis ){ ProductBeamDCA     =vdis; }

  void SetCDHSeg(const int &seg ) { CDHseg.push_back(seg); }

  void SetAngleLab( const double &ang ){ AngleLab =ang; }
  void SetOA(       const double &ang ){ OpenAngle=ang; }
  void SetParameters( double *par ) { for(int i=0;i<5;i++) Param[i] =par[i]; }

  int id() const { return trackID; }
  int daughter1() const { return daughterID1; }
  int daughter2() const { return daughterID2; }

  TLorentzVector get_lmom_daughter1() const { return lmom_daughter1; }
  TLorentzVector get_lmom_daughter2() const { return lmom_daughter2; }

  int    comb()		const { return CombID; }
  int    pid()		const { return PID; }
  int    pid1()		const { return daughterID1; }
  int    pid2()		const { return daughterID2; }
  double mom()		const { return Momentum; }
  double rawmom()	const { return RawMomentum; }
  double mass()		const { return Mass ; }
  double mass2()	const { return Mass2 ; }
  double pdgmass()	const { return pdg::Mass(pid()); }
  double beta()		const { return Beta; } 
  double gamma()        const { return Gamma; }
  double tof()		const { return TOF; };
  double fl()		const { return FlightLength; }
  double ft()		const { return FlightTime; }
  double dt()		const { return Dt; }
  double cdt()		const { return CDt; }
  double dt1()		const { return Dt; }
  double dt2()		const { return CDt; }
  double chi()		const { return Chi; }
  double decdh()	const { return DeCDH; }
  double decdc()	const { return DeCDC; }

  double vdis()		const { return VertexDistance; }
  double vbdis()	const { return VertexBeamDistance; }
  double pbdca()	const { return ProductBeamDCA; }
  double oa()		const { return OpenAngle; }
  double angle()		const { return AngleLab; }

  TVector3 vertex() { return Vertex     ; }
  TVector3 vcdc()   { return VertexCDC  ; }
  TVector3 vbeam()  { return VertexBeam ; }
  TVector3 momdir() { return MomDir.Unit(); }
  int    cdhseg(const int i) const {return CDHseg[i]; }
  bool   ncdh() const { return CDHseg.size(); }

  double* GetParameters() { return Param; }

  TLorentzVector GetLorentzVector() 
  {
    TLorentzVector lv;
    if(PID>=0)
      lv.SetVectM(momdir()*TMath::Abs(mom()),pdgmass());
    else
      lv.SetVectM(momdir()*TMath::Abs(mom()),mass());
    return lv;
  }

  //  ClassDef( pCDS, 1 );
};

class pNeutral //: public TObject
{
 public:
  pNeutral();
  virtual ~pNeutral(){};

 private:
  int PID;
  double Beta;
  double TOF;
  double Mass;

  double CDHTime;
  double CDHDe;
  double FlightLength;

  TVector3 Momentum;
  TVector3 CDHPos;
  TVector3 Vertex;
  TVector3 MomDir;

  std::vector<int> CDHseg;
 public:
  void SetPID(		const int &pid)		{ PID		= pid;		}
  void SetMomentum(	const TVector3 &mom)	{ Momentum	= mom;		}
  void SetMass(		const double &mass)	{ Mass		= mass;		}
  void SetBeta(		const double &beta)	{ Beta		= beta;		}
  void SetTOF(		const double &tof)	{ TOF		= tof;		}
  void SetFL(		const double &fl)	{ FlightLength	= fl;		}
  void SetCDHPos(	const TVector3 vec)	{ CDHPos	= vec;		}
  void SetVertex(	const TVector3 vec)	{ Vertex	= vec;		}
  void SetCDHSeg(	const int &seg )	{ CDHseg.push_back(seg);	}
  void SetCDHTime(	const double &v )	{ CDHTime	= v;		}
  void SetCDHDe(	const double &v )	{ CDHDe		= v;		}

  int      pid()	const { return PID; }
  TVector3 mom()	const { return Momentum; }
  double   mass()	const { return Mass ; }
  double   beta()	const { return Beta; }
  double   tof()	const { return TOF; };
  double   fl()		const { return FlightLength; }

  TVector3 cdhpos() const { return CDHPos; }
  TVector3 vertex() const { return Vertex; }
  int      cdhseg(const int i) const {return CDHseg[i]; }
  double   cdhtime()const { return CDHTime; }
  double   cdhde()  const { return CDHDe; }
  int      ncdh()   const { return CDHseg.size(); }

  TLorentzVector GetLorentzVector() 
  {
    TLorentzVector lv;
    lv.SetVectM(mom(),mass());
    return lv;
  }
  void Calc(pBeam* beam);
  //  ClassDef( pNeutral, 1 );
};

class Particle : public TObject
{
 public:
  Particle();
  virtual ~Particle(){};

 private:
  typedef std::vector<pBeam> pBeamcontainer;
  pBeamcontainer BeamContainer;
  
  typedef std::vector<pNeutral> pNeutralcontainer;
  pNeutralcontainer NeutralContainer;

  typedef std::vector<pCDS> pCDScontainer;
  pCDScontainer CDSContainer;
  pCDScontainer ProductContainer;
  
  std::map<int,std::vector<int>> TrackIDContainer;
  
public:
  const int nProduct() 	const { return ProductContainer.size(); }  
  const int ncds()	const { return CDSContainer.size(); }
  const int nBeam()	const { return BeamContainer.size(); }
  const int nNeutral()	const { return NeutralContainer.size(); }

  const int ncds( const int &pid ) const;
  const int nkaon()	const { return ncds( kKMinus); }
  const int npip()	const { return ncds( kPiPlus); }
  const int npim()	const { return ncds( kPiMinus); }
  const int ndeuteron()	const { return ncds( pdg::kDeuteron); }
  const int ntriton()	const { return ncds( pdg::kTriton); }
  const int nhe3()	const { return ncds( pdg::kHe3); }
  const int nproton()	const { return ncds( kProton); }
  const int nother()	const { return ncds( pdg::kOther); }

  //  pProduct *product(const int &i){ return &ProductContainer[i]; }
  pCDS     *product(const int &i) { return &ProductContainer[i]; }
  pBeam    *beam(   const int &i) { return &BeamContainer[i]; }
  pNeutral *neutral(const int &i) { return &NeutralContainer[i]; }

  pCDS *cds(		     const int &i) { return &CDSContainer[i]; }
  pCDS *cds( const int &pid, const int &i);
  pCDS *kaon(		     const int &i) { return cds( kKMinus   , i ); }
  pCDS *pip(		     const int &i) { return cds( kPiPlus   , i ); }
  pCDS *pim(		     const int &i) { return cds( kPiMinus  , i ); }
  pCDS *proton(		     const int &i) { return cds( kProton   , i ); }
  pCDS *deuteron(	     const int &i) { return cds( pdg::kDeuteron , i ); }
  pCDS *helium3(	     const int &i) { return cds( pdg::kHe3	, i ); }
  pCDS *triton(		     const int &i) { return cds( pdg::kTriton   , i ); }
  pCDS *other(		     const int &i) { return cds( pdg::kOther	, i ); }
  
  //  void AddProduct( const pProduct &pro);
  void AddProduct( const pCDS &pro);
  void AddCDS(     const pCDS &track);
  void AddNeutral( const pNeutral &track);
  void AddBeam(    const pBeam &beam);

  void Clear();
  void CalcAngleCM(const double &targetmass);

  // public:
  //  ClassDef( Particle, 1 );
};

#endif
