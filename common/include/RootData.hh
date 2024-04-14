#ifndef KnuclRootData_h
#define KnuclRootData_h 1

#include <vector>

#include <iostream>
#include <iomanip> 

#include "TObject.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "ComCrossSectionTable.hh"

#define DEFVALI -999
#define DEFVALD -999.0
#define DEFVALT TVector3(DEFVALD,DEFVALD,DEFVALD)

//-----------------------------------------------//
// class DetectorHit
//-----------------------------------------------//
class DetectorHit : public TObject
{
private:
  Int_t nHitID;
  Int_t nDetectorID;
  Int_t nLayerID;
  Int_t nChannelID;
  Double_t dAdc;
  Double_t dTdc;
  TVector3 tPos;
  TVector3 tMomentum;
  Int_t nTrackID;
  Double_t dTime;
  Double_t dDe;
  Double_t dDt;
  Double_t dDx;
  Int_t nPDG;
  Int_t nParentID;

public:
  DetectorHit();
  ~DetectorHit(){}

  void init();
  void setHitID(Int_t val){nHitID=val;}
  void setDetectorID(Int_t val){nDetectorID=val;}
  void setLayerID(Int_t val){nLayerID=val;}
  void setChannelID(Int_t val){nChannelID=val;}
  void setAdc(Double_t val){dAdc=val;}
  void setTdc(Double_t val){dTdc=val;}
  void setPos(const TVector3& val){tPos=val;}
  void setMomentum(const TVector3& val){tMomentum=val;}
  void setTrackID(Int_t val){nTrackID=val;}
  void setTime(Double_t val){dTime=val;}
  void setDe(Double_t val){dDe=val;}
  void setDt(Double_t val){dDt=val;}
  void setDx(Double_t val){dDx=val;}
  void setPDG(Int_t val){nPDG=val;}
  void setParentID(Int_t val){nParentID=val;}

  Int_t hitID() const {return nHitID;}
  Int_t detectorID() const {return nDetectorID;}
  Int_t layerID() const {return nLayerID;}
  Int_t channelID() const {return nChannelID;}
  Double_t adc() const {return dAdc;}
  Double_t tdc() const {return dTdc;}
  TVector3 pos() const {return tPos;}
  TVector3 momentum() const {return tMomentum;}
  Int_t trackID() const {return nTrackID;}
  Double_t time() const {return dTime;}
  Double_t de() const {return dDe;}
  Double_t dt() const {return dDt;}
  Double_t dx() const {return dDx;}
  Int_t pdg() const {return nPDG;}
  Int_t parentID() const {return nParentID;}

  ClassDef(DetectorHit,1);
};

class CounterHit : public TObject 
{
public:
  CounterHit();
  ~CounterHit(){}
  void init();
  void Append(double ch, double t, int pid, int tid, TVector3 pos);
  double getCharge(){return t_charge;     }
  double getTime()  {return t_time;       }
  double getPID()   {return t_pid;        }
  double getTID()   {return t_track_id;   }
  double getPhi();  
  double getX()     {return t_pos.x();}
  double getY()     {return t_pos.y();}
  double getZ()     {return t_pos.z();}
  void   printStatus();
private:
  double t_charge;
  double t_time;
  int    t_track_id;
  int    t_pid;
  TVector3 t_pos;
  int    hit_particle[20];

  ClassDef(CounterHit,1);
};




//-----------------------------------------------//
// class Track
//-----------------------------------------------//
class Track : public TObject
{
private:
  Int_t    nTrackID;
  Int_t    nParentTrackID;
  Int_t    nPdgID;
  TVector3 tVertex;
  TVector3 tMomentum;
  Double_t fFlightLength;
  Double_t fFlightTime;
  std::vector <Int_t>    vnDetectorHitLink;
  Double_t fTime;
  TString  CreatorProcess; 

public:
  Track();
  ~Track(){}

  void init();
  void setTrackID(Int_t val){nTrackID=val;}
  void setParentTrackID(Int_t val){nParentTrackID=val;}
  void setPdgID(Int_t val){nPdgID=val;}
  void setFlightLength(Double_t val){fFlightLength=val;}
  void setFlightTime(Double_t val){fFlightTime=val;}
  void setVertex(const TVector3& val){tVertex=val;}
  void setMomentum(const TVector3& val){tMomentum=val;}
  void setDetectorHitLink(const std::vector <Int_t>& val){vnDetectorHitLink=val;}
  void setDetectorHitLink(Int_t val){vnDetectorHitLink.push_back(val);}
  void setDetectorHitLink(Int_t i, Int_t val){vnDetectorHitLink[i]=val;}
  void setCreatorProcess(TString val){CreatorProcess=val;}
  void setTime(Double_t val){fTime=val;}

  Int_t trackID() const {return nTrackID;}
  Int_t parentTrackID() const {return nParentTrackID;}
  Int_t pdgID() const {return nPdgID;}
  TVector3 vertex() const {return tVertex;}
  TVector3 momentum() const {return tMomentum;}
  Double_t FlightLength() const {return fFlightLength;}
  Double_t FlightTime() const {return fFlightTime;}
  Int_t detectorHitLinkSize() const {return vnDetectorHitLink.size();}
  Int_t detectorHitLink(Int_t i) {return vnDetectorHitLink[i];}
  std::vector <Int_t>* detectorHitLink() {return &vnDetectorHitLink;}
  TString creatorProcess() const {return CreatorProcess; }
  Double_t time() const {return fTime; }

  ClassDef(Track,3);
};


//###############################################//


//-----------------------------------------------//
// class RunHeaderMC
//-----------------------------------------------//
class RunHeaderMC : public TObject
{
private:
  Int_t nSeed;
  Int_t nNumEvent;
  Int_t nNumGenerated;
  std::vector <Int_t> vnCounterList;
  CrossSectionTable fCStable;

public:
  RunHeaderMC();
  ~RunHeaderMC(){}

  void init();
  void setSeed(Int_t val){nSeed=val;}
  void setNumEvent(Int_t val){nNumEvent=val;}
  void setNumGenerated(Int_t val){nNumGenerated=val;}
  void setCounterList(std::vector <Int_t> val){vnCounterList=val;}
  void setCounterList(Int_t val){vnCounterList.push_back(val);}
  void setCStable(const CrossSectionTable& table){fCStable = table;}

  Int_t seed() const {return nSeed;}
  Int_t numEvent() const {return nNumEvent;}
  Int_t numGenerated() const {return nNumGenerated;}
  Int_t counterListSize() const {return vnCounterList.size();}
  Int_t counterList(Int_t i) const {return vnCounterList[i];}
  std::vector <Int_t> counterList() const {return vnCounterList;}
  CrossSectionTable CStable() const {return fCStable;}

  ClassDef(RunHeaderMC,1);
};

//-----------------------------------------------//
// class EventHeaderMC
//-----------------------------------------------//
class EventHeaderMC : public TObject
{
private:
  Int_t nEventID;

public:
  EventHeaderMC();
  ~EventHeaderMC(){}

  void init();
  void setEventID(Int_t val){nEventID=val;}
  Int_t eventID() const {return nEventID;}

  ClassDef(EventHeaderMC,1);
};


//-----------------------------------------------//
// class DetectorData
//-----------------------------------------------//
class DetectorData : public TObject
{
private:
  std::vector <DetectorHit> vcDetectorHit;

public:
  DetectorData();
  ~DetectorData(){}
  
  void setDetectorHit(const std::vector <DetectorHit>& val){vcDetectorHit=val;}
  void setDetectorHit(const DetectorHit& val){vcDetectorHit.push_back(val);}
  void setDetectorHit(Int_t i, const DetectorHit& val){vcDetectorHit[i]=val;}

  Int_t detectorHitSize() const {return vcDetectorHit.size();}
  DetectorHit* detectorHit(Int_t i){return &vcDetectorHit[i];}
  std::vector <DetectorHit>* detectorHit(){return &vcDetectorHit;}

  ClassDef(DetectorData,1);
};


//-----------------------------------------------//
// class MCData
//-----------------------------------------------//
class MCData : public TObject
{
private:
  std::vector <Track> vcTrack;

public:
  MCData();
  ~MCData(){}
  
  void setTrack(const std::vector <Track>& val){vcTrack=val;}
  void setTrack(const Track& val){vcTrack.push_back(val);}
  void setTrack(Int_t i, const Track& val){vcTrack[i]=val;}

  Int_t trackSize() const {return vcTrack.size();}
  Track* track(Int_t i){return &vcTrack[i];}
  std::vector <Track>* track(){return &vcTrack;}
 
  ClassDef(MCData,1);
};


class Stop : public TObject
{
 private:
  TVector3 posStop;
  TString  proName;
  TString  volName;
  TString  matName;
  Int_t  volID;
  Int_t nPDG;

public:
  Stop();
  ~Stop(){}
  
  void SetStopPos(const double x, const double y,const double z)
  { posStop.SetXYZ(x,y,z);}
  void SetStopPos(const TVector3 &pos)      { posStop=pos;}
  void SetProcessName(const TString &name)  { proName=name;}
  void SetVolumeName(const TString &name)   { volName=name;}
  void SetMaterialName(const TString &name) { matName=name;}
  void SetVolumeID(const Int_t &id)         { volID=id;}
  void setPDG(Int_t val)                    {nPDG=val;}
  
  TVector3 stoppos() const {return posStop;}
  TString proname() const {return proName;}
  TString volname() const {return volName;}
  TString matname() const {return matName;}
  Int_t volid() const {return volID;}
  Int_t pdg() const {return nPDG;}

  ClassDef(Stop,1);
};
//-----------------------------------------------//
// class StopData
//-----------------------------------------------//
class StopData : public TObject
{
 private:
  std::vector <Stop> vcStop;

 public:
  StopData();
  ~StopData(){}

  void setStop(const std::vector <Stop>& val){vcStop=val;}
  void setStop(const Stop& val){vcStop.push_back(val);}
  void setStop(Int_t i, const Stop& val){vcStop[i]=val;}
  Int_t trackSize() const {return vcStop.size();}
  Stop* track(Int_t i){return &vcStop[i];}
  std::vector <Stop>* track(){return &vcStop;}
  
  ClassDef(StopData,1);
};


//###############################################//


//-----------------------------------------------//
// class ReactionData
//-----------------------------------------------//
class ReactionData : public TObject
{
 private:
  typedef std::vector<TLorentzVector> TLorentzVectorContainer;
  TLorentzVectorContainer fCMOutParticleContainer;
  TLorentzVectorContainer fOutParticleContainer;
  TLorentzVectorContainer fIntParticleContainer;
  TLorentzVectorContainer fInitParticleContainer;
  
  int fReactionID;
  std::vector <int> fPDG;
  std::vector <int> fIntPDG;
  std::vector <int> fInitPDG;
  int fNParticle[2]; //0:nFin 1:nSpec
  
  double fFermiMom[2];
  double fTmpVal[8]; // temporal container for studies
  
 public:
  ReactionData();
  ~ReactionData(){}
  
  void Init();
  
  void SetCMParticle(Double_t px, Double_t py, Double_t pz, Double_t mass);
  void SetParticle(Double_t px, Double_t py, Double_t pz, Double_t mass);
  void SetIntParticle(Double_t px, Double_t py, Double_t pz, Double_t mass);
  void SetInitParticle(Double_t px, Double_t py, Double_t pz, Double_t mass);
  Int_t CMparticleSize(){ return fCMOutParticleContainer.size(); }
  Int_t ParticleSize(){ return fOutParticleContainer.size(); }
  Int_t IntParticleSize(){ return fIntParticleContainer.size(); }
  Int_t InitParticleSize(){ return fInitParticleContainer.size(); }
  const TLorentzVector& GetCMParticle(Int_t i){ return fCMOutParticleContainer[i]; }
  const TLorentzVector& GetParticle(Int_t i){ return fOutParticleContainer[i]; }
  const TLorentzVector& GetIntParticle(Int_t i){ return fIntParticleContainer[i]; }
  const TLorentzVector& GetInitParticle(Int_t i){ return fInitParticleContainer[i]; }

  void SetReactionID(int val){ fReactionID = val; }
  void SetPDG(int val){ fPDG.push_back(val); }
  void SetIntPDG(int val){ fIntPDG.push_back(val); }
  void SetInitPDG(int val){ fInitPDG.push_back(val); }
  void SetNParticle(int a, int b){ fNParticle[0] = a; fNParticle[1] = b; }
  int  ReactionID() const { return fReactionID; }
  int  PDGSize() const { return fPDG.size(); }
  int  IntPDGSize() const { return fIntPDG.size(); }
  int  InitPDGSize() const { return fInitPDG.size(); }
  const std::vector <int>& PDG() { return fPDG; }
  const std::vector <int>& IntPDG() { return fIntPDG; }
  const std::vector <int>& InitPDG() { return fInitPDG; }
  int  PDG(int i) const { return fPDG[i]; }
  int  IntPDG(int i) const { return fIntPDG[i]; }
  int  InitPDG(int i) const { return fInitPDG[i]; }
  int  NParticle(int i) const { return fNParticle[i]; }

  void   SetFermiMom(int i, double val){ fFermiMom[i] = val; }
  void   SetTmpVal(int i, double val){ fTmpVal[i] = val; }
  double FermiMom(int i) const { return fFermiMom[i]; }
  double TmpVal(int i) const { return fTmpVal[i]; }

  ClassDef(ReactionData,2)
};


#endif
