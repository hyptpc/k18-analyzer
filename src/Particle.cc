#include "Particle.hh"
#include "ELossTools.hh"
#include "GeomTools.hh"
#include "DatabasePDG.hh"

#define DEBUG 0
// ClassImp(Particle);
// ClassImp(pBeam);
// ClassImp(pNeutral);
// ClassImp(pCDS);

Particle::Particle() //: TObject()
{
  Clear();
}

void Particle::Clear(){
  ProductContainer.clear();
  BeamContainer.clear();
  CDSContainer.clear();
  NeutralContainer.clear();
  TrackIDContainer.clear();
}

void Particle::AddBeam( const pBeam &beam )
{
  BeamContainer.push_back(beam);
}

void Particle::AddProduct( const pCDS &pro )
{
  ProductContainer.push_back(pro);
}

void Particle::AddNeutral( const pNeutral &pro )
{
  NeutralContainer.push_back(pro);
}

void Particle::AddCDS( const pCDS &track )
{
  int pid=track.pid();  
  TrackIDContainer[ pid ].push_back(ncds());
  CDSContainer.push_back(track);
  return;
}

const int Particle::ncds(const int &pid) const
{
  std::map< int,std::vector<int> >::const_iterator is = TrackIDContainer.find(pid);
  if( is != TrackIDContainer.end()){
    return is->second.size();
  }
  return 0;  
}
  
pCDS* Particle::cds(const int &pid, const int &id)
{
  if(id<0) return 0;
  unsigned int i;
  std::map< int,std::vector<int> >::iterator is = TrackIDContainer.find(pid);
  if( is != TrackIDContainer.end()){
    i = (id<is->second.size()) ? is->second.at(id) : -1;
    return (0<=i&&i<CDSContainer.size()) ? &CDSContainer[i] : 0;
  }
  return 0;  
}

pBeam::pBeam( void ) 
  : //TObject(),
    PID(-1),T0seg(-1),BHDseg(-1),
    T0time(-999),T0ctime(-999),BHDX(-999),tofbhdt0(-999),Momentum(-999),
    BPCTRACK(false),BLCTRACK(false),MOM(false)
{
}

double pBeam::CalcVertexTime(const TVector3 &vertex)
{
  double momout,tof;
  eloss::CalcElossBeam(t0pos(),vertex,mom(),mass(),momout,tof);
  return t0ctime()+tof;
}
double pBeam::CalcVertexMom(const TVector3 &vertex, const int &sign)
{
  double momout,tof;
  eloss::CalcElossBeam(t0pos(),vertex,mom(),mass(),sign,momout,tof);
  return momout;
}

pNeutral::pNeutral() 
  : //TObject(),
    PID(-1),Beta(-999),TOF(-999),Mass(-999),
    CDHTime(-999),CDHDe(-999),FlightLength(-999)
{
  CDHseg.clear();
}

void pNeutral::Calc(pBeam *beam){
  double time_vertex=beam->CalcVertexTime(this->vertex()); 
  FlightLength=(this->cdhpos()-this->vertex()).Mag();
  TOF=CDHTime-time_vertex;
  double betan=fl()/tof()/(TMath::C()*1e-6);
  double momn= pdg::Mass(kNeutron)/sqrt(1/betan/betan-1); 
  TVector3 mom=(this->cdhpos()-this->vertex()).Unit()*momn;  
  SetBeta(betan);
  SetMomentum(mom);
  SetMass(pdg::Mass(kNeutron));
}

pCDS::pCDS() : //TObject(),
	       trackID(-1),daughterID1(-1),daughterID2(-1),
	       PID(-1),CombID(-1),
	       Momentum(-999),RawMomentum(-999),Beta(-999),Gamma(-999),
	       TOF(-999),Mass(-999),Mass2(-999),
	       VertexDistance(-999),
	       VertexBeamDistance(-999),
	       ProductBeamDCA(-999),
	       OpenAngle(-999),AngleLab(-999),FlightLength(-999),
	       Dt(-999),Chi(-999)
{
  CDHseg.clear();
}
