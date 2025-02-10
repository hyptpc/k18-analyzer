// -*- C++ -*-

#include "RootData.hh"

//-----------------------------------------------//
// class DetectorHit
//-----------------------------------------------//
DetectorHit::DetectorHit()
{
  init();
}

void DetectorHit::init()
{
  nHitID      = DEFVALI;
  nDetectorID = DEFVALI;
  nLayerID    = DEFVALI;
  nChannelID  = DEFVALI;
  dAdc = DEFVALD;
  dTdc = DEFVALD;
  dDe = DEFVALD;
  dDt = DEFVALD;
  dDx = DEFVALD;
  tPos = DEFVALT;
}

CounterHit::CounterHit()
{
  init();
}

void CounterHit::init()
{
  t_charge   = 0.0;
  t_time     = 10000.0;
  t_track_id = DEFVALI;
  t_pid      = DEFVALI;
  t_pos      = DEFVALT;
  for (int i=0; i<20; i++){hit_particle[i] = DEFVALI;}
}

void CounterHit::Append(double ch, double t, int pid, int tid, TVector3 pos)
{
   t_charge += ch;
   if ( ch >0 && t < t_time ) {
     t_time     = t;
     t_pid      = pid ;
     t_track_id = tid ;
     t_pos      = pos ;
   }
}

double CounterHit::getPhi()
{
  double phi = atan2(t_pos.y(),t_pos.x());
  if (phi<0) return phi+2.*TMath::Pi();
  else       return phi;

}

void CounterHit::printStatus()
{
  std::cout << std::fixed << std::setw(5) << std::right << t_track_id << " ";
  std::cout << std::fixed << std::setw(5) << std::right << t_pid      << " ";
  std::cout << std::fixed << std::setw(7) << std::setprecision(2) << std::right << t_pos.x()   << " ";
  std::cout << std::fixed << std::setw(7) << std::setprecision(2) << std::right << t_pos.y()   << " ";
  std::cout << std::fixed << std::setw(7) << std::setprecision(2) << std::right << t_pos.z()   << " ";
  std::cout << std::fixed << std::setw(7) << std::setprecision(2) << std::right << t_charge    << " ";
  std::cout << std::fixed << std::setw(7) << std::setprecision(2) << std::right << t_time      << std::endl;
}

//-----------------------------------------------//
// class Track
//-----------------------------------------------//
Track::Track()
{
  init();
}

void Track::init()
{
  nTrackID       = DEFVALI;
  nParentTrackID = DEFVALI;
  nPdgID         = DEFVALI;
  fFlightLength  = DEFVALD;
  fFlightTime    = DEFVALD;
  fTime    = DEFVALD;
  tVertex        = DEFVALT;
  tMomentum      = DEFVALT;
  vnDetectorHitLink.clear();
}


//-----------------------------------------------//
// class RunHeaderMC
//-----------------------------------------------//
RunHeaderMC::RunHeaderMC()
{
  init();
}

void RunHeaderMC::init()
{
  nSeed         = DEFVALI;
  nNumEvent     = DEFVALI;
  nNumGenerated = DEFVALI;
  vnCounterList.clear();
}


//-----------------------------------------------//
// class EventHeaderMC
//-----------------------------------------------//
EventHeaderMC::EventHeaderMC()
{
  init();
}

void EventHeaderMC::init()
{
  nEventID     = DEFVALI;
}


//-----------------------------------------------//
// class DetectorData
//-----------------------------------------------//
DetectorData::DetectorData()
{
}


//-----------------------------------------------//
// class MCData
//-----------------------------------------------//
MCData::MCData()
{
}

//-----------------------------------------------//
// class Stop
//-----------------------------------------------//
Stop::Stop()
{
  posStop.SetXYZ(-999,-999,-999);
  volID=-1;
  nPDG         = DEFVALI;
}
//-----------------------------------------------//
// class StopData
//-----------------------------------------------//
StopData::StopData()
{
}

//-----------------------------------------------//
// class ReactionData
//-----------------------------------------------//
ReactionData::ReactionData()
{
  Init();
}

void ReactionData::Init()
{
  fCMOutParticleContainer.clear();
  fOutParticleContainer.clear();
  fIntParticleContainer.clear();
  fInitParticleContainer.clear();
  fReactionID = 0;
  fPDG.clear();
  fIntPDG.clear();
  fInitPDG.clear();
  for(int i=0; i<2; i++){ fNParticle[i] = 0; }
  for(int i=0; i<2; i++){ fFermiMom[i] = 0; }
  for(int i=0; i<8; i++){ fTmpVal[i] = 0; }
}

void ReactionData::SetCMParticle(Double_t px, Double_t py, Double_t pz, Double_t mass)
{
  TLorentzVector vec;
  vec.SetVectM( TVector3( px, py, pz), mass );
  fCMOutParticleContainer.push_back( vec );
}

void ReactionData::SetParticle(Double_t px, Double_t py, Double_t pz, Double_t mass)
{
  TLorentzVector vec;
  vec.SetVectM( TVector3( px, py, pz ), mass );
  fOutParticleContainer.push_back( vec );
}

void ReactionData::SetInitParticle(Double_t px, Double_t py, Double_t pz, Double_t mass)
{
  TLorentzVector vec;
  vec.SetVectM( TVector3( px, py, pz ), mass );
  fInitParticleContainer.push_back( vec );
}

void ReactionData::SetIntParticle(Double_t px, Double_t py, Double_t pz, Double_t mass)
{
  TLorentzVector vec;
  vec.SetVectM( TVector3( px, py, pz ), mass );
  fIntParticleContainer.push_back( vec );
}
