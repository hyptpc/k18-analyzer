// TrackTools.h
#ifndef TrackTools_h
#define TrackTools_h 1

#include <string>
#include <iostream>
#include <cmath>

#include "TVector3.h"
#include "TMath.h"

#include "CDCAnalyzer.hh"
#include "CDSTrack.hh"
#include "Particle.hh"
#include "LocalTrack.hh"
#include "HodoAnalyzer.hh"

namespace cds{
  bool FindMass2(CDSTrack *track, pBeam *beam,double tof ,
		 double &beta_calc,double &mass2,double &tmptof);
  bool FindMass2(CDSTrack *track, LocalTrack* bpc,double tof, double beammom, int pid_beam,
		 double &beta_calc,double &mass2,double &tmptof);
  bool FindMass2(CDSTrack *track, TVector3 t0pos,TVector3 vtxbpc, TVector3 vtxcdc, 
		 double tof, double beammom, double beammass,
		 double &beta_calc,double &mass2,double &tmptof);
  
  bool CalcVertex2HelixdE(CDSTrack *cds1, CDSTrack *cds2, const double &mass1, const double &mass2,
			  TVector3 &vtx1, TVector3 &vtx2);
  bool CalcVertex2Helix2(CDSTrack *cds1, CDSTrack *cds2, TVector3 &vtx1, TVector3 &vtx2);
  bool CalcVertex2Helix(CDSTrack *cds1, CDSTrack *cds2, TVector3 &vtx1, TVector3 &vtx2);
  bool CalcLineHelixVertex(LocalTrack *track,CDSTrack *cds, 
			   TVector3 &vtx1, TVector3 &vtx2, double &dis);
  bool CalcLineHelixVertex(pBeam *beam,CDSTrack *cds,
			   TVector3 &vtx1,TVector3 &vtx2,double &dis, bool ELOSS=false);
  void CalcBetaMass(TVector3 vertex,LocalTrack *beam,CDSTrack *cdc, 
		    int pid_beam,double tof, double &beta,double &mass2);
  
  int PID1d(double mom,double mass2);
  int PID2d(double mom,double mass2);
  double Mass2Resol(int pid, double mom, double fl);
  double MomResol(int pid, double mom);
  pCDS *CalcSingleAll(pBeam *beam, CDSTrack *cdc, HodoAnalyzer *hodoAna, bool ELOSS=false,bool REFIT=false);
  pCDS *Calc2HelixAll(pBeam *beam, pCDS* cds1, pCDS* cds2,CDCAnalyzer* cdcAna, bool ELOSS=false, int debug=0);

  void PDFLambda(double *per, double *pdf,bool YAMAGA=true);
}

class CDSAnalysis
{
};

#endif
