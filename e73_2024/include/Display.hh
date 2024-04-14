#ifndef DISPLAY_h
#define DISPLAY_h 1

#include <vector>
#include <iostream>
#include <fstream>

#include <TSystem.h>
#include <TVirtualPad.h>

#include "HodoAnalyzer.hh"
#include "DCAnalyzer.hh"
#include "CDCAnalyzer.hh"
#include "CDSTrack.hh"

enum XYZ {kXY,kZX,kZY,kZPhi,kZR};

namespace disp
{
  // Tools 
  bool Wait();
  bool DrawPLine(const double &xc,const double &yc, 
		 const double &dx, const double &dy, 
		 const double &rot);

  // for hodoscopes (CDH/IH).
  bool DrawSegmentsXY(  TVirtualPad *pad, int cid );
  bool DrawSegmentsZY(  TVirtualPad *pad, int cid );
  bool DrawSegmentsZR(  TVirtualPad *pad, int cid );
  bool DrawSegmentsZX(  TVirtualPad *pad, int cid, bool GLOBAL=false );

  // local display for CDS
  bool DrawCDCLayersXY( TVirtualPad *pad );
  bool DrawCDCLayersZY( TVirtualPad *pad, bool ZR=false );
  bool DrawCDCLayersZX( TVirtualPad *pad );
  bool DrawCDCLayersZR( TVirtualPad *pad );

  bool DrawCDSHits( TVirtualPad *pad, HodoAnalyzer *hodo, int cid, enum XYZ xyz );
  inline bool DrawCDSHitsXY( TVirtualPad *pad, HodoAnalyzer *hodo, int cid )
  { return DrawCDSHits(pad,hodo,cid,kXY); }
  inline bool DrawCDSHitsZY( TVirtualPad *pad, HodoAnalyzer *hodo, int cid )
  { return DrawCDSHits(pad,hodo,cid,kZY); }
  inline bool DrawCDSHitsZX( TVirtualPad *pad, HodoAnalyzer *hodo, int cid )
  { return DrawCDSHits(pad,hodo,cid,kZX); }
  inline bool DrawCDSHitsZR( TVirtualPad *pad, HodoAnalyzer *hodo, int cid )
  { return DrawCDSHits(pad,hodo,cid,kZR); }
  inline bool DrawCDSHitsZPhi( TVirtualPad *pad, HodoAnalyzer *hodo, int cid )
  { return DrawCDSHits(pad,hodo,cid,kZPhi); }

  bool DrawCDCClusterHits(   TVirtualPad *pad, CDCAnalyzer *cdc, enum XYZ xyz );
  inline bool DrawCDCClusterHitsXY( TVirtualPad *pad, CDCAnalyzer *cdc )
  { return DrawCDCClusterHits(pad,cdc,kXY); }
  bool DrawCDCHits(   TVirtualPad *pad, CDCAnalyzer *cdc, enum XYZ xyz );
  inline bool DrawCDCHitsXY( TVirtualPad *pad, CDCAnalyzer *cdc )
  { return DrawCDCHits(pad,cdc,kXY); }
  inline bool DrawCDCHitsZX( TVirtualPad *pad, CDCAnalyzer *cdc )
  { return DrawCDCHits(pad,cdc,kZX); }
  inline bool DrawCDCHitsZY( TVirtualPad *pad, CDCAnalyzer *cdc )
  { return DrawCDCHits(pad,cdc,kZY); }
  inline bool DrawCDCHitsZR( TVirtualPad *pad, CDCAnalyzer *cdc )
  { return DrawCDCHits(pad,cdc,kZR); }
  inline bool DrawCDCHitsZPhi( TVirtualPad *pad, CDCAnalyzer *cdc )
  { return DrawCDCHits(pad,cdc,kZPhi); }

  bool DrawCDCTracks(   TVirtualPad *pad, CDCAnalyzer *cdc, enum XYZ xyz);
  inline bool DrawCDCTracksXY( TVirtualPad *pad, CDCAnalyzer *cdc)
  { return DrawCDCTracks(pad,cdc,kXY); }
  inline bool DrawCDCTracksZX( TVirtualPad *pad, CDCAnalyzer *cdc)
  { return DrawCDCTracks(pad,cdc,kZX); }
  inline bool DrawCDCTracksZY( TVirtualPad *pad, CDCAnalyzer *cdc)
  { return DrawCDCTracks(pad,cdc,kZY); }
  inline bool DrawCDCTracksZR( TVirtualPad *pad, CDCAnalyzer *cdc)
  { return DrawCDCTracks(pad,cdc,kZR); }
  inline bool DrawCDCTracksZPhi( TVirtualPad *pad, CDCAnalyzer *cdc)
  { return DrawCDCTracks(pad,cdc,kZPhi); }
  // inline bool DrawCDCTrackYZ( TVirtualPad *pad, CDCAnalyzer *cdc);
  // inline bool DrawCDCTrackXZ( TVirtualPad *pad, CDCAnalyzer *cdc);
#if 0
  bool DrawCDCTrackYZ( TVirtualPad *pad, CDSTrackingMan *trackMan );
  bool DrawCDCTrackXZ( TVirtualPad *pad, CDSTrackingMan *trackMan );

  // global display for BeamLine
  bool DrawBLHitXZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl, int cid );
  bool DrawBLHitYZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl, int cid );


  // display for Beamline DC
  bool DrawBLDCLayersXZ( TVirtualPad *pad, ConfMan *conf,int cid );
  bool DrawBLDCLayersYZ( TVirtualPad *pad, ConfMan *conf,int cid );
  bool DrawBLDCHit( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl, int cid,int xy );
  bool DrawBLDCTrack( TVirtualPad *pad, ConfMan *conf,BeamLineHitMan* blMan, BeamLineTrackMan *blTrackMan, int cid,int xy );

  bool DrawBLDCXZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
		   BeamLineTrackMan *track,int cid);
  bool DrawBLDCYZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
		   BeamLineTrackMan *track,int cid);
  bool DrawTrackBLDCXZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
		   BeamLineTrackMan *track,int cid);
  bool DrawTrackBLDCYZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
		   BeamLineTrackMan *track,int cid);

  bool DrawFDC(TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl );
  bool DrawFDCWire( TVirtualPad *pad, ConfMan *conf );
  bool DrawFDCHitWire( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl );

  bool DrawBLDCTrackHitXZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *bl, const int &cid);
  bool DrawBLDCTrackHitYZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bl, const int &cid);

  bool DrawBLDCTrackXZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *bl, const int &cid, const int &col=-1);
  bool DrawBLDCTrackYZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *bl, const int &cid, const int &col=-1);

  bool DrawClusterTime( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bl, const int &cid, const int &col=4);

  bool DrawBLC2TrackfromBLC1XZ( TVirtualPad *pad, ConfMan *conf, LocalTrack *track,const double &mom,const int &col);
  bool DrawBLC2TrackfromBLC1YZ( TVirtualPad *pad, ConfMan *conf, LocalTrack *track,const double &mom,const int &col);

#endif
};
#endif
