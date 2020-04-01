/**
 *  file: EventDisplay.hh
 *  date: 2017.04.10
 *
 */

#ifndef EVENT_DISPLAY_HH
#define EVENT_DISPLAY_HH

#include "ThreeVector.hh"
#include "DetectorID.hh"

class DCLocalTrack;

class TApplication;
class TArc;
class TBox;
class TBRIK;
class TCanvas;
class TFile;
class TGeometry;
class TGraph;
class TGraphErrors;
class TH1;
class TH2;
class TF1;
class TLatex;
class TMarker;
class TMarker3DBox;
class TMixture;
class TNode;
class TPad;
class TPave;
class TPaveLabel;
class TPaveText;
class TPolyLine;
class TPolyLine3D;
class TPolyMarker;
class TPolyMarker3D;
class TROOT;
class TRint;
class TRotMatrix;
class TStyle;
class TTRD1;
class TTRD2;
class TTUBS;
class TView;
class TLine;

typedef std::vector <TLine *> TLineContainer;

//______________________________________________________________________________
class EventDisplay //: public TObject
{
public:
  static EventDisplay&      GetInstance( void );
  static const std::string& ClassName( void );
  ~EventDisplay( void );

private:
  EventDisplay( void );
  EventDisplay( const EventDisplay& );
  EventDisplay& operator =( const EventDisplay& );

private:
  bool                       m_is_ready;
  TApplication              *m_theApp;
  TGeometry                 *m_geometry;
  TNode                     *m_node;
  TCanvas                   *m_canvas;
  TCanvas                   *m_canvas_vertex;
  TCanvas                   *m_canvas_hist;
  TCanvas                   *m_canvas_hist2;
  TCanvas                   *m_canvas_hist3;
  TCanvas                   *m_canvas_hist4;
  TCanvas                   *m_canvas_hist5;
  TCanvas                   *m_canvas_hist6;
  TCanvas                   *m_canvas_hist7;
  TCanvas                   *m_canvas_hist8;
  TCanvas                   *m_canvas_hist9;
  TCanvas                   *m_canvas_catch;
  TH1                       *m_hist_vertex_x;
  TH1                       *m_hist_vertex_y;
  TH1                       *m_hist_p;
  TH1                       *m_hist_m2;
  TH1                       *m_hist_missmass;
  TH2                       *m_hist_bh1;
  TH2                       *m_hist_bft;
  TH2                       *m_hist_bcIn;
  std::vector<TBox*>         m_BH1box_cont;
  std::vector<TLine*>        m_BcInTrack;
  TH2                       *m_hist_bft_p;
  TH2                       *m_hist_bcOut;
  std::vector<TBox*>         m_BH2box_cont;
  std::vector<TLine*>        m_BcOutTrack2;
  TH2                       *m_hist_bh2;
  TH2                       *m_hist_bcOut_sdcIn;
  TH2                       *m_hist_sdcIn_predict;
  TH2                       *m_hist_sdcIn_predict2;
  std::vector<TBox*>         m_BH2box_cont2;
  std::vector<TBox*>         m_SCHbox_cont;
  TBox                      *m_TargetXZ_box2;
  TBox                      *m_TargetYZ_box2;
  std::vector<TLine*>        m_BcOutTrack3;
  std::vector<TLine*>        m_SdcInTrack2;
  TH2                       *m_hist_sft_x;
  TH2                       *m_hist_sft_u;
  TH2                       *m_hist_sft_v;
  TH2                       *m_hist_sch;
  TH2                       *m_hist_tof;
  TH2                       *m_hist_sdc1;
  TH2                       *m_hist_sdc1p;
  TH2                       *m_hist_sdc2_l;
  TH2                       *m_hist_sdc2_t;
  TH2                       *m_hist_sdc2p_l;
  TH2                       *m_hist_sdc2p_t;
  TH2                       *m_hist_sdc2y_l;
  TH2                       *m_hist_sdc2y_t;
  TH2                       *m_hist_sdc2yp_l;
  TH2                       *m_hist_sdc2yp_t;
  TH2                       *m_hist_sdc3_l;
  TH2                       *m_hist_sdc3_t;
  TH2                       *m_hist_sdc3p_l;
  TH2                       *m_hist_sdc3p_t;
  TH2                       *m_hist_sdc3y_l;
  TH2                       *m_hist_sdc3y_t;
  TH2                       *m_hist_sdc3yp_l;
  TH2                       *m_hist_sdc3yp_t;

  TH2                       *m_hist_bc3;
  TH2                       *m_hist_bc3p;
  TH2                       *m_hist_bc3u;
  TH2                       *m_hist_bc3up;
  TH2                       *m_hist_bc3v;
  TH2                       *m_hist_bc3vp;
  TH2                       *m_hist_bc4;
  TH2                       *m_hist_bc4p;
  TH2                       *m_hist_bc4u;
  TH2                       *m_hist_bc4up;
  TH2                       *m_hist_bc4v;
  TH2                       *m_hist_bc4vp;
  TH1                       *m_hist_bc3_time;
  TH1                       *m_hist_bc3p_time;
  TH1                       *m_hist_bc4_time;
  TH1                       *m_hist_bc4p_time;

  TH2                       *m_hist_fbt1u;
  TH2                       *m_hist_fbt1up;
  TH2                       *m_hist_fbt1d;
  TH2                       *m_hist_fbt1dp;
  TH2                       *m_hist_fbt2u;
  TH2                       *m_hist_fbt2up;
  TH2                       *m_hist_fbt2d;
  TH2                       *m_hist_fbt2dp;

  TH2                       *m_hist_cft1_l;
  TH2                       *m_hist_cft1_t;
  TH2                       *m_hist_cft1_hi;
  TH2                       *m_hist_cft1_lo;
  TH2                       *m_hist_cft2_l;
  TH2                       *m_hist_cft2_t;
  TH2                       *m_hist_cft2_hi;
  TH2                       *m_hist_cft2_lo;
  TH2                       *m_hist_cft3_l;
  TH2                       *m_hist_cft3_t;
  TH2                       *m_hist_cft3_hi;
  TH2                       *m_hist_cft3_lo;
  TH2                       *m_hist_cft4_l;
  TH2                       *m_hist_cft4_t;
  TH2                       *m_hist_cft4_hi;
  TH2                       *m_hist_cft4_lo;
  TH2                       *m_hist_cft5_l;
  TH2                       *m_hist_cft5_t;
  TH2                       *m_hist_cft5_hi;
  TH2                       *m_hist_cft5_lo;
  TH2                       *m_hist_cft6_l;
  TH2                       *m_hist_cft6_t;
  TH2                       *m_hist_cft6_hi;
  TH2                       *m_hist_cft6_lo;
  TH2                       *m_hist_cft7_l;
  TH2                       *m_hist_cft7_t;
  TH2                       *m_hist_cft7_hi;
  TH2                       *m_hist_cft7_lo;
  TH2                       *m_hist_cft8_l;
  TH2                       *m_hist_cft8_t;
  TH2                       *m_hist_cft8_hi;
  TH2                       *m_hist_cft8_lo;
  TH2                       *m_hist_bgo;
  TH2                       *m_hist_piid_l;
  TH2                       *m_hist_piid_t;

  TNode                     *m_target_node;
  TNode                     *m_kurama_inner_node;
  TNode                     *m_kurama_outer_node;
  TNode                     *m_BH2wall_node;
  std::vector<TNode*>        m_BH2seg_node;
  std::vector<TNode*>        m_BC3x1_node;
  std::vector<TNode*>        m_BC3x2_node;
  std::vector<TNode*>        m_BC3v1_node;
  std::vector<TNode*>        m_BC3v2_node;
  std::vector<TNode*>        m_BC3u1_node;
  std::vector<TNode*>        m_BC3u2_node;
  std::vector<TNode*>        m_BC4u1_node;
  std::vector<TNode*>        m_BC4u2_node;
  std::vector<TNode*>        m_BC4v1_node;
  std::vector<TNode*>        m_BC4v2_node;
  std::vector<TNode*>        m_BC4x1_node;
  std::vector<TNode*>        m_BC4x2_node;
  std::vector<TNode*>        m_SFTu_node;
  std::vector<TNode*>        m_SFTv_node;
  std::vector<TNode*>        m_SFTx_node;
  std::vector<TNode*>        m_SDC1v1_node;
  std::vector<TNode*>        m_SDC1v2_node;
  std::vector<TNode*>        m_SDC1x1_node;
  std::vector<TNode*>        m_SDC1x2_node;
  std::vector<TNode*>        m_SDC1u1_node;
  std::vector<TNode*>        m_SDC1u2_node;
  std::vector<TNode*>        m_SDC2x1_node;
  std::vector<TNode*>        m_SDC2x2_node;
  std::vector<TNode*>        m_SDC2y1_node;
  std::vector<TNode*>        m_SDC2y2_node;
  std::vector<TNode*>        m_SDC3y1_node;
  std::vector<TNode*>        m_SDC3y2_node;
  std::vector<TNode*>        m_SDC3x1_node;
  std::vector<TNode*>        m_SDC3x2_node;
  std::vector<TNode*>        m_SSD1y1_node;
  TNode                     *m_FBHwall_node;
  std::vector<TNode*>        m_FBHseg_node;
  TNode                     *m_SCHwall_node;
  std::vector<TNode*>        m_SCHseg_node;
  TNode                     *m_FBT1d1wall_node;
  std::vector<TNode*>        m_FBT1d1seg_node;
  TNode                     *m_FBT1u1wall_node;
  std::vector<TNode*>        m_FBT1u1seg_node;
  TNode                     *m_FBT1d2wall_node;
  std::vector<TNode*>        m_FBT1d2seg_node;
  TNode                     *m_FBT1u2wall_node;
  std::vector<TNode*>        m_FBT1u2seg_node;
  TNode                     *m_FBT2d1wall_node;
  std::vector<TNode*>        m_FBT2d1seg_node;
  TNode                     *m_FBT2u1wall_node;
  std::vector<TNode*>        m_FBT2u1seg_node;
  TNode                     *m_FBT2d2wall_node;
  std::vector<TNode*>        m_FBT2d2seg_node;
  TNode                     *m_FBT2u2wall_node;
  std::vector<TNode*>        m_FBT2u2seg_node;
  TNode                     *m_TOFwall_node;
  std::vector<TNode*>        m_TOFseg_node;
  TPolyMarker3D             *m_init_step_mark;
  std::vector<TPolyLine3D*>  m_BcOutTrack;
  std::vector<TPolyLine3D*>  m_SdcInTrack;
  std::vector<TPolyLine3D*>  m_SdcOutTrack;
  TPolyMarker3D             *m_kurama_step_mark;
  // vertex
  TBox                      *m_TargetXZ_box;
  TBox                      *m_TargetYZ_box;
  TBox                      *m_EmulsionXZ_box;
  TBox                      *m_EmulsionYZ_box;
  std::vector<TBox*>         m_FBHseg_box;
  TH1                       *m_SSD_x_hist;
  TH1                       *m_SSD_y_hist;
  TMarker                   *m_VertexPointXZ;
  TMarker                   *m_VertexPointYZ;
  std::vector<TPolyLine*>    m_BcOutXZ_line;
  std::vector<TPolyLine*>    m_BcOutYZ_line;
  std::vector<TPolyLine*>    m_SdcInXZ_line;
  std::vector<TPolyLine*>    m_SdcInYZ_line;
  TPolyMarker               *m_KuramaMarkVertexX;
  TPolyMarker               *m_KuramaMarkVertexY;
  TPolyLine                 *m_MissMomXZ_line;
  TPolyLine                 *m_MissMomYZ_line;
  //CATCH
  TArc                      *m_Tgt_Arc;
  TArc                      *m_CFRP_Arc;
  TH2                       *m_hbase_catch;
  TH2                       *m_hbase_catch_zx;
  TH2                       *m_hbase_catch_zy;
  TGeometry                 *m_geometry_catch;
  TNode                     *m_node_catch;
  TCanvas                   *m_canvas_catch3d;
  std::vector<TArc*>         m_CFT_Arc_cont[NumOfPlaneCFT];
  TLineContainer             m_BGO_Line_cont[NumOfSegBGO];
  TLineContainer             m_PiID_Line_cont[NumOfSegPiID];
  std::vector<TNode*>        m_CFT_node_cont[NumOfPlaneCFT/2];
  std::vector<TPolyMarker3D*> m_CFT_UV_cont;
  std::vector<TNode*>        m_BGOseg_node_cont;
  std::vector<TNode*>        m_PiIDseg_node_cont;
  std::vector<TPolyLine3D*>  m_BcOutTrack_Catch_cont;
  std::vector<TPolyLine3D*>  m_SdcInTrack_Catch_cont;
  std::vector<TPolyLine3D*>  m_CFTTrack_cont;
  std::vector<TPolyLine*>    m_CFTTrack_xy_cont;
  std::vector<TPolyLine*>    m_CFTTrack_zx_cont;
  std::vector<TPolyLine*>    m_CFTTrack_zy_cont;
  std::vector<TPolyLine*>    m_BcOutTrack_Catch_xy_cont;
  std::vector<TPolyLine*>    m_BcOutTrack_Catch_zx_cont;
  std::vector<TPolyLine*>    m_BcOutTrack_Catch_zy_cont;
  std::vector<TPolyLine*>    m_SdcInTrack_Catch_xy_cont;
  std::vector<TPolyLine*>    m_SdcInTrack_Catch_zx_cont;
  std::vector<TPolyLine*>    m_SdcInTrack_Catch_zy_cont;
public:
  bool Initialize( void );
  bool IsReady( void ) const { return m_is_ready; }
  bool ConstructEmulsion( void );
  bool ConstructTarget( void );
  bool ConstructBH2( void );
  bool ConstructKURAMA( void );
  bool ConstructCollimator( void );
  bool ConstructBcOut( void );
  bool ConstructSdcIn( void );
  bool ConstructSdcOut( void );
  bool ConstructSSD( void );
  bool ConstructFBH( void );
  bool ConstructSCH( void );
  bool ConstructFBT( void );
  bool ConstructTOF( void );
  bool ConstructCATCH(void);
  bool ConstructCFT(void);
  bool ConstructBGO(void);
  bool ConstructPiID(void);
  bool ConstructCATCH3d(void);
  void DrawInitTrack( int nStep, ThreeVector *StepPoint );
  void DrawInitTrack( void );
  void DrawHitWire( int lid, int hit_wire,
		    bool range_check=true, bool tdc_check=true );
  void DrawTrackWire( int lid, int hit_wire, int it );
  void DrawText( double xpos, double ypos, const std::string& arg );
  void DrawHitHodoscope( int lid, int seg, int Tu=1, int Td=1 );
  void DrawBcOutLocalTrack( DCLocalTrack *tp );
  void DrawBcOutLocalTrack( double x0, double y0, double u0, double v0 );
  void DrawSdcInLocalTrack( DCLocalTrack *tp );
  void DrawSdcOutLocalTrack( DCLocalTrack *tp );
  void DrawSsdHit( int lid, int seg, double de );
  void DrawVertex( const ThreeVector& vertex );
  void DrawKuramaTrack( int nStep, ThreeVector *StepPoint, int Polarity );
  void DrawTarget( void );
  void DrawMissingMomentum( const ThreeVector& mom,
			    const ThreeVector& pos );
  void DrawMomentum( double momentum );
  void DrawMassSquare( double mass_square );
  void DrawMissMass( double mass_square );
  void DrawBH1( int seg, int tdc );
  void SetCorrectTimeBH1( int seg, double de );
  void DrawBFT( int layer, int seg, int tdc );
  void SetCorrectTimeBFT( double pos );
  void DrawBcInTrack( double x0, double u0 );
  void DrawBH2( int seg, int tdc );
  void SetCorrectTimeBH2( int seg, double de );
  void SetCorrectTimeBcOut( int layer, double pos );
  void DrawBcOutTrack( double x0, double u0, double y0, double v0, bool flagGoodForTracking = true );
  void DrawSdcInTrack( double x0, double u0, double y0, double v0, bool flagKurama = false, bool flagBeam = false );
  void DrawSFT( int layer, int seg, int tdc );
  void DrawSFT_X( int seg, int tdc );
  void DrawSFT_U( int seg, int tdc );
  void DrawSFT_V( int seg, int tdc );
  void DrawSCH( int seg, int tdc );
  void SetCorrectTimeSCH( int seg, double de );
  void SetCorrectTimeSdcIn( int layer, double pos );
  void DrawFBT( int detId, int layer, int UorD, int seg, int tdc );
  void DrawCFT( int layer, int seg, int LorT, int tdc );
  void DrawCFT_Adc( int layer, int seg, int HorL, int adc );
  void DrawCFT_AdcCor( int layer, int seg, int HorL, int adc );
  void DrawPiID(int seg, int LorT, int tdc );
  void DrawBGO(int seg, int tdc );
  void SetBGOWaveformCanvas(int nhit );
  void DrawBGOWaveform(int nc, int ngraph, int seg, TGraphErrors* gr );
  void DrawBGOFitFunc(int nc, int ngraph, int seg, TF1* func );
  
  void DrawTOF( int seg, int tdc );
  void DrawBcOutHit(int layer,  int wire, int tdc );
  void DrawBC3( int wire, int tdc );
  void DrawBC3p( int wire, int tdc );
  void DrawBC4( int wire, int tdc );
  void DrawBC4p( int wire, int tdc );
  void DrawSDC1( int wire, int tdc );
  void DrawSDC1p( int wire, int tdc );
  void DrawSdcOutHit(int layer,  int wire, int LorT, int tdc );
  void DrawSDC2_Leading( int wire, int tdc );
  void DrawSDC2_Trailing( int wire, int tdc );
  void DrawSDC2p_Leading( int wire, int tdc );
  void DrawSDC2p_Trailing( int wire, int tdc );
  void DrawSDC3_Leading( int wire, int tdc );
  void DrawSDC3_Trailing( int wire, int tdc );
  void DrawSDC3p_Leading( int wire, int tdc );
  void DrawSDC3p_Trailing( int wire, int tdc );
  void ShowHitFiber(int layer, int segment, double pe);// const;
  void ShowHitBGO(int segment, double de) const;
  void ShowHitPiID(int segment);
  void DrawCFTLocalTrack( DCLocalTrack *tp );
  void UpdateHist( void );
  void UpdateCATCH( void );
  void EndOfEvent( void );
  void ResetVisibility( void );
  void FiberPosPhi(int layer, int seg, double *x, double *y) const;
  void BGOPos(int seg, double *x, double *y) const;
  void CalcRotMatrix( double TA, double RA1, double RA2, double *rotMat );
  int  GetCommand( void ) const;
  void Run( bool flag=kTRUE );

private:
  void ResetVisibility( TNode *& node, Color_t c=kWhite );
  void ResetVisibility( std::vector<TNode*>& node, Color_t c=kWhite );
  void ResetHist( void );
  void ResetCATCH( void );
};

//______________________________________________________________________________
inline EventDisplay&
EventDisplay::GetInstance( void )
{
  static EventDisplay g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const std::string&
EventDisplay::ClassName( void )
{
  static std::string g_name("EventDisplay");
  return g_name;
}

#endif
